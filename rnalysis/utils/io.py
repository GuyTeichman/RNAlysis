import concurrent.futures
import json
import os
import warnings
from functools import lru_cache
from itertools import chain
from pathlib import Path
from typing import List, Set, Union, Iterable, Tuple, Dict, Any

import numpy as np
import pandas as pd
import requests

from rnalysis.utils import parsing, validation, ontology, __path__


def load_csv(filename: str, index_col: int = None, drop_columns: Union[str, List[str]] = False, squeeze=False,
             comment: str = None):
    """
    loads a csv df into a pandas dataframe.

    :type filename: str or pathlib.Path
    :param filename: name of the csv file to be loaded
    :type index_col: int, default None
    :param index_col: number of column to be used as index. default is None, meaning no column will be used as index.
    :type drop_columns: str, list of str, or False (default False)
    :param drop_columns: if a string or list of strings are specified, \
    the columns of the same name/s will be dropped from the loaded DataFrame.
    :type squeeze: bool, default False
    :param squeeze: If the parsed data only contains one column then return a Series.
    :type comment: str (optional)
    :param comment: Indicates remainder of line should not be parsed. \
    If found at the beginning of a line, the line will be ignored altogether. This parameter must be a single character.
    :return: a pandas dataframe of the csv file
    """
    assert isinstance(filename,
                      (str, Path)), f"Filename must be of type str or pathlib.Path, is instead {type(filename)}."
    encoding = 'ISO-8859-1'
    kwargs = dict(sep=None, engine='python', encoding=encoding, squeeze=squeeze, comment=comment, skipinitialspace=True)
    if index_col is not None:
        kwargs['index_col'] = index_col
    df = pd.read_csv(filename, **kwargs)
    df.index = [ind.strip() if isinstance(ind, str) else ind for ind in df.index]
    if isinstance(df, pd.DataFrame):
        df.columns = [col.strip() if isinstance(col, str) else col for col in df.columns]

        for col in df.columns:
            # check if the columns contains string data
            if pd.api.types.is_string_dtype(df[col]):
                df[col] = df[col].str.strip()
    else:
        if pd.api.types.is_string_dtype(df):
            df = df.str.strip()
    # if there remained only empty string "", change to Nan
    df = df.replace({"": np.nan})

    if drop_columns:
        drop_columns_lst = parsing.data_to_list(drop_columns)
        assert validation.isinstanceiter(drop_columns_lst,
                                         str), f"'drop_columns' must be str, list of str, or False; " \
                                               f"is instead {type(drop_columns)}."
        for col in drop_columns_lst:
            col_stripped = col.strip()
            if col_stripped in df:
                df.drop(col_stripped, axis=1, inplace=True)
            else:
                raise IndexError(f"The argument {col} in 'drop_columns' is not a column in the loaded csv file!")
    return df


def save_csv(df: pd.DataFrame, filename: str, suffix: str = None, index: bool = True):
    """
    save a pandas DataFrame to csv.

    :param df: pandas DataFrame to be saved
    :param filename: a string or pathlib.Path object stating the original name of the file
    :type suffix: str, default None
    :param suffix: A suffix to be added to the original name of the file. If None, no suffix will be added.
    :param index: if True, saves the DataFrame with the indices. If false, ignores the index.
    """
    fname = Path(filename)
    if suffix is None:
        suffix = ''
    else:
        assert isinstance(suffix, str), "'suffix' must be either str or None!"
    new_fname = os.path.join(fname.parent.absolute(), f"{fname.stem}{suffix}{fname.suffix}")
    df.to_csv(new_fname, header=True, index=index)


class GOlrAnnotationIterator:
    """
    A class that fetches GO annotations from the GOlr (Gene Ontology Solr) server. \
    This class can be used as an iterable.


    **Attributes**

    n_annotations: int
        The number of annotations found on the server that match the user's query.
    """
    __slots__ = {'taxon_id': 'NCBI Taxon ID for which to fetch GO Annotations',
                 'iter_size': 'number of annotations to be fetched per request',
                 'aspects': 'the GO Aspects for which GO Annotations should be fetched',
                 'qualifiers': 'the evidence types for which GO Annotations should be fetched',
                 'excluded_qualifiers': 'the evidence types for which GO Annotations should NOT be fetched',
                 'databases': 'the ontology databases from which GO Annotations should be fetched',
                 'excluded_databases': 'the ontology databases from which GO Annotations should NOT be fetched',
                 'evidence_types': 'the evidence types for which GO Annotations should be fetched',
                 'excluded_evidence_types': 'the evidence types for which GO Annotations should NOT be fetched',
                 'default_params': 'the default parameters for GET requests',
                 'n_annotations': 'number of annotations matching the filtering criteria'}
    URL = 'http://golr-aux.geneontology.io/solr/select?'

    _EXPERIMENTAL_EVIDENCE = {'EXP', 'IDA', 'IPI', 'IMP', 'IGI', 'IEP', 'HTP', 'HDA', 'HMP', 'HGI', 'HEP'}
    _PHYLOGENETIC_EVIDENCE = {'IBA', 'IBD', 'IKR', 'IRD'}
    _COMPUTATIONAL_EVIDENCE = {'ISS', 'ISO', 'ISA', 'ISM', 'IGC', 'RCA'}
    _AUTHOR_EVIDENCE = {'TAS', 'NAS'}
    _CURATOR_EVIDENCE = {'IC', 'ND'}
    _ELECTRONIC_EVIDENCE = {'IEA'}
    _EVIDENCE_TYPE_DICT = {'experimental': _EXPERIMENTAL_EVIDENCE, 'phylogenetic': _PHYLOGENETIC_EVIDENCE,
                           'computational': _COMPUTATIONAL_EVIDENCE, 'author': _AUTHOR_EVIDENCE,
                           'curator': _CURATOR_EVIDENCE, 'electronic': _ELECTRONIC_EVIDENCE}

    _ASPECTS_DICT = {'biological_process': 'P', 'molecular_function': 'F', 'cellular_component': 'C',
                     'biological process': 'P', 'molecular function': 'F', 'cellular component': 'C'}

    LEGAL_ASPECTS = {'P', 'F', 'C'}
    LEGAL_EVIDENCES = set.union(*[parsing.data_to_set(s) for s in _EVIDENCE_TYPE_DICT.values()])
    LEGAL_QUALIFIERS = {'not', 'contributes_to', 'colocalizes_with'}

    def __init__(self, taxon_id: int, aspects: Union[str, Iterable[str]] = 'any',
                 evidence_types: Union[str, Iterable[str]] = 'any',
                 excluded_evidence_types: Union[str, Iterable[str]] = None,
                 databases: Union[str, Iterable[str]] = 'any',
                 excluded_databases: Union[str, Iterable[str]] = None,
                 qualifiers: Union[str, Iterable[str]] = 'any',
                 excluded_qualifiers: Union[str, Iterable[str]] = None,
                 iter_size: int = 10000):
        """
        :param taxon_id: NCBI Taxon ID to fetch annotations for.
        :type taxon_id: int
        :param aspects: only annotations from the specified GO aspects will be included in the analysis. \
        Legal aspects are 'biological_process' (P), 'molecular_function' (F), and 'cellular_component' (C).
        :type aspects: str, Iterable of str, 'biological_process', 'molecular_function', 'cellular_component', \
        or 'any' (default='any')
        :param evidence_types: only annotations with the specified evidence types will be included in the analysis. \
        For a full list of legal evidence codes and evidence code categories see the GO Consortium website: \
        http://geneontology.org/docs/guide-go-evidence-codes/
        :type evidence_types: str, Iterable of str, 'experimental', 'phylogenetic' ,'computational', 'author', \
        'curator', 'electronic', or 'any' (default='any')
        :param excluded_evidence_types: annotations with the specified evidence types will be \
        excluded from the analysis. \
        For a full list of legal evidence codes and evidence code categories see the GO Consortium website: \
        http://geneontology.org/docs/guide-go-evidence-codes/
        :type excluded_evidence_types: str, Iterable of str, 'experimental', 'phylogenetic' ,'computational', \
        'author', 'curator', 'electronic', or None (default=None)
        :param databases: only annotations from the specified databases will be included in the analysis. \
        For a full list of legal databases see the GO Consortium website:
        http://amigo.geneontology.org/xrefs
        :type databases: str, Iterable of str, or 'any' (default)
        :param excluded_databases: annotations from the specified databases will be excluded from the analysis. \
        For a full list of legal databases see the GO Consortium website:
        http://amigo.geneontology.org/xrefs
        :type excluded_databases: str, Iterable of str, or None (default)
        :param qualifiers: only annotations with the speficied qualifiers will be included in the analysis. \
        Legal qualifiers are 'not', 'contributes_to', and/or 'colocalizes_with'.
        :type qualifiers: str, Iterable of str, or 'any' (default)
        :param excluded_qualifiers: annotations with the speficied qualifiers will be excluded from the analysis. \
        Legal qualifiers are 'not', 'contributes_to', and/or 'colocalizes_with'.
        :type excluded_qualifiers: str, Iterable of str, or None (default 'not')
        :type excluded_qualifiers: str, iterable of str, or None (default)
        :param iter_size: if the number of fetched annotations is larger than iter_size, the request will be \
        split into multiple requests of size iter_size.
        :type iter_size: int (default 10000)
        """
        self.taxon_id: int = taxon_id
        self.iter_size: int = iter_size
        # parse aspects
        self.aspects: Set[str] = self._parse_go_aspects(aspects)
        # parse qualifiers
        self.qualifiers: Set[str] = set() if qualifiers == 'any' else parsing.data_to_set(qualifiers)
        self.excluded_qualifiers: Set[str] = set() if excluded_qualifiers is None else \
            parsing.data_to_set(excluded_qualifiers)
        # parse databases
        self.databases: Set[str] = parsing.data_to_set(databases)
        self.excluded_databases: Set[str] = set() if excluded_databases is None else \
            parsing.data_to_set(excluded_databases)
        # parse evidence types
        self.evidence_types: Set[str] = self._parse_evidence_types(evidence_types)
        self.excluded_evidence_types: Set[str] = self._parse_evidence_types(excluded_evidence_types)
        # validate parameters
        self._validate_parameters()
        # generate default request parameter dictionary
        self.default_params: dict = {
            "q": "*:*",
            "wt": "json",  # return format
            "rows": 0,  # how many rows to return
            # how many annotations to fetch (fetch 0 to find n_annotations, then fetch in iter_size increments
            "start": None,  # from which annotation number to start fetching
            "fq": self._generate_query(),  # search query
            "fl": "source,bioentity_internal_id,annotation_class"  # fields
        }

        self.n_annotations = self._get_n_annotations()

    def _get_n_annotations(self) -> int:
        """
        Check and return the number of annotations on the GOlr server matching the user's query.
        """
        return json.loads(self._golr_request(self.default_params))['response']['numFound']

    def _validate_parameters(self):
        """
        Validate the type and legality of the user's inputs.
        """
        assert isinstance(self.taxon_id, int), f"'taxon_id' must be an integer. Instead got type {type(self.taxon_id)}."
        assert isinstance(self.iter_size, int), \
            f"'iter_size' must be an integer. Instead got type {type(self.iter_size)}."
        assert self.iter_size > 0, f"Invalid value for 'iter_size': {self.iter_size}."
        for field, legals in zip((self.aspects, chain(self.evidence_types, self.excluded_evidence_types),
                                  chain(self.qualifiers, self.excluded_qualifiers)),
                                 (self.LEGAL_ASPECTS, self.LEGAL_EVIDENCES, self.LEGAL_QUALIFIERS)):
            for item in field:
                assert item in legals, f"Illegal item {item}. Legal items are {legals}."

    @staticmethod
    def _golr_request(params: dict) -> str:
        """
        Run a get request to the GOlr server with the specified parameters, and return the server's text response.
        :param params: the get request's parameters.
        :type params: dict
        """
        req = requests.get(GOlrAnnotationIterator.URL, params=params)
        if not req.ok:
            req.raise_for_status()
        return req.text

    @staticmethod
    def _parse_evidence_types(evidence_types: Union[str, Iterable[str]]) -> Set[str]:
        """
        Parse the user's specified evidence types and excluded evidence types into a set of evidence type codes \
        which are supported by GOlr.
        :param evidence_types: evidence types to be parsed
        :type evidence_types: str, Iterable of str, 'experimental', 'phylogenetic' ,'computational', 'author', \
        'curator', 'electronic', or 'any'
        """
        if evidence_types == 'any':
            return set.union(*[parsing.data_to_set(s) for s in GOlrAnnotationIterator._EVIDENCE_TYPE_DICT.values()])

        elif evidence_types is None:
            return set()

        elif isinstance(evidence_types, str) and evidence_types.lower() in GOlrAnnotationIterator._EVIDENCE_TYPE_DICT:
            return parsing.data_to_set(GOlrAnnotationIterator._EVIDENCE_TYPE_DICT[evidence_types.lower()])

        elif validation.isiterable(evidence_types) and any(
            [isinstance(ev_type, str) and ev_type.lower() in GOlrAnnotationIterator._EVIDENCE_TYPE_DICT for ev_type in
             evidence_types]):
            return set.union(*[parsing.data_to_set(GOlrAnnotationIterator._EVIDENCE_TYPE_DICT[ev_type.lower()])
                               if ev_type.lower() in GOlrAnnotationIterator._EVIDENCE_TYPE_DICT else
                               parsing.data_to_set(ev_type) for ev_type in evidence_types])

        else:
            return parsing.data_to_set(evidence_types)

    @staticmethod
    def _parse_go_aspects(aspects: Union[str, Iterable[str]]) -> Set[str]:
        """
        Parse the user's specified GO aspects (namespaces) into a set of GO aspect codes \
        which are supported by GOlr.
        :param aspects: evidence types to be parsed
        :type aspects: str, Iterable of str, 'biological_process', 'molecular_function', 'cellular_component', \
        or 'any'
        """
        aspects = parsing.data_to_set(aspects)

        if aspects == {'any'}:
            return set.union(*[parsing.data_to_set(s) for s in GOlrAnnotationIterator._ASPECTS_DICT.values()])

        elif any(
            [isinstance(aspect, str) and aspect.lower() in GOlrAnnotationIterator._ASPECTS_DICT
             for aspect in aspects]):
            return {GOlrAnnotationIterator._ASPECTS_DICT[aspect.lower()]
                    if aspect.lower() in GOlrAnnotationIterator._ASPECTS_DICT else aspect for aspect in aspects}

        else:
            return aspects

    def _generate_query(self) -> List[str]:
        """
        Generate a Solr filter query (fq=...) to filter annotations based on the user's input.
        """
        # add fields with known legal inputs and cardinality >= 1 to query (taxon ID, aspect, evidence type)
        query = [f'document_category:"annotation"',
                 f'taxon:"NCBITaxon:{self.taxon_id}"',
                 ' OR '.join([f'aspect:"{aspect}"' for aspect in self.aspects]),
                 ' OR '.join([f'evidence_type:"{evidence_type}"' for evidence_type in self.evidence_types])]
        # exclude all 'excluded' items from query
        query.extend([f'-evidence_type:"{evidence_type}"' for evidence_type in self.excluded_evidence_types])
        query.extend([f'-source:"{db}"' for db in self.excluded_databases])
        query.extend([f'-qualifier:"{qual}"' for qual in self.excluded_qualifiers])
        # add union of all requested databases to query
        if not self.databases == {'any'}:
            query.append(' OR '.join(f'source:"{db}"' for db in self.databases))
        # add union of all requested qualifiers to query
        if len(self.qualifiers) > 0:
            query.append(' OR '.join([f'qualifier:"{qual}"' for qual in self.qualifiers]))
        return query

    def _annotation_generator_func(self):
        """
        Generator function that fetches all annotations from GOlr that match the user's input and yields them.
        """
        max_iters = int(np.ceil(self.n_annotations / self.iter_size))
        params = self.default_params.copy()
        params['omitHeader'] = "true"  # omit the header from the json response
        param_dicts_list = []
        start = 0

        for i in range(max_iters):
            params['start'] = start
            start += self.iter_size
            params['rows'] = self.iter_size if i <= max_iters - 1 else self.n_annotations % self.iter_size
            param_dicts_list.append(params.copy())

        processes = []
        with concurrent.futures.ThreadPoolExecutor() as executor:
            for param_dict in param_dicts_list:
                processes.append(executor.submit(self._golr_request, param_dict))
        for task in concurrent.futures.as_completed(processes):
            for record in json.loads(task.result())['response']['docs']:
                yield record

    def __iter__(self):
        return self._annotation_generator_func()


@lru_cache(maxsize=32, typed=False)
def _ensmbl_lookup_post_request(gene_ids: Tuple[str]) -> Dict[str, Dict[str, Any]]:
    """
    Perform an Ensembl 'lookup/id' POST request to find the species and database for several identifiers. \
    See full POST API at https://rest.ensembl.org/documentation/info/lookup_post

    :param gene_ids: a tuple of gene IDs to look up
    :type gene_ids: tuple of str
    :return: a dictionary with gene IDs as keys and dictionaries of attributes as values
    """
    url = 'https://rest.ensembl.org/lookup/id'
    headers = {"Content-Type": "application/json", "Accept": "application/json"}
    # split the gene IDs into chunks of 1000 (the maximum allowed POST request size)
    data_chunks = parsing.partition_list(gene_ids, 1000)
    output = {}
    for chunk in data_chunks:
        data = {"ids": parsing.data_to_list(chunk)}
        req = requests.post(url, headers=headers, data=data.__repr__().replace("'", '"'))
        if not req.ok:
            req.raise_for_status()
        output.update(req.json())
    return output


def infer_sources_from_gene_ids(gene_ids: Iterable[str]) -> Dict[str, Set[str]]:
    """
    Infer the
    :param gene_ids:
    :type gene_ids:
    :return:
    :rtype:
    """
    output = _ensmbl_lookup_post_request(parsing.data_to_tuple(gene_ids))
    sources = {}
    for gene_id in output:
        if output[gene_id] is not None:
            source = output[gene_id]['source']
            if source not in sources:
                sources[source] = set()
            sources[source].add(gene_id)
    return sources


def infer_taxon_from_gene_ids(gene_ids: Iterable[str]) -> Tuple[int, str]:
    """
    Infer the NCBI Taxon ID and name of a collection of gene IDs. \
    In cases where not all gene IDs map to the same taxon, the best-fitting taxon will be picked by a majority vote.

    :param gene_ids: a collection of gene IDs
    :type gene_ids: Iterable of str
    :return: a tuple of the best-matching taxon's NCBI Taxon ID and full scientific name.
    :rtype: Tuple[int ,str]
    """
    output = _ensmbl_lookup_post_request(parsing.data_to_tuple(gene_ids))
    species = dict()
    for gene_id in output:
        if output[gene_id] is not None:
            species[output[gene_id]['species']] = species.setdefault(output[gene_id]['species'], 0) + 1
    if len(species) == 0:
        raise ValueError("No taxon ID could be matched to any of the given gene IDs.")
    chosen_species = list(species.keys())[0]
    chosen_species_n = species[chosen_species]
    if len(species) > 1:
        warnings.warn(f"The given gene IDs match more than one species. "
                      f"Picking the species that fits the majority of gene IDs.")
        for s in species:
            if species[s] > chosen_species_n:
                chosen_species = s
                chosen_species_n = species[s]
    return map_taxon_id(chosen_species.replace('_', ' '))


@lru_cache(maxsize=32, typed=False)
def map_taxon_id(taxon_name: Union[str, int]) -> Tuple[int, str]:
    """
    Maps a given query (taxon name or NCBI Taxon ID) to the best-matching taxon from the NCBI taxonomy database. \
    Mapping is done through UniProt Taxonomy: https://www.uniprot.org/taxonomy/?

    :param taxon_name: a partial/full taxon name (str) or NCBI Taxon ID (int) to map
    :type taxon_name: int or str
    :return: a tuple of the best-matching taxon's NCBI Taxon ID and full scientific name.
    :rtype: Tuple[int ,str]
    """
    url = 'https://www.uniprot.org/taxonomy/?'

    params = {
        'format': 'tab',
        'query': taxon_name,
        'sort': 'score'
    }
    req = requests.get(url, params=params)
    if not req.ok:
        req.raise_for_status()
    res = req.text.splitlines()
    if len(res) == 0:
        raise ValueError(f"No taxons match the search query '{taxon_name}'.")
    header = res[0].split('\t')

    if isinstance(taxon_name, int):
        matched_taxon = None
        for line in res[1::]:
            split_line = line.split('\t')
            if int(split_line[0]) == taxon_name:
                matched_taxon = split_line
                break
        if matched_taxon is None:
            matched_taxon = res[1].split('\t')
    else:
        matched_taxon = res[1].split('\t')

    taxon_id = int(matched_taxon[header.index('Taxon')])
    scientific_name = matched_taxon[header.index('Scientific name')]
    if len(res) > 2 and not (taxon_name == taxon_id or taxon_name == scientific_name):
        warnings.warn(
            f"Found {len(res) - 1} taxons matching the search query '{taxon_name}'. "
            f"Picking the match with the highest score: {scientific_name} (taxonID {taxon_id}).")

    return taxon_id, scientific_name


class GeneIDTranslator:
    """
    A dictionary-like class used for mapping gene IDs from one type to another \
    (for example, from UniProtKB Accession to Entrez Gene ID), or from one type to itself.


    **Attributes**

    mapping_dict: dict or None
        The underlying dictionary that contains mapping from one gene ID type to another. \
        If mapping_dict is None, the GeneIDTranslator will automatically map any given gene ID to itself.
    """
    __slots__ = {'mapping_dict': 'dictionary mapping gene IDs from one type to another'}

    def __init__(self, mapping_dict: Union[dict, None] = None):
        """
        :param mapping_dict: a dictionary mapping gene IDs from one type to another. \
        If mappping_dict is None, gene IDs will be automatically mapped to themselves.
        :type mapping_dict: dict or None (default None)
        """
        self.mapping_dict = mapping_dict

    def __getitem__(self, key):
        if self.mapping_dict is None:
            return key
        return self.mapping_dict[key]

    def __contains__(self, item):
        try:
            _ = self[item]
            return True
        except KeyError:
            return False


def map_gene_ids(ids: Union[str, Iterable[str]], map_from: str, map_to: str = 'UniProtKB AC') -> GeneIDTranslator:
    """
    Map gene IDs from one identifier type to another using the UniProt ID Mapping service. \
    If some IDs cannot be mapped uniquely, duplicate mappings will be resolved by their UniProtKB Annotation Score. \
    Gene IDs that could not be mapped or were not recognized will be dropped from the output.

    :param ids: gene IDs to be mapped
    :type ids: str or an Iterable of strings
    :param map_from: identifier type to map from (for example 'UniProtKB AC' or 'WormBase')
    :type map_from: str
    :param map_to: identifier type to map to (for example 'UniProtKB AC' or 'WormBase'). \
    can be identical to 'map_from'
    :type map_to: str
    :return:a GeneIDTranslator object that uniquely maps each given gene ID in 'map_from' identifier type \
    to its matching gene ID in 'map_to' identifier type.
    :rtype: GeneIDTranslator
    """
    url = 'https://www.uniprot.org/uploadlists/'
    id_dict = _load_id_abbreviation_dict()
    validation.validate_uniprot_dataset_name(id_dict, map_to, map_from)
    ids = parsing.data_to_list(ids)
    n_queries = len(ids)
    # if map_from and map_to are the same, return an empty GeneIDTranslator (which will map any given gene ID to itself)
    if id_dict[map_to] == id_dict[map_from]:
        return GeneIDTranslator()
    # since the Uniprot service can only translate to or from 'UniProtKB AC' identifier type,
    # if we need to map gene IDs between two other identifier types,
    # then we will map from 'map_from' to 'UniProtKB AC' and then from 'UniProtKB AC' to 'map_to'.
    print(f"Mapping {n_queries} entries from '{map_from}' to '{map_to}'...")
    if id_dict[map_to] != id_dict['UniProtKB AC'] and id_dict[map_from] != id_dict['UniProtKB AC']:
        to_uniprot = map_gene_ids(ids, map_from, 'UniProtKB AC').mapping_dict
        from_uniprot = map_gene_ids(to_uniprot.values(), 'UniProtKB AC', map_to).mapping_dict
        output = {key: from_uniprot[val] for key, val in zip(to_uniprot.keys(), to_uniprot.values()) if
                  val in from_uniprot}
    else:
        output = {}
        duplicates = []
        # make sure that 'map_from' and 'map_to' are recognized identifier types
        if id_dict[map_to] != 'Null' and id_dict[map_from] != 'Null':
            # split ids into chunks to keep the Get Request size from getting too big
            for chunk in _format_ids_iter(ids):
                params = {
                    'from': id_dict[map_from],
                    'to': id_dict[map_to],
                    'format': 'tab',
                    'query': chunk,
                    # get 'annotation_score' data if possible (only when mapping to 'UniProtKB AC')
                    'columns': 'id,annotation_score' if id_dict[map_to] == id_dict['UniProtKB AC'] else 'id'}
                req = requests.get(url, params=params)
                if not req.ok:
                    req.raise_for_status()
                # if 'annotation_score' is available (only when mapping to 'UniProtKB AC'), parse it accordingly
                # using 'uniprot_tab_with_score_to_dict' which automatically keeps the best match for duplicates
                if id_dict[map_to] == id_dict['UniProtKB AC']:
                    output.update(parsing.uniprot_tab_with_score_to_dict(req.text))
                # if 'annotation_score' is not available, get a list of duplicate matches to process later
                else:
                    this_output, this_duplicates = parsing.uniprot_tab_to_dict(req.text)
                    output.update(this_output)
                    duplicates.extend(this_duplicates)
        # if there are unprocessed duplicates, map them in reverse (to 'UniProtKB AC' instead of from 'UniProtKB AC')
        # and then use the resulting 'annotation_score' to automatically keep the best match for duplicates
        if len(duplicates) > 0:
            for chunk in _format_ids_iter(duplicates):
                params = {
                    'from': id_dict[map_to],
                    'to': id_dict['UniProtKB AC'],
                    'format': 'tab',
                    'query': chunk,
                    'columns': 'id,annotation_score'}
                req = requests.get(url, params=params)
                if not req.ok:
                    req.raise_for_status()
                output.update(parsing.uniprot_tab_with_score_to_dict(req.text))

    if len(output) < n_queries:
        warnings.warn(f"Failed to map {n_queries - len(output)} entries from '{map_from}' to '{map_to}'. "
                      f"Returning {len(output)} successfully-mapped entries.")
    return GeneIDTranslator(output)


def _format_ids_iter(ids: Union[str, int, list, set], chunk_size: int = 500):
    if isinstance(ids, str):
        yield ids
    elif isinstance(ids, int):
        yield str(ids)
    else:
        for i in range(0, len(ids), chunk_size):
            j = min(chunk_size, len(ids) - i)
            yield " ".join((str(item) for item in ids[i:i + j]))


def _load_id_abbreviation_dict(dict_path: str = os.path.join(__path__[0], 'uniprot_dataset_abbreviation_dict.json')):
    with open(dict_path) as f:
        return json.load(f)


@lru_cache(maxsize=2)
def fetch_go_basic() -> ontology.DAGTree:
    """
    Fetches the basic Gene Ontology OBO file from the geneontology.org website ('go-basic.obo') and parses it into a \
    DAGTree data structure.
    :return: a parsed DAGTree for gene ontology propagation and visualization.
    :rtype: utils.ontology.DAGTree
    """
    url = 'http://current.geneontology.org/ontology/go-basic.obo'
    with requests.get(url, stream=True) as obo_stream:
        return ontology.DAGTree(obo_stream.iter_lines())
