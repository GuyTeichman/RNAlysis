import concurrent.futures
import os
import warnings
from pathlib import Path
from typing import List, Set, Union, Iterable, Tuple, Dict
import requests
import pandas as pd
import json
from itertools import chain
from functools import lru_cache

from rnalysis.utils import parsing, validation, __path__


def load_csv(filename: str, idx_col: int = None, drop_columns: Union[str, List[str]] = False, squeeze=False,
             comment: str = None):
    """
    loads a csv df into a pandas dataframe.

    :type filename: str or pathlib.Path
    :param filename: name of the csv file to be loaded
    :type idx_col: int, default None
    :param idx_col: number of column to be used as index. default is None, meaning no column will be used as index.
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
    if idx_col is not None:
        df = pd.read_csv(filename, index_col=idx_col, encoding=encoding, squeeze=squeeze, comment=comment)
    else:
        df = pd.read_csv(filename, encoding=encoding, squeeze=squeeze, comment=comment)
    if drop_columns:
        if isinstance(drop_columns, str):
            drop_columns = [drop_columns]
        assert isinstance(drop_columns,
                          list), f"'drop_columns' must be str, list, or False; is instead {type(drop_columns)}."
        for col in drop_columns:
            assert isinstance(col, str), f"'drop_columns' must contain strings only. " \
                                         f"Member {col} is of type {type(col)}."
            if col in df:
                df.drop(col, axis=1, inplace=True)
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


def golr_annotations_iterator(taxon_id: int, aspects: Union[str, Iterable[str]] = 'any',
                              evidence_types: Union[str, Iterable[str]] = 'any',
                              excluded_evidence_types: Union[str, Iterable[str]] = None,
                              databases: Union[str, Iterable[str]] = 'any',
                              excluded_databases: Union[str, Iterable[str]] = None,
                              qualifiers: Union[str, Iterable[str]] = 'any',
                              excluded_qualifiers: Union[str, Iterable[str]] = None,
                              iter_size: int = 10000):
    # define constants
    legal_aspects = {'P', 'F', 'C'}
    aspects_dict = {'biological_process': 'P', 'molecular_function': 'F', 'cellular_component': 'C',
                    'biological process': 'P', 'molecular function': 'F', 'cellular component': 'C'}
    experimental_evidence = {'EXP', 'IDA', 'IPI', 'IMP', 'IGI', 'IEP', 'HTP', 'HDA', 'HMP', 'HGI', 'HEP'}
    phylogenetic_evidence = {'IBA', 'IBD', 'IKR', 'IRD'}
    computational_evidence = {'ISS', 'ISO', 'ISA', 'ISM', 'IGC', 'RCA'}
    author_evidence = {'TAS', 'NAS'}
    curator_evidence = {'IC', 'ND'}
    electronic_evidence = {'IEA'}
    evidence_type_dict = {'experimental': experimental_evidence, 'phylogenetic': phylogenetic_evidence,
                          'computational': computational_evidence, 'author': author_evidence,
                          'curator': curator_evidence,
                          'electronic': electronic_evidence}
    legal_evidence = set.union(*[parsing.data_to_set(s) for s in evidence_type_dict.values()])
    legal_qualifiers = {'not', 'contributes_to', 'colocalizes_with'}
    # parse aspects, databases, qualifiers
    aspects = parsing.parse_go_aspects(aspects, aspects_dict)
    databases = parsing.data_to_set(databases)
    qualifiers = () if qualifiers == 'any' else parsing.data_to_set(qualifiers)
    excluded_qualifiers = set() if excluded_qualifiers is None else parsing.data_to_set(excluded_qualifiers)
    excluded_databases = set() if excluded_databases is None else parsing.data_to_set(excluded_databases)
    # parse evidence types
    evidence_types = parsing.parse_evidence_types(evidence_types, evidence_type_dict)
    excluded_evidence_types = parsing.parse_evidence_types(excluded_evidence_types, evidence_type_dict)
    # assert legality of inputs
    for field, legals in zip((aspects, chain(evidence_types, excluded_evidence_types),
                              chain(qualifiers, excluded_qualifiers)),
                             (legal_aspects, legal_evidence, legal_qualifiers)):
        for item in field:
            assert item in legals, f"Illegal item {item}. Legal items are {legals}."
    # add fields with known legal inputs and cardinality >= 1 to query (taxon ID, aspect, evidence type)
    query = [f'document_category:"annotation"',
             f'taxon:"NCBITaxon:{taxon_id}"',
             ' OR '.join([f'aspect:"{aspect}"' for aspect in aspects]),
             ' OR '.join([f'evidence_type:"{evidence_type}"' for evidence_type in evidence_types])]
    # exclude all 'excluded' items from query
    query.extend([f'-evidence_type:"{evidence_type}"' for evidence_type in excluded_evidence_types])
    query.extend([f'-source:"{db}"' for db in excluded_databases])
    query.extend([f'-qualifier:"{qual}"' for qual in excluded_qualifiers])
    # add union of all requested databases to query
    if not databases == {'any'}:
        query.append(' OR '.join(f'source:"{db}"' for db in databases))
    # add union of all requested qualifiers to query
    if len(qualifiers) > 0:
        query.append(' OR '.join([f'qualifier:"{qual}"' for qual in qualifiers]))

    params = {
        "q": "*:*",
        "wt": "json",  # return format
        "rows": 0,  # how many annotations to fetch (fetch 0 to find n_annotations, then fetch in iter_size increments
        "start": None,  # from which annotation number to start fetching
        "fq": query,  # search query
        "fl": "source,bioentity_internal_id,annotation_class"  # fields
    }
    # get number of annotations found in the query
    output = _golr_request(params)
    n_annotations = json.loads(output)['response']['numFound']
    print(f"Fetching {n_annotations} annotations...")

    # fetch all annotations in batches of size iter_size, and yield them one-by-one
    start = 0
    max_iters = n_annotations // iter_size + 1
    params['omitHeader'] = "true"  # omit the header from the json response
    param_dicts_list = []

    for i in range(max_iters):
        params['start'] = start
        start += iter_size
        params['rows'] = iter_size if i < max_iters - 1 else n_annotations % iter_size
        param_dicts_list.append({key: val for key, val in zip(params.keys(), params.values())})

    processes = []
    with concurrent.futures.ThreadPoolExecutor() as executor:
        for param_dict in param_dicts_list:
            processes.append(executor.submit(_golr_request, param_dict))
    for task in concurrent.futures.as_completed(processes):
        for record in json.loads(task.result())['response']['docs']:
            yield record


def _golr_request(params: dict) -> str:
    url = 'http://golr-aux.geneontology.io/solr/select?'
    req = requests.get(url, params=params)
    if not req.ok:
        req.raise_for_status()
    return req.text


@lru_cache(maxsize=32, typed=False)
def _ensmbl_lookup_post(gene_ids: Tuple[str]) -> dict:
    url = 'https://rest.ensembl.org/lookup/id'
    headers = {"Content-Type": "application/json", "Accept": "application/json"}
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
    output = _ensmbl_lookup_post(parsing.data_to_tuple(gene_ids))
    sources = {}
    for gene_id in output:
        if output[gene_id] is not None:
            source = output[gene_id]['source']
            if source not in sources:
                sources[source] = set()
            sources[source].add(gene_id)
    return sources


def infer_taxon_id_from_gene_ids(gene_ids: Iterable[str]) -> Tuple[int, str]:
    output = _ensmbl_lookup_post(parsing.data_to_tuple(gene_ids))
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
    Mapping is done through UniProt Taxonomy.

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
    url = 'https://www.uniprot.org/uploadlists/'
    id_dict = _load_id_abbreviation_dict()
    validation.validate_uniprot_dataset_name(id_dict, map_to, map_from)
    ids = parsing.data_to_list(ids)
    n_queries = len(ids)
    if id_dict[map_to] == id_dict[map_from]:
        return GeneIDTranslator()
    print(f"Mapping {n_queries} entries from '{map_from}' to '{map_to}'...")
    if id_dict[map_to] != id_dict['UniProtKB AC'] and id_dict[map_from] != id_dict['UniProtKB AC']:
        to_uniprot = map_gene_ids(ids, map_from, 'UniProtKB AC').mapping_dict
        from_uniprot = map_gene_ids(to_uniprot.values(), 'UniProtKB AC', map_to).mapping_dict
        output = {key: from_uniprot[val] for key, val in zip(to_uniprot.keys(), to_uniprot.values()) if
                  val in from_uniprot}
    else:
        output = {}
        duplicates = []
        if id_dict[map_to] != 'Null' and id_dict[map_from] != 'Null':
            for chunk in _format_ids_iter(ids):
                params = {
                    'from': id_dict[map_from],
                    'to': id_dict[map_to],
                    'format': 'tab',
                    'query': chunk,
                    'columns': 'id,annotation_score' if id_dict[map_to] == id_dict['UniProtKB AC'] else 'id'}
                req = requests.get(url, params=params)
                if not req.ok:
                    req.raise_for_status()
                if id_dict[map_to] == id_dict['UniProtKB AC']:
                    output.update(parsing.uniprot_tab_with_score_to_dict(req.text))
                else:
                    this_output, this_duplicates = parsing.uniprot_tab_to_dict(req.text)
                    output.update(this_output)
                    duplicates.extend(this_duplicates)
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
                output.update(parsing.uniprot_tab_with_score_to_dict(req.text, reverse_key_value=True))

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
def fetch_go_basic() -> parsing.DAGTree:
    """
    Fetches the basic Gene Ontology OBO file from the geneontology.org website ('go-basic.obo') and parses it into a \
    DAGTree data structure.
    :return: a parsed DAGTree for gene ontology propagation and visualization.
    :rtype: parsing.DAGTree
    """
    url = 'http://current.geneontology.org/ontology/go-basic.obo'
    with requests.get(url, stream=True) as obo_stream:
        return parsing.DAGTree(obo_stream.iter_lines())
