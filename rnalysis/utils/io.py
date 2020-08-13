import os
import warnings
from pathlib import Path
from typing import List, Set, Union, Iterable, Tuple
import requests
import pandas as pd
import json
from itertools import chain

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
                              iter_size: int = 200000):
    url = 'http://golr-aux.geneontology.io/solr/select?'
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
    legal_evidence = set.union(*[set(s) for s in evidence_type_dict.values()])

    legal_qualifiers = {'not', 'contributes_to', 'colocalizes_with'}

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
            assert item in legals, f"Illegal item {item}. Legal items are {legals}.."
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
    req = requests.get(url, params=params)
    if not req.ok:
        req.raise_for_status()
    n_annotations = json.loads(req.text)['response']['numFound']

    print(f"Fetching {n_annotations} annotations...")

    # fetch all annotations in batches of size iter_size, and yield them one-by-one
    start = 0
    max_iters = n_annotations // iter_size + 1
    params['omitHeader'] = "true"  # omit the header from the json response
    for i in range(max_iters):
        params['start'] = start
        start += iter_size
        params['rows'] = iter_size if i < max_iters - 1 else n_annotations % iter_size
        req = requests.get(url, params=params)
        if not req.ok:
            req.raise_for_status()
        for record in json.loads(req.text)['response']['docs']:
            yield record


def map_taxon_id(taxon_name: Union[str, int]) -> Tuple[int, str]:
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
    header = res[0].split('\t')

    matched_taxon = res[1].split('\t')
    taxon_id = int(matched_taxon[header.index('Taxon')])
    scientific_name = matched_taxon[header.index('Scientific name')]
    if len(res) > 2 and not (taxon_name == taxon_id or taxon_name == scientific_name):
        warnings.warn(
            f"Found {len(res) - 1} taxons matching the search term '{taxon_name}'. "
            f"Picking the match with the highest score.")

    return taxon_id, scientific_name


class GeneIDTranslator:
    __slots__ = {'mapping_dict': 'dictionary mapping gene IDs from one type to another'}

    def __init__(self, mapping_dict: Union[dict, None] = None):
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


def map_gene_ids(ids: Union[str, Set[str], List[str]], map_from: str, map_to: str = 'UniProtKB AC') -> GeneIDTranslator:
    url = 'https://www.uniprot.org/uploadlists/'
    id_dict = _load_id_abbreviation_dict()
    validation.validate_uniprot_dataset_name(id_dict, map_to, map_from)

    if id_dict[map_to] == id_dict[map_from]:
        return GeneIDTranslator()

    ids = parsing.data_to_list(ids)
    n_queries = len(ids)
    print(f"Mapping {n_queries} entries from '{map_from}' to '{map_to}'...")
    output = {}
    if id_dict[map_to] != 'Null' and id_dict[map_from] != 'Null':
        for chunk in _format_ids_iter(ids):
            params = {
                'from': id_dict[map_from],
                'to': id_dict[map_to],
                'format': 'tab',
                'query': chunk, }
            req = requests.get(url, params=params)
            if not req.ok:
                req.raise_for_status()
            output.update(parsing.uniprot_tab_to_dict(req.text))

    if len(output) < n_queries:
        warnings.warn(f"Failed to map {n_queries - len(output)} entries from '{map_from}' to '{map_to}'. "
                      f"Returning {len(output)} successfully-mapped entries.")
    return GeneIDTranslator(output)


def _format_ids(ids: Union[str, int, list, set]):
    if isinstance(ids, str):
        return ids
    elif isinstance(ids, int):
        return str(ids)
    return " ".join((str(item) for item in ids))


def _format_ids_iter(ids: Union[str, int, list, set], chunk_size: int = 500):
    if isinstance(ids, str):
        return ids
    elif isinstance(ids, int):
        return str(ids)
    for i in range(0, len(ids), chunk_size):
        j = min(chunk_size, len(ids) - i)
        yield " ".join((str(item) for item in ids[i:j]))


def _load_id_abbreviation_dict(dict_path: str = os.path.join(__path__[0], 'uniprot_dataset_abbreviation_dict.json')):
    with open(dict_path) as f:
        return json.load(f)


def fetch_go_basic():
    url = 'http://current.geneontology.org/ontology/go-basic.obo'
    with requests.get(url, stream=True) as obo_stream:
        return parsing.DAGTreeParser(obo_stream.iter_lines())
