import asyncio
import concurrent.futures
import functools
import gzip
import hashlib
import inspect
import json
import os
import queue
import re
import shlex
import shutil
import subprocess
import time
import typing
import warnings
from datetime import date, datetime
from functools import lru_cache
from io import StringIO
from itertools import chain
from pathlib import Path
from sys import executable
from typing import List, Set, Union, Iterable, Tuple, Dict, Any, Callable
from urllib.parse import urlparse, parse_qs, urlencode

import aiohttp
import aiolimiter
import appdirs
import numpy as np
import pandas as pd
import requests
import yaml
from defusedxml import ElementTree
from requests.adapters import HTTPAdapter, Retry
from tqdm.auto import tqdm

try:
    from typing import Literal
except ImportError:
    from typing_extensions import Literal
from rnalysis.utils import parsing, validation, __path__
from rnalysis import __version__


def get_gui_cache_dir() -> Path:
    cache_dir = Path(appdirs.user_cache_dir('RNAlysis'))
    return cache_dir.joinpath('rnalysis_gui')


def get_data_dir() -> Path:
    data_dir = Path(appdirs.user_data_dir('RNAlysis', roaming=True))
    return data_dir


def get_tutorial_videos_dir() -> Path:
    data_dir = get_data_dir()
    return data_dir.joinpath('videos')


def get_todays_cache_dir() -> Path:
    today = date.today().strftime('%Y_%m_%d')
    cache_dir = Path(appdirs.user_cache_dir('RNAlysis'))
    return cache_dir.joinpath(today)


def load_cached_file(filename: str):
    directory = get_todays_cache_dir()
    file_path = directory.joinpath(filename)
    if file_path.exists():
        with open(file_path) as f:
            return f.read()
    else:
        return None


def cache_file(content: str, filename: str):
    directory = get_todays_cache_dir()
    if not directory.exists():
        directory.mkdir(parents=True)
    file_path = directory.joinpath(filename)
    with open(file_path, 'w') as f:
        f.write(content)


def clear_directory(directory: Union[str, Path]):
    directory = Path(directory)
    if not directory.exists():
        return

    for item in directory.iterdir():
        if item.is_file():
            item.unlink()
        elif item.is_dir():
            shutil.rmtree(item, ignore_errors=True)


def clear_cache():
    cache_dir = Path(appdirs.user_cache_dir('RNAlysis'))
    clear_directory(cache_dir)


def clear_gui_cache():
    directory = get_gui_cache_dir()
    clear_directory(directory)


def load_cached_gui_file(filename: str):
    directory = get_gui_cache_dir()
    file_path = directory.joinpath(filename)
    if file_path.exists():
        if file_path.suffix in {'.csv', '.tsv'}:
            return load_csv(file_path, index_col=0)
        elif file_path.suffix in {'.txt'}:
            with open(file_path) as f:
                return {item.strip() for item in f.readlines()}

        else:
            with open(file_path) as f:
                return f.read()
    else:
        return None


def cache_gui_file(item: Union[pd.DataFrame, set, str], filename: str):
    directory = get_gui_cache_dir()
    if not directory.exists():
        directory.mkdir(parents=True)
    file_path = directory.joinpath(filename)
    if isinstance(item, (pd.DataFrame, pd.Series)):
        save_csv(item, file_path, index=True)
    elif isinstance(item, set):
        save_gene_set(item, file_path)
    elif isinstance(item, str):
        with open(file_path, 'w') as f:
            f.write(item)
    else:
        raise TypeError(type(item))


def save_gui_session(session_filename: Union[str, Path], file_names: List[str], item_names: List[str], item_types: list,
                     item_properties: list, pipeline_names: List[str], pipeline_files: List[str]):
    session_filename = Path(session_filename)
    session_folder = session_filename
    if session_folder.exists():
        if session_folder.is_dir():
            shutil.rmtree(session_folder)
        else:
            session_folder.unlink()
    session_folder.mkdir(parents=True)

    session_data = dict(files=dict(), pipelines=dict(), metadata=dict())
    for file_name, item_name, item_type, item_property in zip(file_names, item_names, item_types, item_properties):
        shutil.move(Path(get_gui_cache_dir().joinpath(file_name)), session_folder.joinpath(file_name))
        session_data['files'][file_name] = (item_name, item_type.__name__, item_property)

    for i, (pipeline_name, pipeline_file) in enumerate(zip(pipeline_names, pipeline_files)):
        pipeline_filename = session_folder.joinpath(f"pipeline_{i}.yaml")
        Path(pipeline_filename).write_text(pipeline_file)
        session_data['pipelines'][pipeline_filename.name] = pipeline_name

    session_data['metadata']['creation_time'] = get_datetime()
    session_data['metadata']['name'] = Path(session_filename).stem
    session_data['metadata']['n_tabs'] = len(session_data['files'])
    session_data['metadata']['n_pipelines'] = len(session_data['pipelines'])
    session_data['metadata']['tab_order'] = file_names

    with open(session_folder.joinpath('session_data.yaml'), 'w') as f:
        yaml.safe_dump(session_data, f)
    shutil.make_archive(session_folder.with_suffix(''), 'zip', session_folder)
    shutil.rmtree(session_folder)
    session_filename.with_suffix('.zip').replace(session_filename.with_suffix('.rnal'))


def load_gui_session(session_filename: Union[str, Path]):
    session_filename = Path(session_filename)
    try:
        session_filename.with_suffix('.rnal').rename(session_filename.with_suffix('.rnal.zip'))
        shutil.unpack_archive(session_filename.with_suffix('.rnal.zip'),
                              get_gui_cache_dir().joinpath(session_filename.name))
    finally:
        session_filename.with_suffix('.rnal.zip').rename(session_filename.with_suffix('.rnal'))

    session_dir = get_gui_cache_dir().joinpath(session_filename.name)
    assert session_dir.exists()

    items = []
    item_names = []
    item_types = []
    item_properties = []
    pipeline_files = []
    pipeline_names = []
    with open(session_dir.joinpath('session_data.yaml')) as f:
        session_data = yaml.safe_load(f)
    if 'tab_order' in session_data['metadata'] and len(session_data['metadata']['tab_order']) == len(
        session_data['files']):
        filenames = session_data['metadata']['tab_order']
    else:
        filenames = session_data['files'].keys()

    for file_name in filenames:
        file_path = session_dir.joinpath(file_name)
        assert file_path.exists() and file_path.is_file()
        item = load_cached_gui_file(Path(session_filename.name).joinpath(file_name))
        item_name, item_type, item_property = session_data['files'][file_name]
        items.append(item)
        item_names.append(item_name)
        item_types.append(item_type)
        item_properties.append(item_property)

    for pipeline_filename in session_data['pipelines'].keys():
        pipeline_path = session_dir.joinpath(pipeline_filename)
        assert pipeline_path.exists() and pipeline_path.is_file()
        pipeline_files.append(pipeline_path.read_text())
        pipeline_names.append(session_data['pipelines'][pipeline_filename])

    shutil.rmtree(session_dir)
    return items, item_names, item_types, item_properties, pipeline_names, pipeline_files


def load_csv(filename: Union[str, Path], index_col: int = None, drop_columns: Union[str, List[str]] = False,
             squeeze=False, comment: str = None):
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
    kwargs = dict(sep=None, engine='python', encoding=encoding, comment=comment, skipinitialspace=True)
    if index_col is not None:
        kwargs['index_col'] = index_col
    df = pd.read_csv(filename, **kwargs)
    if squeeze:
        df = df.squeeze("columns")
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


def save_csv(df: pd.DataFrame, filename: Union[str, Path], suffix: str = None, index: bool = True):
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


class KEGGAnnotationIterator:
    URL = 'https://rest.kegg.jp/'
    TAXON_MAPPING_URL = 'https://www.genome.jp/kegg-bin/download_htext?htext=br08610'
    REQUEST_DELAY_MILLIS = 250
    REQ_MAX_ENTRIES = 10
    TAXON_TREE_CACHED_FILENAME = 'kegg_taxon_tree.json'
    PATHWAY_NAMES_CACHED_FILENAME = 'kegg_pathway_list.txt'
    COMPOUND_LIST_CACHED_FILENAME = 'kegg_compound_list.txt'
    GLYCAN_LIST_CACHED_FILENAME = 'kegg_glycan_list.txt'

    def __init__(self, taxon_id: int, pathways: Union[str, List[str], Literal['all']] = 'all'):
        self.pathway_names = {}
        self.taxon_id = taxon_id
        self.organism_code = self.get_kegg_organism_code(taxon_id)
        if pathways == 'all':
            self.pathway_names, self.n_annotations = self.get_pathways()
        else:
            pathways = parsing.data_to_list(pathways)
            assert len(pathways) > 0, f"No KEGG pathway IDs were given!"
            self.pathway_names = {pathway: None for pathway in pathways}
            self.n_annotations = len(self.pathway_names)
        self.pathway_annotations = None

    @staticmethod
    def _kegg_request(operation: str, arguments: Union[str, List[str]], cached_filename: Union[str, None] = None
                      ) -> Tuple[str, bool]:
        if cached_filename is not None:
            cached_file = load_cached_file(cached_filename)
            if cached_file is not None:
                is_cached = True
                return cached_file, is_cached

        is_cached = False
        address = KEGGAnnotationIterator.URL + operation + '/' + '/'.join(parsing.data_to_list(arguments))
        response = requests.get(address)
        if not response.ok:
            response.raise_for_status()
        if cached_filename is not None:
            cache_file(response.text, cached_filename)
        return response.text, is_cached

    def _generate_cached_filename(self, pathways: Tuple[str, ...]) -> str:
        fname = f'{self.taxon_id}' + ''.join(pathways).replace('path:', '') + '.json'
        return fname

    @staticmethod
    def _get_taxon_tree():
        cached_filename = KEGGAnnotationIterator.TAXON_TREE_CACHED_FILENAME
        cached_file = load_cached_file(cached_filename)
        if cached_file is not None:
            try:
                taxon_tree = json.loads(cached_file)
                return taxon_tree
            except json.decoder.JSONDecodeError:
                pass
        with requests.get(KEGGAnnotationIterator.TAXON_MAPPING_URL, params=dict(format='json')) as req:
            content = req.content.decode('utf8')
            cache_file(content, cached_filename)
            taxon_tree = json.loads(content)
        return taxon_tree

    @staticmethod
    def get_compounds() -> Dict[str, str]:
        compounds = {}
        data = KEGGAnnotationIterator._kegg_request('list', ['compound'],
                                                    KEGGAnnotationIterator.COMPOUND_LIST_CACHED_FILENAME)[0] + '\n' + \
               KEGGAnnotationIterator._kegg_request('list', ['glycan'],
                                                    KEGGAnnotationIterator.GLYCAN_LIST_CACHED_FILENAME)[0]
        data = data.split('\n')
        for line in data:
            split = line.split('\t')
            if len(split) == 2:
                pathway_code, compound_names = split
                main_name = compound_names.split(';')[0]
                compounds[pathway_code] = main_name

        return compounds

    @staticmethod
    @functools.lru_cache(1024)
    def get_kegg_organism_code(taxon_id: int) -> str:
        taxon_tree = KEGGAnnotationIterator._get_taxon_tree()
        q = queue.Queue()
        q.put(taxon_tree)
        while not q.empty():
            this_item = q.get()
            if f"[TAX:{taxon_id}]" in this_item['name']:
                child = this_item['children'][0]
                organism_code = child['name'].split(" ")[0]
                return organism_code
            else:
                children = this_item.get('children', tuple())
                for child in children:
                    q.put(child)
        raise ValueError(f"Could not find organism code for taxon ID {taxon_id}. ")

    def get_pathways(self) -> Tuple[Dict[str, str], int]:
        pathway_names = {}
        data, _ = self._kegg_request('list', ['pathway', self.organism_code], self.PATHWAY_NAMES_CACHED_FILENAME)
        data = data.split('\n')
        for line in data:
            split = line.split('\t')
            if len(split) == 2:
                pathway_code, pathway_name = split
                pathway_names[pathway_code] = pathway_name

        n_annotations = len(pathway_names)
        return pathway_names, n_annotations

    @staticmethod
    def get_pathway_kgml(pathway_id: str) -> ElementTree:
        cached_filename = f'kgml_{pathway_id}.xml'
        data, _ = KEGGAnnotationIterator._kegg_request('get', [pathway_id, 'kgml'], cached_filename)
        cache_file(data, cached_filename)
        return ElementTree.parse(StringIO(data))

    def get_pathway_annotations(self):
        if self.pathway_annotations is not None:
            for pathway, annotations in self.pathway_annotations.items():
                name = self.pathway_names[pathway]
                yield pathway, name, annotations
        else:
            pathway_annotations = {}
            partitioned_pathways = parsing.partition_list(list(self.pathway_names.keys()), self.REQ_MAX_ENTRIES)
            for chunk in partitioned_pathways:
                prev_time = time.time()
                data, was_cached = self._kegg_request('get', '+'.join(chunk), self._generate_cached_filename(chunk))
                entries = data.split('ENTRY')[1:]
                for entry in entries:
                    entry_split = entry.split('\n')
                    pathway = entry_split[0].split()[0]
                    if pathway not in chunk:
                        print(f'Could not find pathway {pathway} in requested chunk: {chunk}')
                    pathway_name = self.pathway_names[pathway]
                    pathway_annotations[pathway] = set()
                    genes_startline = 0
                    genes_endline = 0
                    for i, line in enumerate(entry_split):
                        if line.startswith('GENE'):
                            genes_startline = i
                        elif line.startswith('COMPOUND'):
                            genes_endline = i
                            break
                    for line_num in range(genes_startline, genes_endline):
                        line = entry_split[line_num]
                        if line.startswith('GENE'):
                            line = line.replace('GENE', '', 1)
                        gene_id = f"{self.organism_code}:{line.strip().split(' ')[0]}"
                        pathway_annotations[pathway].add(gene_id)
                        yield pathway, pathway_name, pathway_annotations[pathway]
                if not was_cached:
                    delay = max((self.REQUEST_DELAY_MILLIS / 1000) - (time.time() - prev_time), 0)
                    time.sleep(delay)

            self.pathway_annotations = pathway_annotations

    def __iter__(self):
        return self.get_pathway_annotations()


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
        return json.loads(self._golr_request(self.default_params, self._generate_cached_filename(None)))['response'][
            'numFound']

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
    def _golr_request(params: dict, cached_filename: Union[str, None] = None) -> str:
        """
        Run a get request to the GOlr server with the specified parameters, and return the server's text response.
        :param params: the get request's parameters.
        :type params: dict
        """
        if cached_filename is not None:
            cached_file = load_cached_file(cached_filename)
            if cached_file is not None:
                return cached_file

        response = requests.get(GOlrAnnotationIterator.URL, params=params)
        if not response.ok:
            response.raise_for_status()
        if cached_filename is not None:
            cache_file(response.text, cached_filename)
        return response.text

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

    def _generate_cached_filename(self, start: Union[int, None]) -> str:
        fname = f'{self.taxon_id}' + ''.join(
            [f'{aspect}' for aspect in parsing.data_to_tuple(self.aspects, True)]) + ''.join(
            [f'{evidence_type}' for evidence_type in parsing.data_to_tuple(self.evidence_types, True)])
        # add union of all requested databases to query
        if not self.databases == {'any'}:
            fname += ''.join(f'{db}' for db in parsing.data_to_tuple(self.databases, True))
        # add union of all requested qualifiers to query
        if len(self.qualifiers) > 0:
            fname += ''.join([f'{qual}' for qual in parsing.data_to_tuple(self.qualifiers, True)])

        # exclude all 'excluded' items from query
        fname += 'exc' + ''.join(
            [f'{evidence_type}' for evidence_type in parsing.data_to_tuple(self.excluded_evidence_types, True)])
        fname += ''.join([f'{db}' for db in parsing.data_to_tuple(self.excluded_databases, True)])
        fname += ''.join([f'{qual}' for qual in parsing.data_to_tuple(self.excluded_qualifiers, True)])

        return fname + str(start) + '.json'

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
                processes.append(executor.submit(self._golr_request, param_dict,
                                                 self._generate_cached_filename(param_dict['start'])))
        for task in concurrent.futures.as_completed(processes):
            for record in json.loads(task.result())['response']['docs']:
                yield record

    def __iter__(self):
        return self._annotation_generator_func()


def get_obo_basic_stream():
    url = 'http://current.geneontology.org/ontology/go-basic.obo'
    with requests.get(url, stream=True) as obo_stream:
        content = obo_stream.content.decode('utf8')
        return content


# TODO: cache this! save and load gene IDs individually
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
    processes = []
    with concurrent.futures.ThreadPoolExecutor() as executor:
        for chunk in data_chunks:
            data = {"ids": parsing.data_to_list(chunk)}
            processes.append(
                executor.submit(requests.post, url, headers=headers, data=data.__repr__().replace("'", '"')))

    # req = requests.post(url, headers=headers, data=data.__repr__().replace("'", '"'))
    for task in concurrent.futures.as_completed(processes):
        req = task.result()
        if not req.ok:
            req.raise_for_status()
        output.update(req.json())
    return output


def infer_sources_from_gene_ids(gene_ids: Iterable[str]) -> Dict[str, Set[str]]:
    """
    #TODO
    Infer the
    :param gene_ids:
    :type gene_ids:
    :return:
    :rtype:
    """
    translator, map_from, _ = find_best_gene_mapping(parsing.data_to_tuple(gene_ids), map_from_options=None,
                                                     map_to_options=('Ensembl', 'Ensembl Genomes'))
    output = _ensmbl_lookup_post_request(parsing.data_to_tuple(translator.mapping_dict.values()))
    sources = {}
    for gene_id in output:
        if output[gene_id] is not None:
            source = output[gene_id]['source']
            if source not in sources:
                sources[source] = set()
            sources[source].add(gene_id)
    return sources


def infer_taxon_from_gene_ids(gene_ids: Iterable[str], gene_id_type: str = None) -> Tuple[Tuple[int, str], typing.Any]:
    """
    Infer the NCBI Taxon ID and name of a collection of gene IDs. \
    In cases where not all gene IDs map to the same taxon, the best-fitting taxon will be picked by a majority vote.

    :param gene_ids: a collection of gene IDs
    :type gene_ids: Iterable of str
    :return: a tuple of the best-matching taxon's NCBI Taxon ID and full scientific name.
    :rtype: Tuple[int ,str]
    """
    if gene_id_type is not None:
        gene_id_type = parsing.data_to_tuple(gene_id_type)
    translator, map_from, _ = find_best_gene_mapping(parsing.data_to_tuple(gene_ids), map_from_options=gene_id_type,
                                                     map_to_options=('Ensembl', 'Ensembl Genomes'))
    output = _ensmbl_lookup_post_request(parsing.data_to_tuple(translator.mapping_dict.values()))
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
    return map_taxon_id(chosen_species.replace('_', ' ')), map_from


@lru_cache(maxsize=32, typed=False)
def map_taxon_id(taxon_name: Union[str, int]) -> Tuple[int, str]:
    """
    Maps a given query (taxon name or NCBI Taxon ID) to the best-matching taxon from the NCBI taxonomy database. \
    Mapping is done through UniProt Taxonomy: https://rest.uniprot.org/taxonomy/search?

    :param taxon_name: a partial/full taxon name (str) or NCBI Taxon ID (int) to map
    :type taxon_name: int or str
    :return: a tuple of the best-matching taxon's NCBI Taxon ID and full scientific name.
    :rtype: Tuple[int ,str]
    """
    url = 'https://rest.uniprot.org/taxonomy/search?'

    params = {
        'format': 'tsv',
        'query': taxon_name,
    }
    req = requests.get(url, params=params)
    if not req.ok:
        req.raise_for_status()
    res = pd.read_csv(StringIO(req.text), sep='\t').sort_values(by='Taxon Id', ascending=True)
    if res.shape[0] == 0:
        raise ValueError(f"No taxons match the search query '{taxon_name}'.")

    taxon_id = int(res['Taxon Id'].iloc[0])
    scientific_name = res['Scientific name'].iloc[0]

    if res.shape[0] > 2 and not (taxon_name == taxon_id or taxon_name == scientific_name):
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

    def __len__(self):
        if self.mapping_dict is None:
            return 0
        return len(self.mapping_dict)

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


def submit_id_mapping(url: str, from_db: str, to_db: str, ids: List[str]):
    req = requests.post(f"{url}/idmapping/run", data={"from": from_db, "to": to_db, "ids": ",".join(ids)})
    req.raise_for_status()
    return req.json()["jobId"]


def get_next_link(headers):
    re_next_link = re.compile(r'<(.+)>; rel="next"')
    if "Link" in headers:
        match = re_next_link.match(headers["Link"])
        if match:
            return match.group(1)


def check_id_mapping_results_ready(session, url: str, job_id, polling_interval: float, verbose: bool = True):
    while True:
        r = session.get(f"{url}/idmapping/status/{job_id}")
        r.raise_for_status()
        j = r.json()
        if "jobStatus" in j:
            if j["jobStatus"] in {"RUNNING", "NEW"}:
                if verbose:
                    print(f"Retrying in {polling_interval}s")
                time.sleep(polling_interval)
            else:
                raise Exception(j["jobStatus"])
        else:
            return bool(j.get("results", False) or j.get("failedIds", False))


def get_batch(session, batch_response, file_format):
    batch_url = get_next_link(batch_response.headers)
    while batch_url:
        batch_response = session.get(batch_url)
        batch_response.raise_for_status()
        yield decode_results(batch_response, file_format)
        batch_url = get_next_link(batch_response.headers)


def combine_batches(all_results, batch_results, file_format):
    if file_format == "json":
        for key in ("results", "failedIds"):
            if batch_results[key]:
                all_results[key] += batch_results[key]
    elif file_format == "tsv":
        return all_results + batch_results[1:]  # dump the table header line
    else:
        return all_results + batch_results
    return all_results


def get_id_mapping_results_link(session, url, job_id):
    url = f"{url}/idmapping/details/{job_id}"
    r = session.get(url)
    r.raise_for_status()
    return r.json()["redirectURL"]


def decode_results(response, file_format):
    if file_format == "json":
        return response.json()
    elif file_format == "tsv":
        return [line for line in response.text.split("\n") if line]
    return response.text


def get_id_mapping_results_search(session, url, verbose: bool = True):
    parsed = urlparse(url)
    query = parse_qs(parsed.query)
    query["fields"] = 'accession,annotation_score'
    query["format"] = "tsv"
    file_format = 'tsv'
    # file_format = query["format"][0] if "format" in query else "json"
    if "size" in query:
        size = int(query["size"][0])
    else:
        size = 500
        query["size"] = size

    parsed = parsed._replace(query=urlencode(query, doseq=True))
    url = parsed.geturl()
    r = session.get(url)
    r.raise_for_status()
    results = decode_results(r, file_format)
    try:
        total = int(r.headers["x-total-results"])
    except KeyError:
        return ''
    if verbose:
        print_progress_batches(0, size, total)
    for i, batch in enumerate(get_batch(session, r, file_format)):
        results = combine_batches(results, batch, file_format)
        if verbose:
            print_progress_batches(i + 1, size, total)
    return results


def print_progress_batches(batch_index, size, total):
    n = min((batch_index + 1) * size, total)
    print(f"Fetched: {n} / {total}")


def get_mapping_results(api_url: str, from_db: str, to_db: str, ids: List[str], polling_interval: float, session,
                        verbose: bool = True):
    job_id = submit_id_mapping(api_url, from_db=from_db, to_db=to_db, ids=ids)
    if check_id_mapping_results_ready(session, api_url, job_id, polling_interval, verbose=verbose):
        link = get_id_mapping_results_link(session, api_url, job_id)
        results = get_id_mapping_results_search(session, link, verbose)
        return results


def map_gene_ids(ids: Union[str, Iterable[str]], map_from: str, map_to: str = 'UniProtKB AC',
                 verbose: bool = True) -> GeneIDTranslator:
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
    :param verbose: if False, verbose reports of mapping success/failure will be suppressed.
    :type verbose: bool (default=True)
    :return:a GeneIDTranslator object that uniquely maps each given gene ID in 'map_from' identifier type \
    to its matching gene ID in 'map_to' identifier type.
    :rtype: GeneIDTranslator
    """
    UNIPROTKB_FROM = "UniProtKB_from"
    UNIPROTKB_TO = "UniProtKB_to"

    if len(ids) == 0:
        if verbose:
            warnings.warn('No IDs were given')
        return GeneIDTranslator({})

    # if map_from and map_to are the same, return an empty GeneIDTranslator (which will map any given gene ID to itself)
    id_dicts = _get_id_abbreviation_dicts()
    validation.validate_uniprot_dataset_name(id_dicts, parsing.data_to_list(map_to), parsing.data_to_list(map_from))
    id_dict_to, id_dict_from = id_dicts
    if id_dict_to[map_to] == id_dict_from[map_from]:
        return GeneIDTranslator()

    if map_to in ["UniProtKB", "UniProtKB AC/ID"]:
        map_to = UNIPROTKB_TO
    if map_from in ["UniProtKB", "UniProtKB AC/ID"]:
        map_from = UNIPROTKB_FROM

    ids = parsing.data_to_list(ids)
    n_queries = len(ids)

    # since the Uniprot service can only translate to or from 'UniProtKB' identifier type,
    # if we need to map gene IDs between two other identifier types,
    # then we will map from 'map_from' to 'UniProtKB' and then from 'UniProtKB' to 'map_to'.
    if verbose:
        print(f"Mapping {n_queries} entries from '{map_from}' to '{map_to}'...")

    if id_dict_to[map_to] != id_dict_to[UNIPROTKB_TO] and id_dict_from[map_from] != id_dict_from[UNIPROTKB_FROM]:
        to_uniprot = map_gene_ids(ids, map_from, UNIPROTKB_TO, verbose=verbose).mapping_dict
        if to_uniprot is None:
            to_uniprot = {}
        from_uniprot = map_gene_ids(to_uniprot.values(), UNIPROTKB_FROM, map_to, verbose=verbose).mapping_dict
        if from_uniprot is None:
            from_uniprot = {}
        output_dict = {key: from_uniprot[val] for key, val in zip(to_uniprot.keys(), to_uniprot.values()) if
                       val in from_uniprot}

    # make sure that 'map_from' and 'map_to' are recognized identifier types
    elif id_dict_to[map_to] != 'Null' and id_dict_from[map_from] != 'Null':
        POLLING_INTERVAL = 3
        API_URL = "https://rest.uniprot.org"

        retries = Retry(total=5, backoff_factor=0.25, status_forcelist=[500, 502, 503, 504])
        session = requests.Session()
        session.mount("https://", HTTPAdapter(max_retries=retries))

        results = get_mapping_results(api_url=API_URL, from_db=id_dict_from[map_from], to_db=id_dict_to[map_to],
                                      ids=ids,
                                      polling_interval=POLLING_INTERVAL, session=session, verbose=verbose)

        if results is None or len(results) <= 1:
            if verbose:
                warnings.warn(f"No entries were mapped successfully.")
            return GeneIDTranslator({})

        df = pd.DataFrame([line.split('\t') for line in results[1:]], columns=results[0].split('\t'))
        # sort annotations by decreasing annotation score, so that the most relevant annotations are at the top
        if 'Annotation' in df.columns:
            df['Annotation'][df['Annotation'] == ''] = '0'
            df['Annotation'] = (df['Annotation']).astype(float)
            df = df.sort_values('Annotation', ascending=False)
        output_dict = {}
        duplicates = {}

        # sort duplicates from one-to-one mappings
        for match in df.iterrows():
            match_from = match[1][0]
            match_to = match[1][1]
            if match_from in output_dict or match_from in duplicates:
                if match_from not in duplicates:
                    duplicates[match_from] = [output_dict.pop(match_from)]
                duplicates[match_from].append(match_to)
            else:
                output_dict[match_from] = match_to

        # handle duplicates
        if len(duplicates) > 0:
            if map_to == UNIPROTKB_TO:
                for match_from, match_to_options in duplicates.items():
                    output_dict[match_from] = match_to_options[0]
                duplicates_chosen = {match_from: match_to[0] for match_from, match_to in duplicates.items()}

            # if there are unproccessed duplicates, map them in reverse and sort then by annotation score
            else:
                ids_to_rev_map = parsing.flatten(parsing.data_to_list(duplicates.values()))

                rev_results = get_mapping_results(api_url=API_URL, from_db=id_dict_to[map_to],
                                                  to_db=id_dict_from[UNIPROTKB_TO], ids=ids_to_rev_map,
                                                  polling_interval=POLLING_INTERVAL, session=session, verbose=verbose)
                # TODO: if job fails?
                rev_df = pd.DataFrame([line.split('\t') for line in rev_results[1:]],
                                      columns=rev_results[0].split('\t'))
                rev_df['Annotation'] = (rev_df['Annotation']).astype(float)
                rev_df = rev_df.sort_values('Annotation', ascending=False)
                duplicates_chosen = {}
                for match in rev_df.iterrows():
                    match_from_rev = match[1][0]
                    match_to_rev = match[1][1]
                    if match_to_rev not in output_dict:
                        output_dict[match_to_rev] = match_from_rev
                        duplicates_chosen[match_to_rev] = match_from_rev
            if verbose:
                warnings.warn(f"Duplicate mappings were found for {len(duplicates)} genes.  The following mapping "
                              f"was chosen for them based on their annotation score: {duplicates_chosen}")
    else:
        output_dict = {}
    if len(output_dict) < n_queries and verbose:
        warnings.warn(f"Failed to map {n_queries - len(output_dict)} entries from '{map_from}' to '{map_to}'. "
                      f"Returning {len(output_dict)} successfully-mapped entries.")

    if map_to == 'Ensembl':
        for key, val in output_dict.items():
            output_dict[key] = re.sub('(\.\d+)$', '', val)

    return GeneIDTranslator(output_dict)


def _format_ids_iter(ids: Union[str, int, list, set], chunk_size: int = 250):
    if isinstance(ids, str):
        yield ids
    elif isinstance(ids, int):
        yield str(ids)
    else:
        for i in range(0, len(ids), chunk_size):
            j = min(chunk_size, len(ids) - i)
            yield " ".join((str(item) for item in ids[i:i + j]))


@functools.lru_cache(maxsize=2048)
def find_best_gene_mapping(ids: Tuple[str, ...], map_from_options: Union[Tuple[str, ...], None],
                           map_to_options: Union[Tuple[str, ...], None]):
    all_map_to_options, all_map_from_options = _get_id_abbreviation_dicts()
    if map_to_options is None:
        map_to_options = all_map_to_options
    if map_from_options is None:
        map_from_options = all_map_from_options

    def _key_func(items: Tuple[int, str, str]):
        key = [items[0]]
        key.append(list(map_from_options).index(items[1]))
        key.append(list(map_to_options).index(items[2]))
        return key

    def map_gene_ids_ignore_httpexception(ids: Tuple[str], map_from: str, map_to: str):
        try:
            return map_gene_ids(ids, map_from, map_to, verbose=False), map_from, map_to
        except requests.exceptions.HTTPError:
            return GeneIDTranslator({}), map_from, map_to

    with concurrent.futures.ThreadPoolExecutor(max_workers=5) as executor:
        processes = []
        mapping_combs_from = []
        mapping_combs_to = []
        for map_from in map_from_options:
            for map_to in parsing.data_to_tuple(map_to_options):
                if map_to == map_from:
                    continue
                mapping_combs_from.append(map_from)
                mapping_combs_to.append(map_to)
                # processes.append(executor.submit(map_gene_ids_ignore_httpexception, ids, map_from, map_to))

        results = list(tqdm(executor.map(functools.partial(map_gene_ids_ignore_httpexception, ids), mapping_combs_from,
                                         mapping_combs_to),
                            total=len(mapping_combs_from), desc='Submitting jobs...', unit='jobs'))
    # best_result = None
    # best_result_len = -1
    parsed_results = {}
    with tqdm(desc='Finding the best-matching gene ID type...', total=len(processes) + 1) as pbar:
        pbar.update()
        for result in results:
            mapping_dict, map_from, map_to = result
            result_len = len(mapping_dict)
            parsed_results[(result_len, map_from, map_to)] = result
            pbar.update()

    sorted_keys = sorted(parsed_results.keys(), key=_key_func, reverse=True)
    return parsed_results[sorted_keys[0]]


@functools.lru_cache(maxsize=2)
def get_legal_gene_id_types():
    URL = 'https://rest.uniprot.org/configure/idmapping/fields'
    GROUP_PRIORITIZATION = ['UniProt',
                            'Genome annotation databases',
                            'Organism-specific databases',
                            'Phylogenomic databases',
                            'Sequence databases',
                            'Miscellaneous',
                            'Gene expression databases',
                            'Enzyme and pathway databases',
                            'Proteomic databases']
    abbrev_dict_to = {}
    abbrev_dict_from = {}

    req = requests.get(URL)
    req.raise_for_status()
    entries = json.loads(req.text)['groups']
    entries_filtered = []
    for entry in entries:
        if entry['groupName'] in GROUP_PRIORITIZATION:
            entries_filtered.append(entry)
    entries_sorted = sorted(entries_filtered, key=lambda grp: GROUP_PRIORITIZATION.index(grp['groupName']))
    for grp in entries_sorted:
        for entry in grp['items']:
            name = entry['displayName']
            abbrev = entry['name']
            is_from = entry['from']
            is_to = entry['to']
            if is_to:
                abbrev_dict_to[name] = abbrev
            if is_from:
                abbrev_dict_from[name] = abbrev
    return abbrev_dict_from, abbrev_dict_to


@functools.lru_cache(maxsize=2)
def _get_id_abbreviation_dicts(dict_path: str = os.path.join(Path(__path__[0]).parent,
                                                             'data_files/uniprot_dataset_abbreviation_dict.json')):
    with open(dict_path) as f:
        abbrev_dict_to = json.load(f)
        abbrev_dict_from = abbrev_dict_to.copy()

    legal_from, legal_to = get_legal_gene_id_types()

    abbrev_dict_from.update(legal_from)
    abbrev_dict_to.update(legal_to)

    for val in parsing.data_to_list(abbrev_dict_to.values()):
        abbrev_dict_to[val] = val
    for val in parsing.data_to_list(abbrev_dict_from.values()):
        abbrev_dict_from[val] = val

    return abbrev_dict_to, abbrev_dict_from


def get_datetime():
    now = datetime.now()
    now_str = now.strftime('%H:%M:%S %Y/%m/%d')
    return now_str


def save_gene_set(gene_set: set, path):
    with open(path, 'w') as f:
        f.writelines(
            [f"{item}\n" if (i + 1) < len(gene_set) else f"{item}" for i, item in enumerate(gene_set)])


def calculate_checksum(filename: Union[str, Path]):  # pragma: no cover
    assert Path(filename).exists(), f"file '{filename}' does not exist!"
    with open(filename, 'rb') as file_to_check:
        # read contents of the file
        data = file_to_check.read()
        # pipe contents of the file through
        md5_checksum = hashlib.md5(data).hexdigest()
        return md5_checksum


async def get_video_remote_checksum(video_name: str, session, semaphore, limiter):  # pragma: no cover
    url = 'https://github.com/GuyTeichman/RNAlysis/blob/master/rnalysis/gui/videos/checksums/' \
          f'{Path(video_name).stem}.txt'
    await semaphore.acquire()
    async with limiter:
        async with session.get(url, params=dict(raw='True')) as response:
            content = await response.text()
            semaphore.release()
            return content


async def download_video(video_path: Path, session, semaphore, limiter) -> None:  # pragma: no cover
    content = await _get_video_content(video_path, session, semaphore, limiter)
    await _write_video_to_file(video_path, content)


async def _get_video_content(video_path: Path, session, semaphore, limiter) -> bytes:  # pragma: no cover
    url = 'https://github.com/GuyTeichman/RNAlysis/blob/master/rnalysis/gui/videos/' + video_path.name
    await semaphore.acquire()
    async with limiter:
        async with session.get(url, params=dict(raw='True')) as response:
            content = await response.read()
            semaphore.release()
            return content


async def _write_video_to_file(video_path: Path, content: bytes) -> None:  # pragma: no cover
    with open(video_path, 'wb') as file:
        file.write(content)


async def get_gui_videos(video_filenames: Tuple[str, ...]):  # pragma: no cover
    video_dir_pth = get_tutorial_videos_dir()
    if not video_dir_pth.exists():
        video_dir_pth.mkdir(parents=True)

    videos_to_download = []
    local_checksums = []
    remote_checksums = []
    i = 0
    try:
        semaphore = asyncio.Semaphore(value=10)
        limiter = aiolimiter.AsyncLimiter(10, 1)
        async with aiohttp.ClientSession() as session:
            for video_name in video_filenames:
                this_video_path = video_dir_pth.joinpath(video_name)
                if this_video_path.exists():
                    remote_checksums.append(get_video_remote_checksum(video_name, session, semaphore, limiter))
                    local_checksums.append(calculate_checksum(this_video_path))
                else:
                    videos_to_download.append(this_video_path)

            for this_local, this_remote in zip(local_checksums, await asyncio.gather(*remote_checksums)):
                if this_local == this_remote:
                    yield i
                    i += 1
                else:
                    this_video_path = video_dir_pth.joinpath(video_name)
                    videos_to_download.append(this_video_path)

            tasks = []
            for this_video_path in videos_to_download:
                tasks.append(download_video(this_video_path, session, semaphore, limiter))
            for task in tasks:
                yield i
                i += 1
                await asyncio.gather(task)

    except aiohttp.ClientConnectorError:
        pass


def run_r_script(script_path: Union[str, Path], r_installation_folder: Union[str, Path, Literal['auto']] = 'auto'):
    if r_installation_folder == 'auto':
        prefix = "Rscript"
    else:
        prefix = f'{Path(r_installation_folder).as_posix()}/bin/Rscript'
    script_path = Path(script_path).as_posix()
    assert Path(script_path).exists() and Path(
        script_path).is_file(), f"Could not find the requested R script: {script_path}"

    return_code = run_subprocess([prefix, "--help"], False, False)
    if return_code:
        raise FileNotFoundError(f"Failed to find R executable (return code {return_code}). "
                                "Please make sure your R installation folder is correct. ")

    return_code = run_subprocess([prefix, script_path])

    assert not return_code, f"R script failed to execute (return code {return_code}). "


def run_subprocess(args: List[str], print_stdout: bool = True, print_stderr: bool = True,
                   log_filename: Union[str, None] = None, shell: bool = False):
    # join List of args into a string of args when running in shell mode
    if shell:
        try:
            args = shlex.join(args)
        except AttributeError:
            args = ' '.join([shlex.quote(arg) for arg in args])

    if print_stdout or print_stderr:
        if print_stdout and print_stderr:
            stdout = subprocess.PIPE
            stderr = subprocess.STDOUT
        elif print_stdout:
            stdout = subprocess.PIPE
            stderr = subprocess.DEVNULL
        else:
            stdout = subprocess.DEVNULL
            stderr = subprocess.PIPE

        with subprocess.Popen(args, stdout=stdout, stderr=stderr, shell=shell) as process:
            stream = process.stdout if print_stdout else process.stderr
            if log_filename is not None:
                with open(log_filename, 'w') as logfile:
                    for line in stream:
                        text = line.decode('utf8', errors="ignore")
                        print(text)
                        logfile.write(text)
            else:
                for line in stream:
                    text = line.decode('utf8', errors="ignore")
                    print(text)
        return_code = process.returncode
    else:
        return_code = subprocess.run(args, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, shell=shell).returncode

    return return_code


def is_rnalysis_outdated():
    installed_version = parsing.parse_version(__version__)
    pypi_link = 'https://pypi.python.org/pypi/RNAlysis/json'
    try:
        req = requests.get(pypi_link)
    except (ConnectionError, requests.exceptions.ConnectionError, requests.exceptions.HTTPError):
        return False
    if req.status_code != 200:
        return False
    newest_version = parsing.parse_version(json.loads(req.text)['info']['version'])
    if installed_version < newest_version:
        return True
    return False


def update_rnalysis():
    run_subprocess([executable, '-m', 'pip', 'install', '--upgrade', 'RNAlysis[all]'])


def get_method_docstring(method: Union[str, Callable], obj: object = None) -> Tuple[str, dict]:
    try:
        if isinstance(method, str):
            func = getattr(obj, method)
        else:
            func = method
        raw_docstring = inspect.cleandoc(inspect.getdoc(func))
        return parsing.parse_docstring(raw_docstring)
    except AttributeError:
        return '', {}


def generate_base_call(command: str, installation_folder: Union[str, Path, Literal['auto']],
                       version_command: str = '--version', args=tuple(), shell: bool = False):
    if installation_folder == 'auto':
        call = [command] + parsing.data_to_list(args)
    else:
        installation_folder = Path(installation_folder)
        call = [installation_folder.joinpath(command).as_posix()] + parsing.data_to_list(args)

    try:
        exit_code = run_subprocess(call + [version_command], shell=shell)
        assert exit_code == 0, f"call to {call[0]} exited with exit status {exit_code}."
    except FileNotFoundError:
        raise FileNotFoundError(f"RNAlysis could not find '{command}'. "
                                'Please ensure that your installation folder is correct, or add it to PATH. ')

    return call


def get_gunzip_size(fn):
    size = 0
    with gzip.open(fn) as f:
        while True:
            data = f.read(8192)
            size += len(data)
            if not data:
                break
    return size
