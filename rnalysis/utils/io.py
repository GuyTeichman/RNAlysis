import asyncio
import concurrent.futures
import contextlib
import ftplib
import functools
import gzip
import hashlib
import inspect
import json
import os
import queue
import random
import re
import shlex
import shutil
import subprocess
import sys
import threading
import time
import typing
import warnings
from datetime import date, datetime
from functools import lru_cache
from io import StringIO
from itertools import chain
from pathlib import Path
from sys import executable
from typing import List, Set, Union, Iterable, Tuple, Dict, Any, Callable, Literal, NamedTuple
from urllib.parse import urlparse, parse_qs, urlencode

import aiohttp
import aiolimiter
import appdirs
import matplotlib.pyplot as plt
import nest_asyncio
import numpy as np
import pandas as pd
import requests
import tenacity
import yaml
from defusedxml import ElementTree
from requests.adapters import HTTPAdapter, Retry
from tqdm import tqdm

from rnalysis import __version__
from rnalysis.utils import parsing, validation


def get_datafile_dir() -> Path:
    bundle_dir = getattr(sys, '_MEIPASS', '')
    return Path(bundle_dir)


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
    todays_dir = cache_dir.joinpath(today)
    todays_dir.mkdir(parents=True, exist_ok=True)
    return todays_dir


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


def clear_directory(directory: Union[str, Path], skip_ok=False):
    directory = Path(directory)
    if not directory.exists():
        return

    for item in directory.iterdir():
        try:
            if item.is_file():
                item.unlink()
            elif item.is_dir():
                shutil.rmtree(item, ignore_errors=True)
        except PermissionError as e:
            if not skip_ok:
                raise e


def clear_cache():
    cache_dir = Path(appdirs.user_cache_dir('RNAlysis'))
    clear_directory(cache_dir)


def clear_gui_cache():
    directory = get_gui_cache_dir()
    clear_directory(directory, skip_ok=True)


def load_cached_gui_file(filename: Union[str, Path], load_as_obj: bool = True) -> Union[
    str, set, pd.DataFrame, bytes, None]:
    """
    Load a cached file from the GUI cache directory.

    :param filename: The name of the file to load.
    :type filename: str
    :param load_as_obj: Whether to load the file as an object or raw content as string. Defaults to True.
    :type load_as_obj: bool
    :return: The contents of the file, loaded as either a string, a Pandas DataFrame or a set.
    :rtype: str or pandas.DataFrame or set or Bytes or None
    """
    directory = get_gui_cache_dir()
    file_path = directory.joinpath(filename)
    if file_path.exists():
        if file_path.suffix in {'.csv', '.tsv', '.parquet'} and load_as_obj:
            return load_table(file_path, index_col=0)
        elif file_path.suffix in {'.txt'} and load_as_obj:
            with open(file_path) as f:
                return {item.strip() for item in f.readlines()}

        elif file_path.suffix in {'.R'} and load_as_obj:
            with open(file_path) as f:
                return f.read()
        else:
            with open(file_path, 'rb') as f:
                return f.read()
    else:
        return None


def cache_gui_file(item: Union[pd.DataFrame, set, str], filename: str):
    directory = get_gui_cache_dir()
    if not directory.exists():
        directory.mkdir(parents=True)
    file_path = directory.joinpath(filename)
    if isinstance(item, (pd.DataFrame, pd.Series)):
        save_table(item, file_path, index=True)
    elif isinstance(item, set):
        save_gene_set(item, file_path)
    elif isinstance(item, Path) and item.suffix == '.R':
        assert item.exists(), f"File '{item}' flagged for caching does not exist!"
        shutil.copy(item, file_path)
    elif isinstance(item, str):
        with open(file_path, 'w') as f:
            f.write(item)
    elif isinstance(item, plt.Figure):
        item.savefig(file_path)
    else:
        raise TypeError(type(item))


def check_changed_version():  # pragma: no cover
    data_dir = get_data_dir()
    data_dir.mkdir(parents=True, exist_ok=True)
    filename = data_dir.joinpath('latest_version.txt')
    if not filename.exists():
        ver = ''
    else:
        with open(filename) as f:
            ver = f.read()
    current_ver = __version__
    # update latest version to current version
    with open(filename, 'w') as f:
        f.write(current_ver)

    return ver != current_ver


class FileData(NamedTuple):
    filename: str
    item_name: str
    item_type: str
    item_property: dict
    item_id: int = None
    obj: Union[set, pd.DataFrame] = None


class PipelineData(NamedTuple):
    name: str
    content: str


class GUISessionManager:
    def __init__(self, session_filename: Union[str, Path]):
        self.session_filename = Path(session_filename)

    def save_session(self, file_data: List[FileData], pipeline_data: List[PipelineData],
                     report: dict, report_file_paths: Dict[str, Path]):
        self._prepare_session_folder()
        self._save_report_files_to_session(report_file_paths)
        session_data = self._create_session_data(file_data, pipeline_data, report, report_file_paths)
        self._save_files_to_session(file_data)
        self._save_pipelines_to_session(pipeline_data, session_data)
        self._write_session_data_to_file(session_data)
        self._archive_session_folder()

    def load_session(self) -> Tuple[List[FileData], List[PipelineData], Any]:
        self._unpack_session_archive()
        session_dir = self._get_session_dir()
        with open(session_dir.joinpath('session_data.yaml')) as f:
            session_data = yaml.safe_load(f)

        file_data, pipeline_data = self._process_session_data(session_data)
        shutil.rmtree(session_dir)
        return file_data, pipeline_data, session_data.get('session_report_data', {'report': None})[
            'report']

    def _prepare_session_folder(self):
        if self.session_filename.exists():
            shutil.rmtree(self.session_filename) if self.session_filename.is_dir() else self.session_filename.unlink()
        self.session_filename.mkdir(parents=True)

    def _create_session_data(self, file_data: List[FileData], pipeline_data: List[PipelineData], report: dict,
                             report_item_paths: dict) -> dict:
        return {
            'files': {file.filename: (file.item_name, file.item_type, file.item_property) for file in file_data},
            'pipelines': {},
            'metadata': self._create_metadata(len(file_data), len(pipeline_data),
                                              [file.filename for file in file_data]),
            'session_report_data': {'report': report,
                                    'item_paths': {k: v.as_posix() for k, v in report_item_paths.items()}}}

    def _create_metadata(self, n_files: int, n_pipelines: int, file_names: List[str]) -> dict:
        return {
            'creation_time': get_datetime(),
            'name': self.session_filename.stem,
            'n_tabs': n_files,
            'n_pipelines': n_pipelines,
            'tab_order': file_names,
            'rnalysis_version': __version__
        }

    def _save_files_to_session(self, file_data: List[FileData]):
        for file in file_data:
            shutil.move(Path(get_gui_cache_dir().joinpath(file.filename)),
                        self.session_filename.joinpath(file.filename))

    def _save_report_files_to_session(self, report_file_paths: Dict[str, Path]) -> Dict[str, Path]:
        cache_dir = get_gui_cache_dir()
        for item_id, file_path in report_file_paths.items():
            shutil.copy(cache_dir.joinpath(file_path), self.session_filename.joinpath(file_path))

    def _save_pipelines_to_session(self, pipeline_data: List[PipelineData], session_data: dict):
        for i, pipeline in enumerate(pipeline_data):
            pipeline_filename = self.session_filename.joinpath(f"pipeline_{i}.yaml")
            Path(pipeline_filename).write_text(pipeline.content)
            session_data['pipelines'][pipeline_filename.name] = pipeline.name

    def _write_session_data_to_file(self, session_data: dict):
        with open(self.session_filename.joinpath('session_data.yaml'), 'w') as f:
            yaml.safe_dump(session_data, f)

    def _archive_session_folder(self):
        shutil.make_archive(self.session_filename.with_suffix('').as_posix(), 'zip', self.session_filename)
        shutil.rmtree(self.session_filename)
        self.session_filename.with_suffix('.zip').replace(self.session_filename.with_suffix('.rnal'))

    def _unpack_session_archive(self):
        try:
            self.session_filename.with_suffix('.rnal').rename(self.session_filename.with_suffix('.rnal.zip'))
            shutil.unpack_archive(self.session_filename.with_suffix('.rnal.zip'), self._get_session_dir())
        finally:
            self.session_filename.with_suffix('.rnal.zip').rename(self.session_filename.with_suffix('.rnal'))

    def _get_session_dir(self) -> Path:
        session_dir = get_gui_cache_dir().joinpath(self.session_filename.stem)
        return session_dir

    def _process_session_data(self, session_data: dict) -> Tuple[List[FileData], List[PipelineData]]:
        file_data = []
        pipeline_data = []

        session_dir = self._get_session_dir()
        filenames = session_data['metadata'].get('tab_order', list(session_data['files'].keys()))
        for file_name in filenames:
            file_path = session_dir.joinpath(file_name)
            assert file_path.exists() and file_path.is_file()
            if len(session_data['files'][file_name]) == 3:  # support for legacy session files
                item_name, item_type, item_property = session_data['files'][file_name]
                item_id = None
            else:
                item_name, item_type, item_property, item_id = session_data['files'][file_name]

            obj = load_cached_gui_file(session_dir.joinpath(file_name))

            file_data.append(FileData(filename=file_name,
                                      item_name=item_name,
                                      item_type=item_type,
                                      item_property=item_property,
                                      item_id=item_id,
                                      obj=obj))

        for pipeline_filename, pipeline_name in session_data['pipelines'].items():
            pipeline_path = session_dir.joinpath(pipeline_filename)
            assert pipeline_path.exists() and pipeline_path.is_file()
            pipeline_data.append(PipelineData(name=pipeline_name, content=pipeline_path.read_text()))
        # handle report files, if any by copying them to the cache directory
        cache_dir = get_gui_cache_dir()
        for item_id, file_path in session_data.get('session_report_data', {'item_paths': {}})['item_paths'].items():
            new_pth = cache_dir.joinpath(file_path)
            shutil.copy(session_dir.joinpath(file_path), new_pth)

        return file_data, pipeline_data


def load_table(filename: Union[str, Path], index_col: int = None, drop_columns: Union[str, List[str]] = False,
               squeeze=False, comment: str = None, engine: Literal['pyarrow', 'auto'] = 'auto'):
    """
    loads a csv/parquet table into a pandas dataframe.

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
    filename = Path(filename)
    assert filename.exists() and filename.is_file(), f"File '{filename.as_posix()}' does not exist!"
    assert filename.suffix.lower() in {'.csv', '.tsv', '.txt', '.parquet'}, \
        f"RNAlysis cannot load files of type '{filename.suffix}'. " \
        f"Please convert your file to a .csv, .tsv, .txt, or .parquet file and try again."

    if filename.suffix.lower() == '.parquet':
        df = pd.read_parquet(filename, engine=engine)
    else:
        kwargs = dict(sep=None, engine='python' if engine == 'auto' else engine, encoding='ISO-8859-1', comment=comment,
                      skipinitialspace=True)
        if index_col is not None:
            kwargs['index_col'] = index_col
        df = pd.read_csv(filename, **kwargs)

    if squeeze:
        df = df.squeeze("columns")

    if index_col is not None:
        df.index = df.index.astype('str')
    df.index = [ind.strip() if isinstance(ind, str) else ind for ind in df.index]
    if isinstance(df, pd.DataFrame):
        df.columns = [col.strip() if isinstance(col, str) else str(col).strip() for col in df.columns]

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


def save_table(df: pd.DataFrame, filename: Union[str, Path], postfix: str = None, index: bool = True):
    """
    save a pandas DataFrame to csv/parquet file.

    :param df: pandas DataFrame to be saved
    :param filename: a string or pathlib.Path object stating the original name of the file
    :type postfix: str, default None
    :param postfix: A postfix to be added to the original name of the file. If None, no postfix will be added.
    :param index: if True, saves the DataFrame with the indices. If false, ignores the index.
    """
    fname = Path(filename)
    if postfix is None:
        postfix = ''
    else:
        assert isinstance(postfix, str), "'postfix' must be either str or None!"
    new_fname = os.path.join(fname.parent.absolute(), f"{fname.stem}{postfix}{fname.suffix}")
    if fname.suffix.lower() == '.parquet':
        if isinstance(df, pd.Series):
            df = df.to_frame()
        df.to_parquet(new_fname, index=index)
    else:
        df.to_csv(new_fname, header=True, index=index)


def get_session(retries: Retry):
    session = requests.Session()
    adapter = HTTPAdapter(max_retries=retries)
    session.mount("http://", adapter)
    session.mount("https://", adapter)
    return session


class KEGGAnnotationIterator:
    URL = 'https://rest.kegg.jp/'
    TAXON_MAPPING_URL = 'https://www.genome.jp/kegg-bin/download_htext?htext=br08610'
    REQUEST_DELAY_MILLIS = 250
    REQ_MAX_ENTRIES = 10
    RETRIES = Retry(total=5, backoff_factor=0.25, status_forcelist=[500, 502, 503, 504])
    TAXON_TREE_CACHED_FILENAME = 'kegg_taxon_tree.json'
    PATHWAY_NAMES_CACHED_FILENAME = 'kegg_pathway_list.txt'
    COMPOUND_LIST_CACHED_FILENAME = 'kegg_compound_list.txt'
    GLYCAN_LIST_CACHED_FILENAME = 'kegg_glycan_list.txt'

    def __init__(self, taxon_id: int, pathways: Union[str, List[str], Literal['all']] = 'all'):
        self.pathway_names = {}
        self.taxon_id = taxon_id
        self.session = get_session(self.RETRIES)
        self.organism_code = self.get_kegg_organism_code(taxon_id, self.session)
        if pathways == 'all':
            self.pathway_names, self.n_annotations = self.get_pathways()
        else:
            pathways = parsing.data_to_list(pathways)
            assert len(pathways) > 0, "No KEGG pathway IDs were given!"
            self.pathway_names = {pathway: None for pathway in pathways}
            self.n_annotations = len(self.pathway_names)
        self.pathway_annotations = None

    @staticmethod
    def _kegg_request(session: requests.Session, operation: str, arguments: Union[str, List[str]],
                      cached_filename: Union[str, None] = None, ) -> Tuple[str, bool]:
        if cached_filename is not None:
            cached_file = load_cached_file(cached_filename)
            if cached_file is not None:
                is_cached = True
                return cached_file, is_cached

        is_cached = False
        address = KEGGAnnotationIterator.URL + f'{operation}/' + '/'.join(parsing.data_to_list(arguments))
        response = session.get(address)
        if not response.ok:
            response.raise_for_status()
        if cached_filename is not None:
            cache_file(response.text, cached_filename)
        return response.text, is_cached

    def _generate_cached_filename(self, pathways: Tuple[str, ...]) -> str:
        fname = f'{self.taxon_id}' + ''.join(pathways).replace('path:', '') + '.json'
        return fname

    @staticmethod
    def _get_taxon_tree(session):
        cached_filename = KEGGAnnotationIterator.TAXON_TREE_CACHED_FILENAME
        cached_file = load_cached_file(cached_filename)
        if cached_file is not None:
            try:
                taxon_tree = json.loads(cached_file)
                return taxon_tree
            except json.decoder.JSONDecodeError:
                pass
        with session.get(KEGGAnnotationIterator.TAXON_MAPPING_URL, params=dict(format='json')) as req:
            content = req.content.decode('utf8')
            cache_file(content, cached_filename)
            taxon_tree = json.loads(content)
        return taxon_tree

    @staticmethod
    def get_compounds() -> Dict[str, str]:
        compounds = {}
        session = get_session(KEGGAnnotationIterator.RETRIES)
        data = KEGGAnnotationIterator._kegg_request(session, 'list', ['compound'],
                                                    KEGGAnnotationIterator.COMPOUND_LIST_CACHED_FILENAME)[0] + '\n' + \
               KEGGAnnotationIterator._kegg_request(session, 'list', ['glycan'],
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
    def get_kegg_organism_code(taxon_id: int, session: requests.Session) -> str:
        taxon_tree = KEGGAnnotationIterator._get_taxon_tree(session)
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
        data, _ = self._kegg_request(self.session, 'list', ['pathway', self.organism_code],
                                     self.PATHWAY_NAMES_CACHED_FILENAME)
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
        session = get_session(KEGGAnnotationIterator.RETRIES)
        data, _ = KEGGAnnotationIterator._kegg_request(session, 'get', [pathway_id, 'kgml'], cached_filename)
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
                data, was_cached = self._kegg_request(self.session, 'get', '+'.join(chunk),
                                                      self._generate_cached_filename(chunk))
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
                 'n_annotations': 'number of annotations matching the filtering criteria',
                 'session': 'session'}
    URL = 'http://golr-aux.geneontology.io/solr/select?'
    RETRIES = Retry(total=5, backoff_factor=0.25, status_forcelist=[500, 502, 503, 504])

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

        self.session = get_session(self.RETRIES)
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

    def _golr_request(self, params: dict, cached_filename: Union[str, None] = None) -> str:
        """
        Run a get request to the GOlr server with the specified parameters, and return the server's text response.
        :param params: the get request's parameters.
        :type params: dict
        """
        if cached_filename is not None:
            cached_file = load_cached_file(cached_filename)
            if cached_file is not None:
                return cached_file

        response = self.session.get(GOlrAnnotationIterator.URL, params=params)
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
        query = ['document_category:"annotation"',
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
def _ensembl_lookup_post_request(gene_ids: Tuple[str]) -> Dict[str, Dict[str, Any]]:
    """
    Perform an Ensembl 'lookup/id' POST request to find the species and database for several identifiers. \
    See full POST API at https://rest.ensembl.org/documentation/info/lookup_post

    :param gene_ids: a tuple of gene IDs to look up
    :type gene_ids: tuple of str
    :return: a dictionary with gene IDs as keys and dictionaries of attributes as values
    """
    endpoint = 'lookup/id'
    req_type = 'post'
    # split the gene IDs into chunks of 1000 (the maximum allowed POST request size)
    data_chunks = parsing.partition_list(gene_ids, 1000)
    output = {}
    client = EnsemblRestClient()
    for chunk in tqdm(data_chunks, desc='Submitting jobs...', unit='jobs'):
        data = {"ids": parsing.data_to_list(chunk)}
        client.queue_action(req_type, endpoint, params=repr(data).replace("'", '"'))

    with tqdm('Finding the best-matching species...', total=client.queue.qsize() + 1) as pbar:
        pbar.update()
        for json_res in client.run():
            output.update(json_res)
            pbar.update()
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
    output = _ensembl_lookup_post_request(parsing.data_to_tuple(translator.mapping_dict.values()))
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
    output = _ensembl_lookup_post_request(parsing.data_to_tuple(translator.mapping_dict.values()))
    species = dict()
    for gene_id in output:
        if output[gene_id] is not None:
            species[output[gene_id]['species']] = species.setdefault(output[gene_id]['species'], 0) + 1
    if len(species) == 0:
        raise ValueError("No taxon ID could be matched to any of the given gene IDs.")
    chosen_species = list(species.keys())[0]
    chosen_species_n = species[chosen_species]
    if len(species) > 1:
        warnings.warn("The given gene IDs match more than one species. "
                      "Picking the species that fits the majority of gene IDs.")
        for s in species:
            if species[s] > chosen_species_n:
                chosen_species = s
                chosen_species_n = species[s]
    return map_taxon_id(chosen_species.replace('_', ' ')), map_from


def get_taxon_and_id_type(organism, gene_id_type, gene_set, map_to_options=('Ensembl', 'Ensembl Genomes')):
    """
    Determines the taxon and gene id type for a given organism and set of genes. \
    If 'auto' is provided as the organism, the taxon will be inferred from the gene set. If 'auto' is provided as the \
    gene id type, the function will either use the inferred type (if available) or find the best gene mapping.

    :param organism: The organism to get the taxon and id type for. If 'auto', the taxon will \
    be inferred from the gene set. This can be of type str, int, or Literal['auto'].
    :type organism: Union[str, int, Literal['auto']]
    :param gene_id_type: The gene id type. If 'auto', the id type will be inferred or the best \
    gene mapping will be found. This can be of type str or Literal['auto'].
    :type gene_id_type: Union[str, Literal['auto']]
    :param gene_set: A set of gene ids.
    :type gene_set: Set[str]
    :param map_to_options: The options for mapping to. Defaults to ('Ensembl', 'Ensembl Genomes'). \
    This can be a single string or a tuple of strings.
    :type map_to_options: Union[str, Tuple[str]], optional
    :return: A tuple containing the taxon and the gene id type.
    :rtype: tuple
    """
    # Convert map_to_options to tuple if necessary
    map_to_options = parsing.data_to_tuple(map_to_options)

    # If organism is 'auto', infer taxon from gene IDs
    if isinstance(organism, str) and organism.lower() == 'auto':
        id_type = None if gene_id_type.lower() == 'auto' else gene_id_type
        taxon, map_from = infer_taxon_from_gene_ids(gene_set, id_type)

        # If gene_id_type is also 'auto', use inferred type
        if gene_id_type.lower() == 'auto':
            gene_id_type = map_from

    else:
        # If organism is known and gene_id_type is 'auto', find the best gene mapping
        if gene_id_type.lower() == 'auto':
            _, gene_id_type, _ = find_best_gene_mapping(
                parsing.data_to_tuple(gene_set),
                map_from_options=None,
                map_to_options=map_to_options)

        # Map taxon ID
        taxon = map_taxon_id(organism)

    return taxon, gene_id_type


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


class OrthologDict:
    """
    A dictionary-like class used for mapping genes to their ortholog in a different organism.
    """

    __slots__ = {'mapping_dict': 'dictionary mapping gene IDs from one species to another based on orthology'}

    def __init__(self, mapping_dict: Union[dict, None] = None):
        """
        :param mapping_dict: a dictionary mapping gene IDs from one type to another. \
        If mappping_dict is None, gene IDs will be automatically mapped to themselves.
        :type mapping_dict: dict or None (default None)
        """
        if mapping_dict is None:
            mapping_dict = {}
        self.mapping_dict = mapping_dict

    def __len__(self):
        return len(self.mapping_dict)

    def __getitem__(self, key):
        return self.mapping_dict[key]

    def __contains__(self, item):
        try:
            _ = self[item]
            return True
        except KeyError:
            return False


class GeneIDDict:
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
        if key in self.mapping_dict:
            return self.mapping_dict[key]
        raise KeyError(f"Gene ID '{key}' not found in mapping dictionary.")

    def __contains__(self, item):
        try:
            _ = self[item]
            return True
        except KeyError:
            return False


class EnsemblRestClient:
    SERVER = 'https://rest.ensembl.org/'
    REQS_PER_SEQ = 10
    HEADERS = {"Content-Type": "application/json", "Accept": "application/json"}
    LIMITER = aiolimiter.AsyncLimiter(REQS_PER_SEQ, 1)

    def __init__(self):
        self.queue = queue.Queue()
        self.session = None
        self.semaphore = None

    def queue_action(self, req_type: Literal['get', 'post'], endpoint: str, hdrs=None, params=None):
        self.queue.put((req_type, endpoint, hdrs, params))

    def run(self):
        try:
            res = asyncio.run(self._run())
        except RuntimeError:
            nest_asyncio.apply()
            res = asyncio.run(self._run())
        return res

    async def _run(self):
        tasks = []
        async with aiohttp.ClientSession(raise_for_status=True) as self.session:
            while not self.queue.empty():
                req_type, endpoint, hdrs, params = self.queue.get()
                tasks.append(self.perform_api_action(req_type, endpoint, hdrs, params))

            return await asyncio.gather(*tasks)

    @tenacity.retry(stop=tenacity.stop_after_attempt(5),
                    wait=tenacity.wait_random_exponential(multiplier=1, max=10),
                    retry=tenacity.retry_if_exception_type(
                        (aiohttp.ClientConnectorError, aiohttp.ClientResponseError, asyncio.TimeoutError)))
    async def perform_api_action(self, req_type: Literal['get', 'post'], endpoint: str, hdrs=None, params=None):
        if self.semaphore is None:
            self.semaphore = asyncio.Semaphore(value=self.REQS_PER_SEQ)

        if hdrs is None:
            hdrs = self.HEADERS

        if 'Content-Type' not in hdrs:
            hdrs['Content-Type'] = 'application/json'

        content = None
        try:
            await self.semaphore.acquire()
            async with self.LIMITER:
                if req_type == 'get':
                    request_func = self.session.get
                    kwargs = dict(params=params) if params else {}
                elif req_type == 'post':
                    request_func = self.session.post
                    kwargs = dict(data=params) if params else {}
                else:
                    raise ValueError(f"Invalid request type '{req_type}'. ")
                async with request_func(self.SERVER + endpoint, headers=hdrs, **kwargs) as response:
                    content = await response.json()
                    self.semaphore.release()

        except (aiohttp.ClientConnectorError, asyncio.TimeoutError) as e:
            # check if we are being rate limited by the server
            if e.errno == 429:
                if 'Retry-After' in e.request.headers:
                    retry = e.request.headers['Retry-After']
                    await asyncio.sleep(float(retry))
                    raise e  # wait the minimum amount of time required by 'retry-after',
                    # and then trigger tenacity.retry for random exponential to avoid thundering herd
            else:
                raise e
        return content


def translate_mappings(ids: list, translated_ids: list, mapping_one2one: dict,
                       mapping_one2many: dict):
    requested = parsing.data_to_set(translated_ids)
    mapping_one2one_translated = {}
    mapping_one2many_translated = {}
    for from_id, trans_from_id in zip(ids, translated_ids):
        if trans_from_id in mapping_one2one and trans_from_id in requested:
            mapping_one2one_translated[from_id] = mapping_one2one[trans_from_id]

        if trans_from_id in mapping_one2many and trans_from_id in requested:
            mapping_one2many_translated[from_id] = mapping_one2many[trans_from_id]

    return mapping_one2one_translated, mapping_one2many_translated


class PhylomeDBOrthologMapper:
    URL = 'ftp.phylomedb.org'

    def __init__(self, map_to_organism, map_from_organism='auto', gene_id_type='auto'):
        legal_species = self.get_legal_species()
        assert map_from_organism in legal_species.index, f"organism with taxon id {map_from_organism} is not supported by PhylomeDB. "
        assert map_to_organism in legal_species.index, f"organism with taxon id {map_to_organism} is not supported by PhylomeDB. "
        self.gene_id_type = gene_id_type
        self.map_from_organism = map_from_organism
        self.map_to_organism = map_to_organism

    def translate_ids(self, ids: Tuple[str, ...]) -> Tuple[List[str], List[str]]:
        if self.gene_id_type == 'auto':
            translator, self.gene_id_type, _ = find_best_gene_mapping(ids, None, map_to_options=('UniProtKB AC/ID',))
        else:
            translator = GeneIDTranslator(self.gene_id_type, 'UniProtKB AC/ID').run(ids)

        ids = [this_id for this_id in ids if this_id in translator]
        translated_ids = [translator[this_id] for this_id in ids]

        return ids, translated_ids

    def get_orthologs(self, ids: Tuple[str, ...], non_unique_mode: str, consistency_score_threshold: float,
                      filter_consistency_score: bool = True):
        mapping_one2one = {}
        mapping_one2many = {}

        taxon_table = self._get_taxon_file(self.map_from_organism)
        taxon_table = taxon_table[taxon_table['taxid2'] == self.map_to_organism]
        map_fwd, map_rev = self._get_id_conversion_maps()
        ids, translated_ids = self.translate_ids(ids)

        n_mapped = 0
        for from_id in tqdm(translated_ids, 'Mapping orthologs', unit='genes'):
            if from_id not in map_fwd.index:
                continue
            from_id_conv = map_fwd.loc[from_id]
            if from_id_conv not in taxon_table.index:
                continue
            to_id_conv = taxon_table.loc[from_id_conv, 'protid2']
            if isinstance(to_id_conv, pd.Series):
                for i, this_to_id_conv in enumerate(to_id_conv):
                    if this_to_id_conv not in map_rev.index:
                        continue

                    if from_id not in mapping_one2many:
                        mapping_one2many[from_id] = []

                    to_id = map_rev.loc[this_to_id_conv]
                    score = taxon_table.loc[from_id_conv, 'CS'].values[0]
                    # filter by consistency score
                    if score < consistency_score_threshold:
                        continue

                    mapping_one2many[from_id] = [(to_id, score)]

                # skip genes where all orthologs were filtered out due to consistency score threshold
                if len(mapping_one2many.get(from_id, list())) == 0:
                    continue

                if filter_consistency_score:
                    mapping_one2one[from_id] = max(mapping_one2many[from_id], key=lambda x: x[1])[0]
                else:
                    if non_unique_mode == 'first':
                        mapping_one2one[from_id] = mapping_one2many[from_id][0][0]
                    elif non_unique_mode == 'last':
                        mapping_one2one[from_id] = mapping_one2many[from_id][-1][0]
                    elif non_unique_mode == 'random':
                        mapping_one2one[from_id] = random.choice(mapping_one2many[from_id])[0]
                n_mapped += 1
            else:
                if to_id_conv not in map_rev.index:
                    continue
                to_id = map_rev.loc[to_id_conv]
                score = taxon_table.loc[from_id_conv, 'CS']
                # filter by consistency score
                if score < consistency_score_threshold:
                    continue

                mapping_one2one[from_id] = to_id
                mapping_one2many[from_id] = [(to_id, score)]
                n_mapped += 1

        mapping_one2one, mapping_one2many = translate_mappings(ids, translated_ids, mapping_one2one, mapping_one2many)

        if n_mapped < len(translated_ids):
            warnings.warn(f"Ortholog mapping found for only {n_mapped} out of {len(translated_ids)} gene IDs.")

        return OrthologDict(mapping_one2one), OrthologDict(
            {k: [this_v[0] for this_v in v] for k, v in mapping_one2many.items()})

    @staticmethod
    def _get_taxon_file(taxon_id: int):
        cache_dir = get_todays_cache_dir()
        cache_file = cache_dir.joinpath(f'phylomedb_{taxon_id}.parquet')
        file_path = f"/metaphors/latest/orthologs/{taxon_id}.txt.gz"
        local_zip_file = cache_dir.joinpath(f"{taxon_id}.txt.gz")
        if cache_file.exists():
            df = load_table(cache_file, squeeze=True, index_col=0, engine="pyarrow")
        else:
            ftp = PhylomeDBOrthologMapper._connect()
            # Download the file
            with open(local_zip_file, 'wb') as fp:
                ftp.retrbinary('RETR ' + file_path, fp.write)
            ftp.quit()

            # Unzip the file
            with gzip.open(local_zip_file, 'rb') as f_in:
                with open(cache_file.with_suffix('.txt'), 'wb') as f_out:
                    shutil.copyfileobj(f_in, f_out)

            # Remove the zipped file
            os.remove(local_zip_file)

            # Load into pandas dataframe
            df = pd.read_csv(cache_file.with_suffix('.txt'), sep='\t', index_col=1, engine="pyarrow")

            # cache file locally
            save_table(df, cache_file)
        return df

    @staticmethod
    def _get_id_conversion_maps() -> Tuple[pd.DataFrame, pd.DataFrame]:
        cache_dir = get_todays_cache_dir()
        cache_file = cache_dir.joinpath(f'phylomedb_id_conversion.parquet')
        file_path = "/metaphors/latest/id_conversion.txt.gz"
        local_zip_file = cache_dir.joinpath("id_conversion.txt.gz")

        if cache_file.exists():
            df = pd.read_parquet(cache_file, engine='auto', dtype_backend='pyarrow')
        else:
            ftp = PhylomeDBOrthologMapper._connect()
            # Download the file
            with open(local_zip_file, 'wb') as fp:
                ftp.retrbinary('RETR ' + file_path, fp.write)
            ftp.quit()

            # Unzip the file
            with gzip.open(local_zip_file, 'rb') as f_in:
                with open(cache_file.with_suffix('.txt'), 'wb') as f_out:
                    shutil.copyfileobj(f_in, f_out)

            # Remove the zipped file
            os.remove(local_zip_file)

            # Load into pandas dataframe
            df = pd.read_csv(cache_file.with_suffix('.txt'), index_col=0, sep='\t', dtype_backend='pyarrow',
                             usecols=[0, 2], names=['#extid', 'protid'], header=None, skiprows=1)
            # cache file locally
            save_table(df, cache_file)

        df_inv = df.reset_index().set_index('protid')
        return df.squeeze(), df_inv.squeeze()

    @staticmethod
    def _connect():
        ftp = ftplib.FTP(PhylomeDBOrthologMapper.URL)
        ftp.login()
        ftp.cwd('/metaphors/latest')
        return ftp

    @staticmethod
    @lru_cache(maxsize=2)
    def get_legal_species():
        cache_dir = get_todays_cache_dir()
        file_path = "/metaphors/latest/species.txt.gz"
        local_zip_file = cache_dir.joinpath("species.txt.gz")
        cache_file = cache_dir.joinpath("phylomedb_species.parquet")

        if cache_file.exists():
            df = load_table(cache_file, squeeze=True, index_col=0, engine='pyarrow')
            df.index = df.index.astype(int)
        else:
            ftp = PhylomeDBOrthologMapper._connect()
            # Download the file
            with open(local_zip_file, 'wb') as fp:
                ftp.retrbinary('RETR ' + file_path, fp.write)
            ftp.quit()

            # Unzip the file
            with gzip.open(local_zip_file, 'rb') as f_in:
                with open(cache_file.with_suffix('.txt'), 'wb') as f_out:
                    shutil.copyfileobj(f_in, f_out)

            # Remove the zipped file
            os.remove(local_zip_file)

            # Load into pandas dataframe
            df = pd.read_csv(cache_file.with_suffix('.txt'), sep='\t', index_col=0, engine='pyarrow').dropna()

            # Extract necessary columns and sort by name
            df.index = df.index.astype(int)
            df = df['name'].sort_index()

            # cache file locally
            save_table(df, cache_file)
        return df


class OrthoInspectorOrthologMapper:
    API_URL = 'https://lbgi.fr/api/orthoinspector'
    RETRIES = Retry(total=5, backoff_factor=0.25, status_forcelist=[500, 502, 503, 504])
    HEADERS = {'accept': 'application/json'}

    def __init__(self, map_to_organism, map_from_organism, gene_id_type):
        self.gene_id_type = gene_id_type
        self.map_from_organism = map_from_organism
        self.map_to_organism = map_to_organism

    def translate_ids(self, ids: Tuple[str, ...], session=None) -> Tuple[List[str], List[str]]:
        if self.gene_id_type == 'auto':
            translator, self.gene_id_type, _ = find_best_gene_mapping(ids, None, map_to_options=('UniProtKB AC/ID',))
        else:
            translator = GeneIDTranslator(self.gene_id_type, 'UniProtKB AC/ID', session=session).run(ids)

        ids = [this_id for this_id in ids if this_id in translator]
        translated_ids = [translator[this_id] for this_id in ids]

        return ids, translated_ids

    @staticmethod
    def get_databases(session: requests.Session = None):
        if session is None:
            session = get_session(OrthoInspectorOrthologMapper.RETRIES)
        url = f'{OrthoInspectorOrthologMapper.API_URL}/databases'
        req = session.get(url)
        req.raise_for_status()
        databases = frozenset(req.json()['data'])
        return databases

    @staticmethod
    def get_database_organisms(session: requests.Session = None):
        if session is None:
            session = get_session(OrthoInspectorOrthologMapper.RETRIES)
        # get list of databases
        databases = OrthoInspectorOrthologMapper.get_databases(session)
        db_organisms = {}
        for i, database in enumerate(databases):
            url = f'{OrthoInspectorOrthologMapper.API_URL}/{database}/species'
            req = session.get(url)
            req.raise_for_status()
            content = req.json()
            assert content['meta']['status'] == 'success'
            db_organisms[database] = frozenset({d['id'] for d in content['data']})
            assert len(db_organisms[database]) == int(content['meta']['nbResults'])

        return db_organisms

    def get_cache_filename(self):
        return f'orthoinspector_{self.map_from_organism}_{self.map_to_organism}.json'

    def get_orthologs(self, ids: Tuple[str, ...], non_unique_mode: str, database: str = 'auto'):
        session = get_session(self.RETRIES)
        # find a valid database, or ensure that given database is valid
        if database == 'auto':
            database_organisms = OrthoInspectorOrthologMapper.get_database_organisms()
            dbs_by_size = sorted(database_organisms.keys(), key=lambda x: len(database_organisms[x]), reverse=True)
            valid_dbs = []

            for db in dbs_by_size:
                if self.map_from_organism in database_organisms[db] and self.map_to_organism in database_organisms[db]:
                    valid_dbs.append(db)
            if len(valid_dbs) == 0:
                raise ValueError('No database found that supports mapping from '
                                 f'{self.map_from_organism} to {self.map_to_organism}. ')

        else:
            databases = self.get_databases()
            assert database in databases, f"Invalid database: {database}. Valid databases are: {databases}."
            valid_dbs = [database]

        mapping_one2one = {}
        mapping_one2many = {}
        ids, translated_ids = self.translate_ids(ids, session=session)
        # TODO: progress bars on both download and parsing
        # if a large number of genes is requested, download the entire pairwise dataset and filter it
        cached = load_cached_file(self.get_cache_filename())
        if cached is None:
            content = {}
            for database in valid_dbs:
                url = f'{self.API_URL}/{database}/species/{self.map_from_organism}/orthologs/{self.map_to_organism}'
                req = session.get(url)
                req.raise_for_status()
                content = req.json()['data']
                if len(content) >= 0:
                    cache_file(json.dumps(content), self.get_cache_filename())
                    break
        else:
            content = json.loads(cached)

        for annotation in content:
            rel_type = annotation['type'].lower()
            if rel_type == 'one-to-one':
                to_id = annotation['orthologs'][0]
                from_id = annotation['inparalogs'][0]
                mapping_one2many[from_id] = [to_id]
                mapping_one2one[from_id] = to_id
            elif rel_type == 'one-to-many':
                to_ids = annotation['orthologs']
                from_id = annotation['inparalogs'][0]
                mapping_one2many[from_id] = to_ids
                if non_unique_mode == 'first':
                    mapping_one2one[from_id] = to_ids[0]
                elif non_unique_mode == 'last':
                    mapping_one2one[from_id] = to_ids[-1]
                elif non_unique_mode == 'random':
                    mapping_one2one[from_id] = random.choice(to_ids)

            elif rel_type == 'many-to-one':
                to_id = annotation['orthologs'][0]
                from_ids = annotation['inparalogs']
                for from_id in from_ids:
                    mapping_one2many[from_id] = [to_id]
                    mapping_one2one[from_id] = to_id
            elif rel_type == 'many-to-many':
                to_ids = annotation['orthologs']
                from_ids = annotation['inparalogs']
                for from_id in from_ids:
                    mapping_one2many[from_id] = to_ids
                    if non_unique_mode == 'first':
                        mapping_one2one[from_id] = to_ids[0]
                    elif non_unique_mode == 'last':
                        mapping_one2one[from_id] = to_ids[-1]
                    elif non_unique_mode == 'random':
                        mapping_one2one[from_id] = random.choice(to_ids)

            else:
                raise ValueError(f'Unknown ortholog type: "{annotation["type"]}"')

        mapping_one2one, mapping_one2many = translate_mappings(ids, translated_ids, mapping_one2one, mapping_one2many)
        n_mapped = len(mapping_one2one)
        if n_mapped < len(translated_ids):
            warnings.warn(f"Ortholog mapping found for only {n_mapped} out of {len(translated_ids)} gene IDs.")

        return OrthologDict(mapping_one2one), OrthologDict(mapping_one2many)


class EnsemblOrthologMapper:
    ENDPOINT = 'homology/id/'

    def __init__(self, map_to_organism, map_from_organism, gene_id_type):
        self.gene_id_type = gene_id_type
        self.map_from_organism = map_from_organism
        self.map_to_organism = map_to_organism

    def get_species_name(self):
        return map_taxon_id(self.map_from_organism)[1].replace(' ', '_')

    def translate_ids(self, ids: Tuple[str, ...]) -> Tuple[List[str], List[str]]:
        if self.gene_id_type == 'auto':
            translator, self.gene_id_type, _ = find_best_gene_mapping(ids, None, map_to_options=('Ensembl Genomes',))
        else:
            translator = GeneIDTranslator(self.gene_id_type, 'Ensembl Genomes').run(ids)

        ids = [this_id for this_id in ids if this_id in translator]
        translated_ids = [translator[this_id] for this_id in ids]

        return ids, translated_ids

    def get_paralogs(self, ids: Tuple[str, ...], filter_percent_identity: bool = True):
        ids, translated_ids = self.translate_ids(ids)
        client = EnsemblRestClient()
        mapping_one2many = {}
        species_name = self.get_species_name()

        for gene_id in tqdm(translated_ids, 'Submitting requests', unit='requests'):
            client.queue_action('get', f'{self.ENDPOINT}{species_name}/{gene_id}',
                                params=dict(target_taxon=self.map_to_organism, type='paralogues',
                                            sequence='none', cigar_line=0))

        with tqdm('Mapping paralogs...', total=client.queue.qsize() + 1) as pbar:
            pbar.update()
            for json_res in client.run():
                req_output = json_res['data'][0]
                if len(req_output['homologies']) == 0:
                    continue

                this_id = req_output['id']

                if filter_percent_identity:
                    max_score = 0
                    best_match = None
                    matches = []
                    for homology in req_output['homologies']:
                        current_score = homology['source']['perc_id']
                        if current_score > max_score:
                            max_score = current_score
                            best_match = homology['target']['id']
                            matches = [best_match]
                        elif current_score == max_score:
                            matches.append(homology['target']['id'])

                    mapping_one2many[this_id] = best_match
                else:
                    mapping_one2many[this_id] = [homology['target']['id'] for homology in req_output['homologies']]
                pbar.update()

        _, mapping_one2many = translate_mappings(ids, translated_ids, {}, mapping_one2many)

        n_mapped = len(mapping_one2many)
        if n_mapped < len(translated_ids):
            warnings.warn(f"Paralog mapping found for only {n_mapped} out of {len(translated_ids)} gene IDs.")

        return OrthologDict(mapping_one2many)

    def get_orthologs(self, ids: Tuple[str, ...], non_unique_mode: str, filter_percent_identity: bool = True) -> \
        Tuple[OrthologDict, OrthologDict]:
        ids, translated_ids = self.translate_ids(ids)
        species_name = self.get_species_name()
        client = EnsemblRestClient()
        mapping_one2one = {}
        mapping_one2many = {}

        for gene_id in tqdm(translated_ids, 'Submitting requests', unit='requests'):
            client.queue_action('get', f'{self.ENDPOINT}{species_name}/{gene_id}',
                                params=dict(target_taxon=self.map_to_organism, type='orthologues',
                                            sequence='none', cigar_line=0))

        with tqdm('Mapping orthologs...', total=client.queue.qsize() + 1) as pbar:
            pbar.update()
            for json_res in client.run():
                req_output = json_res['data'][0]
                print(req_output)
                if len(req_output['homologies']) == 0:
                    continue

                this_id = req_output['id']
                mapping_one2many[this_id] = [(homology['target']['id'], homology['source']['perc_id']) for homology in
                                             req_output['homologies']]

                if filter_percent_identity:
                    max_score = 0
                    best_match = None
                    matches = []
                    for homology in req_output['homologies']:
                        current_score = homology['source']['perc_id']
                        if current_score > max_score:
                            max_score = current_score
                            best_match = homology['target']['id']
                            matches = [best_match]
                        elif current_score == max_score:
                            matches.append(homology['target']['id'])

                    if len(matches) > 1:
                        if non_unique_mode == 'first':
                            mapping_one2one[this_id] = matches[0]
                        elif non_unique_mode == 'last':
                            mapping_one2one[this_id] = matches[-1]
                        elif non_unique_mode == 'random':
                            mapping_one2one[this_id] = random.choice(matches)
                    else:
                        mapping_one2one[this_id] = best_match
                else:
                    mapping_one2one[this_id] = req_output['homologies'][0]['target']['id']
                pbar.update()

        mapping_one2one, mapping_one2many = translate_mappings(ids, translated_ids, mapping_one2one, mapping_one2many)

        n_mapped = len(mapping_one2many)
        if n_mapped < len(translated_ids):
            warnings.warn(f"Ortholog mapping found for only {n_mapped} out of {len(translated_ids)} gene IDs.")

        return OrthologDict(mapping_one2one), OrthologDict(
            {k: [this_v[0] for this_v in v] for k, v in mapping_one2many.items()})


class PantherOrthologMapper:
    API_URL = 'http://www.pantherdb.org'
    RETRIES = Retry(total=5, backoff_factor=0.25, status_forcelist=[500, 502, 503, 504])
    HEADERS = {'accept': 'application/json'}

    def __init__(self, map_to_organism, map_from_organism='auto', gene_id_type='auto'):
        self.gene_id_type = gene_id_type
        self.map_from_organism = map_from_organism
        self.map_to_organism = map_to_organism

    def translate_ids(self, ids: Tuple[str, ...], session=None) -> Tuple[List[str], List[str]]:
        if self.gene_id_type == 'auto':
            translator, self.gene_id_type, _ = find_best_gene_mapping(ids, None, map_to_options=('UniProtKB AC/ID',))
        else:
            translator = GeneIDTranslator(self.gene_id_type, 'UniProtKB AC/ID', session=session).run(ids)

        ids = [this_id for this_id in ids if this_id in translator]
        translated_ids = [translator[this_id] for this_id in ids]

        return ids, translated_ids

    def get_orthologs(self, ids: Tuple[str, ...], non_unique_mode: str, filter_least_diverged: bool = True):
        url = f'{self.API_URL}/services/oai/pantherdb/ortholog/matchortho'
        session = get_session(self.RETRIES)
        mapping_one2one = {}
        mapping_one2many = {}

        ids, translated_ids = self.translate_ids(ids, session)
        n_mapped = 0
        for from_id in tqdm(translated_ids, 'Mapping orthologs', unit='genes'):
            req_data = dict(geneInputList=from_id, organism=str(self.map_from_organism),
                            targetOrganism=str(self.map_to_organism),
                            orthologType='all' if filter_least_diverged else 'LDO')
            req = session.post(url, headers=self.HEADERS, params=req_data)
            req.raise_for_status()
            try:
                req_output = req.json()['search']['mapping']['mapped']
                for mapping in parsing.data_to_list(req_output):
                    if len(mapping) <= 1:
                        continue

                    to_id = parsing.parse_uniprot_id(mapping['target_gene'])
                    if to_id is None:
                        continue

                    if from_id not in mapping_one2many:
                        n_mapped += 1
                        mapping_one2many[from_id] = []
                    mapping_one2many[from_id].append((to_id, mapping['ortholog']))

            except KeyError:
                pass

        for from_id, mappings in mapping_one2many.items():
            if filter_least_diverged:
                mappings_filtered = [this_mapping for this_mapping in mappings if this_mapping[1] == 'LDO']
                mappings = mappings_filtered if len(mappings_filtered) > 0 else mappings
            if len(mappings) == 1:
                mapping_one2one[from_id] = mappings[0][0]
            else:
                if non_unique_mode == 'first':
                    mapping_one2one[from_id] = mappings[0][0]
                elif non_unique_mode == 'last':
                    mapping_one2one[from_id] = mappings[-1][0]
                elif non_unique_mode == 'random':
                    mapping_one2one[from_id] = random.choice(mappings)[0]

        mapping_one2one, mapping_one2many = translate_mappings(ids, translated_ids, mapping_one2one, mapping_one2many)
        if n_mapped < len(translated_ids):
            warnings.warn(f"Ortholog mapping found for only {n_mapped} out of {len(translated_ids)} gene IDs.")
            print({k: [this_v[0] for this_v in v] for k, v in mapping_one2many.items()})
        return OrthologDict(mapping_one2one), OrthologDict(
            {k: [this_v[0] for this_v in v] for k, v in mapping_one2many.items()})

    def get_paralogs(self, ids: Tuple[str, ...]):
        session = get_session(self.RETRIES)
        url = f'{self.API_URL}/services/oai/pantherdb/ortholog/homologOther'

        mapping_one2many = {}
        ids, translated_ids = self.translate_ids(ids, session)

        n_mapped = 0
        for from_id in tqdm(translated_ids, 'Mapping paralogs', unit='genes'):
            req_data = dict(geneInputList=from_id, organism=str(self.map_from_organism), homologType='P')
            req = session.post(url, params=req_data)
            req.raise_for_status()

            try:
                req_output = req.json()['search']['mapping']['mapped']
                for mapping in parsing.data_to_list(req_output):
                    if len(mapping) <= 1:
                        continue

                    to_id = parsing.parse_uniprot_id(mapping['target_gene'])
                    if to_id is None:
                        continue

                    if from_id not in mapping_one2many:
                        n_mapped += 1
                        mapping_one2many[from_id] = []
                    mapping_one2many[from_id].append(to_id)

            except KeyError:
                pass

        _, mapping_one2many = translate_mappings(ids, translated_ids, {}, mapping_one2many)

        if n_mapped < len(translated_ids):
            warnings.warn(f"Paralob mapping found for only {n_mapped} out of {len(translated_ids)} gene IDs.")

        return OrthologDict(mapping_one2many)


class GeneIDTranslator:
    UNIPROTKB_FROM = "UniProtKB_from"
    UNIPROTKB_TO = "UniProtKB_to"
    API_URL = "https://rest.uniprot.org"
    POLLING_INTERVAL = 3
    REQUEST_DELAY_MILLIS = 250
    REQ_MAX_ENTRIES = 10
    RETRIES = Retry(total=5, backoff_factor=0.25, status_forcelist=[500, 502, 503, 504])

    def __init__(self, map_from: str, map_to: str = 'UniProtKB AC', verbose: bool = True,
                 session: Union[requests.Session, None] = None):
        """
    :param map_from: identifier type to map from (for example 'UniProtKB AC' or 'WormBase')
    :type map_from: str
    :param map_to: identifier type to map to (for example 'UniProtKB AC' or 'WormBase'). \
    can be identical to 'map_from'
    :type map_to: str
    :param verbose: if False, verbose reports of mapping success/failure will be suppressed.
    :type verbose: bool (default=True)
    """
        self.verbose = verbose
        self.map_from = map_from
        self.map_to = map_to
        self.id_dicts = _get_id_abbreviation_dicts()

        id_dict_to, id_dict_from = self.id_dicts
        if id_dict_to[self.map_to] == id_dict_from[self.map_from]:
            self.to_from_identical = True
        else:
            self.to_from_identical = False

        validation.validate_uniprot_dataset_name(self.id_dicts, parsing.data_to_list(self.map_to),
                                                 parsing.data_to_list(self.map_from))
        if self.map_to in ["UniProtKB", "UniProtKB AC/ID"]:
            self.map_to = self.UNIPROTKB_TO
        if self.map_from in ["UniProtKB", "UniProtKB AC/ID"]:
            self.map_from = self.UNIPROTKB_FROM

        self.session = get_session(self.RETRIES) if session is None else session

    def run(self, ids: Union[str, Iterable[str]]) -> GeneIDDict:
        """
        Map gene IDs from one identifier type to another using the UniProt ID Mapping service. \
        If some IDs cannot be mapped uniquely, duplicate mappings will be resolved by their UniProtKB Annotation Score.\
         Gene IDs that could not be mapped or were not recognized will be dropped from the output.

        :param ids: gene IDs to be mapped
        :type ids: str or an Iterable of strings
        :return:a GeneIDTranslator object that uniquely maps each given gene ID in 'map_from' identifier type \
        to its matching gene ID in 'map_to' identifier type.
        :rtype: GeneIDDict
        """
        if len(ids) == 0:
            if self.verbose:
                warnings.warn('No IDs were given')
            return GeneIDDict({})

        # if self.map_from and self.map_to are the same, return an empty GeneIDTranslator (which will map any given gene ID to itself)
        if self.to_from_identical:
            return GeneIDDict()

        ids = parsing.data_to_list(ids)
        n_queries = len(ids)

        # since the Uniprot service can only translate to or from 'UniProtKB' identifier type,
        # if we need to map gene IDs between two other identifier types,
        # then we will map from 'self.map_from' to 'UniProtKB' and then from 'UniProtKB' to 'self.map_to'.
        if self.verbose:
            print(f"Mapping {n_queries} entries from '{self.map_from}' to '{self.map_to}'...")

        id_dict_to, id_dict_from = self.id_dicts
        if id_dict_to[self.map_to] != id_dict_to[self.UNIPROTKB_TO] and \
            id_dict_from[self.map_from] != id_dict_from[self.UNIPROTKB_FROM]:

            to_translator = GeneIDTranslator(self.map_from, self.UNIPROTKB_TO, self.verbose, self.session)
            to_uniprot = to_translator.run(ids).mapping_dict
            if to_uniprot is None:
                to_uniprot = {}

            from_translator = GeneIDTranslator(self.UNIPROTKB_FROM, self.map_to, self.verbose, self.session)
            from_uniprot = from_translator.run(to_uniprot.values()).mapping_dict
            if from_uniprot is None:
                from_uniprot = {}

            output_dict = {key: from_uniprot[val] for key, val in zip(to_uniprot.keys(), to_uniprot.values()) if
                           val in from_uniprot}

        # make sure that 'self.map_from' and 'self.map_to' are recognized identifier types
        elif id_dict_to[self.map_to] != 'Null' and id_dict_from[self.map_from] != 'Null':
            session = self.session
            results = self.get_mapping_results(self.map_to, self.map_from, ids, session)

            if results is None or len(results) <= 1:
                if self.verbose:
                    warnings.warn("No entries were mapped successfully.")
                return GeneIDDict({})

            output_dict, duplicates = self.format_annotations(results)
            self.handle_duplicates(output_dict, duplicates, session)
        else:
            output_dict = {}

        if len(output_dict) < n_queries and self.verbose:
            warnings.warn(
                f"Failed to map {n_queries - len(output_dict)} entries from '{self.map_from}' to '{self.map_to}'. "
                f"Returning {len(output_dict)} successfully-mapped entries.")

        self.reformat_ids(output_dict)
        return GeneIDDict(output_dict)

    def get_mapping_results(self, map_to: str, map_from: str, ids: Tuple[str, ...], session: requests.Session):
        id_dict_to, id_dict_from = self.id_dicts
        to_db = id_dict_to[map_to]
        from_db = id_dict_from[map_from]

        job_id = self.submit_id_mapping(to_db, from_db, session, ids)

        if self.check_id_mapping_results_ready(session, job_id, self.POLLING_INTERVAL):
            link = self.get_id_mapping_results_link(session, job_id)
            results = self.get_id_mapping_results_search(session, link)
            return results

    @staticmethod
    def submit_id_mapping(to_db: str, from_db: str, session: requests.Session, ids: Tuple[str, ...]):
        req = session.post(f"{GeneIDTranslator.API_URL}/idmapping/run",
                           data={"to": to_db, "from": from_db, "ids": ",".join(ids)})
        req.raise_for_status()
        return req.json()["jobId"]

    @staticmethod
    def get_next_link(headers: dict):
        re_next_link = re.compile(r'<(.+)>; rel="next"')
        if "Link" in headers:
            match = re_next_link.match(headers["Link"])
            if match:
                return match.group(1)

    @staticmethod
    def check_id_mapping_results_ready(session, job_id: str, polling_interval: float, verbose: bool = True):
        while True:
            r = session.get(f"{GeneIDTranslator.API_URL}/idmapping/status/{job_id}")
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

    @staticmethod
    def get_batch(session: requests.Session, batch_response: requests.Response, file_format: Literal['json', 'tsv']):
        batch_url = GeneIDTranslator.get_next_link(batch_response.headers)
        while batch_url:
            batch_response = session.get(batch_url)
            batch_response.raise_for_status()
            yield GeneIDTranslator.decode_results(batch_response, file_format)
            batch_url = GeneIDTranslator.get_next_link(batch_response.headers)

    @staticmethod
    def combine_batches(all_results, batch_results, file_format: Literal['json', 'tsv']):
        if file_format == "json":
            for key in ("results", "failedIds"):
                if batch_results[key]:
                    all_results[key] += batch_results[key]
        elif file_format == "tsv":
            return all_results + batch_results[1:]  # dump the table header line
        else:
            return all_results + batch_results
        return all_results

    @staticmethod
    def get_id_mapping_results_link(session: requests.Session, job_id):
        url = f"{GeneIDTranslator.API_URL}/idmapping/details/{job_id}"
        r = session.get(url)
        r.raise_for_status()
        return r.json()["redirectURL"]

    @staticmethod
    def decode_results(response: requests.Response, file_format: Literal['json', 'tsv']):
        if file_format == "json":
            return response.json()
        elif file_format == "tsv":
            return [line for line in response.text.split("\n") if line]
        return response.text

    def get_id_mapping_results_search(self, session: requests.Session, link):
        parsed = urlparse(link)
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
        results = self.decode_results(r, file_format)
        try:
            total = int(r.headers["x-total-results"])
        except KeyError:
            return ''
        if self.verbose:
            self.print_progress_batches(0, size, total)

        for i, batch in enumerate(self.get_batch(session, r, file_format)):
            results = self.combine_batches(results, batch, file_format)
            if self.verbose:
                self.print_progress_batches(i + 1, size, total)
        return results

    @staticmethod
    def print_progress_batches(batch_index: int, size: int, total: int):
        n = min((batch_index + 1) * size, total)
        print(f"Fetched: {n} / {total}")

    @staticmethod
    def format_annotations(results):
        df = pd.DataFrame([line.split('\t') for line in results[1:]], columns=results[0].split('\t'))
        # sort annotations by decreasing annotation score, so that the most relevant annotations are at the top
        if 'Annotation' in df.columns:
            df.loc[df['Annotation'] == '', 'Annotation'] = '0'
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
        return output_dict, duplicates

    def handle_duplicates(self, output_dict: dict, duplicates: dict, session: requests.Session):
        # handle duplicates
        if len(duplicates) > 0:
            if self.map_to == self.UNIPROTKB_TO:
                for match_from, match_to_options in duplicates.items():
                    output_dict[match_from] = match_to_options[0]
                duplicates_chosen = {match_from: match_to[0] for match_from, match_to in duplicates.items()}

            # if there are unproccessed duplicates, map them in reverse and sort then by annotation score
            else:
                ids_to_rev_map = parsing.flatten(parsing.data_to_list(duplicates.values()))

                rev_results = self.get_mapping_results(self.UNIPROTKB_TO, self.map_to, ids_to_rev_map, session)
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
            if self.verbose:
                warnings.warn(f"Duplicate mappings were found for {len(duplicates)} genes.  The following mapping "
                              f"was chosen for them based on their annotation score: {duplicates_chosen}")

    def reformat_ids(self, output_dict: dict):
        if self.map_to == 'Ensembl':
            for key, val in output_dict.items():
                output_dict[key] = re.sub('(\.\d+)$', '', val)


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
        key = [items[0], list(map_from_options).index(items[1]), list(map_to_options).index(items[2])]
        return key

    def map_gene_ids_ignore_httpexception(ids: Tuple[str], map_from: str, map_to: str):
        try:
            return GeneIDTranslator(map_from, map_to, verbose=False).run(ids), map_from, map_to
        except requests.exceptions.HTTPError:
            return GeneIDDict({}), map_from, map_to

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
def get_legal_panther_taxons():
    URL = 'http://www.pantherdb.org/services/oai/pantherdb/supportedgenomes'
    req = requests.post(URL)
    req.raise_for_status()
    entries = json.loads(req.text)['search']['output']['genomes']['genome']
    taxons = tuple(sorted(d['long_name'] for d in entries))
    return taxons


@functools.lru_cache(maxsize=2)
def get_legal_ensembl_taxons():
    endpoint = 'info/species'
    client = EnsemblRestClient()
    client.queue_action('get', endpoint, params=dict(hide_strain_info=1))
    results = client.run()[0]['species']
    return tuple(sorted(res['name'].replace('_', ' ').capitalize() for res in results))


@functools.lru_cache(maxsize=2)
def get_legal_phylomedb_taxons():
    entries = PhylomeDBOrthologMapper.get_legal_species()
    taxons = tuple(name for name in entries['name'].unique() if name[0].isupper())
    return taxons


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
def _get_id_abbreviation_dicts(
    dict_path: str = Path(__file__).parent.parent.joinpath('data_files/uniprot_dataset_abbreviation_dict.json')):
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
    video_dir_pth.mkdir(parents=True, exist_ok=True)

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

    except (aiohttp.ClientConnectorError, asyncio.TimeoutError):
        pass


def run_r_script(script_path: Union[str, Path], r_installation_folder: Union[str, Path, Literal['auto']] = 'auto'):
    if r_installation_folder == 'auto':
        prefix = "Rscript"
    else:
        prefix = f'{Path(r_installation_folder).as_posix()}/bin/Rscript'
    script_path = Path(script_path).as_posix()
    assert Path(script_path).exists() and Path(
        script_path).is_file(), f"Could not find the requested R script: {script_path}"

    try:
        return_code, _ = run_subprocess([prefix, "--help"], False, False)
        if return_code:
            raise FileNotFoundError

    except FileNotFoundError:
        raise FileNotFoundError("Failed to find R executable. "
                                "Please make sure your R installation folder is correct. \n"
                                "(For example: 'C:/Program Files/R/R-4.2.3')")

    return_code, stderr = run_subprocess([prefix, script_path])
    if return_code:
        full_err = 'See R log below. \n' + '\n'.join([s.rstrip() for s in stderr])
        short_err = stderr[-1].rstrip()
        for i, s in enumerate(stderr):
            if s.startswith('Error'):
                short_err = s.rstrip() + stderr[i + 1].rstrip()
                break
        raise ChildProcessError(f"R script failed to execute: '{short_err}'. See full error report below.") \
            from RuntimeError(full_err)


def stdout_reader(pipe, log_filename, lock, print_output: bool = True):
    with open(log_filename, 'a') if log_filename is not None else contextlib.nullcontext() as logfile:
        for line in (pipe if isinstance(pipe, list) else iter(pipe.readline, b'')):
            decoded_line = line.decode('utf8', errors="ignore")

            if print_output:
                print(decoded_line)

            if log_filename is not None:
                with lock:
                    logfile.write(decoded_line)
    if not isinstance(pipe, list):
        pipe.close()


def stderr_reader(pipe, stderr_record, log_filename, lock, print_output: bool = True):
    with open(log_filename, 'a') if log_filename is not None else contextlib.nullcontext() as logfile:
        for line in (pipe if isinstance(pipe, list) else iter(pipe.readline, b'')):
            decoded_line = line.decode('utf8', errors="ignore")
            stderr_record.append(decoded_line)

            if print_output:
                print(decoded_line)

            if log_filename is not None:
                with lock:
                    logfile.write(decoded_line)
    if not isinstance(pipe, list):
        pipe.close()


def run_subprocess(args: List[str], print_stdout: bool = True, print_stderr: bool = True,
                   log_filename: Union[str, None] = None, shell: bool = False) -> Tuple[int, List[str]]:
    # join List of args into a string of args when running in shell mode
    if shell:
        try:
            args = shlex.join(args)
        except AttributeError:
            args = ' '.join([shlex.quote(arg) for arg in args])

    stderr_record = []
    lock = threading.Lock()
    stdout = subprocess.PIPE
    stderr = subprocess.PIPE

    process = subprocess.Popen(args, stdout=stdout, stderr=stderr, shell=shell)

    stdout_thread = threading.Thread(target=stdout_reader, args=(process.stdout, log_filename, lock, print_stdout))
    stderr_thread = threading.Thread(target=stderr_reader,
                                     args=(process.stderr, stderr_record, log_filename, lock, print_stderr))

    stdout_thread.start()
    stderr_thread.start()
    stdout_thread.join()
    stderr_thread.join()
    process.wait()

    return process.returncode, stderr_record


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
        if shell:
            cmd = parsing.quote_path(installation_folder.joinpath(command))
        else:
            cmd = installation_folder.joinpath(command).as_posix()
        call = [cmd] + parsing.data_to_list(args)

    try:
        exit_code, stderr = run_subprocess(call + [version_command], shell=shell)
        assert exit_code == 0 or (len(stderr) >= 1 and 'Version' in stderr[0]), \
            f"call to {call[0]} exited with exit status {exit_code}: \n'{''.join(stderr)}'"
    except FileNotFoundError:
        raise FileNotFoundError(f"RNAlysis could not find '{command}'. "
                                'Please ensure that your installation folder is correct, or add it to PATH. ')
    except OSError:
        raise OSError(f"RNAlysis could not run '{call + [version_command]}'. \n"
                      f"Please ensure that the installed version of {command} "
                      f"matches your operating system and try again. ")

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
