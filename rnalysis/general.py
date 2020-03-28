"""
This module contains general-purpose functions, such as loading and saving files, \
reading and updating the settings file, identifying type of input variables, etc. \
This module is used mainly by other modules.
"""

import pandas as pd
from pathlib import Path
import os
import re
import time
import subprocess
import yaml
from typing import Union, List, Set, Dict, Tuple
from rnalysis import __attr_file_key__, __biotype_file_key__


def _start_ipcluster(n_engines: int = 'default'):
    """
    Start an ipyparallel ipcluster in order to perform parallelized computation.

    :type n_engines: int or 'default'
    :param n_engines: if 'default', will initiate the default amount of engines. \
    Otherwise, will initiate n_engines engines.
    """
    if n_engines == 'default':
        return subprocess.Popen("ipcluster start", stderr=subprocess.PIPE)
    else:
        return subprocess.Popen(["ipcluster", "start", "-n={:d}".format(n_engines)], stderr=subprocess.PIPE)


def _stop_ipcluster():
    """
    Stop a previously started ipyparallel ipcluster.

    """
    subprocess.Popen("ipcluster stop", stderr=subprocess.PIPE)


def start_parallel_session(n_engines: int = 'default'):
    """
    Stop previous ipyparallel ipcluster and start a new one in order to perform parallelized computation.

    :type n_engines: int or 'default'
    :param n_engines: if 'default', will initiate the default amount of engines. \
    Otherwise, will initiate n_engines engines.

    :Examples:
    >>> from rnalysis import general
    >>> general.start_parallel_session()
    Parallel session started successfully
    """
    _stop_ipcluster()
    time.sleep(1)
    stream = _start_ipcluster(n_engines)
    time.sleep(1)
    while True:
        line = stream.stderr.readline()
        if 'Engines appear to have started successfully' in str(line):
            break
    print('Parallel session started successfully')


def parse_wbgene_string(string):
    """
    Receives a string that contains WBGene indices. Parses the string into a set of WBGene indices. \
    The format of a WBGene index is 'WBGene' and exactly 8 digits.
    :type string: str
    :param string: The string to be parsed. Can be any format of string.
    :return:
    a set of the WBGene indices that appear in the given string.

    :Examples:
    >>> from rnalysis import general
    >>> string =  '''WBGene WBGenes WBGene12345678, WBGene98765432WBGene00000000& the geneWBGene44444444daf-16A5gHB.5
    ... WBGene55555555'''
    >>> parsed = general.parse_wbgene_string(string)
    >>> print(parsed)
    {'WBGene12345678', 'WBGene44444444', 'WBGene98765432', 'WBGene55555555', 'WBGene00000000'}
    """
    return set(re.findall('WBGene[0-9]{8}', string))


def parse_sequence_name_string(string):
    """
    Receives a string that contains sequence names (such as 'Y55D5A.5'). \
    Parses the string into a set of WBGene indices. \
    The format of a sequence name is a sequence consisting of the expression '[A-Z,0-9]{5,6}', \
    the character '.', and a digit.
    :type string: str
    :param string: The string to be parsed. Can be any format of string.
    :return:
    a set of the WBGene indices that appear in the given string.

    :Examples:
    >>> from rnalysis import general
    >>> string = 'CELE_Y55D5A.5T23G5.6WBGene00000000 daf-16^^ZK662.4 '
    >>> parsed = general.parse_sequence_name_string(string)
    >>> print(parsed)
    {'Y55D5A.5', 'T23G5.6', 'ZK662.4'}
    """
    return set(re.findall('[A-Z,0-9]{5,8}\.\d{1,2}', string))


def parse_gene_name_string(string):
    """
    Receives a string that contains gene names (like 'daf-2' or 'lin15B'). \
    Parses the string into a set of gene names. \
    The format of a gene name is a sequence consisting of the expression \
    '[a-z]{3,4}', the character '-', and the expression '[A-Z,0-9]{1,4}'.
    :type string: str
    :param string: The string to be parsed. Can be any format of string.
    :return:
    a set of the WBGene indices that appear in the given string.

    :Examples:
    >>> from rnalysis import general
    >>> string = 'saeg-2 lin-15B cyp-23A1lin-15A WBGene12345678%GHF5H.3'
    >>> parsed = general.parse_gene_name_string(string)
    >>> print(parsed)
    {'saeg-2', 'lin-15B', 'cyp-23A1', 'lin-15A'}
    """
    return set(re.findall('[a-z]{3,4}-[A-Z,0-9,.]{1,4}', string))


def _get_settings_file_path():
    """
    Generates the full path of the settings.yaml file.
    :returns: the path of the settings.yaml file.
    :rtype: pathlib.Path
    """
    return Path(os.path.join(os.path.dirname(__file__), 'settings.yaml'))


def _load_settings_file():
    """
    loads and parses the settings.yaml file into a dictionary.
    :rtype: dict
    """
    settings_pth = _get_settings_file_path()
    if not settings_pth.exists():
        return dict()
    with open(settings_pth) as f:
        settings = yaml.safe_load(f)
        if settings is None:
            settings = dict()
        return settings


def _update_settings_file(value: str, key: str):
    """
    Receives a key and a value, and updates/adds the key and value to the settings.yaml file.
    :param value: the value to be added/updated (such as Reference Table path)
    :type value: str
    :param key: the key to be added/updated (such as __attr_file_key__)
    :type key: str
    """
    settings_pth = _get_settings_file_path()
    out = _load_settings_file()
    out[key] = value
    with open(settings_pth, 'w') as f:
        yaml.safe_dump(out, f)


def _read_value_from_settings(key):
    """
    Attempt to read the value corresponding to a given key from the settings.yaml file. \
    If the key was not previously defined, the user will be prompted to define it.

    :type key: str
    :param key: the key in the settings file whose value to read.

    :return:
    The path of the reference table.
    """
    settings = _load_settings_file()
    if key not in settings:
        _update_settings_file(input(f'Please insert the full path of {key}:\n'), key)
        settings = _load_settings_file()
    return settings[key]


def read_biotype_ref_table_path():
    """
    Reads the Biotype Reference Table path from the settings file.

    :returns: the path of the Biotype Reference Table that is saved in the settings file.
    :rtype: str

    :Examples:
    >>> from rnalysis import general
    >>> general.read_biotype_ref_table_path()
    Biotype Reference Table used: the_biotype_reference_table_path_that_was_saved_in_the_settings_file
    """
    pth = _read_value_from_settings(__biotype_file_key__)
    print(f'Biotype Reference Table used: {pth}')
    return pth


def read_attr_ref_table_path():
    """
    Reads the Attribute Reference Table path from the settings file.

    :returns: the path of the Attribute Reference Table that is saved in the settings file.
    :rtype: str

    :Examples:
    >>> from rnalysis import general
    >>> path = general.read_attr_ref_table_path()
    Attribute Reference Table used: the_attribute_reference_table_path_that_was_saved_in_the_settings_file
    """
    pth = _read_value_from_settings(__attr_file_key__)
    print(f'Attribute Reference Table used: {pth}')
    return pth


def set_attr_ref_table_path(path: str = None):
    """
    Defines/updates the Attribute Reference Table path in the settings file.
    :param path: the path you wish to set as the Attribute Reference Table path
    :type path: str

    :Examples:
    >>> from rnalysis import general
    >>> path="the_new_attribute_reference_table_path"
    >>> general.set_attr_ref_table_path(path)
    Attribute Reference Table path set as: the_new_attribute_reference_table_path
    """
    if path is None:
        path = input("Please write the new Attribute Reference Table Path:\n")
    _update_settings_file(path, __attr_file_key__)
    print(f'Attribute Reference Table path set as: {path}')


def set_biotype_ref_table_path(path: str = None):
    """
    Defines/updates the Biotype Reference Table path in the settings file.
    :param path: the path you wish to set as the Biotype Reference Table path
    :type path: str

    :Examples:
    >>> from rnalysis import general
    >>> path="the_new_biotype_reference_table_path"
    >>> general.set_biotype_ref_table_path(path)
    Attribute Reference Table path set as: the_new_biotype_reference_table_path
    """
    if path is None:
        path = input("Please write the new Attribute Reference Table Path:\n")
    _update_settings_file(path, __biotype_file_key__)
    print(f'Biotype Reference Table path set as: {path}')


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
            for i in drop_columns:
                assert isinstance(i, str), f"'drop_columns' must contain strings only. Member {i} is of type {type(i)}."
                if i in df:
                    df.drop('genes', axis=1, inplace=True)
                else:
                    raise IndexError(f"The argument {i} in 'drop_columns' is not a column in the loaded csv file!")
    return df


def _remove_unindexed_rows(df: pd.DataFrame):
    """
    removes rows which have no WBGene index.

    :param df: a DataFrame with WBGene indices
    :return: a new DataFrame in which all rows with no WBGene index are removed
    """
    return df[[not i for i in df.index.isna()]]


def _check_is_df(inp):
    """
    checks whether an input file is a pandas DataFrame, a string that represent a path of a .csv file, a Path object \
    of a .csv file, or an invalid input.

    :param inp: the input we wish to test
    :return: True if pandas DataFrame, False if a string/Path object that leads to a .csv file.\
    Raises ValueError otherwise.
    """
    if isinstance(inp, pd.DataFrame):
        return True
    elif isinstance(inp, str):
        if inp[-4::] == '.csv':
            return False
    elif isinstance(inp, Path):
        if inp.suffix == '.csv':
            return False
    raise ValueError("The input is neither a pandas DataFrame or a csv file")


def save_to_csv(df: pd.DataFrame, filename: str, suffix: str = None, index: bool = True):
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
    new_fname = Path(f"{str(fname.parent)}\\{fname.stem}{suffix}{fname.suffix}")
    df.to_csv(new_fname, header=True)


def _get_biotype_ref_path(ref: Union[str, Path]):
    """
    Returns the predefined Biotype Reference Table path from the settings file if ref='predefined', \
    otherwise returns 'ref' unchanged.

    :param ref: the 'ref' argument from a filtering module/enrichment module function
    :type ref: str, pathlib.Path or 'predefined'
    :returns: if ref is 'predefined', returns the predefined Biotype Reference Table path from the same file. \
    Otherwise, returns 'ref'.
    :rtype: str
    """
    if ref == 'predefined':
        return read_biotype_ref_table_path()
    else:
        return ref


def _get_attr_ref_path(ref):
    """
    Returns the predefined Attribute Reference Table path from the settings file if ref='predefined', \
    otherwise returns 'ref' unchanged.

    :param ref: the 'ref' argument from a filtering module/enrichment module function
    :type ref: str, pathlib.Path or 'predefined'
    :returns: if ref is 'predefined', returns the predefined Attribute Reference Table path from the same file. \
    Otherwise, returns 'ref'.
    :rtype: str
    """
    if ref == 'predefined':
        return read_attr_ref_table_path()
    else:
        return ref
