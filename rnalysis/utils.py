"""
This module contains various utility functions. \
This module is used mainly by other modules, and is meant for internal use only.
"""

import os
import subprocess
from pathlib import Path
from typing import List, Union
from sklearn.preprocessing import PowerTransformer, StandardScaler

import pandas as pd
import numpy as np
import yaml

from rnalysis import __attr_file_key__, __biotype_file_key__, __path__


def start_ipcluster(n_engines: int = 'default'):
    """
    Start an ipyparallel ipcluster in order to perform parallelized computation.

    :type n_engines: int or 'default'
    :param n_engines: if 'default', will initiate the default amount of engines. \
    Otherwise, will initiate n_engines engines.
    """
    assert (isinstance(n_engines,
                       int) and n_engines > 0) or n_engines == 'default', f"Invalid number of engines {n_engines}"
    if n_engines == 'default':
        return subprocess.Popen("ipcluster start", stderr=subprocess.PIPE, shell=True)
    else:
        return subprocess.Popen(["ipcluster", "start", "-n={:d}".format(n_engines)], stderr=subprocess.PIPE, shell=True)


def stop_ipcluster():
    """
    Stop a previously started ipyparallel ipcluster.

    """
    subprocess.Popen("ipcluster stop", stderr=subprocess.PIPE, shell=True)


def get_settings_file_path():
    """
    Generates the full path of the settings.yaml file.
    :returns: the path of the settings.yaml file.
    :rtype: pathlib.Path
    """
    # return Path(os.path.join(os.path.dirname(__file__), 'settings.yaml'))
    return Path(os.path.join(__path__[0], 'settings.yaml'))


def load_settings_file():
    """
    loads and parses the settings.yaml file into a dictionary.
    :rtype: dict
    """
    settings_pth = get_settings_file_path()
    if not settings_pth.exists():
        return dict()
    with settings_pth.open() as f:
        settings = yaml.safe_load(f)
        if settings is None:
            settings = dict()
        return settings


def update_settings_file(value: str, key: str):
    """
    Receives a key and a value, and updates/adds the key and value to the settings.yaml file.

    :param value: the value to be added/updated (such as Reference Table path)
    :type value: str
    :param key: the key to be added/updated (such as __attr_file_key__)
    :type key: str
    """
    settings_pth = get_settings_file_path()
    out = load_settings_file()
    out[key] = value
    with settings_pth.open('w') as f:
        yaml.safe_dump(out, f)


def read_value_from_settings(key):
    """
    Attempt to read the value corresponding to a given key from the settings.yaml file. \
    If the key was not previously defined, the user will be prompted to define it.

    :type key: str
    :param key: the key in the settings file whose value to read.

    :return:
    The path of the reference table.
    """
    settings = load_settings_file()
    if key not in settings:
        update_settings_file(input(f'Please insert the full path of {key}:\n'), key)
        settings = load_settings_file()
    return settings[key]


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


def remove_unindexed_rows(df: pd.DataFrame):
    """
    removes rows which have no WBGene index.

    :param df: a DataFrame with WBGene indices
    :return: a new DataFrame in which all rows with no WBGene index are removed
    """
    return df[[not i for i in df.index.isna()]]


def check_is_df_like(inp):
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


def get_biotype_ref_path(ref: Union[str, Path]):
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
        pth = read_value_from_settings(__biotype_file_key__)
        print(f'Biotype Reference Table used: {pth}')

        return pth
    else:
        print(f'Biotype Reference Table used: {ref}')
        return ref


def get_attr_ref_path(ref):
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
        pth = read_value_from_settings(__attr_file_key__)
        print(f'Attribute Reference Table used: {pth}')
        return pth
    else:
        print(f'Attribute Reference Table used: {ref}')
        return ref


def biotype_table_assertions(ref_df: pd.DataFrame):
    """
    Assert legality of Biotype Reference Table, and rename column names to standard names ('gene' and 'biotype').
    :param ref_df: the loaded Biotype Reference Table
    :type ref_df: pandas DataFrame

    """
    assert ref_df.shape[
               1] == 2, f"Invalid number of columns in Biotype Reference Table: found {ref_df.shape[1]} columns instead of 2!"
    assert ref_df.shape[
               0] >= 2, f"Biotype Reference Table must have at least two rows, found only  {ref_df.shape[0]}!"
    ref_df.rename(columns={ref_df.columns[0]: 'gene', ref_df.columns[1]: 'biotype'}, inplace=True)


def attr_table_assertions(ref_df: pd.DataFrame):
    """
    Assert legality of Attribute Reference Table, and renames the first column to standard name ('gene').
    :param ref_df:
    :type ref_df: pandas DataFrame

    """
    assert ref_df.shape[
               1] >= 2, f"Attribute Reference Table must have at least two columns, found only  {ref_df.shape[1]}!"
    assert ref_df.shape[
               0] >= 2, f"Attribute Reference Table must have at least two rows, found only  {ref_df.shape[0]}!"
    ref_df.rename(columns={ref_df.columns[0]: 'gene'}, inplace=True)


def from_string(msg: str = '', del_spaces: bool = False, delimiter: str = '\n'):
    """
    Takes a manual string input from the user, and then splits it using a comma delimiter into a list of values. \
    Called when an FeatureSet instance is created without input, \
    or when FeatureSet.enrich_randomization is called without input.

    :param msg: a promprt to be printed to the user
    :param del_spaces: if True, will delete all spaces in each delimited value.
    :param delimiter: the delimiter used to separate the values. Default is '\n'
    :return: A list of the comma-seperated values the user inserted.

    """
    string = input(msg)
    split = string.split(sep=delimiter)
    if del_spaces:
        for i in range(len(split)):
            split[i] = split[i].replace(' ', '')
    if split[-1] == '':
        split = split[:-1]
    return split


def isinstanceinh(obj, parent_class):
    return True if issubclass(obj.__class__, parent_class) else False


def make_temp_copy_of_settings_file():
    pth = get_settings_file_path()
    try:
        remove_temp_copy_of_settings_file()
    except FileNotFoundError:
        print("no previous temporary test file existed")
    if not pth.exists():
        print("no previous settings file exists")
        return
    with open(os.path.join(str(pth.parent), 'temp_settings.yaml'), 'w') as tempfile, pth.open() as originfile:
        tempfile.writelines(originfile.readlines())


def set_temp_copy_of_settings_file_as_default():
    pth = get_settings_file_path()
    if pth.exists():
        pth.unlink()
    if not Path(os.path.join(str(pth.parent), 'temp_settings.yaml')).exists():
        print("no temporary settings file exists")
        return
    with open(os.path.join(str(pth.parent), 'temp_settings.yaml')) as tempfile, pth.open('w') as originfile:
        originfile.writelines(tempfile.readlines())


def remove_temp_copy_of_settings_file():
    pth = get_settings_file_path()
    tmp_pth = Path(os.path.join(str(pth.parent), 'temp_settings.yaml'))
    if tmp_pth.exists():
        tmp_pth.unlink()


def standard_box_cox(data: np.ndarray):
    """

    :param data:
    :type data:
    :return:
    :rtype:
    """
    return StandardScaler().fit_transform(PowerTransformer(method='box-cox').fit_transform(data + 1))


def standardize(data: np.ndarray):
    """

    :param data:
    :type data:
    :return:
    :rtype:
    """
    return StandardScaler().fit_transform(data)
