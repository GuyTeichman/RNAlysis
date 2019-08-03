"""
This module contains general-purpose functions, such as loading and saving files, removing table rows with no WBGene \
 index, identifying type of input variables, etc. This module is used mainly by other modules.
"""

import pandas as pd
from pathlib import Path
import os
import re
from rnalysis import __settings_start_phrase__


def parse_wbgene_string(string):
    """
    Receives a string that contains WBGene indices. Parses the string into a set of WBGene indices. \
    The format of a WBGene index is 'WBGene' and exactly 8 digits.
    :type string: str
    :param string: The string to be parsed. Can be any format of string.
    :return:
    a set of the WBGene indices that appear in the given string.
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
    """
    return set(re.findall('[A-Z,0-9]{5,6}\.\d', string))


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
    """
    return set(re.findall('[a-z]{3,4}-[A-Z,0-9]{1,4}', string))


def is_reference_table_defined(settings_start_phrase=__settings_start_phrase__):
    """
    Check whether a reference table path is defined.

    :param settings_start_phrase: The key phrase to look for in settings.ini ('big_table_path=').
    :return:
    True if defined, False if not defined.
    """
    settings_pth = Path(os.path.join(os.path.dirname(__file__), 'settings.ini'))

    if settings_pth.exists():
        with open(settings_pth) as f:
            setting = f.readline()
            if setting.startswith(settings_start_phrase):
                return True

    return False


def set_reference_table_path(path: str = None):
    """
    Set the path for a user-defined reference table for filtering and enrichment analyses. \
    The path will be saved in a settings.ini file that will persist through sessions.

    :type path: str
    :param path: The full path and name of the .csv reference file.
    """
    settings_pth = Path(os.path.join(os.path.dirname(__file__), 'settings.ini'))
    with open(settings_pth, 'w') as f:
        pth = path if path is not None else input(
            'Please insert the full path of your reference data table ("big table"):\n')
        f.write("big_table_path=" + pth)


def read_reference_table_path():
    """
    Attempt to read the reference table path from settings.ini. \
    If the path was not previously defined, will prompt user to define it.

    :return:
    The path of the reference table.
    """
    settings_pth = Path(os.path.join(os.path.dirname(__file__), 'settings.ini'))
    settings_start_phrase = 'big_table_path='
    if not is_reference_table_defined():
        set_reference_table_path()

    with open(settings_pth) as f:
        return f.readline()[len(settings_start_phrase)::]


def load_csv(filename: str, idx_col: int = None, drop_gene_names: bool = True):
    """
    loads a csv df into a pandas dataframe.

    :param filename: name of the csv file to be loaded
    :param idx_col: number of column to be used as index. default is None, meaning no column will be used as index.
    :param drop_gene_names: bool. If True, tries to drop the 'genes' column from the DataFrame (useful for filtering).
    :return: a pandas dataframe of the csv file
    """
    assert isinstance(filename, (str, Path))
    # with open(filename) as f:
    encoding = 'ISO-8859-1'
    if idx_col is not None:
        df = pd.read_csv(filename, index_col=idx_col, encoding=encoding)
    else:
        df = pd.read_csv(filename, encoding=encoding)
    if drop_gene_names:
        if 'genes' in df:
            df.drop('genes', axis=1, inplace=True)
    return df


def norm_reads_to_rpm(htcount: pd.DataFrame, features: pd.DataFrame):
    """
    Receives a dataframe of reads/counts and a dataframe of feature counts (ambiguous, no feature, not aligned, etc). \
    Divides each column in the read dataframe by (total reads + ambiguous + no feature)*10^-6 in order to normalize the \
    reads to RPM (reads per million).

    :param htcount: DataFrame of reads/counts
    :param features: the corresponding DataFrame of
    :return:
    a DataFrame normalized to RPM
    """
    for column in htcount.columns:
        norm_factor = (htcount[column].sum() + features.loc[r'__ambiguous', column] + features.loc[
            r'__no_feature', column]) / (10 ** 6)
        htcount[column] /= norm_factor
    return htcount


def remove_unindexed_rows(df: pd.DataFrame):
    """
    removes rows which have no WBGene index.

    :param df: a DataFrame with WBGene indices
    :return: a new DataFrame in which all rows with no WBGene index are removed
    """
    return df[[not i for i in df.index.isna()]]


def check_is_df(inp):
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


def save_to_csv(df: pd.DataFrame, filename: str, suffix: str):
    """
    save a pandas DataFrame to csv.

    :param df: pandas DataFrame to be saved
    :param filename: a string or pathlib.Path object stating the original name of the file
    :param suffix: str, a suffix to be added to the original name of the file
    """
    fname = Path(filename)
    new_fname = Path(f"{str(fname.parent)}\\{fname.stem}{suffix}{fname.suffix}")
    df.to_csv(new_fname)


def filter_low_rpm(df: pd.DataFrame, threshold: float = 5):
    """
    remove all features which have less then 'threshold' reads per million in all conditions.

    :param df: pandas DataFrame to be filtered
    :param threshold: float. The minimal rpm a feature needs to have in at least one sample in order not to be filtered out.
    :return:
    A filtered DataFrame
    """
    return df.loc[[True if max(vals) > threshold else False for gene, vals in df.iterrows()]]


def normalize_file(htcount_filename: str, feature_filename: str):
    """
    receives an htcount filename and a feature filename, normalizes the htcount file to rpm, removes unindexed rows, \
    and saves the file with the suffix '_rpm'.

    :param htcount_filename: str or pathlib.Path object. Name of the htcount file to be normalized
    :param feature_filename: str or pathlib.path object. Name of the corresponding feature file to be used for normalizing
    """
    save_to_csv(remove_unindexed_rows(
        norm_reads_to_rpm(load_csv(htcount_filename, idx_col=0), load_csv(feature_filename, idx_col=0))), '_rpm')
