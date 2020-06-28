from pathlib import Path
from typing import Union, List, Set, Iterable
import pandas as pd


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


def isinstanceinh(obj, parent_class):
    return True if issubclass(obj.__class__, parent_class) else False


def validate_biotype_table(ref_df: pd.DataFrame):
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


def validate_attr_table(ref_df: pd.DataFrame):
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


def validate_uniprot_dataset_name(dataset_dict: dict, *names: str):
    for name in names:
        assert name in dataset_dict, f"Dataset '{name}' is not a valid Uniprot Dataset for mapping gene names/IDs. " \
                                     f"Valid Uniprot Datasets are: {', '.join(dataset_dict.keys())}."
