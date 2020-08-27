"""
This module contains general-purpose functions. Those include saving Filter objects and result tables, \
reading and updating the settings file, parsing common types of genomic feature indices, etc.
"""

import re
from typing import Union

import pandas as pd

from rnalysis import __attr_file_key__, __biotype_file_key__
from rnalysis.utils import io, validation, ref_tables
from rnalysis.filtering import Filter


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
    return set(re.findall(r'WBGene[0-9]{8}', string))


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
    return set(re.findall(r'[A-Z,0-9]{5,8}\.\d{1,2}', string))


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
    return set(re.findall(r'[a-z]{3,4}-[A-Z,0-9.]{1,4}', string))


def reset_settings_file():
    """
    Resets the local settings by deleting the local settings file. Warning: this action is irreversible!
    """
    settings_pth = ref_tables.get_settings_file_path()
    if not settings_pth.exists():
        print("No local settings file exists. ")
    else:
        settings_pth.unlink()
        print(f"Local settings file was deleted. ")


def set_attr_ref_table_path(path: str = None):
    """
    Defines/updates the Attribute Reference Table path in the settings file.
    :param path: the path you wish to set as the Attribute Reference Table path
    :type path: str

    :Examples:
    >>> from rnalysis import general
    >>> path="my_attribute_reference_table_path"
    >>> general.set_attr_ref_table_path(path)
    Attribute Reference Table path set as: my_attribute_reference_table_path
    """
    if path is None:
        path = input("Please write the new Attribute Reference Table Path:\n")
    ref_tables.update_settings_file(path, __attr_file_key__)
    print(f'Attribute Reference Table path set as: {path}')


def set_biotype_ref_table_path(path: str = None):
    """
    Defines/updates the Biotype Reference Table path in the settings file.
    :param path: the path you wish to set as the Biotype Reference Table path
    :type path: str

    :Examples:
    >>> from rnalysis import general
    >>> path="my_biotype_reference_table_path"
    >>> general.set_biotype_ref_table_path(path)
    Biotype Reference Table path set as: my_biotype_reference_table_path
    """
    if path is None:
        path = input("Please write the new Attribute Reference Table Path:\n")
    ref_tables.update_settings_file(path, __biotype_file_key__)
    print(f'Biotype Reference Table path set as: {path}')


def print_settings_file():
    """
    Print the current setting file configuration.

    :Examples:
        >>> from rnalysis import general
        >>> general.print_settings_file()
        Attribute Reference Table used: my_attribute_reference_table_path
        Biotype Reference Table used: my_biotype_reference_table_path

    """
    ref_tables.get_attr_ref_path('predefined')
    ref_tables.get_biotype_ref_path('predefined')


def save_to_csv(df: Union[pd.DataFrame, Filter], filename: str):
    """
    save a pandas DataFrame or Filter object to csv.
    :type df: Filter object or pandas DataFrame
    :param df: object to be saved
    :type filename: str
    :param filename: name for the saved file. Specify full path to control the directory where the file will be saved.
    """
    if isinstance(df, pd.DataFrame):
        io.save_csv(df, filename)
    elif validation.isinstanceinh(df, Filter):
        io.save_csv(df.df, filename)
    else:
        raise TypeError(f"Object of type {type(df)} cannot be saved to csv")
