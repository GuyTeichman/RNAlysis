import itertools
import re
from tqdm.auto import tqdm
import warnings
from itertools import islice
from typing import Any, Dict, Union, List, Tuple

import numpy as np
import pandas as pd


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


def uniprot_tab_to_dict(tab_input: str) -> Tuple[Dict[str, str], list]:
    split_list = tab_input.split()
    if len(split_list) == 0:
        return {}, []
    split_list = split_list[2:]

    parsed = {}
    duplicates = []
    for key, val in zip(islice(split_list, 0, len(split_list), 2), islice(split_list, 1, len(split_list), 2)):
        if key in parsed:
            parsed[key].append(val)
        else:
            parsed[key] = [val]

    for key in list(parsed.keys()):
        if len(parsed[key]) > 1:
            duplicates.extend(parsed.pop(key))
        else:
            parsed[key] = parsed[key][0]

    return parsed, duplicates


def uniprot_tab_with_score_to_dict(tab_input: str, reverse_key_value: bool = False) -> Dict[str, str]:
    split_list = re.split('[\t\n]', tab_input)
    if len(split_list) == 0:
        return {}
    split_list = split_list[3:]

    parsed: Dict[str, Tuple[List[str], List[int]]] = {}
    for gene_id, rank, key in zip(islice(split_list, 0, len(split_list), 3), islice(split_list, 1, len(split_list), 3),
                                  islice(split_list, 2, len(split_list), 3)):
        if reverse_key_value:
            gene_id, key = key, gene_id
        numeric_rank = int(rank[0])
        if key in parsed:
            parsed[key][0].append(gene_id)
            parsed[key][1].append(numeric_rank)
        else:
            parsed[key] = ([gene_id], [numeric_rank])

    key_to_id = {}
    duplicates = {}
    for key, (gene_ids, ranks) in zip(parsed.keys(), parsed.values()):
        if len(gene_ids) == 1:
            key_to_id[key] = gene_ids[0]
        else:
            best_id = gene_ids[max(range(len(ranks)), key=lambda i: ranks[i])]
            key_to_id[key] = best_id
            duplicates[key] = best_id

    if len(duplicates) > 0:
        warnings.warn(f"Duplicate mappings were found for {len(duplicates)} genes. "
                      f"The following mapping was chosen for them based on their annotation score: {duplicates}")
    return key_to_id


def data_to_list(data: Any, sort: bool = False) -> list:
    if isinstance(data, list):
        lst = data
    elif isinstance(data, (set, tuple, np.ndarray)):
        lst = list(data)
    elif isinstance(data, (int, float, bool, str)):
        lst = [data]
    elif data is None:
        lst = [None]

    elif callable(data):
        lst = [data]
    else:
        try:
            lst = list(data)
        except TypeError:
            raise TypeError(f"Invalid type {type(data)}.")

    if sort:
        lst.sort()
    return lst


def data_to_tuple(data: Any, sort: bool = False) -> tuple:
    if isinstance(data, tuple):
        tpl = data
    elif isinstance(data, (set, list, np.ndarray)):
        tpl = tuple(data)
    elif isinstance(data, (int, float, bool, str)):
        tpl = data,
    elif data is None:
        tpl = None,
    elif callable(data):
        tpl = data,
    else:
        try:
            tpl = tuple(data)
        except TypeError:
            raise TypeError(f"Invalid type {type(data)}.")
    if sort:
        tpl = tuple(sorted(tpl))
    return tpl


def data_to_set(data: Any) -> set:
    if isinstance(data, set):
        return data
    elif isinstance(data, (list, tuple, np.ndarray)):
        return set(data)
    elif isinstance(data, (int, float, bool, str)):
        return {data}
    elif data is None:
        return {None}
    elif callable(data):
        return {data}
    else:
        try:
            return set(data)
        except TypeError:
            raise TypeError(f"Invalid type {type(data)}.")


def sparse_dict_to_bool_df(sparse_dict: Dict[str, set], progress_bar_desc: str = '') -> pd.DataFrame:
    fmt = '{desc}: {percentage:3.0f}%|{bar}| [{elapsed}<{remaining}]'

    rows = list(sparse_dict.keys())
    columns_set = set()

    for val in sparse_dict.values():
        columns_set.update(val)
    columns = data_to_list(columns_set)
    df = pd.DataFrame(np.zeros((len(rows), len(columns)), dtype=bool), columns=columns, index=rows)
    for row, col in zip(tqdm(sparse_dict.keys(), desc=progress_bar_desc, bar_format=fmt), sparse_dict.values()):
        df.loc[data_to_list(row), data_to_list(col)] = True
    return df


def partition_list(lst: Union[list, tuple], chunk_size: int) -> Union[List[tuple], List[list]]:
    assert isinstance(lst, (list, tuple)), f"'lst' must be a list or tuple; instead got type {type(lst)}."
    assert isinstance(chunk_size, int), f"'chunk_size' must be an integer; instead got type {type(chunk_size)}."
    assert chunk_size > 0, f"'chunk_size' must be >0; instead got {chunk_size}."
    if len(lst) == 0:
        return [type(lst)()]
    return [lst[i: i + chunk_size] for i in range(0, len(lst), chunk_size)]


def flatten(lst: list) -> list:
    """
    Flatten a list of arbitrary depth.

    :param lst: the list to be flattened
    :type lst: list
    :return: a flattened list
    :rtype: list
    """
    output = []
    for item in lst:
        if isinstance(item, list):
            output.extend(flatten(item))
        else:
            output.append(item)
    return output


def parse_docstring(docstring: str) -> Tuple[str, Dict[str, str]]:
    """
    Parse a given docstring (str) to retreive the description text, as well as a dictionary of parameter descriptions.

    :param docstring: the docstring to be parsed
    :type docstring: str
    :return: a string matching the description, and a dictionary containing the parameter descriptions
    """
    docstring = re.sub(' +', ' ', docstring)
    split = docstring.split('\n\n')
    desc = split[0]
    params_str = split[1]
    free_text_match = '[\w\s\.\(\)\-\:%\*=@!\?\+\_,/' + "'" + ']'
    params_matches = list(
        re.finditer('^:param ([a-zA-Z_0-9]+):(' + free_text_match + '+?)(?=^\:.*\:)', params_str, re.MULTILINE))
    params = {match.group(1): match.group(2).replace('. ', '. \n') for match in params_matches}
    return desc, params


def generate_upset_series(objs: dict):
    """
    Receives a dictionary of sets from enrichment._fetch_sets(), \
    and reformats it as a pandas Series to be used by the python package 'upsetplot'.

    :param objs: the output of the enrichment._fetch_sets() function.
    :type objs: dict of sets
    :return: a pandas Series in the format requested by the 'upsetplot' package.

    """
    names = list(objs.keys())
    multi_ind = pd.MultiIndex.from_product([[True, False] for _ in range(len(names))], names=names)[:-1]
    srs = pd.Series(index=multi_ind, dtype='uint32')
    for ind in multi_ind:
        intersection_sets = list(itertools.compress(names, ind))
        difference_sets = list(itertools.compress(names, (not i for i in ind)))
        group = set.intersection(*[objs[s] for s in intersection_sets]).difference(*[objs[s] for s in difference_sets])
        group_size = len(group)
        srs.loc[ind] = group_size
    return srs
