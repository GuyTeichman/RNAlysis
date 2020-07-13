import warnings
from itertools import islice
import numpy as np
from typing import Union, Any, Set, Iterable
import pandas as pd
from rnalysis.utils import  validation


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


def uniprot_tab_to_dict(tab_input: str) -> dict:
    split_list = tab_input.split()
    if len(split_list) == 0:
        return {}
    split_list = split_list[2:]

    parsed = {}
    duplicates = []
    for key, val in zip(islice(split_list, 0, len(split_list), 2), islice(split_list, 1, len(split_list), 2)):
        if key in parsed:
            duplicates.append((key, val))
        else:
            parsed[key] = val

    if len(duplicates) > 0:
        warnings.warn(f"{len(duplicates)} duplicate mappings were found and ignored: {duplicates}")
    return parsed


def data_to_list(data: Any) -> list:
    if isinstance(data, list):
        return data
    elif isinstance(data, (set, tuple, np.ndarray)):
        return list(data)
    elif isinstance(data, (int, float, bool, str)):
        return [data]
    else:
        try:
            return list(data)
        except TypeError:
            raise TypeError(f"Invalid type {type(data)}.")


def data_to_tuple(data: Any) -> tuple:
    if isinstance(data, tuple):
        return data
    elif isinstance(data, (set, list, np.ndarray)):
        return tuple(data)
    elif isinstance(data, (int, float, bool, str)):
        return (data,)
    else:
        try:
            return tuple(data)
        except TypeError:
            raise TypeError(f"Invalid type {type(data)}.")


def data_to_set(data: Any) -> set:
    if isinstance(data, set):
        return data
    elif isinstance(data, (list, tuple, np.ndarray)):
        return set(data)
    elif isinstance(data, (int, float, bool, str)):
        return {data}
    else:
        try:
            return set(data)
        except TypeError:
            raise TypeError(f"Invalid type {type(data)}.")


def sparse_dict_to_bool_df(sparse_dict: dict) -> pd.DataFrame:
    rows = list(sparse_dict.keys())
    columns = set()
    for val in sparse_dict.values():
        for item in val:
            columns.add(item)
    df = pd.DataFrame(np.zeros((len(rows), len(columns)), dtype=bool), columns=columns, index=rows)
    for key in sparse_dict:
        for item in sparse_dict[key]:
            df.loc[key, item] = True
    return df


def parse_evidence_types(evidence_types: Union[str, Iterable[str]], evidence_type_dict: dict) -> Set[str]:
    if evidence_types == 'any':
        return set.union(*[set(s) for s in evidence_type_dict.values()])
    elif isinstance(evidence_types, str) and evidence_types.lower() in evidence_type_dict:
        return evidence_type_dict[evidence_types.lower()]
    elif validation.isiterable(evidence_types) and \
        any([isinstance(ev_type, str) and ev_type.lower() in evidence_type_dict for ev_type in evidence_types]):
        return set.union(
            *[evidence_type_dict[ev_type.lower()] if ev_type.lower() in evidence_type_dict else {ev_type} for ev_type in
              evidence_types])
    elif evidence_types is None:
        return set()
    else:
        return data_to_set(evidence_types)


def parse_go_aspects(aspects: Union[str, Iterable[str]], aspects_dict: dict) -> Set[str]:
    if aspects == 'any':
        return set.union(*[set(s) for s in aspects_dict.values()])
    elif any([isinstance(aspect, str) and aspect.lower() in aspects_dict for aspect in aspects]):
        return {aspects_dict[aspect.lower()] if aspect.lower() in aspects_dict else aspect for aspect in aspects}
    else:
        return data_to_set(aspects)
