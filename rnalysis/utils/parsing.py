import queue
import re
import warnings
from itertools import islice
from typing import Any, Dict, Iterable, Set, Union, List, Tuple

import numpy as np
import pandas as pd

from rnalysis.utils import validation


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

    for key in parsed:
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
        columns.update(val)
    df = pd.DataFrame(np.zeros((len(rows), len(columns)), dtype=bool), columns=columns, index=rows)
    for key in sparse_dict:
        df.loc[key, sparse_dict[key]] = True
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


class GOTerm:
    __slots__ = {'id': 'GO ID', 'name': 'GO Term name', 'level': "GO Term's level in the DAG Tree",
                 'relationships': 'direct parent relationships of the GO Term',
                 'children_relationships': 'direct children relationships of the GO Term'}

    def __init__(self):
        self.id = None
        self.name = None
        self.level = None
        self.relationships: Dict[str, List[str]] = {'is_a': [], 'part_of': []}
        self.children_relationships: Dict[str, List[str]] = {'is_a': [], 'part_of': []}

    def get_parents(self, relationships: Union[str, list] = ('is_a', 'part_of')) -> List[str]:
        go_ids = []
        _ = [go_ids.extend(self.relationships[relationship]) for relationship in data_to_list(relationships) if
             relationship in self.relationships]
        return go_ids

    def get_children(self, relationships: Union[str, list] = ('is_a', 'part_of')) -> List[str]:
        go_ids = []
        _ = [go_ids.extend(self.children_relationships[relationship]) for relationship in data_to_list(relationships) if
             relationship in self.children_relationships]
        return go_ids


def parse_go_id(byte_sequence: bytes) -> str:
    return re.findall(b"GO:[0-9]{7}", byte_sequence)[0].decode('utf8')


class DAGTreeParser:
    __slots__ = {'data_version': 'version of the go-basic.obo file',
                 'go_terms': 'dictionary of GO Terms in the DAG Tree',
                 'alt_ids': 'mapping of alternagive GO IDs to their main GO ID',
                 'levels': 'list of levels in the DAG Tree',
                 'parent_relationship_types': 'the types of relationships that constitute parenthood in the DAG Tree'}

    def __init__(self, line_iterator: Iterable[bytes],
                 parent_relationship_types: Union[str, Iterable[str]] = ('is_a', 'part_of')):
        self.data_version = None
        self.go_terms: Dict[str, GOTerm] = {}
        self.alt_ids: Dict[str, str] = {}
        self.levels: List[dict] = []
        self.parent_relationship_types: list = data_to_list(parent_relationship_types)

        self._parse_file(line_iterator)
        self._populate_levels()
        self._populate_children()

    def __getitem__(self, key) -> 'GOTerm':
        if key in self.go_terms:
            return self.go_terms[key]
        elif key in self.alt_ids:
            return self.go_terms[self.alt_ids[key]]
        raise KeyError(key)

    def __contains__(self, item):
        try:
            _ = self[item]
            return True
        except KeyError:
            return False

    def _parse_file(self, line_iterator: Iterable[bytes]):
        current_term = None
        in_frame = False
        for line in line_iterator:
            if in_frame:
                if line.startswith(b'id: '):
                    current_term.id = parse_go_id(line)
                elif line.startswith(b'name: '):
                    current_term.name = line[6:].decode('utf8')
                elif line.startswith(b'alt_id: '):
                    self.alt_ids[parse_go_id(line)] = current_term.id
                elif line.startswith(b'is_a: '):
                    current_term.relationships['is_a'].append(parse_go_id(line))
                elif line.startswith(b'relationship: '):
                    relationship_type = line.split(b' ')[1].decode('utf8')
                    if relationship_type not in current_term.relationships:
                        current_term.relationships[relationship_type] = []
                    current_term.relationships[relationship_type].append(parse_go_id(line))
                elif line.startswith(b'is_obsolete: true'):
                    in_frame = False
                elif line == b'' or line == b'\n':
                    self.go_terms[current_term.id] = current_term
                    in_frame = False
            else:
                if line.startswith(b'[Term]'):
                    current_term = GOTerm()
                    in_frame = True
                elif line.startswith(b'data-version:'):
                    self.data_version = line[14:].decode('utf8').replace('\n', '')

    def _populate_levels(self):
        levels_dict = {}
        for go_term in self.go_terms.values():
            if go_term.level is None:
                go_term.level = self._get_term_level_rec(go_term)
            if go_term.level not in levels_dict:
                levels_dict[go_term.level] = {}
            levels_dict[go_term.level][go_term.id] = go_term

        self.levels = [levels_dict[i] for i in range(0, max(levels_dict.keys()) + 1)]

    def _get_term_level_rec(self, go_term: GOTerm):
        if go_term.level is not None:
            pass
        elif len(go_term.get_parents(self.parent_relationship_types)) == 0:
            go_term.level = 0
        else:
            go_term.level = 1 + max(
                [self._get_term_level_rec(self[parent_id]) for parent_id in
                 go_term.get_parents(self.parent_relationship_types)])
        return go_term.level

    def _populate_children(self):
        for go_id in self.level_iter():
            for rel_type in self.parent_relationship_types:
                for parent_id in self[go_id].get_parents(rel_type):
                    if rel_type not in self[parent_id].children_relationships:
                        self[parent_id].children_relationships[rel_type] = []
                    self[parent_id].children_relationships[rel_type].append(go_id)

    def level_iter(self):
        for level in self.levels[::-1]:
            for go_id in level:
                yield go_id

    def upper_induced_graph_iter(self, go_id: str):
        node_queue = queue.SimpleQueue()
        processed_nodes = set()
        for parent in self[go_id].get_parents(self.parent_relationship_types):
            node_queue.put(parent)
        while not node_queue.empty():
            this_node = node_queue.get()
            parents = self[this_node].get_parents(self.parent_relationship_types)
            for parent in parents:
                if parent not in processed_nodes:
                    node_queue.put(parent)
                processed_nodes.update(parents)
            yield this_node
