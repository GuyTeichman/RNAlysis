import functools
import queue
import re
from typing import Dict, List, Union, Tuple, Iterable, Set

from rnalysis.utils import parsing


class GOTerm:
    __slots__ = {'_id': 'GO ID', '_name': 'GO Term name',
                 '_namespace': 'biological_process, cellular_component or molecular_function',
                 '_level': "GO Term's level in the DAG Tree",
                 'relationships': 'direct parent relationships of the GO Term',
                 'children_relationships': 'direct children relationships of the GO Term'}

    def __init__(self):
        self._id = None
        self._name = None
        self._namespace = None
        self._level = None
        self.relationships: Dict[str, List[str]] = {'is_a': [], 'part_of': []}
        self.children_relationships: Dict[str, List[str]] = {'is_a': [], 'part_of': []}

    @classmethod
    def with_properties(cls, go_id: str, name: str, namespace: str, level: int):
        go_term = cls()
        go_term.set_id(go_id)
        go_term.set_name(name)
        go_term.set_namespace(namespace)
        go_term.set_level(level)
        return go_term

    @property
    def id(self):
        return self._id

    @property
    def name(self):
        return self._name

    @property
    def namespace(self):
        return self._namespace

    @property
    def level(self):
        return self._level

    def set_id(self, go_id: str):
        self._id = go_id

    def set_name(self, name: str):
        self._name = name

    def set_namespace(self, namespace: str):
        self._namespace = namespace

    def set_level(self, level: int):
        self._level = level

    @functools.lru_cache(maxsize=2)
    def get_parents(self, relationships: Union[str, tuple] = ('is_a', 'part_of')) -> List[str]:
        relationships_filt = [rel for rel in parsing.data_to_list(relationships) if rel in self.relationships]
        go_ids = [go_id for rel in relationships_filt for go_id in self.relationships[rel]]
        return go_ids

    @functools.lru_cache(maxsize=2)
    def get_children(self, relationships: Union[str, Tuple[str]] = ('is_a', 'part_of')) -> List[str]:
        relationships_filt = [rel for rel in parsing.data_to_list(relationships) if rel in self.children_relationships]
        go_ids = [go_id for rel in relationships_filt for go_id in self.children_relationships[rel]]
        return go_ids


def parse_go_id(byte_sequence: bytes) -> str:
    return re.findall(b"GO:[0-9]{7}", byte_sequence)[0].decode('utf8')


class DAGTree:
    __slots__ = {'data_version': 'version of the go-basic.obo file',
                 'go_terms': 'dictionary of GO Terms in the DAG Tree',
                 'alt_ids': 'mapping of alternagive GO IDs to their main GO ID',
                 'namespaces': "namespaces included in the DAGTree",
                 'levels': 'list of levels in the DAG Tree',
                 'parent_relationship_types': 'the types of relationships that constitute parenthood in the DAG Tree',
                 '_upper_induced_graphs': 'memoized upper-induced graphs'}

    def __init__(self, line_iterator: Iterable[bytes],
                 parent_relationship_types: Union[str, Iterable[str]] = ('is_a', 'part_of')):
        self.data_version = None
        self.go_terms: Dict[str, GOTerm] = {}
        self.alt_ids: Dict[str, str] = {}
        self.namespaces: Set[str] = set()
        self.levels: List[dict] = []
        self.parent_relationship_types: tuple = parsing.data_to_tuple(parent_relationship_types)

        self._upper_induced_graphs: Dict[str, Set[str]] = {}

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
                    current_term.set_id(parse_go_id(line))
                elif line.startswith(b'namespace: '):
                    current_term.set_namespace(line[11:].decode('utf8').replace('\n', ''))
                    if current_term.namespace not in self.namespaces:
                        self.namespaces.add(current_term.namespace)
                elif line.startswith(b'name: '):
                    current_term.set_name(line[6:].decode('utf8').replace('\n', ''))
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
                elif line in {b'', b'\n', b'\r\n'}:
                    self.go_terms[current_term.id] = current_term
                    in_frame = False
            else:
                if line.startswith(b'[Term]'):
                    current_term = GOTerm()
                    in_frame = True
                elif line.startswith(b'data-version:'):
                    self.data_version = line[14:].decode('utf8').replace('\n', '')

        if in_frame:  # add last go term to the set, if it was not already added
            self.go_terms[current_term.id] = current_term

    def _populate_levels(self):
        levels_dict = {}
        for go_term in self.go_terms.values():
            if go_term.level is None:
                go_term.set_level(self._get_term_level_rec(go_term))
            if go_term.level not in levels_dict:
                levels_dict[go_term.level] = {}
            levels_dict[go_term.level][go_term.id] = go_term
        if len(levels_dict) == 0:
            self.levels = [{}]
        else:
            self.levels = [levels_dict[i] for i in range(0, max(levels_dict.keys()) + 1)]

    def _get_term_level_rec(self, go_term: GOTerm):
        if go_term.level is not None:
            pass
        elif len(go_term.get_parents(self.parent_relationship_types)) == 0:
            go_term.set_level(0)
        else:
            go_term.set_level(1 + max([self._get_term_level_rec(self[parent_id]) for parent_id in
                                       go_term.get_parents(self.parent_relationship_types)]))
        return go_term.level

    def _populate_children(self):
        for go_id in self.level_iter():
            for rel_type in self.parent_relationship_types:
                for parent_id in self[go_id].get_parents(rel_type):
                    if rel_type not in self[parent_id].children_relationships:
                        self[parent_id].children_relationships[rel_type] = []
                    self[parent_id].children_relationships[rel_type].append(go_id)

    def level_iter(self, namespace: str = 'all'):
        if namespace == 'all':
            for level in self.levels[::-1]:
                for go_id in level:
                    yield go_id
        else:
            for level in self.levels[::-1]:
                for go_id in level:
                    if self[go_id].namespace == namespace:
                        yield go_id

    def upper_induced_graph_iter(self, go_id: str):
        if go_id in self._upper_induced_graphs:
            for upper_induced_node in self._upper_induced_graphs[go_id]:
                yield upper_induced_node

        else:
            # put go_id's parents into the queue
            node_queue = queue.SimpleQueue()
            processed_nodes = set()
            parents = self[go_id].get_parents(self.parent_relationship_types)
            for parent in parents:
                node_queue.put(parent)
            processed_nodes.update(parents)
            # iterate over the queue until it is empty (meaning we reached the top of the graph)
            while not node_queue.empty():
                this_node = node_queue.get()
                yield this_node
                # if this_node's upper-induced graph was already calculated, yield those unprocessed nodes
                if this_node in self._upper_induced_graphs:
                    for upper_induced_node in self._upper_induced_graphs[this_node]:
                        if upper_induced_node not in processed_nodes:
                            yield upper_induced_node
                    processed_nodes.update(self._upper_induced_graphs[this_node])
                # if this_node's upper-induced graph was yet to be calculated, add its unprocessed parents to the queue
                else:
                    parents = self[this_node].get_parents(self.parent_relationship_types)
                    for parent in parents:
                        if parent not in processed_nodes:
                            node_queue.put(parent)
                    processed_nodes.update(parents)
            # memoize the function's output for go_id
            self._upper_induced_graphs[go_id] = processed_nodes
