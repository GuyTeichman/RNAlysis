import typing
from pathlib import Path

import networkx
from pyvis.network import Network

from rnalysis.utils import parsing


class Node:
    __slots__ = ('_node_id', '_node_name', '_predecessors', '_is_active', '_popup_element')

    def __init__(self, node_id: int, node_name: str, predecessors: list, popup_element: str):
        self._node_id = node_id
        self._node_name = node_name
        self._predecessors = parsing.data_to_set(predecessors)
        self._popup_element = popup_element
        self._is_active = True

    @property
    def node_id(self):
        return self._node_id

    @property
    def node_name(self):
        return self._node_name

    @property
    def predecessors(self):
        return self._predecessors

    @property
    def popup_element(self):
        return self._popup_element

    @property
    def is_active(self):
        return self._is_active

    def set_active(self, is_active: bool):
        self._is_active = is_active

    def add_predecessor(self, pred: int):
        self._predecessors.add(pred)


class ReportGenerator:
    def __init__(self):
        self.graph = networkx.DiGraph()
        self.nodes: typing.Dict[int, Node] = {}

        self.add_node('Started RNAlysis session', 0, [])

    def add_node(self, name: str, node_id: int, predecessors: typing.List[int] = (0,), popup_element: str = ''):
        if node_id in self.nodes:
            if self.nodes[node_id].is_active:
                return
            self.nodes[node_id].set_active(True)
            predecessors = self.nodes[node_id].predecessors
            name = self.nodes[node_id].node_name
            popup_element = self.nodes[node_id].popup_element
        else:
            self.nodes[node_id] = Node(node_id, name, predecessors)

        self.graph.add_node(node_id, label=name, title=popup_element)
        for pred in predecessors:
            self.graph.add_edge(pred, node_id)

    def trim_node(self, node_id: int):
        if node_id in self.nodes and self.graph.out_degree(node_id) == 0:
            self.graph.remove_node(node_id)
            self.nodes[node_id].set_active(False)

    def generate_report(self, save_path: Path):
        assert save_path.exists() and save_path.is_dir()
        save_file = save_path.joinpath('report.html').as_posix()
        vis_report = Network(directed=True, layout=True)
        vis_report.from_nx(self.graph)
        vis_report.show(save_file, notebook=False)
