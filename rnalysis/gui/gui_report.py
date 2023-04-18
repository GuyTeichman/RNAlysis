import typing
from pathlib import Path

from pyvis.network import Network


class ReportGenerator:
    def __init__(self):
        self.graph = Network(directed=True, layout=True)
        self.nodes = {}

        self.add_node('Started RNAlysis session', 0, [])

    def add_node(self, name: str, node_id: int, predecessors: typing.List[int] = (0,), popup_element: str = ''):
        assert node_id not in self.nodes, f"Node {node_id} already exists!"
        self.nodes[node_id] = name
        self.graph.add_node(node_id, name, title=popup_element)
        for pred in predecessors:
            self.graph.add_edge(pred, node_id)

    def generate_report(self, save_path: Path):
        assert save_path.exists() and save_path.is_dir()
        save_file = save_path.joinpath('report.html').as_posix()
        self.graph.show(save_file, notebook=False)
