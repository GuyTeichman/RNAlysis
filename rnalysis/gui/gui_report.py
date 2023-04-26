import re
import shutil
import typing
import webbrowser
from pathlib import Path

import networkx
from pyvis.network import Network
from typing_extensions import Literal

from rnalysis import __version__
from rnalysis.utils import parsing, io


class Node:
    __slots__ = ('_node_id', '_node_name', '_predecessors', '_is_active', '_popup_element', '_node_type', '_filename')
    DATA_TYPES = {'Count matrix', 'Differential expression', 'Fold change', 'Other table', 'Gene set', 'dataframe'}

    def __init__(self, node_id: int, node_name: str, predecessors: list, popup_element: str, node_type: str,
                 filename: str = None):
        self._node_id = node_id
        self._node_name = node_name
        self._predecessors = parsing.data_to_set(predecessors)
        self._popup_element = popup_element
        self._node_type = node_type
        self._is_active = True
        self._filename = None if filename is None else Path(filename)
        self._filename = filename

        if filename is not None:
            href = Path('data').joinpath(filename).as_posix()
            self._popup_element += f'<br><a href="{href}" target="_blank" rel="noopener noreferrer">Open file</a>'

        if node_type in self.DATA_TYPES:
            self._node_name += f' (#{node_id})'

    @property
    def node_id(self) -> int:
        return self._node_id

    @property
    def node_name(self) -> str:
        return self._node_name

    @property
    def predecessors(self) -> typing.Set[int]:
        return self._predecessors

    @property
    def popup_element(self) -> str:
        return self._popup_element

    @property
    def node_type(self) -> str:
        return self._node_type

    @property
    def is_active(self) -> bool:
        return self._is_active

    @property
    def filename(self) -> str:
        return self._filename

    def set_active(self, is_active: bool):
        self._is_active = is_active

    def add_predecessor(self, pred: int):
        self._predecessors.add(pred)


class ReportGenerator:
    CSS_TEMPLATE_PATH = Path(__file__).parent.parent.joinpath('data_files/report_templates/vis-network.min.css')
    JS_TEMPLATE_PATH = Path(__file__).parent.parent.joinpath('data_files/report_templates/vis-network.min.js')
    NODE_STYLES = {'root': dict(shape='box', color='#00D4D8'),
                   'function': dict(shape='triangleDown', color='#00D4D8'),
                   'dataframe': dict(shape='square', color='#228B22'),
                   'Count matrix': dict(color='#0D47A1'),
                   'Differential expression': dict(color='#BF360C'),
                   'Fold change': dict(color='#00838F'),
                   'Other table': dict(color='#F7B30A'),
                   'Gene set': dict(color='#BA68C8')}

    def __init__(self):
        self.graph = networkx.DiGraph()
        self.nodes: typing.Dict[int, Node] = {}
        self.create_legend()
        self.add_node('Started RNAlysis session', 0, [], node_type='root', filename='session.rnal')

    def create_legend(self):
        x = -30
        y = -20
        step = 10
        level = 1
        for node_type, kwargs in self.NODE_STYLES.items():
            if node_type in {'root'}:
                continue
            self.graph.add_node(node_type, label=node_type.capitalize(), fixed=True, physics=False, **kwargs)
            y += step
            level += 1

    def add_node(self, name: str, node_id: int, predecessors: typing.List[int] = (0,), popup_element: str = '',
                 node_type: Literal[tuple(NODE_STYLES)] = 'Other table', filename: str = None):
        if node_id in self.nodes:
            if self.nodes[node_id].is_active:
                return
            node = self.nodes[node_id]
            node.set_active(True)

            for pred in predecessors:
                if not self.nodes[pred].is_active:
                    self.add_node('', pred)
        else:
            node = Node(node_id, name, predecessors, popup_element, node_type, filename)
            self.nodes[node_id] = node
        kwargs = self.NODE_STYLES[node_type]
        self.graph.add_node(node.node_id, label=node.node_name, title=node.popup_element, **kwargs)
        for pred in predecessors:
            self.graph.add_edge(pred, node_id)

    def trim_node(self, node_id: int):
        if node_id in self.graph and self.graph.out_degree(node_id) == 0:
            predecessors = self.graph.predecessors(node_id)
            self.graph.remove_node(node_id)
            self.nodes[node_id].set_active(False)
            for pred in predecessors:
                if self.nodes[pred].node_type == 'function':
                    self.trim_node(pred)

    def _modify_html(self, html: str, title: str) -> str:
        if html.count(title) > 1:
            html = re.sub(r'<center>.+?<\/h1>\s+<\/center>', '', html, 1, re.DOTALL)

        css_line = f'<link rel = "stylesheet" href="{self.CSS_TEMPLATE_PATH.name}"/>'
        js_line = f'<script src="{self.JS_TEMPLATE_PATH.name}"></script>'

        html = re.sub(r'<link\s+rel="stylesheet"\s+href\s*=\s*"([^"]+)"[^>]*>', css_line, html, 1, re.DOTALL)
        html = re.sub(r'<script\s+src\s*=\s*"(https?:\/\/[^"]+\.js)"[^>]*><\/script>', js_line, html, 1, re.DOTALL)

        return html

    def generate_report(self, save_path: Path, show_buttons: bool = True):
        assert save_path.exists() and save_path.is_dir()
        save_file = save_path.joinpath('report.html').as_posix()
        title = f"Data analysis report (<i>RNAlysis</i> version {__version__})"
        vis_report = Network(directed=True, layout=True, heading=title)
        vis_report.from_nx(self.graph)
        enabled_str = 'true' if show_buttons else 'false'

        vis_report.set_options("""const options = {
    "configure": {"""
                               f'"enabled": {enabled_str}'
                               """
    },
    "layout": {
        "hierarchical": {
            "enabled": true,
            "levelSeparation": 250,
            "nodeSpacing": 250,
            "treeSpacing": 250,
            "direction": "LR",
            "sortMethod": "directed"
        }
    },
    "physics": {
        "hierarchicalRepulsion": {
            "centralGravity": 0,
            "avoidOverlap": null
        },
        "minVelocity": 0.75,
        "solver": "hierarchicalRepulsion"
    },
    "interaction": {
    "navigationButtons": true
    }
}""")
        html = self._modify_html(vis_report.generate_html(save_file), title)

        with open(save_file, 'w') as f:
            f.write(html)

        for item in [self.CSS_TEMPLATE_PATH, self.JS_TEMPLATE_PATH]:
            with open(item, encoding="utf-8") as f:
                content = f.read()
            with open(save_path.joinpath(item.name), 'w', encoding="utf-8") as outfile:
                outfile.write(content)

        tables_path = save_path.joinpath('data')
        if tables_path.exists():
            shutil.rmtree(tables_path)
        tables_path.mkdir()

        for ind, node in self.nodes.items():
            if ind == 0:  # skip the root node
                continue
            if node.is_active and node.filename is not None:
                content = io.load_cached_gui_file(node.filename, load_as_obj=False)
                if content is not None:
                    with open(tables_path.joinpath(node.filename), 'w') as f:
                        f.write(content)
        webbrowser.open(save_file)
