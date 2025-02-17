import itertools
import json
import re
import shutil
import typing
import webbrowser
from pathlib import Path
from typing import Literal, Union

import networkx
from pyvis.network import Network

from rnalysis import __version__
from rnalysis.gui import gui_windows
from rnalysis.utils import io, parsing


class Node:
    __slots__ = ('_node_id', '_node_name', '_predecessors', '_is_active', '_popup_element', '_node_type', '_filename')
    DATA_TYPES = {'Count matrix', 'Differential expression', 'Fold change', 'Other table', 'Gene set', 'Other output',
                  'Pipeline'}

    def __init__(self, node_id: int, node_name: str, predecessors: list, popup_element: str, node_type: str,
                 filename: str = None):
        self._node_id = node_id
        self._node_name = node_name
        self._predecessors = parsing.data_to_set(predecessors)
        self._popup_element = popup_element
        self._node_type = node_type
        self._is_active = True
        self._filename = None if filename is None else Path(filename)
        if node_type in self.DATA_TYPES:
            self._node_name += f' (#{node_id})'

    def to_json(self):
        # Convert the Node object to a JSON-compatible dictionary
        data = {
            "node_id": self.node_id,
            "node_name": self.node_name,
            "predecessors": parsing.data_to_list(self.predecessors),
            "popup_element": self.popup_element,
            "node_type": self.node_type,
            "filename": str(self.filename) if self.filename else None
        }
        return json.dumps(data)

    @classmethod
    def from_json(cls, json_str):
        # Convert the JSON string back to a dictionary
        data = json.loads(json_str)
        # Create a new Node object using the dictionary
        return cls(
            node_id=data["node_id"],
            node_name=data["node_name"],
            predecessors=parsing.data_to_set(data["predecessors"]),
            popup_element=data["popup_element"],
            node_type=data["node_type"],
            filename=Path(data["filename"]) if data["filename"] else None
        )

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
    def filename(self) -> typing.Union[None, Path]:
        return self._filename

    def set_active(self, is_active: bool):
        self._is_active = is_active

    def add_predecessor(self, pred: int):
        self._predecessors.add(pred)

    def update_filename(self, filename: typing.Union[None, Path]):
        self._filename = filename


class ReportGenerator:
    __slots__ = ('graph', 'nodes')
    CSS_TEMPLATE_PATHS = [Path(__file__).parent.parent.joinpath('data_files/report_templates/vis-network.min.css'),
                          Path(__file__).parent.parent.joinpath('data_files/report_templates/bootstrap.min.css')]
    JS_TEMPLATE_PATHS = [Path(__file__).parent.parent.joinpath('data_files/report_templates/vis-network.min.js'),
                         Path(__file__).parent.parent.joinpath('data_files/report_templates/bootstrap.bundle.min.js')]
    OTHER_PATHS = [Path(__file__).parent.parent.joinpath('data_files/report_templates/question-circle.svg')]
    NODE_GROUPS = {'root': 'root',
                   'Count matrix': 'count',
                   'Differential expression': 'diffexp',
                   'Fold change': 'foldchange',
                   'Other table': 'table',
                   'Gene set': 'geneset',
                   'Function': 'function',
                   'Other output': 'other',
                   'Pipeline': 'pipeline'}
    ROOT_FNAME = 'session.rnal'
    TITLE = f"Data analysis report (<i>RNAlysis</i> version {__version__})"

    def __init__(self):
        self.graph = networkx.DiGraph()
        self.nodes: typing.Dict[int, Node] = {}
        self.create_legend()
        href = Path('data').joinpath(self.ROOT_FNAME).as_posix()
        root_desc = (f'<i>RNAlysis</i> version {__version__}<br>'
                     f'<a href="{href}" target="_blank" rel="noopener noreferrer">Open RNAlysis session</a>')
        self.add_node('Started RNAlysis session', 0, [], root_desc, node_type='root', filename=self.ROOT_FNAME)

    def create_legend(self):
        x = -750
        y = -350
        step = 75
        for node_type, group_id in self.NODE_GROUPS.items():
            if node_type in {'root'}:
                continue
            self.graph.add_node(node_type, group=group_id, is_legend=True, label=node_type.capitalize(), fixed=True,
                                physics=False, x=x, y=y, font={'size': 16}, widthConstraint=100, shape=None)
            y += step

    def add_node(self, name: str, node_id: int, predecessors: typing.List[int] = tuple(), popup_element: str = '',
                 node_type: Literal[tuple(NODE_GROUPS)] = 'Other table', filename: str = None):
        # parentless nodes should be attached to the root node
        if len(predecessors) == 0 and node_id > 0:
            predecessors = [0]

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
        group_id = self.NODE_GROUPS[node.node_type]
        self.graph.add_node(node.node_id, label=node.node_name, title=node.popup_element, group=group_id, shape=None)
        for pred in predecessors:
            self.graph.add_edge(pred, node_id)

    def trim_node(self, node_id: int):
        """
        Removes a node from the report graph if it is a leaf node with no outgoing edges, \
        and trims any of its predecessors that are of type 'Function' and are also now leaf nodes.

        :param node_id: The ID of the node to be trimmed.
        :type node_id: int

        This method checks if the specified node is present in the report graph, and if it is a leaf node \
        with no outgoing edges. If both conditions are met, the node is removed from the graph, \
        its corresponding `ReportNode` object is set to inactive, and any of its predecessors that are of type \
        'Function' are recursively checked to see if they are now leaf nodes themselves (i.e., have no outgoing edges).\
         If a predecessor is a leaf node, it is also trimmed from the graph and marked as inactive.

        This method is used for pruning the report graph of unnecessary nodes when a tab is closed \
        or an action is undone.
            """
        if node_id in self.graph and self.graph.out_degree(node_id) == 0:
            predecessors = self.graph.predecessors(node_id)
            self.graph.remove_node(node_id)
            self.nodes[node_id].set_active(False)
            for pred in predecessors:
                if self.nodes[pred].node_type == 'Function':
                    self.trim_node(pred)

    def trim_function_nodes(self):
        """
        Trims all nodes of type 'Function' that are now leaf nodes from the report graph.
        """
        for node in self.nodes:
            if self.nodes[node].node_type == 'Function' and node in self.graph and self.graph.out_degree(node) == 0:
                self.trim_node(node)

    def _modify_html(self, html: str, title: str, title_fontsize: int) -> str:
        # remove duplicate title
        title = self.TITLE if title == 'auto' else title
        if html.count(title) > 1:
            html = re.sub(r'<center>.+?<\/h1>\s+<\/center>', '', html, 1, re.DOTALL)
        # set title font size
        fontsize_em = title_fontsize / 16
        html = html.replace('<h1>', f'<h1 style="font-size:{fontsize_em}em;">')
        # add "how to use" button
        with open(Path(__file__).parent.parent.joinpath('data_files/report_misc/report_howto_link.html')) as f:
            howto_link = f.read()
            html = re.sub('</center>', howto_link + '\n</center>', html, count=1)
        # remove comments from file
        comment_regex = r"<!--[\s\S]*?-->"
        html = re.sub(comment_regex, "", html)
        # set CSS templates to the correct paths (local version under "assets")
        for css_pth in self.CSS_TEMPLATE_PATHS:
            css_line = f'<link rel="stylesheet" href="assets/{css_pth.name}"/>'
            html = re.sub(r'<link(?:\s+[\w-]+="[^"]*")*\s+href="[^"]+"\s+(?:[\w-]+="[^"]*"\s+)*?\/>', css_line, html, 1,
                          re.DOTALL)
        # set JavaScript templates to the correct paths (local version under "assets")
        for js_pth in self.JS_TEMPLATE_PATHS:
            # change suffix from .js to .jscript, so that services such as Gmail do not block the file
            jscript_name = js_pth.name.replace('.js', '.jscript')
            js_line = f'<script src="assets/{jscript_name}"></script>'
            html = re.sub(r'<script\s+src\s*=\s*"(https?:\/\/[^"]+\.js)"[^>]*><\/script>', js_line, html, 1, re.DOTALL)

        # add path highlighting script to the HTML file
        with open(Path(__file__).parent.parent.joinpath('data_files/report_misc/globalVars.js')) as f:
            global_vars = f.read()
            html = re.sub(r'(var\s+[^;]+;\s*)+', r'\g<0>' + global_vars, html, count=1)
        with open(Path(__file__).parent.parent.joinpath('data_files/report_misc/listeners.js')) as f:
            listeners = f.read()
        with open(Path(__file__).parent.parent.joinpath('data_files/report_misc/drawPath.js')) as f:
            draw_path_func = f.read()
        merged_code = listeners + '\n' + draw_path_func
        html = re.sub(r'(network\s*=\s*new\s+vis\.Network\([^;]+;\s*)', r'\g<0>' + merged_code, html, count=1)
        # point all table links to their .csv version instead of the .parquet version
        html = html.replace('.parquet', '.csv')
        return html

    def _report_from_nx(self, show_settings: bool, title: Union[str, Literal['auto']],
                        hierarchical_layout: bool) -> Network:
        vis_report = Network(directed=True, layout=False, heading=self.TITLE if title == 'auto' else title)
        vis_report.from_nx(self.graph)

        with open(Path(__file__).parent.parent.joinpath('data_files/report_misc/options.json')) as f:
            options = json.load(f)
            options['configure']['enabled'] = show_settings
            options['layout']['hierarchical']['enabled'] = hierarchical_layout
            vis_report.set_options(json.dumps(options))
        return vis_report

    @staticmethod
    def generate_report_dialog(parent=None):
        dialog = ConfigureReportWindow(parent)
        return dialog

    def generate_report(self, output_folder: Path, title: Union[str, Literal['auto']] = 'auto',
                        title_fontsize: int = 24, show_settings_menu: bool = False, hierarchical_layout: bool = True):
        output_folder = Path(output_folder)
        assert output_folder.exists() and output_folder.is_dir()
        save_file = output_folder.joinpath('report.html').as_posix()
        self.trim_function_nodes()
        vis_report = self._report_from_nx(show_settings_menu, title, hierarchical_layout)
        html = self._modify_html(vis_report.generate_html(save_file), title, title_fontsize)

        with open(save_file, 'w') as f:
            f.write(html)

        assets_path = output_folder.joinpath('assets')
        if assets_path.exists():
            shutil.rmtree(assets_path)
        assets_path.mkdir()
        for item in itertools.chain(self.CSS_TEMPLATE_PATHS, self.JS_TEMPLATE_PATHS, self.OTHER_PATHS):
            with open(item, encoding="utf-8") as f:
                content = f.read()
            # change suffix from .js to .jscript, so that services such as Gmail do not block the file
            outfile_path = assets_path.joinpath(item.name.replace('.js', '.jscript'))
            with open(outfile_path, 'w', encoding="utf-8") as outfile:
                outfile.write(content)

        data_path = output_folder.joinpath('data')
        if data_path.exists():
            shutil.rmtree(data_path)
        data_path.mkdir()

        for ind, node in self.nodes.items():
            if ind == 0:  # skip the root node
                continue
            if not node.is_active or node.filename is None:  # skip inactive nodes, or nodes with no data file
                continue

            outfile_path = data_path.joinpath(node.filename)
            suffix = node.filename.suffix.lower()
            if suffix == '.parquet':  # convert parquet files to csv files so users can open them with Excel and such
                df = io.load_cached_gui_file(node.filename)
                io.save_table(df, outfile_path.with_suffix('.csv'))

            content: bytes = io.load_cached_gui_file(node.filename, load_as_obj=False)
            if content is not None:
                with open(outfile_path, 'wb') as f:
                    f.write(content)
        webbrowser.open(save_file)

    def serialize(self):
        self.trim_function_nodes()
        data = {'graph': networkx.node_link_data(self.graph),
                'nodes': {ind: node.to_json() for ind, node in self.nodes.items()}}
        current_file_paths = {ind: node.filename for ind, node in self.nodes.items() if node.filename is not None}
        current_file_paths.pop(0)  # do not reference session file to avoid infinite recursion
        return data, current_file_paths

    @classmethod
    def deserialize(cls, data: dict):
        obj = cls.__new__(cls)
        obj.graph = networkx.node_link_graph(data['graph'])
        obj.nodes = {ind: Node.from_json(node_json) for ind, node_json in data['nodes'].items()}
        return obj


class ConfigureReportWindow(gui_windows.FuncExternalWindow):
    __slots__ = {}
    EXCLUDED_PARAMS = {'self'}

    def __init__(self, parent=None):
        func = ReportGenerator.generate_report
        super().__init__('Generate analysis report', func, None, self.EXCLUDED_PARAMS, threaded=True,
                         parent=parent)
        self.init_ui()
        self.setWindowTitle('Generate analysis report')
