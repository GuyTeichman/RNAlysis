import pytest

from rnalysis.gui.gui_report import *


def test_Node():
    # Test creating a node object
    predecessors = [1, 2, 3]
    node = Node(4, 'Node 4', predecessors, 'popup', 'Count matrix', 'filename.txt')
    assert node.node_id == 4
    assert node.node_name == 'Node 4 (#4)'
    assert node.predecessors == set(predecessors)
    assert node.popup_element == 'popup'
    assert node.node_type == 'Count matrix'
    assert node.is_active is True
    assert node.filename == 'filename.txt'

    # Test setting and getting is_active property
    node.set_active(False)
    assert node.is_active is False
    node.set_active(True)
    assert node.is_active is True

    # Test adding a predecessor
    node.add_predecessor(5)
    assert node.predecessors == {1, 2, 3, 5}


@pytest.fixture
def report_generator():
    return ReportGenerator()


def test_create_legend(report_generator):
    assert len(report_generator.graph.nodes) == 9


def test_add_node(report_generator):
    report_generator.add_node("test", 1)
    assert len(report_generator.graph.nodes) == 10


def test_trim_node(report_generator):
    report_generator.add_node("test", 1)
    report_generator.trim_node(1)
    assert len(report_generator.graph.nodes) == 9


def test_modify_html(report_generator):
    html = report_generator._report_from_nx(True).generate_html()
    modified_html = report_generator._modify_html(html)
    assert 'vis-network.min.css' in modified_html
    assert 'vis-network.min.js' in modified_html
    assert modified_html.count(report_generator.TITLE) == 1


def test_generate_report(report_generator, tmp_path):
    report_generator.generate_report(tmp_path)
    report_file = tmp_path / 'report.html'
    assert report_file.exists()
    assert len(report_generator.graph.nodes) == 9
