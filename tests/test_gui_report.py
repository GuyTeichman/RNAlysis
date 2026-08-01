import polars as pl
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
    assert node.filename == Path('filename.txt')

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


def test_add_duplicate_node(report_generator):
    report_generator.add_node("Test Node", 1, [0], "Popup Element", "Other table")
    report_generator.add_node("Test Node", 1, [0], "Popup Element", "Other table")
    assert report_generator.graph.number_of_nodes() == 10
    assert report_generator.graph.number_of_edges() == 1


def test_add_inactive_node(report_generator):
    report_generator.add_node("Test Node", 1, [0], "Popup Element", "Other table")
    assert report_generator.graph.number_of_nodes() == 10
    assert report_generator.graph.number_of_edges() == 1
    report_generator.trim_node(1)
    assert report_generator.graph.number_of_nodes() == 9
    assert report_generator.graph.number_of_edges() == 0
    report_generator.add_node("Test Node", 1, [0], "Popup Element", "Other table")
    assert report_generator.graph.number_of_nodes() == 10
    assert report_generator.graph.number_of_edges() == 1


def test_add_inactive_predecessor(report_generator):
    report_generator.nodes[0].set_active(False)
    report_generator.add_node("Test Node", 1, [0], "Popup Element", "Other table")
    assert report_generator.graph.number_of_nodes() == 10
    assert report_generator.graph.number_of_edges() == 1


def test_trim_node(report_generator):
    report_generator.add_node("test", 1)
    report_generator.trim_node(1)
    assert len(report_generator.graph.nodes) == 9


def test_trim_function(report_generator):
    report_generator.add_node("Test Node", 1, [0], "Popup Element", "Function")
    report_generator.add_node("Test Node 2", 2, [1], "Popup Element", "Other table")
    report_generator.trim_node(2)
    assert report_generator.graph.number_of_nodes() == 9
    assert report_generator.graph.number_of_edges() == 0


@pytest.mark.parametrize('show_settings,title,fontsize,hierarchical_layout',
                         [(True, 'Test Title', 18, True),
                          (False, 'auto', 32, True),
                          (True, 'auto', 16, False),
                          (True, 'Test Title2', 24, True)])
def test_modify_html(report_generator, show_settings, title, fontsize, hierarchical_layout):
    html = report_generator._report_from_nx(show_settings, title, hierarchical_layout).generate_html()
    modified_html = report_generator._modify_html(html, title, fontsize)
    for filename in ['vis-network.min.css', 'bootstrap.min.css', 'vis-network.min.js', 'bootstrap.bundle.min.js']:
        assert f'assets/{filename}' in modified_html
    # validate title
    if title == 'auto':
        title = report_generator.TITLE
    assert modified_html.count(title) == 1


def test_generate_report(report_generator, tmp_path):
    report_generator.generate_report(tmp_path)
    report_file = tmp_path / 'report.html'
    assert report_file.exists()
    assert len(report_generator.graph.nodes) == 9


def test_generate_report_streams_parquet_to_csv(report_generator, tmp_path, monkeypatch):
    monkeypatch.setattr(webbrowser, 'open', lambda *a, **k: None)
    calls = []
    real = io.stream_cached_parquet_to_csv

    def spy(cached_filename, dst_path):
        calls.append((cached_filename, dst_path))
        return real(cached_filename, dst_path)

    monkeypatch.setattr(io, 'stream_cached_parquet_to_csv', spy)

    df = pl.DataFrame({'': ['g1', 'g2', 'g3'], 'a': [1, 2, 3], 'b': [4.5, None, 6.5]})
    filename = Path('999_report_stream.parquet')
    io.cache_gui_file(df, str(filename))
    io.flush_gui_cache_writes()
    cache_path = io.get_gui_cache_dir().joinpath(filename)
    try:
        report_generator.add_node('Streamed table', 999, [0], 'desc', 'Count matrix', filename)
        report_generator.generate_report(tmp_path)

        out_csv = tmp_path / 'data' / '999_report_stream.csv'
        assert calls, 'report export must transcode parquet->csv via the streaming helper'
        assert out_csv.exists()
        # the streamed csv must be byte-identical to the eager read_parquet -> write_csv path
        expected = tmp_path / 'expected.csv'
        pl.read_parquet(cache_path).write_csv(expected)
        assert out_csv.read_bytes() == expected.read_bytes()
    finally:
        if cache_path.exists():
            cache_path.unlink()
        io._reset_gui_cache_writer()
