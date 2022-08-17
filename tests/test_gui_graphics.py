import pytest
from rnalysis.gui.gui_graphics import *

LEFT_CLICK = QtCore.Qt.LeftButton
RIGHT_CLICK = QtCore.Qt.RightButton


@pytest.fixture
def two_gene_sets():
    return {'first': {'a', 'b', 'c'}, 'second': {'a', 'c', 'e', 'f'}}


@pytest.fixture
def three_gene_sets():
    return {'first': {'a', 'b', 'c'}, 'second': {'c', 'd', 'e'}, 'third': {'a', 'c', 'e', 'f'}}


@pytest.fixture
def three_gene_sets_with_disjoint():
    return {'first': {'a', 'b', 'c'}, 'second': {'d', 'e', 'f'}, 'third': {'a', 'c', 'e', 'f'}}


@pytest.fixture
def four_gene_sets():
    return {'first': {'a', 'b', 'c'}, 'second': {'c', 'd', 'e'}, 'third': {'a', 'c', 'e', 'g'}, 'four': {'c', 'f', 'a'}}


def widget_setup(qtbot, widget_class, *args, **kwargs):
    widget = widget_class(*args, **kwargs)
    widget.show()
    qtbot.add_widget(widget)
    return qtbot, widget


def test_get_icon_invalid(qtbot):
    icon = get_icon('something invalid')
    assert icon is None


def test_get_icon_blank(qtbot):
    pixmap = QtGui.QPixmap(32, 32)
    pixmap.fill(QtCore.Qt.transparent)
    truth = pixmap.toImage()

    icon = get_icon('blank')

    assert icon.pixmap(32, 32).toImage() == truth


@pytest.mark.parametrize("name,path", [
    ('yellow', 'yellow_icon.png'),
    ('Filter', 'filter_icon.png')])
def test_get_icon(qtbot, name, path):
    full_path = 'gui/icons/' + path
    icon = get_icon(name)
    assert isinstance(icon, QtGui.QIcon)
    assert not icon.isNull()
    assert icon.name() == QtGui.QIcon(full_path).name()


def test_EmptyCanvas_init(qtbot):
    qtbot, widget = widget_setup(qtbot, EmptyCanvas, 'text')


def test_VennInteractiveCanvas_init(qtbot, three_gene_sets, three_gene_sets_with_disjoint, four_gene_sets):
    qtbot, widget = widget_setup(qtbot, VennInteractiveCanvas, three_gene_sets)
    qtbot, widget2 = widget_setup(qtbot, VennInteractiveCanvas, three_gene_sets_with_disjoint)
    with pytest.raises(ValueError):
        qtbot, widget = widget_setup(qtbot, VennInteractiveCanvas, four_gene_sets)


def test_VennInteractiveCanvas_clear_selection(qtbot, three_gene_sets):
    qtbot, widget = widget_setup(qtbot, VennInteractiveCanvas, three_gene_sets)
    widget.clear_selection()
    assert widget.get_custom_selection() == set()


def test_VennInteractiveCanvas_union(qtbot, three_gene_sets, three_gene_sets_with_disjoint):
    qtbot, widget = widget_setup(qtbot, VennInteractiveCanvas, three_gene_sets)
    widget.union()
    assert widget.get_custom_selection() == {'a', 'b', 'c', 'd', 'e', 'f'}
    qtbot, widget = widget_setup(qtbot, VennInteractiveCanvas, three_gene_sets_with_disjoint)
    widget.union()
    assert widget.get_custom_selection() == {'a', 'b', 'c', 'd', 'e', 'f'}


def test_VennInteractiveCanvas_intersection(qtbot, three_gene_sets, three_gene_sets_with_disjoint):
    qtbot, widget = widget_setup(qtbot, VennInteractiveCanvas, three_gene_sets)
    widget.intersection()
    assert widget.get_custom_selection() == {'c'}
    qtbot, widget = widget_setup(qtbot, VennInteractiveCanvas, three_gene_sets_with_disjoint)
    widget.intersection()
    assert widget.get_custom_selection() == set()


@pytest.mark.parametrize('primary_set,expected,expected_disjoint', [
    ('first', {'b'}, {'b'}),
    ('second', {'d'}, {'d'}),
    ('third', {'f'}, set())
])
def test_VennInteractiveCanvas_difference(qtbot, three_gene_sets, three_gene_sets_with_disjoint, primary_set, expected,
                                          expected_disjoint):
    qtbot, widget = widget_setup(qtbot, VennInteractiveCanvas, three_gene_sets)
    widget.difference(primary_set)
    assert widget.get_custom_selection() == expected
    qtbot, widget = widget_setup(qtbot, VennInteractiveCanvas, three_gene_sets_with_disjoint)
    widget.difference(primary_set)
    assert widget.get_custom_selection() == expected_disjoint


def test_VennInteractiveCanvas_symmetric_difference(qtbot, two_gene_sets):
    qtbot, widget = widget_setup(qtbot, VennInteractiveCanvas, two_gene_sets)
    widget.symmetric_difference()
    assert widget.get_custom_selection() == {'b', 'e', 'f'}


@pytest.mark.parametrize('threshold,expected,expected_disjoint', [
    (0, {'a', 'b', 'c', 'd', 'e', 'f'}, {'a', 'b', 'c', 'd', 'e', 'f'}),
    (1, {'c'}, set()),
    (0.3, {'a', 'b', 'c', 'd', 'e', 'f'}, {'a', 'b', 'c', 'd', 'e', 'f'}),
    (0.5, {'a', 'c', 'e'}, {'a', 'c', 'e', 'f'}),
    (0.8, {'c'}, set()),
])
def test_VennInteractiveCanvas_majority_vote_intersection(qtbot, three_gene_sets, three_gene_sets_with_disjoint,
                                                          threshold, expected, expected_disjoint):
    qtbot, widget = widget_setup(qtbot, VennInteractiveCanvas, three_gene_sets)
    widget.majority_vote_intersection(threshold)
    assert widget.get_custom_selection() == expected
    qtbot, widget = widget_setup(qtbot, VennInteractiveCanvas, three_gene_sets_with_disjoint)
    widget.majority_vote_intersection(threshold)
    assert widget.get_custom_selection() == expected_disjoint


def test_UpSetInteractiveCanvas_init(qtbot, four_gene_sets, three_gene_sets_with_disjoint):
    qtbot, widget = widget_setup(qtbot, UpSetInteractiveCanvas, four_gene_sets)
    qtbot, widget2 = widget_setup(qtbot, UpSetInteractiveCanvas, three_gene_sets_with_disjoint)


def test_UpSetInteractiveCanvas_clear_selection(qtbot, four_gene_sets):
    qtbot, widget = widget_setup(qtbot, UpSetInteractiveCanvas, four_gene_sets)
    widget.clear_selection()
    assert widget.get_custom_selection() == set()


def test_UpSetInteractiveCanvas_union(qtbot, four_gene_sets, three_gene_sets_with_disjoint):
    qtbot, widget = widget_setup(qtbot, UpSetInteractiveCanvas, four_gene_sets)
    widget.union()
    assert widget.get_custom_selection() == {'a', 'b', 'c', 'd', 'e', 'f', 'g'}
    qtbot, widget = widget_setup(qtbot, UpSetInteractiveCanvas, three_gene_sets_with_disjoint)
    widget.union()
    assert widget.get_custom_selection() == {'a', 'b', 'c', 'd', 'e', 'f'}


def test_UpSetInteractiveCanvas_intersection(qtbot, four_gene_sets, three_gene_sets_with_disjoint):
    qtbot, widget = widget_setup(qtbot, UpSetInteractiveCanvas, four_gene_sets)
    widget.intersection()
    assert widget.get_custom_selection() == {'c'}
    qtbot, widget = widget_setup(qtbot, UpSetInteractiveCanvas, three_gene_sets_with_disjoint)
    widget.intersection()
    assert widget.get_custom_selection() == set()


@pytest.mark.parametrize('primary_set,expected,expected_disjoint', [
    ('first', {'b'}, {'b'}),
    ('second', {'d'}, {'d'}),
    ('third', {'g'}, set())
])
def test_UpSetInteractiveCanvas_difference(qtbot, four_gene_sets, three_gene_sets_with_disjoint, primary_set, expected,
                                           expected_disjoint):
    qtbot, widget = widget_setup(qtbot, UpSetInteractiveCanvas, four_gene_sets)
    widget.difference(primary_set)
    assert widget.get_custom_selection() == expected
    qtbot, widget = widget_setup(qtbot, UpSetInteractiveCanvas, three_gene_sets_with_disjoint)
    widget.difference(primary_set)
    assert widget.get_custom_selection() == expected_disjoint


@pytest.mark.parametrize('threshold,expected,expected_disjoint', [
    (0, {'a', 'b', 'c', 'd', 'e', 'f', 'g'}, {'a', 'b', 'c', 'd', 'e', 'f'}),
    (1, {'c'}, set()),
    (0.3, {'a', 'c', 'e'}, {'a', 'b', 'c', 'd', 'e', 'f'}),
    (0.55, {'a', 'c'}, {'a', 'c', 'e', 'f'}),
    (0.8, {'c'}, set()),
])
def test_UpSetInteractiveCanvas_majority_vote_intersection(qtbot, four_gene_sets, three_gene_sets_with_disjoint,
                                                           threshold, expected, expected_disjoint):
    qtbot, widget = widget_setup(qtbot, UpSetInteractiveCanvas, four_gene_sets)
    widget.majority_vote_intersection(threshold)
    assert widget.get_custom_selection() == expected
    qtbot, widget = widget_setup(qtbot, UpSetInteractiveCanvas, three_gene_sets_with_disjoint)
    widget.majority_vote_intersection(threshold)
    assert widget.get_custom_selection() == expected_disjoint
