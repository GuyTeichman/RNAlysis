import pytest
from rnalysis.gui.gui_graphics import *

LEFT_CLICK = QtCore.Qt.LeftButton
RIGHT_CLICK = QtCore.Qt.RightButton


@pytest.fixture
def gene_sets():
    return {}


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
