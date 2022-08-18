import pytest
from rnalysis.gui.gui import *

LEFT_CLICK = QtCore.Qt.LeftButton
RIGHT_CLICK = QtCore.Qt.RightButton


@pytest.fixture
def blank_icon():
    pixmap = QtGui.QPixmap(32, 32)
    pixmap.fill(QtCore.Qt.transparent)
    return QtGui.QIcon(pixmap)


@pytest.fixture
def red_icon():
    pixmap = QtGui.QPixmap(32, 32)
    pixmap.fill(QtCore.Qt.red)
    return QtGui.QIcon(pixmap)


@pytest.fixture
def green_icon():
    pixmap = QtGui.QPixmap(32, 32)
    pixmap.fill(QtCore.Qt.green)
    return QtGui.QIcon(pixmap)


@pytest.fixture
def use_temp_settings_file(request):
    settings.make_temp_copy_of_settings_file()
    request.addfinalizer(settings.remove_temp_copy_of_settings_file)
    request.addfinalizer(settings.set_temp_copy_of_settings_file_as_default)


@pytest.fixture
def available_objects_no_tabpages(blank_icon, red_icon, green_icon):
    return {'first tab': (None, blank_icon), 'second tab': (None, red_icon), 'third tab': (None, green_icon),
            'fourth tab': (None, red_icon)}


def widget_setup(qtbot, widget_class, *args, **kwargs):
    widget = widget_class(*args, **kwargs)
    widget.show()
    qtbot.add_widget(widget)
    return qtbot, widget


def test_ApplyPipelineWindow_init(qtbot, available_objects_no_tabpages):
    qtbot, window = widget_setup(qtbot, ApplyPipelineWindow, available_objects_no_tabpages)


def test_ApplyPipelineWindow_select_all(qtbot, available_objects_no_tabpages):
    qtbot, window = widget_setup(qtbot, ApplyPipelineWindow, available_objects_no_tabpages)
    qtbot.mouseClick(window.list.select_all_button, LEFT_CLICK)
    assert window.result() == list(available_objects_no_tabpages.keys())


def test_ApplyPipelineWindow_clear_all(qtbot, available_objects_no_tabpages):
    qtbot, window = widget_setup(qtbot, ApplyPipelineWindow, available_objects_no_tabpages)
    qtbot.mouseClick(window.list.select_all_button, LEFT_CLICK)
    qtbot.mouseClick(window.list.clear_all_button, LEFT_CLICK)
    assert window.result() == []


def test_ClicomWindow_init(qtbot):
    # qtbot, window = widget_setup(qtbot, ClicomWindow)
    assert False


def test_EnrichmentWindow_init(qtbot):
    assert False


def test_SetOperationWindowinit(qtbot):
    assert False


def test_SetVisualizationWindowinit(qtbot):
    assert False


def test_FilterTabPage_init(qtbot):
    assert False


def test_SetTabPage_init(qtbot):
    assert False


def test_FuncTypeStack_init(qtbot):
    assert False


def test_CreatePipelineWindow_init(qtbot):
    assert False


def test_MultiKeepWindow_init(qtbot):
    assert False


def test_MultiOpenWindow_init(qtbot):
    assert False


def test_ReactiveTabWidget_init(qtbot):
    assert False


def test_MainWindow_init(qtbot, use_temp_settings_file):
    # qtbot, window = widget_setup(qtbot, MainWindow)

    assert False
