import pytest
from rnalysis.gui.gui import *

LEFT_CLICK = QtCore.Qt.LeftButton
RIGHT_CLICK = QtCore.Qt.RightButton


@pytest.fixture
def use_temp_settings_file(request):
    settings.make_temp_copy_of_settings_file()
    request.addfinalizer(settings.remove_temp_copy_of_settings_file)
    request.addfinalizer(settings.set_temp_copy_of_settings_file_as_default)


def widget_setup(qtbot, widget_class, *args, **kwargs):
    widget = widget_class(*args, **kwargs)
    widget.show()
    qtbot.add_widget(widget)
    return qtbot, widget


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


def test_ApplyPipelineWindow_init(qtbot):
    assert False


def test_MainWindow_init(qtbot, use_temp_settings_file):
    # qtbot, window = widget_setup(qtbot, MainWindow)

    assert False
