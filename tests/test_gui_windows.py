from rnalysis import __version__
from rnalysis.gui.gui_windows import *

LEFT_CLICK = QtCore.Qt.LeftButton
RIGHT_CLICK = QtCore.Qt.RightButton


def widget_setup(qtbot, widget_class, *args, **kwargs):
    widget = widget_class(*args, **kwargs)
    widget.show()
    qtbot.add_widget(widget)
    return qtbot, widget


def test_AboutWindow(qtbot, monkeypatch):
    exit_calls = []

    def mock_close(*args, **kwargs):
        exit_calls.append(1)

    monkeypatch.setattr(AboutWindow, 'close', mock_close)

    qtbot, window = widget_setup(qtbot, AboutWindow)
    qtbot.mouseClick(window.buttons()[0], LEFT_CLICK)
    assert exit_calls == [1]


def test_ErrorMessage_message(qtbot):
    err_text = 'my error text'
    try:
        raise ValueError(err_text)
    except ValueError as e:
        err_tb = e.__traceback__
        err_value = e
    qtbot, dialog = widget_setup(qtbot, ErrorMessage, ValueError, err_value, err_tb)
    assert 'ValueError' in dialog.widgets['error_text'].toPlainText()
    assert err_text in dialog.widgets['error_text'].toPlainText()


def test_ErrorMessage_close(qtbot, monkeypatch):
    closed = []

    def mock_close(*args, **kwargs):
        closed.append(1)

    monkeypatch.setattr(ErrorMessage, 'close', mock_close)
    err_text = 'my error text'
    try:
        raise ValueError(err_text)
    except ValueError as e:
        err_tb = e.__traceback__
        err_value = e
    qtbot, dialog = widget_setup(qtbot, ErrorMessage, ValueError, err_value, err_tb)
    qtbot.mouseClick(dialog.widgets['ok_button'], LEFT_CLICK)

    assert closed == [1]


def test_ErrorMessage_copy_to_clipboard(qtbot, monkeypatch):
    err_text = 'my error text'
    try:
        raise ValueError(err_text)
    except ValueError as e:
        err_tb = e.__traceback__
        err_value = e
    qtbot, dialog = widget_setup(qtbot, ErrorMessage, ValueError, err_value, err_tb)
    qtbot.mouseClick(dialog.widgets['copy_button'], LEFT_CLICK)

    assert 'ValueError' in QtWidgets.QApplication.clipboard().text()
    assert err_text in QtWidgets.QApplication.clipboard().text()


def test_HowToCiteWindow(qtbot, monkeypatch):
    exit_calls = []

    def mock_close(*args, **kwargs):
        exit_calls.append(1)

    monkeypatch.setattr(HowToCiteWindow, 'close', mock_close)

    qtbot, window = widget_setup(qtbot, HowToCiteWindow)
    qtbot.mouseClick(window.ok_button, LEFT_CLICK)
    assert exit_calls == [1]


def test_HowToCiteWindow_copy_to_clipboard(qtbot, monkeypatch):
    qtbot, dialog = widget_setup(qtbot, HowToCiteWindow)
    qtbot.mouseClick(dialog.copy_button, LEFT_CLICK)

    cb = QtWidgets.QApplication.clipboard().text()
    assert 'RNAlysis' in cb
    assert 'version' in cb
    assert str(__version__) in cb


def test_MultiFileSelectionDialog_init(qtbot):
    assert False


def test_MultiFileSelectionDialog_selection(qtbot):
    assert False


def test_MultiFileSelectionDialog_no_selection(qtbot):
    assert False


def test_DataFramePreviewModel(qtbot):
    assert False


def test_GeneSetView_init(qtbot):
    assert False


def test_GeneSetView_save(qtbot):
    assert False


def test_DataFrameView_init(qtbot):
    assert False


def test_DataFrameView_save(qtbot):
    assert False


def test_SettingsWindow_init(qtbot):
    assert False


def test_SettingsWindow_set_defaults(qtbot):
    assert False


def test_SettingsWindow_reset(qtbot):
    assert False


def test_SettingsWindow_save(qtbot):
    assert False
