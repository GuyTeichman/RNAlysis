import os

import pytest
from rnalysis.gui.gui_windows import *
from utils.io import save_gene_set

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
