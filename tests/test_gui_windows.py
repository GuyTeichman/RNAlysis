import numpy as np
import pytest

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


@pytest.mark.parametrize('df,shape_truth', [
    (pd.DataFrame([[1, 2, 3], [4, 5, 6], [7, 8, 9]]), (3, 3)),
    (pd.DataFrame([[1, 2, 3, 0], [4, 5, 6, 0]]), (2, 4)),
    (pd.DataFrame([[1, 2, 3], [4, 5, 6], [7, 8, 9], [10, 11, 12], [13, 14, 15]]), (3, 3)),
    (pd.Series([1, 2]), (2, 1)),
    (pd.Series([1, 2, 3, 4, 5]), (3, 1)),
    (pd.Series(), (0, 1)),
    (pd.DataFrame(), (0, 0))
])
def test_DataFramePreviewModel(qtbot, df, shape_truth):
    model = DataFramePreviewModel(df)

    assert model._dataframe.shape == shape_truth
    if shape_truth[0] > 2:
        assert np.all(model._dataframe.iloc[-1, :] == "...")
    if shape_truth[1] > 3:
        assert np.all(model._dataframe.iloc[:, -1] == "...")


@pytest.mark.parametrize('gene_set', [{1, 2, 3}, {'a', 'b', 'c', 'd'}, set()])
def test_GeneSetView_init(qtbot, gene_set):
    qtbot, dialog = widget_setup(qtbot, GeneSetView, gene_set, 'my set name')
    assert 'my set name' in dialog.label.text()
    assert str(len(gene_set)) in dialog.label.text()
    assert dialog.data_view.count() == len(gene_set)


@pytest.mark.parametrize('gene_set,truth', [
    ({1, 2, 3}, '1\n2\n3'),
    ({'a', 'b', 'c', 'd'}, 'a\nb\nc\nd'),
    (set(), '')])
def test_GeneSetView_save(qtbot, gene_set, truth, monkeypatch):
    pth = 'tests/test_files/my_gene_set_saved_file.txt'

    def get_savepath(*args, **kwargs):
        return pth, '.txt'

    monkeypatch.setattr(QtWidgets.QFileDialog, 'getSaveFileName', get_savepath)
    qtbot, dialog = widget_setup(qtbot, GeneSetView, gene_set, 'my set name')
    try:
        qtbot.mouseClick(dialog.save_button, LEFT_CLICK)
        assert Path(pth).exists()
        with open(pth) as f:
            assert sorted(f.read().split('\n')) == truth.split('\n')

    finally:
        if Path(pth).exists():
            Path(pth).unlink()


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
