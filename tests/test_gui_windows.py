import os.path

import numpy as np
import pytest

from rnalysis import __version__
from rnalysis.gui.gui_windows import *

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


@pytest.mark.parametrize('df,shape_truth', [
    (pd.DataFrame([[1, 2, 3], [4, 5, 6], [7, 8, 9]]), (3, 3)),
    (pd.DataFrame([[1, 2, 3, 0], [4, 5, 6, 0]]), (2, 4)),
    (pd.DataFrame([[1, 2, 3], [4, 5, 6], [7, 8, 9], [10, 11, 12], [13, 14, 15]]), (5, 3)),
    (pd.Series([1, 2]), (2, 1)),
    (pd.Series([1, 2, 3, 4, 5]), (5, 1)),
    (pd.Series(), (0, 1)),
    (pd.DataFrame(), (0, 0))
])
def test_DataFrameView_init(qtbot, df, shape_truth):
    qtbot, dialog = widget_setup(qtbot, DataFrameView, df, 'my df name')
    assert 'my df name' in dialog.label.text()
    assert str(shape_truth[0]) in dialog.label.text()
    assert str(shape_truth[1]) in dialog.label.text()
    assert dialog.data_view.model().rowCount() == shape_truth[0]
    assert dialog.data_view.model().columnCount() == shape_truth[1]


@pytest.mark.parametrize('df,shape_truth', [
    (pd.DataFrame([[1, 2, 3], [4, 5, 6], [7, 8, 9]], columns=['a', 'b', 'c']), (3, 3)),
    (pd.DataFrame([[1, 2, 3, 0], [4, 5, 6, 0]], columns=['a', 'b', 'c', 'd']), (2, 4)),
    (pd.DataFrame([[1, 2, 3], [4, 5, 6], [7, 8, 9], [10, 11, 12], [13, 14, 15]], columns=['a', 'b', 'c']), (5, 3)),
    (pd.Series([1, 2], name='a'), (2, 1)),
    (pd.Series([1, 2, 3, 4, 5], name='a'), (5, 1)),
    (pd.Series(), (0, 1)),
    (pd.DataFrame(), (0, 0))
])
def test_DataFrameView_save(qtbot, monkeypatch, df, shape_truth):
    pth = 'tests/test_files/my_dataframe_saved_file.csv'

    def get_savepath(*args, **kwargs):
        return pth, '.csv'

    monkeypatch.setattr(QtWidgets.QFileDialog, 'getSaveFileName', get_savepath)
    qtbot, dialog = widget_setup(qtbot, DataFrameView, df, 'my df name')
    try:
        qtbot.mouseClick(dialog.save_button, LEFT_CLICK)
        assert Path(pth).exists()
        if isinstance(df, pd.DataFrame):
            assert df.equals(pd.read_csv(pth, index_col=0))
        else:
            if len(df) == 0:
                assert df.to_frame().shape == pd.read_csv(pth, index_col=0).shape
            else:
                assert df.to_frame().equals(pd.read_csv(pth, index_col=0))
    finally:
        if Path(pth).exists():
            Path(pth).unlink()


def test_SettingsWindow_init(qtbot, use_temp_settings_file):
    qtbot, dialog = widget_setup(qtbot, SettingsWindow)


def test_SettingsWindow_get_defaults(qtbot, use_temp_settings_file):
    font_truth = 'Arial'
    font_size_truth = '16'
    theme_truth = 'Light'
    show_tutorial_truth = False
    attr_truth = os.path.abspath('tests/test_files/counted.tsv')
    biotype_truth = os.path.abspath('tests/test_files/counted.csv')

    settings.set_gui_settings(font_truth, int(font_size_truth), theme_truth.lower(), show_tutorial_truth)
    settings.set_table_settings(attr_truth, biotype_truth)

    qtbot, dialog = widget_setup(qtbot, SettingsWindow)

    assert dialog.appearance_widgets['app_font'].currentText() == font_truth
    assert dialog.appearance_widgets['app_font_size'].currentText() == font_size_truth
    assert dialog.appearance_widgets['app_theme'].currentText() == theme_truth
    assert dialog.appearance_widgets['show_tutorial'].isChecked() == show_tutorial_truth
    assert dialog.tables_widgets['attr_ref_path'].text() == attr_truth
    assert dialog.tables_widgets['biotype_ref_path'].text() == biotype_truth


def test_SettingsWindow_reset_settings(qtbot, monkeypatch, use_temp_settings_file):
    reset_done = []

    def mock_reset():
        reset_done.append(True)

    monkeypatch.setattr(settings, 'reset_settings', mock_reset)
    qtbot, dialog = widget_setup(qtbot, SettingsWindow)
    qtbot.mouseClick(dialog.button_box.button(QtWidgets.QDialogButtonBox.RestoreDefaults), LEFT_CLICK)
    assert len(reset_done) == 1
    assert reset_done[0]


def test_SettingsWindow_save_settings(qtbot, monkeypatch, use_temp_settings_file):
    settings_saved = [False, False]
    font_truth = 'Arial'
    font_size_truth = '16'
    theme_truth = 'Light'
    show_tutorial_truth = False
    attr_truth = os.path.abspath('tests/test_files/counted.tsv')
    biotype_truth = os.path.abspath('tests/test_files/counted.csv')

    def mock_save_gui(font, font_size, theme, show_tutorial):
        assert font == font_truth
        assert font_size == int(font_size_truth)
        assert theme == theme_truth.lower()
        assert show_tutorial == show_tutorial_truth
        settings_saved[0] = True

    def mock_save_tables(attr, biotype):
        assert attr == attr_truth
        assert biotype == biotype_truth
        settings_saved[1] = True

    monkeypatch.setattr(settings, 'set_gui_settings', mock_save_gui)
    monkeypatch.setattr(settings, 'set_table_settings', mock_save_tables)

    qtbot, dialog = widget_setup(qtbot, SettingsWindow)

    # set gui selections to truth values
    for widget in ['app_font']:
        dialog.appearance_widgets[widget].clear()
    for widget in ['attr_ref_path', 'biotype_ref_path']:
        dialog.tables_widgets[widget].clear()
    qtbot.keyClicks(dialog.appearance_widgets['app_font'], font_truth)
    qtbot.keyClicks(dialog.appearance_widgets['app_font_size'], font_size_truth)
    qtbot.keyClicks(dialog.appearance_widgets['app_theme'], theme_truth)
    dialog.appearance_widgets['show_tutorial'].setChecked(False)

    qtbot.keyClicks(dialog.tables_widgets['attr_ref_path'].file_path, attr_truth)
    qtbot.keyClicks(dialog.tables_widgets['biotype_ref_path'].file_path, biotype_truth)

    qtbot.mouseClick(dialog.button_box.button(QtWidgets.QDialogButtonBox.Apply), LEFT_CLICK)
    assert settings_saved[0]
    assert settings_saved[1]


def test_SettingsWindow_cancel(qtbot, monkeypatch, use_temp_settings_file):
    save_done = []

    def mock_save():
        save_done.append(True)

    monkeypatch.setattr(settings, 'set_gui_settings', mock_save)
    monkeypatch.setattr(settings, 'set_table_settings', mock_save)
    qtbot, dialog = widget_setup(qtbot, SettingsWindow)
    dialog.appearance_widgets['app_font'].setCurrentText("David")
    qtbot.mouseClick(dialog.button_box.button(QtWidgets.QDialogButtonBox.Cancel), LEFT_CLICK)

    assert len(save_done) == 0


def test_MultiFileSelectionDialog_init(qtbot, use_temp_settings_file):
    _, _ = widget_setup(qtbot, MultiFileSelectionDialog)


@pytest.mark.parametrize('pth', ['tests/test_files/test_deseq.csv', 'tests/test_files/test_fastqs'])
def test_MultiFileSelectionDialog_selection_log(qtbot, use_temp_settings_file, pth):
    qtbot, dialog = widget_setup(qtbot, MultiFileSelectionDialog)
    model = dialog.tree_mycomputer.model()

    for parent in reversed(list(Path(pth).absolute().parents)):
        dialog.tree_mycomputer.expand(model.index(str(parent)))

    index = model.index(str(Path(pth).absolute()))
    model.setCheckState(index, True)

    dialog.update_log()
    txt = dialog.logger.toPlainText()
    assert str(Path(pth).name) in txt

    if Path(pth).is_dir():
        for item in Path(pth).rglob('*'):
            assert str(item.absolute()) in txt


@pytest.mark.parametrize('pth', ['tests/test_files/test_deseq.csv', 'tests/test_files/counted.tsv'])
def test_MultiFileSelectionDialog_selection_files(qtbot, use_temp_settings_file, pth):
    qtbot, dialog = widget_setup(qtbot, MultiFileSelectionDialog)
    model = dialog.tree_mycomputer.model()

    for parent in reversed(list(Path(pth).absolute().parents)):
        dialog.tree_mycomputer.expand(model.index(str(parent)))

    index = model.index(str(Path(pth).absolute()))
    model.setCheckState(index, True)
    assert dialog.result() == [str(Path(pth).absolute())]


@pytest.mark.parametrize('pth', ['tests/test_files/test_count_from_folder', 'tests/test_files'])
def test_MultiFileSelectionDialog_selection_folders(qtbot, use_temp_settings_file, pth):
    qtbot, dialog = widget_setup(qtbot, MultiFileSelectionDialog)
    model = dialog.tree_mycomputer.model()

    for parent in reversed(list(Path(pth).absolute().parents)):
        dialog.tree_mycomputer.expand(model.index(str(parent)))

    index = model.index(str(Path(pth).absolute()))
    model.setCheckState(index, True)
    assert sorted(dialog.result()) == sorted(
        [str(Path(child).absolute()) for child in Path(pth).rglob('*') if child.is_file()])


def test_MultiFileSelectionDialog_no_selection(qtbot, use_temp_settings_file):
    qtbot, dialog = widget_setup(qtbot, MultiFileSelectionDialog)

    assert len(dialog.result()) == 0


def test_splash_screen(qtbot):
    splash = splash_screen()
    splash.show()
    qtbot.add_widget(splash)


@pytest.mark.parametrize('pth', ['tests/test_files/test_deseq.csv', 'tests/test_files/fc_1_nan.csv'])
def test_dataframe_model(qtmodeltester, pth):
    model = DataFrameModel(pd.read_csv(pth, index_col=0))
    qtmodeltester.check(model)


@pytest.mark.parametrize('pth', ['tests/test_files/test_deseq.csv', 'tests/test_files/fc_1_nan.csv'])
def test_dataframe_preview_model(qtmodeltester, pth):
    model = DataFramePreviewModel(pd.read_csv(pth, index_col=0))
    qtmodeltester.check(model)

#
# def test_checkable_file_system_model(qtmodeltester):
#     model = CheckableFileSystemModel()
#     model.setRootPath('')
#     qtmodeltester.check(model, force_py=True)
#
#
# def test_proxy_model(qtmodeltester):
#     model = FilterProxy(False)
#     file_system = CheckableFileSystemModel()
#     file_system.setRootPath('')
#     model.setSourceModel(file_system)
#     model.setFilterRegularExpression(r'^(?![.])(?!.*[-_.]$).+')
#     qtmodeltester.check(model, force_py=True)
