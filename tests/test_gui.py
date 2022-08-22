import pytest

from rnalysis import filtering, enrichment
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


@pytest.fixture
def filtertabpage(qtbot):
    qtbot, window = widget_setup(qtbot, FilterTabPage)
    window.start_from_filter_obj(filtering.DESeqFilter('tests/test_files/test_deseq.csv'))
    return window


@pytest.fixture
def filtertabpage_with_undo_stack(qtbot):
    stack = QtWidgets.QUndoStack()
    qtbot, window = widget_setup(qtbot, FilterTabPage, undo_stack=stack)
    window.start_from_filter_obj(filtering.DESeqFilter('tests/test_files/test_deseq_sig.csv'))
    return window, stack


@pytest.fixture
def settabpage_with_undo_stack(qtbot):
    stack = QtWidgets.QUndoStack()
    qtbot, window = widget_setup(qtbot, SetTabPage, 'my set name', {'a', 'b', 'c', 'd'}, undo_stack=stack)
    return window, stack


@pytest.fixture
def pipeline():
    pipeline = filtering.Pipeline('DESeqFilter')
    pipeline.add_function(filtering.DESeqFilter.describe)
    pipeline.add_function(filtering.DESeqFilter.filter_significant, 0.2)
    pipeline.add_function(filtering.DESeqFilter.filter_top_n, 'padj', n=1)
    return pipeline


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
    qtbot, window = widget_setup(qtbot, FilterTabPage)


def test_FilterTabPage_load_file(qtbot):
    obj_truth = filtering.CountFilter('tests/test_files/counted.csv')
    qtbot, window = widget_setup(qtbot, FilterTabPage)
    assert window.is_empty()
    assert not window.basic_widgets['start_button'].isEnabled()

    window.basic_widgets['file_path'].clear()
    qtbot.keyClicks(window.basic_widgets['file_path'].file_path, str(Path('tests/test_files/counted.csv').absolute()))
    qtbot.keyClicks(window.basic_widgets['table_type_combo'], 'Count matrix')
    qtbot.mouseClick(window.basic_widgets['start_button'], LEFT_CLICK)

    assert not window.is_empty()
    assert window.obj() == obj_truth
    assert window.obj_type() == filtering.CountFilter
    assert window.name == 'counted'
    assert window.get_table_type() == 'Count matrix'


def test_FilterTabPage_from_obj(qtbot):
    table_name = 'table name'
    qtbot, window = widget_setup(qtbot, FilterTabPage)
    obj = filtering.DESeqFilter('tests/test_files/test_deseq.csv')
    assert window.is_empty()
    window.start_from_filter_obj(obj, table_name)

    assert not window.is_empty()
    assert window.obj() == obj
    assert window.obj_type() == filtering.DESeqFilter
    assert window.name == table_name
    assert window.get_table_type() == 'Differential expression'


def test_FilterTabPage_cache(qtbot, monkeypatch):
    qtbot, window = widget_setup(qtbot, FilterTabPage)
    filt = filtering.DESeqFilter('tests/test_files/test_deseq.csv')
    window.start_from_filter_obj(filt, 'table name')
    cached = []

    def mock_cache(obj, filename):
        assert isinstance(obj, pd.DataFrame)
        assert obj.equals(filt.df)
        cached.append(True)

    monkeypatch.setattr(io, 'cache_gui_file', mock_cache)

    fname = window.cache()
    assert fname.endswith('.csv')
    assert len(fname) == 44

    time.sleep(0.01)
    fname2 = window.cache()
    assert cached == [True, True]
    assert fname != fname2


def test_FilterTabPage_obj_properties(qtbot):
    log2fc_col = 'my log2fc col'
    padj_col = 'my padj col'
    qtbot, window = widget_setup(qtbot, FilterTabPage)
    filt = filtering.DESeqFilter('tests/test_files/test_deseq.csv', log2fc_col=log2fc_col, padj_col=padj_col)
    window.start_from_filter_obj(filt, 'table name')

    assert window.obj_properties() == {'log2fc_col': log2fc_col, 'padj_col': padj_col}


def test_FilterTabPage_rename(qtbot, filtertabpage_with_undo_stack):
    new_name = 'my new table name'
    window, stack = filtertabpage_with_undo_stack
    qtbot.keyClicks(window.overview_widgets['table_name'], new_name)
    with qtbot.waitSignal(window.tabNameChange) as blocker:
        qtbot.mouseClick(window.overview_widgets['rename_button'], LEFT_CLICK)
    assert blocker.args[0] == new_name
    assert str(window.obj().fname.stem) == new_name
    assert new_name in window.overview_widgets['table_name_label'].text()
    assert 'test_deseq' not in window.overview_widgets['table_name_label'].text()
    assert window.name == new_name


def test_FilterTabPage_undo_rename(qtbot, filtertabpage_with_undo_stack):
    new_name = 'my new table name'
    window, stack = filtertabpage_with_undo_stack
    prev_name = window.name
    qtbot.keyClicks(window.overview_widgets['table_name'], new_name)
    with qtbot.waitSignal(window.tabNameChange) as blocker:
        qtbot.mouseClick(window.overview_widgets['rename_button'], LEFT_CLICK)
    assert blocker.args[0] == new_name

    with qtbot.waitSignal(window.tabNameChange) as blocker:
        stack.undo()
    assert blocker.args[0] == prev_name
    assert str(window.obj().fname.stem) == prev_name
    assert prev_name in window.overview_widgets['table_name_label'].text()
    assert new_name not in window.overview_widgets['table_name_label'].text()
    assert window.name == prev_name

    with qtbot.waitSignal(window.tabNameChange) as blocker:
        stack.redo()
    assert blocker.args[0] == new_name
    assert str(window.obj().fname.stem) == new_name
    assert new_name in window.overview_widgets['table_name_label'].text()
    assert prev_name not in window.overview_widgets['table_name_label'].text()
    assert window.name == new_name


def test_FilterTabPage_save_table(qtbot, monkeypatch):
    fname = 'my filename.tsv'
    saved = []

    def mock_get_save_name(*args, **kwargs):
        saved.append('got name')
        return fname, '.tsv'

    def mock_save_csv(self, filename):
        saved.append(filename)

    monkeypatch.setattr(QtWidgets.QFileDialog, 'getSaveFileName', mock_get_save_name)
    monkeypatch.setattr(filtering.Filter, 'save_csv', mock_save_csv)

    qtbot, window = widget_setup(qtbot, FilterTabPage)
    filt = filtering.DESeqFilter('tests/test_files/test_deseq.csv')
    window.start_from_filter_obj(filt, 'table name')
    qtbot.mouseClick(window.overview_widgets['save_button'], LEFT_CLICK)

    assert saved == ['got name', fname]


def test_FilterTabPage_view_full_table(qtbot, filtertabpage):
    assert False


def test_FilterTabPage_apply_function(qtbot, filtertabpage_with_undo_stack):
    window, stack = filtertabpage_with_undo_stack
    assert False


def test_FilterTabPage_apply_function_inplace(qtbot, filtertabpage_with_undo_stack):
    window, stack = filtertabpage_with_undo_stack
    truth = window.obj().filter_significant(0.01, opposite=True, inplace=False)
    window.stack_buttons[0].click()
    qtbot.keyClicks(window.stack.currentWidget().func_combo, 'filter_significant')
    window.stack.currentWidget().parameter_widgets['alpha'].clear()
    qtbot.keyClicks(window.stack.currentWidget().parameter_widgets['alpha'], '0.01')
    qtbot.mouseClick(window.stack.currentWidget().parameter_widgets['opposite'].switch, LEFT_CLICK)
    qtbot.mouseClick(window.basic_widgets['apply_button'], LEFT_CLICK)
    assert window.obj() == truth


def test_FilterTabPage_undo_function(qtbot, filtertabpage_with_undo_stack):
    window, stack = filtertabpage_with_undo_stack
    truth = window.obj().filter_significant(0.01, opposite=True, inplace=False)
    orig = window.obj().__copy__()
    window.stack_buttons[0].click()
    qtbot.keyClicks(window.stack.currentWidget().func_combo, 'filter_significant')
    window.stack.currentWidget().parameter_widgets['alpha'].clear()
    qtbot.keyClicks(window.stack.currentWidget().parameter_widgets['alpha'], '0.01')
    qtbot.mouseClick(window.stack.currentWidget().parameter_widgets['opposite'].switch, LEFT_CLICK)
    qtbot.mouseClick(window.basic_widgets['apply_button'], LEFT_CLICK)
    assert window.obj() == truth

    stack.undo()

    assert window.obj() == orig

    stack.redo()

    assert window.obj() == truth


def test_FilterTabPage_apply_pipeline(qtbot, filtertabpage_with_undo_stack, pipeline):
    window, stack = filtertabpage_with_undo_stack
    filter_obj_orig = window.obj().__copy__()

    filter_obj_truth = pipeline.apply_to(window.obj(), inplace=False)[0]

    with qtbot.waitSignal(window.filterObjectCreated) as blocker:
        window.apply_pipeline(pipeline, pipeline_name='my pipeline', inplace=False)
    assert blocker.args[0] == filter_obj_truth
    assert window.obj() == filter_obj_orig

    window.apply_pipeline(pipeline, pipeline_name='my pipeline', inplace=True)
    assert window.obj() == filter_obj_truth


def test_FilterTabPage_undo_pipeline(qtbot, filtertabpage_with_undo_stack, pipeline):
    window, stack = filtertabpage_with_undo_stack
    orig_name = window.name
    filter_obj_orig = window.obj().__copy__()
    filter_obj_truth = pipeline.apply_to(window.obj(), inplace=False)[0]
    window.apply_pipeline(pipeline, pipeline_name='my pipeline', inplace=True)
    assert window.obj() == filter_obj_truth

    stack.undo()
    assert window.obj() == filter_obj_orig
    assert window.name == orig_name

    stack.redo()
    assert window.obj() == filter_obj_truth
    assert window.name != orig_name


def test_FilterTabPage_get_all_actions(qtbot):
    assert False


def test_SetTabPage_init(qtbot):
    _, _ = widget_setup(qtbot, SetTabPage, 'set name')
    _, _ = widget_setup(qtbot, SetTabPage, 'set name', {'aa', 'bb', 'cc', 'dd'})


def test_SetTabPage_from_set(qtbot):
    set_name = 'table name'
    qtbot, window = widget_setup(qtbot, SetTabPage, set_name)
    assert window.is_empty()

    obj = {'abc', 'def', 'ghi', 'jkl'}
    window.update_gene_set(obj)

    assert not window.is_empty()
    assert window.obj() == obj
    assert window.obj_type() == set
    assert window.name == set_name


def test_SetTabPage_cache(qtbot, monkeypatch):
    s = {'abc', 'def', 'ghi', '123'}
    qtbot, window = widget_setup(qtbot, SetTabPage, 'set name', s)
    cached = []

    def mock_cache(obj, filename):
        assert isinstance(obj, set)
        assert obj == s
        cached.append(True)

    monkeypatch.setattr(io, 'cache_gui_file', mock_cache)

    fname = window.cache()
    assert fname.endswith('.txt')
    assert len(fname) == 44

    time.sleep(0.01)
    fname2 = window.cache()
    assert cached == [True, True]
    assert fname != fname2


def test_SetTabPage_obj_properties(qtbot):
    s = {'abc', 'def', 'ghi', '123'}
    qtbot, window = widget_setup(qtbot, SetTabPage, 'set name', s)
    assert window.obj_properties() == {}


def test_SetTabPage_rename(qtbot, settabpage_with_undo_stack):
    window, stack = settabpage_with_undo_stack
    new_name = 'my new set name'
    prev_name = window.name
    qtbot.keyClicks(window.overview_widgets['table_name'], new_name)
    with qtbot.waitSignal(window.tabNameChange) as blocker:
        qtbot.mouseClick(window.overview_widgets['rename_button'], LEFT_CLICK)
    assert blocker.args[0] == new_name
    assert str(window.gene_set.set_name) == new_name
    assert new_name in window.overview_widgets['table_name_label'].text()
    assert prev_name not in window.overview_widgets['table_name_label'].text()
    assert window.name == new_name


def test_SetTabPage_undo_rename(qtbot, settabpage_with_undo_stack):
    window, stack = settabpage_with_undo_stack
    new_name = 'my new set name'
    prev_name = window.name
    qtbot.keyClicks(window.overview_widgets['table_name'], new_name)
    with qtbot.waitSignal(window.tabNameChange) as blocker:
        qtbot.mouseClick(window.overview_widgets['rename_button'], LEFT_CLICK)
    assert window.name == new_name

    with qtbot.waitSignal(window.tabNameChange) as blocker:
        stack.undo()
    assert blocker.args[0] == prev_name
    assert str(window.gene_set.set_name) == prev_name
    assert prev_name in window.overview_widgets['table_name_label'].text()
    assert new_name not in window.overview_widgets['table_name_label'].text()
    assert window.name == prev_name

    with qtbot.waitSignal(window.tabNameChange) as blocker:
        stack.redo()
    assert blocker.args[0] == new_name
    assert str(window.gene_set.set_name) == new_name
    assert new_name in window.overview_widgets['table_name_label'].text()
    assert prev_name not in window.overview_widgets['table_name_label'].text()
    assert window.name == new_name


def test_SetTabPage_save_gene_set(qtbot, monkeypatch):
    fname = 'my filename.txt'
    saved = []

    def mock_get_save_name(*args, **kwargs):
        saved.append('got name')
        return fname, '.txt'

    def mock_save_txt(self, filename):
        saved.append(filename)

    monkeypatch.setattr(QtWidgets.QFileDialog, 'getSaveFileName', mock_get_save_name)
    monkeypatch.setattr(enrichment.FeatureSet, 'save_txt', mock_save_txt)

    qtbot, window = widget_setup(qtbot, SetTabPage, 'set name', {'1', '2', '3'})
    qtbot.mouseClick(window.overview_widgets['save_button'], LEFT_CLICK)

    assert saved == ['got name', fname]


def test_SetTabPage_view_full_set(qtbot):
    assert False


def test_FuncTypeStack_init(qtbot):
    assert False


def test_CreatePipelineWindow_init(qtbot):
    _, _ = widget_setup(qtbot, CreatePipelineWindow)


def test_CreatePipelineWindow_create_pipeline(qtbot, monkeypatch):
    monkeypatch.setattr(QtWidgets.QMessageBox, "question", lambda *args: QtWidgets.QMessageBox.Yes)
    pipeline_name = 'my pipeline name'
    qtbot, window = widget_setup(qtbot, CreatePipelineWindow)
    window.basic_widgets['pipeline_name'].clear()
    qtbot.keyClicks(window.basic_widgets['pipeline_name'], pipeline_name)
    qtbot.keyClicks(window.basic_widgets['table_type_combo'], 'Differential expression')
    qtbot.mouseClick(window.basic_widgets['start_button'], LEFT_CLICK)

    assert not window.basic_group.isVisible()
    assert isinstance(window.pipeline, filtering.Pipeline)
    assert window.pipeline.filter_type == filtering.DESeqFilter
    assert window._get_pipeline_name() == pipeline_name


def test_CreatePipelineWindow_add_function(qtbot, monkeypatch):
    pipeline_truth = filtering.Pipeline('DESeqFilter')
    pipeline_truth.add_function('split_fold_change_direction')

    monkeypatch.setattr(QtWidgets.QMessageBox, "question", lambda *args: QtWidgets.QMessageBox.Yes)

    qtbot, window = widget_setup(qtbot, CreatePipelineWindow)
    window.basic_widgets['pipeline_name'].clear()
    qtbot.keyClicks(window.basic_widgets['pipeline_name'], 'pipeline_name')
    qtbot.keyClicks(window.basic_widgets['table_type_combo'], 'Differential expression')
    qtbot.mouseClick(window.basic_widgets['start_button'], LEFT_CLICK)

    qtbot.mouseClick(window.stack_buttons[0], LEFT_CLICK)
    qtbot.keyClicks(window.stack.currentWidget().func_combo, 'split_fold_change_direction')
    qtbot.mouseClick(window.basic_widgets['apply_button'], LEFT_CLICK)

    assert window.pipeline == pipeline_truth


def test_CreatePipelineWindow_add_function_with_args(qtbot, monkeypatch):
    pipeline_truth = filtering.Pipeline('DESeqFilter')
    pipeline_truth.add_function('filter_significant', alpha=0.01, opposite=True)

    monkeypatch.setattr(QtWidgets.QMessageBox, "question", lambda *args: QtWidgets.QMessageBox.Yes)

    qtbot, window = widget_setup(qtbot, CreatePipelineWindow)
    window.basic_widgets['pipeline_name'].clear()
    qtbot.keyClicks(window.basic_widgets['pipeline_name'], 'pipeline_name')
    qtbot.keyClicks(window.basic_widgets['table_type_combo'], 'Differential expression')
    qtbot.mouseClick(window.basic_widgets['start_button'], LEFT_CLICK)

    qtbot.mouseClick(window.stack_buttons[0], LEFT_CLICK)
    qtbot.keyClicks(window.stack.currentWidget().func_combo, 'filter_significant')
    window.stack.currentWidget().parameter_widgets['alpha'].clear()
    qtbot.keyClicks(window.stack.currentWidget().parameter_widgets['alpha'], '0.01')
    qtbot.mouseClick(window.stack.currentWidget().parameter_widgets['opposite'].switch, LEFT_CLICK)
    qtbot.mouseClick(window.basic_widgets['apply_button'], LEFT_CLICK)
    assert window.pipeline == pipeline_truth

    # add a second function
    pipeline_truth.add_function('split_fold_change_direction')

    window.stack.currentWidget().func_combo.setCurrentText('split_fold_change_direction')
    qtbot.mouseClick(window.basic_widgets['apply_button'], LEFT_CLICK)
    assert window.pipeline == pipeline_truth


def test_CreatePipelineWindow_save_pipeline(qtbot, monkeypatch):
    pipeline_truth = filtering.Pipeline('DESeqFilter')
    pipeline_truth.add_function('describe', percentiles=[0.01, 0.25, 0.5, 0.75, 0.99])
    pipeline_name = 'my pipeline name'

    monkeypatch.setattr(QtWidgets.QMessageBox, "question", lambda *args: QtWidgets.QMessageBox.Yes)

    qtbot, window = widget_setup(qtbot, CreatePipelineWindow)
    window.basic_widgets['pipeline_name'].clear()
    qtbot.keyClicks(window.basic_widgets['pipeline_name'], pipeline_name)
    qtbot.keyClicks(window.basic_widgets['table_type_combo'], 'Differential expression')
    qtbot.mouseClick(window.basic_widgets['start_button'], LEFT_CLICK)

    qtbot.mouseClick(window.stack_buttons[2], LEFT_CLICK)
    qtbot.keyClicks(window.stack.currentWidget().func_combo, 'describe')
    qtbot.mouseClick(window.basic_widgets['apply_button'], LEFT_CLICK)

    with qtbot.waitSignal(window.pipelineSaved) as blocker:
        qtbot.mouseClick(window.overview_widgets['save_button'], LEFT_CLICK)
    assert blocker.args[0] == pipeline_name
    assert blocker.args[1] == pipeline_truth


def test_CreatePipelineWindow_export_pipeline(qtbot, monkeypatch):
    monkeypatch.setattr(QtWidgets.QMessageBox, "question", lambda *args: QtWidgets.QMessageBox.Yes)
    qtbot, window = widget_setup(qtbot, CreatePipelineWindow)
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
