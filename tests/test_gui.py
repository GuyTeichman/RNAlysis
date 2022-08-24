import pytest
import matplotlib

matplotlib.use('Agg')
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
def available_objects(qtbot, red_icon, green_icon):
    qtbot, first = widget_setup(qtbot, SetTabPage, 'first tab',
                                {'WBGene00000002', 'WBGene00000006', 'WBGene00000015', 'WBGene00000017'})

    qtbot, second = widget_setup(qtbot, FilterTabPage, undo_stack=QtWidgets.QUndoStack())
    second.start_from_filter_obj(filtering.DESeqFilter('tests/test_files/test_deseq.csv'))
    second.rename('second tab')

    qtbot, third = widget_setup(qtbot, FilterTabPage, undo_stack=QtWidgets.QUndoStack())
    third.start_from_filter_obj(filtering.CountFilter('tests/test_files/counted.tsv'))
    third.rename('third tab')

    return {'first tab': (first, red_icon), 'second tab': (second, red_icon), 'third tab': (third, green_icon)}


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
def countfiltertabpage_with_undo_stack(qtbot):
    stack = QtWidgets.QUndoStack()
    qtbot, window = widget_setup(qtbot, FilterTabPage, undo_stack=stack)
    window.start_from_filter_obj(filtering.CountFilter('tests/test_files/counted.csv'))
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


@pytest.fixture
def clicom_window(qtbot):
    funcs = {'split_kmeans': 'K-Means', 'split_kmedoids': 'K-Medoids',
             'split_hierarchical': 'Hierarchical (Agglomerative)', 'split_hdbscan': 'HDBSCAN'}
    qtbot, window = widget_setup(qtbot, ClicomWindow, funcs, filtering.CountFilter('tests/test_files/counted.csv'))
    return window


@pytest.fixture
def enrichment_window(qtbot, available_objects):
    qtbot, window = widget_setup(qtbot, EnrichmentWindow, available_objects)
    return window


@pytest.fixture
def set_op_window(qtbot, available_objects):
    qtbot, window = widget_setup(qtbot, SetOperationWindow, available_objects)
    return window


@pytest.fixture
def set_vis_window(qtbot, available_objects):
    qtbot, window = widget_setup(qtbot, SetVisualizationWindow, available_objects)
    return window


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


def test_ClicomWindow_init(qtbot, clicom_window):
    _ = clicom_window


def test_ClicomWindow_add_setup(qtbot, clicom_window):
    truth = dict(method='kmeans', n_clusters=3, n_init=3, max_iter=300, random_seed=None,
                 max_n_clusters_estimate='auto')

    qtbot.keyClicks(clicom_window.stack.func_combo, 'split_kmeans')
    clicom_window.stack.parameter_widgets['n_clusters'].other.set_defaults(3)
    qtbot.mouseClick(clicom_window.setups_widgets['add_button'], LEFT_CLICK)
    assert len(clicom_window.parameter_dicts) == 1
    assert clicom_window.parameter_dicts[0] == truth


def test_ClicomWindow_remove_setup(qtbot, monkeypatch, clicom_window):
    monkeypatch.setattr(QtWidgets.QMessageBox, 'question', lambda *args: QtWidgets.QMessageBox.Yes)
    qtbot.keyClicks(clicom_window.stack.func_combo, 'split_kmeans')
    clicom_window.stack.parameter_widgets['n_clusters'].other.set_defaults(3)
    qtbot.mouseClick(clicom_window.setups_widgets['add_button'], LEFT_CLICK)
    assert len(clicom_window.parameter_dicts) == 1

    qtbot.mouseClick(clicom_window.setups_widgets['list'].delete_all_button, LEFT_CLICK)

    assert len(clicom_window.parameter_dicts) == 0


def test_ClicomWindow_get_analysis_params(qtbot, clicom_window):
    truth = dict(power_transform=[True, False], evidence_threshold=0.35, cluster_unclustered_features=True,
                 min_cluster_size=15, plot_style='all', split_plots=False)

    qtbot.mouseClick(clicom_window.param_widgets['power_transform'].false_button, LEFT_CLICK)
    qtbot.mouseClick(clicom_window.param_widgets['cluster_unclustered_features'].switch, LEFT_CLICK)
    clicom_window.param_widgets['evidence_threshold'].clear()
    qtbot.keyClicks(clicom_window.param_widgets['evidence_threshold'], '0.35')

    assert clicom_window.get_analysis_params() == truth


def test_ClicomWindow_start_clustering(qtbot, clicom_window):
    truth_setups = [dict(method='kmeans', n_clusters=3, n_init=3, max_iter=300, random_seed=None,
                         max_n_clusters_estimate='auto'),
                    dict(method='hierarchical', n_clusters='silhouette', metric='Euclidean', linkage='Average',
                         distance_threshold=None, max_n_clusters_estimate='auto')]
    truth_params = dict(power_transform=[True, False], evidence_threshold=0.35, cluster_unclustered_features=True,
                        min_cluster_size=15, plot_style='all', split_plots=False)

    qtbot.keyClicks(clicom_window.stack.func_combo, 'split_kmeans')
    clicom_window.stack.parameter_widgets['n_clusters'].other.set_defaults(3)
    qtbot.mouseClick(clicom_window.setups_widgets['add_button'], LEFT_CLICK)

    clicom_window.stack.func_combo.setCurrentText('split_hierarchical')
    qtbot.keyClicks(clicom_window.stack.parameter_widgets['n_clusters'].combo, 'silhouette')
    qtbot.mouseClick(clicom_window.setups_widgets['add_button'], LEFT_CLICK)

    qtbot.mouseClick(clicom_window.param_widgets['power_transform'].false_button, LEFT_CLICK)
    qtbot.mouseClick(clicom_window.param_widgets['cluster_unclustered_features'].switch, LEFT_CLICK)
    clicom_window.param_widgets['evidence_threshold'].clear()
    qtbot.keyClicks(clicom_window.param_widgets['evidence_threshold'], '0.35')

    with qtbot.waitSignal(clicom_window.paramsAccepted) as blocker:
        qtbot.mouseClick(clicom_window.start_button, LEFT_CLICK)
    assert blocker.args[0] == truth_setups
    assert blocker.args[1] == truth_params


def test_EnrichmentWindow_init(qtbot, enrichment_window):
    _ = enrichment_window


@pytest.mark.parametrize('button_name,truth', [
    ('Gene Ontology (GO)', 'go'),
    ('Kyoto Encyclopedia of Genes and Genomes (KEGG)', 'kegg'),
    ('Categorical attributes', 'user_defined'),
    ('Non-categorical attributes', 'non_categorical')
])
def test_EnrichmentWindow_get_analysis_type(qtbot, enrichment_window, button_name, truth):
    enrichment_window.widgets['dataset_radiobox'].radio_buttons[button_name].click()
    assert enrichment_window.get_current_analysis_type() == truth


@pytest.mark.parametrize('button_name', [
    'Gene Ontology (GO)',
    'Kyoto Encyclopedia of Genes and Genomes (KEGG)',
    'Categorical attributes',
])
@pytest.mark.parametrize('test_name,truth', [
    ("Fisher's Exact test", False),
    ('Hypergeometric test', False),
    ('Randomization test', False),
    ('Single-set enrichment (XL-mHG test)', True)
])
def test_EnrichmentWindow_is_single_set(qtbot, enrichment_window, button_name, test_name, truth):
    enrichment_window.widgets['dataset_radiobox'].radio_buttons[button_name].click()
    enrichment_window.stats_widgets['stats_radiobox'].radio_buttons[test_name].click()

    assert enrichment_window.is_single_set() == truth


@pytest.mark.parametrize('test_name,truth', [
    ("One-sample T-test (parametric)", False),
    ('Sign test (non-parametric)', False)])
def test_EnrichmentWindow_is_single_set_non_categorical(qtbot, enrichment_window, test_name, truth):
    enrichment_window.widgets['dataset_radiobox'].radio_buttons['Non-categorical attributes'].click()
    enrichment_window.stats_widgets['stats_radiobox'].radio_buttons[test_name].click()

    assert enrichment_window.is_single_set() == truth


@pytest.mark.parametrize('button_name,test_name,func_truth', [
    ('Gene Ontology (GO)', "Fisher's Exact test", enrichment.FeatureSet.go_enrichment),
    ('Gene Ontology (GO)', 'Hypergeometric test', enrichment.FeatureSet.go_enrichment),
    ('Gene Ontology (GO)', 'Randomization test', enrichment.FeatureSet.go_enrichment),
    ('Gene Ontology (GO)', 'Single-set enrichment (XL-mHG test)', enrichment.RankedSet.single_set_go_enrichment),
    ('Kyoto Encyclopedia of Genes and Genomes (KEGG)', "Fisher's Exact test", enrichment.FeatureSet.kegg_enrichment),
    ('Kyoto Encyclopedia of Genes and Genomes (KEGG)', 'Hypergeometric test', enrichment.FeatureSet.kegg_enrichment),
    ('Kyoto Encyclopedia of Genes and Genomes (KEGG)', 'Randomization test', enrichment.FeatureSet.kegg_enrichment),
    ('Kyoto Encyclopedia of Genes and Genomes (KEGG)', 'Single-set enrichment (XL-mHG test)',
     enrichment.RankedSet.single_set_kegg_enrichment),
    ('Categorical attributes', "Fisher's Exact test", enrichment.FeatureSet.user_defined_enrichment),
    ('Categorical attributes', 'Hypergeometric test', enrichment.FeatureSet.user_defined_enrichment),
    ('Categorical attributes', 'Randomization test', enrichment.FeatureSet.user_defined_enrichment),
    ('Categorical attributes', 'Single-set enrichment (XL-mHG test)', enrichment.RankedSet.single_set_enrichment),
    ('Non-categorical attributes', "One-sample T-test (parametric)", enrichment.FeatureSet.non_categorical_enrichment),
    ('Non-categorical attributes', "Sign test (non-parametric)", enrichment.FeatureSet.non_categorical_enrichment)
])
def test_EnrichmentWindow_get_func(qtbot, enrichment_window, button_name, test_name, func_truth):
    enrichment_window.widgets['dataset_radiobox'].radio_buttons[button_name].click()
    enrichment_window.stats_widgets['stats_radiobox'].radio_buttons[test_name].click()

    assert enrichment_window.get_current_func() == func_truth


@pytest.mark.parametrize('button_name,truth', [
    ('Gene Ontology (GO)', True),
    ('Kyoto Encyclopedia of Genes and Genomes (KEGG)', True),
    ('Categorical attributes', True),
    ('Non-categorical attributes', False)
])
def test_EnrichmentWindow_is_categorical(qtbot, enrichment_window, button_name, truth):
    enrichment_window.widgets['dataset_radiobox'].radio_buttons[button_name].click()
    assert enrichment_window.is_categorical() == truth


@pytest.mark.parametrize('en_set,bg_set,en_set_truth,bg_set_truth,', [
    ('first tab', 'second tab', 'first tab', 'second tab'),
    ('third tab', 'first tab', 'third tab', 'first tab'),
    ('second tab', 'third tab', 'second tab', 'third tab'),
])
@pytest.mark.parametrize('button_name,dataset_kwargs', [
    ('Gene Ontology (GO)',
     dict(plot_horizontal=True, plot_ontology_graph=False, organism='auto', excluded_evidence_types='experimental')),
    ('Kyoto Encyclopedia of Genes and Genomes (KEGG)',
     dict(plot_horizontal=True, plot_pathway_graphs=True, gene_id_type='auto')),
    ('Categorical attributes', dict(attributes='all', plot_horizontal=False))
])
@pytest.mark.parametrize('test_name,is_single_set,test_arg_truth,stats_kwargs', [
    ("Fisher's Exact test", False, 'fisher', dict(alpha=0.05)),
    ('Hypergeometric test', False, 'hypergeometric', dict(alpha=0.5)),
    ('Randomization test', False, 'randomization', dict(alpha=0.13, random_seed=42)),
    ('Single-set enrichment (XL-mHG test)', True, 'single_set', dict(alpha=0.01))
])
def test_EnrichmentWindow_get_analysis_params(qtbot, enrichment_window, button_name, test_name, test_arg_truth, en_set,
                                              en_set_truth, bg_set, bg_set_truth, is_single_set, available_objects,
                                              stats_kwargs, dataset_kwargs):
    kwargs_truth = dict()
    kwargs_truth.update(stats_kwargs)
    kwargs_truth.update(dataset_kwargs)

    set_name_truth = available_objects[en_set_truth][0].name
    enrichment_window.widgets['dataset_radiobox'].radio_buttons[button_name].click()
    enrichment_window.stats_widgets['stats_radiobox'].radio_buttons[test_name].click()

    qtbot.keyClicks(enrichment_window.widgets['enrichment_list'], en_set)

    enrichment_window.stats_widgets['alpha'].clear()
    qtbot.keyClicks(enrichment_window.stats_widgets['alpha'], str(kwargs_truth['alpha']))

    for key in stats_kwargs:
        if key not in {'alpha'}:
            enrichment_window.stats_widgets[key].setValue(stats_kwargs[key])

    for key in dataset_kwargs:
        if key in enrichment_window.parameter_widgets:
            qtbot.keyClicks(enrichment_window.parameter_widgets[key].combo, dataset_kwargs[key])
        elif key in enrichment_window.plot_widgets:
            if not dataset_kwargs[key]:
                enrichment_window.plot_widgets[key].switch.click()

    if not is_single_set:
        qtbot.keyClicks(enrichment_window.widgets['bg_list'], bg_set)

    gene_set, bg_set, gene_set_name, kwargs = enrichment_window.get_analysis_params()

    assert gene_set == available_objects[en_set_truth][0].obj()
    if is_single_set:
        assert bg_set is None
        assert 'statistical_test' not in kwargs
    else:
        assert bg_set == available_objects[bg_set_truth][0].obj()
        assert kwargs['statistical_test'] == test_arg_truth

    for key in kwargs_truth:
        assert kwargs[key] == kwargs_truth[key]

    assert gene_set_name == set_name_truth


@pytest.mark.parametrize('en_set,bg_set,en_set_truth,bg_set_truth,', [
    ('first tab', 'second tab', 'first tab', 'second tab'),
    ('third tab', 'first tab', 'third tab', 'first tab'),
    ('second tab', 'third tab', 'second tab', 'third tab'),
])
@pytest.mark.parametrize('button_name,dataset_kwargs', [
    ('Non-categorical attributes', dict(plot_log_scale=False, attributes='all')),
])
@pytest.mark.parametrize('test_name,test_arg_truth,stats_kwargs', [
    ("One-sample T-test (parametric)", True, dict(alpha=0.08)),
    ('Sign test (non-parametric)', False, dict(alpha=0.5))
])
def test_EnrichmentWindow_get_analysis_params_single_set(qtbot, enrichment_window, button_name, test_name,
                                                         test_arg_truth, en_set, en_set_truth, bg_set, bg_set_truth,
                                                         available_objects, stats_kwargs, dataset_kwargs):
    kwargs_truth = dict()
    kwargs_truth.update(stats_kwargs)
    kwargs_truth.update(dataset_kwargs)

    set_name_truth = available_objects[en_set_truth][0].name
    enrichment_window.widgets['dataset_radiobox'].radio_buttons[button_name].click()
    enrichment_window.stats_widgets['stats_radiobox'].radio_buttons[test_name].click()

    qtbot.keyClicks(enrichment_window.widgets['enrichment_list'], en_set)

    enrichment_window.stats_widgets['alpha'].clear()
    qtbot.keyClicks(enrichment_window.stats_widgets['alpha'], str(kwargs_truth['alpha']))

    for key in stats_kwargs:
        if key not in {'alpha'}:
            enrichment_window.stats_widgets[key].setValue(stats_kwargs[key])

    for key in dataset_kwargs:
        if key in enrichment_window.parameter_widgets:
            qtbot.keyClicks(enrichment_window.parameter_widgets[key].combo, dataset_kwargs[key])
        elif key in enrichment_window.plot_widgets:
            if not dataset_kwargs[key]:
                enrichment_window.plot_widgets[key].switch.click()

    qtbot.keyClicks(enrichment_window.widgets['bg_list'], bg_set)

    gene_set, bg_set, gene_set_name, kwargs = enrichment_window.get_analysis_params()

    assert gene_set == available_objects[en_set_truth][0].obj()

    assert bg_set == available_objects[bg_set_truth][0].obj()
    assert kwargs['parametric_test'] == test_arg_truth

    assert gene_set_name == set_name_truth


@pytest.mark.parametrize('en_set,bg_set,en_set_truth,bg_set_truth,', [
    ('third tab', 'first tab', 'third tab', 'first tab'),
    ('second tab', 'third tab', 'second tab', 'third tab'),
])
@pytest.mark.parametrize('button_name,dataset_name_truth,dataset_kwargs', [
    ('Gene Ontology (GO)', 'go',
     dict(plot_horizontal=True, plot_ontology_graph=False, organism='auto', excluded_evidence_types='experimental')),
    ('Kyoto Encyclopedia of Genes and Genomes (KEGG)', 'kegg',
     dict(plot_horizontal=True, plot_pathway_graphs=True, gene_id_type='auto')),
    ('Categorical attributes', 'user_defined', dict(attributes='all', plot_horizontal=False))
])
@pytest.mark.parametrize('test_name,is_single_set,test_arg_truth,stats_kwargs', [
    ("Fisher's Exact test", False, 'fisher', dict(alpha=0.05)),
    ('Hypergeometric test', False, 'hypergeometric', dict(alpha=0.5)),
    ('Randomization test', False, 'randomization', dict(alpha=0.13, random_seed=42)),
    ('Single-set enrichment (XL-mHG test)', True, 'single_set', dict(alpha=0.01))
])
def test_EnrichmentWindow_run_analysis(qtbot, enrichment_window, button_name, test_name,
                                       test_arg_truth, en_set, en_set_truth, bg_set, bg_set_truth, dataset_name_truth,
                                       available_objects, stats_kwargs, dataset_kwargs, is_single_set):
    func_truth = {
        ('go', False): enrichment.FeatureSet.go_enrichment,
        ('go', True): enrichment.RankedSet.single_set_go_enrichment,
        ('kegg', False): enrichment.FeatureSet.kegg_enrichment,
        ('kegg', True): enrichment.RankedSet.single_set_kegg_enrichment,
        ('user_defined', False): enrichment.FeatureSet.user_defined_enrichment,
        ('user_defined', True): enrichment.RankedSet.single_set_enrichment}

    set_name_truth = available_objects[en_set_truth][0].name

    if is_single_set:
        gene_set_truth = enrichment.RankedSet(available_objects[en_set_truth][0].obj(), set_name_truth)
    else:
        gene_set_truth = enrichment.FeatureSet(available_objects[en_set_truth][0].obj() if isinstance(
            available_objects[en_set_truth][0].obj(), set) else available_objects[en_set_truth][0].obj().index_set,
                                               set_name_truth)

    kwargs_truth = dict()
    kwargs_truth.update(stats_kwargs)
    kwargs_truth.update(dataset_kwargs)

    enrichment_window.widgets['dataset_radiobox'].radio_buttons[button_name].click()
    enrichment_window.stats_widgets['stats_radiobox'].radio_buttons[test_name].click()

    qtbot.keyClicks(enrichment_window.widgets['enrichment_list'], en_set)

    enrichment_window.stats_widgets['alpha'].clear()
    qtbot.keyClicks(enrichment_window.stats_widgets['alpha'], str(kwargs_truth['alpha']))

    for key in stats_kwargs:
        if key not in {'alpha'}:
            enrichment_window.stats_widgets[key].setValue(stats_kwargs[key])

    for key in dataset_kwargs:
        if key in enrichment_window.parameter_widgets:
            qtbot.keyClicks(enrichment_window.parameter_widgets[key].combo, dataset_kwargs[key])
        elif key in enrichment_window.plot_widgets:
            if not dataset_kwargs[key]:
                enrichment_window.plot_widgets[key].switch.click()

    if not is_single_set:
        qtbot.keyClicks(enrichment_window.widgets['bg_list'], bg_set)

    with qtbot.waitSignal(enrichment_window.enrichmentStarted) as blocker:
        enrichment_window.widgets['run_button'].click()

    assert callable(blocker.args[0])
    assert blocker.args[1] == set_name_truth

    if is_single_set:
        pass
    else:
        assert blocker.args[0].keywords['statistical_test'] == test_arg_truth
        assert blocker.args[0].keywords['background_genes'].gene_set == available_objects[bg_set_truth][
            0].obj() if isinstance(available_objects[bg_set_truth][0].obj(), set) else available_objects[bg_set_truth][
            0].obj().index_set

    for kw in kwargs_truth:
        assert blocker.args[0].keywords[kw] == kwargs_truth[kw]

    assert blocker.args[0].func == func_truth[(dataset_name_truth, is_single_set)]
    assert blocker.args[0].args[0] == gene_set_truth


@pytest.mark.parametrize('en_set,bg_set,en_set_truth,bg_set_truth,', [
    ('first tab', 'second tab', 'first tab', 'second tab'),
    ('third tab', 'first tab', 'third tab', 'first tab'),
    ('second tab', 'third tab', 'second tab', 'third tab'),
])
@pytest.mark.parametrize('button_name,dataset_kwargs', [
    ('Non-categorical attributes', dict(plot_log_scale=False, attributes='all')),
])
@pytest.mark.parametrize('test_name,test_arg_truth,stats_kwargs', [
    ("One-sample T-test (parametric)", True, dict(alpha=0.08)),
    ('Sign test (non-parametric)', False, dict(alpha=0.5))
])
def test_EnrichmentWindow_run_analysis_non_categorical(qtbot, enrichment_window, button_name, test_name,
                                                       test_arg_truth, en_set, en_set_truth, bg_set, bg_set_truth,
                                                       available_objects, stats_kwargs, dataset_kwargs):
    set_name_truth = available_objects[en_set_truth][0].name
    gene_set_truth = enrichment.FeatureSet(available_objects[en_set_truth][0].obj() if isinstance(
        available_objects[en_set_truth][0].obj(), set) else available_objects[en_set_truth][0].obj().index_set,
                                           set_name_truth)

    kwargs_truth = dict()
    kwargs_truth.update(stats_kwargs)
    kwargs_truth.update(dataset_kwargs)

    enrichment_window.widgets['dataset_radiobox'].radio_buttons[button_name].click()
    enrichment_window.stats_widgets['stats_radiobox'].radio_buttons[test_name].click()

    qtbot.keyClicks(enrichment_window.widgets['enrichment_list'], en_set)

    enrichment_window.stats_widgets['alpha'].clear()
    qtbot.keyClicks(enrichment_window.stats_widgets['alpha'], str(kwargs_truth['alpha']))

    for key in stats_kwargs:
        if key not in {'alpha'}:
            enrichment_window.stats_widgets[key].setValue(stats_kwargs[key])

    for key in dataset_kwargs:
        if key in enrichment_window.parameter_widgets:
            qtbot.keyClicks(enrichment_window.parameter_widgets[key].combo, dataset_kwargs[key])
        elif key in enrichment_window.plot_widgets:
            if not dataset_kwargs[key]:
                enrichment_window.plot_widgets[key].switch.click()

    qtbot.keyClicks(enrichment_window.widgets['bg_list'], bg_set)

    with qtbot.waitSignal(enrichment_window.enrichmentStarted) as blocker:
        enrichment_window.widgets['run_button'].click()

    assert callable(blocker.args[0])
    assert blocker.args[1] == set_name_truth

    assert blocker.args[0].keywords['parametric_test'] == test_arg_truth
    for kw in kwargs_truth:
        assert blocker.args[0].keywords[kw] == kwargs_truth[kw]

    assert blocker.args[0].func == enrichment.FeatureSet.non_categorical_enrichment
    assert blocker.args[0].args[0] == gene_set_truth

    assert blocker.args[0].keywords['background_genes'].gene_set == available_objects[bg_set_truth][
        0].obj() if isinstance(available_objects[bg_set_truth][0].obj(), set) else available_objects[bg_set_truth][
        0].obj().index_set


def test_SetOperationWindow_init(qtbot, set_op_window):
    _ = set_op_window


def test_SetOperationWindow_get_current_func_name(qtbot, set_op_window):
    assert False


def test_SetOperationWindow_canvas_types(qtbot, set_op_window):
    assert False


def test_SetOperationWindow_function_change_canvas(qtbot, set_op_window):
    assert False


def test_SetOperationWindow_primary_set_change(qtbot, set_op_window):
    assert False


def test_SetOperationWindow_apply_set_op(qtbot, set_op_window):
    assert False


def test_SetOperationWindow_apply_set_op_other(qtbot, set_op_window):
    assert False


def test_SetOperationWindow_apply_set_op_inplace(qtbot, set_op_window):
    assert False


def test_SetVisualizationWindow_init(qtbot, set_vis_window):
    _ = set_vis_window


def test_SetVisualizationWindow_get_current_func_name(qtbot, set_vis_window):
    assert False


def test_SetVisualizationWindow_canvas_types(qtbot, set_vis_window):
    assert False


def test_SetVisualizationWindow_function_change_canvas(qtbot, set_vis_window):
    assert False


def test_SetVisualizationWindow_parameter_change_canvas(qtbot, set_vis_window):
    assert False


def test_SetVisualizationWindow_generate_graph(qtbot, set_vis_window):
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


@pytest.mark.parametrize('filter_obj,truth', [
    (filtering.DESeqFilter('tests/test_files/test_deseq.csv', log2fc_col='my log2fc col', padj_col='my padj col'),
     {'log2fc_col': 'my log2fc col', 'padj_col': 'my padj col'}),
    (filtering.CountFilter('tests/test_files/counted.tsv'), {'is_normalized': False}),
    (filtering.FoldChangeFilter('tests/test_files/fc_1.csv', 'num_name', 'denom_name'),
     {'numerator_name': 'num_name', 'denominator_name': 'denom_name'}),
    (filtering.Filter('tests/test_files/test_deseq.csv'), {})

])
def test_FilterTabPage_obj_properties(qtbot, filter_obj, truth):
    qtbot, window = widget_setup(qtbot, FilterTabPage)
    window.start_from_filter_obj(filter_obj, 'table name')

    assert window.obj_properties() == truth


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
    qtbot.mouseClick(filtertabpage.overview_widgets['view_button'], LEFT_CLICK)
    assert isinstance(filtertabpage.overview_widgets['full_table_view'], gui_windows.DataFrameView)
    assert filtertabpage.overview_widgets['full_table_view'].data_view.model()._dataframe.equals(filtertabpage.obj().df)


def test_FilterTabPage_apply_function(qtbot, filtertabpage_with_undo_stack):
    window, stack = filtertabpage_with_undo_stack
    orig = window.obj().__copy__()
    truth = window.obj().filter_significant(0.01, opposite=True, inplace=False)
    window.stack_buttons[0].click()
    qtbot.keyClicks(window.stack.currentWidget().func_combo, 'filter_significant')
    window.stack.currentWidget().parameter_widgets['alpha'].clear()
    qtbot.keyClicks(window.stack.currentWidget().parameter_widgets['alpha'], '0.01')
    qtbot.mouseClick(window.stack.currentWidget().parameter_widgets['opposite'].switch, LEFT_CLICK)
    qtbot.mouseClick(window.stack.currentWidget().parameter_widgets['inplace'].switch, LEFT_CLICK)
    with qtbot.waitSignal(window.filterObjectCreated, timeout=10000) as blocker:
        qtbot.mouseClick(window.basic_widgets['apply_button'], LEFT_CLICK)
    assert blocker.args[0] == truth
    assert window.obj() == orig


def test_FilterTabPage_apply_split_clustering_function(qtbot, monkeypatch, countfiltertabpage_with_undo_stack):
    def mock_show_multikeep(self):
        self.select_all.setChecked(True)
        self.change_all()
        self.accept()
        self.button_box.button(QtWidgets.QDialogButtonBox.Ok).click()

    monkeypatch.setattr(MultiKeepWindow, 'exec', mock_show_multikeep)

    window, stack = countfiltertabpage_with_undo_stack

    def my_slot(partial, func_name):
        window._proccess_outputs(partial()[0], func_name)

    window.startedClustering.connect(my_slot)

    orig = window.obj().__copy__()
    orig_renamed = orig.__copy__()
    orig_renamed.fname = Path('counted')
    truth = orig_renamed.split_kmeans(n_clusters=3, random_seed=0)
    truth = sorted(truth, key=lambda obj: obj.shape[0])

    window.stack_buttons[4].click()
    qtbot.keyClicks(window.stack.currentWidget().func_combo, 'split_kmeans')
    window.stack.currentWidget().parameter_widgets['n_clusters'].other.set_defaults(3)
    qtbot.mouseClick(window.stack.currentWidget().parameter_widgets['random_seed'].checkbox, LEFT_CLICK)
    with qtbot.waitSignals([window.filterObjectCreated, window.filterObjectCreated, window.filterObjectCreated],
                           timeout=15000) as blocker:
        qtbot.mouseClick(window.basic_widgets['apply_button'], LEFT_CLICK)

    res = sorted([sig.args[0] for sig in blocker.all_signals_and_args], key=lambda obj: obj.shape[0])
    assert res == truth
    assert window.obj() == orig


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


def test_FilterTabPage_open_clicom(qtbot):
    assert False


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
    qtbot, window = widget_setup(qtbot, SetTabPage, 'set name', {'a', 'b', 'c', 'd'})
    qtbot.mouseClick(window.overview_widgets['view_button'], LEFT_CLICK)
    assert isinstance(window.overview_widgets['full_table_view'], gui_windows.GeneSetView)

    view = window.overview_widgets['full_table_view'].data_view
    genes_in_view = {view.item(i).text() for i in range(view.count())}
    assert genes_in_view == window.obj()


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


def test_MainMenu_add_new_tab(qtbot, use_temp_settings_file):
    assert False


def test_MainWindow_close_tab(qtbot, use_temp_settings_file):
    # qtbot, window = widget_setup(qtbot, MainWindow)
    assert False


def test_MainWindow_close_tab_undo(qtbot, use_temp_settings_file):
    # qtbot, window = widget_setup(qtbot, MainWindow)
    assert False


def test_MainWindow_sort_tabs_by_type(qtbot, use_temp_settings_file):
    # qtbot, window = widget_setup(qtbot, MainWindow)
    assert False


def test_MainWindow_sort_tabs_by_n_features(qtbot, use_temp_settings_file):
    # qtbot, window = widget_setup(qtbot, MainWindow)
    assert False


def test_MainWindow_sort_tabs_by_creation_time(qtbot, use_temp_settings_file):
    # qtbot, window = widget_setup(qtbot, MainWindow)
    assert False


def test_MainWindow_sort_tabs_by_name(qtbot, use_temp_settings_file):
    # qtbot, window = widget_setup(qtbot, MainWindow)
    assert False


def test_MainWindow_reverse_order(qtbot, use_temp_settings_file):
    # qtbot, window = widget_setup(qtbot, MainWindow)
    assert False


def test_MainWindow_change_tab_icon(qtbot, use_temp_settings_file):
    # qtbot, window = widget_setup(qtbot, MainWindow)
    assert False


def test_MainWindow_clear_history(qtbot, use_temp_settings_file):
    # qtbot, window = widget_setup(qtbot, MainWindow)
    assert False


def test_MainWindow_rename_tab(qtbot, use_temp_settings_file):
    # qtbot, window = widget_setup(qtbot, MainWindow)
    assert False


def test_MainWindow_new_table_from_folder(qtbot, use_temp_settings_file):
    # qtbot, window = widget_setup(qtbot, MainWindow)
    assert False


def test_MainWindow_multiple_new_tables(qtbot, use_temp_settings_file):
    # qtbot, window = widget_setup(qtbot, MainWindow)
    assert False


def test_MainWindow_export_pipeline(qtbot, use_temp_settings_file):
    # qtbot, window = widget_setup(qtbot, MainWindow)
    assert False


def test_MainWindow_import_pipeline(qtbot, use_temp_settings_file):
    # qtbot, window = widget_setup(qtbot, MainWindow)
    assert False


def test_MainWindow_import_gene_set(qtbot, use_temp_settings_file):
    # qtbot, window = widget_setup(qtbot, MainWindow)
    assert False


def test_MainWindow_copy_gene_set(qtbot, use_temp_settings_file):
    # qtbot, window = widget_setup(qtbot, MainWindow)
    assert False


def test_MainWindow_append_to_current_console(qtbot, use_temp_settings_file):
    # qtbot, window = widget_setup(qtbot, MainWindow)
    assert False


def test_MainWindow_add_pipeline(qtbot, use_temp_settings_file):
    # qtbot, window = widget_setup(qtbot, MainWindow)
    assert False


def test_MainWindow_export_gene_set(qtbot, use_temp_settings_file):
    # qtbot, window = widget_setup(qtbot, MainWindow)
    assert False


def test_MainWindow_get_available_objects(qtbot, use_temp_settings_file):
    # qtbot, window = widget_setup(qtbot, MainWindow)
    assert False


def test_MainWindow_choose_set_op(qtbot, use_temp_settings_file):
    # qtbot, window = widget_setup(qtbot, MainWindow)
    assert False


def test_MainWindow_visualize_gene_sets(qtbot, use_temp_settings_file):
    # qtbot, window = widget_setup(qtbot, MainWindow)
    assert False


def test_MainWindow_open_enrichment_analysis(qtbot, use_temp_settings_file):
    # qtbot, window = widget_setup(qtbot, MainWindow)
    assert False


def test_MainWindow_save_session(qtbot, use_temp_settings_file):
    # qtbot, window = widget_setup(qtbot, MainWindow)
    assert False


def test_MainWindow_load_session(qtbot, use_temp_settings_file):
    # qtbot, window = widget_setup(qtbot, MainWindow)
    assert False


def test_MainWindow_about(qtbot, use_temp_settings_file):
    # qtbot, window = widget_setup(qtbot, MainWindow)
    assert False


def test_MainWindow_cite(qtbot, use_temp_settings_file):
    # qtbot, window = widget_setup(qtbot, MainWindow)
    assert False
