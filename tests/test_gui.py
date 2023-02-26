import matplotlib
import pytest

matplotlib.use('Agg')
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
def four_available_objects_and_empty(qtbot, red_icon, green_icon, blank_icon):
    qtbot, first = widget_setup(qtbot, SetTabPage, 'first tab',
                                {'WBGene00008447', 'WBGene00044258', 'WBGene00045410', 'WBGene00010100'})

    qtbot, second = widget_setup(qtbot, FilterTabPage, undo_stack=QtWidgets.QUndoStack())
    second.start_from_filter_obj(filtering.DESeqFilter('tests/test_files/test_deseq_set_ops_1.csv'))
    second.rename('second tab')

    qtbot, third = widget_setup(qtbot, FilterTabPage, undo_stack=QtWidgets.QUndoStack())
    third.start_from_filter_obj(filtering.DESeqFilter('tests/test_files/test_deseq_set_ops_2.csv'))
    third.rename('third tab')

    qtbot, fourth = widget_setup(qtbot, FilterTabPage, undo_stack=QtWidgets.QUndoStack())
    fourth.start_from_filter_obj(filtering.CountFilter('tests/test_files/counted.tsv'))
    fourth.rename('fourth tab')

    qtbot, empty = widget_setup(qtbot, FilterTabPage, undo_stack=QtWidgets.QUndoStack())

    return {'first tab': (first, red_icon), 'second tab': (second, red_icon), 'third tab': (third, red_icon),
            'fourth tab': (fourth, green_icon), 'empty tab': (empty, blank_icon)}


@pytest.fixture
def main_window(qtbot, monkeypatch, use_temp_settings_file):
    monkeypatch.setattr(QtWidgets.QMessageBox, 'question', lambda *args, **kwargs: QtWidgets.QMessageBox.Yes)
    monkeypatch.setattr(gui_widgets.ThreadStdOutStreamTextQueueReceiver, 'run', lambda self: None)
    monkeypatch.setattr(gui_quickstart.QuickStartWizard, '__init__', lambda *args, **kwargs: None)
    settings.set_show_tutorial_settings(False)
    qtbot, window = widget_setup(qtbot, MainWindow)
    warnings.showwarning = customwarn
    sys.excepthook = window.excepthook
    builtins.input = window.input
    return window


@pytest.fixture
def main_window_with_tabs(main_window, monkeypatch):
    monkeypatch.setattr(QtWidgets.QFileDialog, 'getOpenFileName',
                        lambda *args, **kwargs: ('tests/test_files/test_session.rnal', '.rnal'))
    main_window.load_session_action.trigger()
    return main_window


@pytest.fixture
def tab_widget(qtbot):
    qtbot, window = widget_setup(qtbot, ReactiveTabWidget)
    return window


@pytest.fixture
def multi_keep_window(qtbot):
    objs = [filtering.DESeqFilter('tests/test_files/test_deseq.csv'),
            filtering.CountFilter('tests/test_files/counted.tsv'),
            filtering.Filter('tests/test_files/test_deseq_biotype.csv')]
    qtbot, window = widget_setup(qtbot, MultiKeepWindow, objs)
    return window


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
def deseq_window(qtbot) -> DESeqWindow:
    qtbot, window = widget_setup(qtbot, DESeqWindow)
    return window


@pytest.fixture
def cutadapt_single_window(qtbot) -> CutAdaptSingleWindow:
    qtbot, window = widget_setup(qtbot, CutAdaptSingleWindow)
    return window


@pytest.fixture
def cutadapt_paired_window(qtbot) -> CutAdaptPairedWindow:
    qtbot, window = widget_setup(qtbot, CutAdaptPairedWindow)
    return window


@pytest.fixture
def kallisto_index_window(qtbot) -> KallistoIndexWindow:
    qtbot, window = widget_setup(qtbot, KallistoIndexWindow)
    return window


@pytest.fixture
def kallisto_single_window(qtbot) -> KallistoSingleWindow:
    qtbot, window = widget_setup(qtbot, KallistoSingleWindow)
    return window


@pytest.fixture
def kallisto_paired_window(qtbot) -> KallistoPairedWindow:
    qtbot, window = widget_setup(qtbot, KallistoPairedWindow)
    return window


def update_gene_sets_widget(widget: gui_widgets.GeneSetComboBox, objs):
    widget.update_gene_sets(objs)


@pytest.fixture
def enrichment_window(qtbot, available_objects):
    qtbot, window = widget_setup(qtbot, EnrichmentWindow)
    window.geneSetsRequested.connect(functools.partial(update_gene_sets_widget, objs=available_objects))
    return window


@pytest.fixture
def set_op_window(qtbot, four_available_objects_and_empty):
    qtbot, window = widget_setup(qtbot, SetOperationWindow, four_available_objects_and_empty)
    return window


@pytest.fixture
def set_vis_window(qtbot, four_available_objects_and_empty):
    qtbot, window = widget_setup(qtbot, SetVisualizationWindow, four_available_objects_and_empty)
    return window


multi_open_window_files = ['tests/counted.csv', 'tests/test_deseq.csv', 'tests/counted.tsv']


@pytest.fixture
def multi_open_window(qtbot):
    qtbot, window = widget_setup(qtbot, MultiOpenWindow, multi_open_window_files)
    return window


@pytest.fixture
def monkeypatch_create_canvas(monkeypatch):
    canvas_created = []

    def mock_create_canvas(self):
        canvas_created.append(True)

    monkeypatch.setattr(SetVisualizationWindow, 'create_canvas', mock_create_canvas)
    return canvas_created


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


def test_KallistoIndexWindow_init(qtbot, kallisto_index_window):
    _ = kallisto_index_window


def test_KallistoSingleWindow_init(qtbot, kallisto_single_window):
    _ = kallisto_single_window


def test_KallistoPairedWindow_init(qtbot, kallisto_paired_window):
    _ = kallisto_paired_window


def test_KallistoIndexWindow_start_analysis(qtbot, kallisto_index_window):
    truth_args = []
    fa_file = 'path/to/fa/file'
    truth_kwargs = dict(transcriptome_fasta=fa_file,
                        kallisto_installation_folder='auto',
                        kmer_length=31,
                        make_unique=False)

    kallisto_index_window.param_widgets['transcriptome_fasta'].setText(fa_file)

    with qtbot.waitSignal(kallisto_index_window.paramsAccepted) as blocker:
        qtbot.mouseClick(kallisto_index_window.start_button, LEFT_CLICK)
    assert blocker.args[0] == truth_args
    assert blocker.args[1] == truth_kwargs


def test_KallistoSingleWindow_start_analysis(qtbot, kallisto_single_window):
    truth_args = []
    fq_folder = 'path/to/fq/folder'
    out_folder = 'path/to/out/dir'
    index_file = 'path/to/index/file.idx'
    gtf_file = 'path/to/gtf.gtf'
    average_fragment_length = 175
    stdev_fragment_length = 25.5
    truth_kwargs = dict(fastq_folder=fq_folder, output_folder=out_folder,
                        gtf_file=gtf_file,
                        index_file=index_file,
                        average_fragment_length=average_fragment_length,
                        stdev_fragment_length=stdev_fragment_length,
                        kallisto_installation_folder='auto',
                        new_sample_names='auto',
                        stranded='no',
                        learn_bias=False,
                        bootstrap_samples=None, seek_fusion_genes=False)

    kallisto_single_window.param_widgets['fastq_folder'].setText(fq_folder)
    kallisto_single_window.param_widgets['output_folder'].setText(out_folder)
    kallisto_single_window.param_widgets['index_file'].setText(index_file)
    kallisto_single_window.param_widgets['gtf_file'].setText(gtf_file)
    kallisto_single_window.param_widgets['average_fragment_length'].setValue(average_fragment_length)
    kallisto_single_window.param_widgets['stdev_fragment_length'].setValue(stdev_fragment_length)

    with qtbot.waitSignal(kallisto_single_window.paramsAccepted) as blocker:
        qtbot.mouseClick(kallisto_single_window.start_button, LEFT_CLICK)
    assert blocker.args[0] == truth_args
    assert blocker.args[1] == truth_kwargs


def test_KallistoPairedWindow_start_analysis(qtbot, kallisto_paired_window):
    r1_files = ['file1.fq', 'path/file2.fq']
    r2_files = ['file3.fq.gz', 'path/to/file4.fastq.gz']
    out_folder = 'path/to/out/dir'
    index_file = 'path/to/index/file.idx'
    gtf_file = 'path/to/gtf.gtf'
    truth_args = []
    truth_kwargs = dict(r1_files=r1_files, r2_files=r2_files,
                        output_folder=out_folder,
                        gtf_file=gtf_file,
                        index_file=index_file,
                        kallisto_installation_folder='auto',
                        new_sample_names='auto',
                        stranded='no',
                        learn_bias=False,
                        bootstrap_samples=None, seek_fusion_genes=False)

    kallisto_paired_window.pairs_widgets['r1_list'].add_items(r1_files)
    kallisto_paired_window.pairs_widgets['r2_list'].add_items(r2_files)
    kallisto_paired_window.param_widgets['output_folder'].setText(out_folder)
    kallisto_paired_window.param_widgets['index_file'].setText(index_file)
    kallisto_paired_window.param_widgets['gtf_file'].setText(gtf_file)

    with qtbot.waitSignal(kallisto_paired_window.paramsAccepted) as blocker:
        qtbot.mouseClick(kallisto_paired_window.start_button, LEFT_CLICK)
    assert blocker.args[0] == truth_args
    assert blocker.args[1] == truth_kwargs


def test_CutAdaptSingleWindow_init(qtbot, cutadapt_single_window):
    _ = cutadapt_single_window


def test_CutAdaptPairedWindow_init(qtbot, cutadapt_paired_window):
    _ = cutadapt_paired_window


def test_CutAdaptSingleWindow_start_analysis(qtbot, cutadapt_single_window):
    truth_args = []
    fq_folder = 'path/to/fq/folder'
    out_folder = 'path/to/out/dir'
    adapter = 'ATGGA'
    truth_kwargs = dict(fastq_folder=fq_folder, output_folder=out_folder,
                        three_prime_adapters=adapter,
                        five_prime_adapters=None,
                        any_position_adapters=None,
                        quality_trimming=20, trim_n=True,
                        minimum_read_length=10, maximum_read_length=None,
                        discard_untrimmed_reads=True, error_tolerance=0.1,
                        minimum_overlap=3, allow_indels=True, parallel=True)

    cutadapt_single_window.param_widgets['fastq_folder'].setText(fq_folder)
    cutadapt_single_window.param_widgets['output_folder'].setText(out_folder)
    cutadapt_single_window.param_widgets['three_prime_adapters'].setValue(adapter)

    with qtbot.waitSignal(cutadapt_single_window.paramsAccepted) as blocker:
        qtbot.mouseClick(cutadapt_single_window.start_button, LEFT_CLICK)
    assert blocker.args[0] == truth_args
    assert blocker.args[1] == truth_kwargs


def test_CutAdaptPairedWindow_start_analysis(qtbot, cutadapt_paired_window):
    r1_files = ['file1.fq', 'path/file2.fq']
    r2_files = ['file3.fq.gz', 'path/to/file4.fastq.gz']
    out_folder = 'path/to/out/dir'
    adapter1 = 'ATGGA'
    adapter2 = 'CATC'
    truth_args = []
    truth_kwargs = dict(r1_files=r1_files, r2_files=r2_files,
                        output_folder=out_folder,
                        three_prime_adapters_r1=adapter1, three_prime_adapters_r2=adapter2,
                        five_prime_adapters_r1=None, five_prime_adapters_r2=None,
                        any_position_adapters_r1=None, any_position_adapters_r2=None,
                        quality_trimming=20, trim_n=True, minimum_read_length=10, maximum_read_length=None,
                        discard_untrimmed_reads=True, pair_filter_if='both',
                        error_tolerance=0.1, minimum_overlap=3, allow_indels=True, parallel=True)
    cutadapt_paired_window.pairs_widgets['r1_list'].add_items(r1_files)
    cutadapt_paired_window.pairs_widgets['r2_list'].add_items(r2_files)
    cutadapt_paired_window.param_widgets['output_folder'].setText(out_folder)
    cutadapt_paired_window.param_widgets['three_prime_adapters_r1'].setValue(adapter1)
    cutadapt_paired_window.param_widgets['three_prime_adapters_r2'].setValue(adapter2)

    with qtbot.waitSignal(cutadapt_paired_window.paramsAccepted) as blocker:
        qtbot.mouseClick(cutadapt_paired_window.start_button, LEFT_CLICK)
    assert blocker.args[0] == truth_args
    assert blocker.args[1] == truth_kwargs


def test_DESeqWindow_init(qtbot, deseq_window):
    _ = deseq_window


def test_DESeqWindow_load_design_mat(qtbot, deseq_window):
    design_mat_path = 'tests/test_files/test_design_matrix.csv'
    design_mat_truth = io.load_csv(design_mat_path, 0)
    deseq_window.param_widgets['design_matrix'].setText(design_mat_path)
    qtbot.mouseClick(deseq_window.param_widgets['load_design'], LEFT_CLICK)
    assert deseq_window.design_mat.equals(design_mat_truth)
    assert deseq_window.comparisons_widgets['picker'].design_mat.equals(design_mat_truth)


def test_DESeqWindow_get_analysis_params(qtbot, deseq_window):
    design_mat_path = 'tests/test_files/test_design_matrix.csv'
    truth = dict(r_installation_folder='auto', design_matrix=design_mat_path, output_folder=None,
                 comparisons=[('replicate', 'rep3', 'rep2'), ('condition', 'cond1', 'cond1')])

    deseq_window.param_widgets['design_matrix'].setText(design_mat_path)
    qtbot.mouseClick(deseq_window.param_widgets['load_design'], LEFT_CLICK)

    deseq_window.comparisons_widgets['picker'].add_comparison_widget()

    deseq_window.comparisons_widgets['picker'].inputs[0].factor.setCurrentText('replicate')
    deseq_window.comparisons_widgets['picker'].inputs[0].numerator.setCurrentText('rep3')
    deseq_window.comparisons_widgets['picker'].inputs[0].denominator.setCurrentText('rep2')

    assert deseq_window.get_analysis_kwargs() == truth


def test_DESeqWindow_start_analysis(qtbot, deseq_window):
    design_mat_path = 'tests/test_files/test_design_matrix.csv'

    truth_args = []
    truth_kwargs = dict(r_installation_folder='auto', design_matrix=design_mat_path, output_folder=None,
                        comparisons=[('replicate', 'rep3', 'rep2'), ('condition', 'cond1', 'cond1')])

    deseq_window.param_widgets['design_matrix'].setText(design_mat_path)
    qtbot.mouseClick(deseq_window.param_widgets['load_design'], LEFT_CLICK)

    deseq_window.comparisons_widgets['picker'].add_comparison_widget()

    deseq_window.comparisons_widgets['picker'].inputs[0].factor.setCurrentText('replicate')
    deseq_window.comparisons_widgets['picker'].inputs[0].numerator.setCurrentText('rep3')
    deseq_window.comparisons_widgets['picker'].inputs[0].denominator.setCurrentText('rep2')

    with qtbot.waitSignal(deseq_window.paramsAccepted) as blocker:
        qtbot.mouseClick(deseq_window.start_button, LEFT_CLICK)
    assert blocker.args[0] == truth_args
    assert blocker.args[1] == truth_kwargs


def test_ClicomWindow_init(qtbot, clicom_window):
    _ = clicom_window


def test_ClicomWindow_add_setup(qtbot, clicom_window):
    truth = dict(method='kmeans', n_clusters=3, n_init=3, max_iter=300, random_seed=None,
                 max_n_clusters_estimate='auto')

    qtbot.keyClicks(clicom_window.stack.func_combo, filtering.CountFilter.split_kmeans.readable_name)
    clicom_window.stack.parameter_widgets['n_clusters'].other.setValue(3)
    qtbot.mouseClick(clicom_window.setups_widgets['add_button'], LEFT_CLICK)
    assert len(clicom_window.parameter_dicts) == 1
    assert clicom_window.parameter_dicts[0] == truth


def test_ClicomWindow_remove_setup(qtbot, monkeypatch, clicom_window):
    monkeypatch.setattr(QtWidgets.QMessageBox, 'question', lambda *args: QtWidgets.QMessageBox.Yes)
    qtbot.keyClicks(clicom_window.stack.func_combo, filtering.CountFilter.split_kmeans.readable_name)
    clicom_window.stack.parameter_widgets['n_clusters'].other.setValue(3)
    qtbot.mouseClick(clicom_window.setups_widgets['add_button'], LEFT_CLICK)
    assert len(clicom_window.parameter_dicts) == 1

    qtbot.mouseClick(clicom_window.setups_widgets['list'].delete_all_button, LEFT_CLICK)

    assert len(clicom_window.parameter_dicts) == 0


def test_ClicomWindow_get_analysis_params(qtbot, clicom_window):
    truth = dict(replicate_grouping='ungrouped', power_transform=[True, False], evidence_threshold=0.35,
                 cluster_unclustered_features=True, parallel_backend='loky',
                 min_cluster_size=15, plot_style='all', split_plots=False)

    qtbot.mouseClick(clicom_window.param_widgets['power_transform'].false_button, LEFT_CLICK)
    qtbot.mouseClick(clicom_window.param_widgets['cluster_unclustered_features'].switch, LEFT_CLICK)
    clicom_window.param_widgets['evidence_threshold'].clear()
    qtbot.keyClicks(clicom_window.param_widgets['evidence_threshold'], '0.35')

    assert clicom_window.get_analysis_kwargs() == truth


def test_ClicomWindow_start_analysis(qtbot, clicom_window):
    truth_setups = [dict(method='kmeans', n_clusters=3, n_init=3, max_iter=300, random_seed=None,
                         max_n_clusters_estimate='auto'),
                    dict(method='hierarchical', n_clusters='silhouette', metric='Euclidean', linkage='Average',
                         distance_threshold=None, max_n_clusters_estimate='auto')]
    truth_params = dict(replicate_grouping='ungrouped', power_transform=[True, False], evidence_threshold=0.35,
                        cluster_unclustered_features=True, min_cluster_size=15, plot_style='all', split_plots=False,
                        parallel_backend='loky')

    clicom_window.stack.func_combo.setCurrentText(filtering.CountFilter.split_kmeans.readable_name)
    clicom_window.stack.parameter_widgets['n_clusters'].other.setValue(3)
    qtbot.mouseClick(clicom_window.setups_widgets['add_button'], LEFT_CLICK)

    clicom_window.stack.func_combo.setCurrentText(filtering.CountFilter.split_hierarchical.readable_name)
    qtbot.keyClicks(clicom_window.stack.parameter_widgets['n_clusters'].combo, 'silhouette')
    qtbot.mouseClick(clicom_window.setups_widgets['add_button'], LEFT_CLICK)

    qtbot.mouseClick(clicom_window.param_widgets['power_transform'].false_button, LEFT_CLICK)
    qtbot.mouseClick(clicom_window.param_widgets['cluster_unclustered_features'].switch, LEFT_CLICK)
    clicom_window.param_widgets['evidence_threshold'].clear()
    qtbot.keyClicks(clicom_window.param_widgets['evidence_threshold'], '0.35')

    with qtbot.waitSignal(clicom_window.paramsAccepted) as blocker:
        clicom_window.start_button.click()
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
     dict(plot_horizontal=True, gene_id_type='auto')),
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

    enrichment_window.widgets['enrichment_list'].showPopup()
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
        enrichment_window.widgets['bg_list'].showPopup()
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

    enrichment_window.widgets['enrichment_list'].showPopup()
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

    enrichment_window.widgets['bg_list'].showPopup()
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
     dict(plot_horizontal=True, gene_id_type='auto')),
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

    enrichment_window.widgets['enrichment_list'].showPopup()
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
        enrichment_window.widgets['bg_list'].showPopup()
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

    enrichment_window.widgets['enrichment_list'].showPopup()
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

    enrichment_window.widgets['bg_list'].showPopup()
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


@pytest.mark.parametrize('op_name,truth', [
    ('Intersection', 'intersection'),
    ('Union', 'union'),
    ('Symmetric Difference', 'symmetric_difference'),
    ('Majority-Vote Intersection', 'majority_vote_intersection'),
    ('Other', 'other')
])
@pytest.mark.parametrize('second_op_name,second_truth', [
    ('Intersection', 'intersection'),
    ('Union', 'union'),
    ('Symmetric Difference', 'symmetric_difference'),
    ('Majority-Vote Intersection', 'majority_vote_intersection'),
    ('Other', 'other')
])
def test_SetOperationWindow_get_current_func_name(qtbot, set_op_window, op_name, truth, second_op_name, second_truth):
    assert set_op_window.get_current_func_name() is None
    set_op_window.widgets['radio_button_box'].radio_buttons[op_name].click()
    assert set_op_window.get_current_func_name() == truth
    set_op_window.widgets['radio_button_box'].radio_buttons[second_op_name].click()
    assert set_op_window.get_current_func_name() == second_truth


def test_SetOperationWindow_canvas_types(qtbot, set_op_window):
    assert isinstance(set_op_window.widgets['canvas'], gui_graphics.EmptyCanvas)

    set_op_window.widgets['set_list'].list_items[0].setSelected(True)
    assert isinstance(set_op_window.widgets['canvas'], gui_graphics.EmptyCanvas)

    set_op_window.widgets['set_list'].list_items[1].setSelected(True)
    assert isinstance(set_op_window.widgets['canvas'], gui_graphics.VennInteractiveCanvas)

    set_op_window.widgets['set_list'].list_items[2].setSelected(True)
    assert isinstance(set_op_window.widgets['canvas'], gui_graphics.VennInteractiveCanvas)

    set_op_window.widgets['set_list'].select_all_button.click()
    assert isinstance(set_op_window.widgets['canvas'], gui_graphics.UpSetInteractiveCanvas)

    set_op_window.widgets['set_list'].list_items[0].setSelected(False)
    assert isinstance(set_op_window.widgets['canvas'], gui_graphics.VennInteractiveCanvas)

    set_op_window.widgets['set_list'].clear_all_button.click()
    assert isinstance(set_op_window.widgets['canvas'], gui_graphics.EmptyCanvas)


@pytest.mark.parametrize('n_selected', [3, 4])
def test_SetOperationWindow_primary_set_change(qtbot, set_op_window, n_selected):
    for i in range(n_selected):
        set_op_window.widgets['set_list'].list_items[i].setSelected(True)

    with qtbot.waitSignal(set_op_window.primarySetChangedDifference):
        set_op_window.widgets['radio_button_box'].radio_buttons['Difference'].click()

    for tab in ['first tab', 'second tab', 'third tab']:
        with qtbot.waitSignal(set_op_window.primarySetChangedDifference) as blocker:
            set_op_window.widgets['choose_primary_set'].setCurrentText(tab)
            print(qtbot.screenshot(set_op_window))
        assert blocker.args[0] == tab

    with qtbot.waitSignal(set_op_window.primarySetChangedIntersection):
        set_op_window.widgets['radio_button_box'].radio_buttons['Intersection'].click()

    for tab in ['first tab', 'second tab', 'third tab']:
        with qtbot.waitSignal(set_op_window.primarySetChangedIntersection):
            set_op_window.widgets['choose_primary_set'].setCurrentText(tab)


apply_set_ops_parametrize = [
    ('Union', [0, 2],
     {'WBGene00018199', 'WBGene00020407', 'WBGene00045366', 'WBGene00044258', 'WBGene00010100', 'WBGene00018193',
      'WBGene00219307', 'WBGene00021019', 'WBGene00045410', 'WBGene00194708', 'WBGene00021589', 'WBGene00219304',
      'WBGene00023036', 'WBGene00021375', 'WBGene00008447', 'WBGene00044799', 'WBGene00001118', 'WBGene00077437',
      'WBGene00010755', 'WBGene00012919', 'WBGene00021654', 'WBGene00013816', 'WBGene00022486', 'WBGene00019174',
      'WBGene00007674', 'WBGene00012648', 'WBGene00021605'}
     ),
    ('Union', [1, 2, 3],
     {'WBGene00018199', 'WBGene00007064', 'WBGene00020407', 'WBGene00007079', 'WBGene00044478', 'WBGene00045366',
      'WBGene00043989', 'WBGene00007075', 'WBGene00044258', 'WBGene00010100', 'WBGene00043987', 'WBGene00007066',
      'WBGene00018193', 'WBGene00022730', 'WBGene00044022', 'WBGene00077504', 'WBGene00219307', 'WBGene00014997',
      'WBGene00021019', 'WBGene00043990', 'WBGene00045410', 'WBGene00021018', 'WBGene00194708', 'WBGene00007078',
      'WBGene00021589', 'WBGene00219304', 'WBGene00023036', 'WBGene00007069', 'WBGene00021375', 'WBGene00007076',
      'WBGene00008447', 'WBGene00044799', 'WBGene00001118', 'WBGene00077502', 'WBGene00007067', 'WBGene00077503',
      'WBGene00007071', 'WBGene00012961', 'WBGene00077437', 'WBGene00022438', 'WBGene00010755', 'WBGene00007063',
      'WBGene00012919', 'WBGene00021654', 'WBGene00013816', 'WBGene00007074', 'WBGene00010507', 'WBGene00016635',
      'WBGene00022486', 'WBGene00043988', 'WBGene00007077', 'WBGene00019174', 'WBGene00012452', 'WBGene00007674',
      'WBGene00012648', 'WBGene00044951', 'WBGene00021605'}
     ),
    ('Intersection', [0, 2], {'WBGene00044258', 'WBGene00045410', 'WBGene00010100'}),
    ('Intersection', [1, 2, 3], set()),
    ('Difference', [0, 2], {'WBGene00008447'}),
    ('Difference', [1, 2, 3],
     {'WBGene00044478', 'WBGene00008447', 'WBGene00021018', 'WBGene00010507', 'WBGene00016635', 'WBGene00012452',
      'WBGene00022730', 'WBGene00012961', 'WBGene00022438'}
     ),
    ('Symmetric Difference', [0, 1],
     {'WBGene00018199', 'WBGene00044478', 'WBGene00045366', 'WBGene00022730', 'WBGene00219307', 'WBGene00021019',
      'WBGene00021018', 'WBGene00194708', 'WBGene00219304', 'WBGene00023036', 'WBGene00021375', 'WBGene00012961',
      'WBGene00077437', 'WBGene00022438', 'WBGene00013816', 'WBGene00010507', 'WBGene00016635', 'WBGene00022486',
      'WBGene00019174', 'WBGene00012452', 'WBGene00007674', 'WBGene00012648'}
     ),
    ('Symmetric Difference', [2, 3],
     {'WBGene00018199', 'WBGene00045366', 'WBGene00043989', 'WBGene00043987', 'WBGene00007066', 'WBGene00219307',
      'WBGene00021019', 'WBGene00043990', 'WBGene00007078', 'WBGene00219304', 'WBGene00023036', 'WBGene00044799',
      'WBGene00077502', 'WBGene00001118', 'WBGene00007067', 'WBGene00077437', 'WBGene00010755', 'WBGene00007063',
      'WBGene00021654', 'WBGene00013816', 'WBGene00007674', 'WBGene00012648', 'WBGene00007064', 'WBGene00020407',
      'WBGene00007079', 'WBGene00007075', 'WBGene00044258', 'WBGene00010100', 'WBGene00021605', 'WBGene00018193',
      'WBGene00044022', 'WBGene00077504', 'WBGene00045410', 'WBGene00194708', 'WBGene00021589', 'WBGene00007069',
      'WBGene00021375', 'WBGene00007076', 'WBGene00077503', 'WBGene00007071', 'WBGene00012919', 'WBGene00007074',
      'WBGene00043988', 'WBGene00007077', 'WBGene00022486', 'WBGene00019174', 'WBGene00044951', 'WBGene00014997'}
     )
]


@pytest.mark.parametrize('operation,set_indices,truth', apply_set_ops_parametrize)
def test_SetOperationWindow_apply_set_op(qtbot, set_op_window, operation, set_indices, truth):
    for ind in set_indices:
        set_op_window.widgets['set_list'].list_items[ind].setSelected(True)
    set_op_window.widgets['radio_button_box'].radio_buttons[operation].click()
    if operation in ['Difference', 'Intersection']:
        set_op_window.widgets['choose_primary_set'].setCurrentText(
            set_op_window.widgets['set_list'].items[set_indices[0]])
    with qtbot.waitSignal(set_op_window.geneSetReturned) as blocker:
        set_op_window.widgets['apply_button'].click()
    assert blocker.args[0] == truth


@pytest.mark.parametrize('operation,set_indices,truth', apply_set_ops_parametrize)
def test_SetOperationWindow_apply_set_op_other(qtbot, set_op_window, operation, set_indices, truth):
    for ind in set_indices:
        set_op_window.widgets['set_list'].list_items[ind].setSelected(True)
    set_op_window.widgets['radio_button_box'].radio_buttons[operation].click()
    if operation in ['Difference', 'Intersection']:
        set_op_window.widgets['choose_primary_set'].setCurrentText(
            set_op_window.widgets['set_list'].items[set_indices[0]])

    set_op_window.widgets['radio_button_box'].radio_buttons['Other'].click()
    with qtbot.waitSignal(set_op_window.geneSetReturned) as blocker:
        set_op_window.widgets['apply_button'].click()
    assert blocker.args[0] == truth


@pytest.mark.parametrize('operation,set_indices,primary_set,truth',
                         [('Intersection', [0, 2], 2, {'WBGene00044258', 'WBGene00045410', 'WBGene00010100'}),
                          ('Intersection', [1, 2, 3], 1, set()),
                          ('Difference', [0, 2], 2,
                           {'WBGene00001118', 'WBGene00007674', 'WBGene00010755', 'WBGene00012648', 'WBGene00012919',
                            'WBGene00013816', 'WBGene00018193', 'WBGene00018199', 'WBGene00019174', 'WBGene00020407',
                            'WBGene00021019', 'WBGene00021375', 'WBGene00021589', 'WBGene00021605', 'WBGene00021654',
                            'WBGene00022486', 'WBGene00023036', 'WBGene00044799', 'WBGene00045366', 'WBGene00077437',
                            'WBGene00194708', 'WBGene00219304', 'WBGene00219307'}),
                          ('Difference', [1, 2, 3], 1,
                           {'WBGene00044478', 'WBGene00008447', 'WBGene00021018', 'WBGene00010507', 'WBGene00016635',
                            'WBGene00012452',
                            'WBGene00022730', 'WBGene00012961', 'WBGene00022438'})])
def test_SetOperationWindow_apply_set_op_inplace(qtbot, four_available_objects_and_empty, set_op_window, operation,
                                                 set_indices, primary_set, truth):
    primary_set_name = set_op_window.widgets['set_list'].items[primary_set]

    for ind in set_indices:
        set_op_window.widgets['set_list'].list_items[ind].setSelected(True)
    set_op_window.widgets['radio_button_box'].radio_buttons[operation].click()
    set_op_window.widgets['choose_primary_set'].setCurrentText(
        primary_set_name)

    with qtbot.waitSignal(set_op_window.geneSetReturned) as blocker:
        set_op_window.widgets['apply_button'].click()
    assert blocker.args[0] == truth

    inplace_truth = four_available_objects_and_empty[primary_set_name][0].obj().__copy__()
    obj_names = [set_op_window.widgets['set_list'].items[ind] for ind in set_indices if ind != primary_set]
    objs_for_operation = [four_available_objects_and_empty[name][0].obj() for name in obj_names]
    if operation == 'Difference':
        inplace_truth.difference(
            *objs_for_operation, inplace=True)
    else:
        inplace_truth.intersection(*objs_for_operation, inplace=True)
    set_op_window.parameter_widgets['inplace'].switch.click()
    with qtbot.assertNotEmitted(set_op_window.geneSetReturned) as blocker:
        set_op_window.widgets['apply_button'].click()
    assert four_available_objects_and_empty[primary_set_name][0].obj() == inplace_truth


@pytest.mark.parametrize('threshold,truth', [
    (0, {'WBGene00194708', 'WBGene00044951', 'WBGene00018193', 'WBGene00022730', 'WBGene00012919', 'WBGene00044022',
         'WBGene00044799', 'WBGene00001118', 'WBGene00007069', 'WBGene00021375', 'WBGene00021654', 'WBGene00077437',
         'WBGene00010507', 'WBGene00043987', 'WBGene00010755', 'WBGene00012648', 'WBGene00077503', 'WBGene00007079',
         'WBGene00010100', 'WBGene00012452', 'WBGene00013816', 'WBGene00022438', 'WBGene00012961', 'WBGene00016635',
         'WBGene00007064', 'WBGene00219307', 'WBGene00043989', 'WBGene00007063', 'WBGene00023036', 'WBGene00007078',
         'WBGene00043988', 'WBGene00077504', 'WBGene00007066', 'WBGene00007674', 'WBGene00044258', 'WBGene00021589',
         'WBGene00021605', 'WBGene00021019', 'WBGene00007071', 'WBGene00219304', 'WBGene00043990', 'WBGene00014997',
         'WBGene00045410', 'WBGene00077502', 'WBGene00020407', 'WBGene00007075', 'WBGene00018199', 'WBGene00045366',
         'WBGene00007067', 'WBGene00044478', 'WBGene00022486', 'WBGene00007074', 'WBGene00007076', 'WBGene00007077',
         'WBGene00008447', 'WBGene00019174', 'WBGene00021018'}),
    (0.25, {'WBGene00194708', 'WBGene00044951', 'WBGene00018193', 'WBGene00022730', 'WBGene00012919', 'WBGene00044022',
            'WBGene00044799', 'WBGene00001118', 'WBGene00007069', 'WBGene00021375', 'WBGene00021654', 'WBGene00077437',
            'WBGene00010507', 'WBGene00043987', 'WBGene00010755', 'WBGene00012648', 'WBGene00077503', 'WBGene00007079',
            'WBGene00010100', 'WBGene00012452', 'WBGene00013816', 'WBGene00022438', 'WBGene00012961', 'WBGene00016635',
            'WBGene00007064', 'WBGene00219307', 'WBGene00043989', 'WBGene00007063', 'WBGene00023036', 'WBGene00007078',
            'WBGene00043988', 'WBGene00077504', 'WBGene00007066', 'WBGene00007674', 'WBGene00044258', 'WBGene00021589',
            'WBGene00021605', 'WBGene00021019', 'WBGene00007071', 'WBGene00219304', 'WBGene00043990', 'WBGene00014997',
            'WBGene00045410', 'WBGene00077502', 'WBGene00020407', 'WBGene00007075', 'WBGene00018199', 'WBGene00045366',
            'WBGene00007067', 'WBGene00044478', 'WBGene00022486', 'WBGene00007074', 'WBGene00007076', 'WBGene00007077',
            'WBGene00008447', 'WBGene00019174', 'WBGene00021018'}),
    (0.57, {'WBGene00044258', 'WBGene00010100', 'WBGene00045410'}),
    (0.99, set()),
    (1, set())
])
def test_SetOperationWindow_apply_set_op_majority_vote(qtbot, set_op_window, threshold, truth):
    for ind in range(4):
        set_op_window.widgets['set_list'].list_items[ind].setSelected(True)
    set_op_window.widgets['radio_button_box'].radio_buttons['Majority-Vote Intersection'].click()
    set_op_window.parameter_widgets['majority_threshold'].setValue(threshold)
    with qtbot.waitSignal(set_op_window.geneSetReturned) as blocker:
        set_op_window.widgets['apply_button'].click()
    assert blocker.args[0] == truth


def test_SetVisualizationWindow_init(qtbot, set_vis_window):
    _ = set_vis_window


@pytest.mark.parametrize('op_name,truth', [
    ('Venn Diagram', 'venn_diagram'),
    ('UpSet Plot', 'upset_plot')
])
@pytest.mark.parametrize('second_op_name,second_truth', [
    ('Venn Diagram', 'venn_diagram'),
    ('UpSet Plot', 'upset_plot')
])
def test_SetVisualizationWindow_get_current_func_name(qtbot, set_vis_window, op_name, truth, second_op_name,
                                                      second_truth):
    assert set_vis_window.get_current_func_name() is None
    set_vis_window.widgets['radio_button_box'].radio_buttons[op_name].click()
    assert set_vis_window.get_current_func_name() == truth
    set_vis_window.widgets['radio_button_box'].radio_buttons[second_op_name].click()
    assert set_vis_window.get_current_func_name() == second_truth


@pytest.mark.parametrize('is_func_selected', ['Venn Diagram', 'UpSet Plot', False])
def test_SetVisualizationWindow_canvas_types(qtbot, set_vis_window, is_func_selected):
    expected_canvas = gui_graphics.BasePreviewCanvas if is_func_selected else gui_graphics.EmptyCanvas
    if is_func_selected:
        qtbot.mouseClick(set_vis_window.widgets['radio_button_box'].radio_buttons[is_func_selected], LEFT_CLICK)
    assert isinstance(set_vis_window.widgets['canvas'], gui_graphics.EmptyCanvas)

    set_vis_window.widgets['set_list'].list_items[0].setSelected(True)
    assert isinstance(set_vis_window.widgets['canvas'], gui_graphics.EmptyCanvas)

    set_vis_window.widgets['set_list'].list_items[1].setSelected(True)
    assert isinstance(set_vis_window.widgets['canvas'], expected_canvas)

    set_vis_window.widgets['set_list'].list_items[2].setSelected(True)
    assert isinstance(set_vis_window.widgets['canvas'], expected_canvas)

    set_vis_window.widgets['set_list'].select_all_button.click()
    assert isinstance(set_vis_window.widgets['canvas'], expected_canvas)

    set_vis_window.widgets['set_list'].list_items[0].setSelected(False)
    assert isinstance(set_vis_window.widgets['canvas'], expected_canvas)

    set_vis_window.widgets['set_list'].clear_all_button.click()
    assert isinstance(set_vis_window.widgets['canvas'], gui_graphics.EmptyCanvas)


@pytest.mark.parametrize('op_name', [
    'Venn Diagram',
    'UpSet Plot'
])
@pytest.mark.parametrize('second_op_name', [
    'Venn Diagram',
    'UpSet Plot'
])
def test_SetVisualizationWindow_function_change_canvas(monkeypatch_create_canvas, set_vis_window, qtbot, op_name,
                                                       second_op_name):
    n_sets = 3
    for i in range(n_sets):
        set_vis_window.widgets['set_list'].list_items[i].setSelected(True)

    while len(monkeypatch_create_canvas) > 0:
        monkeypatch_create_canvas.pop(-1)

    set_vis_window.widgets['radio_button_box'].radio_buttons[op_name].click()
    assert monkeypatch_create_canvas == [True]
    set_vis_window.widgets['radio_button_box'].radio_buttons[second_op_name].click()
    assert monkeypatch_create_canvas == [True, True]


@pytest.mark.parametrize('op_name,n_sets,sample_bool_param', [
    ('Venn Diagram', 2, 'weighted'),
    ('UpSet Plot', 4, 'show_percentages')
])
def test_SetVisualizationWindow_parameter_change_canvas(monkeypatch, qtbot, set_vis_window, op_name, n_sets,
                                                        sample_bool_param):
    canvas_created = []

    def mock_create_canvas(self):
        canvas_created.append(True)

    monkeypatch.setattr(SetVisualizationWindow, 'create_canvas', mock_create_canvas)

    set_vis_window.widgets['radio_button_box'].radio_buttons[op_name].click()
    for i in range(n_sets):
        set_vis_window.widgets['set_list'].list_items[i].setSelected(True)

    set_vis_window.parameter_widgets['title_fontsize'].setValue(27)

    assert canvas_created == [True]

    qtbot.mouseClick(set_vis_window.parameter_widgets[sample_bool_param].switch, LEFT_CLICK)

    assert canvas_created == [True, True]


@pytest.mark.parametrize('func_name,op_name,n_sets,kwargs_truth', [
    ('venn_diagram', 'Venn Diagram', 2, {'title': 'default', 'weighted': True, 'transparency': 0.4}),
    ('venn_diagram', 'Venn Diagram', 3, {'title': 'default', 'weighted': True, 'linestyle': 'solid'}),
    ('upset_plot', 'UpSet Plot', 2, {'title': 'UpSet Plot', 'title_fontsize': 20}),
    ('upset_plot', 'UpSet Plot', 4, {'title': 'UpSet Plot', 'show_percentages': True}),

])
def test_SetVisualizationWindow_generate_graph(qtbot, set_vis_window, monkeypatch, func_name, op_name, n_sets,
                                               kwargs_truth, four_available_objects_and_empty):
    called = []

    def mock_func(*args, **kwargs):
        for key in kwargs_truth:
            assert kwargs[key] == kwargs_truth[key]
        assert 'fig' not in kwargs
        assert len(args) == 1
        objs_truth = {
            list(four_available_objects_and_empty.keys())[i]:
                four_available_objects_and_empty[list(four_available_objects_and_empty.keys())[i]][
                    0].obj() for i in range(n_sets)}
        assert args[0] == objs_truth

        called.append(True)

    set_vis_window.widgets['radio_button_box'].radio_buttons[op_name].click()
    for i in range(n_sets):
        set_vis_window.widgets['set_list'].list_items[i].setSelected(True)

    monkeypatch.setattr(enrichment, func_name, mock_func)

    qtbot.mouseClick(set_vis_window.widgets['generate_button'], LEFT_CLICK)
    assert called == [True]


def test_FilterTabPage_init(qtbot):
    _, _ = widget_setup(qtbot, FilterTabPage)


def test_FilterTabPage_load_file(qtbot):
    obj_truth = filtering.CountFilter('tests/test_files/counted.csv')
    qtbot, window = widget_setup(qtbot, FilterTabPage)
    assert window.is_empty()
    assert not window.basic_widgets['start_button'].isEnabled()

    window.basic_widgets['file_path'].clear()
    qtbot.keyClicks(window.basic_widgets['file_path'].file_path, str(Path('tests/test_files/counted.csv').absolute()))
    window.basic_widgets['table_type_combo'].setCurrentText('Count matrix')
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
    window.overview_widgets['table_name'].setText(new_name)
    with qtbot.waitSignal(window.tabNameChange) as blocker:
        window.overview_widgets['rename_button'].click()
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
        window.overview_widgets['rename_button'].click()
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


@pytest.mark.parametrize('tab_type,args', [(FilterTabPage, [])])
def test_save_table_empty_tab(qtbot, tab_type, args):
    qtbot, window = widget_setup(qtbot, tab_type, *args)
    with qtbot.assertNotEmitted(window.tabSaved):
        window.save_file()


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
    filtertabpage.overview_widgets['view_button'].click()
    assert isinstance(filtertabpage.overview_widgets['full_table_view'], gui_windows.DataFrameView)
    assert filtertabpage.overview_widgets['full_table_view'].data_view.model()._dataframe.equals(filtertabpage.obj().df)


def test_FilterTabPage_apply_function(qtbot, filtertabpage_with_undo_stack):
    window, stack = filtertabpage_with_undo_stack
    orig = window.obj().__copy__()
    truth = window.obj().filter_significant(0.01, opposite=True, inplace=False)
    window.stack_buttons[0].click()
    window.stack.currentWidget().func_combo.setCurrentText(filtering.DESeqFilter.filter_significant.readable_name)
    window.stack.currentWidget().parameter_widgets['alpha'].clear()
    qtbot.keyClicks(window.stack.currentWidget().parameter_widgets['alpha'], '0.01')
    qtbot.mouseClick(window.stack.currentWidget().parameter_widgets['opposite'].switch, LEFT_CLICK)
    qtbot.mouseClick(window.stack.currentWidget().parameter_widgets['inplace'].switch, LEFT_CLICK)
    with qtbot.waitSignal(window.filterObjectCreated, timeout=10000) as blocker:
        qtbot.mouseClick(window.apply_button, LEFT_CLICK)
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
        window.process_outputs(partial()[0], func_name)

    window.startedClustering.connect(my_slot)

    orig = window.obj().__copy__()

    truth = orig.split_kmeans(n_clusters=3, random_seed=0)
    truth = sorted(truth, key=lambda obj: sorted(obj.df.index)[0])

    window.stack_buttons[4].click()
    qtbot.keyClicks(window.stack.currentWidget().func_combo, filtering.CountFilter.split_kmeans.readable_name)
    window.stack.currentWidget().parameter_widgets['n_clusters'].other.setValue(3)
    qtbot.mouseClick(window.stack.currentWidget().parameter_widgets['random_seed'].checkbox, LEFT_CLICK)
    with qtbot.waitSignals([window.filterObjectCreated, window.filterObjectCreated, window.filterObjectCreated],
                           timeout=15000) as blocker:
        qtbot.mouseClick(window.apply_button, LEFT_CLICK)

    res = sorted([sig.args[0] for sig in blocker.all_signals_and_args], key=lambda obj: sorted(obj.df.index)[0])
    for i in range(3):
        res[i].df.sort_index(inplace=True)
        truth[i].df.sort_index(inplace=True)
        assert np.all(np.isclose(res[i].df, truth[i].df, equal_nan=True))

    assert window.obj() == orig


def test_FilterTabPage_apply_function_inplace(qtbot, filtertabpage_with_undo_stack):
    window, stack = filtertabpage_with_undo_stack
    truth = window.obj().filter_significant(0.01, opposite=True, inplace=False)
    window.stack_buttons[0].click()
    window.stack.currentWidget().func_combo.setCurrentText(filtering.DESeqFilter.filter_significant.readable_name)
    window.stack.currentWidget().parameter_widgets['alpha'].clear()
    qtbot.keyClicks(window.stack.currentWidget().parameter_widgets['alpha'], '0.01')
    qtbot.mouseClick(window.stack.currentWidget().parameter_widgets['opposite'].switch, LEFT_CLICK)
    qtbot.mouseClick(window.apply_button, LEFT_CLICK)
    assert window.obj() == truth


def test_FilterTabPage_undo_function(qtbot, filtertabpage_with_undo_stack):
    window, stack = filtertabpage_with_undo_stack
    truth = window.obj().filter_significant(0.01, opposite=True, inplace=False)
    orig = window.obj().__copy__()
    window.stack_buttons[0].click()
    window.stack.currentWidget().func_combo.setCurrentText(filtering.DESeqFilter.filter_significant.readable_name)
    window.stack.currentWidget().parameter_widgets['alpha'].clear()
    qtbot.keyClicks(window.stack.currentWidget().parameter_widgets['alpha'], '0.01')
    qtbot.mouseClick(window.stack.currentWidget().parameter_widgets['opposite'].switch, LEFT_CLICK)
    qtbot.mouseClick(window.apply_button, LEFT_CLICK)
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


def test_FilterTabPage_open_clicom(countfiltertabpage_with_undo_stack, qtbot, monkeypatch):
    tabpage = countfiltertabpage_with_undo_stack[0]
    opened = []

    def mock_show(*args, **kwargs):
        opened.append(True)

    monkeypatch.setattr(ClicomWindow, 'show', mock_show)

    tabpage.stack_buttons[4].click()
    tabpage.stack_widgets['Cluster'].func_combo.setCurrentText(filtering.CountFilter.split_clicom.readable_name)

    assert opened == [True]


def test_FilterTabPage_open_deseq(countfiltertabpage_with_undo_stack, qtbot, monkeypatch):
    tabpage = countfiltertabpage_with_undo_stack[0]
    opened = []

    def mock_show(*args, **kwargs):
        opened.append(True)

    monkeypatch.setattr(DESeqWindow, 'show', mock_show)

    tabpage.stack_buttons[5].click()
    tabpage.stack_widgets['General'].func_combo.setCurrentText(
        filtering.CountFilter.differential_expression_deseq2.readable_name)

    assert opened == [True]


def test_FilterTabPage_get_all_actions(qtbot, countfiltertabpage_with_undo_stack, filtertabpage_with_undo_stack):
    countfilter = countfiltertabpage_with_undo_stack[0]
    deseqfilter = filtertabpage_with_undo_stack[0]
    truth_counts = {'Filter': [], 'Summarize': [], 'Visualize': [], 'Normalize': [], 'Cluster': [], 'General': []}
    truth_deseq = {'Filter': [], 'Summarize': [], 'Visualize': [], 'Normalize': [], 'Cluster': [], 'General': []}

    counts_res = countfilter.get_all_actions()
    deseq_res = deseqfilter.get_all_actions()

    assert sorted(counts_res.keys()) == sorted(truth_counts.keys())
    assert sorted(deseq_res.keys()) == sorted(truth_deseq.keys())

    assert len(counts_res['Cluster']) >= 5
    assert len(deseq_res['Cluster']) == 0

    for res in (counts_res, deseq_res):
        for action in res['Filter']:
            assert 'filter' in action or action.startswith('split')
        for action in res['Normalize']:
            assert action.startswith('normalize')
        for action in ['sort', 'transform']:
            assert action in res['General']
        for action in itertools.chain(res['General'], res['Visualize'], res['Summarize']):
            for keyword in ['split', 'filter', 'normalize']:
                assert keyword not in action


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
    assert window.obj().gene_set == obj
    assert window.obj_type() == enrichment.FeatureSet
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
    assert genes_in_view == window.obj().gene_set


@pytest.mark.parametrize('exc_params', [None, ['self', 'other']])
@pytest.mark.parametrize('pipeline_mode', [True, False])
def test_FuncTypeStack_init(qtbot, pipeline_mode, exc_params):
    qtbot, stack = widget_setup(qtbot, FuncTypeStack, ['filter_biotype_from_ref_table', 'number_filters', 'describe'],
                                filtering.Filter('tests/test_files/test_deseq.csv'),
                                additional_excluded_params=exc_params, pipeline_mode=pipeline_mode)


def test_CreatePipelineWindow_init(qtbot):
    _, _ = widget_setup(qtbot, CreatePipelineWindow)


def test_CreatePipelineWindow_from_pipeline(qtbot):
    name = 'my pipeline name'
    pipeline = filtering.Pipeline.import_pipeline('tests/test_files/test_pipeline.yaml')
    qtbot, window = widget_setup(qtbot, CreatePipelineWindow.start_from_pipeline, pipeline, name)
    assert window.pipeline == pipeline
    assert window._get_pipeline_name() == name


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
    qtbot.keyClicks(window.stack.currentWidget().func_combo,
                    filtering.DESeqFilter.split_fold_change_direction.readable_name)
    qtbot.mouseClick(window.apply_button, LEFT_CLICK)

    assert window.pipeline == pipeline_truth


def test_CreatePipelineWindow_remove_function(qtbot, monkeypatch):
    warned = []
    monkeypatch.setattr(QtWidgets.QMessageBox, "question", lambda *args: QtWidgets.QMessageBox.Yes)
    monkeypatch.setattr(QtWidgets.QMessageBox, 'exec', lambda *args, **kwargs: warned.append(True))

    qtbot, window = widget_setup(qtbot, CreatePipelineWindow)
    window.basic_widgets['pipeline_name'].clear()
    qtbot.keyClicks(window.basic_widgets['pipeline_name'], 'pipeline_name')
    qtbot.keyClicks(window.basic_widgets['table_type_combo'], 'Differential expression')
    qtbot.mouseClick(window.basic_widgets['start_button'], LEFT_CLICK)

    qtbot.mouseClick(window.stack_buttons[0], LEFT_CLICK)
    qtbot.keyClicks(window.stack.currentWidget().func_combo,
                    filtering.DESeqFilter.split_fold_change_direction.readable_name)
    qtbot.mouseClick(window.apply_button, LEFT_CLICK)

    assert len(window.pipeline) == 1

    qtbot.mouseClick(window.overview_widgets['remove_button'], LEFT_CLICK)
    assert len(window.pipeline) == 0
    assert warned == []

    qtbot.mouseClick(window.overview_widgets['remove_button'], LEFT_CLICK)
    assert len(window.pipeline) == 0
    assert warned == [True]


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
    window.stack.currentWidget().func_combo.setCurrentText(filtering.DESeqFilter.filter_significant.readable_name)
    window.stack.currentWidget().parameter_widgets['alpha'].clear()
    qtbot.keyClicks(window.stack.currentWidget().parameter_widgets['alpha'], '0.01')
    qtbot.mouseClick(window.stack.currentWidget().parameter_widgets['opposite'].switch, LEFT_CLICK)
    qtbot.mouseClick(window.apply_button, LEFT_CLICK)
    assert window.pipeline == pipeline_truth

    # add a second function
    pipeline_truth.add_function('split_fold_change_direction')

    window.stack.currentWidget().func_combo.setCurrentText(
        filtering.DESeqFilter.split_fold_change_direction.readable_name)
    qtbot.mouseClick(window.apply_button, LEFT_CLICK)
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
    qtbot.keyClicks(window.stack.currentWidget().func_combo, filtering.Filter.describe.readable_name)
    qtbot.mouseClick(window.apply_button, LEFT_CLICK)

    with qtbot.waitSignal(window.pipelineSaved) as blocker:
        qtbot.mouseClick(window.overview_widgets['save_button'], LEFT_CLICK)
    assert blocker.args[0] == pipeline_name
    assert blocker.args[1] == pipeline_truth


def test_CreatePipelineWindow_export_pipeline(qtbot, monkeypatch):
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
    qtbot.keyClicks(window.stack.currentWidget().func_combo, filtering.Filter.describe.readable_name)
    qtbot.mouseClick(window.apply_button, LEFT_CLICK)

    with qtbot.waitSignal(window.pipelineExported) as blocker:
        qtbot.mouseClick(window.overview_widgets['export_button'], LEFT_CLICK)
    assert blocker.args[0] == pipeline_name
    assert blocker.args[1] == pipeline_truth


def test_MultiKeepWindow_init(qtbot, multi_keep_window):
    _ = multi_keep_window


@pytest.mark.parametrize('keep_ops,name_ops,truth', [
    ({}, {}, []),
    ({'all': True}, {}, [filtering.DESeqFilter('tests/test_files/test_deseq.csv'),
                         filtering.CountFilter('tests/test_files/counted.tsv'),
                         filtering.Filter('tests/test_files/test_deseq_biotype.csv')]),
    ({'test_deseq': True, 'test_deseq_biotype': True}, {},
     [filtering.DESeqFilter('tests/test_files/test_deseq.csv'),
      filtering.Filter('tests/test_files/test_deseq_biotype.csv')]),
    ({'test_deseq': True, 'test_deseq_biotype': True},
     {'test_deseq': 'new name1', 'test_deseq_biotype': 'new name 2', 'counted': 'new name 3'},
     [filtering.DESeqFilter('tests/test_files/test_deseq.csv'),
      filtering.Filter('tests/test_files/test_deseq_biotype.csv')])
])
def test_MultiKeepWindow_result(qtbot, multi_keep_window, keep_ops, name_ops, truth):
    for ind, op in keep_ops.items():
        if ind == 'all':
            qtbot.mouseClick(multi_keep_window.select_all, LEFT_CLICK)
        else:
            qtbot.mouseClick(multi_keep_window.keep_marks[ind], LEFT_CLICK)
    for ind, op in name_ops.items():
        qtbot.keyClicks(multi_keep_window.names[ind], op)

        for item in truth:
            if item.fname.stem == ind:
                item.fname = Path(op)

    assert multi_keep_window.result() == truth


def test_MultiOpenWindow_init(qtbot, multi_open_window):
    _ = multi_open_window


@pytest.mark.parametrize('path_ops,type_ops,name_ops,truth', [
    ({}, {}, {}, ({'tests/counted.csv': 'tests/counted.csv', 'tests/test_deseq.csv': 'tests/test_deseq.csv',
                   'tests/counted.tsv': 'tests/counted.tsv'},
                  {'tests/counted.csv': 'Other', 'tests/test_deseq.csv': 'Other', 'tests/counted.tsv': 'Other'},
                  {'tests/counted.csv': '', 'tests/test_deseq.csv': '', 'tests/counted.tsv': ''},
                  {'tests/counted.csv': {'drop_columns': []},
                   'tests/counted.tsv': {'drop_columns': []},
                   'tests/test_deseq.csv': {'drop_columns': []}})),
    (
        {}, {}, {1: 'new name', 2: 'second new name'},
        ({'tests/counted.csv': 'tests/counted.csv', 'tests/test_deseq.csv': 'tests/test_deseq.csv',
          'tests/counted.tsv': 'tests/counted.tsv'},
         {'tests/counted.csv': 'Other', 'tests/test_deseq.csv': 'Other',
          'tests/counted.tsv': 'Other'},
         {'tests/counted.csv': '', 'tests/test_deseq.csv': 'new name', 'tests/counted.tsv': 'second new name'},
         {'tests/counted.csv': {'drop_columns': []},
          'tests/counted.tsv': {'drop_columns': []},
          'tests/test_deseq.csv': {'drop_columns': []}})),

    ({}, {0: 'Count matrix', 1: 'Differential expression'}, {},
     ({'tests/counted.csv': 'tests/counted.csv', 'tests/test_deseq.csv': 'tests/test_deseq.csv',
       'tests/counted.tsv': 'tests/counted.tsv'},
      {'tests/counted.csv': 'Count matrix', 'tests/test_deseq.csv': 'Differential expression',
       'tests/counted.tsv': 'Other'},
      {'tests/counted.csv': '', 'tests/test_deseq.csv': '', 'tests/counted.tsv': ''},
      {'tests/counted.csv': {'drop_columns': [], 'is_normalized': False},
       'tests/counted.tsv': {'drop_columns': []},
       'tests/test_deseq.csv': {'drop_columns': [],
                                'log2fc_col': 'log2FoldChange',
                                'padj_col': 'padj'}})),

    ({1: 'tests/big_counted.csv'}, {}, {},
     ({'tests/counted.csv': 'tests/counted.csv', 'tests/test_deseq.csv': 'tests/big_counted.csv',
       'tests/counted.tsv': 'tests/counted.tsv'},
      {'tests/counted.csv': 'Other', 'tests/test_deseq.csv': 'Other', 'tests/counted.tsv': 'Other'},
      {'tests/counted.csv': '', 'tests/test_deseq.csv': '', 'tests/counted.tsv': ''},
      {'tests/counted.csv': {'drop_columns': []},
       'tests/counted.tsv': {'drop_columns': []},
       'tests/test_deseq.csv': {'drop_columns': []}})),

    ({}, {0: 'Count matrix', 1: 'Differential expression'}, {1: 'new name', 2: 'second new name'},
     ({'tests/counted.csv': 'tests/counted.csv', 'tests/test_deseq.csv': 'tests/test_deseq.csv',
       'tests/counted.tsv': 'tests/counted.tsv'},
      {'tests/counted.csv': 'Count matrix', 'tests/test_deseq.csv': 'Differential expression',
       'tests/counted.tsv': 'Other'},
      {'tests/counted.csv': '', 'tests/test_deseq.csv': 'new name', 'tests/counted.tsv': 'second new name'},
      {'tests/counted.csv': {'drop_columns': [], 'is_normalized': False},
       'tests/counted.tsv': {'drop_columns': []},
       'tests/test_deseq.csv': {'drop_columns': [],
                                'log2fc_col': 'log2FoldChange',
                                'padj_col': 'padj'}})),
])
def test_MultiOpenWindow_result(qtbot, multi_open_window, path_ops, type_ops, name_ops, truth):
    files = multi_open_window_files
    for ind, op in path_ops.items():
        multi_open_window.paths[files[ind]].clear()
        qtbot.keyClicks(multi_open_window.paths[files[ind]].file_path, op)
    for ind, op in type_ops.items():
        qtbot.keyClicks(multi_open_window.table_types[files[ind]], op)
    for ind, op in name_ops.items():
        qtbot.keyClicks(multi_open_window.names[files[ind]], op)

    assert multi_open_window.result() == truth


def test_ReactiveTabWidget_init(qtbot, tab_widget):
    _ = tab_widget


def test_ReactiveTabWidget_new_tab_from_item(qtbot, tab_widget):
    with qtbot.waitSignal(tab_widget.newTabFromSet) as blocker:
        tab_widget.new_tab_from_item({'a', 'b', 'c', 'd'}, 'my name')
    assert blocker.args == [{'a', 'b', 'c', 'd'}, 'my name']

    with qtbot.waitSignal(tab_widget.newTabFromFilter) as blocker:
        tab_widget.new_tab_from_item(filtering.CountFilter('tests/test_files/counted.tsv'), 'my name')
    assert blocker.args == [filtering.CountFilter('tests/test_files/counted.tsv'), 'my name']

    with pytest.raises(TypeError):
        with qtbot.assertNotEmitted(tab_widget.newTabFromSet):
            with qtbot.assertNotEmitted(tab_widget.newTabFromFilter):
                tab_widget.new_tab_from_item(['invalid item type'], 'my name')


def test_ReactiveTabWidget_remove_tab(qtbot, tab_widget):
    qtbot, widget = widget_setup(qtbot, QtWidgets.QWidget)
    tab_widget.setCornerWidget(widget, QtCore.Qt.TopRightCorner)
    for i in range(3):
        qtbot, widget = widget_setup(qtbot, QtWidgets.QWidget)
        tab_widget.addTab(widget, 'name')

    for i in range(3):
        tab_widget.removeTab(0)

    assert tab_widget.count() == 0


def test_ReactiveTabWidget_left_click(qtbot, tab_widget, monkeypatch):
    super_triggered = []

    def mock_mouse_press_event(self, event):
        super_triggered.append(True)

    monkeypatch.setattr(QtWidgets.QTabWidget, 'mousePressEvent', mock_mouse_press_event)

    for i in range(3):
        qtbot.mouseClick(tab_widget, LEFT_CLICK)
        assert super_triggered == [True] * (i + 1)


def test_ReactiveTabWidget_right_click(qtbot, tab_widget, monkeypatch):
    super_triggered = []

    def mock_mouse_press_event(self, event):
        super_triggered.append(True)

    monkeypatch.setattr(QtWidgets.QTabWidget, 'mousePressEvent', mock_mouse_press_event)

    for i in range(3):
        # with qtbot.waitSignal(tab_widget.tabRightClicked) as blocker:
        qtbot.mouseClick(tab_widget, RIGHT_CLICK)
        assert super_triggered == []

    for i in range(3):
        qtbot, widget = widget_setup(qtbot, QtWidgets.QWidget)
        tab_widget.addTab(widget, 'name')

        for j in range(i + 1):
            with qtbot.waitSignal(tab_widget.tabRightClicked) as blocker:
                qtbot.mouseClick(tab_widget, RIGHT_CLICK, pos=tab_widget.tabBar().tabRect(j).center())
            assert super_triggered == []
            assert blocker.args[0] == j


def test_MainWindow_init(qtbot, main_window):
    _ = main_window


def test_MainWindow_add_new_tab(qtbot, main_window):
    main_window.new_table_action.trigger()
    assert main_window.tabs.count() == 2


@pytest.mark.parametrize('ind', range(6))
def test_MainWindow_add_new_tab_at(qtbot, main_window_with_tabs, ind):
    main_window_with_tabs.add_new_tab_at(ind)
    assert main_window_with_tabs.tabs.count() == 6
    assert main_window_with_tabs.tabs.currentIndex() == ind
    assert main_window_with_tabs.tabs.currentWidget().is_empty()


def test_MainWindow_close_current_tab(qtbot, main_window_with_tabs):
    current_tab_name = main_window_with_tabs.tabs.currentWidget().get_tab_name()
    main_window_with_tabs.close_current_action.trigger()
    assert main_window_with_tabs.tabs.count() == 4
    assert current_tab_name not in main_window_with_tabs.get_available_objects()
    assert main_window_with_tabs.tabs.currentWidget().get_tab_name() != current_tab_name


@pytest.mark.parametrize('ind', range(5))
def test_MainWindow_close_tabs_to_the_right(qtbot, main_window_with_tabs, ind):
    tab_name = main_window_with_tabs.tabs.widget(ind).get_tab_name()
    all_names = [main_window_with_tabs.tabs.widget(i).get_tab_name() for i in range(5)]
    main_window_with_tabs.close_tabs_to_the_right(ind)
    assert main_window_with_tabs.tabs.count() == ind + 1
    assert tab_name in main_window_with_tabs.get_available_objects()
    for i in range(ind, 5):
        if i <= ind:
            assert all_names[i] in main_window_with_tabs.get_available_objects()
        else:
            assert all_names[i] not in main_window_with_tabs.get_available_objects()


@pytest.mark.parametrize('ind', range(5))
def test_MainWindow_close_tabs_to_the_left(qtbot, main_window_with_tabs, ind):
    tab_name = main_window_with_tabs.tabs.widget(ind).get_tab_name()
    all_names = [main_window_with_tabs.tabs.widget(i).get_tab_name() for i in range(5)]
    main_window_with_tabs.close_tabs_to_the_left(ind)
    assert main_window_with_tabs.tabs.count() == 5 - ind
    assert tab_name in main_window_with_tabs.get_available_objects()
    for i in range(ind, 5):
        if i >= ind:
            assert all_names[i] in main_window_with_tabs.get_available_objects()
        else:
            assert all_names[i] not in main_window_with_tabs.get_available_objects()


@pytest.mark.parametrize('ind', range(5))
def test_MainWindow_close_other_tabs(qtbot, main_window_with_tabs, ind):
    tab_name = main_window_with_tabs.tabs.widget(ind).get_tab_name()
    all_names = [main_window_with_tabs.tabs.widget(i).get_tab_name() for i in range(5)]
    main_window_with_tabs.close_other_tabs(ind)
    assert main_window_with_tabs.tabs.count() == 1
    assert tab_name in main_window_with_tabs.get_available_objects()
    for i in range(ind, 5):
        if i != ind:
            assert all_names[i] not in main_window_with_tabs.get_available_objects()


@pytest.mark.parametrize('ind', range(5))
def test_MainWindow_close_tab(qtbot, main_window_with_tabs, ind):
    tab_name = main_window_with_tabs.tabs.widget(ind).get_tab_name()
    main_window_with_tabs.close_tab(ind)
    assert main_window_with_tabs.tabs.count() == 4
    assert tab_name not in main_window_with_tabs.get_available_objects()
    assert main_window_with_tabs.tabs.currentWidget().get_tab_name() != tab_name


def test_MainWindow_close_tab_undo(qtbot, main_window_with_tabs):
    current_tab_name = main_window_with_tabs.tabs.currentWidget().get_tab_name()
    current_obj = copy.copy(main_window_with_tabs.tabs.currentWidget().obj())
    main_window_with_tabs.close_current_action.trigger()
    assert main_window_with_tabs.tabs.count() == 4

    main_window_with_tabs.restore_tab_action.trigger()

    assert main_window_with_tabs.tabs.count() == 5
    assert current_tab_name in main_window_with_tabs.get_available_objects()
    assert main_window_with_tabs.tabs.currentWidget().get_tab_name() == current_tab_name
    assert main_window_with_tabs.tabs.currentWidget().obj() == current_obj


def test_MainWindow_sort_tabs_by_type(qtbot, main_window_with_tabs):
    truth = ['counted', 'test_deseq', 'my table', 'counted_6cols', 'majority_vote_intersection output']
    main_window_with_tabs.sort_tabs_by_type()
    for i, name in enumerate(truth):
        assert main_window_with_tabs.tabs.widget(i).get_tab_name() == name


def test_MainWindow_sort_tabs_by_n_features(qtbot, main_window_with_tabs):
    truth = ['test_deseq', 'my table', 'majority_vote_intersection output', 'counted_6cols', 'counted']
    main_window_with_tabs.sort_tabs_by_name()
    main_window_with_tabs.sort_tabs_by_n_features()
    for i, name in enumerate(truth):
        assert main_window_with_tabs.tabs.widget(i).get_tab_name() == name


def test_MainWindow_sort_tabs_by_creation_time(qtbot, main_window_with_tabs):
    truth = [main_window_with_tabs.tabs.widget(i).get_tab_name() for i in range(5)]
    main_window_with_tabs.sort_tabs_by_name()
    main_window_with_tabs.sort_tabs_by_creation_time()
    for i, name in enumerate(truth):
        assert main_window_with_tabs.tabs.widget(i).get_tab_name() == name


def test_MainWindow_sort_tabs_by_name(qtbot, main_window_with_tabs):
    truth = ['counted', 'counted_6cols', 'majority_vote_intersection output', 'my table', 'test_deseq']
    main_window_with_tabs.sort_tabs_by_name()
    for i, name in enumerate(truth):
        assert main_window_with_tabs.tabs.widget(i).get_tab_name() == name


def test_MainWindow_reverse_order(qtbot, main_window_with_tabs):
    truth = reversed(['counted', 'counted_6cols', 'majority_vote_intersection output', 'my table', 'test_deseq'])
    main_window_with_tabs.sort_tabs_by_name()
    main_window_with_tabs.sort_reverse()
    for i, name in enumerate(truth):
        assert main_window_with_tabs.tabs.widget(i).get_tab_name() == name


def test_MainWindow_rename_tab(qtbot, main_window_with_tabs):
    new_name = 'my new tab name'

    qtbot.keyClicks(main_window_with_tabs.tabs.currentWidget().overview_widgets['table_name'], new_name)
    qtbot.mouseClick(main_window_with_tabs.tabs.currentWidget().overview_widgets['rename_button'], LEFT_CLICK)

    assert main_window_with_tabs.tabs.currentWidget().get_tab_name() == new_name
    assert main_window_with_tabs.tabs.tabText(main_window_with_tabs.tabs.currentIndex()).rstrip('*') == new_name


@pytest.mark.parametrize('normalize', [False, True])
def test_MainWindow_new_table_from_folder(qtbot, main_window_with_tabs, normalize, monkeypatch):
    dir_path = 'tests/test_files/test_count_from_folder'

    def mock_get_dir(*args, **kwargs):
        return dir_path

    def mock_question(*args, **kwargs):
        if args[1] == 'Close program':
            return QtWidgets.QMessageBox.Yes
        return QtWidgets.QMessageBox.Yes if normalize else QtWidgets.QMessageBox.No

    monkeypatch.setattr(QtWidgets.QFileDialog, 'getExistingDirectory', mock_get_dir)
    monkeypatch.setattr(QtWidgets.QMessageBox, 'question', mock_question)

    main_window_with_tabs.new_table_from_folder_action.trigger()
    assert main_window_with_tabs.tabs.count() == 6
    assert main_window_with_tabs.tabs.currentWidget().obj() == filtering.CountFilter.from_folder(dir_path, normalize)


def test_MainWindow_multiple_new_tables(qtbot, main_window, monkeypatch):
    filenames = ['tests/test_files/test_deseq.csv', 'tests/test_files/counted.tsv', 'tests/test_files/fc_1.csv']
    objs_truth = [filtering.DESeqFilter(filenames[0]), filtering.CountFilter(filenames[1], drop_columns='cond2'),
                  filtering.FoldChangeFilter(filenames[2], 'num', 'denom')]
    objs_truth[1].fname = Path('new name')

    def mock_exec(self):
        return True

    def mock_multi_selection_result(self):
        return filenames

    def mock_multi_open_result(self):
        filename_dict = {fname: fname for fname in filenames}
        types_dict = {filenames[0]: 'Differential expression', filenames[1]: 'Count matrix',
                      filenames[2]: 'Fold change'}
        names_dict = {filenames[0]: '', filenames[1]: 'new name', filenames[2]: ''}
        kwargs_dict = {filenames[0]: {}, filenames[1]: {'drop_columns': ['cond2']},
                       filenames[2]: {'numerator_name': 'num', 'denominator_name': 'denom'}}
        return filename_dict, types_dict, names_dict, kwargs_dict

    monkeypatch.setattr(gui_windows.MultiFileSelectionDialog, 'exec', mock_exec)
    monkeypatch.setattr(MultiOpenWindow, 'exec', mock_exec)
    monkeypatch.setattr(gui_windows.MultiFileSelectionDialog, 'result', mock_multi_selection_result)
    monkeypatch.setattr(MultiOpenWindow, 'result', mock_multi_open_result)
    main_window.new_multiple_action.trigger()

    assert main_window.tabs.count() == 3
    for i in range(3):
        assert main_window.tabs.widget(i).obj() == objs_truth[i]


def test_MainWindow_export_pipeline(use_temp_settings_file, main_window, monkeypatch):
    fname = 'path/to/pipeline.yaml'
    pipeline_exported = []
    pipeline_truth = filtering.Pipeline.import_pipeline('tests/test_files/test_pipeline.yaml')

    def mock_export(pipeline, filename):
        assert filename == fname
        assert pipeline == pipeline_truth
        pipeline_exported.append(True)

    monkeypatch.setattr(QtWidgets.QInputDialog, 'getItem', lambda *args, **kwargs: ('test_pipeline', True))
    monkeypatch.setattr(QtWidgets.QFileDialog, 'getSaveFileName', lambda *args, **kwargs: (fname, '.yaml'))
    monkeypatch.setattr(filtering.Pipeline, 'export_pipeline', mock_export)
    main_window.pipelines['test_pipeline'] = pipeline_truth

    main_window.export_pipeline_action.trigger()
    assert pipeline_exported == [True]


def test_MainWindow_import_pipeline(use_temp_settings_file, main_window, monkeypatch):
    fname = 'tests/test_files/test_pipeline.yaml'
    monkeypatch.setattr(QtWidgets.QFileDialog, 'getOpenFileName', lambda *args, **kwargs: (fname, '.yaml'))
    main_window.import_pipeline_action.trigger()
    assert main_window.pipelines == {'test_pipeline': filtering.Pipeline.import_pipeline(fname)}


def test_MainWindow_import_multiple_gene_sets(qtbot, main_window_with_tabs, monkeypatch):
    filenames = ['tests/test_files/counted.tsv', 'tests/test_files/test_deseq.csv',
                 'tests/test_files/test_gene_set.txt']
    truth = []
    for filename in filenames:
        if filename.endswith('.txt'):
            with open(filename) as f:
                truth_set = set(f.read().split())
        else:
            df = io.load_csv(filename, index_col=0)
            truth_set = set(df.index)
        truth_featureset = enrichment.FeatureSet(truth_set, Path(filename).stem)
        truth.append(truth_featureset)

    def mock_get_files(*args, **kwargs):
        return filenames

    monkeypatch.setattr(gui_windows.MultiFileSelectionDialog, 'result', mock_get_files)
    monkeypatch.setattr(gui_windows.MultiFileSelectionDialog, 'exec', lambda *args, **kwargs: True)
    main_window_with_tabs.import_multiple_sets_action.trigger()
    assert main_window_with_tabs.tabs.count() == 5 + len(filenames)
    assert isinstance(main_window_with_tabs.tabs.currentWidget(), SetTabPage)
    for i in range(3):
        assert main_window_with_tabs.tabs.widget(i + 5).obj() == truth[i]


@pytest.mark.parametrize('filename', ['tests/test_files/counted.tsv', 'tests/test_files/test_deseq.csv',
                                      'tests/test_files/test_gene_set.txt'])
def test_MainWindow_import_gene_set(qtbot, main_window_with_tabs, monkeypatch, filename):
    if filename.endswith('.txt'):
        with open(filename) as f:
            truth_set = enrichment.FeatureSet(set(f.read().split()), Path(filename).stem)
    else:
        df = io.load_csv(filename, index_col=0)
        truth_set = enrichment.FeatureSet(set(df.index), Path(filename).stem)

    def mock_get_file(*args, **kwargs):
        return filename, '.csv'

    monkeypatch.setattr(QtWidgets.QFileDialog, 'getOpenFileName', mock_get_file)
    main_window_with_tabs.import_set_action.trigger()
    assert main_window_with_tabs.tabs.count() == 6
    assert isinstance(main_window_with_tabs.tabs.currentWidget(), SetTabPage)
    assert main_window_with_tabs.tabs.currentWidget().obj() == truth_set


@pytest.mark.parametrize('ind', range(5))
def test_MainWindow_export_gene_set(qtbot, use_temp_settings_file, main_window_with_tabs, monkeypatch, ind):
    save_path = 'test/save/path.csv'
    save_called = []

    def mock_save(gene_set, filename):
        assert gene_set == gene_set_truth
        assert filename == save_path
        save_called.append(True)

    def mock_get_file(*args, **kwargs):
        return save_path, '.csv'

    monkeypatch.setattr(io, 'save_gene_set', mock_save)
    monkeypatch.setattr(QtWidgets.QFileDialog, 'getSaveFileName', mock_get_file)
    main_window_with_tabs.tabs.setCurrentIndex(ind)
    obj = main_window_with_tabs.tabs.currentWidget().obj()
    if isinstance(obj, set):
        gene_set_truth = obj
    else:
        gene_set_truth = obj.index_set
    main_window_with_tabs.export_set_action.trigger()
    assert save_called == [True]


@pytest.mark.parametrize('ind,gene_set', [
    (4, {'WBGene00007069', 'WBGene00007064', 'WBGene00007063', 'WBGene00007074', 'WBGene00077502',
         'WBGene00007076', 'WBGene00044951', 'WBGene00007067', 'WBGene00044022', 'WBGene00043990',
         'WBGene00077504', 'WBGene00007066', 'WBGene00043987', 'WBGene00014997', 'WBGene00043989',
         'WBGene00007071', 'WBGene00007075', 'WBGene00007078', 'WBGene00007079', 'WBGene00007077',
         'WBGene00077503', 'WBGene00043988'}),
    (1, {'WBGene00007066', 'WBGene00007076', 'WBGene00044022', 'WBGene00007067', 'WBGene00043987',
         'WBGene00007077', 'WBGene00044951', 'WBGene00007075', 'WBGene00077502', 'WBGene00077504',
         'WBGene00007069', 'WBGene00007079', 'WBGene00043990', 'WBGene00043989', 'WBGene00014997',
         'WBGene00007074', 'WBGene00007071', 'WBGene00077503', 'WBGene00007063', 'WBGene00043988',
         'WBGene00007064', 'WBGene00007078'})
])
def test_MainWindow_copy_gene_set(qtbot, main_window_with_tabs, ind, gene_set):
    main_window_with_tabs.tabs.setCurrentIndex(1)
    main_window_with_tabs.copy_action.trigger()
    txt = QtWidgets.QApplication.clipboard().text()
    assert len(txt.split()) == len(gene_set)
    assert sorted(gene_set) == sorted(txt.split())


def test_MainWindow_add_pipeline(qtbot, main_window, monkeypatch):
    window_opened = []
    monkeypatch.setattr(CreatePipelineWindow, 'exec', functools.partial(window_opened.append, True))

    main_window.new_pipeline_action.trigger()
    assert window_opened == [True]


def test_MainWindow_apply_function(qtbot, main_window_with_tabs):
    main_window_with_tabs.choose_tab_by_name('test_deseq')
    tab = main_window_with_tabs.tabs.currentWidget()
    orig = filtering.DESeqFilter('tests/test_files/test_deseq.csv')

    truth = tab.obj().filter_significant(0.01, opposite=True, inplace=False)
    tab.stack_buttons[0].click()
    tab.stack.currentWidget().func_combo.setCurrentText(filtering.DESeqFilter.filter_significant.readable_name)
    tab.stack.currentWidget().parameter_widgets['alpha'].clear()
    qtbot.keyClicks(tab.stack.currentWidget().parameter_widgets['alpha'], '0.01')
    qtbot.mouseClick(tab.stack.currentWidget().parameter_widgets['opposite'].switch, LEFT_CLICK)
    qtbot.mouseClick(tab.stack.currentWidget().parameter_widgets['inplace'].switch, LEFT_CLICK)
    with qtbot.waitSignal(tab.filterObjectCreated, timeout=10000) as blocker:
        qtbot.mouseClick(tab.apply_button, LEFT_CLICK)
    assert blocker.args[0] == truth
    assert np.all(np.isclose(tab.obj().df, orig.df))

    assert main_window_with_tabs.tabs.count() == 6


def test_MainWindow_get_available_objects(qtbot, use_temp_settings_file, main_window_with_tabs):
    objs_truth = {'my table': filtering.FoldChangeFilter('tests/test_files/fc_1.csv', 'a', 'b'),
                  'counted': filtering.CountFilter('tests/test_files/counted.tsv'),
                  'counted_6cols': filtering.Filter('tests/test_files/counted_6cols.csv'),
                  'test_deseq': filtering.DESeqFilter('tests/test_files/test_deseq.csv'),
                  'majority_vote_intersection output': enrichment.FeatureSet(
                      {'WBGene00007069', 'WBGene00007064', 'WBGene00007063',
                       'WBGene00007074', 'WBGene00077502',
                       'WBGene00007076', 'WBGene00044951', 'WBGene00007067',
                       'WBGene00044022', 'WBGene00043990',
                       'WBGene00077504', 'WBGene00007066', 'WBGene00043987',
                       'WBGene00014997', 'WBGene00043989',
                       'WBGene00007071', 'WBGene00007075', 'WBGene00007078',
                       'WBGene00007079', 'WBGene00007077',
                       'WBGene00077503', 'WBGene00043988'}, 'majority_vote_intersection output')}
    objs_truth['my table'].fname = Path('my table')
    objs_truth['counted'].fname = Path('counted')
    objs_truth['counted_6cols'].fname = Path('counted_6cols')
    objs_truth['test_deseq'].fname = Path('test_deseq')

    res = main_window_with_tabs.get_available_objects()
    assert len(res) == len(objs_truth)
    for name in res.keys():
        assert isinstance(res[name][0], TabPage)
        assert (res[name][0].obj() == objs_truth[name]) or (
            np.all(np.isclose(np.squeeze(res[name][0].obj().df), np.squeeze(objs_truth[name].df))) and (
            res[name][0].obj().fname == objs_truth[name].fname))

        assert isinstance(res[name][1], QtGui.QIcon)


def test_MainWindow_choose_set_op(qtbot, use_temp_settings_file, main_window, monkeypatch):
    def mock_init(self, available_objs, parent=None):
        assert available_objs == 'my available objects'
        QtWidgets.QWidget.__init__(self)

    monkeypatch.setattr(main_window, 'get_available_objects', lambda: 'my available objects')
    window_opened = []
    monkeypatch.setattr(SetOperationWindow, '__init__', mock_init)
    monkeypatch.setattr(SetOperationWindow, 'show', functools.partial(window_opened.append, True))

    main_window.set_op_action.trigger()
    assert window_opened == [True]


def test_MainWindow_visualize_gene_sets(qtbot, use_temp_settings_file, main_window, monkeypatch):
    def mock_init(self, available_objs, parent=None):
        assert available_objs == 'my available objects'
        QtWidgets.QWidget.__init__(self)

    monkeypatch.setattr(main_window, 'get_available_objects', lambda: 'my available objects')
    window_opened = []
    monkeypatch.setattr(SetVisualizationWindow, '__init__', mock_init)
    monkeypatch.setattr(SetVisualizationWindow, 'show', functools.partial(window_opened.append, True))

    main_window.set_vis_action.trigger()
    assert window_opened == [True]


def test_MainWindow_open_enrichment_analysis(qtbot, main_window, monkeypatch):
    window_opened = []
    monkeypatch.setattr(EnrichmentWindow, 'show', functools.partial(window_opened.append, True))

    main_window.enrichment_action.trigger()
    assert window_opened == [True]


def test_MainWindow_choose_tab_by_name(qtbot, use_temp_settings_file, main_window_with_tabs, monkeypatch):
    truth = sorted(['my table', 'counted', 'counted_6cols', 'test_deseq', 'majority_vote_intersection output'])
    keys = truth.copy()
    np.random.shuffle(keys)
    main_window_with_tabs.sort_tabs_by_name()

    for key in keys:
        main_window_with_tabs.choose_tab_by_name(key)
        assert main_window_with_tabs.tabs.currentIndex() == truth.index(key)


def test_MainWindow_save_session(qtbot, use_temp_settings_file, main_window, monkeypatch):
    monkeypatch.setattr(QtWidgets.QFileDialog, 'getOpenFileName',
                        lambda *args, **kwargs: ('tests/test_files/test_session.rnal', '.rnal'))
    main_window.load_session()

    session_fname = 'session filename.rnal'
    n_files = main_window.tabs.count()
    n_pipelines = len(main_window.pipelines)
    item_types_truth = [filtering.FoldChangeFilter, filtering.CountFilter, filtering.Filter, filtering.DESeqFilter,
                        set]
    item_properties_truth = [{'numerator_name': 'a', 'denominator_name': 'b'}, {'is_normalized': False}, {},
                             {'log2fc_col': 'log2FoldChange', 'padj_col': 'padj'}, {}]
    pipeline_names_truth = ['New Pipeline', 'Other Pipeline']
    pipeline_files_truth = [pipeline.export_pipeline(filename=None) for pipeline in main_window.pipelines.values()]
    func_called = []

    def mock_save_session(session_filename: Union[str, Path], file_names: List[str], item_names: List[str],
                          item_types: list, item_properties: list, pipeline_names: List[str],
                          pipeline_files: List[str]):
        assert session_filename == session_fname
        assert len(file_names) == n_files
        assert item_types == item_types_truth
        assert item_properties == item_properties_truth

        assert len(pipeline_names) == n_pipelines
        assert len(pipeline_files) == n_pipelines
        assert pipeline_names == pipeline_names_truth
        assert pipeline_files == pipeline_files_truth

        func_called.append(True)

    monkeypatch.setattr(io, 'save_gui_session', mock_save_session)
    monkeypatch.setattr(QtWidgets.QFileDialog, 'getSaveFileName', lambda *args, **kwargs: (session_fname, '.rnal'))

    main_window.save_session_action.trigger()
    assert func_called == [True]


def test_MainWindow_load_session(qtbot, use_temp_settings_file, main_window, monkeypatch):
    pipelines_truth = {'New Pipeline': filtering.Pipeline('Filter'),
                       'Other Pipeline': filtering.Pipeline('DESeqFilter')}
    pipelines_truth['New Pipeline'].add_function('filter_top_n', by='log2FoldChange', n=99, ascending=True,
                                                 na_position='last', opposite=False)
    pipelines_truth['New Pipeline'].add_function('describe', percentiles=[0.01, 0.25, 0.5, 0.75, 0.99])
    pipelines_truth['Other Pipeline'].add_function('split_fold_change_direction')
    pipelines_truth['Other Pipeline'].add_function('volcano_plot', alpha=0.1)

    objs_truth = [filtering.FoldChangeFilter('tests/test_files/fc_1.csv', 'a', 'b'),
                  filtering.CountFilter('tests/test_files/counted.tsv'),
                  filtering.Filter('tests/test_files/counted_6cols.csv'),
                  filtering.DESeqFilter('tests/test_files/test_deseq.csv'),
                  enrichment.FeatureSet(
                      {'WBGene00007069', 'WBGene00007064', 'WBGene00007063', 'WBGene00007074', 'WBGene00077502',
                       'WBGene00007076', 'WBGene00044951', 'WBGene00007067', 'WBGene00044022', 'WBGene00043990',
                       'WBGene00077504', 'WBGene00007066', 'WBGene00043987', 'WBGene00014997', 'WBGene00043989',
                       'WBGene00007071', 'WBGene00007075', 'WBGene00007078', 'WBGene00007079', 'WBGene00007077',
                       'WBGene00077503', 'WBGene00043988'}, 'majority_vote_intersection output')]
    objs_truth[0].fname = Path('my table')
    objs_truth[1].fname = Path('counted')
    objs_truth[2].fname = Path('counted_6cols')
    objs_truth[3].fname = Path('test_deseq')

    monkeypatch.setattr(QtWidgets.QFileDialog, 'getOpenFileName',
                        lambda *args, **kwargs: ('tests/test_files/test_session.rnal', '.rnal'))
    main_window.load_session_action.trigger()
    assert main_window.tabs.count() == 5
    assert len(main_window.pipelines) == 2
    assert main_window.pipelines == pipelines_truth

    for i in range(1, main_window.tabs.count()):
        assert (main_window.tabs.widget(i).obj() == objs_truth[i]) or (
            np.all(np.isclose(main_window.tabs.widget(i).obj().df, objs_truth[i].df)) and (
            main_window.tabs.widget(i).obj().fname == objs_truth[i].fname))


def test_MainWindow_about(qtbot, main_window, monkeypatch):
    window_opened = []
    monkeypatch.setattr(gui_windows.AboutWindow, 'exec', functools.partial(window_opened.append, True))

    main_window.about_action.trigger()
    assert window_opened == [True]


def test_MainWindow_cite(qtbot, main_window, monkeypatch):
    window_opened = []
    monkeypatch.setattr(gui_windows.HowToCiteWindow, 'exec', functools.partial(window_opened.append, True))

    main_window.cite_action.trigger()
    assert window_opened == [True]


def test_MainWindow_settings(qtbot, main_window, monkeypatch):
    window_opened = []
    monkeypatch.setattr(gui_windows.SettingsWindow, 'exec', functools.partial(window_opened.append, True))

    main_window.settings_action.trigger()
    assert window_opened == [True]


def test_MainWindow_user_guide(qtbot, main_window, monkeypatch):
    window_opened = []

    def mock_open_url(url):
        window_opened.append(True)
        return True

    monkeypatch.setattr(QtGui.QDesktopServices, 'openUrl', mock_open_url)

    main_window.user_guide_action.trigger()
    assert window_opened == [True]


def test_MainWindow_context_menu(qtbot, main_window_with_tabs, monkeypatch):
    opened = []

    def mock_exec(*args, **kwargs):
        opened.append(True)

    monkeypatch.setattr(QtWidgets.QMenu, 'exec', mock_exec)
    qtbot.mouseClick(main_window_with_tabs.tabs.tabBar(), RIGHT_CLICK,
                     pos=main_window_with_tabs.tabs.tabBar().tabRect(1).center())
    assert opened == [True]


def test_MainWindow_clear_history(qtbot, main_window_with_tabs, monkeypatch):
    cleared = [False for i in range(main_window_with_tabs.tabs.count())]
    truth = [True for i in cleared]

    def mock_clear(*args, ind):
        cleared[ind] = True

    for i in range(len(truth)):
        monkeypatch.setattr(main_window_with_tabs.undo_group.stacks()[i], 'clear', functools.partial(mock_clear, ind=i))

    main_window_with_tabs.clear_history()
    assert cleared == truth


def test_MainWindow_clear_session(qtbot, main_window_with_tabs):
    main_window_with_tabs.clear_session(confirm_action=False)
    assert main_window_with_tabs.tabs.count() == 1
    assert main_window_with_tabs.tabs.widget(0).is_empty()


@pytest.mark.parametrize('action_name', ['ontology_graph_action', 'pathway_graph_action', 'featurecounts_single_action',
                                         'featurecounts_paired_action', 'bowtie2_index_action',
                                         'bowtie2_single_action', 'bowtie2_paired_action', 'kallisto_index_action',
                                         'kallisto_single_action', 'kallisto_paired_action', 'cutadapt_single_action',
                                         'cutadapt_paired_action', 'set_op_action', 'enrichment_action',
                                         'set_vis_action', 'bar_plot_action'])
def test_MainWindow_open_windows(qtbot, main_window_with_tabs, action_name):
    action = getattr(main_window_with_tabs, action_name)
    action.trigger()


@pytest.mark.parametrize('action_name, window_attr_name',
                         [('new_pipeline_action', 'pipeline_window'),
                          ('cite_action', 'cite_window'),
                          ('about_action', 'about_window'),
                          ('settings_action', 'settings_window')])
def test_MainWindow_open_dialogs(qtbot, main_window_with_tabs, action_name, window_attr_name, monkeypatch):
    action = getattr(main_window_with_tabs, action_name)

    def win():
        return getattr(main_window_with_tabs, window_attr_name)

    def handle_dialog():
        while win() is None or not win().isVisible():
            QtWidgets.QApplication.processEvents()
        win().close()

    QtCore.QTimer.singleShot(100, handle_dialog)
    action.trigger()
