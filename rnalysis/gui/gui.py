import builtins
import copy
import functools
import hashlib
import importlib
import itertools
import os
import platform
import sys
import time
import typing
import warnings
from collections import OrderedDict
from pathlib import Path
from queue import Queue
from typing import List, Tuple, Union, Callable, Type

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import yaml
from PyQt5 import QtCore, QtWidgets, QtGui
from joblib import parallel_backend

from rnalysis import fastq, filtering, enrichment, __version__
from rnalysis.gui import gui_style, gui_widgets, gui_windows, gui_graphics, gui_quickstart
from rnalysis.utils import io, validation, generic, parsing, settings, clustering

FILTER_OBJ_TYPES = {'Count matrix': filtering.CountFilter, 'Differential expression': filtering.DESeqFilter,
                    'Fold change': filtering.FoldChangeFilter, 'Other table': filtering.Filter}
FILTER_OBJ_TYPES_INV = {val.__name__: key for key, val in FILTER_OBJ_TYPES.items()}
INIT_EXCLUDED_PARAMS = {'self', 'fname', 'suppress_warnings'}

FROZEN_ENV = getattr(sys, 'frozen', False) and hasattr(sys, '_MEIPASS')

JOB_COUNTER = gui_widgets.JobCounter()


def check_run_success(result: gui_widgets.WorkerOutput):
    if result.raised_exception:
        raise result.raised_exception


class BarPlotWindow(gui_windows.FuncExternalWindow):
    EXCLUDED_PARAMS = set()
    __slots__ = {}

    def __init__(self, parent=None):
        func = enrichment.enrichment_bar_plot
        help_link = f"https://guyteichman.github.io/RNAlysis/build/rnalysis.enrichment.{func.__name__}.html"
        super().__init__('Enrichment bar-plot', func, help_link, self.EXCLUDED_PARAMS, threaded=False, parent=parent)
        self.init_ui()
        self.setWindowTitle('Create enrichment bar-plot')


class OntologyGraphWindow(gui_windows.FuncExternalWindow):
    EXCLUDED_PARAMS = set()
    __slots__ = {}

    def __init__(self, parent=None):
        func = enrichment.gene_ontology_graph
        help_link = f"https://guyteichman.github.io/RNAlysis/build/rnalysis.fastq.{func.__name__}.html"
        super().__init__('Gene Ontology graph', func, help_link, self.EXCLUDED_PARAMS, threaded=False, parent=parent)
        self.init_ui()
        self.setWindowTitle('Plot Gene Ontology Graph')


class PathwayGraphWindow(gui_windows.FuncExternalWindow):
    EXCLUDED_PARAMS = set()
    __slots__ = {}

    def __init__(self, parent=None):
        func = enrichment.kegg_pathway_graph
        help_link = f"https://guyteichman.github.io/RNAlysis/build/rnalysis.fastq.{func.__name__}.html"
        super().__init__('KEGG Pathway graph', func, help_link, self.EXCLUDED_PARAMS, threaded=False, parent=parent)
        self.init_ui()
        self.setWindowTitle('Plot KEGG Pathway Graph')


class FeatureCountsSingleWindow(gui_windows.FuncExternalWindow):
    EXCLUDED_PARAMS = set()
    __slots__ = {}

    def __init__(self, parent=None):
        func = fastq.featurecounts_single_end
        help_link = f"https://guyteichman.github.io/RNAlysis/build/rnalysis.fastq.{func.__name__}.html"
        super().__init__(func.readable_name, func, help_link, self.EXCLUDED_PARAMS, parent=parent)
        self.init_ui()
        self.setWindowTitle('featureCounts single-end counting setup')


class SamToFastqSingleWindow(gui_windows.FuncExternalWindow):
    EXCLUDED_PARAMS = set()
    __slots__ = {}

    def __init__(self, parent=None):
        func = fastq.sam_to_fastq_single
        help_link = f"https://guyteichman.github.io/RNAlysis/build/rnalysis.fastq.{func.__name__}.html"
        super().__init__(func.readable_name, func, help_link, self.EXCLUDED_PARAMS, parent=parent)
        self.init_ui()
        self.setWindowTitle('Convert SAM/BAM to FASTQ (single-end)')


class SamToFastqPairedWindow(gui_windows.FuncExternalWindow):
    EXCLUDED_PARAMS = set()
    __slots__ = {}

    def __init__(self, parent=None):
        func = fastq.sam_to_fastq_paired
        help_link = f"https://guyteichman.github.io/RNAlysis/build/rnalysis.fastq.{func.__name__}.html"
        super().__init__(func.readable_name, func, help_link, self.EXCLUDED_PARAMS, parent=parent)
        self.init_ui()
        self.setWindowTitle('Convert SAM/BAM to FASTQ (paired-end)')


class FastqToSamSingleWindow(gui_windows.FuncExternalWindow):
    EXCLUDED_PARAMS = {'self', 'return_new_filenames', 'legacy_args'}
    __slots__ = {}

    def __init__(self, parent=None):
        func = fastq.fastq_to_sam_single
        help_link = f"https://guyteichman.github.io/RNAlysis/build/rnalysis.fastq.{func.__name__}.html"

        super().__init__(func.readable_name, func, help_link, self.EXCLUDED_PARAMS, parent=parent)
        self.init_ui()
        self.setWindowTitle('Convert FASTQ to SAM/BAM (single-end)')


class FastqToSamPairedWindow(gui_windows.PairedFuncExternalWindow):
    __slots__ = {}

    def __init__(self, parent=None):
        func = fastq.fastq_to_sam_paired
        help_link = f"https://guyteichman.github.io/RNAlysis/build/rnalysis.fastq.{func.__name__}.html"

        super().__init__(func.readable_name, func, help_link, self.EXCLUDED_PARAMS, parent=parent)
        self.init_ui()
        self.setWindowTitle('Convert FASTQ to SAM/BAM (paired-end)')


class ConvertSamFormatWindow(gui_windows.FuncExternalWindow):
    EXCLUDED_PARAMS = set()
    __slots__ = {}

    def __init__(self, parent=None):
        func = fastq.convert_sam_format
        help_link = f"https://guyteichman.github.io/RNAlysis/build/rnalysis.fastq.{func.__name__}.html"
        super().__init__(func.readable_name, func, help_link, self.EXCLUDED_PARAMS, parent=parent)
        self.init_ui()
        self.setWindowTitle('Convert SAM/BAM format')


class BamIndexWindow(gui_windows.FuncExternalWindow):
    EXCLUDED_PARAMS = set()
    __slots__ = {}

    def __init__(self, parent=None):
        func = fastq.create_bam_index
        help_link = f"https://guyteichman.github.io/RNAlysis/build/rnalysis.fastq.{func.__name__}.html"

        super().__init__('Create BAM index', func, help_link, self.EXCLUDED_PARAMS, parent=parent)
        self.init_ui()
        self.setWindowTitle('Create BAM index')


class ValidateSamWindow(gui_windows.FuncExternalWindow):
    EXCLUDED_PARAMS = set()
    __slots__ = {}

    def __init__(self, parent=None):
        func = fastq.validate_sam
        help_link = f"https://guyteichman.github.io/RNAlysis/build/rnalysis.fastq.{func.__name__}.html"

        super().__init__('Validate SAM/BAM file', func, help_link, self.EXCLUDED_PARAMS, parent=parent)
        self.init_ui()
        self.setWindowTitle('Validate SAM/BAM file')


class SortSamWindow(gui_windows.FuncExternalWindow):
    EXCLUDED_PARAMS = set()
    __slots__ = {}

    def __init__(self, parent=None):
        func = fastq.sort_sam
        help_link = f"https://guyteichman.github.io/RNAlysis/build/rnalysis.fastq.{func.__name__}.html"

        super().__init__('Sort SAM/BAM file', func, help_link, self.EXCLUDED_PARAMS, parent=parent)
        self.init_ui()
        self.setWindowTitle('Sort SAM/BAM file')


class FindDuplicatesWindow(gui_windows.FuncExternalWindow):
    EXCLUDED_PARAMS = set()
    __slots__ = {}

    def __init__(self, parent=None):
        func = fastq.find_duplicates
        help_link = f"https://guyteichman.github.io/RNAlysis/build/rnalysis.fastq.{func.__name__}.html"

        super().__init__('Find duplicate reads', func, help_link, self.EXCLUDED_PARAMS, parent=parent)
        self.init_ui()
        self.setWindowTitle('Find duplicate reads')


class FeatureCountsPairedWindow(gui_windows.FuncExternalWindow):
    EXCLUDED_PARAMS = set()
    __slots__ = {}

    def __init__(self, parent=None):
        func = fastq.featurecounts_paired_end
        help_link = f"https://guyteichman.github.io/RNAlysis/build/rnalysis.fastq.{func.__name__}.html"
        super().__init__(func.readable_name, func, help_link, self.EXCLUDED_PARAMS, parent=parent)
        self.init_ui()
        self.setWindowTitle('featureCounts paired-end counting setup')


class Bowtie2IndexWindow(gui_windows.FuncExternalWindow):
    EXCLUDED_PARAMS = set()
    __slots__ = {}

    def __init__(self, parent=None):
        func = fastq.bowtie2_create_index
        help_link = f"https://guyteichman.github.io/RNAlysis/build/rnalysis.fastq.{func.__name__}.html"
        super().__init__(func.readable_name, func, help_link, self.EXCLUDED_PARAMS, parent=parent)
        self.init_ui()
        self.setWindowTitle('Bowtie2 - build genome index')


class Bowtie2SingleWindow(gui_windows.FuncExternalWindow):
    EXCLUDED_PARAMS = set()
    __slots__ = {}

    def __init__(self, parent=None):
        func = fastq.bowtie2_align_single_end
        help_link = f"https://guyteichman.github.io/RNAlysis/build/rnalysis.fastq.{func.__name__}.html"
        super().__init__(func.readable_name, func, help_link, self.EXCLUDED_PARAMS, parent=parent)
        self.init_ui()
        self.setWindowTitle('Bowtie2 single-end alignment setup')


class Bowtie2PairedWindow(gui_windows.PairedFuncExternalWindow):
    __slots__ = {}

    def __init__(self, parent=None):
        func = fastq.bowtie2_align_paired_end
        help_link = f"https://guyteichman.github.io/RNAlysis/build/rnalysis.fastq.{func.__name__}.html"
        super().__init__(func.readable_name, func, help_link, self.EXCLUDED_PARAMS, parent)
        self.init_ui()
        self.setWindowTitle('Bowtie2 paired-end alignment setup')


class ShortStackWindow(gui_windows.FuncExternalWindow):
    EXCLUDED_PARAMS = set()
    __slots__ = {}

    def __init__(self, parent=None):
        func = fastq.shortstack_align_smallrna
        help_link = f"https://guyteichman.github.io/RNAlysis/build/rnalysis.fastq.{func.__name__}.html"
        super().__init__(func.readable_name, func, help_link, self.EXCLUDED_PARAMS, parent=parent)
        self.init_ui()
        self.setWindowTitle('ShortStack small RNA alignment setup')


class KallistoIndexWindow(gui_windows.FuncExternalWindow):
    EXCLUDED_PARAMS = set()

    __slots__ = {}

    def __init__(self, parent=None):
        func = fastq.kallisto_create_index
        help_link = f"https://guyteichman.github.io/RNAlysis/build/rnalysis.fastq.{func.__name__}.html"
        super().__init__(func.readable_name, func, help_link, self.EXCLUDED_PARAMS, parent=parent)
        self.init_ui()
        self.setWindowTitle('Kallisto - build transcriptome index')


class KallistoSingleWindow(gui_windows.FuncExternalWindow):
    EXCLUDED_PARAMS = {'legacy_args'}
    __slots__ = {}

    def __init__(self, parent=None):
        func = fastq.kallisto_quantify_single_end
        help_link = f"https://guyteichman.github.io/RNAlysis/build/rnalysis.fastq.{func.__name__}.html"
        super().__init__(func.readable_name, func, help_link, self.EXCLUDED_PARAMS, parent=parent)
        self.init_ui()
        self.setWindowTitle('Kallisto single-end quantification setup')


class KallistoPairedWindow(gui_windows.PairedFuncExternalWindow):
    EXCLUDED_PARAMS = gui_windows.PairedFuncExternalWindow.EXCLUDED_PARAMS.copy()
    EXCLUDED_PARAMS.add('legacy_args')
    __slots__ = {}

    def __init__(self, parent=None):
        func = fastq.kallisto_quantify_paired_end
        help_link = f"https://guyteichman.github.io/RNAlysis/build/rnalysis.fastq.{func.__name__}.html"

        super().__init__(func.readable_name, func, help_link, self.EXCLUDED_PARAMS,
                         parent)
        self.init_ui()
        self.setWindowTitle('Kallisto paired-end quantification setup')


class CutAdaptSingleWindow(gui_windows.FuncExternalWindow):
    EXCLUDED_PARAMS = set()
    __slots__ = {}

    def __init__(self, parent=None):
        func = fastq.trim_adapters_single_end
        help_link = f"https://guyteichman.github.io/RNAlysis/build/rnalysis.fastq.{func.__name__}.html"
        super().__init__(func.readable_name, func, help_link, self.EXCLUDED_PARAMS, parent=parent)
        self.init_ui()
        self.setWindowTitle('CutAdapt single-end adapter trimming setup')


class CutAdaptPairedWindow(gui_windows.PairedFuncExternalWindow):
    EXCLUDED_PARAMS = gui_windows.PairedFuncExternalWindow.EXCLUDED_PARAMS.copy()
    EXCLUDED_PARAMS.add('return_new_filenames')

    __slots__ = {}

    def __init__(self, parent=None):
        func = fastq.trim_adapters_paired_end
        help_link = f"https://guyteichman.github.io/RNAlysis/build/rnalysis.fastq.{func.__name__}.html"
        super().__init__(func.readable_name, func, help_link, self.EXCLUDED_PARAMS, parent)
        self.init_ui()
        self.setWindowTitle('CutAdapt paired-end adapter trimming setup')


class DiffExpWindow(gui_windows.FuncExternalWindow):
    EXCLUDED_PARAMS = {'self', 'comparisons'}
    IGNORED_WIDGETS = gui_windows.FuncExternalWindow.IGNORED_WIDGETS | {'load_design'}

    __slots__ = {'comparisons': 'list of comparisons to make',
                 'design_mat': 'design matrix',
                 'comparisons_group': 'widget group for choosing comparisons',
                 'comparisons_grid': 'layout for choosing comparisons',
                 'comparisons_widgets': 'widgets for choosing comparisons'}

    def __init__(self, func: Callable, name: str, parent=None):
        self.name = name
        help_link = f"https://guyteichman.github.io/RNAlysis/build/rnalysis.filtering.CountFilter.{func.__name__}.html"
        super().__init__(name, func, help_link, self.EXCLUDED_PARAMS, parent=parent)

        self.comparisons = []
        self.design_mat = None

        self.comparisons_group = QtWidgets.QGroupBox(f"2. Choose pairwise comparisons for {name}")
        self.comparisons_grid = QtWidgets.QGridLayout(self.comparisons_group)
        self.comparisons_widgets = {}

        self.init_ui()

    def init_ui(self):
        self.setWindowTitle(f'{self.name} differential expression setup')
        self.scroll_layout.addWidget(self.comparisons_group, 0, 1)
        super().init_ui()
        _, _, width, height = self.scroll.geometry().getRect()
        self.resize(1100, height)

    def init_param_ui(self):
        super().init_param_ui()
        self.param_widgets['load_design'] = QtWidgets.QPushButton('Load design matrix')
        self.param_widgets['load_design'].setEnabled(False)
        self.param_widgets['load_design'].clicked.connect(self.init_comparisons_ui)
        self.param_grid.addWidget(self.param_widgets['load_design'], self.param_grid.rowCount(), 0, 1, 2)
        self.param_widgets['design_matrix'].textChanged.connect(self._change_accept_button_state)

    def _change_accept_button_state(self, is_legal: bool):
        self.param_widgets['load_design'].setEnabled(is_legal)

    def init_comparisons_ui(self):
        design_mat = io.load_table(self.param_widgets['design_matrix'].text(), index_col=0)
        for factor in design_mat.columns:
            assert parsing.slugify(factor) == factor, f"Invalid factor name '{factor}': contains invalid characters." \
                                                      f" \nSuggested alternative name: '{parsing.slugify(factor)}'. "
        self.design_mat = design_mat
        if 'picker' in self.comparisons_widgets:
            self.comparisons_grid.removeWidget(self.comparisons_widgets['picker'])
            self.comparisons_widgets['picker'].deleteLater()

        self.comparisons_widgets['picker'] = gui_widgets.ComparisonPickerGroup(self.design_mat, self)
        self.comparisons_grid.addWidget(self.comparisons_widgets['picker'], 0, 0)

    def get_analysis_kwargs(self):
        kwargs = super().get_analysis_kwargs()
        kwargs['comparisons'] = self.comparisons_widgets['picker'].get_comparison_values()
        return kwargs


class DESeqWindow(DiffExpWindow):
    def __init__(self, parent=None):
        func = filtering.CountFilter.differential_expression_deseq2
        name = 'DESeq2'
        super().__init__(func, name, parent)


class LimmaWindow(DiffExpWindow):
    def __init__(self, parent=None):
        func = filtering.CountFilter.differential_expression_limma_voom
        name = 'Limma-Voom'
        super().__init__(func, name, parent)


class ClicomWindow(gui_windows.FuncExternalWindow):
    EXCLUDED_PARAMS = {'self', 'parameter_dicts', 'gui_mode'}
    ADDITIONAL_EXCLUDED_PARAMS = {'power_transform', 'plot_style', 'split_plots', 'return_probabilities', 'gui_mode',
                                  'parallel_backend'}
    if FROZEN_ENV:
        EXCLUDED_PARAMS.add('parallel_backend')

    def __init__(self, funcs: dict, filter_obj: filtering.Filter, parent=None):
        func = filtering.CountFilter.split_clicom
        help_link = f"https://guyteichman.github.io/RNAlysis/build/rnalysis.filtering.CountFilter.{func.__name__}.html"
        super().__init__('CLICOM', func, help_link, self.EXCLUDED_PARAMS, parent=parent)
        self.parameter_dicts: List[dict] = []
        self.funcs = funcs
        self.setups_counter = {key: 0 for key in self.funcs.keys()}
        self.filter_obj = filter_obj
        self.stack = FuncTypeStack(self.funcs, self.filter_obj, self,
                                   additional_excluded_params=self.ADDITIONAL_EXCLUDED_PARAMS)

        self.setups_group = QtWidgets.QGroupBox("2. Choose clustering setups for CLICOM")
        self.setups_grid = QtWidgets.QGridLayout(self.setups_group)
        self.setups_widgets = {}

        self.init_ui()

    def connect_widget(self, widget: QtWidgets.QWidget):
        super().connect_widget(widget)
        if isinstance(widget, (gui_widgets.TableColumnPicker, gui_widgets.TableColumnPicker)):
            widget.add_columns(self.filter_obj.columns)

    def init_ui(self):
        super().init_ui()
        self.setWindowTitle('CLICOM clustering setup')
        self.setMinimumWidth(1500)
        self.scroll_layout.addWidget(self.setups_group, 0, 1)
        self.init_setups_ui()

    def init_setups_ui(self):
        self.setups_grid.addWidget(self.stack, 1, 0)
        self.setups_widgets['list'] = gui_widgets.MultiChoiceListWithDelete(list(), parent=self.setups_group)
        self.setups_widgets['list'].itemDeleted.connect(self.remove_clustering_setup)
        self.setups_grid.addWidget(QtWidgets.QLabel('<b>Added setups</b>'), 0, 1, QtCore.Qt.AlignCenter)
        self.setups_grid.addWidget(self.setups_widgets['list'], 1, 1, 2, 1)

        self.setups_widgets['add_button'] = QtWidgets.QPushButton('Add setup')
        self.setups_widgets['add_button'].clicked.connect(self.add_clustering_setup)
        self.stack.funcSelected.connect(self.setups_widgets['add_button'].setEnabled)
        self.setups_grid.addWidget(self.setups_widgets['add_button'], 2, 0)
        self.setups_widgets['add_button'].setDisabled(True)
        self.setups_grid.setRowStretch(1, 1)

    @QtCore.pyqtSlot(int)
    def remove_clustering_setup(self, ind: int):
        self.parameter_dicts.pop(ind)

    def add_clustering_setup(self):
        func_name = self.stack.get_function_name()
        func_params = self.stack.get_function_params()
        func_params['method'] = func_name.lstrip('split_').lower()
        self.parameter_dicts.append(func_params)
        self.setups_counter[func_name] += 1
        self.setups_widgets['list'].add_item(f"{func_params['method']}_{self.setups_counter[func_name]}")

    def get_analysis_args(self):
        return self.parameter_dicts


class EnrichmentWindow(gui_widgets.MinMaxDialog):
    EXCLUDED_PARAMS = {'self', 'save_csv', 'fname', 'return_fig', 'biotype', 'background_genes',
                       'statistical_test', 'parametric_test', 'biotype_ref_path', 'gui_mode'}
    if FROZEN_ENV:
        EXCLUDED_PARAMS.add('parallel_backend')

    ANALYSIS_TYPES = {'Gene Ontology (GO)': 'go',
                      'Kyoto Encyclopedia of Genes and Genomes (KEGG)': 'kegg',
                      'Categorical attributes': 'user_defined', 'Non-categorical attributes': 'non_categorical'}

    ANALYSIS_TYPES_BUTTONS = (('External datasets:', ('Gene Ontology (GO)',
                                                      'Kyoto Encyclopedia of Genes and Genomes (KEGG)')),
                              ('Custom dataset:', ('Categorical attributes', 'Non-categorical attributes')))

    ANALYSIS_FUNCS = {('go', False): enrichment.FeatureSet.go_enrichment,
                      ('go', True): enrichment.RankedSet.single_set_go_enrichment,
                      ('kegg', False): enrichment.FeatureSet.kegg_enrichment,
                      ('kegg', True): enrichment.RankedSet.single_set_kegg_enrichment,
                      ('user_defined', False): enrichment.FeatureSet.user_defined_enrichment,
                      ('user_defined', True): enrichment.RankedSet.single_set_enrichment,
                      ('non_categorical', False): enrichment.FeatureSet.non_categorical_enrichment}

    STATISTICAL_TESTS = {'Randomization test': 'randomization', "Fisher's Exact test": 'fisher',
                         'Hypergeometric test': 'hypergeometric', 'Single-set enrichment (XL-mHG test)': 'single_set'}

    ORDINAL_STATISTICAL_TESTS = {'One-sample T-test (parametric)': True, 'Sign test (non-parametric)': False}

    STATISTICAL_TEST_ARGS = {'randomization': {'alpha', 'randomization_reps', 'random_seed'},
                             'fisher': {'alpha'},
                             'hypergeometric': {'alpha'},
                             'single_set': {'alpha', 'min_positive_genes', 'lowest_cutoff'},
                             'non_categorical': {'alpha'},
                             None: {},
                             True: {'alpha'},
                             False: {'alpha'}}

    PLOT_ARGS = {'user_defined': {'plot_horizontal', 'plot_style', 'show_expected'},
                 'go': {'plot_horizontal', 'plot_style', 'show_expected', 'plot_ontology_graph',
                        'ontology_graph_format'},
                 'kegg': {'plot_horizontal', 'plot_style', 'show_expected', 'plot_pathway_graphs',
                          'pathway_graphs_format'},
                 'non_categorical': {'plot_log_scale', 'plot_style', 'n_bins'}}

    enrichmentStarted = QtCore.pyqtSignal(object, object)
    geneSetsRequested = QtCore.pyqtSignal(object)
    functionApplied = QtCore.pyqtSignal(gui_widgets.WorkerOutput)

    def __init__(self, parent=None):
        super().__init__(parent)

        self.parameters_signature = {}
        self.stats_signature = {}
        self.plot_signature = {}

        self.scroll = QtWidgets.QScrollArea()

        self.widgets = {}
        self.list_group = QtWidgets.QGroupBox('Enrichment analysis', self)
        self.list_grid = QtWidgets.QGridLayout(self.list_group)

        self.parameter_group = QtWidgets.QGroupBox('Additional parameters', self)
        self.parameter_grid = QtWidgets.QGridLayout(self.parameter_group)
        self.parameter_widgets = {}

        self.stats_group = QtWidgets.QGroupBox('Statistical test', self)
        self.stats_grid = QtWidgets.QGridLayout(self.stats_group)
        self.stats_widgets = {}

        self.plot_group = QtWidgets.QGroupBox('Configure plot', self)
        self.plot_grid = QtWidgets.QGridLayout(self.plot_group)
        self.plot_widgets = {}

        self.scroll_widget = QtWidgets.QWidget(self.scroll)
        self.main_layout = QtWidgets.QVBoxLayout(self)
        self.scroll_layout = QtWidgets.QVBoxLayout(self.scroll_widget)

        self.init_basic_ui()

    def init_basic_ui(self):
        self.setWindowTitle('Enrichment Analysis')
        self.main_layout.addWidget(self.scroll)

        self.parameter_group.setVisible(False)
        self.plot_group.setVisible(False)
        self.stats_group.setVisible(False)

        self.scroll.setHorizontalScrollBarPolicy(QtCore.Qt.ScrollBarAlwaysOff)
        self.scroll.setWidgetResizable(True)
        self.scroll.setWidget(self.scroll_widget)
        self.scroll_layout.setSizeConstraint(QtWidgets.QLayout.SetMinAndMaxSize)

        self.scroll_layout.addWidget(self.list_group)
        self.scroll_layout.addWidget(self.stats_group)
        self.scroll_layout.addWidget(self.parameter_group)
        self.scroll_layout.addWidget(self.plot_group)

        self.widgets['run_button'] = QtWidgets.QPushButton('Run')
        self.widgets['run_button'].clicked.connect(self.run_analysis)
        self.widgets['run_button'].setVisible(False)
        self.scroll_layout.addWidget(self.widgets['run_button'])

        self.scroll_layout.addStretch(1)

        for lst in ['enrichment_list', 'bg_list']:
            self.widgets[lst] = gui_widgets.GeneSetComboBox(self)
            self.widgets[lst].boxOpened.connect(functools.partial(self.geneSetsRequested.emit, self.widgets[lst]))
            self.widgets[lst].currentTextChanged.connect(self._verify_inputs)

        self.widgets['dataset_radiobox'] = gui_widgets.RadioButtonBox('Choose enrichment dataset',
                                                                      self.ANALYSIS_TYPES_BUTTONS)
        self.widgets['dataset_radiobox'].buttonClicked.connect(self.update_uis)

        self.list_grid.addWidget(self.widgets['dataset_radiobox'], 0, 0, 6, 1)
        self.list_grid.addWidget(self.widgets['enrichment_list'], 1, 1)
        self.list_grid.addWidget(self.widgets['bg_list'], 3, 1)
        self.list_grid.addWidget(QtWidgets.QLabel('<b>Choose enrichment set:</b>', self), 0, 1)
        self.list_grid.addWidget(QtWidgets.QLabel('<b>Choose background set:</b>', self), 2, 1)

        self.scroll.setMinimumWidth(self.scroll_widget.sizeHint().width() + 150)

    def _set_background_select_mode(self, selectable: bool = True):
        if selectable:
            self.widgets['bg_list'].setDisabled(False)
        else:
            self.widgets['bg_list'].setDisabled(True)

    def _get_statistical_test_name(self):
        try:
            button = self.stats_widgets['stats_radiobox'].checkedButton()
        except KeyError:
            button = None

        if button is None:
            return None
        return button.text()

    def _get_statistical_test(self):
        stat_test = self._get_statistical_test_name()
        if stat_test is None:
            return None

        try:
            if self.is_categorical():
                return self.STATISTICAL_TESTS[stat_test]
            return self.ORDINAL_STATISTICAL_TESTS[stat_test]
        except KeyError:
            return None

    def _update_single_set(self):
        if self.is_single_set():
            self._set_background_select_mode(False)
        else:
            self._set_background_select_mode(True)

    def _update_signatures(self, update_params: bool = True, update_stats: bool = True, update_plot: bool = True):
        if update_params:
            self.parameters_signature = {}
        if update_stats:
            self.stats_signature = {}
        if update_plot:
            self.plot_signature = {}

        analysis_type = self.get_current_analysis_type()
        chosen_func = self.get_current_func()
        signature = generic.get_method_signature(chosen_func)
        func_desc, param_desc = io.get_method_docstring(chosen_func)
        for name, param in signature.items():
            this_desc = param_desc.get(name, '')
            if name in self.EXCLUDED_PARAMS:
                continue
            elif update_plot and name in self.PLOT_ARGS[analysis_type]:
                self.plot_signature[name] = (param, this_desc)
            elif update_stats and name in set.union(*self.STATISTICAL_TEST_ARGS.values()):
                self.stats_signature[name] = (param, this_desc)
            elif update_params:
                self.parameters_signature[name] = (param, this_desc)

    def update_uis(self):
        self._update_signatures()

        self.update_parameters_ui()
        self.update_stats_ui()
        self.update_plot_ui()
        self.widgets['run_button'].setVisible(True)
        self.widgets['run_button'].setDisabled(True)
        self._verify_inputs()

        if 'help_link' in self.widgets:
            self.scroll_layout.removeWidget(self.widgets['help_link'])
            self.widgets['help_link'].deleteLater()
        chosen_func_name = self.get_current_func().__name__
        obj_type = enrichment.RankedSet if self.is_single_set() else enrichment.FeatureSet
        help_address = f"https://guyteichman.github.io/RNAlysis/build/rnalysis.enrichment." \
                       f"{obj_type.__name__}.{chosen_func_name}.html"
        self.widgets['help_link'] = QtWidgets.QLabel(f'<a href="{help_address}">Open documentation for function '
                                                     f'<b>{obj_type.__name__}.{chosen_func_name}</b></a>')
        self.widgets['help_link'].setOpenExternalLinks(True)
        self.scroll_layout.insertWidget(4, self.widgets['help_link'])

        _, _, width, height = self.scroll.geometry().getRect()
        self.resize(width, 750)

    @QtCore.pyqtSlot()
    def _verify_inputs(self):
        is_legal = True
        if not self.widgets['enrichment_list'].is_legal():
            is_legal = False
        if self.widgets['bg_list'].isEnabled() and (not self.widgets['bg_list'].is_legal()):
            is_legal = False
        if self._get_statistical_test() is None:
            is_legal = False
        self.widgets['run_button'].setEnabled(is_legal)

    def get_current_analysis_type(self):
        button = self.widgets['dataset_radiobox'].checkedButton()
        if button is None:
            return None
        return self.ANALYSIS_TYPES[self.widgets['dataset_radiobox'].checkedButton().text()]

    def get_current_func(self):
        single_set = self.is_single_set()
        name = self.get_current_analysis_type()
        func = self.ANALYSIS_FUNCS[(name, single_set)]
        return func

    def update_plot_ui(self):
        self.plot_group.setVisible(True)
        # delete previous widgets
        self.plot_widgets = {}
        gui_widgets.clear_layout(self.plot_grid)

        i = 0
        for name, (param, desc) in self.plot_signature.items():
            self.plot_widgets[name] = gui_widgets.param_to_widget(param, name)
            label = QtWidgets.QLabel(f'{name}:', self.plot_widgets[name])
            label.setToolTip(desc)
            help_button = gui_widgets.HelpButton()
            self.plot_grid.addWidget(label, i, 0)
            self.plot_grid.addWidget(self.plot_widgets[name], i, 1)
            self.plot_grid.addWidget(help_button, i, 2)
            help_button.set_param_help(name, desc)

            i += 1

    def update_parameters_ui(self):
        self.parameter_group.setVisible(True)
        # delete previous widgets
        self.parameter_widgets = {}
        gui_widgets.clear_layout(self.parameter_grid)

        i = 0
        for name, (param, desc) in self.parameters_signature.items():
            self.parameter_widgets[name] = gui_widgets.param_to_widget(param, name)
            label = QtWidgets.QLabel(f'{name}:', self.parameter_widgets[name])
            label.setToolTip(desc)
            help_button = gui_widgets.HelpButton()
            self.parameter_grid.addWidget(help_button, i, 2)
            self.parameter_grid.addWidget(label, i, 0)
            self.parameter_grid.addWidget(self.parameter_widgets[name], i, 1)
            help_button.set_param_help(name, desc)
            i += 1

    def update_stats_ui(self):
        self._update_signatures(False, True, False)
        self.stats_group.setVisible(True)
        prev_test = self._get_statistical_test()
        prev_test_name = self._get_statistical_test_name()
        self._update_single_set()

        # delete previous widgets
        self.stats_widgets = {}
        gui_widgets.clear_layout(self.stats_grid)

        radio_options = parsing.data_to_list(self.STATISTICAL_TESTS.keys() if self.is_categorical() else \
                                                 self.ORDINAL_STATISTICAL_TESTS.keys())
        self.stats_widgets['stats_radiobox'] = gui_widgets.RadioButtonBox('Choose statistical test:', radio_options)
        self.stats_widgets['stats_radiobox'].selectionChanged.connect(self._verify_inputs)
        self.stats_widgets['stats_radiobox'].buttonClicked.connect(self._verify_inputs)
        if prev_test_name is not None:
            self.stats_widgets['stats_radiobox'].set_selection(prev_test_name)
        self.stats_widgets['stats_radiobox'].buttonClicked.connect(self.update_stats_ui)

        self.stats_grid.addWidget(self.stats_widgets['stats_radiobox'], 0, 0, 3, 2)

        i = 0
        for name, (param, desc) in self.stats_signature.items():
            if name in self.STATISTICAL_TEST_ARGS[prev_test]:
                self.stats_widgets[name] = gui_widgets.param_to_widget(param, name)
                label = QtWidgets.QLabel(f'{name}:', self.stats_widgets[name])
                label.setToolTip(desc)
                help_button = gui_widgets.HelpButton()
                self.stats_grid.addWidget(help_button, i, 4)
                self.stats_grid.addWidget(label, i, 2)
                self.stats_grid.addWidget(self.stats_widgets[name], i, 3)
                help_button.set_param_help(name, desc)
                i += 1

    def is_single_set(self):
        stat_test = self._get_statistical_test()
        return stat_test == 'single_set'

    def is_categorical(self):
        analysis_type = self.get_current_analysis_type()
        return analysis_type != 'non_categorical'

    def get_analysis_params(self):
        kwargs = dict()
        pred_ids = []
        stat_test = self._get_statistical_test()

        if not self.is_single_set():
            if self.is_categorical():
                kwargs['statistical_test'] = stat_test
            else:
                kwargs['parametric_test'] = stat_test
        gene_set = self.widgets['enrichment_list'].value()
        pred_ids.append(self.widgets['enrichment_list'].current_id())
        gene_set_name = self.widgets['enrichment_list'].currentText()

        for param_name, widget in itertools.chain(self.parameter_widgets.items(), self.plot_widgets.items(),
                                                  self.stats_widgets.items()):
            if param_name in {'help_link', 'dataset_radiobox', 'stats_radiobox'}:
                continue
            val = gui_widgets.get_val_from_widget(widget)
            kwargs[param_name] = val

        if not self.is_single_set():
            bg_set = self.widgets['bg_list'].value()
            pred_ids.append(self.widgets['bg_list'].current_id())

        else:
            bg_set = None
        return gene_set, bg_set, gene_set_name.rstrip('*'), kwargs, pred_ids

    @QtCore.pyqtSlot()
    def run_analysis(self):
        func = self.get_current_func()
        gene_set, bg_set, set_name, kwargs, pred_ids = self.get_analysis_params()
        kwargs['gui_mode'] = True
        print("Enrichment analysis started")
        self.showMinimized()
        try:
            is_single_set = self.is_single_set()
            if is_single_set:
                bg_set_obj = None
            else:
                bg_set_obj = enrichment.FeatureSet(bg_set, 'background_set')
            feature_set_obj = enrichment.RankedSet(gene_set, set_name) if is_single_set \
                else enrichment.FeatureSet(gene_set, set_name)

            if is_single_set:
                partial = functools.partial(func, feature_set_obj, **kwargs)
            else:
                partial = functools.partial(func, feature_set_obj, background_genes=bg_set_obj, **kwargs)

            worker = gui_widgets.Worker(partial, JOB_COUNTER.get_id(), pred_ids, set_name)
            worker.finished.connect(self.functionApplied.emit)

            self.enrichmentStarted.emit(worker, self.showNormal)

        except Exception as e:
            self.showNormal()
            raise e


class SetOperationWindow(gui_widgets.MinMaxDialog):
    SET_OPERATIONS = {'Union': 'union', 'Majority-Vote Intersection': 'majority_vote_intersection',
                      'Intersection': 'intersection', 'Difference': 'difference',
                      'Symmetric Difference': 'symmetric_difference', 'Other': 'other'}
    EXCLUDED_PARAMS = {'self', 'other', 'others', 'return_type'}

    geneSetReturned = QtCore.pyqtSignal(set, str, list, dict)
    primarySetUsed = QtCore.pyqtSignal(str)
    primarySetChangedDifference = QtCore.pyqtSignal(str)
    primarySetChangedIntersection = QtCore.pyqtSignal()

    def __init__(self, available_objects: dict, parent=None):
        super().__init__(parent)
        self.available_objects = available_objects
        self.widgets = {}
        self.list_group = QtWidgets.QGroupBox('Choose gene sets', self)
        self.list_grid = QtWidgets.QGridLayout(self.list_group)
        self.parameter_group = QtWidgets.QGroupBox('Additional parameters', self)
        self.parameter_grid = QtWidgets.QGridLayout(self.parameter_group)
        self.parameter_widgets = {}
        self.operations_group = QtWidgets.QGroupBox('Set operation')
        self.operations_grid = QtWidgets.QGridLayout(self.operations_group)
        self.layout = QtWidgets.QHBoxLayout(self)

        self.init_ui()

    def create_canvas(self):
        set_names = [item.text() for item in self.widgets['set_list'].get_sorted_selection()]
        sets = [self.available_objects[name][0].obj() for name in set_names]
        ind = 0
        while ind < len(set_names):
            if sets[ind] is None or self.available_objects[set_names[ind]][0].is_empty():
                set_names.pop(ind)
                sets.pop(ind)
            else:
                ind += 1

        if len(set_names) < 2:
            canvas = gui_graphics.EmptyCanvas('Please select 2 or more gene sets to continue', self)
        else:
            items = {}
            for s, s_name in zip(sets, set_names):
                if validation.isinstanceinh(s, filtering.Filter):
                    s_set = s.index_set
                elif validation.isinstanceinh(s, enrichment.FeatureSet):
                    s_set = s.gene_set
                elif isinstance(s, set):
                    s_set = s
                else:
                    raise TypeError(type(s))
                items[s_name] = s_set
            if 2 <= len(set_names) <= 3:
                canvas = gui_graphics.VennInteractiveCanvas(items, self)

            else:
                canvas = gui_graphics.UpSetInteractiveCanvas(items, self)
        if 'canvas' in self.widgets:
            self.widgets['canvas'].deleteLater()
            self.widgets['toolbar'].deleteLater()
            self.operations_grid.removeWidget(self.widgets['canvas'])
            self.operations_grid.removeWidget(self.widgets['toolbar'])

        self.widgets['canvas'] = canvas
        self.widgets['toolbar'] = gui_graphics.CleanPlotToolBar(self.widgets['canvas'], self)
        self.operations_grid.addWidget(self.widgets['canvas'], 1, 2, 3, 3)
        self.operations_grid.addWidget(self.widgets['toolbar'], 0, 2, 1, 3)

        for col in range(2, 5):
            self.operations_grid.setColumnStretch(col, 1)
        for row in range(1, 4):
            self.operations_grid.setRowStretch(row, 1)

        if not isinstance(canvas, gui_graphics.EmptyCanvas):
            self._connect_canvas(canvas)

    def _connect_canvas(self, canvas: gui_graphics.BaseInteractiveCanvas):
        canvas.manualChoice.connect(self._set_op_other)
        self.widgets['radio_button_box'].radio_buttons['Union'].clicked.connect(canvas.union)
        self.widgets['radio_button_box'].radio_buttons['Intersection'].clicked.connect(canvas.intersection)
        self.primarySetChangedDifference.connect(canvas.difference)
        self.primarySetChangedIntersection.connect(canvas.intersection)

        if isinstance(canvas, gui_graphics.VennInteractiveCanvas):
            self.widgets['radio_button_box'].radio_buttons['Symmetric Difference'].clicked.connect(
                canvas.symmetric_difference)
        this_button = self.widgets['radio_button_box'].checkedButton()
        if this_button is not None:
            self.widgets['radio_button_box'].set_selection(this_button.text())

    def init_ui(self):
        self.setWindowTitle('Set Operations')
        self.setGeometry(600, 50, 1050, 800)
        self.setLayout(self.layout)
        self.widgets['splitter'] = QtWidgets.QSplitter(QtCore.Qt.Horizontal)
        self.layout.addWidget(self.widgets['splitter'])
        self.widgets['splitter'].addWidget(self.list_group)
        self.widgets['splitter'].addWidget(self.operations_group)
        self.widgets['splitter'].setSizes([int(self.width() * 0.15), int(self.width() * 0.85)])
        self.parameter_group.setVisible(False)

        self.init_sets_ui()
        self.init_operations_ui()

    def init_sets_ui(self):
        self.widgets['set_list'] = gui_widgets.MultiChoiceListWithReorder(self.available_objects,
                                                                          [val[1] for val in
                                                                           self.available_objects.values()],
                                                                          self)

        for func in [self.create_canvas, self._check_legal_operations, self._validate_input,
                     self._toggle_choose_primary_set]:
            self.widgets['set_list'].itemSelectionChanged.connect(func)
            self.widgets['set_list'].itemOrderChanged.connect(func)
        self.list_grid.addWidget(self.widgets['set_list'], 0, 0)

    def _toggle_choose_primary_set(self):
        if self.get_current_func_name() in ['difference', 'intersection']:
            self.widgets['choose_primary_set'].setVisible(True)
            self.widgets['choose_primary_set_label'].setVisible(True)

            self.widgets['choose_primary_set'].clear()
            self.widgets['choose_primary_set'].addItems(
                [item.text() for item in self.widgets['set_list'].get_sorted_selection()])
            self.widgets['canvas'].clear_selection()
        else:
            self.widgets['choose_primary_set'].setVisible(False)
            self.widgets['choose_primary_set_label'].setVisible(False)

    def init_operations_ui(self):
        self.widgets['set_op_box'] = QtWidgets.QWidget(self)
        self.widgets['set_op_box_layout'] = QtWidgets.QVBoxLayout(self.widgets['set_op_box'])
        self.operations_grid.addWidget((self.widgets['set_op_box']), 1, 0, 3, 1)
        self.widgets['radio_button_box'] = gui_widgets.RadioButtonBox('Choose set operation',
                                                                      self.SET_OPERATIONS.keys())

        for func in [self.update_paremeter_ui, self._validate_input, self._toggle_choose_primary_set]:
            self.widgets['radio_button_box'].buttonClicked.connect(func)
            self.widgets['radio_button_box'].selectionChanged.connect(func)
        self.widgets['radio_button_box'].radio_buttons['Majority-Vote Intersection'].clicked.connect(
            self._majority_vote_intersection)
        self.widgets['set_op_box_layout'].addWidget(self.widgets['radio_button_box'], stretch=1)

        self.widgets['choose_primary_set'] = gui_widgets.MandatoryComboBox('Choose primary set...', self)
        self.widgets['choose_primary_set'].currentTextChanged.connect(self._primary_set_changed)
        self.widgets['choose_primary_set_label'] = QtWidgets.QLabel('Primary set for operation:')
        self.widgets['set_op_box_layout'].addWidget(self.widgets['choose_primary_set_label'])
        self.widgets['set_op_box_layout'].addWidget(self.widgets['choose_primary_set'])

        self.widgets['set_op_box_layout'].addWidget(self.parameter_group)

        self.widgets['set_op_box_layout'].addStretch(1)

        self._toggle_choose_primary_set()

        self.create_canvas()

        self.widgets['apply_button'] = QtWidgets.QPushButton('Apply')
        self.widgets['apply_button'].clicked.connect(self.apply_set_op)
        self.widgets['apply_button'].setEnabled(False)
        self.operations_grid.addWidget(self.widgets['apply_button'], 4, 0, 1, 6)

    def _majority_vote_intersection(self):
        if not isinstance(self.widgets['canvas'], gui_graphics.EmptyCanvas):
            if 'majority_threshold' in self.parameter_widgets:
                threshold = gui_widgets.get_val_from_widget(self.parameter_widgets['majority_threshold'])
            else:
                threshold = 0.5
            self.widgets['canvas'].majority_vote_intersection(threshold)

    @QtCore.pyqtSlot(str)
    def _primary_set_changed(self, set_name: str):
        func_name = self.get_current_func_name()
        if func_name == 'difference':
            self.primarySetChangedDifference.emit(set_name)
        elif func_name == 'intersection':
            self.primarySetChangedIntersection.emit()
        self._validate_input()

    def update_paremeter_ui(self):
        # delete previous widgets
        try:
            self.parameter_widgets['help_link'].deleteLater()
            self.operations_grid.removeWidget(self.parameter_widgets['help_link'])
        except KeyError:
            pass
        self.parameter_widgets = {}
        gui_widgets.clear_layout(self.parameter_grid)

        chosen_func_name = self.get_current_func_name()
        signature = generic.get_method_signature(chosen_func_name, filtering.Filter)
        i = 0
        for name, param in signature.items():
            if name in self.EXCLUDED_PARAMS:
                continue
            self.parameter_widgets[name] = gui_widgets.param_to_widget(param, name)
            self.parameter_grid.addWidget(QtWidgets.QLabel(f'{name}:', self.parameter_widgets[name]), i, 0)
            self.parameter_grid.addWidget(self.parameter_widgets[name], i, 1)
            if chosen_func_name == 'majority_vote_intersection':
                self.parameter_widgets[name].valueChanged.connect(self._majority_vote_intersection)
                self._majority_vote_intersection()
            i += 1

        self.parameter_group.setVisible(i > 0)

        if chosen_func_name != 'other':
            help_address = f"https://guyteichman.github.io/RNAlysis/build/rnalysis.filtering." \
                           f"{filtering.Filter.__name__}.{chosen_func_name}.html"
            self.parameter_widgets['help_link'] = QtWidgets.QLabel(
                f'<a href="{help_address}">Open documentation for function '
                f'<b>{filtering.Filter.__name__}.{chosen_func_name}</b></a>', self)
            self.parameter_widgets['help_link'].setOpenExternalLinks(True)
            self.operations_grid.addWidget(self.parameter_widgets['help_link'], 5, 0, 1, 6)

    def get_current_func_name(self):
        button = self.widgets['radio_button_box'].checkedButton()
        if button is None:
            return None
        return self.SET_OPERATIONS[button.text()]

    def _check_legal_operations(self):
        n_items = len(self.widgets['set_list'].get_sorted_selection())
        if self.get_current_func_name() == 'symmetric_difference' and n_items > 2:
            self.widgets['radio_button_box'].set_selection('Other')
        sym_diff_button = self.widgets['radio_button_box'].radio_buttons['Symmetric Difference']
        sym_diff_button.setEnabled(n_items <= 2)

    def _validate_input(self):
        is_legal = True

        if isinstance(self.widgets['canvas'], gui_graphics.EmptyCanvas):
            is_legal = False

        if self.get_current_func_name() is None:
            is_legal = False

        if self.widgets['choose_primary_set'].isVisible() and not self.widgets['choose_primary_set'].is_legal():
            is_legal = False

        self.widgets['apply_button'].setEnabled(is_legal)

    def _set_op_other(self):
        self.widgets['radio_button_box'].set_selection('Other')

    def _get_set_names(self):
        return [item.text() for item in self.widgets['set_list'].get_sorted_selection()]

    def _get_function_params(self):
        set_names = self._get_set_names()
        if self.get_current_func_name() in ['intersection', 'difference']:
            primary_set_name = self.widgets['choose_primary_set'].currentText()
            self.primarySetUsed.emit(primary_set_name)
        else:
            primary_set_name = set_names[0]
        kwargs = {}
        for param_name, widget in self.parameter_widgets.items():
            if param_name in {'apply_button', 'help_link'}:
                continue
            val = gui_widgets.get_val_from_widget(widget)

            kwargs[param_name] = val
        return set_names, primary_set_name, kwargs

    @QtCore.pyqtSlot()
    def apply_set_op(self):
        func_name = self.get_current_func_name()
        ancestor_ids = [self.available_objects[name][0].tab_id for name in self._get_set_names()]
        if func_name == 'other':
            output_set = self.widgets['canvas'].get_custom_selection()
            output_name = "Other set operation output"
            kwargs = {}
        else:
            set_names, primary_set_name, kwargs = self._get_function_params()

            first_obj = self.available_objects[primary_set_name][0].obj()
            if isinstance(first_obj, set):
                first_obj = filtering.Filter.from_dataframe(pd.DataFrame(index=list(first_obj)), 'placeholder')
            other_objs = []
            for name in set_names:
                if name != primary_set_name:
                    other_objs.append(self.available_objects[name][0].obj())
            output_name = f"{func_name} output"

            if kwargs.get('inplace', False):
                command = SetOpInplacCommand(self.available_objects[primary_set_name][0], func_name, other_objs, kwargs,
                                             f'Apply "{func_name}"', ancestor_ids)
                self.available_objects[primary_set_name][0].undo_stack.push(command)
                output_set = None
            else:
                output_set = getattr(first_obj, func_name)(*other_objs, **kwargs)

        if output_set is not None:
            self.geneSetReturned.emit(output_set, output_name, ancestor_ids, kwargs)
        self.close()


class SetVisualizationWindow(gui_widgets.MinMaxDialog):
    VISUALIZATION_FUNCS = {'Venn Diagram': 'venn_diagram', 'UpSet Plot': 'upset_plot'}
    EXCLUDED_PARAMS = {'objs', 'attr_ref_table_path', 'fig'}
    figureGenerated = QtCore.pyqtSignal(plt.Figure, str, list, dict)

    def __init__(self, available_objects: dict, parent=None):
        super().__init__(parent)
        self.available_objects = available_objects
        self.layout = QtWidgets.QHBoxLayout(self)
        self.widgets = {}

        self.list_group = QtWidgets.QGroupBox('Choose gene sets', self)
        self.list_grid = QtWidgets.QGridLayout(self.list_group)

        self.visualization_group = QtWidgets.QGroupBox('Gene set visualization')
        self.visualization_grid = QtWidgets.QGridLayout(self.visualization_group)

        self.parameter_widgets = {}
        self.parameter_group = QtWidgets.QGroupBox('Additional parameters')
        self.parameter_grid = QtWidgets.QGridLayout(self.parameter_group)

        self.init_ui()

    def init_ui(self):
        self.setWindowTitle('Gene Set Visualization')
        self.setGeometry(600, 50, 1050, 800)
        self.setLayout(self.layout)
        self.widgets['splitter'] = QtWidgets.QSplitter(QtCore.Qt.Horizontal)
        self.layout.addWidget(self.widgets['splitter'])
        self.widgets['splitter'].addWidget(self.list_group)
        self.widgets['splitter'].addWidget(self.visualization_group)
        self.widgets['splitter'].setSizes([int(self.width() * 0.1), int(self.width() * 0.9)])

        self.parameter_group.setVisible(False)

        self.init_list_ui()
        self.init_visualization_ui()

    def init_list_ui(self):
        self.widgets['set_list'] = gui_widgets.MultiChoiceListWithReorder(self.available_objects,
                                                                          [val[1] for val in
                                                                           self.available_objects.values()],
                                                                          self)

        for func in [self._check_legal_operations, self._validate_input, self.create_canvas]:
            self.widgets['set_list'].itemSelectionChanged.connect(func)
            self.widgets['set_list'].itemOrderChanged.connect(func)

        self.list_grid.addWidget(self.widgets['set_list'], 0, 0)

    def init_visualization_ui(self):
        self.widgets['radio_button_box'] = gui_widgets.RadioButtonBox('Choose visualization type:',
                                                                      self.VISUALIZATION_FUNCS, parent=self)
        for func in [self.update_parameter_ui, self._validate_input, self.create_canvas]:
            self.widgets['radio_button_box'].buttonClicked.connect(func)
            self.widgets['radio_button_box'].selectionChanged.connect(func)

        self.visualization_grid.addWidget(self.widgets['radio_button_box'], 0, 0, 2, 1)
        self.visualization_grid.addWidget(self.parameter_group, 2, 0, 2, 1)
        self.visualization_grid.setRowStretch(self.visualization_grid.count(), 1)

        self.widgets['generate_button'] = QtWidgets.QPushButton('Generate graph')
        self.widgets['generate_button'].clicked.connect(self.generate_graph)
        self.widgets['generate_button'].setEnabled(False)
        self.visualization_grid.addWidget(self.widgets['generate_button'], 4, 0, 1, 5)

        self.create_canvas()

    def create_canvas(self):
        set_names = [item.text() for item in self.widgets['set_list'].get_sorted_selection()]
        sets = [self.available_objects[name][0].obj() for name in set_names]
        ind = 0
        while ind < len(set_names):
            if sets[ind] is None or self.available_objects[set_names[ind]][0].is_empty():
                set_names.pop(ind)
                sets.pop(ind)
            else:
                ind += 1

        func_name = self.get_current_func_name()

        if len(set_names) < 2:
            canvas = gui_graphics.EmptyCanvas('Please select 2 or more gene sets to continue', self)
        elif func_name is None:
            canvas = gui_graphics.EmptyCanvas("Please choose a visualization function to continue", self)
        else:
            objs_to_plot, kwargs, _ = self._get_function_params()
            try:
                canvas = gui_graphics.BasePreviewCanvas(getattr(enrichment, func_name), self, objs=objs_to_plot,
                                                        **kwargs)
            except Exception:
                canvas = gui_graphics.EmptyCanvas("Invalid input; please change one or more of your parameters")
        if 'canvas' in self.widgets:
            self.widgets['canvas'].deleteLater()
            self.visualization_grid.removeWidget(self.widgets['canvas'])

        self.widgets['canvas'] = canvas
        self.visualization_grid.addWidget(self.widgets['canvas'], 0, 2, 4, 3)

        for col in range(2, self.visualization_grid.columnCount()):
            self.visualization_grid.setColumnStretch(col, 2)
        for row in range(0, 4):
            self.visualization_grid.setRowStretch(row, 1)

    def _validate_input(self):
        is_legal = True

        if self.get_current_func_name() is None:
            is_legal = False

        if len(self.widgets['set_list'].get_sorted_selection()) < 2:
            is_legal = False

        self.widgets['generate_button'].setEnabled(is_legal)

    def _check_legal_operations(self):
        n_items = len(self.widgets['set_list'].get_sorted_selection())
        if self.get_current_func_name() == 'venn_diagram' and n_items > 3:
            self.widgets['radio_button_box'].set_selection('UpSet Plot')
        venn_button = self.widgets['radio_button_box'].radio_buttons['Venn Diagram']
        venn_button.setEnabled(n_items <= 3)

    def update_parameter_ui(self):
        # delete previous widgets
        try:
            self.parameter_widgets['help_link'].deleteLater()
            self.visualization_grid.removeWidget(self.parameter_widgets['help_link'])
        except KeyError:
            pass
        self.parameter_widgets = {}
        gui_widgets.clear_layout(self.parameter_grid)

        chosen_func_name = self.get_current_func_name()
        signature = generic.get_method_signature(chosen_func_name, enrichment)
        i = 0
        for name, param in signature.items():
            if name in self.EXCLUDED_PARAMS:
                continue
            self.parameter_widgets[name] = gui_widgets.param_to_widget(param, name,
                                                                       actions_to_connect=self.create_canvas)
            self.parameter_grid.addWidget(QtWidgets.QLabel(f'{name}:', self.parameter_widgets[name]), i, 0)
            self.parameter_grid.addWidget(self.parameter_widgets[name], i, 1)
            i += 1

        help_address = f"https://guyteichman.github.io/RNAlysis/build/rnalysis.enrichment.{chosen_func_name}.html"
        self.parameter_widgets['help_link'] = QtWidgets.QLabel(
            f'<a href="{help_address}">Open documentation for function '
            f'<b>enrichment.{chosen_func_name}</b></a>', self)
        self.parameter_widgets['help_link'].setOpenExternalLinks(True)
        self.visualization_grid.addWidget(self.parameter_widgets['help_link'], 5, 0, 1, 4)

        self.parameter_group.setVisible(i > 0)

    def get_current_func_name(self):
        button = self.widgets['radio_button_box'].checkedButton()
        if button is None:
            return None
        return self.VISUALIZATION_FUNCS[button.text()]

    def _get_function_params(self):
        set_names = [item.text() for item in self.widgets['set_list'].get_sorted_selection()]
        objs_to_plot = {name: self.available_objects[name][0].obj() for name in set_names if
                        not self.available_objects[name][0].is_empty()}
        ancestor_ids = [self.available_objects[name][0].tab_id for name in set_names if
                        not self.available_objects[name][0].is_empty()]

        kwargs = {}
        for param_name, widget in self.parameter_widgets.items():
            if param_name in {'apply_button', 'help_link'}:
                continue
            val = gui_widgets.get_val_from_widget(widget)

            kwargs[param_name] = val
        return objs_to_plot, kwargs, ancestor_ids

    @QtCore.pyqtSlot()
    def generate_graph(self):
        func_name = self.get_current_func_name()
        objs_to_plot, kwargs, ancestor_ids = self._get_function_params()
        fig = getattr(enrichment, func_name)(objs_to_plot, **kwargs)
        self.figureGenerated.emit(fig, func_name, ancestor_ids, kwargs)


class TabPage(QtWidgets.QWidget):
    functionApplied = QtCore.pyqtSignal(gui_widgets.WorkerOutput)
    itemSpawned = QtCore.pyqtSignal(str, int, int, object)
    filterObjectCreated = QtCore.pyqtSignal(object, int)
    featureSetCreated = QtCore.pyqtSignal(object, int)
    startedJob = QtCore.pyqtSignal(object, object, object)
    tabNameChange = QtCore.pyqtSignal(str, bool, int, int)
    tabSaved = QtCore.pyqtSignal()
    changeIcon = QtCore.pyqtSignal(str)
    geneSetsRequested = QtCore.pyqtSignal(object)
    tabLoaded = QtCore.pyqtSignal(int, str, object)
    tabReverted = QtCore.pyqtSignal(int)

    EXCLUDED_FUNCS = set()
    SUMMARY_FUNCS = set()
    CLUSTERING_FUNCS = ()
    GENERAL_FUNCS = ()
    THREADED_FUNCS = set()

    def __init__(self, parent=None, undo_stack: QtWidgets.QUndoStack = None, tab_id: int = None):
        super().__init__(parent)
        self.tab_id = tab_id
        self.undo_stack = undo_stack
        self.sup_layout = QtWidgets.QVBoxLayout(self)
        self.container = QtWidgets.QWidget(self)
        self.layout = QtWidgets.QVBoxLayout(self.container)
        self.scroll = QtWidgets.QScrollArea()
        self.scroll.setWidget(self.container)
        self.scroll.setVerticalScrollBarPolicy(QtCore.Qt.ScrollBarAlwaysOn)
        self.scroll.setHorizontalScrollBarPolicy(QtCore.Qt.ScrollBarAlwaysOff)
        self.scroll.setWidgetResizable(True)

        self.name = None
        self.creation_time = time.time()
        self.object_views = []

        # initiate function stack
        self.stack = QtWidgets.QStackedWidget(self)
        self.button_box = QtWidgets.QButtonGroup(self)
        self.stack_buttons = []
        self.stack_widgets = {}

        # initiate standard tab groups and their widget containers
        self.overview_group = QtWidgets.QGroupBox('Data overview')
        self.overview_grid = QtWidgets.QGridLayout(self.overview_group)
        self.overview_widgets = {}

        self.function_group = QtWidgets.QGroupBox('Apply functions')
        self.function_grid = QtWidgets.QGridLayout(self.function_group)
        self.function_widgets = {}
        self.layout.insertWidget(2, self.function_group)
        self.function_group.setVisible(False)

        self.stdout_group = QtWidgets.QGroupBox('Log')
        self.stdout_grid = QtWidgets.QGridLayout(self.stdout_group)
        self.stdout_widgets = {}

        # initiate apply button
        self.apply_button = QtWidgets.QPushButton('Apply')
        self.apply_button.clicked.connect(self.apply_function)
        self.layout.insertWidget(3, self.apply_button)
        self.apply_button.setVisible(False)

        # initiate the splitter layout for the tab
        self.splitter = QtWidgets.QSplitter(QtCore.Qt.Vertical)
        self.sup_layout.addWidget(self.splitter)
        self.splitter.addWidget(self.scroll)
        self.splitter.addWidget(self.stdout_group)
        self.splitter.setStretchFactor(0, 1)

        self.init_stdout_ui()

    def obj(self):
        raise NotImplementedError

    def obj_name(self):
        raise NotImplementedError

    def update_obj(self, obj):
        raise NotImplementedError

    def update_tab(self, is_unsaved: bool = True):
        raise NotImplementedError

    def obj_type(self):
        raise NotImplementedError

    def obj_properties(self) -> dict:
        return {}

    def is_empty(self):
        return True

    def init_overview_ui(self):
        raise NotImplementedError

    def init_function_ui(self):
        self.function_group.setVisible(True)
        sorted_actions = self.get_all_actions()

        self.stack.addWidget(QtWidgets.QWidget(self))  # start the stack empty
        i = 0
        for i, action_type in enumerate(sorted_actions):
            bttn = QtWidgets.QPushButton(action_type)
            bttn.setCheckable(True)
            bttn.setStyleSheet('''QPushButton::checked {background-color : purple;
                                                        color: white;
                                                        border: 1px solid #ba32ba;
                                                        border-radius: 4px;}''')
            bttn.setEnabled(len(sorted_actions[action_type]) > 0)  # disable the button if there are no relevant actions

            self.stack_widgets[action_type] = FuncTypeStack(sorted_actions[action_type], self.obj())
            self.stack_widgets[action_type].funcSelected.connect(self.apply_button.setVisible)
            self.stack_widgets[action_type].funcSelected.connect(self._check_for_special_functions)
            self.stack_widgets[action_type].geneSetsRequested.connect(self.geneSetsRequested)
            self.stack.addWidget(self.stack_widgets[action_type])
            bttn.clicked.connect(functools.partial(self.stack.setCurrentIndex, i + 1))
            self.button_box.addButton(bttn)
            self.stack_buttons.append(bttn)
            self.function_grid.addWidget(bttn, 0, i)
        self.stack.currentChanged.connect(self._update_stack_status)

        self.function_grid.addWidget(self.stack, 1, 0, 1, i + 1)

    def _check_for_special_functions(self):
        pass

    def _update_stack_status(self, ind: int):
        self.stack.widget(ind).check_selection_status()

    def get_all_actions(self):
        assert self.obj() is not None, "No object was loaded!"
        all_methods = dir(self.obj())
        public_methods = [mthd for mthd in all_methods if
                          (not mthd.startswith('_')) and (callable(getattr(type(self.obj()), mthd))) and (
                              mthd not in self.EXCLUDED_FUNCS)]
        sorted_methods = {'Filter': [], 'Normalize': [], 'Summarize': [], 'Visualize': [], 'Cluster': [], 'General': []}

        for method in public_methods:
            if method in self.SUMMARY_FUNCS:
                sorted_methods['Summarize'].append(method)
            elif method in self.CLUSTERING_FUNCS:
                sorted_methods['Cluster'].append(method)
            elif method in self.GENERAL_FUNCS:
                sorted_methods['General'].append(method)
            elif 'normalize' in method:
                sorted_methods['Normalize'].append(method)
            elif 'filter' in method or 'split' in method:
                sorted_methods['Filter'].append(method)
            else:
                sorted_methods['Visualize'].append(method)

        return sorted_methods

    def apply_function(self):
        this_stack: FuncTypeStack = self.stack.currentWidget()
        func_name = this_stack.get_function_name()
        func_params = this_stack.get_function_params()
        predecessors = this_stack.get_function_predecessors()
        if func_params.get('inplace', False):
            if func_name in self.THREADED_FUNCS:
                command = InplaceCachedCommand(self, func_name, args=[], kwargs=func_params,
                                               description=f'Apply "{func_name}"', predecessors=predecessors)
            else:
                command = InplaceCommand(self, func_name, args=[], kwargs=func_params,
                                         description=f'Apply "{func_name}"', predecessors=predecessors)
            self.undo_stack.push(command)
        else:
            self._apply_function_from_params(func_name, args=[], kwargs=func_params, predecessors=predecessors)

    def _apply_function_from_params(self, func_name, args: list, kwargs: dict, finish_slot=None, job_id: int = None,
                                    predecessors: list = None):
        partial = functools.partial(getattr(self.obj(), func_name), *args, **kwargs)
        source_name = generic.get_method_readable_name(partial.func)
        job_id = JOB_COUNTER.get_id() if job_id is None else job_id
        predecessors = predecessors if isinstance(predecessors, list) else []
        worker = gui_widgets.Worker(partial, job_id, predecessors + [self.tab_id], source_name)

        if func_name in self.THREADED_FUNCS:
            self.startedJob.emit(self, worker, [self.functionApplied.emit, finish_slot])
            return

        prev_name = self.get_tab_name()
        worker.finished.connect(check_run_success)
        worker.finished.connect(self.functionApplied.emit)
        result = worker.run()
        if kwargs.get('inplace', False):
            self.tab_id = JOB_COUNTER.get_id()
            self.itemSpawned.emit(f"'{source_name}'\noutput", self.tab_id, job_id, self.obj())

        self.update_tab(prev_name != self.obj_name())
        self.process_outputs(result, job_id, source_name)

    def process_outputs(self, outputs, job_id: int, source_name: str = ''):
        if isinstance(outputs, (list, tuple)) and len(outputs) == 0:
            return
        elif validation.isinstanceinh(outputs, (filtering.Filter, fastq.filtering.Filter)):
            new_id = JOB_COUNTER.get_id()
            self.itemSpawned.emit(f"'{source_name}'\noutput", new_id, job_id, outputs)
            self.filterObjectCreated.emit(outputs, new_id)
        elif validation.isinstanceinh(outputs, enrichment.FeatureSet):
            new_id = JOB_COUNTER.get_id()
            self.itemSpawned.emit(f"'{source_name}'\noutput", new_id, job_id, outputs)
            self.featureSetCreated.emit(outputs, new_id)
        elif isinstance(outputs, pd.DataFrame):
            new_id = JOB_COUNTER.get_id()
            self.itemSpawned.emit(f"'{source_name}'\noutput", new_id, job_id, outputs)
            self.object_views.append(gui_windows.DataFrameView(outputs, source_name, self))
            self.object_views[-1].show()
        elif isinstance(outputs, np.ndarray):
            df = pd.DataFrame(outputs)
            self.process_outputs(df, job_id, source_name)
        elif isinstance(outputs, plt.Figure):
            new_id = JOB_COUNTER.get_id()
            self.itemSpawned.emit(f"'{source_name}'\ngraph", new_id, job_id, outputs)
        elif isinstance(outputs, (tuple, list)):
            if validation.isinstanceiter_inh(outputs,
                                             (filtering.Filter, fastq.filtering.Filter, enrichment.FeatureSet)):
                dialog = MultiKeepWindow(outputs, job_id, self)
                dialog.accepted.connect(functools.partial(self._multi_keep_window_accepted, dialog, source_name))
                dialog.exec()
            else:
                for output in outputs:
                    self.process_outputs(output, job_id, source_name)
        elif isinstance(outputs, dict):
            tab_name = self.get_tab_name()
            for this_src_name, output in outputs.items():
                self.process_outputs(output, job_id, f"{this_src_name} {tab_name}")

    def _multi_keep_window_accepted(self, dialog: 'MultiKeepWindow', source_name: str):
        kept_outputs = dialog.result()
        job_id = dialog.job_id
        for i, output in enumerate(kept_outputs):
            self.process_outputs(output, job_id, source_name)

    def init_stdout_ui(self):
        self.stdout_widgets['text_edit_stdout'] = gui_widgets.StdOutTextEdit(self)
        self.stdout_widgets['text_edit_stdout'].setStyleSheet("""QTextEdit {background: #dddddd;}""")
        self.stdout_grid.addWidget(self.stdout_widgets['text_edit_stdout'], 0, 0, 3, 4)

    def get_console(self):
        return self.stdout_widgets['text_edit_stdout']

    @QtCore.pyqtSlot()
    @generic.readable_name('Rename tab')
    def rename(self, new_name: str = None):
        if new_name is None:
            new_name = self.overview_widgets['table_name'].text()
        prev_name = self.get_tab_name()
        command = RenameCommand(prev_name, new_name, self, f'Rename "{prev_name}" to "{new_name}"')
        self.undo_stack.push(command)

    def _rename(self, new_name: str = None, job_id: int = None):
        prev_id = self.tab_id
        self.tab_id = JOB_COUNTER.get_id() if job_id is not None else prev_id

        self.tabNameChange.emit(new_name, True, self.tab_id, prev_id)
        self.overview_widgets['table_name_label'].setText(f"Table name: '<b>{new_name}</b>'")
        self.overview_widgets['table_name'].setText('')
        self.name = new_name.rstrip('*')
        if job_id is not None:
            worker_output = gui_widgets.WorkerOutput(self.obj(), functools.partial(self.rename, new_name=new_name),
                                                     job_id, [prev_id])
            self.functionApplied.emit(worker_output)
            self.itemSpawned.emit(self.name, self.tab_id, job_id, self.obj())

    def cache(self):
        raise NotImplementedError

    def get_tab_name(self):
        return self.name.rstrip("*")

    # custom method to write anything printed out to console/terminal to my QTextEdit widget via append function.
    def output_terminal_written(self, text):
        self.stdout_widgets['stdout'].append(text)


class SetTabPage(TabPage):
    EXCLUDED_FUNCS = {'union', 'intersection', 'majority_vote_intersection', 'difference', 'symmetric_difference'}
    SUMMARY_FUNCS = {'biotypes_from_ref_table', 'biotypes_from_gtf'}
    GENERAL_FUNCS = {'translate_gene_ids', 'map_orthologs_phylomedb', 'map_orthologs_orthoinspector',
                     'map_orthologs_ensembl', 'map_orthologs_panther', 'find_paralogs_panther', 'find_paralogs_ensembl'}
    THREADED_FUNCS = {'translate_gene_ids', 'filter_by_kegg_annotations', 'filter_by_go_annotations',
                      'map_orthologs_phylomedb', 'map_orthologs_orthoinspector',
                      'map_orthologs_ensembl', 'map_orthologs_panther', 'find_paralogs_panther', 'find_paralogs_ensembl'
                      }

    def __init__(self, set_name: str, gene_set: typing.Union[set, enrichment.FeatureSet] = None, parent=None,
                 undo_stack: QtWidgets.QUndoStack = None, tab_id: int = None):
        super().__init__(parent, undo_stack, tab_id)
        if gene_set is None:
            gene_set = enrichment.FeatureSet(set(), set_name)
        elif isinstance(gene_set, set):
            gene_set = enrichment.FeatureSet(gene_set, set_name)
        self.gene_set = gene_set
        self.name = set_name

        self.init_overview_ui()
        self.init_function_ui()

    def obj(self):
        return self.gene_set

    def obj_name(self):
        return self.gene_set.set_name

    def update_obj(self, obj: set):
        self.update_gene_set(obj)

    def obj_type(self):
        return type(self.gene_set)

    def get_index_string(self):
        return "\n".join(self.obj())

    def start_from_gene_set(self, tab_id: int, gene_set: set):
        self.tab_id = tab_id
        self.update_obj(gene_set)

    def init_overview_ui(self):
        this_row = 0
        self.layout.insertWidget(0, self.overview_group)
        self.overview_widgets['table_name_label'] = QtWidgets.QLabel(f"Gene set name: '<b>{self.obj().set_name}</b>'")
        self.overview_widgets['table_name_label'].setWordWrap(True)

        self.overview_widgets['preview'] = gui_widgets.ReactiveListWidget()
        self.update_set_preview()
        self.overview_grid.addWidget(self.overview_widgets['table_name_label'], this_row, 0, 1, 4)
        this_row += 1
        self.overview_widgets['table_name'] = QtWidgets.QLineEdit()
        self.overview_widgets['rename_label'] = QtWidgets.QLabel('Rename your gene set (optional):')
        self.overview_widgets['rename_button'] = QtWidgets.QPushButton('Rename')
        self.overview_widgets['rename_button'].clicked.connect(self.rename)

        self.overview_grid.addWidget(self.overview_widgets['rename_label'], this_row, 0)
        self.overview_grid.addWidget(self.overview_widgets['table_name'], this_row, 1)
        self.overview_grid.addWidget(self.overview_widgets['rename_button'], this_row, 2)
        this_row += 1
        self.overview_grid.addWidget(self.overview_widgets['preview'], this_row, 0, 1, 4)
        this_row += 1

        self.overview_widgets['shape'] = QtWidgets.QLabel()
        self.overview_grid.addWidget(self.overview_widgets['shape'], this_row, 0, 1, 2)
        self.update_set_shape()

        self.overview_widgets['view_button'] = QtWidgets.QPushButton('View full gene set')
        self.overview_widgets['view_button'].clicked.connect(self.view_full_gene_set)
        self.overview_grid.addWidget(self.overview_widgets['view_button'], this_row, 2, 2, 1)

        self.overview_widgets['save_button'] = QtWidgets.QPushButton('Save gene set')
        self.overview_widgets['save_button'].clicked.connect(self.save_file)
        self.overview_grid.addWidget(self.overview_widgets['save_button'], this_row, 3, 2, 1)
        this_row += 1
        self.overview_grid.addWidget(QtWidgets.QLabel(''), this_row, 0, 1, 1)
        this_row += 1

    def view_full_gene_set(self):
        set_window = gui_windows.GeneSetView(self.gene_set.gene_set, self.get_tab_name(), self)
        self.overview_widgets['full_table_view'] = set_window
        set_window.show()

    def save_file(self):
        default_name = parsing.slugify(self.get_tab_name().rstrip("*")) + '.txt'
        filename, _ = QtWidgets.QFileDialog.getSaveFileName(self, "Save gene set",
                                                            str(Path.home().joinpath(default_name)),
                                                            "Text document (*.txt);;"
                                                            "All Files (*)")
        if filename:
            self.gene_set.save_txt(filename)
            print(f"Successfully saved at {io.get_datetime()} under {filename}")
            self.tabSaved.emit()

    @QtCore.pyqtSlot()
    def _rename(self, new_name: str = None, job_id: int = None):
        self.gene_set.change_set_name(new_name)
        super()._rename(new_name, job_id)

    def cache(self) -> str:
        base_str = str(time.time_ns()) + self.get_tab_name() + str(len(self.gene_set))
        hex_hash = hashlib.sha1(base_str.encode('utf-8')).hexdigest()
        filename = f"{hex_hash}.txt"
        io.cache_gui_file(self.gene_set.gene_set, filename)
        return filename

    def is_empty(self):
        return self.gene_set is None or len(self.gene_set) == 0

    def update_set_shape(self):
        if self.gene_set is not None:
            shape = len(self.gene_set)
            self.overview_widgets['shape'].setText(f'This gene set contains {shape} features')

    def update_set_preview(self):
        if self.gene_set is not None:
            self.overview_widgets['preview'].addItems([str(item) for item in self.gene_set])
        if self.overview_widgets['preview'].count() > 0:
            item_height = self.overview_widgets['preview'].sizeHintForRow(0)
            self.overview_widgets['preview'].setFixedHeight((item_height + 8) * 4)

    def update_gene_set(self, gene_set: set, set_name: str = None):
        if isinstance(gene_set, enrichment.FeatureSet):
            self.gene_set = gene_set
        else:
            set_name = self.gene_set.set_name if set_name is None else set_name
            self.gene_set = enrichment.FeatureSet(gene_set, set_name)
        self.update_tab()

    def update_tab(self, is_unsaved: bool = True):
        self.update_set_shape()
        self.update_set_preview()
        self.changeIcon.emit('set')

    def get_all_actions(self):
        sorted_methods = super().get_all_actions()
        for discarded_category in ('Normalize', 'Cluster', 'Visualize'):
            sorted_methods[discarded_category] = []
        return sorted_methods


class FuncTypeStack(QtWidgets.QWidget):
    EXCLUDED_PARAMS = {'self', 'backend', 'gui_mode'}
    if FROZEN_ENV:
        EXCLUDED_PARAMS.add('parallel_backend')

    NO_FUNC_CHOSEN_TEXT = "Choose a function..."
    funcSelected = QtCore.pyqtSignal(bool)
    geneSetsRequested = QtCore.pyqtSignal(object)
    __slots__ = {'parameter_widgets': 'widgets for function parameters',
                 'layout': 'layout',
                 'parameter_grid': 'layout for function parameters',
                 'func_combo': 'combo box for choosing functions',
                 'func_help_button': 'help button for function combo box',
                 'func_combo_layout': 'layout for function combo box',
                 'func': 'dict of functions',
                 'filter_obj': 'filtering.Filter object to apply functions to',
                 'excluded_params': 'parameters to be excluded from functions',
                 'pipeline_mode': 'indicating if in the function selector is in Pipeline mode'}

    def __init__(self, funcs: typing.Iterable, filter_obj: filtering.Filter, parent=None,
                 additional_excluded_params: set = None, pipeline_mode: bool = False):
        super().__init__(parent)
        self.parameter_widgets = {}
        self.layout = QtWidgets.QVBoxLayout(self)
        self.parameter_grid = QtWidgets.QGridLayout()
        self.func_combo = QtWidgets.QComboBox(self)
        self.func_help_button = gui_widgets.HelpButton(self)
        self.func_combo_layout = QtWidgets.QHBoxLayout()
        self.funcs = {}
        for func in funcs:
            self.funcs[generic.get_method_readable_name(func, filter_obj)] = func
        self.filter_obj = filter_obj
        self.excluded_params = self.EXCLUDED_PARAMS.copy()
        if additional_excluded_params is not None:
            self.excluded_params.update(additional_excluded_params)
        self.pipeline_mode = pipeline_mode
        self.init_ui()

    def init_ui(self):

        self.layout.addLayout(self.func_combo_layout)
        self.layout.addLayout(self.parameter_grid)
        self.func_combo_layout.addWidget(self.func_combo)
        self.func_combo_layout.addWidget(self.func_help_button)
        self._set_empty_tooltip()
        self.layout.addStretch(1)
        self.func_combo.addItem(self.NO_FUNC_CHOSEN_TEXT)
        self.func_combo.addItems(sorted(self.funcs.keys()))
        self.func_combo.currentTextChanged.connect(self.update_parameter_ui)

    def _set_empty_tooltip(self):
        txt = "Choose a function from this list to read its description. "
        self.func_combo.setToolTip(txt)
        self.func_help_button.set_desc_help(txt)

    def deselect(self):
        self.func_combo.setCurrentIndex(0)

    def check_selection_status(self):
        self.funcSelected.emit(self.get_function_name() != self.NO_FUNC_CHOSEN_TEXT)

    def update_parameter_ui(self):
        # delete previous widgets
        gui_widgets.clear_layout(self.parameter_grid)
        self.parameter_widgets = {}

        chosen_func_name = self.get_function_name()
        if chosen_func_name == self.NO_FUNC_CHOSEN_TEXT:
            self._set_empty_tooltip()
            self.funcSelected.emit(False)
            return

        signature = generic.get_method_signature(chosen_func_name, self.filter_obj)
        desc, param_desc = io.get_method_docstring(chosen_func_name, self.filter_obj)
        self.func_combo.setToolTip(desc)
        self.func_help_button.set_param_help(self.get_function_readable_name(), desc)

        i = 1
        for name, param in signature.items():
            if name in self.excluded_params:
                continue
            self.parameter_widgets[name] = gui_widgets.param_to_widget(param, name, pipeline_mode=self.pipeline_mode)
            self.connect_widget(self.parameter_widgets[name])

            label = QtWidgets.QLabel(f'{name}:', self.parameter_widgets[name])
            if name in param_desc:
                label.setToolTip(param_desc[name])
                help_button = gui_widgets.HelpButton()
                self.parameter_grid.addWidget(help_button, i, 2)
                help_button.set_param_help(name, param_desc[name])

            self.parameter_grid.addWidget(label, i, 0)
            self.parameter_grid.addWidget(self.parameter_widgets[name], i, 1)

            i += 1

        help_address = f"https://guyteichman.github.io/RNAlysis/build/rnalysis.filtering." \
                       f"{type(self.filter_obj).__name__}.{chosen_func_name}.html"
        self.parameter_widgets['help_link'] = QtWidgets.QLabel(
            text=f'<a href="{help_address}">Open documentation for function '
                 f'<b>{type(self.filter_obj).__name__}.{chosen_func_name}</b></a>')
        self.parameter_widgets['help_link'].setOpenExternalLinks(True)
        self.parameter_grid.addWidget(self.parameter_widgets['help_link'], i + 1, 0, 1, 2)
        self.parameter_grid.setColumnStretch(1, 1)
        self.funcSelected.emit(True)

    def connect_widget(self, widget: QtWidgets.QWidget):
        if self.pipeline_mode:
            return

        if isinstance(widget, (gui_widgets.ComboBoxOrOtherWidget, gui_widgets.OptionalWidget)):
            self.connect_widget(widget.other)
        elif isinstance(widget, (gui_widgets.TableColumnPicker, gui_widgets.TableColumnPicker)):
            widget.add_columns(self.filter_obj.columns)
        elif isinstance(widget, gui_widgets.GeneSetComboBox):
            widget.boxOpened.connect(functools.partial(self.geneSetsRequested.emit, widget))

    def get_function_params(self):
        func_params = {}
        for param_name, widget in self.parameter_widgets.items():
            if param_name in {'function_combo', 'help_link'}:
                continue
            val = gui_widgets.get_val_from_widget(widget)

            func_params[param_name] = val
        return func_params

    def get_function_predecessors(self):
        predecessors = []
        for param_name, widget in self.parameter_widgets.items():
            if isinstance(widget, (gui_widgets.OptionalWidget, gui_widgets.ComboBoxOrOtherWidget)):
                if widget.other.isEnabled():
                    widget = widget.other

            if isinstance(widget, gui_widgets.GeneSetComboBox):
                predecessors.append(widget.current_id())
        return predecessors

    def get_function_readable_name(self):
        return self.func_combo.currentText()

    def get_function_name(self):
        readable_name = self.get_function_readable_name()
        if readable_name == self.NO_FUNC_CHOSEN_TEXT:
            return self.NO_FUNC_CHOSEN_TEXT
        name = self.funcs[readable_name]
        return name


class FilterTabPage(TabPage):
    EXCLUDED_FUNCS = {'union', 'intersection', 'majority_vote_intersection', 'difference', 'symmetric_difference',
                      'from_folder', 'from_folder_htseqcount', 'save_txt', 'save_csv', 'save_table', 'save_parquet',
                      'from_dataframe', 'print_features'}
    CLUSTERING_FUNCS = {'split_kmeans': 'K-Means', 'split_kmedoids': 'K-Medoids',
                        'split_hierarchical': 'Hierarchical (Agglomerative)', 'split_hdbscan': 'HDBSCAN',
                        'split_clicom': 'CLICOM (Ensemble)'}
    SUMMARY_FUNCS = {'describe', 'head', 'tail', 'biotypes_from_ref_table', 'biotypes_from_gtf', 'print_features'}
    GENERAL_FUNCS = {'sort', 'sort_by_principal_component', 'transform', 'translate_gene_ids',
                     'differential_expression_deseq2', 'fold_change', 'average_replicate_samples', 'drop_columns',
                     'differential_expression_limma_voom', 'map_orthologs_phylomedb', 'map_orthologs_orthoinspector',
                     'map_orthologs_ensembl', 'map_orthologs_panther', 'find_paralogs_panther'}
    THREADED_FUNCS = {'translate_gene_ids', 'differential_expression_deseq2', 'filter_by_kegg_annotations',
                      'filter_by_go_annotations', 'differential_expression_limma_voom', 'map_orthologs_phylomedb',
                      'map_orthologs_orthoinspector',
                      'map_orthologs_ensembl', 'map_orthologs_panther', 'find_paralogs_panther',
                      'find_paralogs_ensembl'}
    startedClustering = QtCore.pyqtSignal(object, object, object)
    widthChanged = QtCore.pyqtSignal()

    def __init__(self, parent=None, undo_stack: QtWidgets.QUndoStack = None, tab_id: int = None):
        super().__init__(parent, undo_stack, tab_id)
        self.filter_obj = None

        self.basic_group = QtWidgets.QGroupBox('Load a data table')
        self.basic_grid = QtWidgets.QGridLayout(self.basic_group)
        self.basic_widgets = {}
        self.basic_param_container = QtWidgets.QWidget(self)
        self.basic_param_widgets = {}
        self.basic_param_grid = QtWidgets.QGridLayout(self.basic_param_container)

        self.clicom_window = None
        self.deseq_window = None
        self.limma_window = None

        self.widthChanged.connect(self.update_table_preview_width)

        self.init_basic_ui()

    def obj(self):
        return self.filter_obj

    def obj_name(self):
        name = self.filter_obj.fname.name
        for suffix in ['.csv', '.tsv', '.txt', '.parquet']:
            if self.filter_obj.fname.suffix == suffix:
                name = self.filter_obj.fname.stem
                break
        return name

    def update_obj(self, obj: filtering.Filter):
        self.filter_obj = obj
        self.update_tab()

    def obj_type(self):
        return type(self.filter_obj)

    def obj_properties(self) -> dict:
        if self.obj_type() == filtering.CountFilter:
            return dict(is_normalized=self.filter_obj.is_normalized)
        elif self.obj_type() == filtering.DESeqFilter:
            return dict(log2fc_col=self.filter_obj.log2fc_col, padj_col=self.filter_obj.padj_col)
        elif self.obj_type() == filtering.FoldChangeFilter:
            return dict(numerator_name=self.filter_obj.numerator, denominator_name=self.filter_obj.denominator)
        else:
            return {}

    @QtCore.pyqtSlot()
    def _rename(self, new_name: str = None, job_id: int = None):
        self.filter_obj.fname = Path(
            os.path.join(str(self.filter_obj.fname.parent), f"{new_name}{self.filter_obj.fname.suffix}"))
        super()._rename(new_name, job_id)

    def cache(self):
        base_str = str(time.time_ns()) + str(self.filter_obj.fname) + str(len(self.filter_obj.shape))
        hex_hash = hashlib.sha1(base_str.encode('utf-8')).hexdigest()
        filename = f"{hex_hash}.parquet"
        io.cache_gui_file(self.filter_obj.df, filename)
        return filename

    def is_empty(self):
        return self.filter_obj is None

    def get_table_type(self):
        return FILTER_OBJ_TYPES_INV[type(self.filter_obj).__name__]

    def update_table_name_label(self):
        self.overview_widgets['table_name_label'].setText(f"Table name: '<b>{self.get_tab_name()}</b>'")
        self.overview_widgets['table_name_label'].setWordWrap(True)

    def init_overview_ui(self):
        this_row = 0
        self.layout.insertWidget(1, self.overview_group)
        self.overview_widgets['table_type_label'] = QtWidgets.QLabel(
            f"Table type: {self.get_table_type()}")
        self.overview_widgets['table_name_label'] = QtWidgets.QLabel()
        self.overview_widgets['table_name_label'].setWordWrap(True)

        self.overview_widgets['preview'] = gui_widgets.ReactiveTableView()
        self.overview_widgets['preview'].setVerticalScrollBarPolicy(QtCore.Qt.ScrollBarAlwaysOff)
        self.overview_widgets['preview'].setHorizontalScrollBarPolicy(QtCore.Qt.ScrollBarAlwaysOff)

        self.overview_grid.addWidget(self.overview_widgets['table_name_label'], this_row, 0, 1, 4)
        this_row += 1
        self.overview_widgets['table_name'] = QtWidgets.QLineEdit()
        self.overview_widgets['rename_label'] = QtWidgets.QLabel('Rename your table (optional):')
        self.overview_widgets['rename_button'] = QtWidgets.QPushButton('Rename')
        self.overview_widgets['rename_button'].clicked.connect(self.rename)

        self.overview_grid.addWidget(self.overview_widgets['rename_label'], this_row, 0)
        self.overview_grid.addWidget(self.overview_widgets['table_name'], this_row, 1)
        self.overview_grid.addWidget(self.overview_widgets['rename_button'], this_row, 2)
        this_row += 1
        self.overview_grid.addWidget(self.overview_widgets['preview'], this_row, 0, 1, 4)
        this_row += 1

        self.overview_widgets['save_button'] = QtWidgets.QPushButton('Save table')
        self.overview_widgets['save_button'].clicked.connect(self.save_file)
        self.overview_grid.addWidget(self.overview_widgets['save_button'], this_row, 3, 2, 1)

        self.overview_widgets['shape'] = QtWidgets.QLabel()
        self.overview_grid.addWidget(self.overview_widgets['shape'], this_row, 0, 1, 2)

        self.overview_widgets['view_button'] = QtWidgets.QPushButton('View full table')
        self.overview_widgets['view_button'].clicked.connect(self.view_full_dataframe)
        self.overview_grid.addWidget(self.overview_widgets['view_button'], this_row, 2, 2, 1)

        this_row += 1
        self.overview_grid.addWidget(self.overview_widgets['table_type_label'], this_row, 0, 1, 1)
        this_row += 1

        self.update_tab(False)

    def update_filter_obj_shape(self):
        if self.is_empty():
            return
        shape = self.filter_obj.shape
        if len(shape) == 1:
            shape = (shape[0], 1)
        self.overview_widgets['shape'].setText(f'This table contains {shape[0]} rows, {shape[1]} columns')

    def _change_start_button_state(self, is_legal: bool):
        self.basic_widgets['start_button'].setEnabled(is_legal)

    def update_basic_ui(self):
        # clear previous layout
        gui_widgets.clear_layout(self.basic_param_grid)
        self.basic_param_widgets = {}
        func_name = '__init__'
        filter_obj_type = FILTER_OBJ_TYPES[self.basic_widgets['table_type_combo'].currentText()]
        signature = generic.get_method_signature(func_name, filter_obj_type)
        desc, param_desc = io.get_method_docstring(func_name, filter_obj_type)
        self.basic_widgets['table_type_combo'].setToolTip(desc)
        i = 1
        for name, param in signature.items():
            if name in INIT_EXCLUDED_PARAMS:
                continue
            self.basic_param_widgets[name] = gui_widgets.param_to_widget(param, name)
            label = QtWidgets.QLabel(f'{name}:', self.basic_param_widgets[name])
            if name in param_desc:
                label.setToolTip(param_desc[name])
                help_button = gui_widgets.HelpButton()
                self.basic_param_grid.addWidget(help_button, i, 2)
                help_button.set_param_help(name, param_desc[name])

            self.basic_param_grid.addWidget(label, i, 0)
            self.basic_param_grid.addWidget(self.basic_param_widgets[name], i, 1)
            i += 1

    def init_basic_ui(self):
        self.layout.insertWidget(0, self.basic_group)
        self.basic_widgets['table_type_combo'] = QtWidgets.QComboBox()
        self.basic_widgets['table_type_combo'].addItems(FILTER_OBJ_TYPES.keys())
        self.basic_widgets['table_type_combo'].currentIndexChanged.connect(self.update_basic_ui)
        self.basic_widgets['table_type_combo'].setCurrentText('Other table')

        self.basic_widgets['start_button'] = QtWidgets.QPushButton('Start')
        self.basic_widgets['start_button'].clicked.connect(self.start)
        self.basic_widgets['start_button'].setEnabled(False)

        self.basic_widgets['file_path'] = gui_widgets.PathLineEdit(file_types="Data table "
                                                                              "(*.csv;*.tsv;*.txt;*.parquet);;"
                                                                              "All Files (*)")
        self.basic_widgets['file_path'].textChanged.connect(self._change_start_button_state)

        self.basic_widgets['table_name'] = QtWidgets.QLineEdit()

        self.basic_widgets['file_label'] = QtWidgets.QLabel('Choose a file:')
        self.basic_widgets['type_label'] = QtWidgets.QLabel('Choose table type:')
        self.basic_widgets['name_label'] = QtWidgets.QLabel('Name your table (optional):')

        self.basic_grid.addWidget(self.basic_widgets['file_label'], 0, 0, 1, 2)
        self.basic_grid.addWidget(self.basic_widgets['type_label'], 0, 2)
        self.basic_grid.addWidget(self.basic_widgets['name_label'], 0, 3)

        self.basic_grid.addWidget(self.basic_widgets['file_path'], 1, 0, 1, 2)
        self.basic_grid.addWidget(self.basic_widgets['table_type_combo'], 1, 2)
        self.basic_grid.addWidget(self.basic_widgets['table_name'], 1, 3)
        self.basic_grid.addWidget(self.basic_param_container, 2, 0, 1, 4)
        self.basic_grid.addWidget(self.basic_widgets['start_button'], 3, 0, 1, 4)

        self.basic_grid.addWidget(QtWidgets.QWidget(self), 4, 0, 1, 4)
        self.basic_grid.addWidget(QtWidgets.QWidget(self), 0, 4, 4, 1)
        self.basic_grid.setRowStretch(4, 1)
        self.basic_grid.setColumnStretch(4, 1)

    def _check_for_special_functions(self, is_selected: bool):
        if not is_selected:
            return
        this_stack: FuncTypeStack = self.stack.currentWidget()
        func_name = this_stack.get_function_name()

        if func_name == 'split_clicom':
            other_clustering_funcs = self.CLUSTERING_FUNCS.copy()
            other_clustering_funcs.pop('split_clicom')
            this_stack.deselect()
            self.clicom_window = ClicomWindow(other_clustering_funcs, self.filter_obj, self)
            self.clicom_window.paramsAccepted.connect(
                functools.partial(self._apply_function_from_params, func_name))
            self.clicom_window.show()
        elif func_name == 'differential_expression_deseq2':
            this_stack.deselect()
            self.deseq_window = DESeqWindow(self)
            self.deseq_window.paramsAccepted.connect(functools.partial(self._apply_function_from_params, func_name))
            self.deseq_window.show()
        elif func_name == 'differential_expression_limma_voom':
            this_stack.deselect()
            self.limma_window = LimmaWindow(self)
            self.limma_window.paramsAccepted.connect(functools.partial(self._apply_function_from_params, func_name))
            self.limma_window.show()

    def view_full_dataframe(self):
        df_window = gui_windows.DataFrameView(self.filter_obj.df, self.name, self)
        self.overview_widgets['full_table_view'] = df_window
        df_window.show()

    def update_table_preview(self):
        if self.is_empty():
            return
        model = gui_windows.DataFramePreviewModel(self.filter_obj.df)
        self.overview_widgets['preview'].setModel(model)
        self.update_table_preview_width()

    def update_table_preview_width(self):
        if self.is_empty():
            return
        self.overview_widgets['preview'].resizeColumnsToContents()
        self.overview_widgets['preview'].resizeRowsToContents()

        new_width = min(self.overview_widgets['preview'].horizontalHeader().length() + self.overview_widgets[
            'preview'].verticalHeader().width(), self.width() - 100)
        self.overview_widgets['preview'].setFixedWidth(new_width)

        self.overview_widgets['preview'].setFixedHeight(
            self.overview_widgets['preview'].verticalHeader().length() + self.overview_widgets[
                'preview'].horizontalHeader().height())

    def resizeEvent(self, a0: QtGui.QResizeEvent) -> None:
        if a0.size().width() != a0.oldSize().width():
            self.widthChanged.emit()

    def update_tab(self, is_unsaved: bool = True):
        if self.is_empty():
            return
        self.update_filter_obj_shape()
        self.update_table_preview()
        self.name = str(self.obj_name())
        self.tabNameChange.emit(self.name, is_unsaved, -1, -1)
        self.update_table_name_label()

    def _apply_function_from_params(self, func_name, args: list, kwargs: dict, finish_slot=None, job_id: int = None,
                                    predecessors: list = None):
        # since clustering functions can be computationally intensive, start them in another thread.
        # furthermode, make sure they use the 'multiprocessing' backend instead of the 'loky' backend -
        # otherwise this could cause issues in Pyinstaller-frozen versions of RNAlysis
        if func_name in self.CLUSTERING_FUNCS:
            kwargs['gui_mode'] = True
            if FROZEN_ENV:
                kwargs['parallel_backend'] = 'multiprocessing'
            partial = functools.partial(getattr(self.filter_obj, func_name), *args, **kwargs)
            predecessors = predecessors if isinstance(predecessors, list) else []
            worker = gui_widgets.Worker(partial, JOB_COUNTER.get_id(), predecessors + [self.tab_id], func_name)
            worker.finished.connect(self.functionApplied.emit)
            self.startedClustering.emit(self, worker, finish_slot)
            return

        return super()._apply_function_from_params(func_name, args, kwargs, finish_slot, job_id, predecessors)

    def apply_pipeline(self, pipeline: filtering.Pipeline, pipeline_name: str, pipeline_id: int, inplace: bool):
        predecessors = [pipeline_id]
        if inplace:
            command = PipelineInplaceCommand(self, pipeline, pipeline_name, f'Apply Pipeline "{pipeline_name}"',
                                             predecessors)
            self.undo_stack.push(command)
        else:
            self._apply_pipeline(pipeline, pipeline_name, False, JOB_COUNTER.get_id(), predecessors)

    def _apply_pipeline(self, pipeline, pipeline_name: str, inplace: bool, job_id: int = None,
                        predecessors: list = None):

        job_name = f'Pipeline {pipeline_name} on {self.get_tab_name()}'
        partial = functools.partial(pipeline.apply_to, self.filter_obj, inplace)
        job_id = JOB_COUNTER.get_id() if job_id is None else job_id
        predecessors = predecessors if isinstance(predecessors, list) else []
        worker = gui_widgets.Worker(partial, job_id, predecessors + [self.tab_id], f"Pipeline '{pipeline_name}'")
        worker.finished.connect(check_run_success)
        worker.finished.connect(self.functionApplied.emit)
        prev_name = self.get_tab_name()
        result = worker.run()
        self.update_tab(prev_name != self.filter_obj.fname.name)
        self.process_outputs(result, job_id, job_name)

    def get_index_string(self):
        if self.is_empty():
            return ''
        return self.filter_obj.index_string

    def save_file(self):
        if self.filter_obj is None:
            warnings.warn("Cannot save an empty tab!")
            return
        default_name = parsing.slugify(str(self.filter_obj.fname.stem).rstrip("*")) + '.csv'
        filename, _ = QtWidgets.QFileDialog.getSaveFileName(self, "Save filtering result",
                                                            str(Path.home().joinpath(default_name)),
                                                            "Comma-Separated Values (*.csv);;"
                                                            "Tab-Separated Values (*.tsv);;"
                                                            "Parquet file (*.parquet);;"
                                                            "All Files (*)")
        if filename:
            suffix = Path(filename).suffix
            self.filter_obj.save_table(suffix, filename)
            print(f"Successfully saved at {io.get_datetime()} under {filename}")
            self.tabSaved.emit()

    def start(self):
        filter_obj_type = FILTER_OBJ_TYPES[self.basic_widgets['table_type_combo'].currentText()]
        file_path = self.basic_widgets['file_path'].text()
        kwargs = {}
        for name, widget in self.basic_param_widgets.items():
            kwargs[name] = gui_widgets.get_val_from_widget(widget)
        self.filter_obj = filter_obj_type(file_path, **kwargs)

        print(self.filter_obj)

        table_name_user_input = self.basic_widgets['table_name'].text()
        if table_name_user_input != '':
            new_name = table_name_user_input
            self.filter_obj._update(fname=Path(new_name).with_suffix('.csv'))
        else:
            new_name = self.filter_obj.fname.stem
        self.tabNameChange.emit(new_name, False, -1, -1)

        self.init_overview_ui()
        self.init_function_ui()

        gui_widgets.clear_layout(self.basic_grid)
        self.layout.removeWidget(self.basic_group)
        self.basic_group.deleteLater()

        self.changeIcon.emit(type(self.filter_obj).__name__)

        self.tabLoaded.emit(self.tab_id, self.name, self.obj())

    def start_from_filter_obj(self, filter_obj: filtering.Filter, tab_id: int, name: str = None):
        self.tab_id = tab_id
        self.filter_obj = filter_obj
        self.basic_group.setVisible(False)
        self.init_overview_ui()
        self.init_function_ui()
        if name is not None:
            self._rename(name)
        print(self.filter_obj)
        self.changeIcon.emit(type(self.filter_obj).__name__)


class CreatePipelineWindow(gui_widgets.MinMaxDialog, FilterTabPage):
    pipelineSaved = QtCore.pyqtSignal(str, generic.GenericPipeline)
    pipelineExported = QtCore.pyqtSignal(str, generic.GenericPipeline)
    widthChanged = QtCore.pyqtSignal()
    geneSetsRequested = QtCore.pyqtSignal()
    PIPELINE_TYPES = {name: filtering.Pipeline for name in FILTER_OBJ_TYPES.keys()}
    PIPELINE_TYPES.update({'Sequence files (paired-end)': fastq.PairedEndPipeline,
                           'Sequence files (single-end)': fastq.SingleEndPipeline})
    INV_PIPELINE_TYPES = {val.__name__: key for key, val in PIPELINE_TYPES.items()}
    __slots__ = {'pipeline': 'Pipeline object',
                 'is_unsaved': 'indicates whether the Pipeline was saved since changes were last made'}

    def __init__(self, parent=None):
        super().__init__(parent=parent)
        self.setLayout(self.layout)
        self.setWindowTitle('Create new Pipeline')
        self.setGeometry(500, 200, 900, 800)
        self.pipeline = None
        self.is_unsaved = False

    @classmethod
    def start_from_pipeline(cls, pipeline: generic.GenericPipeline, pipeline_name: str, parent=None):
        window = cls(parent)
        window.basic_widgets['pipeline_name'].setText(pipeline_name)
        window.pipeline = pipeline
        pipeline_type = cls.INV_PIPELINE_TYPES[type(pipeline).__name__]
        if isinstance(pipeline, filtering.Pipeline):
            window.filter_obj = pipeline.filter_type.__new__(pipeline.filter_type)
            pipeline_type = FILTER_OBJ_TYPES_INV[pipeline.filter_type.__name__]

        window.basic_widgets['table_type_combo'].setCurrentText(pipeline_type)
        window.init_overview_ui()
        window.init_function_ui()
        window.basic_group.setVisible(False)
        window.is_unsaved = False
        return window

    def init_basic_ui(self):
        self.layout.insertWidget(0, self.basic_group)
        self.basic_group.setTitle("Choose data table type for Pipeline")

        self.basic_widgets['table_type_combo'] = QtWidgets.QComboBox()
        self.basic_widgets['table_type_combo'].addItems(self.PIPELINE_TYPES.keys())
        self.basic_widgets['table_type_combo'].setCurrentText('Other table')

        self.basic_widgets['pipeline_name'] = QtWidgets.QLineEdit()
        self.basic_widgets['pipeline_name'].setText('New Pipeline')

        self.basic_widgets['name_label'] = QtWidgets.QLabel('Name your Pipeline:')

        self.basic_widgets['start_button'] = QtWidgets.QPushButton('Start')
        self.basic_widgets['start_button'].clicked.connect(self.start)

        self.apply_button = QtWidgets.QPushButton('Add to Pipeline')
        self.apply_button.clicked.connect(self.apply_function)
        self.layout.insertWidget(2, self.apply_button)
        self.apply_button.setVisible(False)

        self.basic_widgets['type_label'] = QtWidgets.QLabel('Choose Pipeline type:')

        self.basic_grid.addWidget(self.basic_widgets['pipeline_name'], 1, 1)
        self.basic_grid.addWidget(self.basic_widgets['name_label'], 0, 0, 1, 2)
        self.basic_grid.addWidget(self.basic_widgets['type_label'], 0, 2)
        self.basic_grid.addWidget(self.basic_widgets['table_type_combo'], 1, 2)
        self.basic_grid.addWidget(self.basic_widgets['start_button'], 1, 3)

        self.basic_grid.addWidget(QtWidgets.QWidget(), 2, 0)
        self.basic_grid.addWidget(QtWidgets.QWidget(), 0, 4)
        self.basic_grid.setRowStretch(2, 1)
        self.basic_grid.setColumnStretch(4, 1)

    def get_all_fastq_actions(self):
        func_type = 'single' if isinstance(self.pipeline, fastq.SingleEndPipeline) else 'paired'
        all_funcs = dir(fastq)
        public_funcs = []
        for func in all_funcs:
            func_obj = getattr(fastq, func)
            if func.startswith('_'):
                continue
            if not hasattr(func_obj, 'func_type'):
                continue
            if getattr(func_obj, 'func_type') not in [func_type, 'both']:
                continue
            public_funcs.append(func)

        return public_funcs

    def init_function_ui(self):
        if isinstance(self.pipeline, (fastq.SingleEndPipeline, fastq.PairedEndPipeline)):
            self.function_group.setVisible(True)
            self.funcs = {}
            actions = self.get_all_fastq_actions()
            for func in actions:
                self.funcs[generic.get_method_readable_name(func, fastq)] = func
            self.stack_widgets['main'] = FuncTypeStack(actions, fastq, self,
                                                       {'input_folder', 'fastq_folder', 'output_folder', 'r1_files',
                                                        'r2_files', 'return_new_filenames'}, True)
            self.stack_widgets['main'].funcSelected.connect(self.apply_button.setVisible)

            self.stack.addWidget(self.stack_widgets['main'])
            self.function_grid.addWidget(self.stack, 1, 0, 1, 1)
        else:
            super().init_function_ui()
            sorted_actions = self.get_all_actions()
            for action_type in sorted_actions:
                self.stack_widgets[action_type].pipeline_mode = True
                self.stack_widgets[action_type].excluded_params.add('inplace')

        self.function_group.setTitle("Add functions to Pipeline")

    def apply_function(self):
        this_stack: FuncTypeStack = self.stack.currentWidget()
        func_name = this_stack.get_function_name()
        func_params = this_stack.get_function_params()
        self.pipeline.add_function(func_name, **func_params)
        self.update_pipeline_preview()
        self.is_unsaved = True

    def _apply_function_from_params(self, func_name, args: list, kwargs: dict, finish_slot=None, job_id: int = None):
        raise NotImplementedError

    def update_table_preview_width(self):
        return

    def update_pipeline_preview(self):
        self.overview_widgets['preview'].setPlainText(str(self.pipeline))

    def update_table_preview(self):
        raise NotImplementedError

    def update_filter_obj_shape(self):
        raise NotImplementedError

    def set_file_path_bg_color(self):
        raise NotImplementedError

    def start(self):
        table_type = self.basic_widgets['table_type_combo'].currentText()
        pipeline_type = self.PIPELINE_TYPES[table_type]
        if pipeline_type == filtering.Pipeline:
            filt_obj_type = FILTER_OBJ_TYPES[self.basic_widgets['table_type_combo'].currentText()]
            self.filter_obj = filt_obj_type.__new__(filt_obj_type)
            self.pipeline = filtering.Pipeline(filt_obj_type)
        else:
            self.pipeline = pipeline_type()

        self.init_overview_ui()
        self.init_function_ui()
        self.basic_group.setVisible(False)
        self.is_unsaved = True

    def init_overview_ui(self):
        self.overview_group.setTitle("Pipeline preview")
        self.layout.insertWidget(1, self.overview_group)
        self.overview_widgets['preview'] = QtWidgets.QPlainTextEdit()
        self.overview_widgets['preview'].setReadOnly(True)
        self.update_pipeline_preview()

        self.overview_grid.addWidget(QtWidgets.QLabel(f"Pipeline name: "
                                                      f"'<b>{self._get_pipeline_name()}</b>'"), 0, 0, 1, 6)
        self.overview_grid.addWidget(self.overview_widgets['preview'], 2, 0, 2, 6)

        self.overview_grid.addWidget(QtWidgets.QLabel(f"Pipeline type: "
                                                      f"{self.basic_widgets['table_type_combo'].currentText()}"),
                                     5, 0, 1, 3)

        self.overview_widgets['remove_button'] = QtWidgets.QPushButton('Remove last function')
        self.overview_widgets['remove_button'].clicked.connect(self.remove_last_function)
        self.overview_grid.addWidget(self.overview_widgets['remove_button'], 5, 3)

        self.overview_widgets['export_button'] = QtWidgets.QPushButton('Export Pipeline')
        self.overview_widgets['export_button'].clicked.connect(self.export_pipeline)
        self.overview_grid.addWidget(self.overview_widgets['export_button'], 5, 4)

        self.overview_widgets['save_button'] = QtWidgets.QPushButton('Save Pipeline')
        self.overview_widgets['save_button'].clicked.connect(self.save_file)
        self.overview_grid.addWidget(self.overview_widgets['save_button'], 5, 5)

    def remove_last_function(self):
        try:
            self.pipeline.remove_last_function()
            self.update_pipeline_preview()
        except AssertionError:
            err = QtWidgets.QMessageBox(self)
            err.setWindowTitle('Pipeline is already empty!')
            err.setText('Cannot remove functions from the Pipeline - it is already empty!')
            err.setIcon(err.Warning)
            err.exec()

    def _get_pipeline_name(self):
        return self.basic_widgets['pipeline_name'].text()

    def export_pipeline(self):
        self.pipelineExported.emit(self._get_pipeline_name(), self.pipeline)

    def save_file(self):
        try:
            self.pipelineSaved.emit(self._get_pipeline_name(), self.pipeline)
            self.is_unsaved = False
        except Exception as e:
            print("Failed to save Pipeline")
            raise e

    def closeEvent(self, event):  # pragma: no cover
        if self.is_unsaved:
            quit_msg = "Are you sure you want to close this window without saving your Pipeline?\n" \
                       "All unsaved progress will be lost"

            reply = QtWidgets.QMessageBox.question(self, "Close 'Create Pipeline' window?",
                                                   quit_msg, QtWidgets.QMessageBox.No, QtWidgets.QMessageBox.Yes)

            if reply == QtWidgets.QMessageBox.Yes:
                event.accept()
            else:
                event.ignore()
        else:
            event.accept()


class MultiKeepWindow(gui_widgets.MinMaxDialog):
    __slots__ = {'objs': 'objects to keep or discard',
                 'job_id': 'job ID',
                 'files': 'filenames of the objects to keep',
                 'button_box': 'button box for selecting objects to keep',
                 'labels': 'labels for the objects',
                 'keep_marks': 'check boxes for the objects',
                 'names': 'potential new names for the objects',
                 'select_all': 'select all checkbox',
                 'scroll': 'scroll area widget',
                 'scroll_layout': 'layout for scroll area',
                 'scroll_widget': 'widget containing scroll area',
                 'main_layout': 'main layout for the window'}

    def __init__(self, objs: List[Union[filtering.Filter, enrichment.FeatureSet]], job_id: int, parent=None):
        super().__init__(parent)
        self.job_id = job_id
        # make sure there are no two objects with the same name
        self.objs = {}
        for obj in objs:
            key = str(obj.fname.stem) if validation.isinstanceinh(obj, (
                filtering.Filter, fastq.filtering.Filter)) else obj.set_name
            if key in self.objs:
                i = 1
                key += '_{i}'
                while key.format(i=i) in self.objs:
                    i += 1
                key = key.format(i=i)
            self.objs[key] = obj

        self.files = list(self.objs.keys())
        self.button_box = QtWidgets.QDialogButtonBox(QtWidgets.QDialogButtonBox.Ok | QtWidgets.QDialogButtonBox.Cancel)
        self.labels = dict()
        self.keep_marks = dict()
        self.names = dict()
        self.select_all = QtWidgets.QCheckBox("Select all")

        self.scroll = QtWidgets.QScrollArea()
        self.scroll_layout = QtWidgets.QGridLayout(self.scroll)
        self.scroll_widget = QtWidgets.QWidget(self.scroll)
        self.main_layout = QtWidgets.QVBoxLayout(self)

        self.init_ui()

    def init_ui(self):
        self.setLayout(self.main_layout)
        self.main_layout.addWidget(self.scroll)
        self.scroll.setWidgetResizable(True)
        self.scroll.setWidget(self.scroll_widget)
        self.scroll_widget.setLayout(self.scroll_layout)
        self.scroll_layout.setSizeConstraint(QtWidgets.QLayout.SetMinAndMaxSize)

        self.select_all.clicked.connect(self.change_all)
        self.button_box.accepted.connect(self.accept)
        self.button_box.rejected.connect(self.reject)
        self.setWindowTitle('Choose tables to keep')
        self.scroll_layout.addWidget(
            QtWidgets.QLabel(
                'Please choose which tables to keep (mandatory) out of the tables generated by the last function call, '
                'and rename them (optional):\n\n'), 0, 0, 1, 3)
        self.scroll_layout.addWidget(QtWidgets.QLabel('<b>Tables:</b>'), 1, 0)
        self.scroll_layout.addWidget(QtWidgets.QLabel('<b>Choose tables to keep</b>'), 1, 1)
        self.scroll_layout.addWidget(QtWidgets.QLabel('<b>Rename tables (optional):</b>'), 1, 2)
        self.scroll_layout.addWidget(self.select_all, 2, 1)
        for i, file in enumerate(self.files, start=3):
            self.labels[file] = QtWidgets.QLabel(file)
            self.keep_marks[file] = QtWidgets.QCheckBox("Keep table?")
            self.names[file] = QtWidgets.QLineEdit()

            self.scroll_layout.addWidget(self.labels[file], i, 0)
            self.scroll_layout.addWidget(self.keep_marks[file], i, 1)
            self.scroll_layout.addWidget(self.names[file], i, 2)
        self.main_layout.addWidget(self.button_box)

        self.scroll.setMinimumWidth(self.scroll_widget.sizeHint().width() + 150)

    def change_all(self):
        for widget in self.keep_marks.values():
            widget.setChecked(self.select_all.isChecked())

    def result(self) -> list:
        keep_tables = {file: widget.isChecked() for file, widget in self.keep_marks.items()}
        new_names = {file: gui_widgets.get_val_from_widget(widget) for file, widget in self.names.items()}
        out = []
        for file in keep_tables:
            if keep_tables[file]:
                obj = self.objs[file]
                if new_names[file] != '':
                    if validation.isinstanceinh(obj, (filtering.Filter, fastq.filtering.Filter)):
                        obj.fname = Path(f"{new_names[file]}.csv")
                    else:
                        obj.set_name = new_names[file]
                out.append(obj)
        return out


class MultiOpenWindow(QtWidgets.QDialog):
    __slots__ = {'files': 'filenames of the objects to open',
                 'button_box': 'button box for selecting objects to open',
                 'paths': 'paths for the objects to open',
                 'table_types': 'table types for the objects',
                 'names': 'potential new names for the objects',
                 'kwargs': 'kwargs for the objects to open',
                 'kwargs_widgets': 'widgets representing kwargs for the objects',
                 'scroll': 'scroll area widget',
                 'scroll_layout': 'layout for scroll area',
                 'scroll_widget': 'widget containing scroll area',
                 'main_layout': 'main layout for the window'}

    def __init__(self, files: List[str], parent=None):
        super().__init__(parent)
        self.files = files
        self.button_box = QtWidgets.QDialogButtonBox(QtWidgets.QDialogButtonBox.Ok | QtWidgets.QDialogButtonBox.Cancel)
        self.all_types_combo = QtWidgets.QComboBox(self)
        self.paths = dict()
        self.table_types = dict()
        self.names = dict()
        self.kwargs = dict()
        self.kwargs_widgets = dict()

        self.scroll = QtWidgets.QScrollArea()
        self.scroll_layout = QtWidgets.QGridLayout(self.scroll)
        self.scroll_widget = QtWidgets.QWidget(self.scroll)
        self.main_layout = QtWidgets.QVBoxLayout(self)

        self.init_ui()

    def init_ui(self):
        self.setLayout(self.main_layout)
        self.main_layout.addWidget(self.scroll)
        self.scroll.setWidgetResizable(True)
        self.scroll.setWidget(self.scroll_widget)
        self.scroll_widget.setLayout(self.scroll_layout)
        self.scroll_layout.setSizeConstraint(QtWidgets.QLayout.SetMinAndMaxSize)

        self.button_box.accepted.connect(self.accept)
        self.button_box.rejected.connect(self.reject)
        self.setWindowTitle('Choose table types and names')
        self.scroll_layout.addWidget(
            QtWidgets.QLabel('Please choose a table type (mandatory) and table name (optional) '
                             'for each loaded table\n\n'), 0, 0, 1, 3)
        self.scroll_layout.addWidget(QtWidgets.QLabel('Table paths:'), 1, 0)
        self.scroll_layout.addWidget(QtWidgets.QLabel('Table types:'), 1, 1)
        self.scroll_layout.addWidget(QtWidgets.QLabel('Table names (optional):'), 1, 2)
        self.scroll_layout.addWidget(QtWidgets.QLabel('Additional parameters:'), 1, 3)
        self.scroll_layout.addWidget(self.all_types_combo, 2, 1)

        self.all_types_combo.addItems(list(FILTER_OBJ_TYPES.keys()))
        self.all_types_combo.setCurrentText('Other table')
        self.all_types_combo.currentTextChanged.connect(self.change_all_table_types)

        for i, file in enumerate(self.files, 3):
            self.paths[file] = gui_widgets.PathLineEdit(file, parent=self)
            self.table_types[file] = QtWidgets.QComboBox(self)
            self.table_types[file].addItems(list(FILTER_OBJ_TYPES.keys()))
            self.names[file] = QtWidgets.QLineEdit(self)
            kwargs_widget = QtWidgets.QWidget(self)
            self.kwargs[file] = QtWidgets.QGridLayout(kwargs_widget)

            self.table_types[file].currentTextChanged.connect(functools.partial(self.update_args_ui, file))
            self.table_types[file].setCurrentText('Other table')

            self.scroll_layout.addWidget(self.paths[file], i, 0)
            self.scroll_layout.addWidget(self.table_types[file], i, 1)
            self.scroll_layout.addWidget(self.names[file], i, 2)
            self.scroll_layout.addWidget(kwargs_widget, i, 3)
        self.main_layout.addWidget(self.button_box)

        self.scroll.setMinimumWidth(self.scroll_widget.sizeHint().width() + 150)

    @QtCore.pyqtSlot(str)
    def change_all_table_types(self, new_type: str):
        for file in self.files:
            self.table_types[file].setCurrentText(new_type)

    @QtCore.pyqtSlot(str)
    def update_args_ui(self, file: str):
        # clear previous layout
        gui_widgets.clear_layout(self.kwargs[file])
        self.kwargs_widgets[file] = {}
        func_name = '__init__'
        filter_obj_type = FILTER_OBJ_TYPES[self.table_types[file].currentText()]
        signature = generic.get_method_signature(func_name, filter_obj_type)
        desc, param_desc = io.get_method_docstring(func_name, filter_obj_type)
        self.table_types[file].setToolTip(desc)
        i = 1
        for name, param in signature.items():
            if name in INIT_EXCLUDED_PARAMS:
                continue
            self.kwargs_widgets[file][name] = gui_widgets.param_to_widget(param, name)
            label = QtWidgets.QLabel(f'{name}:', self.kwargs_widgets[file][name])
            if name in param_desc:
                label.setToolTip(param_desc[name])
                help_button = gui_widgets.HelpButton()
                self.kwargs[file].addWidget(help_button, i, 2)
                help_button.set_param_help(name, param_desc[name])

            self.kwargs[file].addWidget(label, i, 0)
            self.kwargs[file].addWidget(self.kwargs_widgets[file][name], i, 1)
            i += 1

    def result(self):
        paths = {file: gui_widgets.get_val_from_widget(widget) for file, widget in self.paths.items()}
        types = {file: gui_widgets.get_val_from_widget(widget) for file, widget in self.table_types.items()}
        names = {file: gui_widgets.get_val_from_widget(widget) for file, widget in self.names.items()}
        kwargs = {file: {key: gui_widgets.get_val_from_widget(val) for key, val in widget_dict.items()} for
                  file, widget_dict in self.kwargs_widgets.items()}
        return paths, types, names, kwargs


class ReactiveTabWidget(QtWidgets.QTabWidget):
    tabRightClicked = QtCore.pyqtSignal(int)
    newTabFromSet = QtCore.pyqtSignal(set, int, str)
    newTabFromFilter = QtCore.pyqtSignal(filtering.Filter, int, str)
    tabClosed = QtCore.pyqtSignal(int)

    def __init__(self, parent=None):
        super().__init__(parent)
        self.tabBar().setDocumentMode(True)
        self.tabBar().setMovable(True)
        self.setTabsClosable(True)
        self.setElideMode(QtCore.Qt.TextElideMode.ElideMiddle)

    def mousePressEvent(self, event: QtGui.QMouseEvent):
        if event.button() == QtCore.Qt.LeftButton:
            super().mousePressEvent(event)
        elif event.button() == QtCore.Qt.RightButton:
            point = event.pos()
            if point.isNull():
                return
            ind = self.tabBar().tabAt(point)
            if ind == -1:
                return
            self.tabRightClicked.emit(ind)

    def new_tab_from_item(self, item: Union[filtering.Filter, set], name: str, tab_id: int):
        if isinstance(item, set):
            self.newTabFromSet.emit(item, tab_id, name)
        elif validation.isinstanceinh(item, filtering.Filter):
            self.newTabFromFilter.emit(item, tab_id, name)
        else:
            raise TypeError(type(item))

    def removeTab(self, index):
        h = self.cornerWidget().height()
        tab_id = self.widget(index).tab_id

        super().removeTab(index)
        self.update()
        if self.count() == 0:
            self.cornerWidget().setMinimumHeight(h)
            self.setMinimumHeight(h)
        self.tabClosed.emit(tab_id)

    def addTab(self, widget: TabPage, a1: str) -> int:
        return super().addTab(widget, a1)

    def currentWidget(self) -> Union[FilterTabPage, SetTabPage]:
        return super().currentWidget()

    def widget(self, index: int) -> Union[FilterTabPage, SetTabPage]:
        return super().widget(index)


class RenameCommand(QtWidgets.QUndoCommand):
    __slots__ = {'prev_name': 'previous name of the tab',
                 'new_name': 'new name of the tab',
                 'prev_id': 'previous ID',
                 'job_id': 'job ID of the rename command',
                 'tab': 'tab widget'}

    def __init__(self, prev_name: str, new_name: str, tab: TabPage, description: str):
        super().__init__(description)
        self.prev_name = prev_name
        self.new_name = new_name
        self.prev_id = tab.tab_id
        self.job_id = JOB_COUNTER.get_id()
        self.tab = tab

    def undo(self):
        self.tab.tabReverted.emit(self.tab.tab_id)
        self.tab._rename(self.prev_name)
        self.tab.tab_id = self.prev_id

    def redo(self):
        self.job_id = JOB_COUNTER.get_id() if self.job_id is None else self.job_id
        self.tab._rename(self.new_name, self.job_id)


class CloseTabCommand(QtWidgets.QUndoCommand):
    __slots__ = {'tab_container': 'ReactiveTabWidget containing the tabs',
                 'tab_index': 'index of the tab to be closed',
                 'tab_icon': 'icon of the tab',
                 'tab_name': 'name of the tab',
                 'obj_type': "object type of the tab's underlying object",
                 'filename': "filename to cache the tab's underlying object under",
                 'kwargs': 'kwargs for the underlying object'}

    def __init__(self, tab_container: ReactiveTabWidget, tab_index: int, description: str):
        super().__init__(description)
        self.tab_container = tab_container
        self.tab_index = tab_index
        self.tab_id = self.tab_container.widget(self.tab_index).tab_id
        self.tab_icon = tab_container.tabIcon(self.tab_index)
        self.tab_name = tab_container.tabText(self.tab_index).rstrip('*')
        self.obj_type = self.tab_container.widget(self.tab_index).obj_type()
        self.filename = self.tab_container.widget(self.tab_index).cache()

        self.kwargs = self.tab_container.widget(self.tab_index).obj_properties()

    def undo(self):
        item = io.load_cached_gui_file(self.filename)
        if isinstance(item, pd.DataFrame):
            item = self.obj_type.from_dataframe(item, self.tab_name, **self.kwargs)

        self.tab_container.new_tab_from_item(item, self.tab_name, self.tab_id)
        self.tab_container.setTabIcon(self.tab_container.currentIndex(), self.tab_icon)

    def redo(self):
        self.tab_container.removeTab(self.tab_index)


class InplaceCommand(QtWidgets.QUndoCommand):
    __slots__ = {'tab': 'tab widget',
                 'prev_job_id': 'previous job ID',
                 'new_job_id': 'new job ID',
                 'func_name': 'name of the function to apply/undo',
                 'args': 'function args',
                 'kwargs': 'function kwargs',
                 'obj_copy': 'copy of the original object',
                 'predecessors': 'ancestor job IDs for the action'}

    def __init__(self, tab: TabPage, func_name: str, args: list, kwargs: dict, description: str,
                 predecessors: List[int]):
        super().__init__(description)
        self.tab = tab
        self.prev_job_id = tab.tab_id
        self.new_job_id = None
        self.func_name = func_name
        self.args = args
        self.kwargs = kwargs
        self.obj_copy = copy.copy(self.tab.obj())
        self.predecessors = predecessors if isinstance(predecessors, list) else []

    def undo(self):
        obj = self.tab.obj()
        del obj
        self.tab.update_obj(copy.copy(self.obj_copy))
        self.tab.update_tab(is_unsaved=True)
        self.tab.tabReverted.emit(self.tab.tab_id)
        self.tab.tab_id = self.prev_job_id

    def redo(self):
        self.new_job_id = JOB_COUNTER.get_id() if self.new_job_id is None else self.new_job_id
        self.tab._apply_function_from_params(self.func_name, self.args, self.kwargs, job_id=self.new_job_id,
                                             predecessors=self.predecessors)


class InplaceCachedCommand(InplaceCommand):
    __slots__ = {'first_pass': 'indicates whether the command was already applied once',
                 'processed_obj': 'object after application of the function',
                 'new_spawn_id': 'tab id of the tab after running redo()'}

    def __init__(self, tab: TabPage, func_name: str, args: list, kwargs: dict, description: str, predecessors: list):
        super().__init__(tab, func_name, args, kwargs, description, predecessors)
        self.first_pass = True
        self.processed_obj = None
        self.new_spawn_id = None

    def redo(self):
        if self.first_pass:
            self.first_pass = False
            super().redo()
        else:
            assert self.new_spawn_id is not None
            source_name = generic.get_method_readable_name(getattr(self.tab.obj(), self.func_name))
            self.new_job_id = JOB_COUNTER.get_id() if self.new_job_id is None else self.new_job_id
            self.tab.update_obj(copy.copy(self.processed_obj))
            self.tab.update_tab(is_unsaved=True)
            self.tab.tab_id = self.new_spawn_id
            mock_partial = functools.partial(lambda *args, **kwargs: None, self.args, self.kwargs)
            worker_output = gui_widgets.WorkerOutput(self.tab.obj(), mock_partial, self.new_job_id,
                                                     self.predecessors + [self.prev_job_id])
            self.tab.functionApplied.emit(worker_output)
            self.tab.itemSpawned.emit(f"'{source_name}'\noutput", self.new_spawn_id, self.new_job_id,
                                      self.processed_obj)

    def undo(self):
        processed_obj = self.tab.obj()
        self.processed_obj = copy.copy(processed_obj)
        self.tab.update_obj(copy.copy(self.obj_copy))
        self.tab.update_tab(is_unsaved=True)
        self.new_spawn_id = self.tab.tab_id
        self.tab.tabReverted.emit(self.tab.tab_id)
        self.tab.tab_id = self.prev_job_id


class SetOpInplacCommand(InplaceCommand):
    def redo(self):
        source_name = generic.get_method_readable_name(getattr(self.tab.obj(), self.func_name))
        self.new_job_id = JOB_COUNTER.get_id() if self.new_job_id is None else self.new_job_id
        is_filter_obj = validation.isinstanceinh(self.tab.obj(), filtering.Filter)
        first_obj = self.tab.obj()
        if not is_filter_obj:
            first_obj = filtering.Filter.from_dataframe(pd.DataFrame(index=first_obj), 'placeholder')
        partial = functools.partial(getattr(first_obj, self.func_name), *self.args, **self.kwargs)
        worker = gui_widgets.Worker(partial, self.new_job_id, self.predecessors + [self.prev_job_id])
        worker.finished.connect(check_run_success)
        worker.finished.connect(self.tab.functionApplied.emit)
        worker.run()

        if not is_filter_obj:
            self.tab.update_obj(first_obj.index_set)
        self.tab.update_tab()
        new_spawn_id = JOB_COUNTER.get_id()
        self.tab.tab_id = new_spawn_id
        self.tab.itemSpawned.emit(f"'{source_name}'\noutput", new_spawn_id, self.new_job_id, self.tab.obj())


class PipelineInplaceCommand(QtWidgets.QUndoCommand):
    __slots__ = {'tab': 'tab object',
                 'pipeline': 'Pipeline to apply',
                 'pipeline_name': 'Pipeline name',
                 'obj_copy': 'copy of the original Filter object of the tab',
                 'prev_job_id': 'previous job ID',
                 'new_job_id': 'new job ID',
                 'predecessors': 'ancestor job IDs for the action'}

    def __init__(self, tab: FilterTabPage, pipeline: filtering.Pipeline, pipeline_name: str, description: str,
                 predecessors: List[int]):
        super().__init__(description)
        self.tab = tab
        self.prev_job_id = tab.tab_id
        self.new_job_id = None
        self.pipeline = pipeline
        self.pipeline_name = pipeline_name
        self.predecessors = predecessors if isinstance(predecessors, list) else []
        self.obj_copy = copy.copy(self.tab.filter_obj)

    def undo(self):
        del self.tab.filter_obj
        self.tab.filter_obj = copy.copy(self.obj_copy)
        self.tab.update_tab(is_unsaved=True)
        self.tab.tabReverted.emit(self.tab.tab_id)
        self.tab.tab_id = self.prev_job_id

    def redo(self):
        self.new_job_id = JOB_COUNTER.get_id() if self.new_job_id is None else self.new_job_id
        self.tab._apply_pipeline(self.pipeline, self.pipeline_name, True, self.new_job_id, self.predecessors)


class MainWindow(QtWidgets.QMainWindow):
    USER_GUIDE_URL = 'https://guyteichman.github.io/RNAlysis/build/user_guide_gui.html'
    TUTORIAL_URL = 'https://guyteichman.github.io/RNAlysis/build/tutorial.html'
    FAQ_URL = 'https://guyteichman.github.io/RNAlysis/build/faq.html'
    BUGS_URL = 'https://github.com/GuyTeichman/RNAlysis/issues/new?assignees=&labels=bug+report&projects=' \
               '&template=bug_report.yaml&title=Bug+Report%3A+'
    FEATURE_URL = 'https://github.com/GuyTeichman/RNAlysis/issues/new?assignees=&labels=feature+request&projects=' \
                  '&template=feature_request.yaml&title=Feature+Request%3A+'
    QUESTION_URL = 'https://github.com/GuyTeichman/RNAlysis/discussions'

    jobQueued = QtCore.pyqtSignal()

    def __init__(self, gather_stdout: bool = True):
        super().__init__()

        self._generate_report = False
        self.report = None
        self.tabs = ReactiveTabWidget(self)

        self.closed_tabs_stack = QtWidgets.QUndoStack(self)
        self.undo_group = QtWidgets.QUndoGroup(self)
        self.tabs.currentChanged.connect(self._change_undo_stack)

        self.undo_view = QtWidgets.QUndoView(self.undo_group)
        self.command_history_dock = QtWidgets.QDockWidget('Command history', self)

        self.add_tab_button = QtWidgets.QToolButton()
        self.add_tab_button.setToolTip('Add New Tab')
        self.add_tab_button.clicked.connect(functools.partial(self.add_new_tab, name=None))
        self.add_tab_button.setText("+")
        self.status_bar = gui_windows.StatusBar(self)
        self.error_window = None

        self.menu_bar = QtWidgets.QMenuBar(self)
        self.tab_contextmenu = None
        self.pipelines: typing.OrderedDict[str, (generic.GenericPipeline, int)] = OrderedDict()
        self.pipeline_window = None

        self.about_window = gui_windows.AboutWindow(self)
        self.settings_window = gui_windows.SettingsWindow(self)
        self.settings_window.styleSheetUpdated.connect(self.update_style_sheet)
        self.set_op_window = None
        self.set_visualization_window = None
        self.quickstart_window = gui_quickstart.QuickStartWizard(self)
        self.cite_window = gui_windows.HowToCiteWindow(self)
        self.enrichment_window = None
        self.enrichment_results = []
        self.task_queue_window = gui_windows.TaskQueueWindow(self)
        self.external_windows = {}

        self.init_ui()
        self.add_new_tab()
        self.init_actions()
        self.init_menus()

        if gather_stdout:
            self.queue_stdout = Queue()
            # create console text read thread + receiver object
            self.thread_stdout_queue_listener = QtCore.QThread()
            self.stdout_receiver = gui_widgets.ThreadStdOutStreamTextQueueReceiver(self.queue_stdout)
            sys.stdout = gui_widgets.WriteStream(self.queue_stdout)

            # attach console text receiver to console text thread
            self.stdout_receiver.moveToThread(self.thread_stdout_queue_listener)
            # connect receiver object to widget for text update
            self.stdout_receiver.queue_stdout_element_received_signal.connect(self.append_text_to_current_console)
            # attach to start / stop methods
            self.thread_stdout_queue_listener.started.connect(self.stdout_receiver.run)
            self.thread_stdout_queue_listener.start()
        else:
            self.queue_stdout = None
            self.stdout_receiver = None
            self.thread_stdout_queue_listener = None

        # init thread execution attributes
        self.current_worker = None
        self.job_thread = QtCore.QThread()
        self.job_queue = Queue()
        self.job_timer = QtCore.QTimer(self)
        self.job_timer.timeout.connect(self.run_threaded_workers)
        self.job_timer.start(3500)
        self.jobQueued.connect(self.run_threaded_workers)

    def show_tutorial(self):
        if settings.get_show_tutorial_settings():
            self.quickstart_window.show()

    @QtCore.pyqtSlot(bool)
    def _toggle_reporting(self, state: bool):
        if state:
            print("Turning on report generation...")
            try:
                from rnalysis.gui import gui_report
                self._generate_report = True
                self.report = gui_report.ReportGenerator()
                cleared = self.clear_session(not (self.tabs.count() == 1 and self.tabs.currentWidget().is_empty()))
                assert cleared
                self.toggle_report_action.setChecked(True)
                print("Report generation turned on. ")
            except ImportError:
                warnings.warn("The RNAlysis 'reports' module is not installed. Please install it and try again. ")
                self._toggle_reporting(False)
            except AssertionError:
                warnings.warn("You must clear the current session before turning report generation on. "
                              "Please clear your current session and try again. ")
                self._toggle_reporting(False)

        else:
            self._generate_report = False
            self.report = None
            self.toggle_report_action.setChecked(False)
            print("Report generation turned off. ")

    def prompt_auto_report_gen(self):
        preset = settings.get_report_gen_settings()
        if preset is None:
            message_box = gui_windows.ReportGenerationMessageBox(self)
            choice, dont_ask_again = message_box.exec()
            if dont_ask_again:
                new_setting = choice
            else:
                new_setting = None
            settings.set_report_gen_settings(new_setting)
        else:
            choice = preset

        self._toggle_reporting(choice)

    def init_ui(self):
        self.setWindowTitle(f'RNAlysis {__version__}')
        self.setGeometry(600, 50, 1050, 800)
        self.update_style_sheet()

        self.tabs.tabRightClicked.connect(self.init_tab_contextmenu)
        self.tabs.tabCloseRequested.connect(self.close_tab)
        self.tabs.newTabFromSet.connect(self.new_tab_from_gene_set)
        self.tabs.newTabFromFilter.connect(self.new_tab_from_filter_obj)
        self.tabs.tabClosed.connect(self.remove_tab_from_report)

        self.command_history_dock.setWidget(self.undo_view)
        self.command_history_dock.setFloating(False)
        self.addDockWidget(QtCore.Qt.RightDockWidgetArea, self.command_history_dock)

        self.setStatusBar(self.status_bar)
        self.task_queue_window.cancelRequested.connect(self.cancel_job)

        self.tabs.setCornerWidget(self.add_tab_button, QtCore.Qt.TopRightCorner)
        self.setCentralWidget(self.tabs)

    @QtCore.pyqtSlot()
    def clear_history(self, confirm_action: bool = True):
        if confirm_action:
            clear_msg = """Are you sure you want to clear all command history?
            This cannot be undone!"""
            reply = QtWidgets.QMessageBox.question(self, 'Clear history',
                                                   clear_msg, QtWidgets.QMessageBox.No, QtWidgets.QMessageBox.Yes)
        else:
            reply = QtWidgets.QMessageBox.Yes

        if reply == QtWidgets.QMessageBox.Yes:
            for stack in self.undo_group.stacks():
                stack.clear()
            self.closed_tabs_stack.clear()

    @QtCore.pyqtSlot(int)
    def init_tab_contextmenu(self, ind: int):
        self.tab_contextmenu = QtWidgets.QMenu(self)

        new_to_right_action = QtWidgets.QAction("New tab to the right")
        new_to_right_action.triggered.connect(
            functools.partial(self.add_new_tab_at, index=ind + 1, name=None, is_set=False))
        self.tab_contextmenu.addAction(new_to_right_action)
        self.tab_contextmenu.addSeparator()

        color_menu = self.tab_contextmenu.addMenu("Change tab &color")
        actions = []
        for color in gui_graphics.COLOR_ICONS:
            this_action = QtWidgets.QAction(color.capitalize())
            this_action.setIcon(gui_graphics.get_icon(color))
            this_action.triggered.connect(functools.partial(self.set_tab_icon, ind, icon_name=color))
            actions.append(this_action)
            color_menu.addAction(this_action)
        reset_action = QtWidgets.QAction("Reset color")
        reset_action.triggered.connect(functools.partial(self.set_tab_icon, ind, icon_name=None))
        color_menu.addAction(reset_action)
        self.tab_contextmenu.addSeparator()
        # sort_menu = self.tab_contextmenu.addMenu("Sort tabs")
        sort_by_name = QtWidgets.QAction("Sort by tab &name")
        sort_by_name.triggered.connect(self.sort_tabs_by_name)
        sort_by_time = QtWidgets.QAction("Sort by creation &time")
        sort_by_time.triggered.connect(self.sort_tabs_by_creation_time)
        sort_by_type = QtWidgets.QAction("Sort by tab type")
        sort_by_type.triggered.connect(self.sort_tabs_by_type)
        sort_by_size = QtWidgets.QAction("Sort by number of features")
        sort_by_size.triggered.connect(self.sort_tabs_by_n_features)
        reverse = QtWidgets.QAction("Reverse tab order")
        reverse.triggered.connect(self.sort_reverse)
        self.tab_contextmenu.addActions([sort_by_name, sort_by_time, sort_by_type, sort_by_size, reverse])
        self.tab_contextmenu.addSeparator()

        close_this_action = QtWidgets.QAction("Close")
        close_this_action.triggered.connect(functools.partial(self.close_tab, ind))
        close_others_action = QtWidgets.QAction("Close other tabs")
        close_others_action.triggered.connect(functools.partial(self.close_other_tabs, ind))
        close_right_action = QtWidgets.QAction("Close tabs to the right")
        close_right_action.triggered.connect(functools.partial(self.close_tabs_to_the_right, ind))
        close_left_action = QtWidgets.QAction("Close tabs to the left")
        close_left_action.triggered.connect(functools.partial(self.close_tabs_to_the_left, ind))
        self.tab_contextmenu.addActions([close_this_action, close_others_action, close_right_action, close_left_action])

        self.tab_contextmenu.exec(QtGui.QCursor.pos())

    def update_style_sheet(self):
        self.setStyleSheet(gui_style.get_stylesheet())

    @QtCore.pyqtSlot(int)
    def _change_undo_stack(self, ind: int):
        if self.tabs.count() == 0:
            self.add_new_tab()
            self.undo_group.setActiveStack(self.undo_group.stacks()[0])
        else:
            stack = self.tabs.widget(ind).undo_stack
            self.undo_group.setActiveStack(stack)

    def sort_reverse(self):
        prev_order = {i: self.tabs.widget(i) for i in range(self.tabs.count())}
        self._sort_by_map(prev_order, reversed(list(prev_order.keys())))

    def sort_tabs_by_name(self):
        tab_names = {}
        for i in range(self.tabs.count()):
            key = [self.tabs.tabText(i), 0]
            while tuple(key) in tab_names:
                key[1] += 1
            tab_names[tuple(key)] = self.tabs.widget(i)
        sorted_names = sorted(tab_names.keys())
        self._sort_by_map(tab_names, sorted_names)

    def sort_tabs_by_creation_time(self):
        widgets = [self.tabs.widget(i) for i in range(self.tabs.count())]
        tab_times = {}
        for widget in widgets:
            key = [widget.creation_time, 0]
            while tuple(key) in tab_times:
                key[1] += 1
            tab_times[tuple(key)] = widget

        sorted_times = sorted(tab_times.keys())
        self._sort_by_map(tab_times, sorted_times)

    def sort_tabs_by_n_features(self):
        widgets = [self.tabs.widget(i) for i in range(self.tabs.count())]
        tab_count = {}
        for widget in widgets:
            if isinstance(widget, FilterTabPage):
                if widget.filter_obj is None:
                    key = [0, 0]
                else:
                    key = [widget.filter_obj.shape[0], 0]
            elif isinstance(widget, SetTabPage):
                key = [len(widget.gene_set), 0]
            else:
                raise TypeError(type(widget))
            while tuple(key) in tab_count:
                key[1] += 1
            tab_count[tuple(key)] = widget
        self._sort_by_map(tab_count, sorted(tab_count.keys(), reverse=True))

    def sort_tabs_by_type(self):
        widgets = [self.tabs.widget(i) for i in range(self.tabs.count())]
        type_to_order = {filtering.CountFilter: 1, filtering.DESeqFilter: 2, filtering.FoldChangeFilter: 3,
                         filtering.Filter: 4, enrichment.FeatureSet: 5, type(None): 6}
        tab_types = {}
        for widget in widgets:
            if isinstance(widget, FilterTabPage):
                key = [type_to_order[type(widget.filter_obj)], 0]
            elif isinstance(widget, SetTabPage):
                key = [type_to_order[type(widget.gene_set)], 0]
            else:
                raise TypeError(type(widget))
            while tuple(key) in tab_types:
                key[1] += 1
            tab_types[tuple(key)] = widget
        self._sort_by_map(tab_types, sorted(tab_types.keys()))

    def _sort_by_map(self, key_map: dict, sorted_keys: typing.Iterable):
        for to_ind, name in enumerate(sorted_keys):
            widget = key_map[name]
            from_ind = self.tabs.indexOf(widget)

            self.tabs.tabBar().moveTab(from_ind, to_ind)

    def add_new_tab_at(self, index: int, name: str = None, is_set: bool = False):
        self.add_new_tab(name, is_set)
        self.tabs.tabBar().moveTab(self.tabs.currentIndex(), index)

    def add_new_tab(self, name: str = None, is_set: bool = False):
        tab_id = JOB_COUNTER.get_id()
        new_undo_stack = QtWidgets.QUndoStack()
        self.undo_group.addStack(new_undo_stack)
        if name is None:
            name = 'New Table'

        else:
            print(name)

        if is_set:
            tab = SetTabPage(name, parent=self.tabs, undo_stack=new_undo_stack, tab_id=tab_id)
        else:
            tab = FilterTabPage(self.tabs, undo_stack=new_undo_stack, tab_id=tab_id)
            tab.startedClustering.connect(self.start_clustering)

        tab.startedJob.connect(self.start_generic_job)
        tab.filterObjectCreated.connect(self.new_tab_from_filter_obj)
        tab.featureSetCreated.connect(self.new_tab_from_gene_set)
        tab.changeIcon.connect(self.set_current_tab_icon)
        tab.tabNameChange.connect(self.rename_tab)
        tab.geneSetsRequested.connect(self.update_gene_sets_widget)

        if self._generate_report:
            tab.functionApplied.connect(self.update_report_from_worker)
            tab.itemSpawned.connect(self.update_report_spawn)
            tab.tabLoaded.connect(self.add_loaded_item_to_report)
            tab.tabReverted.connect(self.remove_tab_from_report)

        self.tabs.addTab(tab, name)
        self.tabs.setCurrentIndex(self.tabs.count() - 1)

    @QtCore.pyqtSlot(object)
    def update_gene_sets_widget(self, widget: gui_widgets.GeneSetComboBox):
        widget.update_gene_sets(self.get_available_objects())

    @QtCore.pyqtSlot(str, bool, int, int)
    def rename_tab(self, new_name: str, is_unsaved: bool, new_id: int = -1, prev_id: int = -1):
        if is_unsaved:
            new_name += '*'
        else:
            new_name = new_name.rstrip('*')
        self.tabs.setTabText(self.tabs.currentIndex(), new_name)
        self.tabs.setTabToolTip(self.tabs.currentIndex(), new_name)

    @QtCore.pyqtSlot()
    def remove_tab_asterisk(self):
        current_name = self.tabs.tabText(self.tabs.currentIndex())
        self.tabs.setTabText(self.tab.currentIndex(), current_name.rstrip('*'))

    def new_table_from_folder(self):
        folder_name = QtWidgets.QFileDialog.getExistingDirectory(self, "Choose directory", str(Path.home()))
        if folder_name:
            filter_obj = filtering.CountFilter.from_folder(folder_name)
            if self.tabs.currentWidget().is_empty():
                self.tabs.removeTab(self.tabs.currentIndex())
            self.new_tab_from_filter_obj(filter_obj, JOB_COUNTER.get_id())

    def new_table_from_folder_htseqcount(self):
        folder_name = QtWidgets.QFileDialog.getExistingDirectory(self, "Choose directory", str(Path.home()))
        if folder_name:
            normalize_answer = QtWidgets.QMessageBox.question(self, 'Normalize values?',
                                                              "Do you want to normalize your count table to "
                                                              "reads-per-million (RPM)?",
                                                              QtWidgets.QMessageBox.Yes | QtWidgets.QMessageBox.No)
            to_normalize = normalize_answer == QtWidgets.QMessageBox.Yes

            filter_obj = filtering.CountFilter.from_folder_htseqcount(folder_name, norm_to_rpm=to_normalize)
            if self.tabs.currentWidget().is_empty():
                self.tabs.removeTab(self.tabs.currentIndex())
            self.new_tab_from_filter_obj(filter_obj, JOB_COUNTER.get_id())

    def load_multiple_files(self):
        dialog = gui_windows.MultiFileSelectionDialog()
        accepted = dialog.exec()
        if accepted == QtWidgets.QDialog.Accepted:
            filenames = dialog.result()
            if len(filenames) > 0:
                window = MultiOpenWindow(filenames, self)
                accepted = window.exec()
                if accepted:
                    paths, types, names, kwargs = window.result()
                    tabs_to_close = None
                    if self.tabs.currentWidget().is_empty():
                        tabs_to_close = self.tabs.currentIndex()

                    for filename in filenames:
                        path = paths[filename]
                        table_type = FILTER_OBJ_TYPES[types[filename]]
                        name = names[filename]
                        filter_obj = table_type(path, **kwargs[filename])
                        if name == '':
                            self.new_tab_from_filter_obj(filter_obj, JOB_COUNTER.get_id())
                        else:
                            self.new_tab_from_filter_obj(filter_obj, JOB_COUNTER.get_id(), name)
                        QtWidgets.QApplication.processEvents()
                    if tabs_to_close is not None:
                        self.tabs.removeTab(tabs_to_close)

    @QtCore.pyqtSlot(filtering.Filter, int, str)
    @QtCore.pyqtSlot(filtering.Filter, int)
    def new_tab_from_filter_obj(self, filter_obj: filtering.Filter, tab_id: int, name: str = None):
        self.add_new_tab(filter_obj.fname.name if name is None else name)
        self.tabs.currentWidget().start_from_filter_obj(filter_obj, tab_id)
        if self._generate_report:
            self.add_loaded_item_to_report(tab_id, self.tabs.currentWidget().name, self.tabs.currentWidget().obj())

    @QtCore.pyqtSlot(str)
    def set_current_tab_icon(self, icon_name: str = None):
        self.set_tab_icon(self.tabs.currentIndex(), icon_name)

    def set_tab_icon(self, ind: int, icon_name: str = None):
        if icon_name is None:
            if isinstance(self.tabs.currentWidget(), SetTabPage):
                obj_type_str = 'set'
            else:
                obj_type = type(self.tabs.currentWidget().filter_obj)
                obj_type_str = 'blank' if obj_type == type(None) else obj_type.__name__
            icon = gui_graphics.get_icon(obj_type_str)
        else:
            icon = gui_graphics.get_icon(icon_name)
        self.tabs.setTabIcon(ind, icon)

    @QtCore.pyqtSlot(set, int, str)
    @QtCore.pyqtSlot(set, int)
    def new_tab_from_gene_set(self, gene_set: enrichment.FeatureSet, tab_id: int, gene_set_name: str = None):
        if gene_set_name is None and validation.isinstanceinh(gene_set, enrichment.FeatureSet):
            gene_set_name = gene_set.set_name
        self.add_new_tab(gene_set_name, is_set=True)
        if tab_id == -1:
            tab_id = self.tabs.currentWidget().tab_id
        self.tabs.currentWidget().start_from_gene_set(tab_id, gene_set)
        if self._generate_report:
            self.add_loaded_item_to_report(tab_id, self.tabs.currentWidget().name, self.tabs.currentWidget().obj())

    def add_loaded_item_to_report(self, item_id: int, item_name: str,
                                  obj: Union[filtering.Filter, enrichment.FeatureSet, generic.GenericPipeline]):
        if item_id in self.report.nodes:
            if self.report.nodes[item_id].is_active:
                return
            self.report.add_node("reactivate node", item_id, [0], '', '', '')
            return
        obj_type = self._get_spawn_type(obj)
        prefix = f'{item_id}_{item_name}' if obj_type in ('Other output', 'Pipeline') else str(item_id)
        filename = self._cache_spawn(obj, prefix)
        desc = self._format_report_desc(obj, filename, obj_type)
        label = f"Loaded {obj_type}\n'{item_name}'"
        self.report.add_node(label, item_id, [0], desc, obj_type, filename)

    def remove_tab_from_report(self, tab_id: int):
        if self._generate_report:
            self.report.trim_node(tab_id)

    def update_report(self, name: str, job_id: int, predecessor_ids: List[int], description: str,
                      node_type: str = 'Other table', filename: Union[str, Path] = None):
        self.report.add_node(name, job_id, predecessor_ids, description, node_type, filename)

    @QtCore.pyqtSlot(gui_widgets.WorkerOutput)
    def update_report_from_worker(self, worker_output: gui_widgets.WorkerOutput):
        if worker_output.raised_exception:
            return
        partial = worker_output.partial
        kwargs = partial.keywords
        name = generic.get_method_readable_name(partial.func)
        desc = f'{partial.func.__name__}({parsing.format_dict_for_display(kwargs)})'.replace('\n', '<br>')
        self.update_report(name, worker_output.job_id, worker_output.predecessor_ids, desc, 'Function')

    def update_report_spawn(self, name: str, spawn_id: int, predecessor_id: int,
                            spawn: Union[filtering.Filter, enrichment.FeatureSet,
                            pd.DataFrame, plt.Figure, generic.GenericPipeline]):
        if spawn is None:
            return
        assert self.is_valid_spawn(spawn), f"Invalid spawn type '{type(spawn)}'!"
        spawn_type = self._get_spawn_type(spawn)
        prefix = f'{spawn_id}_{name}' if spawn_type in ('Other output', 'Pipeline') else str(spawn_id)
        filename = self._cache_spawn(spawn, prefix)
        desc = self._format_report_desc(spawn, filename, spawn_type)
        self.update_report(name, spawn_id, [predecessor_id], desc, spawn_type, filename)

    @staticmethod
    def is_valid_spawn(spawn: object):
        if validation.isinstanceinh(spawn, filtering.Filter) \
            or validation.isinstanceinh(spawn, enrichment.FeatureSet) \
            or isinstance(spawn, (pd.Series, pd.DataFrame, plt.Figure)) \
            or validation.isinstanceinh(spawn, generic.GenericPipeline):
            return True
        return False

    @staticmethod
    def _get_spawn_type(spawn: Union[filtering.Filter, enrichment.FeatureSet, pd.DataFrame, plt.Figure]):
        if validation.isinstanceinh(spawn, filtering.Filter):
            spawn_type = FILTER_OBJ_TYPES_INV.get(type(spawn).__name__, 'Other table')
        elif validation.isinstanceinh(spawn, enrichment.FeatureSet):
            spawn_type = 'Gene set'
        elif isinstance(spawn, (pd.Series, pd.DataFrame, plt.Figure)):
            spawn_type = 'Other output'
        elif validation.isinstanceinh(spawn, generic.GenericPipeline):
            spawn_type = 'Pipeline'
        else:
            raise TypeError(f"Invalid spawn type '{type(spawn)}' for spawn '{spawn}'!")
        return spawn_type

    @staticmethod
    def _cache_spawn(spawn: Union[filtering.Filter, enrichment.FeatureSet, pd.DataFrame], suffix: str):
        obj = spawn
        if validation.isinstanceinh(spawn, filtering.Filter):
            obj = spawn.df
            filename = parsing.slugify(f'{suffix}_{spawn.fname.stem}') + '.csv'
        elif validation.isinstanceinh(spawn, enrichment.FeatureSet):
            filename = parsing.slugify(f'{suffix}_{spawn.set_name}') + '.txt'
        elif isinstance(spawn, (pd.Series, pd.DataFrame)):
            filename = parsing.slugify(suffix) + '.csv'
        elif isinstance(spawn, plt.Figure):
            filename = parsing.slugify(suffix) + '.svg'
        elif validation.isinstanceinh(spawn, generic.GenericPipeline):
            filename = parsing.slugify(suffix) + '.yaml'
            obj = spawn.export_pipeline(None)
        else:
            raise TypeError(f"Invalid spawn type '{type(spawn)}' for spawn '{spawn}'!")
        io.cache_gui_file(obj, filename)
        return filename

    @staticmethod
    def _format_report_desc(obj: Union[filtering.Filter, enrichment.FeatureSet, pd.DataFrame], filename: str,
                            obj_type: str):
        href = Path('data').joinpath(filename).as_posix()
        if validation.isinstanceinh(obj, filtering.Filter):
            html = parsing.df_to_html(obj.df)
            shape = obj.shape if len(obj.shape) >= 2 else (obj.shape[0], 1)
            desc = f'{obj_type}:<br>"{obj.fname.stem}"<br>{html}{shape[0]} rows, {shape[1]} columns'
        elif validation.isinstanceinh(obj, enrichment.FeatureSet):
            items = []
            for i, item in enumerate(obj.gene_set):
                if i == 5:
                    break
                items.append(item)
            items.append('...')
            desc = f'{obj_type}:<br>"{obj.set_name}"<br>{parsing.items_to_html_table(items)}' \
                   f'{len(obj.gene_set)} features'
        elif isinstance(obj, (pd.DataFrame, pd.Series)):
            desc = parsing.df_to_html(obj)
        elif isinstance(obj, plt.Figure):
            desc = f'<img src="data/{filename}" alt="Figure" height="400">'
        elif validation.isinstanceinh(obj, generic.GenericPipeline):
            desc = str(obj).replace('\n', '<br>')
        else:
            raise TypeError(f"Invalid object type '{type(obj)}' of object '{obj}'.")

        desc += f'<br><a href="{href}" target="_blank" rel="noopener noreferrer">Open file</a>'
        return desc

    def delete_pipeline(self):
        if len(self.pipelines) == 0:
            warnings.warn('No Pipelines to delete')
            return

        pipeline_name, status = QtWidgets.QInputDialog.getItem(
            self, 'Delete Pipeline', 'Choose Pipeline to delete:', self.pipelines.keys())
        if status:
            reply = QtWidgets.QMessageBox.question(self, 'Delete Pipeline?',
                                                   "Are you sure you want to delete this Pipeline? "
                                                   "This action cannot be undone!",
                                                   QtWidgets.QMessageBox.No | QtWidgets.QMessageBox.Yes)
            if reply == QtWidgets.QMessageBox.Yes:
                self.pipelines.pop(pipeline_name)
                print(f"Pipeline '{pipeline_name}' deleted successfully")

    def export_pipeline(self):
        if len(self.pipelines) == 0:
            warnings.warn('No Pipelines to export')
            return

        pipeline_name, status = QtWidgets.QInputDialog.getItem(
            self, 'Export Pipeline', 'Choose Pipeline to export:', self.pipelines.keys())
        if status:
            pipeline = self.pipelines[pipeline_name][0]
            self._export_pipeline_from_obj(pipeline_name, pipeline)

    def _export_pipeline_from_obj(self, pipeline_name: str, pipeline: filtering.Pipeline):
        default_name = parsing.slugify(pipeline_name) + '.yaml'
        filename, _ = QtWidgets.QFileDialog.getSaveFileName(self, "Save Pipeline",
                                                            str(Path.home().joinpath(default_name)),
                                                            "YAML file (*.yaml)")
        if filename:
            pipeline.export_pipeline(filename)
            print(f"Successfully saved at {io.get_datetime()} under {filename}")

    def _import_pipeline_from_str(self, pipeline_name: str, content: str):
        d = yaml.safe_load(content)
        if d.get('filter_type') is not None:
            pipeline = filtering.Pipeline.import_pipeline(content)
        elif d['metadata'].get('pipeline_type') == 'single':
            pipeline = fastq.SingleEndPipeline.import_pipeline(content)
        elif d['metadata'].get('pipeline_type') == 'paired':
            pipeline = fastq.PairedEndPipeline.import_pipeline(content)
        else:
            raise TypeError(f"Pipeline file '{content}' is invalid.")
        pipeline_id = JOB_COUNTER.get_id()
        self.pipelines[pipeline_name] = (pipeline, pipeline_id)
        if self._generate_report:
            self.add_loaded_item_to_report(pipeline_id, pipeline_name, pipeline)

    def import_pipeline(self):
        filename, _ = QtWidgets.QFileDialog.getOpenFileName(self, "Choose a Pipeline file", str(Path.home()),
                                                            "YAML file (*.yaml)")
        if filename:
            pipeline_name = str(Path(filename).stem)
            with open(filename) as f:
                content = f.read()
            self._import_pipeline_from_str(pipeline_name, content)

    def import_multiple_gene_sets(self):
        dialog = gui_windows.MultiFileSelectionDialog()
        accepted = dialog.exec()
        if accepted == QtWidgets.QDialog.Accepted:
            filenames = dialog.result()
            tabs_to_close = None
            if len(filenames) > 0 and self.tabs.currentWidget().is_empty():
                tabs_to_close = self.tabs.currentIndex()
            for filename in filenames:
                gene_set = self._filename_to_gene_set(filename)
                self.new_tab_from_gene_set(gene_set, JOB_COUNTER.get_id(), Path(filename).stem)
            if tabs_to_close is not None:
                self.tabs.removeTab(tabs_to_close)

    def import_gene_set(self):
        filename, _ = QtWidgets.QFileDialog.getOpenFileName(self, "Choose a file", str(Path.home()),
                                                            "Text Document (*.txt);;"
                                                            "Comma-Separated Values (*.csv);;"
                                                            "Tab-Separated Values (*.tsv);;"
                                                            "All Files (*)")
        if filename:
            tabs_to_close = None
            if self.tabs.currentWidget().is_empty():
                tabs_to_close = self.tabs.currentIndex()
            gene_set = self._filename_to_gene_set(filename)
            self.new_tab_from_gene_set(gene_set, -1, Path(filename).stem)
            if tabs_to_close is not None:
                self.tabs.removeTab(tabs_to_close)

    @staticmethod
    def _filename_to_gene_set(filename: str):
        if filename.endswith('.csv'):
            gene_set = set(pd.read_csv(filename, index_col=0).index)
        elif filename.endswith('.tsv'):
            gene_set = set(pd.read_csv(filename, index_col=0, sep='\t').index)
        else:
            with open(filename) as f:
                gene_set = {line.strip() for line in f.readlines()}
        return gene_set

    def close_current_tab(self):
        ind = self.tabs.currentIndex()
        self.close_tab(ind)

    def close_other_tabs(self, index: int):
        self.close_tabs_to_the_right(index)
        self.close_tabs_to_the_left(index)

    def close_tabs_to_the_right(self, index: int):
        while self.tabs.count() > index + 1:
            self.close_tab(self.tabs.count() - 1)

    def close_tabs_to_the_left(self, index: int):
        start_count = self.tabs.count()
        while self.tabs.count() > start_count - index:
            self.close_tab(0)

    def close_tab(self, index: int):
        if self.tabs.widget(index).is_empty():
            self.tabs.removeTab(index)
        else:
            self.undo_group.removeStack(self.tabs.widget(index).undo_stack)
            command = CloseTabCommand(self.tabs, index, '"' + self.tabs.tabText(index).rstrip('*') + '"')
            self.closed_tabs_stack.push(command)

    def copy_gene_set(self):
        gene_set = self.tabs.currentWidget().get_index_string()
        QtWidgets.QApplication.clipboard().setText(gene_set)

    def excepthook(self, exc_type, exc_value, exc_tb):  # pragma: no cover
        sys.__excepthook__(exc_type, exc_value, exc_tb)
        self.error_window = gui_windows.ErrorMessage(exc_type, exc_value, exc_tb, self)
        self.error_window.exec()

    def _get_current_console(self):  # pragma: no cover
        if self.pipeline_window is not None and self.pipeline_window.isVisible():
            current_console = self.pipeline_window.get_console()
        else:
            current_console = self.tabs.currentWidget().get_console()
        return current_console

    @QtCore.pyqtSlot(str)
    def append_text_to_current_console(self, text: str):  # pragma: no cover
        current_console = self._get_current_console()
        current_console.append_text(text)

    def add_pipeline(self):
        self.pipeline_window = CreatePipelineWindow(self)
        self.pipeline_window.pipelineSaved.connect(self.save_pipeline)
        self.pipeline_window.pipelineExported.connect(self._export_pipeline_from_obj)
        self.pipeline_window.exec()
        # self.pipeline_window = None

    @QtCore.pyqtSlot(str)
    def edit_pipeline(self, pipeline_name: str):
        assert pipeline_name in self.pipelines, f"Pipeline {pipeline_name} doesn't exist!"
        pipeline = self.pipelines[pipeline_name][0]
        self.pipeline_window = CreatePipelineWindow.start_from_pipeline(pipeline, pipeline_name, self)
        self.pipeline_window.pipelineSaved.connect(self.save_pipeline)
        self.pipeline_window.pipelineExported.connect(self._export_pipeline_from_obj)
        self.pipeline_window.exec()
        self.pipeline_window = None

    @QtCore.pyqtSlot(str, generic.GenericPipeline)
    def save_pipeline(self, pipeline_name: str, pipeline: filtering.Pipeline):
        if pipeline_name in self.pipelines:
            is_new = False
            response = QtWidgets.QMessageBox.question(self, 'Overwrite Pipeline?',
                                                      'A Pipeline with this name already exists. '
                                                      'Are you sure you want to overwrite it?',
                                                      defaultButton=QtWidgets.QMessageBox.No)

        else:
            is_new = True
            response = QtWidgets.QMessageBox.Yes

        if response == QtWidgets.QMessageBox.Yes:
            new_pipeline_id = JOB_COUNTER.get_id()
            if self._generate_report:
                if is_new:
                    self.add_loaded_item_to_report(new_pipeline_id, pipeline_name, pipeline)
                else:
                    self.update_report_spawn(pipeline_name, new_pipeline_id, self.pipelines[pipeline_name][1], pipeline)

            self.pipelines[pipeline_name] = (pipeline, new_pipeline_id)
            print(f"Successfully saved Pipeline '{pipeline_name}'")

    def settings(self):
        self.settings_window.exec()

    def create_action(self, name, triggered_func, checkable=False, checked=False, enabled=True, shortcut=None):
        action = QtWidgets.QAction(name, self)
        action.triggered.connect(triggered_func)
        action.setCheckable(checkable)
        action.setChecked(checked)
        action.setEnabled(enabled)
        if shortcut:
            action.setShortcut(QtGui.QKeySequence(shortcut))
        return action

    def init_actions(self):
        # Table Actions
        self.new_table_action = self.create_action("&New table", functools.partial(self.add_new_tab, name=None),
                                                   shortcut="Ctrl+N")
        self.new_table_from_folder_action = self.create_action("New table from &folder", self.new_table_from_folder)
        self.new_table_from_folder_htseq_action = self.create_action("New table from folder (HTSeq-count output)",
                                                                     self.new_table_from_folder_htseqcount)
        self.new_multiple_action = self.create_action("&Multiple new tables", self.load_multiple_files,
                                                      shortcut="Ctrl+Shift+N")

        # File Actions
        self.save_action = self.create_action("&Save...", self.save_file, shortcut="Ctrl+S")
        self.load_session_action = self.create_action("&Load session...", self.load_session)
        self.save_session_action = self.create_action("Sa&ve session...", self.save_session)
        self.clear_session_action = self.create_action("Clea&r session...", self.clear_session)

        # Settings and Updates
        self.clear_cache_action = self.create_action("&Clear cache...", self.clear_cache)
        self.settings_action = self.create_action("&Settings...", self.settings)
        self.exit_action = self.create_action("&Exit", self.close)
        self.check_update_action = self.create_action("Check for &updates...", self.check_for_updates)

        # Report Actions
        self.toggle_report_action = self.create_action("&Enable report generation", self._toggle_reporting, True, False)
        self.generate_report_action = self.create_action("&Create session report", self.generate_report)

        # Undo and Redo Actions
        self.undo_action = self.undo_group.createUndoAction(self, "Ctrl+Z")
        self.redo_action = self.undo_group.createRedoAction(self, "Ctrl+Shift+Z")
        self.restore_tab_action = self.create_action("Restore tab", self.closed_tabs_stack.undo,
                                                     shortcut="Ctrl+Shift+T")
        self.close_current_action = self.create_action("&Close current tab", self.close_current_tab, shortcut="Ctrl+W")

        # View and History Actions
        self.close_figs_action = self.create_action("Close all &Figures", functools.partial(plt.close, 'all'))
        self.task_queue_action = self.create_action("Task &queue", self.task_queue_window.show)
        self.status_bar.taskQueueRequested.connect(self.task_queue_action.trigger)
        self.show_history_action = self.create_action("Command &History", self.toggle_history_window, checkable=True,
                                                      checked=True)
        self.clear_history_action = self.create_action("Clea&r command history", self.clear_history)

        # Gene Set Actions
        self.copy_action = self.create_action("&Copy Gene Set", self.copy_gene_set, shortcut="Ctrl+C")
        self.import_set_action = self.create_action("&Import Gene Set...", self.import_gene_set,
                                                    shortcut="Ctrl+Shift+I")
        self.import_multiple_sets_action = self.create_action("Import &Multiple Gene Sets...",
                                                              self.import_multiple_gene_sets, shortcut="Ctrl+Shift+M")
        self.export_set_action = self.create_action("&Export Gene Set...", self.export_gene_set,
                                                    shortcut="Ctrl+Shift+E")
        self.set_op_action = self.create_action("Set &Operations...", self.choose_set_op, shortcut="Ctrl+Shift+O")
        self.enrichment_action = self.create_action("Enrichment &Analysis...", self.open_enrichment_analysis,
                                                    shortcut="Ctrl+Shift+A")
        self.set_vis_action = self.create_action("&Visualize Gene Sets...", self.visualize_gene_sets,
                                                 shortcut="Ctrl+Shift+V")

        # External Tool Actions
        self.cutadapt_single_action = self.create_action("&Single-end adapter trimming...",
                                                         functools.partial(self.start_external_window,
                                                                           CutAdaptSingleWindow))
        self.cutadapt_paired_action = self.create_action("&Paired-end adapter trimming...",
                                                         functools.partial(self.start_external_window,
                                                                           CutAdaptPairedWindow))
        self.kallisto_index_action = self.create_action("kallisto build &index...",
                                                        functools.partial(self.start_external_window,
                                                                          KallistoIndexWindow))
        self.kallisto_single_action = self.create_action("&kallisto Single-end RNA-seq quantification...",
                                                         functools.partial(self.start_external_window,
                                                                           KallistoSingleWindow))
        self.kallisto_paired_action = self.create_action("kallisto &Paired-end RNA-seq quantification...",
                                                         functools.partial(self.start_external_window,
                                                                           KallistoPairedWindow))
        self.bowtie2_index_action = self.create_action("Bowtie2 build &index...",
                                                       functools.partial(self.start_external_window,
                                                                         Bowtie2IndexWindow))
        self.bowtie2_single_action = self.create_action("Bowtie2 &Single-end alignment...",
                                                        functools.partial(self.start_external_window,
                                                                          Bowtie2SingleWindow))
        self.bowtie2_paired_action = self.create_action("Bowtie2 &Paired-end alignment...",
                                                        functools.partial(self.start_external_window,
                                                                          Bowtie2PairedWindow))

        self.featurecounts_single_action = self.create_action("&featureCounts Single-end counting...",
                                                              functools.partial(self.start_external_window,
                                                                                FeatureCountsSingleWindow))
        self.featurecounts_paired_action = self.create_action("featureCounts &Paired-end counting...",
                                                              functools.partial(self.start_external_window,
                                                                                FeatureCountsPairedWindow))
        self.sam2fastq_single_action = self.create_action("Convert SAM/BAM to FASTQ (single-end)...",
                                                          functools.partial(self.start_external_window,
                                                                            SamToFastqSingleWindow))
        self.sam2fastq_paired_action = self.create_action("Convert SAM/BAM to FASTQ (paired-end)...",
                                                          functools.partial(self.start_external_window,
                                                                            SamToFastqPairedWindow))
        self.fastq2sam_single_action = self.create_action("Convert FASTQ to SAM (single-end)...",
                                                          functools.partial(self.start_external_window,
                                                                            FastqToSamSingleWindow))
        self.fastq2sam_paired_action = self.create_action("Convert FASTQ to SAM (paired-end)...",
                                                          functools.partial(self.start_external_window,
                                                                            FastqToSamPairedWindow))
        self.convert_sam_action = self.create_action("Convert SAM/BAM to BAM/SAM...",
                                                     functools.partial(self.start_external_window,
                                                                       ConvertSamFormatWindow))
        self.bam_index_action = self.create_action("Create BAM index...",
                                                   functools.partial(self.start_external_window, BamIndexWindow))
        self.sort_sam_action = self.create_action("Sort SAM/BAM...", functools.partial(self.start_external_window,
                                                                                       SortSamWindow))
        self.validate_sam_action = self.create_action("Validate SAM/BAM...",
                                                      functools.partial(self.start_external_window, ValidateSamWindow))
        self.find_duplicates_action = self.create_action("Find PCR/optical &duplicates...",
                                                         functools.partial(self.start_external_window,
                                                                           FindDuplicatesWindow))

        # Conditional Action for Platform-Specific Use
        action_name = "ShortStack small &RNA alignment..." if platform.system() != 'Windows' else "ShortStack small &RNA alignment (not available on Windows)"
        self.shortstack_action = self.create_action(action_name,
                                                    functools.partial(self.start_external_window, ShortStackWindow),
                                                    enabled=platform.system() != 'Windows')

        # Visualization Actions
        self.bar_plot_action = self.create_action("Create enrichment &bar-plot...",
                                                  functools.partial(self.start_external_window, BarPlotWindow))
        self.ontology_graph_action = self.create_action("Visualize &Gene Ontology...",
                                                        functools.partial(self.start_external_window,
                                                                          OntologyGraphWindow))
        self.pathway_graph_action = self.create_action("Visualize &KEGG Pathway...",
                                                       functools.partial(self.start_external_window,
                                                                         PathwayGraphWindow))

        # Help Actions
        self.quick_start_action = self.create_action("&Quick-start guide", self.quickstart_window.show)
        self.user_guide_action = self.create_action("&User Guide",
                                                    functools.partial(self.open_link, self.USER_GUIDE_URL))
        self.tutorial_action = self.create_action("&Tutorial", functools.partial(self.open_link, self.TUTORIAL_URL))
        self.faq_action = self.create_action("&Frequently Asked Questions",
                                             functools.partial(self.open_link, self.FAQ_URL))
        self.bug_report_action = self.create_action("Submit an &issue",
                                                    functools.partial(self.open_link, self.BUGS_URL))
        self.request_feature_action = self.create_action("&Request a feature",
                                                         functools.partial(self.open_link, self.FEATURE_URL))
        self.ask_question_action = self.create_action("Ask a &question",
                                                      functools.partial(self.open_link, self.QUESTION_URL))
        self.check_update_action = self.create_action("Check for &updates...", self.check_for_updates)
        self.about_action = self.create_action("&About", self.about)
        self.cite_action = self.create_action("How to &cite RNAlysis", self.cite)

        # Pipeline Actions
        self.new_pipeline_action = self.create_action("&New Pipeline...", self.add_pipeline, shortcut="Ctrl+Alt+N")
        self.import_pipeline_action = self.create_action("&Import Pipeline...", self.import_pipeline,
                                                         shortcut="Ctrl+Alt+I")
        self.export_pipeline_action = self.create_action("&Export Pipeline...", self.export_pipeline,
                                                         shortcut="Ctrl+Alt+E")
        self.delete_pipeline_action = self.create_action("&Delete Pipeline...", self.delete_pipeline,
                                                         shortcut="Ctrl+Alt+D")

    @QtCore.pyqtSlot()
    def clear_cache(self):
        reply = QtWidgets.QMessageBox.question(self, 'Clear cache?',
                                               'Are you sure you want to clear the <i>RNAlysis</i> cache? '
                                               'This cannot be undone!',
                                               defaultButton=QtWidgets.QMessageBox.No)
        if reply == QtWidgets.QMessageBox.Yes:
            io.clear_gui_cache()
            io.clear_cache()

    @QtCore.pyqtSlot()
    def check_for_updates(self, confirm_updated: bool = True):  # pragma: no cover
        if io.is_rnalysis_outdated():
            # frozen releases of RNAlysis cannot update using pip. new version must be downloaded manually
            if getattr(sys, 'frozen', False) and hasattr(sys, '_MEIPASS'):
                reply = QtWidgets.QMessageBox.question(self, 'A new version is available',
                                                       'A new version of <i>RNAlysis</i> is available! '
                                                       'Do you wish to download it?')
                if reply == QtWidgets.QMessageBox.Yes:
                    url = QtCore.QUrl('https://github.com/GuyTeichman/RNAlysis/releases/latest')
                    if not QtGui.QDesktopServices.openUrl(url):
                        QtGui.QMessageBox.warning(self, 'Connection failed', 'Could not download new version')
                return

            reply = QtWidgets.QMessageBox.question(self, 'A new version is available',
                                                   'A new version of <i>RNAlysis</i> is available! '
                                                   'Do you wish to update?')
            if reply == QtWidgets.QMessageBox.Yes:
                io.update_rnalysis()
                QtCore.QCoreApplication.quit()
                self.deleteLater()
                QtCore.QProcess.startDetached(
                    Path(sys.executable).parent.joinpath('Scripts', 'rnalysis-gui').as_posix(), sys.argv)

        else:
            if confirm_updated:
                _ = QtWidgets.QMessageBox.information(self, 'You are using the latest version of RNAlysis',
                                                      f'Your version of <i>RNAlysis</i> ({__version__}) is up to date!')

    @QtCore.pyqtSlot(bool)
    def toggle_history_window(self, state: bool):
        if state:
            self.command_history_dock.show()
        else:
            self.command_history_dock.close()

    def start_external_window(self, window_type: Type):
        self.external_windows[window_type] = window_type(self)

        func_name = self.external_windows[window_type].func_name
        func = self.external_windows[window_type].func
        self.external_windows[window_type].paramsAccepted.connect(
            functools.partial(self.start_generic_job_from_params, func_name, func))
        self.external_windows[window_type].geneSetsRequested.connect(self.update_gene_sets_widget)

        self.external_windows[window_type].show()

    def save_file(self):
        self.tabs.currentWidget().save_file()

    def export_gene_set(self):
        this_tab = self.tabs.currentWidget()
        if isinstance(this_tab, FilterTabPage):
            filter_obj = this_tab.obj()
            gene_set = filter_obj.index_set if filter_obj is not None else None
        else:
            gene_set = this_tab.obj()
        if gene_set is None:
            warnings.warn('Cannot export an empty gene set')
            return
        default_name = parsing.slugify(self.tabs.tabText(self.tabs.currentIndex()).rstrip("*")) + '.txt'
        filename, _ = QtWidgets.QFileDialog.getSaveFileName(self, "Save gene set",
                                                            str(Path.home().joinpath(default_name)),
                                                            "Text document (*.txt);;"
                                                            "All Files (*)")
        if filename:
            io.save_gene_set(gene_set, filename)
            print(f"Successfully saved at {io.get_datetime()} under {filename}")

    def open_link(self, link: str):
        url = QtCore.QUrl(link)
        if not QtGui.QDesktopServices.openUrl(url):
            QtGui.QMessageBox.warning(self, 'Connection failed', 'Could not open link. Please try again later. ')

    def get_gene_set_by_ind(self, ind: int):
        gene_set = self.tabs.widget(ind).filter_obj if \
            isinstance(self.tabs.currentWidget(), FilterTabPage) else self.tabs.widget(ind).gene_set.gene_set
        return gene_set

    def get_available_objects(self):
        tab_names = self.get_tab_names()
        checked = {}
        available_objects_unique = {}
        for i, name in enumerate(tab_names):
            if name in checked:
                key = f"{name}_{checked[name]}"

            else:
                checked[name] = 1
                key = name
            available_objects_unique[key] = (self.tabs.widget(i), self.tabs.tabIcon(i))
            checked[name] += 1
        return available_objects_unique

    def get_gene_set_by_name(self, name: str):
        ind = self.get_tab_names().index(name)
        return self.get_gene_set_by_ind(ind)

    def choose_set_op(self):
        available_objs = self.get_available_objects()
        self.set_op_window = SetOperationWindow(available_objs, self)
        self.set_op_window.primarySetUsed.connect(self.choose_tab_by_name)
        self.set_op_window.geneSetReturned.connect(self.resolve_set_op)
        self.set_op_window.show()

    @QtCore.pyqtSlot(set, str, list, dict)
    def resolve_set_op(self, output_set: set, output_name: str, ancestor_ids: list, kwargs: dict):
        set_op_id = JOB_COUNTER.get_id()
        new_tab_id = JOB_COUNTER.get_id()
        if self._generate_report:
            self.update_report(output_name.replace(' output', ''), set_op_id, ancestor_ids,
                               parsing.format_dict_for_display(kwargs), 'Function')
            self.update_report_spawn('Set operation output', new_tab_id, set_op_id,
                                     enrichment.FeatureSet(output_set, output_name))

        self.new_tab_from_gene_set(output_set, new_tab_id, output_name)

    def visualize_gene_sets(self):
        available_objs = self.get_available_objects()
        self.set_visualization_window = SetVisualizationWindow(available_objs, self)
        self.set_visualization_window.show()
        if self._generate_report:
            self.set_visualization_window.figureGenerated.connect(self.resolve_set_vis)

    @QtCore.pyqtSlot(plt.Figure, str, list, dict)
    def resolve_set_vis(self, output_fig: plt.Figure, output_name: str, ancestor_ids: list, kwargs: dict):
        func_id = JOB_COUNTER.get_id()
        fig_id = JOB_COUNTER.get_id()
        if self._generate_report:
            self.update_report(output_name, func_id, ancestor_ids, parsing.format_dict_for_display(kwargs), 'Function')
            self.update_report_spawn('Set visualization output', fig_id, func_id, output_fig)

    @QtCore.pyqtSlot(str)
    def choose_tab_by_name(self, set_name: str):
        available_objs = self.get_available_objects()
        for i, name in enumerate(available_objs.keys()):
            if name == set_name:
                self.tabs.setCurrentIndex(i)
                return

    def display_enrichment_results(self, result: pd.DataFrame, gene_set_name: str):
        df_window = gui_windows.DataFrameView(result, f'Enrichment results for set "{gene_set_name}"', self)
        self.enrichment_results.append(df_window)
        df_window.show()

    def open_enrichment_analysis(self):
        self.enrichment_window = EnrichmentWindow(self)
        self.enrichment_window.geneSetsRequested.connect(self.update_gene_sets_widget)
        self.enrichment_window.enrichmentStarted.connect(self.start_enrichment)
        if self._generate_report:
            self.enrichment_window.functionApplied.connect(self.update_report_from_worker)
        self.enrichment_window.show()

    def get_tab_names(self) -> List[str]:
        return [self.tabs.tabText(i).rstrip('*') for i in range(self.tabs.count())]

    def generate_report(self):
        if not self._generate_report:
            QtWidgets.QMessageBox.warning(self, 'Report generation was not enabled!',
                                          'Cannot generate a report since report generation was not enabled. '
                                          'You can enable it from the Settings menu.')
            return
        outdir = QtWidgets.QFileDialog.getExistingDirectory(self, "Save report")
        if outdir:
            outdir = Path(outdir).joinpath('RNAlysis_report')
            if not outdir.exists():
                outdir.mkdir()
            self.report.generate_report(outdir)
            self._save_session_to(outdir.joinpath('data', self.report.ROOT_FNAME))

    def init_menus(self):
        self.setMenuBar(self.menu_bar)
        file_menu = self.menu_bar.addMenu("&File")
        self.new_menu = file_menu.addMenu("&New...")
        self.new_menu.addActions([self.new_table_action, self.new_multiple_action, self.new_table_from_folder_action,
                                  self.new_table_from_folder_htseq_action])
        file_menu.addActions(
            [self.save_action, self.load_session_action, self.save_session_action, self.clear_session_action,
             self.clear_cache_action, self.toggle_report_action, self.generate_report_action, self.settings_action,
             self.exit_action])

        edit_menu = self.menu_bar.addMenu("&Edit")
        edit_menu.addActions([self.undo_action, self.redo_action, self.restore_tab_action, self.close_current_action,
                              self.clear_history_action])

        view_menu = self.menu_bar.addMenu("&View")
        view_menu.addActions([self.show_history_action, self.task_queue_action, self.close_figs_action])

        fastq_menu = self.menu_bar.addMenu("&FASTQ/SAM")

        self.quality_menu = fastq_menu.addMenu('&Quality control')
        self.quality_menu.addActions([self.validate_sam_action])

        # self.preprocess_menu = fastq_menu.addMenu('&Pre-processing')

        self.trimming_menu = fastq_menu.addMenu('Adapter &trimming')
        self.trimming_menu.addActions([self.cutadapt_single_action, self.cutadapt_paired_action])

        self.kallisto_menu = fastq_menu.addMenu("RNA-sequencing &quantification")
        self.kallisto_menu.addActions(
            [self.kallisto_index_action, self.kallisto_single_action, self.kallisto_paired_action])

        self.alignment_menu = fastq_menu.addMenu("Read &alignment")
        self.alignment_menu.addActions(
            [self.bowtie2_index_action, self.bowtie2_single_action, self.bowtie2_paired_action, self.shortstack_action])

        self.convert_menu = fastq_menu.addMenu('&Conversion')
        self.convert_menu.addActions(
            [self.convert_sam_action, self.sam2fastq_single_action, self.sam2fastq_paired_action,
             self.fastq2sam_single_action, self.fastq2sam_paired_action])

        self.process_menu = fastq_menu.addMenu('P&ost-processing')
        self.process_menu.addActions([self.sort_sam_action, self.find_duplicates_action, self.bam_index_action])

        # self.filtering_menu = fastq_menu.addMenu('&Filtering')

        self.count_menu = fastq_menu.addMenu("&Feature counting")
        self.count_menu.addActions([self.featurecounts_single_action, self.featurecounts_paired_action])

        gene_sets_menu = self.menu_bar.addMenu("&Gene sets")
        gene_sets_menu.addActions(
            [self.copy_action, self.import_set_action, self.import_multiple_sets_action,
             self.export_set_action, self.set_op_action, self.set_vis_action, self.enrichment_action,
             self.ontology_graph_action, self.pathway_graph_action, self.bar_plot_action])

        pipeline_menu = self.menu_bar.addMenu("&Pipelines")
        pipeline_menu.addActions([self.new_pipeline_action, self.import_pipeline_action, self.export_pipeline_action,
                                  self.delete_pipeline_action])
        self.edit_pipeline_menu = pipeline_menu.addMenu("&Edit Pipeline")
        self.edit_pipeline_menu.aboutToShow.connect(
            functools.partial(self._populate_pipelines, self.edit_pipeline_menu, self.edit_pipeline,
                              pipeline_arg=False))
        self.apply_pipeline_menu = pipeline_menu.addMenu("&Apply Pipeline")
        self.apply_pipeline_menu.aboutToShow.connect(
            functools.partial(self._populate_pipelines, self.apply_pipeline_menu, self.apply_pipeline))

        help_menu = self.menu_bar.addMenu("&Help")
        help_menu.addActions(
            [self.quick_start_action, self.tutorial_action, self.user_guide_action, self.faq_action,
             self.bug_report_action, self.request_feature_action, self.ask_question_action,
             self.check_update_action, self.about_action, self.cite_action])

    def _populate_pipelines(self, menu: QtWidgets.QMenu, func: Callable, pipeline_arg: bool = True,
                            name_arg: bool = True):
        # Remove the old options from the menu
        menu.clear()
        # Dynamically create the actions
        actions = []
        for name, (pipeline, pipeline_id) in self.pipelines.items():
            action = QtWidgets.QAction(name, self)
            args = []
            if pipeline_arg:
                args.append(pipeline)
            if name_arg:
                args.append(name)
            action.triggered.connect(functools.partial(func, *args))
            actions.append(action)
        # Step 3. Add the actions to the menu
        menu.addActions(actions)

    def apply_pipeline(self, pipeline: generic.GenericPipeline, pipeline_name: str):
        pipeline_id = self.pipelines[pipeline_name][1]
        if isinstance(pipeline, filtering.Pipeline):
            self._apply_table_pipeline(pipeline, pipeline_name, pipeline_id)
        elif validation.isinstanceinh(pipeline, fastq._FASTQPipeline):
            self._apply_fastq_pipeline(pipeline, pipeline_name, pipeline_id)
        else:
            raise TypeError(f"Invalid Pipeline type: {type(pipeline)}")

    def _apply_fastq_pipeline(self, pipeline: fastq._FASTQPipeline, pipeline_name: str, pipeline_id: int):
        if isinstance(pipeline, fastq.SingleEndPipeline):
            dialog = gui_windows.FuncExternalWindow('Pipeline', pipeline.apply_to, None, {'self'}, parent=self)
        else:
            dialog = gui_windows.PairedFuncExternalWindow('Pipeline', pipeline.apply_to, None,
                                                          {'self', 'r1_files', 'r2_files'}, parent=self)
        dialog.paramsAccepted.connect(
            functools.partial(self.start_generic_job_from_params, pipeline_name, pipeline.apply_to,
                              predecessor=pipeline_id))
        dialog.paramsAccepted.connect(dialog.deleteLater)
        dialog.init_ui()
        self.external_windows['fastq_pipeline'] = dialog
        dialog.exec()

    def _apply_table_pipeline(self, pipeline: filtering.Pipeline, pipeline_name: str, pipeline_id: int):
        apply_msg = f"Do you want to apply Pipeline '{pipeline_name}' inplace?"
        reply = QtWidgets.QMessageBox.question(self, f"Apply Pipeline '{pipeline_name}'",
                                               apply_msg,
                                               QtWidgets.QMessageBox.Yes | QtWidgets.QMessageBox.No |
                                               QtWidgets.QMessageBox.Cancel)
        if reply == QtWidgets.QMessageBox.Cancel:
            return
        inplace = reply == QtWidgets.QMessageBox.Yes

        available_objs = self.get_available_objects()
        filtered_available_objs = {}
        for i, (key, val) in enumerate(available_objs.items()):
            if (self.tabs.widget(i).obj_type() == pipeline.filter_type) or (
                pipeline.filter_type == filtering.Filter and issubclass(self.tabs.widget(i).obj_type(),
                                                                        filtering.Filter)):
                filtered_available_objs[key] = val
        window = gui_windows.ApplyTablePipelineWindow(filtered_available_objs, self)
        accepted = window.exec()
        if accepted:
            current_ind = self.tabs.currentIndex()
            chosen_names = window.result()
            for name in chosen_names:
                self.tabs.setCurrentWidget(filtered_available_objs[name][0])
                filtered_available_objs[name][0].apply_pipeline(pipeline, name, pipeline_id, inplace)
            self.tabs.setCurrentIndex(current_ind)

    @QtCore.pyqtSlot()
    def clear_session(self, confirm_action: bool = True) -> bool:
        if confirm_action:
            response = QtWidgets.QMessageBox.question(self, 'Clear session?',
                                                      'Are you sure you want to clear your session? '
                                                      'All unsaved changes will be lost!',
                                                      defaultButton=QtWidgets.QMessageBox.No)
        else:
            response = QtWidgets.QMessageBox.Yes

        if response == QtWidgets.QMessageBox.Yes:
            while self.tabs.count() > 1:
                self.close_tab(0)
            self.close_tab(0)

            self.pipelines = OrderedDict()
            self.clear_history(confirm_action=False)
            self._change_undo_stack(0)
        return response == QtWidgets.QMessageBox.Yes

    def load_session(self):
        session_filename, _ = QtWidgets.QFileDialog.getOpenFileName(self, "Load session",
                                                                    str(Path.home()),
                                                                    "RNAlysis session files (*.rnal);;"
                                                                    "All Files (*)")
        if session_filename:
            items, item_names, item_types, item_properties, pipeline_names, pipeline_files = io.load_gui_session(
                session_filename)

            for pipeline_name, pipeline_file in zip(pipeline_names, pipeline_files):
                self._import_pipeline_from_str(pipeline_name, pipeline_file)

            tab_to_close = None
            if self.tabs.currentWidget().is_empty():
                tab_to_close = self.tabs.currentIndex()

            for item, item_name, item_type, item_property in zip(items, item_names, item_types, item_properties):
                tab_id = JOB_COUNTER.get_id()
                if item_type == 'set':
                    self.new_tab_from_gene_set(item, tab_id, item_name)
                else:
                    try:
                        cls = getattr(filtering, item_type)
                    except AttributeError:
                        raise TypeError(f"Invalid object type in session file: '{item_type}'")
                    obj = cls.from_dataframe(item, item_name, **item_property)
                    self.new_tab_from_filter_obj(obj, tab_id)

                QtWidgets.QApplication.processEvents()

            if tab_to_close is not None:
                self.tabs.removeTab(tab_to_close)

    def _save_session_to(self, session_filename: Union[str, Path]):
        filenames = []
        item_names = []
        item_types = []
        item_properties = []
        for ind in range(self.tabs.count()):
            tab = self.tabs.widget(ind)
            if not tab.is_empty():
                filenames.append(tab.cache())
                item_names.append(self.tabs.tabText(ind).rstrip('*'))
                if isinstance(tab, SetTabPage):
                    item_types.append(set)
                else:
                    item_types.append(tab.obj_type())
                item_properties.append(tab.obj_properties())

        pipeline_names = []
        pipeline_files = []
        for pipeline_name, (pipeline, pipeline_id) in self.pipelines.items():
            pipeline_files.append(pipeline.export_pipeline(filename=None))
            pipeline_names.append(pipeline_name)

        io.save_gui_session(session_filename, filenames, item_names, item_types, item_properties, pipeline_names,
                            pipeline_files)

    def save_session(self):
        default_name = 'Untitled session.rnal'
        session_filename, _ = QtWidgets.QFileDialog.getSaveFileName(self, "Save session",
                                                                    str(Path.home().joinpath(default_name)),
                                                                    "RNAlysis session files (*.rnal);;"
                                                                    "All Files (*)")
        if session_filename:
            self._save_session_to(session_filename)
            print(f"Session saved successfully at {io.get_datetime()} under {session_filename}")

    def about(self):
        self.about_window.exec()

    def cite(self):
        self.cite_window.exec()

    def input(self, message: str = ''):
        dialog = gui_widgets.PathInputDialog(message, parent=self)
        self.external_windows['input'] = dialog
        accepted = dialog.exec()
        if accepted:
            return dialog.result()
        return None

    def closeEvent(self, event):  # pragma: no cover

        quit_msg = "Are you sure you want to close <i>RNAlysis</i>?\n" \
                   "All unsaved progress will be lost"

        reply = QtWidgets.QMessageBox.question(self, 'Close program',
                                               quit_msg, QtWidgets.QMessageBox.No, QtWidgets.QMessageBox.Yes)

        if reply == QtWidgets.QMessageBox.Yes:
            # quit job and STDOUT listener threads
            try:
                self.job_thread.quit()
            except AttributeError:
                pass
            try:
                self.thread_stdout_queue_listener.quit()
            except AttributeError:
                pass
            # clear cache
            io.clear_gui_cache()
            # close all figures
            plt.close("all")

            # close all external windows
            try:
                for window in itertools.chain(self.external_windows.values(), self.enrichment_results):
                    window.close()
                if self.error_window is not None:
                    self.error_window.close()
            except AttributeError:  # Python 3.8/3.9 error on test suite
                pass
            event.accept()
        else:
            event.ignore()

    @QtCore.pyqtSlot(object, object, object)
    def start_generic_job(self, parent_tab: Union[FilterTabPage, None], worker: gui_widgets.Worker,
                          finish_slots: Union[Callable, List[Callable], None]):
        slots = [functools.partial(self.finish_generic_job, parent_tab=parent_tab)] + parsing.data_to_list(finish_slots)
        self.queue_worker(worker, slots)

    def start_generic_job_from_params(self, func_name, func, args, kwargs,
                                      finish_slots: Union[Callable, List[Callable], None],
                                      predecessor: Union[int, None] = None):
        partial = functools.partial(func, *args, **kwargs)
        worker = gui_widgets.Worker(partial, JOB_COUNTER.get_id(), [] if predecessor is None else [predecessor],
                                    func_name)
        self.start_generic_job(None, worker, finish_slots)

    @QtCore.pyqtSlot(gui_widgets.WorkerOutput, object)
    def finish_generic_job(self, worker_output: gui_widgets.WorkerOutput, parent_tab: TabPage = None):
        if worker_output.raised_exception:
            raise worker_output.raised_exception
        if worker_output.result is None or len(worker_output.result) == 0:
            print("Done")
            return
        assert isinstance(worker_output, gui_widgets.WorkerOutput), f"invalid worker output: {worker_output}"
        func_name: str = worker_output.emit_args[0]
        job_id = worker_output.job_id

        if parent_tab is not None:
            parent_tab.update_tab()
            parent_tab.process_outputs(worker_output.result, job_id, func_name)
        else:
            source_name = generic.get_method_readable_name(worker_output.partial.func)
            if self._generate_report:
                self.update_report_from_worker(worker_output)
            for output in parsing.data_to_list(worker_output.result):
                output_id = JOB_COUNTER.get_id()
                if self._generate_report and self.is_valid_spawn(output):
                    self.update_report_spawn(source_name, output_id, job_id, output)

                if validation.isinstanceinh(output, filtering.Filter):
                    self.new_tab_from_filter_obj(output, output_id)
                elif isinstance(output, pd.DataFrame):
                    self.tabs.currentWidget().object_views.append(
                        gui_windows.DataFrameView(output, source_name, self.tabs.currentWidget()))
                    self.tabs.currentWidget().object_views[-1].show()

    @QtCore.pyqtSlot(object, object, object)
    def start_clustering(self, parent_tab: FilterTabPage, worker: gui_widgets.Worker,
                         finish_slot: Union[Callable, None]):
        slots = (functools.partial(self.finish_clustering, parent_tab=parent_tab), finish_slot)
        self.queue_worker(worker, slots)

    @QtCore.pyqtSlot(gui_widgets.WorkerOutput)
    def finish_clustering(self, worker_output: gui_widgets.WorkerOutput, parent_tab: FilterTabPage):
        if worker_output.raised_exception:
            raise worker_output.raised_exception
        if worker_output.result is None or len(worker_output.result) == 0:
            print("Done")
            return
        assert isinstance(worker_output, gui_widgets.WorkerOutput), f"invalid worker output: {worker_output}"

        func_name: str = generic.get_method_readable_name(worker_output.partial.func)
        clustering_runner: clustering.ClusteringRunner = worker_output.result[1]
        return_val: clustering.ClusteringRunner = worker_output.result[0]
        job_id = worker_output.job_id
        figs = clustering_runner.plot_clustering()
        parent_tab.update_tab()
        parent_tab.process_outputs(figs, job_id, func_name)
        parent_tab.process_outputs(return_val, job_id, func_name)

    @QtCore.pyqtSlot(object, object)
    def start_enrichment(self, worker: gui_widgets.Worker, finish_slot: Union[Callable, None]):
        slots = (self.finish_enrichment, finish_slot)
        self.queue_worker(worker, slots)

    @QtCore.pyqtSlot(gui_widgets.WorkerOutput)
    def finish_enrichment(self, worker_output: gui_widgets.WorkerOutput):
        if worker_output.raised_exception:
            raise worker_output.raised_exception
        if len(worker_output.result) == 0:
            print("Done")
            return
        assert isinstance(worker_output, gui_widgets.WorkerOutput), f"invalid worker output: {worker_output}"

        set_name = worker_output.emit_args[0]
        job_id = worker_output.job_id
        results, enrichment_runner = worker_output.result
        self.show()
        outputs = parsing.data_to_list(enrichment_runner.plot_results()) + [results]
        self.display_enrichment_results(results, set_name)

        if self._generate_report:
            spawn_name = f'Enrichment results for set\n"{set_name}"'
            for item in outputs:
                self.update_report_spawn(spawn_name, JOB_COUNTER.get_id(), job_id, item)

    def queue_worker(self, worker: gui_widgets.Worker,
                     output_slots: Union[Callable, Tuple[Callable, ...], None] = None):
        self.job_queue.put((worker, output_slots))
        self.jobQueued.emit()

    @QtCore.pyqtSlot(int, str)
    def cancel_job(self, index: int, func_name: str):
        func_name = func_name[0:func_name.find(' ')]
        if index == 0:
            QtWidgets.QMessageBox.warning(self, "Can't stop a running job!",
                                          "<i>RNAlysis</i> can't stop a task that is currently running. ")
            return
        print(f'Cancelling job "{func_name}"...')
        index -= 1
        modified_queue = parsing.data_to_list(self.job_queue.queue)
        worker = modified_queue.pop(index)
        worker.deleteLater()
        while not self.job_queue.empty():
            self.job_queue.get()
        for item in modified_queue:
            self.job_queue.put(item)
        self.run_threaded_workers()

    def run_threaded_workers(self):  # pragma: no cover
        job_running = (self.job_thread is not None) and self.job_thread.isRunning()
        self.status_bar.update_n_tasks(self.job_queue.qsize() + job_running)
        # if there are no jobs available, don't proceed
        if self.job_queue.qsize() == 0:
            self._update_queue_window(job_running)
            return
        # if there is a job currently running, don't proceed
        if job_running:
            self._update_queue_window(job_running)
            self.status_bar.update_time()
            return

        worker, output_slots = self.job_queue.get()
        # Create a worker object
        if isinstance(self.current_worker, gui_widgets.Worker):
            try:
                self.current_worker.deleteLater()
            except RuntimeError:
                pass
        self.current_worker = worker

        @QtCore.pyqtSlot()
        def alt_tqdm(iter_obj: typing.Iterable = None, desc: str = '', unit: str = '', bar_format: str = '',
                     total: int = None):
            self.current_worker.startProgBar.emit(
                dict(iter_obj=iter_obj, desc=desc, unit=unit, bar_format=bar_format, total=total))
            obj = gui_widgets.AltTQDM(iter_obj, desc=desc, unit=unit, bar_format=bar_format, total=total)
            obj.barUpdate.connect(self.status_bar.move_progress_bar)
            obj.barFinished.connect(self.status_bar.reset_progress)

            self.current_worker.startProgBar.emit(
                dict(iter_obj=iter_obj, desc=desc, unit=unit, bar_format=bar_format, total=total))
            return obj

        @QtCore.pyqtSlot()
        def alt_parallel(n_jobs: int = -1, desc: str = '', unit: str = '', bar_format: str = '',
                         total: int = None, **kwargs):
            self.current_worker.startProgBar.emit(
                dict(iter_obj=None, desc=desc, unit=unit, bar_format=bar_format, total=total))
            print(f"{desc}: started" + (f" {total} jobs\r" if isinstance(total, int) else "\r"))
            obj = gui_widgets.AltParallel(n_jobs=n_jobs, desc=desc, unit=unit, bar_format=bar_format,
                                          total=total, **kwargs)
            obj.barUpdate.connect(self.status_bar.move_progress_bar)
            obj.barFinished.connect(self.status_bar.reset_progress)
            obj.barTotalUpdate.connect(self.status_bar.update_bar_total)

            self.current_worker.startProgBar.emit(
                dict(iter_obj=None, desc=desc, unit=unit, bar_format=bar_format, total=total))
            return obj

        enrichment.enrichment_runner.generic.ProgressParallel = alt_parallel
        generic.ProgressParallel = alt_parallel
        filtering.clustering.generic.ProgressParallel = alt_parallel

        enrichment.enrichment_runner.parsing.tqdm = alt_tqdm
        enrichment.enrichment_runner.io.tqdm = alt_tqdm
        enrichment.enrichment_runner.tqdm = alt_tqdm
        filtering.clustering.tqdm = alt_tqdm
        filtering.tqdm = alt_tqdm
        fastq.tqdm = alt_tqdm

        #  Move worker to the thread
        self.current_worker.moveToThread(self.job_thread)
        #  Connect signals and slots
        self.current_worker.startProgBar.connect(self.start_progress_bar)
        self.job_thread.started.connect(self.current_worker.run)
        for slot in parsing.data_to_list(output_slots):
            if slot is not None:
                self.current_worker.finished.connect(slot)
        self.current_worker.finished.connect(self.job_thread.quit)
        self.current_worker.finished.connect(self.run_threaded_workers)
        self.current_worker.finished.connect(self.current_worker.deleteLater)

        self.job_thread.start()

        self._update_queue_window(True)

    def _update_queue_window(self, job_running: bool):
        jobs = [self.current_worker.partial.func.__name__ + ' (running)'] if job_running else []
        jobs += [item[0].partial.func.__name__ + ' (queued)' for item in self.job_queue.queue]
        self.task_queue_window.update_tasks(jobs)

    def start_progress_bar(self, arg_dict):  # pragma: no cover
        total = arg_dict.get('total', None)
        iter_obj = arg_dict['iter_obj']
        progbar_desc = arg_dict.get('desc', '')
        if total is not None:
            progbar_total = total
        else:
            try:
                progbar_total = len(iter_obj)
            except TypeError:
                progbar_total = 1

        self.status_bar.start_progress(progbar_total, progbar_desc)


def customwarn(message, category, filename, lineno, file=None, line=None):  # pragma: no cover
    sys.stdout.write(warnings.formatwarning(message, category, filename, lineno))


async def run():  # pragma: no cover
    warnings.showwarning = customwarn
    # close built-in splash screen in frozen app version of RNAlysis
    if '_PYIBoot_SPLASH' in os.environ and importlib.util.find_spec("pyi_splash"):
        import pyi_splash
        pyi_splash.close()
    # when using GUI, the joblib parallel backend should always be multiprocessing (since Loky is not supported in freeze mode)
    if FROZEN_ENV:
        parallel_backend('multiprocessing')
    lockfile = QtCore.QLockFile(QtCore.QDir.tempPath() + '/RNAlysis.lock')
    if lockfile.tryLock(100):
        show_app = True
    else:
        show_app = False

    app = QtWidgets.QApplication([])
    app.setDesktopFileName('RNAlysis')
    icon_pth = str(Path(__file__).parent.parent.joinpath('favicon.ico').absolute())
    app.setWindowIcon(QtGui.QIcon(icon_pth))

    if show_app:
        splash = gui_windows.splash_screen()
        app.processEvents()
        base_message = f"<i>RNAlysis</i> version {__version__}:\t"
        splash.showMessage(base_message + 'loading dependencies', QtCore.Qt.AlignBottom | QtCore.Qt.AlignHCenter)

        if io.check_changed_version():
            video_files = gui_quickstart.QuickStartWizard.VIDEO_FILES
            splash.showMessage(base_message + 'validating tutorial videos',
                               QtCore.Qt.AlignBottom | QtCore.Qt.AlignHCenter)
            async for i in io.get_gui_videos(video_files):
                splash.showMessage(base_message + f'getting tutorial videos {i + 1}/{len(video_files)}',
                                   QtCore.Qt.AlignBottom | QtCore.Qt.AlignHCenter)

        splash.showMessage(base_message + 'loading application', QtCore.Qt.AlignBottom | QtCore.Qt.AlignHCenter)
    matplotlib.use('Qt5Agg')
    window = MainWindow()
    sys.excepthook = window.excepthook
    builtins.input = window.input

    try:
        pass
    except ImportError:
        warnings.warn("RNAlysis can perform faster when package 'numba' is installed. \n"
                      "If you want to improve the performance of slow operations on RNAlysis, "
                      "please install package 'numba'. ")

    if show_app:
        window.show()
        window.check_for_updates(False)
        window.show_tutorial()
        window.prompt_auto_report_gen()
        splash.finish(window)
    sys.exit(app.exec())
