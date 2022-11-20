import builtins
import copy
import functools
import hashlib
import itertools
import os
import sys
import time
import typing
import typing_extensions
import warnings
from collections import OrderedDict
from pathlib import Path
from queue import Queue
from typing import List, Tuple, Union, Callable

import matplotlib
import numpy as np
import pandas as pd
from PyQt5 import QtCore, QtWidgets, QtGui
from joblib.externals.loky import get_reusable_executor

from rnalysis import fastq, filtering, enrichment, __version__
from rnalysis.gui import gui_style, gui_widgets, gui_windows, gui_graphics, gui_tutorial
from rnalysis.utils import io, validation, generic, parsing, settings, enrichment_runner, clustering

FILTER_OBJ_TYPES = {'Count matrix': filtering.CountFilter, 'Differential expression': filtering.DESeqFilter,
                    'Fold change': filtering.FoldChangeFilter, 'Other': filtering.Filter}
FILTER_OBJ_TYPES_INV = {val.__name__: key for key, val in FILTER_OBJ_TYPES.items()}
INIT_EXCLUDED_PARAMS = {'self', 'fname', 'suppress_warnings'}


class FuncExternalWindow(gui_widgets.MinMaxDialog):
    IGNORED_WIDGETS = {'help_link'}
    paramsAccepted = QtCore.pyqtSignal(list, dict, object)

    def __init__(self, func_name: str, func: Callable, excluded_params: set, parent=None):
        super().__init__(parent)
        self.func_name = func_name
        self.func = func
        self.signature = generic.get_method_signature(self.func)
        self.desc, self.param_desc = generic.get_method_docstring(self.func)
        self.excluded_params = excluded_params.copy()

        self.widgets = {}
        self.main_layout = QtWidgets.QVBoxLayout(self)

        self.scroll = QtWidgets.QScrollArea()
        self.scroll_widget = QtWidgets.QWidget(self.scroll)
        self.layout = QtWidgets.QGridLayout(self.scroll_widget)

        self.param_group = QtWidgets.QGroupBox(f"1. Set {func_name} parameters")
        self.param_grid = QtWidgets.QGridLayout(self.param_group)
        self.param_widgets = {}

        self.start_button = QtWidgets.QPushButton(f'Start {self.func_name}')
        self.close_button = QtWidgets.QPushButton('Close')

    def init_ui(self):
        self.scroll.setHorizontalScrollBarPolicy(QtCore.Qt.ScrollBarAlwaysOff)
        self.scroll.setWidgetResizable(True)
        self.scroll.setWidget(self.scroll_widget)
        self.layout.setSizeConstraint(QtWidgets.QLayout.SetMinAndMaxSize)

        self.main_layout.addWidget(self.scroll)

        self.layout.addWidget(self.param_group, 0, 0)
        self.layout.addWidget(self.start_button, 1, 0, 1, 2)
        self.layout.addWidget(self.close_button, 2, 0, 1, 2)

        self.start_button.clicked.connect(self.start_analysis)
        self.close_button.clicked.connect(self.close)

        self.init_param_ui()

    def init_param_ui(self):
        i = 0

        for name, param in self.signature.items():
            if name in self.excluded_params:
                continue
            this_desc = self.param_desc.get(name, '')
            self.param_widgets[name] = gui_widgets.param_to_widget(param, name)
            if isinstance(self.param_widgets[name], (gui_widgets.TableColumnPicker, gui_widgets.TableColumnPicker)):
                self.param_widgets[name].add_columns(self.filter_obj.columns)
            elif isinstance(self.param_widgets[name], gui_widgets.ComboBoxOrOtherWidget) and isinstance(
                self.param_widgets[name].other, (gui_widgets.TableColumnPicker, gui_widgets.TableColumnPicker)):
                self.param_widgets[name].other.add_columns(self.filter_obj.columns)
            label = QtWidgets.QLabel(f'{name}:', self.param_widgets[name])
            label.setToolTip(this_desc)
            help_button = gui_widgets.HelpButton()
            self.param_grid.addWidget(help_button, i, 2)
            self.param_grid.addWidget(label, i, 0)
            self.param_grid.addWidget(self.param_widgets[name], i, 1)
            help_button.connect_param_help(name, this_desc)
            i += 1
        self.param_grid.setRowStretch(i, 1)

    def get_analysis_kwargs(self):
        kwargs = {}
        for param_name, widget in self.param_widgets.items():
            if param_name in self.IGNORED_WIDGETS:
                continue
            val = gui_widgets.get_val_from_widget(widget)
            kwargs[param_name] = val
        return kwargs

    def get_analysis_args(self):
        args = []
        return args

    def start_analysis(self):
        args = self.get_analysis_args()
        kwargs = self.get_analysis_kwargs()
        self.showMinimized()
        try:
            self.paramsAccepted.emit(args, kwargs, self.showNormal)
        except Exception as e:
            self.showNormal()
            raise e


class KallistoIndexWindow(FuncExternalWindow):
    EXCLUDED_PARAMS = set()
    IGNORED_WIDGETS = {'help_link'}

    def __init__(self, parent=None):
        super().__init__('Kallisto create index', fastq.kallisto_create_index,
                         self.EXCLUDED_PARAMS, parent)
        self.init_ui()

    def init_ui(self):
        self.setWindowTitle('Kallisto - create transcriptome index')
        super().init_ui()


class KallistoSingleWindow(FuncExternalWindow):
    EXCLUDED_PARAMS = set()
    IGNORED_WIDGETS = {'help_link'}

    def __init__(self, parent=None):
        super().__init__('Kallisto quantify (single-end reads)', fastq.kallisto_quantify_single_end,
                         self.EXCLUDED_PARAMS, parent)
        self.init_ui()

    def init_ui(self):
        self.setWindowTitle('Kallisto single-end quantification setup')
        super().init_ui()


class KallistoPairedWindow(FuncExternalWindow):
    EXCLUDED_PARAMS = {'r1_files', 'r2_files'}
    IGNORED_WIDGETS = {'help_link'}

    def __init__(self, parent=None):
        super().__init__('Kallisto quantify (paired-end reads)', fastq.kallisto_quantify_paired_end,
                         self.EXCLUDED_PARAMS, parent)

        self.pairs_group = QtWidgets.QGroupBox("2. Choose FASTQ file pairs")
        self.pairs_grid = QtWidgets.QGridLayout(self.pairs_group)
        self.pairs_widgets = {}
        self.init_ui()

    def init_ui(self):
        self.setWindowTitle('Kallisto paired-end quantification setup')
        self.layout.addWidget(self.pairs_group, 0, 1)
        self.setMinimumSize(1250, 650)
        super().init_ui()
        self.init_pairs_ui()

    def init_pairs_ui(self):
        self.pairs_widgets['r1_list'] = gui_widgets.OrderedFileList(self)
        self.pairs_widgets['r2_list'] = gui_widgets.OrderedFileList(self)

        self.pairs_grid.addWidget(self.pairs_widgets['r1_list'], 1, 0)
        self.pairs_grid.addWidget(self.pairs_widgets['r2_list'], 1, 1)
        self.pairs_grid.addWidget(QtWidgets.QLabel("<b>R1 files:</b>"), 0, 0)
        self.pairs_grid.addWidget(QtWidgets.QLabel("<b>R2 files:</b>"), 0, 1)

    def get_analysis_kwargs(self):
        kwargs = super().get_analysis_kwargs()
        kwargs['r1_files'] = self.pairs_widgets['r1_list'].get_sorted_names()
        kwargs['r2_files'] = self.pairs_widgets['r2_list'].get_sorted_names()
        return kwargs


class CutAdaptSingleWindow(FuncExternalWindow):
    EXCLUDED_PARAMS = set()
    IGNORED_WIDGETS = {'help_link'}

    def __init__(self, parent=None):
        super().__init__('CutAdapt (single-end reads)', fastq.trim_adapters_single_end, self.EXCLUDED_PARAMS, parent)
        self.init_ui()

    def init_ui(self):
        self.setWindowTitle('CutAdapt single-end adapter trimming setup')
        super().init_ui()


class CutAdaptPairedWindow(FuncExternalWindow):
    EXCLUDED_PARAMS = {'r1_files', 'r2_files'}
    IGNORED_WIDGETS = {'help_link'}

    def __init__(self, parent=None):
        super().__init__('CutAdapt (paired-end reads)', fastq.trim_adapters_paired_end, self.EXCLUDED_PARAMS, parent)

        self.pairs_group = QtWidgets.QGroupBox("2. Choose FASTQ file pairs")
        self.pairs_grid = QtWidgets.QGridLayout(self.pairs_group)
        self.pairs_widgets = {}
        self.init_ui()

    def init_ui(self):
        self.setWindowTitle('CutAdapt paired-end adapter trimming setup')
        self.layout.addWidget(self.pairs_group, 0, 1)
        self.setMinimumSize(1250, 650)
        super().init_ui()
        self.init_pairs_ui()

    def init_pairs_ui(self):
        self.pairs_widgets['r1_list'] = gui_widgets.OrderedFileList(self)
        self.pairs_widgets['r2_list'] = gui_widgets.OrderedFileList(self)

        self.pairs_grid.addWidget(self.pairs_widgets['r1_list'], 1, 0)
        self.pairs_grid.addWidget(self.pairs_widgets['r2_list'], 1, 1)
        self.pairs_grid.addWidget(QtWidgets.QLabel("<b>R1 files:</b>"), 0, 0)
        self.pairs_grid.addWidget(QtWidgets.QLabel("<b>R2 files:</b>"), 0, 1)

    def get_analysis_kwargs(self):
        kwargs = super().get_analysis_kwargs()
        kwargs['r1_files'] = self.pairs_widgets['r1_list'].get_sorted_names()
        kwargs['r2_files'] = self.pairs_widgets['r2_list'].get_sorted_names()
        return kwargs


class DESeqWindow(FuncExternalWindow):
    EXCLUDED_PARAMS = {'self', 'comparisons'}
    IGNORED_WIDGETS = {'help_link', 'load_design'}

    def __init__(self, parent=None):
        super().__init__('DESeq2', filtering.CountFilter.differential_expression_deseq2, self.EXCLUDED_PARAMS, parent)

        self.comparisons = []
        self.design_mat = None

        self.comparisons_group = QtWidgets.QGroupBox("2. Choose pairwise comparisons for DESeq2")
        self.comparisons_grid = QtWidgets.QGridLayout(self.comparisons_group)
        self.comparisons_widgets = {}

        self.init_ui()

    def init_ui(self):
        self.setWindowTitle('DESeq2 differential expression setup')
        self.layout.addWidget(self.comparisons_group, 0, 1)
        super().init_ui()

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
        self.design_mat = io.load_csv(self.param_widgets['design_matrix'].text(), index_col=0)

        if 'picker' in self.comparisons_widgets:
            self.comparisons_grid.removeWidget(self.comparisons_widgets['picker'])
            self.comparisons_widgets['picker'].deleteLater()

        self.comparisons_widgets['picker'] = gui_widgets.ComparisonPickerGroup(self.design_mat, self)
        self.comparisons_grid.addWidget(self.comparisons_widgets['picker'], 0, 0)

    def get_analysis_kwargs(self):
        kwargs = super().get_analysis_kwargs()
        kwargs['comparisons'] = self.comparisons_widgets['picker'].get_comparison_values()
        return kwargs


class ClicomWindow(FuncExternalWindow):
    EXCLUDED_PARAMS = {'self', 'parameter_dicts', 'gui_mode'}
    ADDITIONAL_EXCLUDED_PARAMS = {'power_transform', 'plot_style', 'split_plots', 'return_probabilities', 'gui_mode'}

    def __init__(self, funcs: dict, filter_obj: filtering.Filter, parent=None):
        super().__init__('CLICOM', filtering.CountFilter.split_clicom, self.EXCLUDED_PARAMS, parent)
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

    def init_ui(self):
        super().init_ui()
        self.setWindowTitle('CLICOM clustering setup')
        self.setMinimumWidth(1500)
        self.layout.addWidget(self.setups_group, 0, 1)
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
                             'single_set': {'alpha'},
                             'non_categorical': {'alpha'},
                             None: {},
                             True: {'alpha'},
                             False: {'alpha'}}

    PLOT_ARGS = {'user_defined': {'plot_horizontal'},
                 'go': {'plot_horizontal', 'plot_ontology_graph', 'ontology_graph_format'},
                 'kegg': {'plot_horizontal', 'plot_pathway_graphs', 'pathway_graphs_format'},
                 'non_categorical': {'plot_log_scale', 'plot_style', 'n_bins'}}

    enrichmentStarted = QtCore.pyqtSignal(object, str, object)

    def __init__(self, available_objects: dict, parent=None):
        super().__init__(parent)

        self.available_objects = available_objects

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
            self.widgets[lst] = gui_widgets.MandatoryComboBox('Choose gene set...', self)
            for obj_name in self.available_objects:
                self.widgets[lst].addItem(self.available_objects[obj_name][1], obj_name)
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

    def update_uis(self):
        self.parameters_signature = {}
        self.stats_signature = {}
        self.plot_signature = {}

        analysis_type = self.get_current_analysis_type()
        chosen_func = self.get_current_func()
        signature = generic.get_method_signature(chosen_func)
        func_desc, param_desc = generic.get_method_docstring(chosen_func)
        for name, param in signature.items():
            this_desc = param_desc.get(name, '')
            if name in self.EXCLUDED_PARAMS:
                continue
            elif name in self.PLOT_ARGS[analysis_type]:
                self.plot_signature[name] = (param, this_desc)
            elif name in set.union(*self.STATISTICAL_TEST_ARGS.values()):
                self.stats_signature[name] = (param, this_desc)
            else:
                self.parameters_signature[name] = (param, this_desc)

        self.update_parameters_ui()
        self.update_stats_ui()
        self.update_plot_ui()
        self.widgets['run_button'].setVisible(True)
        self.widgets['run_button'].setDisabled(True)

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
            help_button.connect_param_help(name, desc)

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
            help_button.connect_param_help(name, desc)
            i += 1

    def update_stats_ui(self):
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
                help_button.connect_param_help(name, desc)
                i += 1

    def is_single_set(self):
        stat_test = self._get_statistical_test()
        return stat_test == 'single_set'

    def is_categorical(self):
        analysis_type = self.get_current_analysis_type()
        return analysis_type != 'non_categorical'

    def get_analysis_params(self):
        kwargs = {}
        stat_test = self._get_statistical_test()

        if not self.is_single_set():
            if self.is_categorical():
                kwargs['statistical_test'] = stat_test
            else:
                kwargs['parametric_test'] = stat_test
        gene_set_name = self.widgets['enrichment_list'].currentText()
        gene_set = self.available_objects[gene_set_name][0].obj()

        for param_name, widget in itertools.chain(self.parameter_widgets.items(), self.plot_widgets.items(),
                                                  self.stats_widgets.items()):
            if param_name in {'help_link', 'dataset_radiobox', 'stats_radiobox'}:
                continue
            val = gui_widgets.get_val_from_widget(widget)
            kwargs[param_name] = val

        if not self.is_single_set():
            bg_set_name = self.widgets['bg_list'].currentText()
            bg_set = self.available_objects[bg_set_name][0].obj()
        else:
            bg_set = None
        return gene_set, bg_set, gene_set_name.rstrip('*'), kwargs

    @QtCore.pyqtSlot()
    def run_analysis(self):
        func = self.get_current_func()
        gene_set, bg_set, set_name, kwargs = self.get_analysis_params()
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

            self.enrichmentStarted.emit(partial, set_name, self.showNormal)

        except Exception as e:
            self.showNormal()
            raise e


class SetOperationWindow(gui_widgets.MinMaxDialog):
    SET_OPERATIONS = {'Union': 'union', 'Majority-Vote Intersection': 'majority_vote_intersection',
                      'Intersection': 'intersection', 'Difference': 'difference',
                      'Symmetric Difference': 'symmetric_difference', 'Other': 'other'}
    EXCLUDED_PARAMS = {'self', 'other', 'others', 'return_type'}

    geneSetReturned = QtCore.pyqtSignal(set, str)
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

    def _get_function_params(self):
        set_names = [item.text() for item in self.widgets['set_list'].get_sorted_selection()]
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
        if func_name == 'other':
            output_set = self.widgets['canvas'].get_custom_selection()
            output_name = f"Other set operation output"
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
                                             f'Apply "{func_name}"')
                self.available_objects[primary_set_name][0].undo_stack.push(command)
                output_set = None
            else:
                output_set = getattr(first_obj, func_name)(*other_objs, **kwargs)

        if isinstance(output_set, set):
            self.geneSetReturned.emit(output_set, output_name)
        self.close()


class SetVisualizationWindow(gui_widgets.MinMaxDialog):
    VISUALIZATION_FUNCS = {'Venn Diagram': 'venn_diagram', 'UpSet Plot': 'upset_plot'}
    EXCLUDED_PARAMS = {'objs', 'attr_ref_table_path', 'fig'}

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
            canvas = gui_graphics.EmptyCanvas("Please choose a visualization function to continue")
        else:
            objs_to_plot, kwargs = self._get_function_params()
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
        kwargs = {}
        for param_name, widget in self.parameter_widgets.items():
            if param_name in {'apply_button', 'help_link'}:
                continue
            val = gui_widgets.get_val_from_widget(widget)

            kwargs[param_name] = val
        return objs_to_plot, kwargs

    @QtCore.pyqtSlot()
    def generate_graph(self):
        func_name = self.get_current_func_name()
        objs_to_plot, kwargs = self._get_function_params()
        _ = getattr(enrichment, func_name)(objs_to_plot, **kwargs)


class TabPage(QtWidgets.QWidget):
    tabNameChange = QtCore.pyqtSignal(str, bool)
    tabSaved = QtCore.pyqtSignal()
    changeIcon = QtCore.pyqtSignal(str)

    def __init__(self, parent=None, undo_stack: QtWidgets.QUndoStack = None):
        super().__init__(parent)
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

        self.stdout_group = QtWidgets.QGroupBox('Log')

        self.splitter = QtWidgets.QSplitter(QtCore.Qt.Vertical)
        self.sup_layout.addWidget(self.splitter)
        self.splitter.addWidget(self.scroll)
        self.splitter.addWidget(self.stdout_group)
        self.splitter.setStretchFactor(0, 1)

        # self.sup_layout.addStretch(1)
        self.stdout_grid = QtWidgets.QGridLayout(self.stdout_group)
        self.stdout_widgets = {}

        self.init_stdout_ui()

    def obj(self):
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

    def init_stdout_ui(self):
        self.stdout_widgets['text_edit_stdout'] = gui_widgets.StdOutTextEdit(self)
        self.stdout_widgets['text_edit_stdout'].setStyleSheet("""QTextEdit {background: #dddddd;}""")
        self.stdout_grid.addWidget(self.stdout_widgets['text_edit_stdout'], 0, 0, 3, 4)

    def get_console(self):
        return self.stdout_widgets['text_edit_stdout']

    @QtCore.pyqtSlot()
    def rename(self, new_name: str = None):
        if new_name is None:
            new_name = self.overview_widgets['table_name'].text()
        prev_name = self.get_tab_name()
        command = RenameCommand(prev_name, new_name, self, f'Rename "{prev_name}" to "{new_name}"')
        self.undo_stack.push(command)

    def _rename(self, new_name: str = None):
        self.tabNameChange.emit(new_name, True)
        self.overview_widgets['table_name_label'].setText(f"Table name: '<b>{new_name}</b>'")
        self.overview_widgets['table_name'].setText('')
        self.name = new_name.rstrip('*')

    def cache(self):
        raise NotImplementedError

    def get_tab_name(self):
        return self.name.rstrip("*")

    # custom method to write anything printed out to console/terminal to my QTextEdit widget via append function.
    def output_terminal_written(self, text):
        self.stdout_widgets['stdout'].append(text)


class SetTabPage(TabPage):
    def __init__(self, set_name: str, gene_set: typing.Union[set, enrichment.FeatureSet] = None, parent=None,
                 undo_stack: QtWidgets.QUndoStack = None):
        super().__init__(parent, undo_stack)
        if gene_set is None:
            gene_set = enrichment.FeatureSet(set(), set_name)
        elif isinstance(gene_set, set):
            gene_set = enrichment.FeatureSet(gene_set, set_name)
        self.gene_set = gene_set
        self.name = set_name
        self.overview_group = QtWidgets.QGroupBox('Data overview')
        self.overview_grid = QtWidgets.QGridLayout(self.overview_group)
        self.overview_widgets = {}
        self.init_overview_ui(set_name)

    def obj(self):
        return self.gene_set.gene_set

    def update_obj(self, obj: set):
        self.update_gene_set(obj)

    def obj_type(self):
        return type(self.gene_set.gene_set)

    def get_index_string(self):
        return "\n".join(self.obj())

    def init_overview_ui(self, set_name: str):
        this_row = 0
        self.layout.insertWidget(0, self.overview_group)
        self.layout.addStretch(1)
        self.overview_widgets['table_name_label'] = QtWidgets.QLabel(f"Gene set name: '<b>{set_name}</b>'")
        self.overview_widgets['table_name_label'].setWordWrap(True)

        self.overview_widgets['preview'] = QtWidgets.QListWidget()
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
        self.overview_grid.addWidget(self.overview_widgets['preview'], this_row, 0, 3, 4)
        this_row += 3

        self.overview_widgets['save_button'] = QtWidgets.QPushButton('Save gene set')
        self.overview_widgets['save_button'].clicked.connect(self.save_file)
        self.overview_grid.addWidget(self.overview_widgets['save_button'], this_row, 3, 2, 1)

        self.overview_widgets['shape'] = QtWidgets.QLabel()
        self.overview_grid.addWidget(self.overview_widgets['shape'], this_row, 0, 1, 2)
        self.update_set_shape()

        self.overview_widgets['view_button'] = QtWidgets.QPushButton('View full gene set')
        self.overview_widgets['view_button'].clicked.connect(self.view_full_gene_set)
        self.overview_grid.addWidget(self.overview_widgets['view_button'], this_row, 2, 2, 1)
        this_row += 2

        self.overview_grid.addWidget(QtWidgets.QWidget(self), this_row, 0)
        self.overview_grid.addWidget(QtWidgets.QWidget(self), 0, 4)
        self.overview_grid.setRowStretch(this_row, 1)
        self.overview_grid.setColumnStretch(4, 1)

    def view_full_gene_set(self):
        set_window = gui_windows.GeneSetView(self.gene_set.gene_set, self.get_tab_name())
        self.overview_widgets['full_table_view'] = set_window
        set_window.show()

    def save_file(self):
        default_name = self.get_tab_name().rstrip("*") + '.txt'
        filename, _ = QtWidgets.QFileDialog.getSaveFileName(self, "Save gene set",
                                                            str(Path.home().joinpath(default_name)),
                                                            "Text document (*.txt);;"
                                                            "All Files (*)")
        if filename:
            self.gene_set.save_txt(filename)
            print(f"Successfully saved at {io.get_datetime()} under {filename}")
            self.tabSaved.emit()

    @QtCore.pyqtSlot()
    def _rename(self, new_name: str = None):
        super()._rename(new_name)
        self.gene_set.change_set_name(new_name)

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

    def update_gene_set(self, gene_set: set):
        self.gene_set.gene_set = gene_set
        self.update_tab()

    def update_tab(self, is_unsaved: bool = True):
        self.update_set_shape()
        self.update_set_preview()
        self.changeIcon.emit('set')


class FuncTypeStack(QtWidgets.QWidget):
    EXCLUDED_PARAMS = {'self', 'gui_mode'}
    NO_FUNC_CHOSEN_TEXT = "Choose a function..."
    funcSelected = QtCore.pyqtSignal(bool)

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
        txt = f"Choose a function from this list to read its description. "
        self.func_combo.setToolTip(txt)
        self.func_help_button.connect_desc_help(txt)

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
        desc, param_desc = generic.get_method_docstring(chosen_func_name, self.filter_obj)
        self.func_combo.setToolTip(desc)
        self.func_help_button.connect_param_help(chosen_func_name, desc)

        i = 1
        for name, param in signature.items():
            if name in self.excluded_params:
                continue
            self.parameter_widgets[name] = gui_widgets.param_to_widget(param, name, pipeline_mode=self.pipeline_mode)
            if not self.pipeline_mode:
                if isinstance(self.parameter_widgets[name],
                              (gui_widgets.TableColumnPicker, gui_widgets.TableColumnPicker)):
                    self.parameter_widgets[name].add_columns(self.filter_obj.columns)
                elif isinstance(self.parameter_widgets[name], gui_widgets.ComboBoxOrOtherWidget) and isinstance(
                    self.parameter_widgets[name].other, (gui_widgets.TableColumnPicker, gui_widgets.TableColumnPicker)):
                    self.parameter_widgets[name].other.add_columns(self.filter_obj.columns)

            label = QtWidgets.QLabel(f'{name}:', self.parameter_widgets[name])
            if name in param_desc:
                label.setToolTip(param_desc[name])
                help_button = gui_widgets.HelpButton()
                self.parameter_grid.addWidget(help_button, i, 2)
                help_button.connect_param_help(name, param_desc[name])

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

    def get_function_params(self):
        func_params = {}
        for param_name, widget in self.parameter_widgets.items():
            if param_name in {'function_combo', 'help_link'}:
                continue
            val = gui_widgets.get_val_from_widget(widget)

            func_params[param_name] = val
        return func_params

    def get_function_name(self):
        readable_name = self.func_combo.currentText()
        if readable_name == self.NO_FUNC_CHOSEN_TEXT:
            return self.NO_FUNC_CHOSEN_TEXT
        name = self.funcs[readable_name]
        return name


class FilterTabPage(TabPage):
    EXCLUDED_FUNCS = {'union', 'intersection', 'majority_vote_intersection', 'difference', 'symmetric_difference',
                      'from_folder', 'save_txt', 'save_csv', 'from_dataframe', 'print_features'}
    CLUSTERING_FUNCS = {'split_kmeans': 'K-Means', 'split_kmedoids': 'K-Medoids',
                        'split_hierarchical': 'Hierarchical (Agglomerative)', 'split_hdbscan': 'HDBSCAN',
                        'split_clicom': 'CLICOM (Ensemble)'}
    SUMMARY_FUNCS = {'describe', 'head', 'tail', 'biotypes', 'print_features'}
    GENERAL_FUNCS = {'sort', 'transform', 'translate_gene_ids', 'differential_expression_deseq2', 'fold_change'}
    THREADED_FUNCS = {'translate_gene_ids', 'differential_expression_deseq2', 'filter_by_kegg_annotations',
                      'filter_by_go_annotations'}
    filterObjectCreated = QtCore.pyqtSignal(object)
    startedClustering = QtCore.pyqtSignal(object, str, object)
    startedJob = QtCore.pyqtSignal(object, str, object)

    def __init__(self, parent=None, undo_stack: QtWidgets.QUndoStack = None):
        super().__init__(parent, undo_stack)
        self.filter_obj = None
        self.df_views = []

        self.basic_group = QtWidgets.QGroupBox('Load a data table')
        self.basic_grid = QtWidgets.QGridLayout(self.basic_group)
        self.basic_widgets = {}
        self.basic_param_container = QtWidgets.QWidget(self)
        self.basic_param_widgets = {}
        self.basic_param_grid = QtWidgets.QGridLayout(self.basic_param_container)

        self.overview_group = QtWidgets.QGroupBox('Data overview')
        self.overview_grid = QtWidgets.QGridLayout(self.overview_group)
        self.overview_widgets = {}

        self.function_group = QtWidgets.QGroupBox('Apply functions')
        self.function_grid = QtWidgets.QGridLayout(self.function_group)
        self.function_widgets = {}

        self.stack = QtWidgets.QStackedWidget(self)
        self.button_box = QtWidgets.QButtonGroup(self)
        self.stack_buttons = []
        self.stack_widgets = {}

        self.clicom_window = None
        self.deseq_window = None

        self.init_basic_ui()

    def obj(self):
        return self.filter_obj

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
    def _rename(self, new_name: str = None):
        super()._rename(new_name)
        self.filter_obj.fname = Path(
            os.path.join(str(self.filter_obj.fname.parent), f"{new_name}{self.filter_obj.fname.suffix}"))

    def cache(self):
        base_str = str(time.time_ns()) + str(self.filter_obj.fname) + str(len(self.filter_obj.shape))
        hex_hash = hashlib.sha1(base_str.encode('utf-8')).hexdigest()
        filename = f"{hex_hash}.csv"
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

        self.overview_widgets['preview'] = QtWidgets.QTableView()
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
        desc, param_desc = generic.get_method_docstring(func_name, filter_obj_type)
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
                help_button.connect_param_help(name, param_desc[name])

            self.basic_param_grid.addWidget(label, i, 0)
            self.basic_param_grid.addWidget(self.basic_param_widgets[name], i, 1)
            i += 1

    def init_basic_ui(self):
        self.layout.insertWidget(0, self.basic_group)
        self.basic_widgets['table_type_combo'] = QtWidgets.QComboBox()
        self.basic_widgets['table_type_combo'].addItems(FILTER_OBJ_TYPES.keys())
        self.basic_widgets['table_type_combo'].currentIndexChanged.connect(self.update_basic_ui)
        self.basic_widgets['table_type_combo'].setCurrentText("Other")

        self.basic_widgets['start_button'] = QtWidgets.QPushButton('Start')
        self.basic_widgets['start_button'].clicked.connect(self.start)
        self.basic_widgets['start_button'].setEnabled(False)

        self.basic_widgets['file_path'] = gui_widgets.PathLineEdit()
        self.basic_widgets['file_path'].textChanged.connect(self._change_start_button_state)

        self.basic_widgets['table_name'] = QtWidgets.QLineEdit()

        self.basic_widgets['file_label'] = QtWidgets.QLabel('Choose a file:')
        self.basic_widgets['type_label'] = QtWidgets.QLabel('Choose table type:')
        self.basic_widgets['name_label'] = QtWidgets.QLabel('Name your table (optional):')

        self.basic_widgets['apply_button'] = QtWidgets.QPushButton('Apply')
        self.basic_widgets['apply_button'].clicked.connect(self.apply_function)
        self.layout.insertWidget(1, self.basic_widgets['apply_button'])
        self.basic_widgets['apply_button'].setVisible(False)

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

    @QtCore.pyqtSlot(int)
    def _update_stack_status(self, ind: int):
        self.stack.widget(ind).check_selection_status()

    def init_function_ui(self):
        self.layout.insertWidget(2, self.function_group)
        sorted_actions = self.get_all_actions()

        self.stack.addWidget(QtWidgets.QWidget(self))  # start the stack empty
        for i, action_type in enumerate(sorted_actions):
            bttn = QtWidgets.QPushButton(action_type)
            bttn.setCheckable(True)
            bttn.setStyleSheet('''QPushButton::checked {background-color : purple;
                                                        color: white;
                                                        border: 1px solid #ba32ba;
                                                        border-radius: 4px;}''')
            self.stack_widgets[action_type] = FuncTypeStack(sorted_actions[action_type], self.filter_obj)
            self.stack_widgets[action_type].funcSelected.connect(self.basic_widgets['apply_button'].setVisible)
            self.stack_widgets[action_type].funcSelected.connect(self._check_for_special_functions)
            self.stack.addWidget(self.stack_widgets[action_type])
            bttn.clicked.connect(functools.partial(self.stack.setCurrentIndex, i + 1))
            self.button_box.addButton(bttn)
            self.stack_buttons.append(bttn)
            self.function_grid.addWidget(bttn, 0, i)
        self.stack.currentChanged.connect(self._update_stack_status)

        self.function_grid.addWidget(self.stack, 1, 0, 1, i + 1)

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

    def view_full_dataframe(self):
        df_window = gui_windows.DataFrameView(self.filter_obj.df, self.name)
        self.overview_widgets['full_table_view'] = df_window
        df_window.show()

    def update_table_preview(self):
        if self.is_empty():
            return
        model = gui_windows.DataFramePreviewModel(self.filter_obj.df)
        self.overview_widgets['preview'].setModel(model)
        self.overview_widgets['preview'].resizeColumnsToContents()
        self.overview_widgets['preview'].resizeRowsToContents()

        self.overview_widgets['preview'].setFixedWidth(
            self.overview_widgets['preview'].horizontalHeader().length() + self.overview_widgets[
                'preview'].verticalHeader().width())

        self.overview_widgets['preview'].setFixedHeight(
            self.overview_widgets['preview'].verticalHeader().length() + self.overview_widgets[
                'preview'].horizontalHeader().height())

    def apply_function(self):
        this_stack: FuncTypeStack = self.stack.currentWidget()
        func_name = this_stack.get_function_name()
        func_params = this_stack.get_function_params()
        if func_params.get('inplace', False):
            if func_name in self.THREADED_FUNCS:
                command = InplaceCachedCommand(self, func_name, args=[], kwargs=func_params,
                                               description=f'Apply "{func_name}"')
            else:
                command = InplaceCommand(self, func_name, args=[], kwargs=func_params,
                                         description=f'Apply "{func_name}"')
            self.undo_stack.push(command)
        else:
            self._apply_function_from_params(func_name, args=[], kwargs=func_params)

    def update_tab(self, is_unsaved: bool = True):
        if self.is_empty():
            return
        self.update_filter_obj_shape()
        self.update_table_preview()
        self.tabNameChange.emit(self.filter_obj.fname.stem, is_unsaved)
        self.name = str(self.filter_obj.fname.stem)
        self.update_table_name_label()

    def _apply_function_from_params(self, func_name, args: list, kwargs: dict, finish_slot=None):
        # since clustering functions can be computationally intensive, start it in another thread
        if func_name in self.CLUSTERING_FUNCS:
            kwargs['gui_mode'] = True
            partial = functools.partial(getattr(self.filter_obj, func_name), *args, **kwargs)
            self.startedClustering.emit(partial, func_name, finish_slot)
            return
        elif func_name in self.THREADED_FUNCS and (not kwargs.get('inplace', False)):

            partial = functools.partial(getattr(self.filter_obj, func_name), *args, **kwargs)
            self.startedJob.emit(partial, func_name, finish_slot)
            return

        prev_name = self.get_tab_name()
        result = getattr(self.filter_obj, func_name)(*args, **kwargs)
        self.update_tab(prev_name != self.filter_obj.fname.name)
        self.process_outputs(result, func_name)

    def process_outputs(self, outputs, source_name: str = ''):
        if validation.isinstanceinh(outputs, (filtering.Filter, fastq.filtering.Filter)):
            self.filterObjectCreated.emit(outputs)
        elif isinstance(outputs, pd.DataFrame):
            self.df_views.append(gui_windows.DataFrameView(outputs, source_name))
            self.df_views[-1].show()
        elif isinstance(outputs, np.ndarray):
            df = pd.DataFrame(outputs)
            self.process_outputs(df, source_name)

        elif isinstance(outputs, (tuple, list)):
            if validation.isinstanceiter_inh(outputs, filtering.Filter):
                dialog = MultiKeepWindow(outputs, self)
                dialog.accepted.connect(functools.partial(self._multi_keep_window_accepted, dialog, source_name))
                dialog.exec()
            else:
                for output in outputs:
                    self.process_outputs(output, source_name)
        elif isinstance(outputs, dict):
            tab_name = self.get_tab_name()
            for this_src_name, output in outputs.items():
                self.process_outputs(output, f"{this_src_name} {tab_name}")

    def _multi_keep_window_accepted(self, dialog: QtWidgets.QDialog, source_name: str):
        kept_outputs = dialog.result()
        for output in kept_outputs:
            self.process_outputs(output, source_name)

    def apply_pipeline(self, pipeline: filtering.Pipeline, pipeline_name: str, inplace: bool):
        if inplace:
            command = PipelineInplaceCommand(self, pipeline, description=f'Apply Pipeline "{pipeline_name}"')
            self.undo_stack.push(command)
        else:
            self._apply_pipeline(pipeline, inplace=False)

    def _apply_pipeline(self, pipeline, inplace: bool):
        prev_name = self.filter_obj.fname.name
        result = pipeline.apply_to(self.filter_obj, inplace)
        self.update_tab(prev_name != self.filter_obj.fname.name)
        self.process_outputs(result, f'Pipeline on {self.get_tab_name()}')

    def get_index_string(self):
        if self.is_empty():
            return ''
        return self.filter_obj.index_string

    def get_all_actions(self):
        assert self.filter_obj is not None, "No table was loaded!"
        all_methods = dir(self.filter_obj)
        public_methods = [mthd for mthd in all_methods if
                          (not mthd.startswith('_')) and (callable(getattr(type(self.filter_obj), mthd))) and (
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

    def save_file(self):
        if self.filter_obj is None:
            warnings.warn("Cannot save an empty tab!")
            return
        default_name = str(self.filter_obj.fname).rstrip("*")
        filename, _ = QtWidgets.QFileDialog.getSaveFileName(self, "Save filtering result",
                                                            str(Path.home().joinpath(default_name)),
                                                            "Comma-Separated Values (*.csv);;"
                                                            "Tab-Separated Values (*.tsv);;"
                                                            "All Files (*)")
        if filename:
            self.filter_obj.save_csv(filename)
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
        self.tabNameChange.emit(new_name, False)

        self.init_overview_ui()
        self.init_function_ui()

        gui_widgets.clear_layout(self.basic_grid)
        self.layout.removeWidget(self.basic_group)
        self.basic_group.deleteLater()

        self.changeIcon.emit(type(self.filter_obj).__name__)

    def start_from_filter_obj(self, filter_obj: filtering.Filter, name: str = None):
        self.filter_obj = filter_obj
        self.basic_group.setVisible(False)
        self.init_overview_ui()
        self.init_function_ui()
        if name is not None:
            self._rename(name)
        print(self.filter_obj)
        self.changeIcon.emit(type(self.filter_obj).__name__)


class CreatePipelineWindow(gui_widgets.MinMaxDialog, FilterTabPage):
    pipelineSaved = QtCore.pyqtSignal(str, filtering.Pipeline)
    pipelineExported = QtCore.pyqtSignal(str, filtering.Pipeline)

    def __init__(self, parent=None):
        super().__init__(parent)
        self.setLayout(self.layout)
        self.setWindowTitle(f'Create new Pipeline')
        self.setGeometry(500, 200, 800, 800)
        self.pipeline = None
        self.is_unsaved = False

    def init_basic_ui(self):
        self.layout.insertWidget(0, self.basic_group)
        self.basic_group.setTitle("Choose data table type for Pipeline")

        self.basic_widgets['table_type_combo'] = QtWidgets.QComboBox()
        self.basic_widgets['table_type_combo'].addItems(FILTER_OBJ_TYPES.keys())
        self.basic_widgets['table_type_combo'].setCurrentText('Other')

        self.basic_widgets['pipeline_name'] = QtWidgets.QLineEdit()
        self.basic_widgets['pipeline_name'].setText('New Pipeline')

        self.basic_widgets['name_label'] = QtWidgets.QLabel('Name your Pipeline:')

        self.basic_widgets['start_button'] = QtWidgets.QPushButton('Start')
        self.basic_widgets['start_button'].clicked.connect(self.start)

        self.basic_widgets['apply_button'] = QtWidgets.QPushButton('Add to Pipeline')
        self.basic_widgets['apply_button'].clicked.connect(self.apply_function)
        self.layout.insertWidget(1, self.basic_widgets['apply_button'])
        self.basic_widgets['apply_button'].setVisible(False)

        self.basic_widgets['type_label'] = QtWidgets.QLabel('Choose table type:')

        self.basic_grid.addWidget(self.basic_widgets['pipeline_name'], 1, 1)
        self.basic_grid.addWidget(self.basic_widgets['name_label'], 0, 0, 1, 2)
        self.basic_grid.addWidget(self.basic_widgets['type_label'], 0, 2)
        self.basic_grid.addWidget(self.basic_widgets['table_type_combo'], 1, 2)
        self.basic_grid.addWidget(self.basic_widgets['start_button'], 1, 3)

        self.basic_grid.addWidget(QtWidgets.QWidget(), 2, 0)
        self.basic_grid.addWidget(QtWidgets.QWidget(), 0, 4)
        self.basic_grid.setRowStretch(2, 1)
        self.basic_grid.setColumnStretch(4, 1)

    @QtCore.pyqtSlot(int)
    def _update_stack_status(self, ind: int):
        self.stack.widget(ind).check_selection_status()

    def init_function_ui(self):
        super().init_function_ui()
        self.function_group.setTitle("Add functions to Pipeline")
        sorted_actions = self.get_all_actions()
        for action_type in sorted_actions:
            self.stack_widgets[action_type].pipeline_mode = True
            self.stack_widgets[action_type].excluded_params.add('inplace')

    def apply_function(self):
        this_stack: FuncTypeStack = self.stack.currentWidget()
        func_name = this_stack.get_function_name()
        func_params = this_stack.get_function_params()
        self.pipeline.add_function(func_name, **func_params)
        self.update_pipeline_preview()
        self.is_unsaved = True

    def _apply_function_from_params(self, func_name, args: list, kwargs: dict):
        raise NotImplementedError

    def update_pipeline_preview(self):
        self.overview_widgets['preview'].setPlainText(str(self.pipeline))

    def update_table_preview(self):
        raise NotImplementedError

    def update_filter_obj_shape(self):
        raise NotImplementedError

    def set_file_path_bg_color(self):
        raise NotImplementedError

    def start(self):
        filt_obj_type = FILTER_OBJ_TYPES[self.basic_widgets['table_type_combo'].currentText()]
        self.filter_obj = filt_obj_type.__new__(filt_obj_type)
        self.pipeline = filtering.Pipeline(filt_obj_type)
        self.init_overview_ui()
        self.init_function_ui()
        self.basic_group.setVisible(False)
        self.is_unsaved = True

    def init_overview_ui(self):
        self.function_group.setTitle("Pipeline preview")
        self.layout.insertWidget(1, self.overview_group)
        self.overview_widgets['preview'] = QtWidgets.QPlainTextEdit()
        self.overview_widgets['preview'].setReadOnly(True)
        self.update_pipeline_preview()

        self.overview_grid.addWidget(QtWidgets.QLabel(f"Pipeline name: "
                                                      f"'<b>{self._get_pipeline_name()}</b>'"), 0, 0, 1, 2)
        self.overview_grid.addWidget(self.overview_widgets['preview'], 2, 0, 2, 4)

        self.overview_grid.addWidget(QtWidgets.QLabel(f"Pipeline table type: "
                                                      f"{self.basic_widgets['table_type_combo'].currentText()}"),
                                     5, 0, 1, 2)

        self.overview_widgets['save_button'] = QtWidgets.QPushButton('Save Pipeline')
        self.overview_widgets['save_button'].clicked.connect(self.save_file)
        self.overview_grid.addWidget(self.overview_widgets['save_button'], 5, 3)

        self.overview_widgets['export_button'] = QtWidgets.QPushButton('Export Pipeline')
        self.overview_widgets['export_button'].clicked.connect(self.export_pipeline)
        self.overview_grid.addWidget(self.overview_widgets['export_button'], 5, 2, 1, 1)

    def _get_pipeline_name(self):
        return self.basic_widgets['pipeline_name'].text()

    def export_pipeline(self):
        self.pipelineExported.emit(self._get_pipeline_name(), self.pipeline)

    def save_file(self):
        try:
            self.pipelineSaved.emit(self._get_pipeline_name(), self.pipeline)
            print(f"Successfully saved Pipeline '{self.basic_widgets['pipeline_name'].text()}'")
            self.is_unsaved = False
        except Exception as e:
            print("Failed to save Pipeline")
            raise e

    def closeEvent(self, event):
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
    def __init__(self, objs: List[filtering.Filter], parent=None):
        super().__init__(parent)
        self.objs = {str(obj.fname.stem): obj for obj in objs}
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

    def result(self):
        keep_tables = {file: widget.isChecked() for file, widget in self.keep_marks.items()}
        new_names = {file: gui_widgets.get_val_from_widget(widget) for file, widget in self.names.items()}
        out = []
        for file in keep_tables:
            if keep_tables[file]:
                obj = self.objs[file]
                if new_names[file] != '':
                    obj.fname = Path(new_names[file])
                out.append(obj)
        return out


class MultiOpenWindow(QtWidgets.QDialog):
    def __init__(self, files: List[str], parent=None):
        super().__init__(parent)
        self.files = files
        self.button_box = QtWidgets.QDialogButtonBox(QtWidgets.QDialogButtonBox.Ok | QtWidgets.QDialogButtonBox.Cancel)
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
            QtWidgets.QLabel(
                'Please choose a table type (mandatory) and table name (optional) for each loaded table\n\n'),
            0, 0, 1, 3)
        self.scroll_layout.addWidget(QtWidgets.QLabel('Table paths:'), 1, 0)
        self.scroll_layout.addWidget(QtWidgets.QLabel('Table types:'), 1, 1)
        self.scroll_layout.addWidget(QtWidgets.QLabel('Table names (optional):'), 1, 2)
        self.scroll_layout.addWidget(QtWidgets.QLabel('Additional parameters:'), 1, 3)
        for i, file in enumerate(self.files):
            self.paths[file] = gui_widgets.PathLineEdit(file, parent=self)
            self.table_types[file] = QtWidgets.QComboBox(self)
            self.table_types[file].addItems(list(FILTER_OBJ_TYPES.keys()))
            self.names[file] = QtWidgets.QLineEdit(self)
            kwargs_widget = QtWidgets.QWidget(self)
            self.kwargs[file] = QtWidgets.QGridLayout(kwargs_widget)

            self.table_types[file].currentTextChanged.connect(functools.partial(self.update_args_ui, file))
            self.table_types[file].setCurrentText('Other')

            self.scroll_layout.addWidget(self.paths[file], i + 2, 0)
            self.scroll_layout.addWidget(self.table_types[file], i + 2, 1)
            self.scroll_layout.addWidget(self.names[file], i + 2, 2)
            self.scroll_layout.addWidget(kwargs_widget, i + 2, 3)
        self.main_layout.addWidget(self.button_box)

        self.scroll.setMinimumWidth(self.scroll_widget.sizeHint().width() + 150)

    @QtCore.pyqtSlot(str)
    def update_args_ui(self, file: str):
        # clear previous layout
        gui_widgets.clear_layout(self.kwargs[file])
        self.kwargs_widgets[file] = {}
        func_name = '__init__'
        filter_obj_type = FILTER_OBJ_TYPES[self.table_types[file].currentText()]
        signature = generic.get_method_signature(func_name, filter_obj_type)
        desc, param_desc = generic.get_method_docstring(func_name, filter_obj_type)
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
                help_button.connect_param_help(name, param_desc[name])

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
    newTabFromSet = QtCore.pyqtSignal(set, str)
    newTabFromFilter = QtCore.pyqtSignal(filtering.Filter, str)

    def __init__(self, parent=None):
        super().__init__(parent)
        self.tabBar().setDocumentMode(True)
        self.tabBar().setMovable(True)
        self.setTabsClosable(True)

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

    def new_tab_from_item(self, item: Union[filtering.Filter, set], name: str):
        if isinstance(item, set):
            self.newTabFromSet.emit(item, name)
        elif validation.isinstanceinh(item, filtering.Filter):
            self.newTabFromFilter.emit(item, name)
        else:
            raise TypeError(type(item))

    def removeTab(self, index):
        h = self.cornerWidget().height()

        super().removeTab(index)
        self.update()
        if self.count() == 0:
            self.cornerWidget().setMinimumHeight(h)
            self.setMinimumHeight(h)


class RenameCommand(QtWidgets.QUndoCommand):
    def __init__(self, prev_name: str, new_name: str, tab: TabPage, description: str):
        super().__init__(description)
        self.prev_name = prev_name
        self.new_name = new_name
        self.tab = tab

    def undo(self):
        self.tab._rename(self.prev_name)

    def redo(self):
        self.tab._rename(self.new_name)


class CloseTabCommand(QtWidgets.QUndoCommand):
    def __init__(self, tab_container: ReactiveTabWidget, tab_index: int, description: str):
        super().__init__(description)
        self.tab_container = tab_container
        self.tab_index = tab_index
        self.tab_icon = tab_container.tabIcon(self.tab_index)
        self.tab_name = tab_container.tabText(self.tab_index).rstrip('*')
        self.obj_type = self.tab_container.widget(self.tab_index).obj_type()
        self.filename = self.tab_container.widget(self.tab_index).cache()

        self.kwargs = self.tab_container.widget(self.tab_index).obj_properties()

    def undo(self):
        item = io.load_cached_gui_file(self.filename)
        if isinstance(item, pd.DataFrame):
            item = self.obj_type.from_dataframe(item, self.tab_name, **self.kwargs)

        self.tab_container.new_tab_from_item(item, self.tab_name)
        self.tab_container.setTabIcon(self.tab_container.currentIndex(), self.tab_icon)

    def redo(self):
        self.tab_container.removeTab(self.tab_index)


class InplaceCommand(QtWidgets.QUndoCommand):
    def __init__(self, tab_page: FilterTabPage, func_name: str, args: list, kwargs: dict, description: str):
        super().__init__(description)
        self.tab_page = tab_page
        self.func_name = func_name
        self.args = args
        self.kwargs = kwargs
        self.obj_copy = copy.copy(self.tab_page.obj())

    def undo(self):
        obj = self.tab_page.obj()
        del obj
        self.tab_page.update_obj(copy.copy(self.obj_copy))
        self.tab_page.update_tab(is_unsaved=True)

    def redo(self):
        self.tab_page._apply_function_from_params(self.func_name, self.args, self.kwargs)


class InplaceCachedCommand(InplaceCommand):
    def __init__(self, tab_page: FilterTabPage, func_name: str, args: list, kwargs: dict, description: str):
        super().__init__(tab_page, func_name, args, kwargs, description)
        self.first_pass = True
        self.processed_obj = None

    def redo(self):
        if self.first_pass:
            self.first_pass = False
            super().redo()
        else:
            self.tab_page.update_obj(copy.copy(self.processed_obj))
            self.tab_page.update_tab(is_unsaved=True)

    def undo(self):
        processed_obj = self.tab_page.obj()
        self.processed_obj = copy.copy(processed_obj)
        self.tab_page.update_obj(copy.copy(self.obj_copy))
        self.tab_page.update_tab(is_unsaved=True)


class SetOpInplacCommand(InplaceCommand):
    def redo(self):
        is_filter_obj = validation.isinstanceinh(self.tab_page.obj(), filtering.Filter)
        first_obj = self.tab_page.obj()
        if not is_filter_obj:
            first_obj = filtering.Filter.from_dataframe(pd.DataFrame(index=first_obj), 'placeholder')
        getattr(first_obj, self.func_name)(*self.args, **self.kwargs)
        if not is_filter_obj:
            self.tab_page.update_obj(first_obj.index_set)
        self.tab_page.update_tab()


class PipelineInplaceCommand(QtWidgets.QUndoCommand):
    def __init__(self, tab_page: FilterTabPage, pipeline: filtering.Pipeline, description: str):
        super().__init__(description)
        self.tab_page = tab_page
        self.pipeline = pipeline
        self.obj_copy = copy.copy(self.tab_page.filter_obj)

    def undo(self):
        del self.tab_page.filter_obj
        self.tab_page.filter_obj = copy.copy(self.obj_copy)
        self.tab_page.update_tab(is_unsaved=True)

    def redo(self):
        self.tab_page._apply_pipeline(self.pipeline, inplace=True)


class ApplyPipelineWindow(gui_widgets.MinMaxDialog):
    def __init__(self, available_objects: dict, parent=None):
        super().__init__(parent)
        self.available_objects = available_objects
        self.layout = QtWidgets.QVBoxLayout(self)
        self.label = QtWidgets.QLabel('Choose the tables you wish to apply your Pipeline to', self)
        self.button_box = QtWidgets.QDialogButtonBox(QtWidgets.QDialogButtonBox.Ok | QtWidgets.QDialogButtonBox.Cancel)
        self.list = gui_widgets.MultipleChoiceList(self.available_objects,
                                                   [val[1] for val in self.available_objects.values()],
                                                   self)

        self.init_ui()

    def init_ui(self):
        self.button_box.accepted.connect(self.accept)
        self.button_box.rejected.connect(self.reject)
        self.layout.addWidget(self.label)
        self.layout.addWidget(self.list)
        self.layout.addWidget(self.button_box, QtCore.Qt.AlignCenter)

    def result(self):
        return [item.text() for item in self.list.get_sorted_selection()]


class MainWindow(QtWidgets.QMainWindow):
    USER_GUIDE_URL = 'https://guyteichman.github.io/RNAlysis/build/user_guide.html'
    jobQueued = QtCore.pyqtSignal()

    def __init__(self):
        super().__init__()
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
        self.error_window = None

        self.menu_bar = QtWidgets.QMenuBar(self)
        self.tab_contextmenu = None
        self.pipelines: typing.OrderedDict[str, filtering.Pipeline] = OrderedDict()
        self.pipeline_window = None

        self.about_window = gui_windows.AboutWindow(self)
        self.settings_window = gui_windows.SettingsWindow(self)
        self.settings_window.styleSheetUpdated.connect(self.update_style_sheet)
        self.set_op_window = None
        self.set_visualization_window = None
        self.tutorial_window = gui_tutorial.WelcomeWizard(self)
        self.cite_window = gui_windows.HowToCiteWindow(self)
        self.enrichment_window = None
        self.enrichment_results = []
        self.cutadapt_window = None
        self.kallisto_window = None

        self.init_ui()
        self.add_new_tab()
        self.init_actions()
        self.init_menu_ui()
        self.init_shortcuts()

        if settings.get_show_tutorial_settings():
            self.tutorial_window.show()

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

        # init progress bar and thread execution attributes
        self.progress_bar = None
        self.progbar_desc = ''
        self.progbar_total = 0
        self.progbar_start_time = 0
        self.progbar_completed_items = 0
        self.worker = None
        self.thread = QtCore.QThread()
        self.job_queue = Queue()
        self.job_timer = QtCore.QTimer(self)
        self.job_timer.timeout.connect(self.run_partial)
        self.job_timer.start(3500)
        self.jobQueued.connect(self.run_partial)

    def init_ui(self):
        self.setWindowTitle(f'RNAlysis {__version__}')
        icon_pth = str(Path(__file__).parent.parent.joinpath('favicon.ico').absolute())
        self.setWindowIcon(QtGui.QIcon(icon_pth))
        self.setGeometry(600, 50, 1050, 800)
        self.update_style_sheet()

        self.tabs.tabRightClicked.connect(self.init_tab_contextmenu)
        self.tabs.tabCloseRequested.connect(self.close_tab)
        self.tabs.newTabFromSet.connect(self.new_tab_from_gene_set)
        self.tabs.newTabFromFilter.connect(self.new_tab_from_filter_obj)

        self.command_history_dock.setWidget(self.undo_view)
        self.command_history_dock.setFloating(False)
        self.addDockWidget(QtCore.Qt.RightDockWidgetArea, self.command_history_dock)

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

        sort_menu = self.tab_contextmenu.addMenu("Sort tabs")
        sort_by_name = QtWidgets.QAction("Sort by tab &name")
        sort_by_name.triggered.connect(self.sort_tabs_by_name)
        sort_by_time = QtWidgets.QAction("Sort by creation &time")
        sort_by_time.triggered.connect(self.sort_tabs_by_creation_time)
        sort_by_type = QtWidgets.QAction("Sort by tab type")
        sort_by_type.triggered.connect(self.sort_tabs_by_type)
        sort_by_size = QtWidgets.QAction("Sort by number of features")
        sort_by_size.triggered.connect(self.sort_tabs_by_n_features)
        reverse = QtWidgets.QAction("Reverse order")
        reverse.triggered.connect(self.sort_reverse)
        sort_menu.addActions([sort_by_name, sort_by_time, sort_by_type, sort_by_size, reverse])

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

    def add_new_tab(self, name: str = None, is_set: bool = False):
        new_undo_stack = QtWidgets.QUndoStack()
        self.undo_group.addStack(new_undo_stack)
        if name is None:
            name = 'New Table'

        else:
            print(name)

        if is_set:
            tab = SetTabPage(name, parent=self.tabs, undo_stack=new_undo_stack)
        else:
            tab = FilterTabPage(self.tabs, undo_stack=new_undo_stack)
            tab.filterObjectCreated.connect(self.new_tab_from_filter_obj)
            tab.startedClustering.connect(self.start_clustering)
            tab.startedJob.connect(self.start_generic_job)
        tab.changeIcon.connect(self.set_current_tab_icon)
        tab.tabNameChange.connect(self.rename_tab)
        self.tabs.addTab(tab, name)
        self.tabs.setCurrentIndex(self.tabs.count() - 1)

    @QtCore.pyqtSlot(str, bool)
    def rename_tab(self, new_name: str, is_unsaved: bool):
        if is_unsaved:
            new_name += '*'
        else:
            new_name = new_name.rstrip('*')
        self.tabs.setTabText(self.tabs.currentIndex(), new_name)

    @QtCore.pyqtSlot()
    def remove_tab_asterisk(self):
        current_name = self.tabs.tabText(self.tabs.currentIndex())
        self.tabs.setTabText(self.tab.currentIndex(), current_name.rstrip('*'))

    def new_table_from_folder(self):
        folder_name = QtWidgets.QFileDialog.getExistingDirectory(self, "Choose directory", str(Path.home()))
        if folder_name:
            normalize_answer = QtWidgets.QMessageBox.question(self, 'Normalize values?',
                                                              "Do you want to normalize your count table to "
                                                              "reads-per-million (RPM)?",
                                                              QtWidgets.QMessageBox.Yes | QtWidgets.QMessageBox.No)
            to_normalize = normalize_answer == QtWidgets.QMessageBox.Yes

            filter_obj = filtering.CountFilter.from_folder(folder_name, norm_to_rpm=to_normalize)
            if self.tabs.currentWidget().is_empty():
                self.tabs.removeTab(self.tabs.currentIndex())
            self.new_tab_from_filter_obj(filter_obj)

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
                            self.new_tab_from_filter_obj(filter_obj)
                        else:
                            self.new_tab_from_filter_obj(filter_obj, name)
                    if tabs_to_close is not None:
                        self.tabs.removeTab(tabs_to_close)

    def new_tab_from_filter_obj(self, filter_obj: filtering.Filter, name: str = None):
        self.add_new_tab(filter_obj.fname.name)
        self.tabs.currentWidget().start_from_filter_obj(filter_obj, name)
        QtWidgets.QApplication.processEvents()

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

    def new_tab_from_gene_set(self, gene_set: set, gene_set_name: str):
        self.add_new_tab(gene_set_name, is_set=True)
        self.tabs.currentWidget().update_gene_set(gene_set)
        QtWidgets.QApplication.processEvents()

    def export_pipeline(self):
        if len(self.pipelines) == 0:
            warnings.warn('No Pipelines to export')
            return

        pipeline_name, status = QtWidgets.QInputDialog.getItem(
            self, 'Export Pipeline', 'Choose Pipeline to export:', self.pipelines.keys())
        if status:
            pipeline = self.pipelines[pipeline_name]
            self._export_pipeline_from_obj(pipeline_name, pipeline)

    def _export_pipeline_from_obj(self, pipeline_name: str, pipeline: filtering.Pipeline):
        default_name = pipeline_name + '.yaml'
        filename, _ = QtWidgets.QFileDialog.getSaveFileName(self, "Save Pipeline",
                                                            str(Path.home().joinpath(default_name)),
                                                            "YAML file (*.yaml)")
        if filename:
            pipeline.export_pipeline(filename)
            print(f"Successfully saved at {io.get_datetime()} under {filename}")

    def import_pipeline(self):
        filename, _ = QtWidgets.QFileDialog.getOpenFileName(self, "Choose a Pipeline file", str(Path.home()),
                                                            "YAML file (*.yaml)")
        if filename:
            pipeline = filtering.Pipeline.import_pipeline(filename)
            self.pipelines[str(Path(filename).stem)] = pipeline

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
                self.new_tab_from_gene_set(gene_set, Path(filename).stem)
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
            self.new_tab_from_gene_set(gene_set, Path(filename).stem)
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

    def close_tab(self, index):
        if self.tabs.widget(index).is_empty():
            self.tabs.removeTab(index)
        else:
            self.undo_group.removeStack(self.tabs.widget(index).undo_stack)
            command = CloseTabCommand(self.tabs, index, '"' + self.tabs.tabText(index).rstrip('*') + '"')
            self.closed_tabs_stack.push(command)

    def copy_gene_set(self):
        gene_set = self.tabs.currentWidget().get_index_string()
        QtWidgets.QApplication.clipboard().setText(gene_set)

    def excepthook(self, exc_type, exc_value, exc_tb):
        sys.__excepthook__(exc_type, exc_value, exc_tb)
        self.error_window = gui_windows.ErrorMessage(exc_type, exc_value, exc_tb, self)
        self.error_window.exec()

    def _get_current_console(self):
        if self.pipeline_window is not None and self.pipeline_window.isVisible():
            current_console = self.pipeline_window.get_console()
        else:
            current_console = self.tabs.currentWidget().get_console()
        return current_console

    @QtCore.pyqtSlot(str)
    def append_text_to_current_console(self, text: str):
        current_console = self._get_current_console()
        current_console.append_text(text)

    def add_pipeline(self):
        self.pipeline_window = CreatePipelineWindow(self)
        self.pipeline_window.pipelineSaved.connect(self.save_pipeline)
        self.pipeline_window.pipelineExported.connect(self._export_pipeline_from_obj)
        self.pipeline_window.exec()
        self.pipeline_window = None

    @QtCore.pyqtSlot(str, filtering.Pipeline)
    def save_pipeline(self, pipeline_name: str, pipeline: filtering.Pipeline):
        self.pipelines[pipeline_name] = pipeline

    def settings(self):
        self.settings_window.exec()

    def init_actions(self):
        self.new_table_action = QtWidgets.QAction("&New table", self)
        self.new_table_action.triggered.connect(functools.partial(self.add_new_tab, name=None))
        self.new_table_from_folder_action = QtWidgets.QAction("New table from &folder")
        self.new_table_from_folder_action.triggered.connect(self.new_table_from_folder)

        self.new_multiple_action = QtWidgets.QAction("&Multiple new tables")
        self.new_multiple_action.triggered.connect(self.load_multiple_files)

        self.save_action = QtWidgets.QAction("&Save...", self)
        self.save_action.triggered.connect(self.save_file)

        self.load_session_action = QtWidgets.QAction("&Load session...")
        self.load_session_action.triggered.connect(self.load_session)
        self.save_session_action = QtWidgets.QAction("Sa&ve session...")
        self.save_session_action.triggered.connect(self.save_session)
        self.clear_session_action = QtWidgets.QAction("Clea&r session...")
        self.clear_session_action.triggered.connect(self.clear_session)

        self.settings_action = QtWidgets.QAction("&Settings...", self)
        self.settings_action.triggered.connect(self.settings)
        self.exit_action = QtWidgets.QAction("&Exit", self)
        self.exit_action.triggered.connect(self.close)
        self.check_update_action = QtWidgets.QAction("Check for &updates...", self)
        self.check_update_action.triggered.connect(self.check_for_updates)

        self.undo_action = self.undo_group.createUndoAction(self)
        self.redo_action = self.undo_group.createRedoAction(self)
        self.restore_tab_action = self.closed_tabs_stack.createUndoAction(self, 'Restore tab')
        self.close_current_action = QtWidgets.QAction("&Close current tab", self)
        self.close_current_action.triggered.connect(self.close_current_tab)

        self.show_history_action = QtWidgets.QAction("Command &History")
        self.show_history_action.setCheckable(True)
        self.show_history_action.setChecked(True)
        self.show_history_action.triggered.connect(self.toggle_history_window)
        self.clear_history_action = QtWidgets.QAction("Clea&r command history")
        self.clear_history_action.triggered.connect(self.clear_history)

        self.copy_action = QtWidgets.QAction("&Copy Gene Set", self)
        self.copy_action.triggered.connect(self.copy_gene_set)
        self.set_op_action = QtWidgets.QAction("Set &Operations...", self)
        self.set_op_action.triggered.connect(self.choose_set_op)
        self.enrichment_action = QtWidgets.QAction("Enrichment &Analysis...", self)
        self.enrichment_action.triggered.connect(self.open_enrichment_analysis)
        self.set_vis_action = QtWidgets.QAction("&Visualize Gene Sets...", self)
        self.set_vis_action.triggered.connect(self.visualize_gene_sets)
        self.import_set_action = QtWidgets.QAction("&Import Gene Set...", self)
        self.import_set_action.triggered.connect(self.import_gene_set)
        self.import_multiple_sets_action = QtWidgets.QAction("Import &Multiple Gene Sets...", self)
        self.import_multiple_sets_action.triggered.connect(self.import_multiple_gene_sets)
        self.export_set_action = QtWidgets.QAction("&Export Gene Set...", self)
        self.export_set_action.triggered.connect(self.export_gene_set)

        self.cutadapt_single_action = QtWidgets.QAction("&Single-end adapter trimming...", self)
        self.cutadapt_single_action.triggered.connect(functools.partial(self.trim_adapters, True))
        self.cutadapt_paired_action = QtWidgets.QAction("&Paired-end adapter trimming...", self)
        self.cutadapt_paired_action.triggered.connect(functools.partial(self.trim_adapters, False))

        self.kallisto_index_action = QtWidgets.QAction("Create kallisto &index...", self)
        self.kallisto_index_action.triggered.connect(functools.partial(self.kallisto, 'index'))
        self.kallisto_single_action = QtWidgets.QAction("&Single-end RNA-seq quantification...", self)
        self.kallisto_single_action.triggered.connect(functools.partial(self.kallisto, 'single'))
        self.kallisto_paired_action = QtWidgets.QAction("&Paired-end RNA-seq quantification...", self)
        self.kallisto_paired_action.triggered.connect(functools.partial(self.kallisto, 'paired'))

        self.tutorial_action = QtWidgets.QAction("&Tutorial", self)
        self.tutorial_action.triggered.connect(self.tutorial_window.show)
        self.user_guide_action = QtWidgets.QAction("&User Guide", self)
        self.user_guide_action.triggered.connect(self.open_user_guide)
        self.about_action = QtWidgets.QAction("&About", self)
        self.about_action.triggered.connect(self.about)
        self.cite_action = QtWidgets.QAction("How to &cite RNAlysis", self)
        self.cite_action.triggered.connect(self.cite)

        self.new_pipeline_action = QtWidgets.QAction("&New Pipeline...", self)
        self.new_pipeline_action.triggered.connect(self.add_pipeline)
        self.import_pipeline_action = QtWidgets.QAction("&Import Pipeline...", self)
        self.import_pipeline_action.triggered.connect(self.import_pipeline)
        self.export_pipeline_action = QtWidgets.QAction("&Export Pipeline...", self)
        self.export_pipeline_action.triggered.connect(self.export_pipeline)

    def init_shortcuts(self):
        self.copy_action.setShortcut(QtGui.QKeySequence("Ctrl+C"))
        self.save_action.setShortcut(QtGui.QKeySequence("Ctrl+S"))
        self.new_multiple_action.setShortcut(QtGui.QKeySequence("Ctrl+Shift+N"))
        self.new_table_action.setShortcut(QtGui.QKeySequence("Ctrl+N"))

        self.undo_action.setShortcut(QtGui.QKeySequence("Ctrl+Z"))
        self.redo_action.setShortcut(QtGui.QKeySequence("Ctrl+Shift+Z"))
        self.restore_tab_action.setShortcut(QtGui.QKeySequence("Ctrl+Shift+T"))
        self.close_current_action.setShortcut(QtGui.QKeySequence("Ctrl+W"))

        self.import_set_action.setShortcut(QtGui.QKeySequence("Ctrl+Shift+I"))
        self.import_multiple_sets_action.setShortcut(QtGui.QKeySequence("Ctrl+Shift+M"))
        self.export_set_action.setShortcut(QtGui.QKeySequence("Ctrl+Shift+E"))
        self.set_vis_action.setShortcut(QtGui.QKeySequence("Ctrl+Shift+V"))
        self.enrichment_action.setShortcut(QtGui.QKeySequence("Ctrl+Shift+A"))
        self.set_op_action.setShortcut(QtGui.QKeySequence("Ctrl+Shift+O"))

        self.new_pipeline_action.setShortcut(QtGui.QKeySequence("Ctrl+Alt+N"))
        self.import_pipeline_action.setShortcut(QtGui.QKeySequence("Ctrl+Alt+I"))
        self.export_pipeline_action.setShortcut(QtGui.QKeySequence("Ctrl+Alt+I"))

    @QtCore.pyqtSlot()
    def check_for_updates(self, confirm_updated: bool = True):
        if io.is_rnalysis_outdated():
            reply = QtWidgets.QMessageBox.question(self, 'A new version is available',
                                                   f'A new version of <i>RNAlysis</i> is available! '
                                                   f'Do you wish to update?')
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

    @QtCore.pyqtSlot(bool)
    def trim_adapters(self, single_end: bool):
        if single_end:
            self.cutadapt_window = CutAdaptSingleWindow(self)
        else:
            self.cutadapt_window = CutAdaptPairedWindow(self)
        func_name = self.cutadapt_window.func_name
        func = self.cutadapt_window.func
        self.cutadapt_window.paramsAccepted.connect(
            functools.partial(self.start_generic_job_from_params, func_name, func))
        self.cutadapt_window.show()

    @QtCore.pyqtSlot(str)
    def kallisto(self, mode: typing_extensions.Literal['index', 'single', 'paired']):
        if mode == 'single':
            self.kallisto_window = KallistoSingleWindow(self)
        elif mode == 'paired':
            self.kallisto_window = KallistoPairedWindow(self)
        elif mode == 'index':
            self.kallisto_window = KallistoIndexWindow(self)
        else:
            raise ValueError(f"invalid mode {mode}")
        func_name = self.kallisto_window.func_name
        func = self.kallisto_window.func
        self.kallisto_window.paramsAccepted.connect(
            functools.partial(self.start_generic_job_from_params, func_name, func))
        self.kallisto_window.show()

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
        default_name = self.tabs.tabText(self.tabs.currentIndex()).rstrip("*") + '.txt'
        filename, _ = QtWidgets.QFileDialog.getSaveFileName(self, "Save gene set",
                                                            str(Path.home().joinpath(default_name)),
                                                            "Text document (*.txt);;"
                                                            "All Files (*)")
        if filename:
            io.save_gene_set(gene_set, filename)
            print(f"Successfully saved at {io.get_datetime()} under {filename}")

    def open_user_guide(self):
        url = QtCore.QUrl(self.USER_GUIDE_URL)
        if not QtGui.QDesktopServices.openUrl(url):
            QtGui.QMessageBox.warning(self, 'User Guide', 'Could not open User Guide')

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
        self.set_op_window.geneSetReturned.connect(self.new_tab_from_gene_set)
        self.set_op_window.show()

    def visualize_gene_sets(self):
        available_objs = self.get_available_objects()
        self.set_visualization_window = SetVisualizationWindow(available_objs, self)
        self.set_visualization_window.show()

    @QtCore.pyqtSlot(str)
    def choose_tab_by_name(self, set_name: str):
        available_objs = self.get_available_objects()
        for i, name in enumerate(available_objs.keys()):
            if name == set_name:
                self.tabs.setCurrentIndex(i)
                return

    def display_enrichment_results(self, result: pd.DataFrame, gene_set_name: str):
        df_window = gui_windows.DataFrameView(result, "Enrichment results for set " + gene_set_name)
        self.enrichment_results.append(df_window)
        df_window.show()

    def open_enrichment_analysis(self):
        self.enrichment_window = EnrichmentWindow(self.get_available_objects(), self)
        self.enrichment_window.enrichmentStarted.connect(self.start_enrichment)
        self.enrichment_window.show()

    def get_tab_names(self) -> List[str]:
        return [self.tabs.tabText(i).rstrip('*') for i in range(self.tabs.count())]

    def init_menu_ui(self):
        self.setMenuBar(self.menu_bar)
        file_menu = self.menu_bar.addMenu("&File")
        self.new_menu = file_menu.addMenu("&New...")
        self.new_menu.addActions([self.new_table_action, self.new_multiple_action, self.new_table_from_folder_action])
        file_menu.addActions(
            [self.save_action, self.load_session_action, self.save_session_action, self.clear_session_action,
             self.settings_action, self.exit_action])

        edit_menu = self.menu_bar.addMenu("&Edit")
        edit_menu.addActions([self.undo_action, self.redo_action, self.restore_tab_action, self.close_current_action,
                              self.clear_history_action])

        view_menu = self.menu_bar.addMenu("&View")
        view_menu.addActions([self.show_history_action])

        fastq_menu = self.menu_bar.addMenu("&FASTQ")
        self.trimming_menu = fastq_menu.addMenu('&Adapter trimming')
        self.kallisto_menu = fastq_menu.addMenu("RNA-sequencing &quantification")
        self.trimming_menu.addActions([self.cutadapt_single_action, self.cutadapt_paired_action])
        self.kallisto_menu.addActions(
            [self.kallisto_index_action, self.kallisto_single_action, self.kallisto_paired_action])

        gene_sets_menu = self.menu_bar.addMenu("&Gene sets")
        gene_sets_menu.addActions(
            [self.copy_action, self.set_op_action, self.enrichment_action, self.set_vis_action, self.import_set_action,
             self.import_multiple_sets_action, self.export_set_action])

        pipeline_menu = self.menu_bar.addMenu("&Pipelines")
        pipeline_menu.addActions([self.new_pipeline_action, self.import_pipeline_action, self.export_pipeline_action])
        self.apply_pipeline_menu = pipeline_menu.addMenu("&Apply Pipeline")
        self.apply_pipeline_menu.aboutToShow.connect(self._populate_pipelines)

        help_menu = self.menu_bar.addMenu("&Help")
        help_menu.addActions([self.tutorial_action, self.user_guide_action, self.check_update_action, self.about_action,
                              self.cite_action])

    def _populate_pipelines(self):
        # Remove the old options from the menu
        self.apply_pipeline_menu.clear()
        # Dynamically create the actions
        actions = []
        for name, pipeline in self.pipelines.items():
            action = QtWidgets.QAction(name, self)
            action.triggered.connect(functools.partial(self.apply_pipeline, pipeline, name))
            actions.append(action)
        # Step 3. Add the actions to the menu
        self.apply_pipeline_menu.addActions(actions)

    def apply_pipeline(self, pipeline: filtering.Pipeline, pipeline_name: str):
        apply_msg = f"Do you want to apply Pipeline '{pipeline_name}' inplace?"
        reply = QtWidgets.QMessageBox.question(self, f"Apply Pipeline '{pipeline_name}'",
                                               apply_msg, QtWidgets.QMessageBox.Yes, QtWidgets.QMessageBox.No)

        inplace = reply == QtWidgets.QMessageBox.Yes

        available_objs = self.get_available_objects()
        filtered_available_objs = {}
        for i, (key, val) in enumerate(available_objs.items()):
            if (self.tabs.widget(i).obj_type() == pipeline.filter_type) or (
                pipeline.filter_type == filtering.Filter and issubclass(self.tabs.widget(i).obj_type(),
                                                                        filtering.Filter)):
                filtered_available_objs[key] = val
        window = ApplyPipelineWindow(filtered_available_objs, self)
        accepted = window.exec()
        if accepted:
            current_ind = self.tabs.currentIndex()
            chosen_names = window.result()
            for name in chosen_names:
                self.tabs.setCurrentWidget(filtered_available_objs[name][0])
                filtered_available_objs[name][0].apply_pipeline(pipeline, name, inplace)
                QtWidgets.QApplication.processEvents()
            self.tabs.setCurrentIndex(current_ind)

    @QtCore.pyqtSlot()
    def clear_session(self, confirm_action: bool = True):
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

            self.pipelines = {}
            self.clear_history(confirm_action=False)

    def load_session(self):
        session_filename, _ = QtWidgets.QFileDialog.getOpenFileName(self, "Load session",
                                                                    str(Path.home()),
                                                                    "RNAlysis session files (*.rnal);;"
                                                                    "All Files (*)")
        if session_filename:
            items, item_names, item_types, item_properties, pipeline_names, pipeline_files = io.load_gui_session(
                session_filename)

            for pipeline_name, pipeline_file in zip(pipeline_names, pipeline_files):
                pipeline = filtering.Pipeline.import_pipeline(pipeline_file)
                self.pipelines[pipeline_name] = pipeline
            QtWidgets.QApplication.processEvents()

            tabs_to_close = None
            if self.tabs.currentWidget().is_empty():
                tabs_to_close = self.tabs.currentIndex()

            for item, item_name, item_type, item_property in zip(items, item_names, item_types, item_properties):
                if item_type == 'set':
                    self.new_tab_from_gene_set(item, item_name)
                else:
                    try:
                        cls = getattr(filtering, item_type)
                    except AttributeError:
                        raise TypeError(f"Invalid object type in session file: '{item_type}'")
                    obj = cls.from_dataframe(item, item_name, **item_property)
                    self.new_tab_from_filter_obj(obj, item_name)
            if tabs_to_close is not None:
                self.tabs.removeTab(tabs_to_close)

    def save_session(self):
        default_name = 'Untitled session.rnal'
        session_filename, _ = QtWidgets.QFileDialog.getSaveFileName(self, "Save session",
                                                                    str(Path.home().joinpath(default_name)),
                                                                    "RNAlysis session files (*.rnal);;"
                                                                    "All Files (*)")
        if session_filename:
            filenames = []
            item_names = []
            item_types = []
            item_properties = []
            for ind in range(self.tabs.count()):
                tab = self.tabs.widget(ind)
                if not tab.is_empty():
                    filenames.append(tab.cache())
                    item_names.append(self.tabs.tabText(ind).rstrip('*'))
                    item_types.append(tab.obj_type())
                    item_properties.append(tab.obj_properties())

            pipeline_names = []
            pipeline_files = []
            for pipeline_name, pipeline in self.pipelines.items():
                pipeline_files.append(pipeline.export_pipeline(filename=None))
                pipeline_names.append(pipeline_name)

            io.save_gui_session(session_filename, filenames, item_names, item_types, item_properties, pipeline_names,
                                pipeline_files)
            print(f"Session saved successfully at {io.get_datetime()} under {session_filename}")

    def about(self):
        self.about_window.exec()

    def cite(self):
        self.cite_window.exec()

    def input(self, message: str):
        dialog = gui_widgets.PathInputDialog(message, parent=self)
        accepted = dialog.exec()
        if accepted:
            return dialog.result()
        return None

    def closeEvent(self, event):

        quit_msg = "Are you sure you want to close <i>RNAlysis</i>?\n" \
                   "All unsaved progress will be lost"

        reply = QtWidgets.QMessageBox.question(self, 'Close program',
                                               quit_msg, QtWidgets.QMessageBox.No, QtWidgets.QMessageBox.Yes)

        if reply == QtWidgets.QMessageBox.Yes:
            self.thread_stdout_queue_listener.exit()

            io.clear_gui_cache()
            get_reusable_executor().shutdown(wait=True)

            event.accept()
        else:
            event.ignore()

    @QtCore.pyqtSlot(object, str, object)
    def start_generic_job(self, partial: Callable, func_name: str, finish_slot: Union[Callable, None]):
        slots = (self.finish_generic_job, finish_slot)
        self.queue_partial(partial, slots, func_name)

    def start_generic_job_from_params(self, func_name, func, args, kwargs, finish_slot: Union[Callable, None]):
        partial = functools.partial(func, *args, **kwargs)
        self.start_generic_job(partial, func_name, finish_slot)

    @QtCore.pyqtSlot(tuple)
    def finish_generic_job(self, worker_output):
        if len(worker_output) == 1:
            if isinstance(worker_output[0], Exception):
                raise worker_output[0]
            print("Done")
        else:
            assert len(worker_output) == 2
            func_name: str = worker_output[-1]
            return_val: tuple = worker_output[0]
            self.tabs.currentWidget().process_outputs(return_val, func_name)
            self.tabs.currentWidget().update_tab()

    @QtCore.pyqtSlot(object, str, object)
    def start_clustering(self, partial: Callable, func_name: str, finish_slot: Union[Callable, None]):
        slots = (self.finish_clustering, finish_slot)
        self.queue_partial(partial, slots, func_name)

    @QtCore.pyqtSlot(tuple)
    def finish_clustering(self, worker_output: tuple):
        if len(worker_output) == 1:
            if isinstance(worker_output[0], Exception):
                raise worker_output[0]
            return
        func_name: str = worker_output[-1]
        clustering_runner: clustering.ClusteringRunner = worker_output[0][1]
        return_val: tuple = worker_output[0][0]
        clustering_runner.plot_clustering()
        self.tabs.currentWidget().process_outputs(return_val, func_name)
        self.tabs.currentWidget().update_tab()

    @QtCore.pyqtSlot(object, str, object)
    def start_enrichment(self, partial: Callable, set_name: str, finish_slot: Union[Callable, None]):
        slots = (self.finish_enrichment, finish_slot)
        self.queue_partial(partial, slots, set_name)

    @QtCore.pyqtSlot(tuple)
    def finish_enrichment(self, worker_output: tuple):
        if len(worker_output) == 1:
            if isinstance(worker_output[0], Exception):
                raise worker_output[0]
            return
        set_name = worker_output[-1]
        results, enrichment_runner = worker_output[0]
        self.show()
        enrichment_runner.plot_results()
        self.display_enrichment_results(results, set_name)

    def queue_partial(self, partial: Callable, output_slots: Union[Callable, Tuple[Callable, ...], None] = None, *args):
        self.job_queue.put((partial, output_slots, args))
        self.jobQueued.emit()

    def run_partial(self):
        # if there are no jobs available, don't proceed
        if self.job_queue.qsize() == 0:
            return
        # if there is a job currently running, don't proceed
        try:
            if (self.thread is not None) and self.thread.isRunning():
                return
        except RuntimeError:
            self.thread = None
        partial, output_slots, args = self.job_queue.get()
        # Create a worker object
        self.worker = gui_widgets.Worker(partial, *args)

        def alt_tqdm(iter_obj: typing.Iterable = None, desc: str = '', unit: str = '', bar_format: str = '',
                     total: int = None):
            self.worker.startProgBar.emit(
                dict(iter_obj=iter_obj, desc=desc, unit=unit, bar_format=bar_format, total=total))
            obj = gui_widgets.AltTQDM(iter_obj, desc=desc, unit=unit, bar_format=bar_format, total=total)
            obj.barUpdate.connect(self.move_progress_bar)
            obj.barFinished.connect(self.end_progress_bar)
            return obj

        def alt_parallel(n_jobs: int = -1, desc: str = '', unit: str = '', bar_format: str = '',
                         total: int = None, **kwargs):
            self.worker.startProgBar.emit(
                dict(iter_obj=None, desc=desc, unit=unit, bar_format=bar_format, total=total))
            print(f"{desc}: started" + (f" {total} jobs\r" if isinstance(total, int) else "\r"))
            obj = gui_widgets.AltParallel(n_jobs=n_jobs, desc=desc, unit=unit, bar_format=bar_format,
                                          total=total, **kwargs)
            obj.barUpdate.connect(self.move_progress_bar)
            obj.barFinished.connect(self.end_progress_bar)
            obj.barTotalUpdate.connect(self.update_bar_total)

            return obj

        enrichment.enrichment_runner.generic.ProgressParallel = alt_parallel
        generic.ProgressParallel = alt_parallel
        filtering.clustering.generic.ProgressParallel = alt_parallel

        enrichment.enrichment_runner.parsing.tqdm = alt_tqdm
        enrichment.enrichment_runner.io.tqdm = alt_tqdm
        enrichment.enrichment_runner.tqdm = alt_tqdm
        filtering.clustering.tqdm = alt_tqdm
        fastq.tqdm = alt_tqdm

        #  Move worker to the thread
        self.worker.moveToThread(self.thread)
        #  Connect signals and slots
        self.worker.startProgBar.connect(self.start_progress_bar)
        self.thread.started.connect(self.worker.run)
        for slot in parsing.data_to_tuple(output_slots):
            if slot is not None:
                self.worker.finished.connect(slot)
        self.worker.finished.connect(self.run_partial)
        self.worker.finished.connect(self.thread.quit)
        self.worker.finished.connect(self.worker.deleteLater)

        self.thread.start()

    def start_progress_bar(self, arg_dict):
        total = arg_dict.get('total', None)
        iter_obj = arg_dict['iter_obj']
        self.progbar_desc = arg_dict.get('desc', '')
        if total is not None:
            self.progbar_total = total
        else:
            try:
                self.progbar_total = len(iter_obj)
            except TypeError:
                self.progbar_total = 2
        if self.progress_bar is not None:
            self.progress_bar.close()

        self.progress_bar = QtWidgets.QProgressDialog(self.progbar_desc, "Hide", 0, self.progbar_total, self)
        self.progbar_start_time = time.time()
        self.progress_bar.setMinimumDuration(0)
        self.progress_bar.setWindowTitle(self.progbar_desc)
        self.progress_bar.setValue(0)
        self.progress_bar.setWindowModality(QtCore.Qt.WindowModal)
        self.progbar_completed_items = 0

    def move_progress_bar(self, n: int):
        self.progbar_completed_items += n
        elapsed_time = time.time() - self.progbar_start_time
        remaining_time = (elapsed_time / self.progbar_completed_items) * abs(
            self.progbar_total - self.progbar_completed_items)
        try:
            self.progress_bar.setLabelText(self.progbar_desc + '\n' + f"Elapsed time: {elapsed_time:.2f} seconds"
                                           + '\n' + f"Remaining time: {remaining_time:.2f} seconds")
            self.progress_bar.setValue(self.progbar_completed_items)
        except AttributeError:
            pass

    def update_bar_total(self, n: int):
        self.progbar_total = n
        try:
            self.progress_bar.setMaximum(n)
        except AttributeError:
            pass

    def end_progress_bar(self):
        if self.progress_bar is not None:
            self.progress_bar.close()
            self.progress_bar = None


def customwarn(message, category, filename, lineno, file=None, line=None):
    sys.stdout.write(warnings.formatwarning(message, category, filename, lineno))


def run():
    app = QtWidgets.QApplication(sys.argv)
    splash = gui_windows.splash_screen()
    app.processEvents()
    base_message = f"<i>RNAlysis</i> version {__version__}:\t"
    splash.showMessage(base_message + 'loading dependencies', QtCore.Qt.AlignBottom | QtCore.Qt.AlignHCenter)
    matplotlib.use('Qt5Agg')

    video_files = gui_tutorial.WelcomeWizard.VIDEO_FILES
    for i in io.get_gui_videos(video_files):
        splash.showMessage(
            base_message + f'getting tutorial videos {i + 1}/{len(video_files)}',
            QtCore.Qt.AlignBottom | QtCore.Qt.AlignHCenter)

    splash.showMessage(base_message + 'loading application', QtCore.Qt.AlignBottom | QtCore.Qt.AlignHCenter)
    window = MainWindow()
    warnings.showwarning = customwarn
    sys.excepthook = window.excepthook
    builtins.input = window.input

    try:
        import numba
    except ImportError:
        warnings.warn("RNAlysis can perform faster when package 'numba' is installed. \n"
                      "If you want to improve the performance of slow operations on RNAlysis, "
                      "please install package 'numba'. ")

    window.show()
    window.check_for_updates(False)
    splash.finish(window)
    sys.exit(app.exec())
