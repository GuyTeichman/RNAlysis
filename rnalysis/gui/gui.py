import builtins
import functools
import itertools
import time
import os
import sys
import typing
import warnings
from collections import OrderedDict
from pathlib import Path
from queue import Queue
from typing import List, Union, Callable

import matplotlib
import numpy as np
import pandas as pd
from PyQt5 import QtCore, QtWidgets, QtGui

from rnalysis import filtering, enrichment, __version__
from rnalysis.gui import gui_style, gui_utils, gui_graphics
from rnalysis.utils import io, validation, generic, parsing

FILTER_OBJ_TYPES = {'Count matrix': filtering.CountFilter, 'Differential expression': filtering.DESeqFilter,
                    'Fold change': filtering.FoldChangeFilter, 'Other': filtering.Filter}
FILTER_OBJ_TYPES_INV = {val: key for key, val in FILTER_OBJ_TYPES.items()}


class ClicomWindow(gui_utils.MinMaxDialog):
    CLICOM_FUNC = filtering.CountFilter.split_clicom
    CLICOM_SIGNATURE = generic.get_method_signature(CLICOM_FUNC)
    CLICOM_DESC, CLICOM_PARAM_DESC = generic.get_method_docstring(CLICOM_FUNC)
    EXCLUDED_PARAMS = {'self', 'parameter_dicts'}
    ADDITIONAL_EXCLUDED_PARAMS = {'power_transform', 'plot_style', 'split_plots', 'return_probabilities'}
    paramsAccepted = QtCore.pyqtSignal(list, dict)

    def __init__(self, funcs: dict, filter_obj: filtering.Filter, parent=None):
        super().__init__(parent)
        self.parameter_dicts: List[dict] = []
        self.funcs = funcs
        self.setups_counter = {key: 0 for key in self.funcs.keys()}
        self.filter_obj = filter_obj
        self.stack = FuncTypeStack(self.funcs, self.filter_obj, self,
                                   additional_excluded_params=self.ADDITIONAL_EXCLUDED_PARAMS)
        self.widgets = {}
        self.layout = QtWidgets.QGridLayout(self)

        self.param_group = QtWidgets.QGroupBox('Choose CLICOM parameters')
        self.param_grid = QtWidgets.QGridLayout(self.param_group)
        self.param_widgets = {}

        self.setups_group = QtWidgets.QGroupBox("Choose clustering setups for CLICOM")
        self.setups_grid = QtWidgets.QGridLayout(self.setups_group)
        self.setups_widgets = {}

        self.start_button = QtWidgets.QPushButton('Accept CLICOM parameters')
        self.close_button = QtWidgets.QPushButton('Close')

        self.init_ui()

    def init_ui(self):
        self.setWindowTitle('CLICOM clustering setup')
        self.layout.addWidget(self.param_group, 0, 0)
        self.layout.addWidget(self.setups_group, 0, 1)
        self.layout.addWidget(self.start_button, 1, 0, 1, 2)
        self.layout.addWidget(self.close_button, 2, 0, 1, 2)

        self.start_button.clicked.connect(self.start_clustering)
        self.close_button.clicked.connect(self.close)

        self.init_param_ui()
        self.init_setups_ui()

    def init_param_ui(self):
        i = 0

        for name, param in self.CLICOM_SIGNATURE.items():
            if name in self.EXCLUDED_PARAMS:
                continue
            this_desc = self.CLICOM_PARAM_DESC.get(name, '')
            self.param_widgets[name] = gui_utils.param_to_widget(param, name)
            if isinstance(self.param_widgets[name], (gui_utils.TableColumnPicker, gui_utils.TableColumnPicker)):
                self.param_widgets[name].add_columns(self.filter_obj.columns)
            elif isinstance(self.param_widgets[name], gui_utils.ComboBoxOrOtherWidget) and isinstance(
                self.param_widgets[name].other, (gui_utils.TableColumnPicker, gui_utils.TableColumnPicker)):
                self.param_widgets[name].other.add_columns(self.filter_obj.columns)
            label = QtWidgets.QLabel(f'{name}:', self.param_widgets[name])
            label.setToolTip(this_desc)
            help_button = gui_utils.HelpButton()
            self.param_grid.addWidget(help_button, i, 2)
            self.param_grid.addWidget(label, i, 0)
            self.param_grid.addWidget(self.param_widgets[name], i, 1)
            help_button.connect_param_help(name, this_desc)
            i += 1
        self.param_grid.setRowStretch(i, 1)

    def init_setups_ui(self):
        self.setups_grid.addWidget(self.stack, 1, 0)
        self.setups_widgets['list'] = gui_utils.MultiChoiceListWithDelete(list(), parent=self.setups_group)
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

    def get_analysis_params(self):
        kwargs = {}
        for param_name, widget in self.param_widgets.items():
            if param_name in {'help_link'}:
                continue
            val = gui_utils.get_val_from_widget(widget)
            kwargs[param_name] = val
        return kwargs

    def start_clustering(self):
        kwargs = self.get_analysis_params()
        self.close()
        try:
            self.paramsAccepted.emit(self.parameter_dicts, kwargs)
        finally:
            self.show()


class EnrichmentWindow(gui_utils.MinMaxDialog):
    EXCLUDED_PARAMS = {'self', 'save_csv', 'fname', 'return_fig', 'biotype', 'background_genes',
                       'statistical_test', 'parametric_test', 'biotype_ref_path'}

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
                         'Hypergeometric test': 'hypergeometric', 'Single-list enrichment (XL-mHG test)': 'single_set'}

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
                 'go': {'plot_horizontal', 'plot_plot_ontology_graph', 'ontology_graph_format'},
                 'kegg': {'plot_horizontal', 'plot_pathway_graphs', 'pathway_graphs_format'},
                 'non_categorical': {'plot_log_scale', 'plot_style', 'n_bins'}}

    enrichmentFinished = QtCore.pyqtSignal(pd.DataFrame, str)

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

        self.scroll_layout = QtWidgets.QVBoxLayout(self.scroll)
        self.scroll_widget = QtWidgets.QWidget(self.scroll)
        self.main_layout = QtWidgets.QVBoxLayout(self)
        self.init_basic_ui()

    def init_basic_ui(self):
        self.setWindowTitle('Enrichment Analysis')
        self.setLayout(self.main_layout)
        self.main_layout.addWidget(self.scroll)

        self.parameter_group.setVisible(False)
        self.plot_group.setVisible(False)
        self.stats_group.setVisible(False)

        # self.scroll.setVerticalScrollBarPolicy(QtCore.Qt.ScrollBarAlwaysOn)
        self.scroll.setHorizontalScrollBarPolicy(QtCore.Qt.ScrollBarAlwaysOff)
        self.scroll.setWidgetResizable(True)
        self.scroll.setWidget(self.scroll_widget)
        self.scroll_widget.setLayout(self.scroll_layout)
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
            self.widgets[lst] = gui_utils.MandatoryComboBox('Choose gene set...', self)
            for obj_name in self.available_objects:
                self.widgets[lst].addItem(self.available_objects[obj_name][1], obj_name)

        self.widgets['dataset_radiobox'] = gui_utils.RadioButtonBox('Choose enrichment dataset',
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

    def _get_analysis_type(self):
        return self.ANALYSIS_TYPES[self.widgets['dataset_radiobox'].checkedButton().text()]

    def _update_single_set(self):
        if self.is_single_set():
            self._set_background_select_mode(False)
        else:
            self._set_background_select_mode(True)

    def update_uis(self):
        self.parameters_signature = {}
        self.stats_signature = {}
        self.plot_signature = {}

        analysis_type = self._get_analysis_type()
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
        gui_utils.clear_layout(self.plot_grid)

        i = 0
        for name, (param, desc) in self.plot_signature.items():
            self.plot_widgets[name] = gui_utils.param_to_widget(param, name)
            label = QtWidgets.QLabel(f'{name}:', self.plot_widgets[name])
            label.setToolTip(desc)
            help_button = gui_utils.HelpButton()
            self.plot_grid.addWidget(label, i, 0)
            self.plot_grid.addWidget(self.plot_widgets[name], i, 1)
            self.plot_grid.addWidget(help_button, i, 2)
            help_button.connect_param_help(name, desc)

            i += 1

    def update_parameters_ui(self):
        self.parameter_group.setVisible(True)
        # delete previous widgets
        self.parameter_widgets = {}
        gui_utils.clear_layout(self.parameter_grid)

        i = 0
        for name, (param, desc) in self.parameters_signature.items():
            self.parameter_widgets[name] = gui_utils.param_to_widget(param, name)
            label = QtWidgets.QLabel(f'{name}:', self.parameter_widgets[name])
            label.setToolTip(desc)
            help_button = gui_utils.HelpButton()
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
        gui_utils.clear_layout(self.stats_grid)

        radio_options = parsing.data_to_list(self.STATISTICAL_TESTS.keys() if self.is_categorical() else \
                                                 self.ORDINAL_STATISTICAL_TESTS.keys())
        self.stats_widgets['stats_radiobox'] = gui_utils.RadioButtonBox('Choose statistical test:', radio_options)
        if prev_test_name is not None:
            self.stats_widgets['stats_radiobox'].set_selection(prev_test_name)
        self.stats_widgets['stats_radiobox'].buttonClicked.connect(self.update_stats_ui)

        self.stats_grid.addWidget(self.stats_widgets['stats_radiobox'], 0, 0, 3, 2)

        i = 0
        for name, (param, desc) in self.stats_signature.items():
            if name in self.STATISTICAL_TEST_ARGS[prev_test]:
                self.stats_widgets[name] = gui_utils.param_to_widget(param, name)
                label = QtWidgets.QLabel(f'{name}:', self.stats_widgets[name])
                label.setToolTip(desc)
                help_button = gui_utils.HelpButton()
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
        gene_set = self.available_objects[gene_set_name][0]

        for param_name, widget in itertools.chain(self.parameter_widgets.items(), self.plot_widgets.items(),
                                                  self.stats_widgets.items()):
            if param_name in {'help_link', 'dataset_radiobox', 'stats_radiobox'}:
                continue
            val = gui_utils.get_val_from_widget(widget)
            kwargs[param_name] = val

        if not self.is_single_set():
            bg_set_name = self.widgets['bg_list'].currentText()
            bg_set = self.available_objects[bg_set_name][0]
        else:
            bg_set = None
        return gene_set, bg_set, gene_set_name.rstrip('*'), kwargs

    @QtCore.pyqtSlot()
    def run_analysis(self):
        func = self.get_current_func()
        gene_set, bg_set, set_name, kwargs = self.get_analysis_params()
        print("Enrichment analysis started")
        self.close()
        try:
            is_single_set = self.is_single_set()
            if is_single_set:
                bg_set_obj = None
            else:
                bg_set_obj = enrichment.FeatureSet(bg_set, 'background_set')
            feature_set_obj = enrichment.RankedSet(gene_set, set_name) if is_single_set \
                else enrichment.FeatureSet(gene_set, set_name)
            if is_single_set:
                result = func(feature_set_obj, **kwargs)
            else:
                result = func(feature_set_obj, background_genes=bg_set_obj, **kwargs)
        finally:
            self.show()
            self.enrichmentFinished.emit(result, set_name)


class SetOperationWindow(gui_utils.MinMaxDialog):
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
        set_names = [item.text() for item in self.widgets['set_list'].selectedItems()]
        sets = [self.available_objects[name][0] for name in set_names]
        ind = 0
        while ind < len(set_names):
            if sets[ind] is None:
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
                self._connect_canvas(canvas)

            else:
                canvas = gui_graphics.UpSetInteractiveCanvas(items, self)
                self._connect_canvas(canvas)
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

    def _connect_canvas(self, canvas: gui_graphics.BaseInteractiveCanvas):
        canvas.manualChoice.connect(self._set_op_other)
        self.widgets['radio_button_box'].radio_buttons['Union'].clicked.connect(canvas.union)
        self.widgets['radio_button_box'].radio_buttons['Intersection'].clicked.connect(canvas.intersection)
        self.primarySetChangedDifference.connect(canvas.difference)
        self.primarySetChangedIntersection.connect(canvas.intersection)

        if isinstance(canvas, gui_graphics.VennInteractiveCanvas):
            self.widgets['radio_button_box'].radio_buttons['Symmetric Difference'].clicked.connect(
                canvas.symmetric_difference)

    def init_ui(self):
        self.setWindowTitle('Set Operations')
        self.setGeometry(200, 200, 1250, 500)
        self.setLayout(self.layout)
        self.widgets['splitter'] = QtWidgets.QSplitter(QtCore.Qt.Horizontal)
        self.layout.addWidget(self.widgets['splitter'])
        self.widgets['splitter'].addWidget(self.list_group)
        self.widgets['splitter'].addWidget(self.operations_group)
        self.widgets['splitter'].setSizes([self.width() * 0.2, self.width() * 0.8])
        self.parameter_group.setVisible(False)

        self.init_sets_ui()
        self.init_operations_ui()

    def init_sets_ui(self):
        self.widgets['set_list'] = gui_utils.MultipleChoiceList(self.available_objects,
                                                                [val[1] for val in self.available_objects.values()],
                                                                self)

        for func in [self.create_canvas, self._check_legal_operations, self._validate_input,
                     self._toggle_choose_primary_set]:
            self.widgets['set_list'].itemSelectionChanged.connect(func)
        self.list_grid.addWidget(self.widgets['set_list'], 0, 0)

    def _toggle_choose_primary_set(self):
        if self.get_current_func_name() in ['difference', 'intersection']:
            self.widgets['choose_primary_set'].setVisible(True)
            self.widgets['choose_primary_set_label'].setVisible(True)

            self.widgets['choose_primary_set'].clear()
            self.widgets['choose_primary_set'].addItems(
                [item.text() for item in self.widgets['set_list'].selectedItems()])
            self.widgets['canvas'].clear_selection()
        else:
            self.widgets['choose_primary_set'].setVisible(False)
            self.widgets['choose_primary_set_label'].setVisible(False)

    def init_operations_ui(self):
        self.widgets['set_op_box'] = QtWidgets.QWidget(self)
        self.widgets['set_op_box_layout'] = QtWidgets.QVBoxLayout(self.widgets['set_op_box'])
        self.operations_grid.addWidget((self.widgets['set_op_box']), 1, 0, 3, 1)
        self.widgets['radio_button_box'] = gui_utils.RadioButtonBox('Choose set operation', self.SET_OPERATIONS.keys())

        for func in [self.update_paremeter_ui, self._validate_input, self._toggle_choose_primary_set]:
            self.widgets['radio_button_box'].buttonClicked.connect(func)
            self.widgets['radio_button_box'].selectionChanged.connect(func)

        self.widgets['set_op_box_layout'].addWidget(self.widgets['radio_button_box'], stretch=1)

        self.widgets['choose_primary_set'] = gui_utils.MandatoryComboBox('Choose primary set...', self)
        self.widgets['choose_primary_set'].currentTextChanged.connect(self._primary_set_changed)
        self.widgets['choose_primary_set_label'] = QtWidgets.QLabel('Primary set for operation:')
        self.widgets['set_op_box_layout'].addWidget(self.widgets['choose_primary_set_label'])
        self.widgets['set_op_box_layout'].addWidget(self.widgets['choose_primary_set'])

        self.widgets['set_op_box_layout'].addWidget(self.parameter_group)

        self.widgets['set_op_box_layout'].addStretch(1)

        self._toggle_choose_primary_set()

        self.create_canvas()

        self.widgets['apply_button'] = QtWidgets.QPushButton(text='Apply')
        self.widgets['apply_button'].clicked.connect(self.apply_set_op)
        self.widgets['apply_button'].setEnabled(False)
        self.operations_grid.addWidget(self.widgets['apply_button'], 4, 0, 1, 6)

    def _majority_vote_intersection(self):
        if not isinstance(self.widgets['canvas'], gui_graphics.EmptyCanvas):
            threshold = gui_utils.get_val_from_widget(self.parameter_widgets['majority_threshold'])
            self.widgets['canvas'].majority_vote_intersection(threshold)

    @QtCore.pyqtSlot(str)
    def _primary_set_changed(self, set_name: str):
        func_name = self.get_current_func_name()
        if func_name == 'difference':
            self.primarySetChangedDifference.emit(set_name)
        elif func_name == 'intersection':
            self.primarySetChangedIntersection.emit()

    def update_paremeter_ui(self):
        # delete previous widgets
        try:
            self.parameter_widgets['help_link'].deleteLater()
            self.operations_grid.removeWidget(self.parameter_widgets['help_link'])
        except KeyError:
            pass
        self.parameter_widgets = {}
        gui_utils.clear_layout(self.parameter_grid)

        chosen_func_name = self.get_current_func_name()
        signature = generic.get_method_signature(chosen_func_name, filtering.Filter)
        i = 0
        for name, param in signature.items():
            if name in self.EXCLUDED_PARAMS:
                continue
            self.parameter_widgets[name] = gui_utils.param_to_widget(param, name)
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
        n_items = len(self.widgets['set_list'].selectedItems())
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
        set_names = [item.text() for item in self.widgets['set_list'].selectedItems()]
        if self.get_current_func_name() in ['intersection', 'difference']:
            primary_set_name = self.widgets['choose_primary_set'].currentText()
            self.primarySetUsed.emit(primary_set_name)
        else:
            primary_set_name = set_names[0]
        kwargs = {}
        for param_name, widget in self.parameter_widgets.items():
            if param_name in {'apply_button', 'help_link'}:
                continue
            val = gui_utils.get_val_from_widget(widget)

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

            first_obj = self.available_objects[primary_set_name][0]
            if isinstance(first_obj, set):
                first_obj = filtering.Filter(('placeholder', pd.DataFrame(index=first_obj)))
            other_objs = []
            for name in set_names:
                if name != primary_set_name:
                    other_objs.append(self.available_objects[name][0])
            output_set = getattr(first_obj, func_name)(*other_objs, **kwargs)
            output_name = f"{func_name} output"
        if isinstance(output_set, set):
            self.geneSetReturned.emit(output_set, output_name)
        self.close()


class SetVisualizationWindow(gui_utils.MinMaxDialog):
    VISUALIZATION_FUNCS = {'Venn Diagram': 'venn_diagram', 'UpSet Plot': 'upset_plot'}
    EXCLUDED_PARAMS = {'objs', 'ref', 'fig'}

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
        self.setGeometry(200, 200, 1250, 500)
        self.setLayout(self.layout)
        self.widgets['splitter'] = QtWidgets.QSplitter(QtCore.Qt.Horizontal)
        self.layout.addWidget(self.widgets['splitter'])
        self.widgets['splitter'].addWidget(self.list_group)
        self.widgets['splitter'].addWidget(self.visualization_group)
        self.widgets['splitter'].setSizes([int(self.width() * 0.2), int(self.width() * 0.8)])

        self.parameter_group.setVisible(False)

        self.init_list_ui()
        self.init_visualization_ui()

    def init_list_ui(self):
        self.widgets['set_list'] = gui_utils.MultipleChoiceList(self.available_objects,
                                                                [val[1] for val in self.available_objects.values()],
                                                                self)

        for func in [self._check_legal_operations, self._validate_input, self.create_canvas]:
            self.widgets['set_list'].itemSelectionChanged.connect(func)

        self.list_grid.addWidget(self.widgets['set_list'], 0, 0)

    def init_visualization_ui(self):
        self.widgets['radio_button_box'] = gui_utils.RadioButtonBox('Choose visualization type:',
                                                                    self.VISUALIZATION_FUNCS, parent=self)
        for func in [self.update_parameter_ui, self._validate_input, self.create_canvas]:
            self.widgets['radio_button_box'].buttonClicked.connect(func)
            self.widgets['radio_button_box'].selectionChanged.connect(func)

        self.visualization_grid.addWidget(self.widgets['radio_button_box'], 0, 0, 2, 1)
        self.visualization_grid.addWidget(self.parameter_group, 2, 0, 2, 1)
        self.visualization_grid.setRowStretch(self.visualization_grid.count(), 1)

        self.widgets['generate_button'] = QtWidgets.QPushButton(text='Generate graph')
        self.widgets['generate_button'].clicked.connect(self.generate_graph)
        self.widgets['generate_button'].setEnabled(False)
        self.visualization_grid.addWidget(self.widgets['generate_button'], 4, 0, 1, 5)

        self.create_canvas()

    def create_canvas(self):
        set_names = [item.text() for item in self.widgets['set_list'].selectedItems()]
        sets = [self.available_objects[name][0] for name in set_names]
        ind = 0
        while ind < len(set_names):
            if sets[ind] is None:
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

    def select_all(self):
        for ind in range(self.widgets['set_list'].count()):
            item = self.widgets['set_list'].item(ind)
            if not item.isSelected():
                item.setSelected(True)

    def clear_all(self):
        for item in self.widgets['set_list'].selectedItems():
            item.setSelected(False)

    def _validate_input(self):
        is_legal = True

        if self.get_current_func_name() is None:
            is_legal = False

        if len(self.widgets['set_list'].selectedItems()) < 2:
            is_legal = False

        self.widgets['generate_button'].setEnabled(is_legal)

    def _check_legal_operations(self):
        n_items = len(self.widgets['set_list'].selectedItems())
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
        gui_utils.clear_layout(self.parameter_grid)

        chosen_func_name = self.get_current_func_name()
        signature = generic.get_method_signature(chosen_func_name, enrichment)
        i = 0
        for name, param in signature.items():
            if name in self.EXCLUDED_PARAMS:
                continue
            self.parameter_widgets[name] = gui_utils.param_to_widget(param, name, actions_to_connect=self.create_canvas)
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
        set_names = [item.text() for item in self.widgets['set_list'].selectedItems()]
        objs_to_plot = {key: val[0] for key, val in self.available_objects.items() if key in set_names}
        kwargs = {}
        for param_name, widget in self.parameter_widgets.items():
            if param_name in {'apply_button', 'help_link'}:
                continue
            val = gui_utils.get_val_from_widget(widget)

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
    commandIssued = QtCore.pyqtSignal(QtWidgets.QUndoCommand)

    def __init__(self, parent=None):
        super().__init__(parent)
        self.sup_layout = QtWidgets.QVBoxLayout(self)
        self.container = QtWidgets.QWidget(self)
        self.layout = QtWidgets.QVBoxLayout(self.container)
        self.name = None
        self.creation_time = time.time()

        self.stdout_group = QtWidgets.QGroupBox('Log')

        self.splitter = QtWidgets.QSplitter(QtCore.Qt.Vertical)
        self.sup_layout.addWidget(self.splitter)
        self.splitter.addWidget(self.container)
        self.splitter.addWidget(self.stdout_group)
        # self.splitter.setSizes([int(self.width() * 0.2), int(self.width() * 0.8)

        # self.layout.addWidget(self.stdout_group)
        self.sup_layout.addStretch(1)
        self.stdout_grid = QtWidgets.QGridLayout(self.stdout_group)
        self.stdout_widgets = {}

        self.init_stdout_ui()

    def is_empty(self):
        return True

    def init_stdout_ui(self):
        self.stdout_widgets['text_edit_stdout'] = gui_utils.StdOutTextEdit(self)
        self.stdout_widgets['text_edit_stdout'].setStyleSheet("""QTextEdit {background: #dddddd;}""")
        self.stdout_grid.addWidget(self.stdout_widgets['text_edit_stdout'], 0, 0, 3, 4)

    def get_console(self):
        return self.stdout_widgets['text_edit_stdout']

    def get_pbar(self):
        return self.stdout_widgets['progress_bar']

    @QtCore.pyqtSlot()
    def rename(self, new_name: str = None):
        if new_name is None:
            new_name = self.overview_widgets['table_name'].text()
        prev_name = self.get_tab_name()
        command = RenameCommand(prev_name, new_name, self, f'Rename "{prev_name}" to "{new_name}"')
        self.commandIssued.emit(command)

    def _rename(self, new_name: str = None):
        self.tabNameChange.emit(new_name, True)
        self.overview_widgets['table_name_label'].setText(f"Table name: '<b>{new_name}</b>'")
        self.overview_widgets['table_name'].setText('')

    def _get_parent_window(self):
        parent = self.parent()
        while not isinstance(parent, MainWindow):
            parent = parent.parent()
        return parent

    def _get_parent_tabwidget(self):
        parent = self.parent()
        while not isinstance(parent, QtWidgets.QTabWidget):
            parent = parent.parent()
        return parent

    def get_tab_name(self):
        parent = self._get_parent_tabwidget()
        return str(parent.tabText(parent.currentIndex())).rstrip("*")

    # custom method to write anything printed out to console/terminal to my QTextEdit widget via append function.
    def output_terminal_written(self, text):
        self.stdout_widgets['stdout'].append(text)


class SetTabPage(TabPage):
    def __init__(self, set_name: str, gene_set: typing.Union[set, enrichment.FeatureSet] = None, parent=None):
        super().__init__(parent)
        if gene_set is None:
            gene_set = enrichment.FeatureSet(set(), set_name)
        elif isinstance(gene_set, set):
            gene_set = enrichment.FeatureSet(gene_set, set_name)
        self.gene_set = gene_set

        self.overview_group = QtWidgets.QGroupBox('Data overview')
        self.overview_grid = QtWidgets.QGridLayout(self.overview_group)
        self.overview_widgets = {}
        self.init_overview_ui(set_name)

    def init_overview_ui(self, set_name: str):
        this_row = 0
        self.layout.insertWidget(0, self.overview_group)
        self.overview_widgets['table_name_label'] = QtWidgets.QLabel(f"Gene set name: '<b>{set_name}</b>'")

        self.overview_widgets['preview'] = QtWidgets.QListWidget()
        self.update_set_preview()
        self.overview_grid.addWidget(self.overview_widgets['table_name_label'], this_row, 0, 1, 1)
        this_row += 1
        self.overview_widgets['table_name'] = QtWidgets.QLineEdit()
        self.overview_widgets['rename_label'] = QtWidgets.QLabel('Rename your gene set (optional):')
        self.overview_widgets['rename'] = QtWidgets.QPushButton('Rename')
        self.overview_widgets['rename'].clicked.connect(self.rename)

        self.overview_grid.addWidget(self.overview_widgets['rename_label'], this_row, 0)
        self.overview_grid.addWidget(self.overview_widgets['table_name'], this_row, 1)
        self.overview_grid.addWidget(self.overview_widgets['rename'], this_row, 2)
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

    def view_full_gene_set(self):
        set_window = gui_utils.GeneSetView(self.gene_set.gene_set, self.get_tab_name())
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

    def is_empty(self):
        return self.gene_set is None

    def update_set_shape(self):
        if self.gene_set is not None:
            shape = len(self.gene_set)
            self.overview_widgets['shape'].setText(f'This gene set contains {shape} features')

    def update_set_preview(self):
        if self.gene_set is not None:
            self.overview_widgets['preview'].addItems([str(item) for item in self.gene_set])

    def update_gene_set(self, gene_set: set):
        self.gene_set.gene_set = gene_set
        self.update_set_shape()
        self.update_set_preview()
        self.changeIcon.emit('set')


class FuncTypeStack(QtWidgets.QWidget):
    EXCLUDED_PARAMS = {'self'}
    NO_FUNC_CHOSEN_TEXT = "Choose a function..."
    funcSelected = QtCore.pyqtSignal(bool)

    def __init__(self, funcs: list, filter_obj: filtering.Filter, parent=None, additional_excluded_params: set = None):
        super().__init__(parent)
        self.parameter_widgets = {}
        self.layout = QtWidgets.QVBoxLayout(self)
        self.parameter_grid = QtWidgets.QGridLayout()
        self.func_combo = QtWidgets.QComboBox(self)
        self.func_help_button = gui_utils.HelpButton(self)
        self.func_combo_layout = QtWidgets.QHBoxLayout()
        self.funcs = funcs
        self.filter_obj = filter_obj
        self.excluded_params = self.EXCLUDED_PARAMS
        if additional_excluded_params is not None:
            self.excluded_params.update(additional_excluded_params)
        self.init_ui()

    def init_ui(self):
        self.layout.addLayout(self.func_combo_layout)
        self.func_combo_layout.addWidget(self.func_combo)
        self.func_combo_layout.addWidget(self.func_help_button)
        self._set_empty_tooltip()
        self.layout.addLayout(self.parameter_grid)
        self.layout.addStretch(1)
        self.func_combo.addItem(self.NO_FUNC_CHOSEN_TEXT)
        self.func_combo.addItems(self.funcs)
        self.func_combo.currentTextChanged.connect(self.update_parameter_ui)

    def _set_empty_tooltip(self):
        txt = f"Choose a function from this list to read its description. "
        self.func_combo.setToolTip(txt)
        self.func_help_button.connect_desc_help(txt)

    def deselect(self):
        self.func_combo.setCurrentIndex(0)

    def update_parameter_ui(self):
        # delete previous widgets
        gui_utils.clear_layout(self.parameter_grid)
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
            self.parameter_widgets[name] = gui_utils.param_to_widget(param, name)
            if isinstance(self.parameter_widgets[name], (gui_utils.TableColumnPicker, gui_utils.TableColumnPicker)):
                self.parameter_widgets[name].add_columns(self.filter_obj.columns)
            elif isinstance(self.parameter_widgets[name], gui_utils.ComboBoxOrOtherWidget) and isinstance(
                self.parameter_widgets[name].other, (gui_utils.TableColumnPicker, gui_utils.TableColumnPicker)):
                self.parameter_widgets[name].other.add_columns(self.filter_obj.columns)
            label = QtWidgets.QLabel(f'{name}:', self.parameter_widgets[name])
            if name in param_desc:
                label.setToolTip(param_desc[name])
                help_button = gui_utils.HelpButton()
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
            val = gui_utils.get_val_from_widget(widget)

            func_params[param_name] = val
        return func_params

    def get_function_name(self):
        name = self.func_combo.currentText()
        return name


class FilterTabPage(TabPage):
    EXCLUDED_FUNCS = {'union', 'intersection', 'majority_vote_intersection', 'difference', 'symmetric_difference',
                      'from_folder', 'save_txt', 'save_csv'}
    CLUSTERING_FUNCS = {'split_kmeans': 'K-Means', 'split_kmedoids': 'K-Medoids',
                        'split_hierarchical': 'Hierarchical (Agglomerative)', 'split_hdbscan': 'HDBSCAN',
                        'split_clicom': 'CLICOM (Ensemble)'}
    SUMMARY_FUNCS = {'describe', 'head', 'tail', 'biotypes', 'print_features'}
    GENERAL_FUNCS = {'sort', 'transform', 'fold_change'}
    filterObjectCreated = QtCore.pyqtSignal(filtering.Filter)

    def __init__(self, parent=None):
        super().__init__(parent)
        self.filter_obj = None
        self.df_views = []

        self.basic_group = QtWidgets.QGroupBox('Load a data table')
        self.basic_grid = QtWidgets.QGridLayout(self.basic_group)
        self.basic_widgets = {}

        self.overview_group = QtWidgets.QGroupBox('Data overview')
        self.overview_grid = QtWidgets.QGridLayout(self.overview_group)
        self.overview_widgets = {}

        self.function_group = QtWidgets.QGroupBox('Apply functions')
        self.function_grid = QtWidgets.QGridLayout(self.function_group)
        self.function_widgets = {}

        self.clicom_window = None

        self.init_basic_ui()

    @QtCore.pyqtSlot()
    def _rename(self, new_name: str = None):
        super()._rename(new_name)
        self.filter_obj.fname = Path(
            os.path.join(str(self.filter_obj.fname.parent), f"{new_name}{self.filter_obj.fname.suffix}"))

    def is_empty(self):
        return self.filter_obj is None

    def get_table_type(self):
        return FILTER_OBJ_TYPES_INV[type(self.filter_obj)]

    def update_table_name_label(self):
        self.overview_widgets['table_name_label'].setText(f"Table name: '<b>{self.get_tab_name()}</b>'")

    def init_overview_ui(self):
        this_row = 0
        self.layout.insertWidget(1, self.overview_group)
        self.overview_widgets['table_type_label'] = QtWidgets.QLabel(
            f"Table type: {self.get_table_type()}")
        self.overview_widgets['table_name_label'] = QtWidgets.QLabel()
        self.update_table_name_label()

        self.overview_widgets['preview'] = QtWidgets.QTableView()
        # self.overview_widgets['preview'].setFixedWidth(550)
        # self.overview_widgets['preview'].setFixedHeight(175)
        self.overview_widgets['preview'].setVerticalScrollBarPolicy(QtCore.Qt.ScrollBarAlwaysOff)
        self.overview_widgets['preview'].setHorizontalScrollBarPolicy(QtCore.Qt.ScrollBarAlwaysOff)
        # self.overview_widgets['preview'].setFlags(self.overview_widgets['preview'].flags() & ~QtCore.Qt.ItemIsEditable)
        self.overview_widgets['preview'].horizontalHeader().setStretchLastSection(True)
        self.overview_widgets['preview'].verticalHeader().setStretchLastSection(True)
        self.update_table_preview()
        self.overview_grid.addWidget(self.overview_widgets['table_name_label'], this_row, 0, 1, 1)
        this_row += 1
        self.overview_widgets['table_name'] = QtWidgets.QLineEdit()
        self.overview_widgets['rename_label'] = QtWidgets.QLabel('Rename your table (optional):')
        self.overview_widgets['rename'] = QtWidgets.QPushButton('Rename')
        self.overview_widgets['rename'].clicked.connect(self.rename)

        self.overview_grid.addWidget(self.overview_widgets['rename_label'], this_row, 0)
        self.overview_grid.addWidget(self.overview_widgets['table_name'], this_row, 1)
        self.overview_grid.addWidget(self.overview_widgets['rename'], this_row, 2)
        this_row += 1
        self.overview_grid.addWidget(self.overview_widgets['preview'], this_row, 0, 4, 4)
        this_row += 4

        self.overview_widgets['save_button'] = QtWidgets.QPushButton('Save table')
        self.overview_widgets['save_button'].clicked.connect(self.save_file)
        self.overview_grid.addWidget(self.overview_widgets['save_button'], this_row, 3, 2, 1)

        self.overview_widgets['shape'] = QtWidgets.QLabel()
        self.overview_grid.addWidget(self.overview_widgets['shape'], this_row, 0, 1, 2)
        self.update_filter_obj_shape()

        self.overview_widgets['view_button'] = QtWidgets.QPushButton('View full table')
        self.overview_widgets['view_button'].clicked.connect(self.view_full_dataframe)
        self.overview_grid.addWidget(self.overview_widgets['view_button'], this_row, 2, 2, 1)

        this_row += 1
        self.overview_grid.addWidget(self.overview_widgets['table_type_label'], this_row, 0, 1, 1)

    def update_filter_obj_shape(self):
        shape = self.filter_obj.shape
        self.overview_widgets['shape'].setText(f'This table contains {shape[0]} rows, {shape[1]} columns')

    def _change_start_button_state(self, is_legal: bool):
        self.basic_widgets['start_button'].setEnabled(is_legal)

    def init_basic_ui(self):
        self.layout.insertWidget(0, self.basic_group)
        self.basic_widgets['table_type_combo'] = QtWidgets.QComboBox()
        self.basic_widgets['table_type_combo'].addItems(FILTER_OBJ_TYPES.keys())

        self.basic_widgets['start_button'] = QtWidgets.QPushButton('Start')
        self.basic_widgets['start_button'].clicked.connect(self.start)
        self.basic_widgets['start_button'].setEnabled(False)

        self.basic_widgets['file_path'] = gui_utils.PathLineEdit()
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
        self.basic_grid.addWidget(self.basic_widgets['start_button'], 2, 0, 1, 4)

    def init_function_ui(self):
        self.layout.insertWidget(2, self.function_group)
        sorted_actions = self.get_all_actions()

        self.stack = QtWidgets.QStackedWidget(self)
        self.button_box = QtWidgets.QButtonGroup(self)
        self.stack_buttons = []
        self.stack_widgets = {}
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

    def view_full_dataframe(self):
        df_window = gui_utils.DataFrameView(self.filter_obj.df, self.filter_obj.fname)
        self.overview_widgets['full_table_view'] = df_window
        df_window.show()

    def update_table_preview(self):
        model = gui_utils.DataFramePreviewModel(self.filter_obj.df)
        self.overview_widgets['preview'].setModel(model)
        self.overview_widgets['preview'].resizeRowsToContents()

    def apply_function(self):
        this_stack: FuncTypeStack = self.stack.currentWidget()
        func_name = this_stack.get_function_name()
        func_params = this_stack.get_function_params()
        self._apply_function_from_params(func_name, args=[], kwargs=func_params)

    def _apply_function_from_params(self, func_name, args: list, kwargs: dict):
        prev_name = self.get_tab_name()
        result = getattr(self.filter_obj, func_name)(*args, **kwargs)
        self.update_filter_obj_shape()
        self.update_table_preview()
        if prev_name != self.filter_obj.fname.name:
            self.tabNameChange.emit(self.filter_obj.fname.stem, True)
            self.update_table_name_label()

        self._proccess_outputs(result, func_name)

    def _proccess_outputs(self, outputs, source_name: str = ''):
        if validation.isinstanceinh(outputs, filtering.Filter):
            self.filterObjectCreated.emit(outputs)
        elif isinstance(outputs, pd.DataFrame):
            self.df_views.append(gui_utils.DataFrameView(outputs, source_name))
            self.df_views[-1].show()
        elif isinstance(outputs, np.ndarray):
            df = pd.DataFrame(outputs)
            self._proccess_outputs(df, source_name)

        elif isinstance(outputs, (tuple, list)):
            if validation.isinstanceiter_inh(outputs, filtering.Filter):
                dialog = MultiKeepWindow(outputs, self)
                dialog.accepted.connect(functools.partial(self._multi_keep_window_accepted, dialog, source_name))
                dialog.show()
            else:
                for output in outputs:
                    self._proccess_outputs(output, source_name)
        elif isinstance(outputs, dict):
            for output, this_src_name in outputs.items():
                self._proccess_outputs(output, this_src_name)

    def _multi_keep_window_accepted(self, dialog: QtWidgets.QDialog, source_name: str):
        kept_outputs = dialog.result()
        for output in kept_outputs:
            self._proccess_outputs(output, source_name)

    def apply_pipeline(self, pipeline: filtering.Pipeline, pipeline_name: str):
        apply_msg = f"Do you want to apply Pipeline '{pipeline_name}' inplace?"
        reply = QtWidgets.QMessageBox.question(self, f"Apply Pipeline '{pipeline_name}'",
                                               apply_msg, QtWidgets.QMessageBox.Yes, QtWidgets.QMessageBox.No)

        if reply == QtWidgets.QMessageBox.Yes:
            inplace = True
        else:
            inplace = False
        prev_name = self.filter_obj.fname.name
        result = pipeline.apply_to(self.filter_obj, inplace)
        self.update_filter_obj_shape()
        self.update_table_preview()
        if prev_name != self.filter_obj.fname.name:
            self.tabNameChange.emit(self.filter_obj.fname.stem, True)

        self._proccess_outputs(result, 'pipeline')

    def get_index_string(self):
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
        self.filter_obj = FILTER_OBJ_TYPES[self.basic_widgets['table_type_combo'].currentText()](
            self.basic_widgets['file_path'].text())
        print(self.filter_obj)
        table_name_user_input = self.basic_widgets['table_name'].text()
        if table_name_user_input != '':
            new_name = table_name_user_input
        else:
            new_name = self.filter_obj.fname.stem
        self.tabNameChange.emit(new_name, False)

        self.init_overview_ui()
        self.init_function_ui()

        gui_utils.clear_layout(self.basic_grid)
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


class CreatePipelineWindow(gui_utils.MinMaxDialog, FilterTabPage):
    pipelineSaved = QtCore.pyqtSignal(str, filtering.Pipeline)
    pipelineExported = QtCore.pyqtSignal(str, filtering.Pipeline)

    def __init__(self, parent=None):
        super().__init__(parent)
        self.setLayout(self.layout)
        self.setWindowTitle(f'Create new Pipeline')
        self.setGeometry(500, 200, 500, 300)
        self.pipeline = None
        self.is_unsaved = False

    def init_basic_ui(self):
        self.layout.insertWidget(0, self.basic_group)
        self.basic_group.setTitle("Choose data table type for Pipeline")

        self.basic_widgets['table_type_combo'] = QtWidgets.QComboBox()
        self.basic_widgets['table_type_combo'].addItems(FILTER_OBJ_TYPES.keys())

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

    def init_function_ui(self):
        super().init_function_ui()
        self.function_group.setTitle("Add functions to Pipeline")

    def _apply_function_from_params(self, func_name, *args, **kwargs):
        self.pipeline.add_function(func_name, *args, **kwargs)
        self.update_pipeline_preview()
        self.is_unsaved = True

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
        self.is_unsaved = True

    def init_overview_ui(self):
        self.function_group.setTitle("Pipeline preview")
        self.layout.insertWidget(1, self.overview_group)
        self.overview_widgets['preview'] = QtWidgets.QPlainTextEdit()
        self.overview_widgets['preview'].setReadOnly(True)
        self.update_pipeline_preview()

        self.overview_grid.addWidget(QtWidgets.QLabel('Pipeline preview', self.overview_widgets['preview']), 0, 0, 2, 4)
        self.overview_grid.addWidget(self.overview_widgets['preview'], 1, 0, 2, 4)

        self.overview_widgets['save_button'] = QtWidgets.QPushButton('Save Pipeline')
        self.overview_widgets['save_button'].clicked.connect(self.save_file)
        self.overview_grid.addWidget(self.overview_widgets['save_button'], 3, 3)

        self.overview_widgets['export_button'] = QtWidgets.QPushButton('Export Pipeline')
        self.overview_widgets['export_button'].clicked.connect(self.export_pipeline)
        self.overview_grid.addWidget(self.overview_widgets['export_button'], 3, 2, 1, 1)

    def _get_pipeline_name(self):
        return self.basic_widgets['pipeline_name'].text()

    def export_pipeline(self):
        self.pipelineExported.emit(self._get_pipeline_name(), self.pipeline)

    def save_file(self):
        try:
            self.pipelineSaved.emit(self._get_pipeline_name(), self.pipeline)
            print(f"Successfully saved Pipeline '{self.basic_widgets['pipeline_name'].text()}'")
            self.is_unsaved = False
        except Exception:
            print("Failed to save Pipeline")

    def closeEvent(self, event):
        if self.is_unsaved:
            quit_msg = "Are you sure you want to close this window without saving your Pipeline?\n" \
                       "All unsaved progress will be lost"

            reply = QtWidgets.QMessageBox.question(self, 'Close program',
                                                   quit_msg, QtWidgets.QMessageBox.No, QtWidgets.QMessageBox.Yes)

            if reply == QtWidgets.QMessageBox.Yes:
                event.accept()
            else:
                event.ignore()
        else:
            event.accept()


class MultiKeepWindow(gui_utils.MinMaxDialog):
    def __init__(self, objs: List[filtering.Filter], parent=None):
        super().__init__(parent)
        self.objs = {str(obj.fname.stem): obj for obj in objs}
        self.files = list(self.objs.keys())
        self.layout = QtWidgets.QGridLayout(self)
        self.button_box = QtWidgets.QDialogButtonBox(QtWidgets.QDialogButtonBox.Ok | QtWidgets.QDialogButtonBox.Cancel)
        self.labels = dict()
        self.keep_marks = dict()
        self.names = dict()
        self.select_all = QtWidgets.QCheckBox("Select all")

        self.init_ui()

    def init_ui(self):
        self.select_all.clicked.connect(self.change_all)
        self.button_box.accepted.connect(self.accept)
        self.button_box.rejected.connect(self.reject)
        self.setWindowTitle('Choose tables to keep')
        self.layout.addWidget(
            QtWidgets.QLabel(
                'Please choose which tables to keep (mandatory) out of the tables generated by the last function call, '
                'and _rename them (optional):\n\n'), 0, 0, 1, 3)
        self.layout.addWidget(QtWidgets.QLabel('<b>Tables:</b>'), 1, 0)
        self.layout.addWidget(QtWidgets.QLabel('<b>Choose tables to keep</b>'), 1, 1)
        self.layout.addWidget(QtWidgets.QLabel('<b>Rename tables (optional):</b>'), 1, 2)
        self.layout.addWidget(self.select_all, 2, 1)
        for i, file in enumerate(self.files, start=3):
            self.labels[file] = QtWidgets.QLabel(file)
            self.keep_marks[file] = QtWidgets.QCheckBox("Keep table?")
            self.names[file] = QtWidgets.QLineEdit()

            self.layout.addWidget(self.labels[file], i, 0)
            self.layout.addWidget(self.keep_marks[file], i, 1)
            self.layout.addWidget(self.names[file], i, 2)
        self.layout.addWidget(self.button_box, i + 1, 0, 1, 3)

    def change_all(self):
        for widget in self.keep_marks.values():
            widget.setChecked(self.select_all.isChecked())

    def result(self):
        keep_tables = {file: widget.isChecked() for file, widget in self.keep_marks.items()}
        new_names = {file: gui_utils.get_val_from_widget(widget) for file, widget in self.names.items()}
        out = []
        for file in keep_tables:
            if keep_tables[file]:
                name = str(file) if new_names[file] == '' else new_names[file]
                obj = self.objs[file]
                obj.fname = Path(name)
                out.append(obj)
        return out


class MultiOpenWindow(QtWidgets.QDialog):
    def __init__(self, files: List[str], parent=None):
        super().__init__(parent)
        self.files = files
        self.layout = QtWidgets.QGridLayout(self)
        self.button_box = QtWidgets.QDialogButtonBox(QtWidgets.QDialogButtonBox.Ok | QtWidgets.QDialogButtonBox.Cancel)
        self.paths = dict()
        self.table_types = dict()
        self.names = dict()

        self.init_ui()

    def init_ui(self):
        self.button_box.accepted.connect(self.accept)
        self.button_box.rejected.connect(self.reject)
        self.setWindowTitle('Choose table types and names')
        self.layout.addWidget(
            QtWidgets.QLabel(
                'Please choose a table type (mandatory) and table name (optional) for each loaded table\n\n'),
            0, 0, 1, 3)
        self.layout.addWidget(QtWidgets.QLabel('Table paths:'), 1, 0)
        self.layout.addWidget(QtWidgets.QLabel('Table types:'), 1, 1)
        self.layout.addWidget(QtWidgets.QLabel('Table names (optional):'), 1, 2)
        for i, file in enumerate(self.files):
            self.paths[file] = gui_utils.PathLineEdit(file, parent=self)
            self.table_types[file] = QtWidgets.QComboBox(self)
            self.table_types[file].addItems(list(FILTER_OBJ_TYPES.keys()))
            self.names[file] = QtWidgets.QLineEdit(self)

            self.layout.addWidget(self.paths[file], i + +2, 0)
            self.layout.addWidget(self.table_types[file], i + +2, 1)
            self.layout.addWidget(self.names[file], i + +2, 2)
        self.layout.addWidget(self.button_box, i + 3, 0, 1, 3)

    def result(self):
        paths = {file: gui_utils.get_val_from_widget(widget) for file, widget in self.paths.items()}
        types = {file: gui_utils.get_val_from_widget(widget) for file, widget in self.table_types.items()}
        names = {file: gui_utils.get_val_from_widget(widget) for file, widget in self.names.items()}
        return paths, types, names


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


class MainWindow(QtWidgets.QMainWindow):
    USER_GUIDE_URL = 'https://guyteichman.github.io/RNAlysis/build/user_guide.html'

    def __init__(self):
        super().__init__()
        self.tabs = ReactiveTabWidget(self)
        self.next_tab_n = 0

        self.undo_stack = QtWidgets.QUndoStack(self)
        self.undo_view = QtWidgets.QUndoView(self.undo_stack)
        self.undo_stack_widget = QtWidgets.QDockWidget('Command history', self)

        self.add_tab_button = QtWidgets.QToolButton()
        self.add_tab_button.setToolTip('Add New Tab')
        self.add_tab_button.clicked.connect(functools.partial(self.add_new_tab, name=None))
        self.add_tab_button.setText("+")
        self.error_window = None

        self.menu_bar = QtWidgets.QMenuBar(self)

        self.pipelines = OrderedDict()
        self.pipeline_window = None

        self.about_window = gui_utils.AboutWindow(self)
        self.settings_window = gui_utils.SettingsWindow(self)
        self.set_op_window = None
        self.set_visualization_window = None
        self.enrichment_window = None
        self.enrichment_results = []

        self.init_ui()
        self.add_new_tab()
        self.init_actions()
        self.init_menu_ui()
        self.init_shortcuts()

        self.queue_stdout = Queue()
        # create console text read thread + receiver object
        self.thread_stdout_queue_listener = QtCore.QThread()
        self.stdout_receiver = gui_utils.ThreadStdOutStreamTextQueueReceiver(self.queue_stdout)
        sys.stdout = gui_utils.WriteStream(self.queue_stdout)

        # connect receiver object to widget for text update
        self.stdout_receiver.queue_stdout_element_received_signal.connect(self.append_text_to_current_console)
        # attach console text receiver to console text thread
        self.stdout_receiver.moveToThread(self.thread_stdout_queue_listener)
        # attach to start / stop methods
        self.thread_stdout_queue_listener.started.connect(self.stdout_receiver.run)
        self.thread_stdout_queue_listener.start()

    def init_ui(self):
        self.setWindowTitle(f'RNAlysis {__version__}')
        self.setWindowIcon(QtGui.QIcon('../../docs/source/favicon.ico'))
        self.setGeometry(600, 50, 1000, 500)
        self.update_style_sheet()

        self.tabs.tabRightClicked.connect(self.init_tab_contextmenu)
        self.tabs.tabCloseRequested.connect(self.delete)

        self.undo_stack_widget.setWidget(self.undo_view)
        self.undo_stack_widget.setFloating(False)
        self.addDockWidget(QtCore.Qt.RightDockWidgetArea, self.undo_stack_widget)

        self.tabs.setCornerWidget(self.add_tab_button, QtCore.Qt.TopRightCorner)
        self.setCentralWidget(self.tabs)

    @QtCore.pyqtSlot(int)
    def init_tab_contextmenu(self, ind: int):
        tab_contextmenu = QtWidgets.QMenu(self)
        color_menu = tab_contextmenu.addMenu("Change tab &color")
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

        sort_menu = tab_contextmenu.addMenu("Sort tabs")
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

        tab_contextmenu.exec(QtGui.QCursor.pos())

    def update_style_sheet(self):
        self.setStyleSheet(gui_style.get_stylesheet())

    def sort_reverse(self):
        prev_order = {i: self.tabs.widget(i) for i in range(self.tabs.count())}
        self._sort_by_map(prev_order, reversed(prev_order.keys()))

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
                         filtering.Filter: 4, set: 5, type(None): 6}
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

    def _sort_by_map(self, key_map: dict, sorted_keys: list):
        for to_ind, name in enumerate(sorted_keys):
            widget = key_map[name]
            from_ind = self.tabs.indexOf(widget)

            self.tabs.tabBar().moveTab(from_ind, to_ind)

    def add_new_tab(self, name: str = None, is_set: bool = False):
        if name is None:
            name = f'New Table {self.next_tab_n + 1}'
            self.next_tab_n += 1

        else:
            print(name)

        if is_set:
            tab = SetTabPage(name, parent=self.tabs)
        else:
            tab = FilterTabPage(self.tabs)
            tab.filterObjectCreated.connect(self.new_tab_from_filter_obj)
        tab.changeIcon.connect(self.set_current_tab_icon)
        tab.tabNameChange.connect(self.rename_tab)
        tab.commandIssued.connect(self.undo_stack.push)
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
            self.new_tab_from_filter_obj(filter_obj)

    def load_multiple_files(self):
        dialog = gui_utils.MultiFileSelectionDialog()
        accepted = dialog.exec()
        if accepted == QtWidgets.QDialog.Accepted:
            filenames = dialog.result()
            if len(filenames) > 0:
                window = MultiOpenWindow(filenames, self)
                accepted = window.exec()
                if accepted:
                    paths, types, names = window.result()
                    if self.tabs.currentWidget().is_empty():
                        self.tabs.removeTab(self.tabs.currentIndex())
                    for filename in filenames:
                        path = paths[filename]
                        table_type = FILTER_OBJ_TYPES[types[filename]]
                        name = names[filename]
                        filter_obj = table_type(path)
                        if name == '':
                            self.new_tab_from_filter_obj(filter_obj)
                        else:
                            self.new_tab_from_filter_obj(filter_obj, name)

    def new_tab_from_filter_obj(self, filter_obj: filtering.Filter, name: str = None):
        self.add_new_tab(filter_obj.fname.name)
        self.tabs.currentWidget().start_from_filter_obj(filter_obj, name)

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

    def import_pipeline(self):
        filename, _ = QtWidgets.QFileDialog.getOpenFileName(self, "Choose a Pipeline file", str(Path.home()),
                                                            "YAML file (*.yaml)")
        if filename:
            pipeline = filtering.Pipeline.import_pipeline(filename)
            self.pipelines[str(Path(filename).stem)] = pipeline

    def import_multiple_gene_sets(self):
        dialog = gui_utils.MultiFileSelectionDialog()
        accepted = dialog.exec()
        if accepted == QtWidgets.QDialog.Accepted:
            filenames = dialog.result()
            if len(filenames) > 0 and self.tabs.currentWidget().is_empty():
                self.tabs.removeTab(self.tabs.currentIndex())
            for filename in filenames:
                gene_set = self._filename_to_gene_set(filename)
                self.new_tab_from_gene_set(gene_set, filename)

    def import_gene_set(self):
        filename, _ = QtWidgets.QFileDialog.getOpenFileName(self, "Choose a file", str(Path.home()),
                                                            "Text Document (*.txt);;"
                                                            "Comma-Separated Values (*.csv);;"
                                                            "Tab-Separated Values (*.tsv);;"
                                                            "All Files (*)")
        if filename:
            if self.tabs.currentWidget().is_empty():
                self.tabs.removeTab(self.tabs.currentIndex())
            gene_set = self._filename_to_gene_set(filename)
            self.new_tab_from_gene_set(gene_set, filename)

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

    def delete(self, index):
        self.tabs.removeTab(index)

    def copy_gene_set(self):
        gene_set = self.tabs.currentWidget().get_index_string()
        QtWidgets.QApplication.clipboard().setText(gene_set)

    def excepthook(self, exc_type, exc_value, exc_tb):
        sys.__excepthook__(exc_type, exc_value, exc_tb)
        self.error_window = gui_utils.ErrorMessage(exc_type, exc_value, exc_tb, self)
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

    def _get_current_pbar(self):

        current_pbar = self.tabs.currentWidget().get_pbar()
        return current_pbar

    @QtCore.pyqtSlot(str)
    def append_text_to_current_pbar(self, text: str):
        current_pbar = self._get_current_pbar()
        current_pbar.set_text(text)

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
        self.settings_action = QtWidgets.QAction("&Settings...", self)
        self.settings_action.triggered.connect(self.settings)
        self.exit_action = QtWidgets.QAction("&Exit", self)
        self.exit_action.triggered.connect(self.closeEvent)

        self.undo_action = self.undo_stack.createUndoAction(self)
        self.redo_action = self.undo_stack.createRedoAction(self)

        self.show_history_action = QtWidgets.QAction("Command &History")
        self.show_history_action.setCheckable(True)
        self.show_history_action.setChecked(True)
        self.show_history_action.triggered.connect(self.toggle_history_window)

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

        self.user_guide_action = QtWidgets.QAction("&User Guide", self)
        self.user_guide_action.triggered.connect(self.open_user_guide)
        self.about_action = QtWidgets.QAction("&About", self)
        self.about_action.triggered.connect(self.about)

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

        self.import_set_action.setShortcut(QtGui.QKeySequence("Ctrl+Shift+I"))
        self.import_multiple_sets_action.setShortcut(QtGui.QKeySequence("Ctrl+Shift+M"))
        self.export_set_action.setShortcut(QtGui.QKeySequence("Ctrl+Shift+E"))
        self.set_vis_action.setShortcut(QtGui.QKeySequence("Ctrl+Shift+V"))
        self.enrichment_action.setShortcut(QtGui.QKeySequence("Ctrl+Shift+A"))
        self.set_op_action.setShortcut(QtGui.QKeySequence("Ctrl+Shift+O"))

        self.new_pipeline_action.setShortcut(QtGui.QKeySequence("Ctrl+Alt+N"))
        self.import_pipeline_action.setShortcut(QtGui.QKeySequence("Ctrl+Alt+I"))
        self.export_pipeline_action.setShortcut(QtGui.QKeySequence("Ctrl+Alt+I"))

    @QtCore.pyqtSlot(bool)
    def toggle_history_window(self, state: bool):
        if state:
            self.undo_stack_widget.show()
        else:
            self.undo_stack_widget.close()

    def save_file(self):
        self.tabs.currentWidget().save_file()

    def export_gene_set(self):
        this_tab = self.tabs.currentWidget()
        if isinstance(this_tab, FilterTabPage):
            filter_obj = this_tab.filter_obj
            gene_set = filter_obj.index_set if filter_obj is not None else None
        else:
            gene_set = this_tab.gene_set
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
        prev_ind = self.tabs.currentIndex()
        try:
            self.tabs.setCurrentIndex(ind)
            gene_set = self.tabs.currentWidget().filter_obj if \
                isinstance(self.tabs.currentWidget(), FilterTabPage) else self.tabs.currentWidget().gene_set.gene_set
        finally:
            self.tabs.setCurrentIndex(prev_ind)
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
            available_objects_unique[key] = (self.get_gene_set_by_ind(i), self.tabs.tabIcon(i))
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
        df_window = gui_utils.DataFrameView(result, "Enrichment results for set " + gene_set_name)
        self.enrichment_results.append(df_window)
        df_window.show()

    def open_enrichment_analysis(self):
        self.enrichment_window = EnrichmentWindow(self.get_available_objects(), self)
        self.enrichment_window.enrichmentFinished.connect(self.display_enrichment_results)
        self.enrichment_window.show()

    def get_tab_names(self) -> List[str]:
        return [self.tabs.tabText(i) for i in range(self.tabs.count())]

    def init_menu_ui(self):
        self.setMenuBar(self.menu_bar)
        file_menu = self.menu_bar.addMenu("&File")
        self.new_menu = file_menu.addMenu("&New...")
        self.new_menu.addActions([self.new_table_action, self.new_multiple_action, self.new_table_from_folder_action])
        file_menu.addActions(
            [self.save_action, self.settings_action, self.exit_action])

        edit_menu = self.menu_bar.addMenu("&Edit")
        edit_menu.addActions([self.undo_action, self.redo_action])

        view_menu = self.menu_bar.addMenu("&View")
        view_menu.addActions([self.show_history_action])

        gene_sets_menu = self.menu_bar.addMenu("&Gene sets")
        gene_sets_menu.addActions(
            [self.copy_action, self.set_op_action, self.enrichment_action, self.set_vis_action, self.import_set_action,
             self.import_multiple_sets_action, self.export_set_action])

        pipeline_menu = self.menu_bar.addMenu("&Pipelines")
        pipeline_menu.addActions([self.new_pipeline_action, self.import_pipeline_action, self.export_pipeline_action])
        self.apply_pipeline_menu = pipeline_menu.addMenu("&Apply Pipeline")
        self.apply_pipeline_menu.aboutToShow.connect(self._populate_pipelines)

        help_menu = self.menu_bar.addMenu("&Help")
        help_menu.addActions([self.user_guide_action, self.about_action])

    def _populate_pipelines(self):
        # Remove the old options from the menu
        self.apply_pipeline_menu.clear()
        # Dynamically create the actions
        actions = []
        for name, pipeline in self.pipelines.items():
            action = QtWidgets.QAction(name, self)
            action.triggered.connect(
                functools.partial(self.tabs.currentWidget().apply_pipeline, pipeline, name))
            actions.append(action)
        # Step 3. Add the actions to the menu
        self.apply_pipeline_menu.addActions(actions)

    def about(self):
        self.about_window.exec()

    def input(self, message: str):
        dialog = gui_utils.PathInputDialog(message, parent=self)
        accepted = dialog.exec()
        if accepted:
            return dialog.result()
        return None

    def closeEvent(self, event):

        quit_msg = "Are you sure you want to close <i>RNAlysis</i>?\n" \
                   "All unsaved progress will be lost"

        reply = QtWidgets.QMessageBox.question(self, 'Close program',
                                               quit_msg, QtWidgets.QMessageBox.No, QtWidgets.QMessageBox.Yes)
        # reply.setWindowIcon(self.style().standardIcon(QtWidgets.QStyle.SP_MessageBoxQuestion))

        if reply == QtWidgets.QMessageBox.Yes:
            event.accept()
        else:
            event.ignore()


def customwarn(message, category, filename, lineno, file=None, line=None):
    sys.stdout.write(warnings.formatwarning(message, category, filename, lineno))


def run():
    matplotlib.use('Qt5Agg')
    app = QtWidgets.QApplication(sys.argv)
    window = MainWindow()
    warnings.showwarning = customwarn
    sys.excepthook = window.excepthook
    builtins.input = window.input

    # enrichment.enrichment_runner.generic.ProgressParallel = functools.partial(gui_utils.ProgressParallelGui,
    #                                                                           parent=window)
    enrichment.enrichment_runner.tqdm = functools.partial(gui_utils.ProgressSerialGui, parent=window)
    # generic.ProgressParallel = functools.partial(gui_utils.ProgressParallelGui, parent=window)
    enrichment.enrichment_runner.parsing.tqdm = functools.partial(gui_utils.ProgressSerialGui, parent=window)
    filtering.clustering.tqdm = functools.partial(gui_utils.ProgressSerialGui, parent=window)
    enrichment.enrichment_runner.io.tqdm = functools.partial(gui_utils.ProgressSerialGui, parent=window)
    # filtering.clustering.generic.ProgressParallel = functools.partial(gui_utils.ProgressParallelGui, parent=window)

    window.show()
    sys.exit(app.exec())


if __name__ == '__main__':
    run()
