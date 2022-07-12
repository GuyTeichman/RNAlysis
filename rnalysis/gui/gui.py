import functools
import sys
import os
from collections import OrderedDict
from pathlib import Path
from queue import Queue
import typing
import matplotlib.pyplot as plt
import matplotlib
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg, NavigationToolbar2QT
import matplotlib_venn
import numpy as np
import pandas as pd
from PyQt5 import QtCore, QtWidgets, QtGui
import itertools

from rnalysis import filtering, enrichment, __version__
from rnalysis.utils import io, validation, generic
from rnalysis.gui import gui_style, gui_utils
from typing import List, Union, Callable


class EnrichmentWindow(gui_utils.MinMaxDialog):
    IGNORED_ARGS = {'self', 'save_csv', 'fname', 'return_fig', 'parallel', 'biotype', 'background_genes',
                    'statistical_test', 'parametric_test'}

    ANALYSIS_TYPES = {'User-defined attributes': 'user_defined', 'Gene Ontology (GO)': 'go',
                      'Kyoto Encyclopedia of Genes and Genomes (KEGG) Pathways': 'kegg',
                      'Non-categorical variables': 'non_categorical'}

    ANALYSIS_FUNCS = {('go', False): enrichment.FeatureSet.go_enrichment,
                      ('go', True): enrichment.RankedSet.single_set_go_enrichment,
                      ('kegg', False): enrichment.FeatureSet.kegg_enrichment,
                      ('kegg', True): enrichment.RankedSet.single_set_kegg_enrichment,
                      ('user_defined', False): enrichment.FeatureSet.enrich_hypergeometric,  # TODO: update me!
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

    def __init__(self, available_objects: list, parent=None):
        super().__init__(parent)

        checked = {}
        available_objects_unique = []
        for obj in available_objects:
            if obj in checked:
                available_objects_unique.append(f"{obj}_{checked[obj]}")
            else:
                checked[obj] = 1
                available_objects_unique.append(obj)
            checked[obj] += 1
        self.available_objects = available_objects_unique

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
        self.setGeometry(200, 200, 1100, 600)
        self.setLayout(self.main_layout)
        self.main_layout.addWidget(self.scroll)

        self.parameter_group.setVisible(False)
        self.plot_group.setVisible(False)
        self.stats_group.setVisible(False)

        # self.scroll.setVerticalScrollBarPolicy(QtCore.Qt.ScrollBarAlwaysOn)
        # self.scroll.setHorizontalScrollBarPolicy(QtCore.Qt.ScrollBarAlwaysOff)
        self.scroll.setWidgetResizable(True)
        self.scroll.setWidget(self.scroll_widget)
        self.scroll_widget.setLayout(self.scroll_layout)

        self.scroll_layout.addWidget(self.list_group)
        self.scroll_layout.addWidget(self.parameter_group)
        for lst in ['list1', 'list2']:
            self.widgets[lst] = QtWidgets.QListWidget(self)
            self.widgets[lst].insertItems(0, self.available_objects)

        self.widgets['radio_button_box'] = gui_utils.RadioButtonBox('Choose analysis type', self.ANALYSIS_TYPES.keys())
        self.widgets['radio_button_box'].buttonClicked.connect(self.update_uis)

        self.widgets['run_button'] = QtWidgets.QPushButton('Run')
        self.widgets['run_button'].clicked.connect(self.run_analysis)

        self.list_grid.addWidget(self.widgets['radio_button_box'], 2, 0, 3, 1)
        self.list_grid.addWidget(self.widgets['list1'], 2, 1, 3, 1)
        self.list_grid.addWidget(self.widgets['list2'], 2, 2, 3, 1)
        self.list_grid.addWidget(QtWidgets.QLabel('Choose the enrichment set:', self), 1, 1, 1, 1)
        self.list_grid.addWidget(QtWidgets.QLabel('Choose background gene set:', self), 1, 2, 1, 1)

    def _set_background_select_mode(self, selectable: bool = True):
        if selectable:
            self.widgets['list2'].setSelectionMode(QtWidgets.QAbstractItemView.SingleSelection)
            self.widgets['list2'].setDisabled(False)
        else:
            self.widgets['list2'].setSelectionMode(QtWidgets.QAbstractItemView.NoSelection)
            self.widgets['list2'].setDisabled(True)
            self.widgets['list2'].clearSelection()

    def _get_statistical_test_name(self):
        try:
            button = self.stats_widgets['radio_button_box'].checkedButton()
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
        return self.ANALYSIS_TYPES[self.widgets['radio_button_box'].checkedButton().text()]

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
        for name, param in signature.items():
            if name in self.IGNORED_ARGS:
                continue
            elif name in self.PLOT_ARGS[analysis_type]:
                self.plot_signature[name] = param
            elif name in set.union(*self.STATISTICAL_TEST_ARGS.values()):
                self.stats_signature[name] = param
            else:
                self.parameters_signature[name] = param

        self.update_parameters_ui()
        self.update_stats_ui()
        self.update_plot_ui()
        self.scroll_layout.insertWidget(4, self.widgets['run_button'])

    def get_current_analysis_type(self):
        button = self.widgets['radio_button_box'].checkedButton()
        if button is None:
            return None
        return self.ANALYSIS_TYPES[self.widgets['radio_button_box'].checkedButton().text()]

    def get_current_func(self):
        single_set = self.is_single_set()
        name = self.get_current_analysis_type()
        func = self.ANALYSIS_FUNCS[(name, single_set)]
        return func

    @staticmethod
    def clear_layout(layout):
        while layout.count():
            child = layout.takeAt(0)
            if child.widget():
                child.widget().deleteLater()

    def update_plot_ui(self):
        self.plot_group.setVisible(True)
        # delete previous widgets
        self.plot_widgets = {}
        self.clear_layout(self.plot_grid)
        self.scroll_layout.insertWidget(3, self.plot_group)

        i = 0
        for name, param in self.plot_signature.items():
            self.plot_widgets[name] = gui_utils.param_to_widget(param, name)
            self.plot_grid.addWidget(QtWidgets.QLabel(f'{name}:', self.plot_widgets[name]), i, 0)
            self.plot_grid.addWidget(self.plot_widgets[name], i, 1)
            i += 1

    def update_parameters_ui(self):
        self.parameter_group.setVisible(True)
        # delete previous widgets
        self.parameter_widgets = {}
        self.clear_layout(self.parameter_grid)
        self.scroll_layout.insertWidget(2, self.parameter_group)

        i = 0
        for name, param in self.parameters_signature.items():
            self.parameter_widgets[name] = gui_utils.param_to_widget(param, name)
            self.parameter_grid.addWidget(QtWidgets.QLabel(f'{name}:', self.parameter_widgets[name]), i, 0)
            self.parameter_grid.addWidget(self.parameter_widgets[name], i, 1)
            i += 1

    def update_stats_ui(self):
        self.stats_group.setVisible(True)
        prev_test = self._get_statistical_test()
        prev_test_name = self._get_statistical_test_name()
        self._update_single_set()

        # delete previous widgets
        self.stats_widgets = {}
        self.clear_layout(self.stats_grid)
        self.scroll_layout.insertWidget(1, self.stats_group)

        radio_options = self.STATISTICAL_TESTS.keys() if self.is_categorical() else \
            self.ORDINAL_STATISTICAL_TESTS.keys()
        self.stats_widgets['radio_button_box'] = gui_utils.RadioButtonBox('Choose statistical test:', radio_options)
        if prev_test_name is not None:
            self.stats_widgets['radio_button_box'].set_selection(prev_test_name)
        self.stats_widgets['radio_button_box'].buttonClicked.connect(self.update_stats_ui)

        self.stats_grid.addWidget(self.stats_widgets['radio_button_box'], 0, 0, 3, 2)

        i = 0
        for name, param in self.stats_signature.items():
            if name in self.STATISTICAL_TEST_ARGS[prev_test]:
                self.stats_widgets[name] = gui_utils.param_to_widget(param, name)
                self.stats_grid.addWidget(QtWidgets.QLabel(f'{name}:', self.stats_widgets[name]), i, 2)
                self.stats_grid.addWidget(self.stats_widgets[name], i, 3)
                i += 1

        # help_address = f"https://guyteichman.github.io/RNAlysis/build/rnalysis.filtering." \
        #                f"{filtering.Filter.__name__}.{chosen_func_name}.html"  # TODO: fix me!
        # self.stats_widgets['help_link'] = QtWidgets.QLabel(
        #     text=f'<a href="{help_address}">Open documentation for function '
        #          f'<b>{filtering.Filter.__name__}.{chosen_func_name}</b></a>')
        # self.stats_widgets['help_link'].setOpenExternalLinks(True)
        # self.parameter_grid.addWidget(self.stats_widgets['help_link'], i + 1, 0, 1, 2)

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
        gene_set_name = self.widgets['list1'].currentItem().text()
        gene_set = self.available_objects.index(gene_set_name)

        for param_name, widget in itertools.chain(self.parameter_widgets.items(), self.plot_widgets.items(),
                                                  self.stats_widgets.items()):
            if param_name in {'help_link', 'radio_button_box'}:
                continue
            val = gui_utils.get_val_from_widget(widget)
            kwargs[param_name] = val

        if not self.is_single_set():
            bg_set_name = self.widgets['list2'].currentItem().text()
            bg_set = self.available_objects.index(bg_set_name)
        else:
            bg_set = None
        return gene_set, bg_set, kwargs

    @QtCore.pyqtSlot()
    def run_analysis(self):
        func = self.get_current_func()
        gene_set, bg_set, kwargs = self.get_analysis_params()

        self.parent().run_enrichment_analysis(func, gene_set, bg_set, kwargs)

        self.close()


# class SetOpWindow(gui_utils.MinMaxDialog):
#     SET_OPERATIONS = {'Intersection': 'intersection',
#                       'Majority-Vote Intersection': 'majority_vote_intersection',
#                       'Union': 'union', 'Difference': 'difference',
#                       'Symmetric Difference': 'symmetric_difference'}
#
#     def __init__(self, available_objects: dict, parent=None):
#         super().__init__(parent)
#         self.available_objects = available_objects
#         self.widgets = {}
#         self.list_group = QtWidgets.QGroupBox('Set operation', self)
#         self.list_grid = QtWidgets.QGridLayout(self.list_group)
#         self.parameter_group = QtWidgets.QGroupBox('Additional parameters', self)
#         self.parameter_grid = QtWidgets.QGridLayout(self.parameter_group)
#         self.parameter_widgets = {}
#         self.layout = QtWidgets.QHBoxLayout(self)
#
#         self.init_ui()
#
#     def init_ui(self):
#         self.setWindowTitle('Set Operations')
#         self.setLayout(self.layout)
#         self.layout.addWidget(self.list_group)
#         self.layout.addWidget(self.parameter_group)
#
#         self.parameter_group.setVisible(False)
#
#         for lst in ['list1', 'list2']:
#             self.widgets[lst] = QtWidgets.QListWidget(self)
#             self.widgets[lst].insertItems(0, self.available_objects)
#         self.widgets['list2'].setSelectionMode(QtWidgets.QAbstractItemView.MultiSelection)
#         self.widgets['radio_button_box'] = gui_utils.RadioButtonBox('Choose set operation', self.SET_OPERATIONS.keys())
#         self.widgets['radio_button_box'].buttonClicked.connect(self.update_paremeter_ui)
#
#         self.list_grid.addWidget(self.widgets['radio_button_box'], 2, 0, 3, 1)
#         self.list_grid.addWidget(self.widgets['list1'], 2, 1, 3, 1)
#         self.list_grid.addWidget(self.widgets['list2'], 2, 2, 3, 1)
#         self.list_grid.addWidget(QtWidgets.QLabel('Choose first gene set:', self), 1, 1, 1, 1)
#         self.list_grid.addWidget(QtWidgets.QLabel('Choose other gene set(s):', self), 1, 2, 1, 1)
#
#     def get_current_func_name(self):
#         return self.SET_OPERATIONS[self.widgets['radio_button_box'].checkedButton().text()]
#
#     @staticmethod
#     def clear_layout(layout):
#         while layout.count():
#             child = layout.takeAt(0)
#             if child.widget():
#                 child.widget().deleteLater()
#
#     def update_paremeter_ui(self):
#         self.parameter_group.setVisible(True)
#
#         # delete previous widgets
#         self.parameter_widgets = {}
#         self.clear_layout(self.parameter_grid)
#
#         chosen_func_name = self.get_current_func_name()
#         signature = generic.get_method_signature(chosen_func_name, filtering.Filter)
#         i = 0
#         for name, param in signature.items():
#             if name in {'self', 'other', 'others', 'return_type'}:
#                 continue
#             self.parameter_widgets[name] = gui_utils.param_to_widget(param, name)
#             self.parameter_grid.addWidget(QtWidgets.QLabel(f'{name}:', self.parameter_widgets[name]), i, 0)
#             self.parameter_grid.addWidget(self.parameter_widgets[name], i, 1)
#             i += 1
#
#         self.parameter_widgets['apply_button'] = QtWidgets.QPushButton(text='Apply')
#         self.parameter_widgets['apply_button'].clicked.connect(self.apply_set_op)
#         self.parameter_grid.addWidget(self.parameter_widgets['apply_button'], i, 0, 1, 2)
#
#         help_address = f"https://guyteichman.github.io/RNAlysis/build/rnalysis.filtering." \
#                        f"{filtering.Filter.__name__}.{chosen_func_name}.html"
#         self.parameter_widgets['help_link'] = QtWidgets.QLabel(
#             text=f'<a href="{help_address}">Open documentation for function '
#                  f'<b>{filtering.Filter.__name__}.{chosen_func_name}</b></a>')
#         self.parameter_widgets['help_link'].setOpenExternalLinks(True)
#         self.parameter_grid.addWidget(self.parameter_widgets['help_link'], i + 1, 0, 1, 2)
#
#     def _get_function_params(self):
#         first_name = self.widgets['list1'].currentItem().text()
#         first = self.available_objects.index(first_name)
#         second_names = [item.text() for item in self.widgets['list2'].selectedItems()]
#         second = [self.available_objects.index(name) for name in second_names]
#         kwargs = {}
#         for param_name, widget in self.parameter_widgets.items():
#             if param_name in {'apply_button', 'help_link'}:
#                 continue
#             val = gui_utils.get_val_from_widget(widget)
#
#             kwargs[param_name] = val
#         return first, second, kwargs
#
#     @QtCore.pyqtSlot()
#     def apply_set_op(self):
#         func_name = self.get_current_func_name()
#         first, seconds, kwargs = self._get_function_params()
#
#         self.parent().apply_set_op(func_name, first, seconds, kwargs)
#
#         self.close()


class EmptyCanvas(FigureCanvasQTAgg):
    def __init__(self, text: str, parent=None):
        self.fig = plt.Figure(constrained_layout=True)
        self.ax = self.fig.add_subplot()
        self.ax.text(0, 0.5, text, fontsize=15)
        super().__init__(self.fig)
        plt.close(self.fig)
        self.parent = parent

        for spine in ['right', 'left', 'top', 'bottom']:
            self.ax.spines[spine].set_visible(False)
        self.ax.set_xticks([])
        self.ax.set_yticks([])


class VennCanvas(FigureCanvasQTAgg):
    DESELECTED_STATE = 0
    SELECTED_STATE = 1
    HOVER_STATE = 2
    HOVER_SELECTED_STATE = 3

    SELECTED_COLOR = (0.15, 0.01, 0.35)
    HOVER_COLOR = (0.5, 0.5, 0.5)
    HOVER_SELECTED_COLOR = (0.4, 0.1, 0.75)

    manualChoice = QtCore.pyqtSignal()

    def __init__(self, gene_sets: dict, parent=None):
        self.parent = parent
        self.gene_sets = gene_sets
        self.fig = plt.Figure(tight_layout=True)
        self.ax = self.fig.add_subplot()
        super().__init__(self.fig)

        if len(gene_sets) == 2:
            funcs = matplotlib_venn.venn2, matplotlib_venn.venn2_circles
            colors = (self.SELECTED_COLOR, self.SELECTED_COLOR)

        elif len(gene_sets) == 3:
            funcs = matplotlib_venn.venn3, matplotlib_venn.venn3_circles
            colors = (self.SELECTED_COLOR, self.SELECTED_COLOR, self.SELECTED_COLOR)
        else:
            raise ValueError("Cannot proccess more than 3 sets!")

        self.venn = funcs[0](gene_sets.values(), gene_sets.keys(), set_colors=colors, ax=self.ax)
        self.venn_circles = funcs[1](gene_sets.values(), linestyle='solid', linewidth=2.0, ax=self.ax)
        self.set_font_size(16, 12)
        self.fig.canvas.mpl_connect('button_press_event', self.on_click)
        self.fig.canvas.mpl_connect('motion_notify_event', self.on_hover)

        self.clear_selection()

    def set_font_size(self, set_label_size: float, subset_label_size: float):
        for label in self.venn.set_labels:
            if label is not None:
                label.set_fontsize(set_label_size)
        for label in self.venn.subset_labels:
            if label is not None:
                label.set_fontsize(subset_label_size)

    def clear_selection(self, draw: bool = True):
        for patch in self.venn.patches:
            if patch is None:
                continue
            patch.set_fill(False)
            patch.state = self.DESELECTED_STATE
        if draw:
            self.draw()

    def select(self, ind, draw: bool = True):
        patch = self.venn.get_patch_by_id(ind)
        if patch is None:
            return
        patch.state = self.SELECTED_STATE
        patch.set_facecolor(self.SELECTED_COLOR)
        patch.set_fill(True)
        if draw:
            self.draw()

    def deselect(self, ind, draw: bool = True):
        patch = self.venn.get_patch_by_id(ind)
        if patch is None:
            return
        patch.state = self.DESELECTED_STATE
        patch.set_fill(False)
        if draw:
            self.draw()

    def flip(self, ind, draw: bool = True):
        patch = self.venn.get_patch_by_id(ind)
        if patch is None:
            return
        if patch.state in [self.SELECTED_STATE, self.HOVER_SELECTED_STATE]:
            patch.state = self.DESELECTED_STATE
            patch.set_fill(False)
        else:
            patch.state = self.SELECTED_STATE
            patch.set_facecolor(self.SELECTED_COLOR)
            patch.set_fill(True)

        if draw:
            self.draw()

    def on_click(self, event):
        for patch in self.venn.patches:
            if patch is None:
                continue

            if patch.contains_point((event.x, event.y)):
                self.manualChoice.emit()
                patch.set_facecolor(self.HOVER_SELECTED_COLOR)
                if patch.state in [self.DESELECTED_STATE, self.HOVER_STATE]:
                    patch.state = self.SELECTED_STATE
                    patch.set_fill(True)
                else:
                    patch.state = self.DESELECTED_STATE
                    patch.set_facecolor(self.HOVER_COLOR)
        self.draw()

    def on_hover(self, event):
        for patch in self.venn.patches:
            if patch is None:
                continue

            if patch.contains_point((event.x, event.y)):
                if patch.state in [self.HOVER_STATE, self.DESELECTED_STATE]:
                    patch.set_facecolor(self.HOVER_COLOR)
                    patch.state = self.HOVER_STATE
                else:
                    patch.set_facecolor(self.HOVER_SELECTED_COLOR)
                    patch.state = self.HOVER_SELECTED_STATE
                patch.set_fill(True)
            else:
                if patch.state in [self.HOVER_STATE, self.DESELECTED_STATE]:
                    patch.set_fill(False)
                    patch.state = self.DESELECTED_STATE
                else:
                    patch.set_facecolor(self.SELECTED_COLOR)
                    patch.state = self.SELECTED_STATE
                    patch.set_fill(True)
        self.draw()

    def get_tuple_patch_ids(self):
        return list(itertools.product([0, 1], repeat=len(self.gene_sets)))[1:]

    @QtCore.pyqtSlot()
    def union(self):
        for patch_id in self.get_tuple_patch_ids():
            str_patch_id = ''.join(str(i) for i in patch_id)
            self.select(str_patch_id, draw=False)
        self.draw()

    @QtCore.pyqtSlot()
    def intersection(self):
        self.clear_selection(draw=False)
        self.select("1" * len(self.gene_sets), draw=False)
        self.draw()

    @QtCore.pyqtSlot()
    def symmetric_difference(self):
        if len(self.gene_sets) != 2:
            return
        self.clear_selection(draw=False)
        self.select("10", draw=False)
        self.select("01", draw=False)
        self.draw()

    @QtCore.pyqtSlot(str)
    def difference(self, primary_set: str):
        self.clear_selection(draw=False)
        if primary_set in self.gene_sets:
            key = ""
            for set_name in self.gene_sets:
                key += "1" if set_name == primary_set else "0"
            self.select(key, draw=False)
        self.draw()

    @QtCore.pyqtSlot()
    def majority_vote_intersection(self, majority_threshold: float):
        if len(self.gene_sets) == 3:
            if majority_threshold <= (1 / 3):
                self.union()
            elif (1 / 3) < majority_threshold <= (2 / 3):
                self.clear_selection(draw=False)
                for ind in ["111", "110", "101", "011"]:
                    self.select(ind, draw=False)
                self.draw()
            else:
                self.intersection()
        else:
            if majority_threshold <= 0.5:
                self.union()
            else:
                self.intersection()

        self.draw()

    def get_custom_selection(self):
        selection = set()
        for patch_id in self.get_tuple_patch_ids():
            str_patch_id = ''.join(str(i) for i in patch_id)
            patch = self.venn.get_patch_by_id(str_patch_id)
            if patch is None:
                continue
            if patch.state in [self.SELECTED_STATE, self.HOVER_SELECTED_STATE]:
                included_sets = [s for s, ind in zip(self.gene_sets.values(), patch_id) if ind]
                excluded_sets = [s for s, ind in zip(self.gene_sets.values(), patch_id) if not ind]
                patch_content = set.intersection(*included_sets).difference(*excluded_sets)
                selection.add(patch_content)
        return selection


class SimpleSetOpWindow(gui_utils.MinMaxDialog):
    SET_OPERATIONS = {'Union': 'union', 'Majority-Vote Intersection': 'majority_vote_intersection',
                      'Intersection': 'intersection', 'Difference': 'difference',
                      'Symmetric Difference': 'symmetric_difference', 'Other': 'other'}

    geneSetReturned = QtCore.pyqtSignal(set, str)
    primarySetUsed = QtCore.pyqtSignal(str)
    primarySetChangedDifference = QtCore.pyqtSignal(str)
    primarySetChangedIntersection = QtCore.pyqtSignal()

    def __init__(self, available_objects: dict, parent=None):
        super().__init__(parent)
        self.available_objects = available_objects
        self.widgets = {}
        self.list_group = QtWidgets.QGroupBox('Set operation', self)
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
        sets = [self.available_objects[name] for name in set_names]
        ind = 0
        while ind < len(set_names):
            if sets[ind] is None:
                set_names.pop(ind)
                sets.pop(ind)
            else:
                ind += 1

        if len(set_names) < 2:
            canvas = EmptyCanvas('Please select 2 or more gene sets to continue', self)
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
                canvas = VennCanvas(items, self)
                self._connect_canvas(canvas)

            else:
                canvas = EmptyCanvas('UpSet not implemented yet!', self)
        if 'venn' in self.widgets:
            self.widgets['canvas'].deleteLater()
            self.widgets['toolbar'].deleteLater()
            self.operations_grid.removeWidget(self.widgets['canvas'])
            self.operations_grid.removeWidget(self.widgets['toolbar'])

        self.widgets['canvas'] = canvas
        self.widgets['toolbar'] = NavigationToolbar2QT(self.widgets['canvas'], self)
        self.operations_grid.addWidget(self.widgets['canvas'], 1, 3, 3, 3)
        self.operations_grid.addWidget(self.widgets['toolbar'], 0, 3, 1, 3)

    def _connect_canvas(self, canvas: VennCanvas):
        canvas.manualChoice.connect(self._set_op_other)
        self.widgets['radio_button_box'].radio_buttons['Union'].clicked.connect(canvas.union)
        self.widgets['radio_button_box'].radio_buttons['Intersection'].clicked.connect(canvas.intersection)
        self.widgets['radio_button_box'].radio_buttons['Symmetric Difference'].clicked.connect(
            canvas.symmetric_difference)
        self.primarySetChangedDifference.connect(canvas.difference)
        self.primarySetChangedIntersection.connect(canvas.intersection)

    def init_ui(self):
        self.list_group.setTitle('Gene sets')
        self.setWindowTitle('Set Operations')

        self.setLayout(self.layout)
        self.layout.addWidget(self.list_group)
        self.layout.addWidget(self.operations_group)
        self.operations_grid.addWidget(self.parameter_group, 1, 1, 3, 1)

        self.init_sets_ui()
        self.init_operations_ui()

    def init_sets_ui(self):
        self.widgets['set_list'] = QtWidgets.QListWidget(self)
        self.widgets['set_list'].insertItems(0, self.available_objects)
        self.widgets['set_list'].setSelectionMode(QtWidgets.QAbstractItemView.MultiSelection)

        for func in [self.create_canvas, self._check_legal_operations, self._validate_input,
                     self._toggle_choose_primary_set]:
            self.widgets['set_list'].itemSelectionChanged.connect(func)
        self.list_grid.addWidget(self.widgets['set_list'], 2, 1, 3, 1)

        self.widgets['select_all'] = QtWidgets.QPushButton('Select all', self)
        self.widgets['select_all'].clicked.connect(self.select_all)
        self.list_grid.addWidget(self.widgets['select_all'], 5, 1)

        self.widgets['clear_all'] = QtWidgets.QPushButton('Clear all', self)
        self.widgets['clear_all'].clicked.connect(self.clear_all)
        self.list_grid.addWidget(self.widgets['clear_all'], 6, 1)

        self.list_grid.addWidget(QtWidgets.QLabel('<b>Choose gene sets:</b>', self), 1, 1, 1, 1)

    def select_all(self):
        for ind in range(self.widgets['set_list'].count()):
            item = self.widgets['set_list'].item(ind)
            if not item.isSelected():
                item.setSelected(True)

    def clear_all(self):
        for item in self.widgets['set_list'].selectedItems():
            item.setSelected(False)

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
        self.widgets['choose_primary_set_label'] = QtWidgets.QLabel('Choose primary set for this operation:')
        self.widgets['set_op_box_layout'].addWidget(self.widgets['choose_primary_set_label'])
        self.widgets['set_op_box_layout'].addWidget(self.widgets['choose_primary_set'])

        self._toggle_choose_primary_set()

        self.create_canvas()

        self.widgets['apply_button'] = QtWidgets.QPushButton(text='Apply')
        self.widgets['apply_button'].clicked.connect(self.apply_set_op)
        self.widgets['apply_button'].setEnabled(False)
        self.operations_grid.addWidget(self.widgets['apply_button'], 4, 0, 1, 6)

    def _majority_vote_intersection(self):
        if not isinstance(self.widgets['canvas'], EmptyCanvas):
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
        self.clear_layout(self.parameter_grid)

        chosen_func_name = self.get_current_func_name()
        signature = generic.get_method_signature(chosen_func_name, filtering.Filter)
        i = 0
        for name, param in signature.items():
            if name in {'self', 'other', 'others', 'return_type'}:
                continue
            self.parameter_widgets[name] = gui_utils.param_to_widget(param, name)
            self.parameter_grid.addWidget(QtWidgets.QLabel(f'{name}:', self.parameter_widgets[name]), i, 0)
            self.parameter_grid.addWidget(self.parameter_widgets[name], i, 1)
            if chosen_func_name == 'majority_vote_intersection':
                self.parameter_widgets[name].valueChanged.connect(self._majority_vote_intersection)
                self._majority_vote_intersection()

        if chosen_func_name != 'other':
            help_address = f"https://guyteichman.github.io/RNAlysis/build/rnalysis.filtering." \
                           f"{filtering.Filter.__name__}.{chosen_func_name}.html"
            self.parameter_widgets['help_link'] = QtWidgets.QLabel(
                f'<a href="{help_address}">Open documentation for function '
                f'<b>{filtering.Filter.__name__}.{chosen_func_name}</b></a>', self)
            self.parameter_widgets['help_link'].setOpenExternalLinks(True)
            self.operations_grid.addWidget(self.parameter_widgets['help_link'], 5, 0, 1, 2)

    def get_current_func_name(self):
        button = self.widgets['radio_button_box'].checkedButton()
        if button is None:
            return None
        return self.SET_OPERATIONS[button.text()]

    @staticmethod
    def clear_layout(layout):
        while layout.count():
            child = layout.takeAt(0)
            if child.widget():
                child.widget().deleteLater()

    def _check_legal_operations(self):
        n_items = len(self.widgets['set_list'].selectedItems())
        if self.get_current_func_name() == 'symmetric_difference' and n_items > 2:
            self.widgets['radio_button_box'].set_selection('Other')
        sym_diff_button = self.widgets['radio_button_box'].radio_buttons['Symmetric Difference']
        sym_diff_button.setEnabled(n_items <= 2)

    def _validate_input(self):
        is_legal = True

        if isinstance(self.widgets['canvas'], EmptyCanvas):
            is_legal = False

        if self.widgets['radio_button_box'].checkedButton() is None:
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

            first_obj = self.available_objects[primary_set_name]
            if isinstance(first_obj, set):
                first_obj = filtering.Filter(('placeholder', pd.DataFrame(index=first_obj)))
            other_objs = []
            for name in set_names:
                if name != primary_set_name:
                    other_objs.append(self.available_objects[name])
            output_set = getattr(first_obj, func_name)(*other_objs, **kwargs)
            output_name = f"{func_name} output"
        if isinstance(output_set, set):
            self.geneSetReturned.emit(output_set, output_name)
        self.close()


class TabPage(QtWidgets.QWidget):
    def __init__(self, parent=None):
        super().__init__(parent)
        self.layout = QtWidgets.QVBoxLayout(self)
        self.name = None

        self.stdout_group = QtWidgets.QGroupBox('Log')
        self.layout.addWidget(self.stdout_group)
        self.stdout_grid = QtWidgets.QGridLayout(self.stdout_group)
        self.stdout_widgets = {}

        self.init_stdout_ui()

    def init_stdout_ui(self):
        self.stdout_widgets['text_edit_stdout'] = gui_utils.StdOutTextEdit(self)
        self.stdout_widgets['text_edit_stdout'].setStyleSheet("""QTextEdit {background: #dddddd;}""")
        # self.stdout_widgets['progress_bar'] = gui_utils.TQDMLineEdit(self)

        self.stdout_grid.addWidget(self.stdout_widgets['text_edit_stdout'], 0, 0, 3, 4)
        # self.stdout_grid.addWidget(self.stdout_widgets['progress_bar'], 3, 0, 1, 4)

    def get_console(self):
        return self.stdout_widgets['text_edit_stdout']

    def get_pbar(self):
        return self.stdout_widgets['progress_bar']

    def rename(self):
        new_name = self.overview_widgets['table_name'].text()
        self.set_tab_name(new_name, is_unsaved=True)
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

    def set_tab_name(self, new_name: str, is_unsaved: bool):
        if is_unsaved:
            new_name += '*'
        parent = self._get_parent_tabwidget()
        parent.setTabText(parent.currentIndex(), new_name)

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
            self.set_tab_name(self.get_tab_name().strip("*"), is_unsaved=False)

    def rename(self):
        new_name = self.overview_widgets['table_name'].text()
        super().rename()
        self.gene_set.change_set_name(new_name)

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


class FilterTabPage(TabPage):
    FILTER_OBJ_TYPES = {'Count matrix': filtering.CountFilter, 'Differential expression': filtering.DESeqFilter,
                        'Fold change': filtering.FoldChangeFilter, 'Other': filtering.Filter}
    EXCLUDED_FUNCS = {'union', 'intersection', 'majority_vote_intersection', 'difference', 'symmetric_difference'}

    def __init__(self, parent=None):
        super().__init__(parent)
        self.filter_obj = None
        self.df_views = []

        self.basic_group = QtWidgets.QGroupBox('Load a data table')
        self.layout.insertWidget(0, self.basic_group)
        self.basic_grid = QtWidgets.QGridLayout(self.basic_group)
        self.basic_widgets = {}

        self.overview_group = QtWidgets.QGroupBox('Data overview')
        self.overview_grid = QtWidgets.QGridLayout(self.overview_group)
        self.overview_widgets = {}

        self.parameter_group = QtWidgets.QGroupBox('Apply functions')
        self.parameter_grid = QtWidgets.QGridLayout(self.parameter_group)
        self.parameter_widgets = {}

        self.init_basic_ui()

        # self.grid.addWidget(QtWidgets.QSpinBox(), 2, 0)
        # self.grid.addWidget(QtWidgets.QPushButton('Clear Text'), 2, 2)

    @staticmethod
    def clear_layout(layout, exceptions: set = frozenset()):
        while layout.count() > len(exceptions):
            child = layout.takeAt(0)
            if child.widget() and child.widget() not in exceptions:
                child.widget().deleteLater()

    def rename(self):
        new_name = self.overview_widgets['table_name'].text()
        super().rename()
        self.filter_obj.fname = Path(
            os.path.join(str(self.filter_obj.fname.parent), f"{new_name}{self.filter_obj.fname.suffix}"))
        print(self.filter_obj.fname)

    def init_overview_ui(self):
        this_row = 0
        self.layout.insertWidget(1, self.overview_group)
        self.overview_widgets['table_type_label'] = QtWidgets.QLabel(
            f"Table type: {self.basic_widgets['table_type_combo'].currentText()}")
        self.overview_widgets['table_name_label'] = QtWidgets.QLabel(f"Table name: '<b>{self.get_tab_name()}</b>'")

        self.overview_widgets['preview'] = QtWidgets.QTableView()
        # self.overview_widgets['preview'].setFixedWidth(550)
        # self.overview_widgets['preview'].setFixedHeight(175)
        self.overview_widgets['preview'].setVerticalScrollBarPolicy(QtCore.Qt.ScrollBarAlwaysOff)
        self.overview_widgets['preview'].setHorizontalScrollBarPolicy(QtCore.Qt.ScrollBarAlwaysOff)
        # self.overview_widgets['preview'].setFlags(self.overview_widgets['preview'].flags() & ~QtCore.Qt.ItemIsEditable)
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
        self.overview_grid.addWidget(self.overview_widgets['preview'], this_row, 0, 3, 4)
        this_row += 3

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
        self.basic_widgets['table_type_combo'] = QtWidgets.QComboBox()
        self.basic_widgets['table_type_combo'].addItems(self.FILTER_OBJ_TYPES.keys())

        self.basic_widgets['start_button'] = QtWidgets.QPushButton('Start')
        self.basic_widgets['start_button'].clicked.connect(self.start)
        self.basic_widgets['start_button'].setEnabled(False)

        self.basic_widgets['file_path'] = gui_utils.PathLineEdit()
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
        self.basic_grid.addWidget(self.basic_widgets['start_button'], 2, 0, 1, 4)

    def init_action_ui(self):
        self.layout.insertWidget(2, self.parameter_group)
        self.parameter_widgets['function_combo'] = QtWidgets.QComboBox()
        self.parameter_widgets['function_combo'].addItem('Choose a function...')
        self.parameter_widgets['function_combo'].addItems(self.get_all_actions())
        self.parameter_widgets['function_combo'].currentTextChanged.connect(self.update_parameter_ui)

        self.parameter_grid.addWidget(self.parameter_widgets['function_combo'], 0, 0, 1, 2)

    def update_parameter_ui(self):
        # delete previous widgets
        self.parameter_grid.removeWidget(self.parameter_widgets['function_combo'])
        self.clear_layout(self.parameter_grid)
        self.parameter_widgets = {'function_combo': self.parameter_widgets['function_combo']}
        self.parameter_grid.addWidget(self.parameter_widgets['function_combo'], 0, 0, 1, 2)

        chosen_func_name = self.parameter_widgets['function_combo'].currentText()
        signature = generic.get_method_signature(chosen_func_name, self.filter_obj)
        i = 1
        for name, param in signature.items():
            if name in {'self'}:
                continue
            self.parameter_widgets[name] = gui_utils.param_to_widget(param, name)
            self.parameter_grid.addWidget(QtWidgets.QLabel(f'{name}:', self.parameter_widgets[name]), i, 0, )
            self.parameter_grid.addWidget(self.parameter_widgets[name], i, 1)
            i += 1

        help_address = f"https://guyteichman.github.io/RNAlysis/build/rnalysis.filtering." \
                       f"{type(self.filter_obj).__name__}.{chosen_func_name}.html"
        self.parameter_widgets['help_link'] = QtWidgets.QLabel(
            text=f'<a href="{help_address}">Open documentation for function '
                 f'<b>{type(self.filter_obj).__name__}.{chosen_func_name}</b></a>')
        self.parameter_widgets['help_link'].setOpenExternalLinks(True)
        self.parameter_grid.addWidget(self.parameter_widgets['help_link'], i + 1, 0, 1, 2)

    def view_full_dataframe(self):
        df_window = gui_utils.DataFrameView(self.filter_obj.df, self.filter_obj.fname)
        self.overview_widgets['full_table_view'] = df_window
        df_window.show()

    def update_table_preview(self):
        model = gui_utils.DataFramePreviewModel(self.filter_obj.df)
        self.overview_widgets['preview'].setModel(model)

    def _get_function_params(self):
        func_params = {}
        for param_name, widget in self.parameter_widgets.items():
            if param_name in {'apply_button', 'help_link'}:
                continue
            val = gui_utils.get_val_from_widget(widget)

            func_params[param_name] = val
        return func_params

    def apply_function(self):
        func_name = self.parameter_widgets['function_combo'].currentText()
        func_params = self._get_function_params()
        prev_name = self.get_tab_name()
        result = getattr(self.filter_obj, func_name)(**func_params)
        self.update_filter_obj_shape()
        self.update_table_preview()
        # TODO: update table name!
        if prev_name != self.filter_obj.fname.name:
            self.set_tab_name(self.filter_obj.fname.name, is_unsaved=True)

        self._proccess_outputs(result, func_name)

    def _proccess_outputs(self, outputs, source_name: str = ''):
        if validation.isinstanceinh(outputs, filtering.Filter):
            self._get_parent_window().new_tab_from_filter_obj(outputs)
        elif isinstance(outputs, pd.DataFrame):
            self.df_views.append(gui_utils.DataFrameView(outputs, source_name))
            self.df_views[-1].show()
        elif isinstance(outputs, np.ndarray):
            df = pd.DataFrame(outputs)
            self._proccess_outputs(df, source_name)

        elif isinstance(outputs, (tuple, list)):
            for output in outputs:
                self._proccess_outputs(output, source_name)
        elif isinstance(outputs, dict):
            for output, this_src_name in outputs.values():
                self._proccess_outputs(output, this_src_name)

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
            self.set_tab_name(self.filter_obj.fname.name + '*')

        self._proccess_outputs(result, 'pipeline')

    def get_index_string(self):
        return self.filter_obj.index_string

    def get_all_actions(self):
        assert self.filter_obj is not None, "No table was loaded!"
        all_methods = dir(self.filter_obj)
        public_methods = [mthd for mthd in all_methods if
                          (not mthd.startswith('_')) and (callable(getattr(type(self.filter_obj), mthd))) and (
                              mthd not in self.EXCLUDED_FUNCS)]
        return public_methods

    # def choose_file(self):
    #     filename, _ = QtWidgets.QFileDialog.getOpenFileName(self, "Choose a file", str(Path.home()),
    #                                                         "Comma-Separated Values (*.csv);;"
    #                                                         "Tab-Separated Values (*.tsv);;"
    #                                                         "All Files (*)")
    #     if filename:
    #         self.basic_widgets['file_path'].setText(filename)

    def save_file(self):
        default_name = str(self.filter_obj.fname).rstrip("*")
        filename, _ = QtWidgets.QFileDialog.getSaveFileName(self, "Save filtering result",
                                                            str(Path.home().joinpath(default_name)),
                                                            "Comma-Separated Values (*.csv);;"
                                                            "Tab-Separated Values (*.tsv);;"
                                                            "All Files (*)")
        if filename:
            self.filter_obj.save_csv(filename)
            print(f"Successfully saved at {io.get_datetime()} under {filename}")
            self.set_tab_name(self.filter_obj.fname.name, is_unsaved=False)

    def start(self):
        self.filter_obj = self.FILTER_OBJ_TYPES[self.basic_widgets['table_type_combo'].currentText()](
            self.basic_widgets['file_path'].text())
        print(self.filter_obj)
        table_name_user_input = self.basic_widgets['table_name'].text()
        if table_name_user_input != '':
            new_name = table_name_user_input
        else:
            new_name = self.filter_obj.fname.name
        self.set_tab_name(new_name, is_unsaved=False)

        self.init_overview_ui()
        self.init_action_ui()

        self.clear_layout(self.basic_grid)
        self.layout.removeWidget(self.basic_group)
        self.basic_group.deleteLater()

    def start_from_filter_obj(self, filter_obj: filtering.Filter):
        self.filter_obj = filter_obj
        print(self.filter_obj)
        self.init_overview_ui()
        self.init_action_ui()


class CreatePipelineWindow(gui_utils.MinMaxDialog, FilterTabPage):
    def __init__(self, parent=None):
        super().__init__(parent)
        self.setLayout(self.layout)
        self.setWindowTitle(f'Create new Pipeline')
        self.setGeometry(500, 200, 500, 300)
        self.pipeline = None

    def init_basic_ui(self):
        self.basic_group.setTitle("Choose data table type for Pipeline")

        self.basic_widgets['table_type_combo'] = QtWidgets.QComboBox()
        self.basic_widgets['table_type_combo'].addItems(self.FILTER_OBJ_TYPES.keys())

        self.basic_widgets['pipeline_name'] = QtWidgets.QLineEdit()
        self.basic_widgets['pipeline_name'].setText('New Pipeline')

        self.basic_widgets['name_label'] = QtWidgets.QLabel('Name your Pipeline:')

        self.basic_widgets['start_button'] = QtWidgets.QPushButton('start_button')
        self.basic_widgets['start_button'].clicked.connect(self.start)

        self.basic_widgets['type_label'] = QtWidgets.QLabel('Choose table type:')

        self.basic_grid.addWidget(self.basic_widgets['pipeline_name'], 1, 1)
        self.basic_grid.addWidget(self.basic_widgets['name_label'], 0, 0, 1, 2)
        self.basic_grid.addWidget(self.basic_widgets['type_label'], 0, 2)
        self.basic_grid.addWidget(self.basic_widgets['table_type_combo'], 1, 2)
        self.basic_grid.addWidget(self.basic_widgets['start_button'], 1, 3)

    def init_action_ui(self):
        super().init_action_ui()
        self.parameter_group.setTitle("Add functions to Pipeline")

    def update_parameter_ui(self):
        super().update_parameter_ui()
        self.parameter_widgets['apply_button'].setText('Add to Pipeline')

    def apply_function(self):
        func_name = self.action_widgets['function_combo'].currentText()
        func_params = self._get_function_params()
        self.pipeline.add_function(func_name, **func_params)
        self.update_pipeline_preview()

    def update_pipeline_preview(self):
        self.overview_widgets['preview'].setPlainText(str(self.pipeline))

    def update_table_preview(self):
        raise NotImplementedError

    def update_filter_obj_shape(self):
        raise NotImplementedError

    def set_file_path_bg_color(self):
        raise NotImplementedError

    def start(self):
        filt_obj_type = self.FILTER_OBJ_TYPES[self.basic_widgets['table_type_combo'].currentText()]
        self.filter_obj = filt_obj_type.__new__(filt_obj_type)
        self.pipeline = filtering.Pipeline(filt_obj_type)
        self.init_overview_ui()
        self.init_action_ui()

    def init_overview_ui(self):
        self.parameter_group.setTitle("Pipeline preview")
        self.layout.insertWidget(1, self.overview_group)
        self.overview_widgets['preview'] = QtWidgets.QPlainTextEdit()
        self.overview_widgets['preview'].setReadOnly(True)
        self.update_pipeline_preview()

        self.overview_grid.addWidget(QtWidgets.QLabel('Pipeline preview', self.overview_widgets['preview']), 0, 0, 2, 4)
        self.overview_grid.addWidget(self.overview_widgets['preview'], 1, 0, 2, 4)

        self.overview_widgets['save_button'] = QtWidgets.QPushButton('Save Pipeline')
        self.overview_widgets['save_button'].clicked.connect(self.save_file)
        self.overview_grid.addWidget(self.overview_widgets['save_button'], 3, 3)

        # self.overview_widgets['shape'] = QtWidgets.QLabel()
        # self.overview_grid.addWidget(self.overview_widgets['shape'], 3, 0, 1, 2)

        self.overview_widgets['export_button'] = QtWidgets.QPushButton('Export Pipeline')
        # self.overview_widgets['export_button'].clicked.connect(self.export_pipeline)
        self.overview_grid.addWidget(self.overview_widgets['export_button'], 3, 2, 1, 1)

    def set_tab_name(self, new_name: str, is_unsaved: bool):
        raise NotImplementedError

    def save_file(self):
        self._get_parent_window().pipelines[self.basic_widgets['pipeline_name'].text()] = self.pipeline
        print(f"Successfully saved Pipeline '{self.basic_widgets['pipeline_name'].text()}'")


class MainWindow(QtWidgets.QMainWindow):
    USER_GUIDE_URL = 'https://guyteichman.github.io/RNAlysis/build/user_guide.html'

    def __init__(self):
        super().__init__()
        self.setWindowTitle(f'RNAlysis {__version__}')
        self.setWindowIcon(QtGui.QIcon('../../docs/source/favicon.ico'))
        self.setGeometry(600, 50, 750, 350)
        self.update_style_sheet()

        self.tabs = QtWidgets.QTabWidget()
        self.next_tab_n = 0
        self.tabs.setTabsClosable(True)
        self.tabs.tabCloseRequested.connect(self.delete)

        self.layout = QtWidgets.QVBoxLayout(self)
        self.layout.addWidget(self.tabs)

        self.add_tab_button = QtWidgets.QToolButton()
        self.add_tab_button.setToolTip('Add New Tab')
        self.add_tab_button.clicked.connect(functools.partial(self.add_new_tab, name=None))
        self.add_tab_button.setText("+")
        self.error_window = None

        self.tabs.setCornerWidget(self.add_tab_button, QtCore.Qt.TopRightCorner)
        self.setCentralWidget(self.tabs)

        self.menu_bar = QtWidgets.QMenuBar(self)

        self.pipelines = OrderedDict()
        self.pipeline_window = CreatePipelineWindow(self)

        self.about_window = gui_utils.AboutWindow(self)
        self.settings_window = gui_utils.SettingsWindow(self)
        self.set_op_window = None
        self.enrichment_window = None
        self.enrichment_results = []

        self.add_new_tab()
        self.init_actions()
        self.init_menu_ui()

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

        # self.queue_tqdm = Queue()
        # self.thread_tqdm_queue_listener = QtCore.QThread()
        # self.tqdm_receiver = gui_utils.ThreadTQDMStreamTextQueueReceiver(self.queue_tqdm)
        # sys.stderr = gui_utils.WriteStream(self.queue_tqdm)
        #
        # # connect receiver object to widget for text update
        # self.tqdm_receiver.queue_tqdm_element_received_signal.connect(self.append_text_to_current_pbar)
        # # attach console text receiver to console text thread
        # self.tqdm_receiver.moveToThread(self.thread_tqdm_queue_listener)
        # # attach to start / stop methods
        # self.thread_tqdm_queue_listener.started.connect(self.tqdm_receiver.run)
        # self.thread_tqdm_queue_listener.start()

    def update_style_sheet(self):
        self.setStyleSheet(gui_style.get_stylesheet())

    def add_new_tab(self, name: str = None, is_set: bool = False):
        if name is None:
            name = f'New Table {self.next_tab_n + 1}'
            self.next_tab_n += 1

        else:
            print(name)

        if is_set:
            self.tabs.addTab(SetTabPage(name, parent=self.tabs), name)
        else:
            self.tabs.addTab(FilterTabPage(parent=self.tabs), name)
        self.tabs.setCurrentIndex(self.tabs.count() - 1)

    def new_tab_from_filter_obj(self, filter_obj: filtering.Filter):
        self.add_new_tab(filter_obj.fname.name)

        self.tabs.currentWidget().start_from_filter_obj(filter_obj)

    def new_tab_from_gene_set(self, gene_set: set, gene_set_name: str):
        self.add_new_tab(gene_set_name, is_set=True)
        self.tabs.currentWidget().update_gene_set(gene_set)

    def import_gene_set(self):
        filename, _ = QtWidgets.QFileDialog.getOpenFileName(self, "Choose a file", str(Path.home()),
                                                            "Text Document (*.txt);;"
                                                            "Comma-Separated Values (*.csv);;"
                                                            "Tab-Separated Values (*.tsv);;"
                                                            "All Files (*)")
        if filename:
            if filename.endswith('.csv'):
                gene_set = set(pd.read_csv(filename, index_col=0).index)
            elif filename.endswith('.tsv'):
                gene_set = set(pd.read_csv(filename, index_col=0, sep='\t').index)
            else:
                with open(filename) as f:
                    gene_set = {line.strip() for line in f.readlines()}

            self.new_tab_from_gene_set(gene_set, filename)

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
        self.pipeline_window.exec()

    def settings(self):
        self.settings_window.exec()

    def init_actions(self):
        self.new_action = QtWidgets.QAction("&New...", self)
        self.new_action.triggered.connect(functools.partial(self.add_new_tab, name=None))
        self.save_action = QtWidgets.QAction("&Save...", self)
        self.save_action.triggered.connect(self.tabs.currentWidget().save_file)
        self.settings_action = QtWidgets.QAction("&Settings...", self)
        self.settings_action.triggered.connect(self.settings)
        self.exit_action = QtWidgets.QAction("&Exit", self)
        self.exit_action.triggered.connect(self.closeEvent)

        self.copy_action = QtWidgets.QAction("&Copy Gene Set", self)
        self.copy_action.triggered.connect(self.copy_gene_set)
        self.set_op_action = QtWidgets.QAction("Set &Operations...")
        self.set_op_action.triggered.connect(self.choose_set_op)
        self.enrichment_action = QtWidgets.QAction("E&nrichment Analysis...")
        self.enrichment_action.triggered.connect(self.open_enrichment_analysis)
        self.set_viz_action = QtWidgets.QAction("Gene Set &Visualization...")
        self.import_set_action = QtWidgets.QAction("&Import Gene Set...", self)
        self.import_set_action.triggered.connect(self.import_gene_set)
        self.export_set_action = QtWidgets.QAction("&Export Gene Set...", self)

        self.user_guide_action = QtWidgets.QAction("&User Guide", self)
        self.user_guide_action.triggered.connect(self.open_user_guide)
        self.help_action = QtWidgets.QAction("&Help", self)
        self.about_action = QtWidgets.QAction("&About", self)
        self.about_action.triggered.connect(self.about)

        self.new_pipeline_action = QtWidgets.QAction("&New Pipeline...", self)
        self.new_pipeline_action.triggered.connect(self.add_pipeline)
        self.import_pipeline_action = QtWidgets.QAction("&Import Pipeline...", self)
        self.export_pipeline_action = QtWidgets.QAction("&Export Pipeline...", self)

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
                available_objects_unique[f"{name}_{checked[name]}"] = self.get_gene_set_by_ind(i)
            else:
                checked[name] = 1
                available_objects_unique[name] = self.get_gene_set_by_ind(i)
            checked[name] += 1
        return available_objects_unique

    def get_gene_set_by_name(self, name: str):
        ind = self.get_tab_names().index(name)
        return self.get_gene_set_by_ind(ind)

    def choose_set_op(self):
        available_objs = self.get_available_objects()
        self.set_op_window = SimpleSetOpWindow(available_objs, self)
        self.set_op_window.primarySetUsed.connect(self.choose_tab_by_name)
        self.set_op_window.geneSetReturned.connect(self.new_tab_from_gene_set)
        self.set_op_window.show()

    @QtCore.pyqtSlot(str)
    def choose_tab_by_name(self, set_name: str):
        available_objs = self.get_available_objects()
        for i, name in enumerate(available_objs):
            if name == set_name:
                self.tabs.setCurrentIndex(i)
                return

    def run_enrichment_analysis(self, func: Callable, gene_set_ind: int, bg_set_ind: Union[int, None], kwargs: dict):
        is_single_set = bg_set_ind is None

        if is_single_set:
            bg_set_obj = None
        else:
            self.tabs.setCurrentIndex(bg_set_ind)
            bg_set_obj = enrichment.FeatureSet(self.tabs.currentWidget().filter_obj) if \
                isinstance(self.tabs.currentWidget(), FilterTabPage) else self.tabs.currentWidget().gene_set

        self.tabs.setCurrentIndex(gene_set_ind)
        set_name = self.tabs.currentWidget().get_tab_name()
        gene_set = self.tabs.currentWidget().filter_obj if \
            isinstance(self.tabs.currentWidget(),
                       FilterTabPage) else self.tabs.currentWidget().gene_set.gene_set

        feature_set_obj = enrichment.RankedSet(gene_set, set_name) if is_single_set \
            else enrichment.FeatureSet(gene_set, set_name)

        if is_single_set:
            result = func(feature_set_obj, **kwargs)
        else:
            result = func(feature_set_obj, background_genes=bg_set_obj, **kwargs)

        df_window = gui_utils.DataFrameView(result, "Enrichment results for set " + feature_set_obj.set_name)
        self.enrichment_results.append(df_window)
        df_window.show()

    def open_enrichment_analysis(self):
        tab_names = self.get_tab_names()
        self.enrichment_window = EnrichmentWindow(tab_names, self)
        self.enrichment_window.show()

    def get_tab_names(self) -> List[str]:
        return [self.tabs.tabText(i) for i in range(self.tabs.count())]

    def init_menu_ui(self):
        self.setMenuBar(self.menu_bar)
        file_menu = self.menu_bar.addMenu("&File")
        file_menu.addActions(
            [self.new_action, self.save_action, self.settings_action, self.exit_action])

        gene_sets_menu = self.menu_bar.addMenu("&Gene sets")
        gene_sets_menu.addActions(
            [self.copy_action, self.set_op_action, self.enrichment_action, self.set_viz_action, self.import_set_action,
             self.export_set_action])

        pipeline_menu = self.menu_bar.addMenu("&Pipelines")
        pipeline_menu.addActions([self.new_pipeline_action, self.import_pipeline_action, self.export_pipeline_action])
        self.apply_pipeline_menu = pipeline_menu.addMenu("&Apply Pipeline")
        self.apply_pipeline_menu.aboutToShow.connect(self._populate_pipelines)

        help_menu = self.menu_bar.addMenu("&Help")
        help_menu.addActions([self.help_action, self.user_guide_action, self.about_action])

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


def run():
    matplotlib.use('Qt5Agg')
    app = QtWidgets.QApplication(sys.argv)
    window = MainWindow()
    sys.excepthook = window.excepthook
    window.show()
    sys.exit(app.exec())


if __name__ == '__main__':
    run()
