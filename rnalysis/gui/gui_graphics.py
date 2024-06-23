import abc
import functools
import itertools
from pathlib import Path
from typing import List, Tuple, Callable
import pandas as pd
import matplotlib
import matplotlib_venn
import upsetplot
from PyQt5 import QtCore, QtGui
from matplotlib import pyplot as plt
from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg, NavigationToolbar2QT

from rnalysis.utils import parsing, generic

COLOR_ICONS = {'black': 'black_icon.png', 'green': 'green_icon.png', 'red': 'red_icon.png', 'white': 'white_icon.png',
               'blue': 'blue_icon.png', 'grey': 'grey_icon.png', 'pink': 'pink_icon.png', 'orange': 'orange_icon.png',
               'yellow': 'yellow_icon.png'}
AVAILABLE_ICONS = {'CountFilter': 'countfilter_icon.png', 'DESeqFilter': 'deseqfilter_icon.png',
                   'set': 'featureset_icon.png', 'Filter': 'filter_icon.png',
                   'FoldChangeFilter': 'foldchangefilter_icon.png', }
AVAILABLE_ICONS.update(COLOR_ICONS)


class CleanPlotToolBar(NavigationToolbar2QT):
    toolitems = [t for t in NavigationToolbar2QT.toolitems if
                 t[0] in ('Home', 'Pan', 'Zoom')]

    def set_message(self, msg):
        pass


class BasePreviewCanvas(FigureCanvasQTAgg):
    def __init__(self, plotting_func: Callable, parent=None, tight_layout: bool = True, **plotting_kwargs):
        self.fig = plt.Figure(tight_layout=tight_layout)
        self.parent = parent
        self.generated_plot = plotting_func(**plotting_kwargs, fig=self.fig)
        super().__init__(figure=self.fig)


class BaseInteractiveCanvas(FigureCanvasQTAgg):
    DESELECTED_STATE = 0
    SELECTED_STATE = 1
    HOVER_STATE = 2
    HOVER_SELECTED_STATE = 3

    DESELECTED_COLOR = (0, 0, 0)
    SELECTED_COLOR = (0.55, 0.1, 0.8)
    HOVER_COLOR = (0.7, 0.7, 0.7)
    HOVER_SELECTED_COLOR = (0.7, 0.3, 1)

    COLORMAP = {DESELECTED_STATE: DESELECTED_COLOR, SELECTED_STATE: SELECTED_COLOR, HOVER_STATE: HOVER_COLOR,
                HOVER_SELECTED_STATE: HOVER_SELECTED_COLOR}

    manualChoice = QtCore.pyqtSignal()

    def __init__(self, gene_sets: dict, parent=None, constrained_layout: bool = False):
        self.parent = parent
        self.gene_sets = gene_sets
        self.fig = plt.Figure(constrained_layout=constrained_layout)
        super().__init__(self.fig)

        self.fig.canvas.mpl_connect('button_press_event', self.on_click)
        self.fig.canvas.mpl_connect('motion_notify_event', self.on_hover)

    @abc.abstractmethod
    def on_click(self, event):
        raise NotImplementedError

    @abc.abstractmethod
    def on_hover(self, event):
        raise NotImplementedError

    @abc.abstractmethod
    def union(self):
        raise NotImplementedError

    @abc.abstractmethod
    def intersection(self):
        raise NotImplementedError

    @abc.abstractmethod
    def difference(self, primary_set: str):
        raise NotImplementedError

    @abc.abstractmethod
    def majority_vote_intersection(self, majority_threshold: float):
        raise NotImplementedError

    @abc.abstractmethod
    def clear_selection(self):
        raise NotImplementedError

    @abc.abstractmethod
    def select(self, ind, draw: bool = True):
        raise NotImplementedError

    @abc.abstractmethod
    def deselect(self, ind, draw: bool = True):
        raise NotImplementedError

    @abc.abstractmethod
    def get_custom_selection(self) -> set:
        raise NotImplementedError

    @abc.abstractmethod
    def get_tuple_patch_ids(self) -> List[Tuple[int, ...]]:
        raise NotImplementedError

    def draw(self):
        self.fig.suptitle(f"{len(self.get_custom_selection())} selected genes")
        super().draw()


class VennInteractiveCanvas(BaseInteractiveCanvas):
    def __init__(self, gene_sets: dict, parent=None):
        super().__init__(gene_sets, parent)
        self.ax = self.fig.add_subplot()

        if len(gene_sets) == 2:
            funcs = matplotlib_venn.venn2, matplotlib_venn.venn2_circles
            colors = (self.SELECTED_COLOR, self.SELECTED_COLOR)

        elif len(gene_sets) == 3:
            funcs = matplotlib_venn.venn3, matplotlib_venn.venn3_circles
            colors = (self.SELECTED_COLOR, self.SELECTED_COLOR, self.SELECTED_COLOR)
        else:
            raise ValueError("Cannot process more than 3 sets!")

        self.venn = funcs[0](gene_sets.values(), gene_sets.keys(), set_colors=colors, ax=self.ax, alpha=1)
        self.venn_circles = funcs[1](gene_sets.values(), linestyle='solid', linewidth=2.0, ax=self.ax)
        self.default_subset_fontsize = 14
        self.states = [self.DESELECTED_STATE for _ in range(len(self.venn.patches))]
        self.set_font_size(16, self.default_subset_fontsize)
        self.clear_selection()

    def set_font_size(self, set_label_size: float, subset_label_size: float):
        for label in self.venn.set_labels:
            if label is not None:
                label.set_fontsize(set_label_size)
        for label in self.venn.subset_labels:
            if label is not None:
                label.set_fontsize(subset_label_size)

    def clear_selection(self, draw: bool = True):
        for subset in range(len(self.states)):
            patch = self.venn.patches[subset]
            if patch is None:
                continue
            self.update_color(subset, self.DESELECTED_STATE)
        if draw:
            self.draw()

    def select(self, ind, draw: bool = True):
        patch = self.venn.get_patch_by_id(ind)
        subset = self.venn.patches.index(patch)
        if patch is None:
            return
        self.update_color(subset, self.SELECTED_STATE)
        if draw:
            self.draw()

    def deselect(self, ind, draw: bool = True):
        patch = self.venn.get_patch_by_id(ind)
        subset = self.venn.patches.index(patch)
        if patch is None:
            return
        self.update_color(subset, self.DESELECTED_STATE)
        if draw:
            self.draw()

    def update_color(self, subset: int, state: int):
        color = self.COLORMAP[state]
        font_color = 'orange' if state != self.DESELECTED_STATE else 'black'
        fontweight = 'bold' if state != self.DESELECTED_STATE else 'regular'
        fontsize = 16 if state != self.DESELECTED_STATE else self.default_subset_fontsize

        self.venn.patches[subset].set_facecolor(color)
        self.venn.patches[subset].set_fill(state != self.DESELECTED_STATE)

        self.venn.subset_labels[subset].set_fontsize(fontsize)
        self.venn.subset_labels[subset].set_color(font_color)
        self.venn.subset_labels[subset].set_fontweight(fontweight)

        self.states[subset] = state

    def on_click(self, event):
        for subset in range(len(self.states)):
            patch = self.venn.patches[subset]
            if patch is None:
                continue

            if patch.contains_point((event.x, event.y)):
                self.manualChoice.emit()
                if self.states[subset] in [self.DESELECTED_STATE, self.HOVER_STATE]:
                    self.update_color(subset, self.HOVER_SELECTED_STATE)
                else:
                    self.update_color(subset, self.HOVER_STATE)
        self.draw()

    def on_hover(self, event):
        for subset in range(len(self.states)):
            patch = self.venn.patches[subset]
            if patch is None:
                continue

            if patch.contains_point((event.x, event.y)):
                if self.states[subset] in [self.HOVER_STATE, self.DESELECTED_STATE]:
                    self.update_color(subset, self.HOVER_STATE)
                else:
                    self.update_color(subset, self.HOVER_SELECTED_STATE)
            else:
                if self.states[subset] in [self.HOVER_STATE, self.DESELECTED_STATE]:
                    self.update_color(subset, self.DESELECTED_STATE)
                else:
                    self.update_color(subset, self.SELECTED_STATE)
        self.draw()

    def get_tuple_patch_ids(self) -> List[Tuple[int, ...]]:
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

    @QtCore.pyqtSlot(float)
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

    def get_custom_selection(self) -> set:
        selection = set()
        for patch_id in self.get_tuple_patch_ids():
            str_patch_id = ''.join(str(i) for i in patch_id)
            patch = self.venn.get_patch_by_id(str_patch_id)
            subset = self.venn.patches.index(patch)
            if patch is None:
                continue
            if self.states[subset] in [self.SELECTED_STATE, self.HOVER_SELECTED_STATE]:
                included_sets = [s for s, ind in zip(self.gene_sets.values(), patch_id) if ind]
                excluded_sets = [s for s, ind in zip(self.gene_sets.values(), patch_id) if not ind]
                patch_content = set.intersection(*included_sets).difference(*excluded_sets)
                selection.update(patch_content)
        return selection


class UpSetInteractiveCanvas(BaseInteractiveCanvas):
    def __init__(self, gene_sets: dict, parent=None):
        super().__init__(gene_sets, parent, constrained_layout=False)
        self.upset_df = parsing.generate_upset_series(gene_sets)
        with pd.option_context("mode.copy_on_write", False):
            self.upset = upsetplot.UpSet(self.upset_df, sort_by='degree', sort_categories_by=None)
            self.subset_states = {i: self.DESELECTED_STATE for i in range(len(self.upset.subset_styles))}

            self.axes = self.upset.plot(self.fig)
            for ax in self.axes.values():
                generic.despine(ax)

            self.bounding_boxes = []
            axis_ymax = self.axes['intersections'].dataLim.bounds[3]
            for patch in self.axes['intersections'].patches:
                x, y = patch.get_xy()
                width = patch.get_width()
                bbox = matplotlib.patches.Rectangle((x, y), width, axis_ymax, visible=False)
                self.bounding_boxes.append(bbox)

            for bbox in self.bounding_boxes:
                self.axes['intersections'].add_patch(bbox)

    def draw(self):
        # update matrix
        with pd.option_context("mode.copy_on_write", False):
            matrix_ax = self.axes['matrix']
            matrix_ax.clear()
            self.upset.plot_matrix(matrix_ax)
            super().draw()

    def on_click(self, event):
        graph_modified = False
        for subset in self.subset_states:
            bounding_box = self.bounding_boxes[subset]

            if bounding_box.contains_point((event.x, event.y)):
                graph_modified = True
                self.manualChoice.emit()
                if self.subset_states[subset] in [self.DESELECTED_STATE, self.HOVER_STATE]:
                    _ = self.update_color(subset, self.HOVER_SELECTED_STATE)

                else:
                    _ = self.update_color(subset, self.HOVER_STATE)
        if graph_modified:
            self.draw()

    def update_color(self, subset: int, state: int) -> bool:
        color = self.COLORMAP[state]
        self.upset.subset_styles[subset]['facecolor'] = color
        self.axes['intersections'].patches[subset].set_facecolor(color)
        if self.subset_states[subset] != state:
            self.subset_states[subset] = state
            graph_modified = True
        else:
            graph_modified = False

        return graph_modified

    @staticmethod
    def _compare_ids(id1: Tuple[int, ...], id2: Tuple[int, ...]):
        if sum(id1) > sum(id2):
            return 1
        elif sum(id1) < sum(id2):
            return -1
        else:
            for i in range(len(id1), 0, -1):
                ind = i - 1
                if id1[ind] > id2[ind]:
                    return 1
                elif id1[ind] < id2[ind]:
                    return -1
            return 0

    def get_tuple_patch_ids(self) -> List[Tuple[int, ...]]:
        unsorted_ids = list(itertools.product([0, 1], repeat=len(self.gene_sets)))[1:]
        sorted_ids = sorted(unsorted_ids, key=functools.cmp_to_key(self._compare_ids))
        return sorted_ids

    def select(self, ind, draw: bool = True):
        graph_modified = self.update_color(ind, self.SELECTED_STATE)
        if graph_modified and draw:
            self.draw()

    def deselect(self, ind, draw: bool = True):
        graph_modified = self.update_color(ind, self.DESELECTED_STATE)
        if graph_modified and draw:
            self.draw()

    def on_hover(self, event):
        graph_modified = False

        for subset in self.subset_states:
            bounding_box = self.bounding_boxes[subset]
            if bounding_box.contains_point((event.x, event.y)):
                if self.subset_states[subset] in [self.HOVER_STATE, self.DESELECTED_STATE]:
                    graph_modified |= self.update_color(subset, self.HOVER_STATE)
                else:
                    graph_modified |= self.update_color(subset, self.HOVER_SELECTED_STATE)

            else:
                if self.subset_states[subset] in [self.HOVER_STATE, self.DESELECTED_STATE]:
                    graph_modified |= self.update_color(subset, self.DESELECTED_STATE)
                else:
                    graph_modified |= self.update_color(subset, self.SELECTED_STATE)

        if graph_modified:
            self.draw()

    def clear_selection(self, draw: bool = True):
        graph_modified = False
        for subset in self.subset_states:
            graph_modified |= self.update_color(subset, self.DESELECTED_STATE)
        if graph_modified and draw:
            self.draw()

    @QtCore.pyqtSlot()
    def union(self):
        for subset in self.subset_states:
            self.select(subset, draw=False)
        self.draw()

    @QtCore.pyqtSlot()
    def intersection(self):
        self.clear_selection(draw=False)
        self.select(len(self.subset_states) - 1, draw=False)
        self.draw()

    @QtCore.pyqtSlot(float)
    def majority_vote_intersection(self, majority_threshold: float):
        thresholds = [-float("inf")] + [(i + 1) * (1 / len(self.gene_sets)) for i in range(len(self.gene_sets))]

        self.clear_selection(draw=False)
        for i in range(len(self.gene_sets)):
            if thresholds[i] < majority_threshold <= thresholds[i + 1]:
                start = sum([generic.combination(len(self.gene_sets), j + 1) for j in range(i)])
                print(i, start)
                for subset in range(start, len(self.subset_states)):
                    self.select(subset, draw=False)
        self.draw()

    @QtCore.pyqtSlot(str)
    def difference(self, primary_set: str):
        self.clear_selection(draw=False)
        if primary_set in self.gene_sets:
            subset_idx = list(self.gene_sets.keys()).index(primary_set)
            self.select(subset_idx, draw=False)
        self.draw()

    def get_custom_selection(self) -> set:
        selection = set()
        patch_ids = self.get_tuple_patch_ids()
        for subset, id in zip(self.subset_states, patch_ids):
            if self.subset_states[subset] in [self.SELECTED_STATE, self.HOVER_SELECTED_STATE]:
                included_sets = [s for s, ind in zip(self.gene_sets.values(), id) if ind]
                excluded_sets = [s for s, ind in zip(self.gene_sets.values(), id) if not ind]
                patch_content = set.intersection(*included_sets).difference(*excluded_sets)
                selection.update(patch_content)
        return selection


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

    def clear_selection(self):
        pass


def get_icon(name: str):
    if name in AVAILABLE_ICONS:
        pth = str(Path(__file__).parent.joinpath('icons', AVAILABLE_ICONS[name]))
        icon = QtGui.QIcon(pth)
        return icon
    elif name == 'blank':
        pixmap = QtGui.QPixmap(32, 32)
        pixmap.fill(QtCore.Qt.transparent)
        return QtGui.QIcon(pixmap)
    return None
