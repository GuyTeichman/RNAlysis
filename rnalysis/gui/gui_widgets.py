import ast
import collections
import functools
import inspect
import json
from pathlib import Path
from queue import Queue
from typing import List, Dict, Tuple, Sequence, Iterable, Union, Callable

import matplotlib
import pandas as pd
from PyQt5 import QtCore, QtWidgets, QtGui
from joblib import Parallel, parallel_backend
from tqdm.auto import tqdm
from typing_extensions import Literal, get_origin, get_args

from rnalysis.utils import parsing, validation, generic, param_typing, settings

EMPTY = inspect._empty


class TableColumnPicker(QtWidgets.QPushButton):
    valueChanged = QtCore.pyqtSignal()
    IS_MULTI_INPUT = True
    __slots__ = {'dialog': 'dialog window of the widget',
                 'dialog_table': 'table containing widgets inside the dialog window',
                 'done_button': 'done button',
                 'select_all_button': 'select all button',
                 'clear_button': 'clear button',
                 'column_labels': 'labels of the differenet columns',
                 'column_checks': 'toggle switches of the different columns'}

    def __init__(self, text: str = 'Choose columns', parent=None):
        super().__init__(text, parent)
        self.columns: list = []

        self.dialog = QtWidgets.QDialog(self)
        self.dialog_table = QtWidgets.QTableWidget()
        self.dialog_layout = QtWidgets.QGridLayout(self.dialog)
        self.done_button = QtWidgets.QPushButton('Done')
        self.select_all_button = QtWidgets.QPushButton('Select all')
        self.clear_button = QtWidgets.QPushButton('Clear selection')
        self.column_labels: List[QtWidgets.QLabel] = []
        self.column_checks: List[ToggleSwitch] = []

        self.init_ui()

    def init_ui(self):
        self.clicked.connect(self.open_dialog)

        self.dialog.setWindowTitle(f"Choose table columns")

        self.select_all_button.clicked.connect(self.select_all)
        self.dialog_layout.addWidget(self.select_all_button, 5, 0, 1, 2)

        self.clear_button.clicked.connect(self.clear_selection)
        self.dialog_layout.addWidget(self.clear_button, 6, 0, 1, 2)

        self.done_button.clicked.connect(self.dialog.close)
        self.done_button.clicked.connect(self.valueChanged.emit)
        self.dialog_layout.addWidget(self.done_button, 7, 0, 1, 2)

        self.dialog_table.setColumnCount(2)
        self.dialog_table.setSortingEnabled(False)
        self.dialog_table.setHorizontalHeaderLabels(['Column names', 'Include column?'])
        self.dialog_layout.addWidget(self.dialog_table, 0, 0, 5, 2)

    def select_all(self):
        for checkbox in self.column_checks:
            checkbox.setChecked(True)

    def clear_selection(self):
        for checkbox in self.column_checks:
            checkbox.setChecked(False)

    def add_columns(self, columns: list):
        self.columns.extend([str(col) for col in columns])
        self.update_table()

    def open_dialog(self):
        self.dialog.exec()

    def value(self):
        picked_cols = []
        for checkbox, col_name in zip(self.column_checks, self.columns):
            if checkbox.isChecked():
                picked_cols.append(col_name)
        return picked_cols

    def _update_window_size(self):
        self.dialog_table.setSizePolicy(QtWidgets.QSizePolicy.Minimum, QtWidgets.QSizePolicy.Minimum)
        # self.dialog_table.setVerticalScrollBarPolicy(QtCore.Qt.ScrollBarAlwaysOff)
        self.dialog_table.setHorizontalScrollBarPolicy(QtCore.Qt.ScrollBarAlwaysOff)
        self.dialog_table.resizeRowsToContents()
        self.dialog_table.resizeColumnsToContents()
        self.dialog_table.horizontalHeader().setSectionResizeMode(QtWidgets.QHeaderView.ResizeToContents)
        self.dialog_table.horizontalHeader().setSectionResizeMode(0, QtWidgets.QHeaderView.Stretch)

        self.dialog_table.resize(self.dialog_table.horizontalHeader().length() +
                                 self.dialog_table.verticalHeader().width(),
                                 self.dialog_table.verticalHeader().length() +
                                 self.dialog_table.horizontalHeader().height())
        self.dialog.resize(self.dialog_table.size())

    def update_table(self):
        if self.dialog_table.rowCount() < len(self.columns):
            start_ind = self.dialog_table.rowCount()
            self.dialog_table.setRowCount(len(self.columns))
            for i in range(start_ind, len(self.columns)):
                col_name = self.columns[i]
                col_label = QtWidgets.QLabel(col_name)
                col_checkbox = ToggleSwitch()

                col_checkbox.setChecked(True)

                self.column_labels.append(col_label)
                self.column_checks.append(col_checkbox)

                self.dialog_table.setCellWidget(i, 0, col_label)
                self.dialog_table.setCellWidget(i, 1, col_checkbox)

            self._update_window_size()


class TableSingleColumnPicker(TableColumnPicker):
    __slots__ = {'button_group': 'group of all of the toggle switches, making sure exactly one is selected'}

    def __init__(self, text: str = 'Choose a column', parent=None):
        super().__init__(text, parent)
        self.button_group = QtWidgets.QButtonGroup(self)
        self.button_group.setExclusive(True)

    def select_all(self):
        raise NotImplementedError

    def init_ui(self):
        super().init_ui()
        self.dialog.setWindowTitle("Choose a table column")
        self.dialog_layout.removeWidget(self.select_all_button)
        self.select_all_button.deleteLater()
        self.dialog_layout.removeWidget(self.clear_button)
        self.clear_button.deleteLater()

    def update_table(self):
        start_ind = self.dialog_table.rowCount()
        super().update_table()
        self.clear_selection()
        for i in range(start_ind, len(self.column_checks)):
            checkbox = self.column_checks[i]
            self.button_group.addButton(checkbox.switch)

    def value(self):
        picked_cols = super().value()
        assert len(picked_cols) > 0, "Not enough columns were picked!"
        assert len(picked_cols) < 2, "Too many columns were picked!"
        return picked_cols[0]


class TableColumnGroupPicker(TableColumnPicker):
    __slots__ = {'column_combos': 'combo boxes for choosing column groups',
                 'colors': 'list of colors for each row',
                 'reset_button': 'button for resetting the selections',
                 '_color_gen': 'generator for choosing row colors'}

    def __init__(self, text: str = 'Choose columns', parent=None):
        self.column_combos: List[QtWidgets.QComboBox] = []

        self._color_gen = generic.color_generator()
        self.colors = []
        self.reset_button = QtWidgets.QPushButton('Reset')
        self.save_button = QtWidgets.QPushButton('Export selection')
        self.load_button = QtWidgets.QPushButton('Import selection')
        super().__init__(text, parent)

    def add_columns(self, columns: list):
        for i in range(len(columns)):
            self.colors.append(matplotlib.colors.to_hex(next(self._color_gen)))
        super().add_columns(columns)

    def set_row_color(self):
        groups = self._get_groups_in_use()
        for row in range(self.dialog_table.rowCount()):
            color = self.colors[groups.index(self.column_combos[row].currentText())] if self.column_checks[
                row].isChecked() else "#FFFFFF"
            for col in range(self.dialog_table.columnCount()):
                self.dialog_table.item(row, col).setBackground(QtGui.QColor(color))
                # self.column_combos[row].setStyleSheet("QComboBox {background-color: " + str(color) + "}")

    def init_ui(self):
        super().init_ui()
        self.dialog_table.setColumnCount(3)
        self.dialog_table.setHorizontalHeaderLabels(['Column names', 'Include column?', 'Column group'])
        self.reset_button.clicked.connect(self.reset)

        self.save_button.clicked.connect(self.export_selection)
        self.dialog_layout.addWidget(self.save_button, 7, 0, 1, 2)

        self.load_button.clicked.connect(self.import_selection)
        self.dialog_layout.addWidget(self.load_button, 8, 0, 1, 2)
        self.dialog_layout.addWidget(self.reset_button, 9, 0, 1, 2)
        self.dialog_layout.addWidget(self.done_button, 10, 0, 1, 2)

    def export_selection(self):
        value = self.value()
        default_name = 'selection.txt'
        output_filename, _ = QtWidgets.QFileDialog.getSaveFileName(self, "Export selection",
                                                                   str(Path.home().joinpath(default_name)),
                                                                   "Text files (*.txt);;"
                                                                   "All Files (*)")
        if output_filename:
            with open(output_filename, 'w') as f:
                json.dump(value, f)

    def import_selection(self):
        input_filename, _ = QtWidgets.QFileDialog.getOpenFileName(self, "Import selection",
                                                                  str(Path.home()),
                                                                  "Text files (*.txt);;"
                                                                  "All Files (*)")
        if input_filename:
            with open(input_filename) as f:
                selection = json.load(f)
            self.set_selection(selection)

    def set_selection(self, selection: list):
        self.clear_selection()
        for group_ind, grp in enumerate(selection):
            for item in grp:
                assert item in self.columns, f"Item '{item}' is not a column in this table!"
                item_ind = self.columns.index(item)
                self.column_checks[item_ind].setChecked(True)
                self.column_combos[item_ind].setCurrentText(str(group_ind + 1))

    def reset(self):
        self.select_all()
        for ind in range(len(self.columns)):
            combo = self.column_combos[ind]
            combo.setCurrentText(str(ind + 1))

    def _get_groups_in_use(self):
        existing_groups = set()
        for check, combo, col_name in zip(self.column_checks, self.column_combos, self.columns):
            if check.isChecked():
                grp = combo.currentText()
                existing_groups.add(grp)
        return sorted(existing_groups)

    def value(self):
        existing_groups = {}
        for combo, col_name in zip(self.column_combos, self.columns):
            if combo.isEnabled():
                grp = combo.currentText()
                if grp not in existing_groups:
                    existing_groups[grp] = []
                existing_groups[grp].append(col_name)
        output = [existing_groups[grp] for grp in sorted(existing_groups.keys())]
        return output

    def update_table(self):
        if self.dialog_table.rowCount() < len(self.columns):
            start_ind = self.dialog_table.rowCount()

            super().update_table()

            for i in range(start_ind, len(self.columns)):
                col_checkbox = self.column_checks[i]
                col_combo = QtWidgets.QComboBox()

                col_checkbox.stateChanged.connect(col_combo.setEnabled)
                col_checkbox.setChecked(True)
                col_checkbox.stateChanged.connect(self.set_row_color)

                col_combo.currentIndexChanged.connect(self.set_row_color)

                self.column_combos.append(col_combo)
                self.dialog_table.setCellWidget(i, 2, col_combo)

                for j in range(self.dialog_table.columnCount()):
                    self.dialog_table.setItem(i, j, QtWidgets.QTableWidgetItem())
            for ind in range(len(self.columns)):
                combo = self.dialog_table.cellWidget(ind, 2)
                combo.clear()
                combo.addItems([str(i + 1) for i in range(len(self.columns))])
                combo.setCurrentText(str(ind + 1))

            self._update_window_size()


class PathInputDialog(QtWidgets.QDialog):
    __slots__ = {'message': 'default message to display',
                 'layout': 'layout of the widget',
                 'button_box': 'button box for accept/reject buttons',
                 'path': 'path line edit widget'}

    def __init__(self, message: str = "No prompt available", parent=None):
        super().__init__(parent)
        self.message = message
        self.layout = QtWidgets.QVBoxLayout(self)
        self.button_box = QtWidgets.QDialogButtonBox(QtWidgets.QDialogButtonBox.Ok)
        self.path = PathLineEdit(parent=self)

        self.init_ui()

    def init_ui(self):
        self.button_box.accepted.connect(self.accept)
        self.button_box.rejected.connect(self.reject)
        self.setWindowTitle('User input required')
        self.layout.addWidget(QtWidgets.QLabel(self.message))
        self.layout.addWidget(self.path)
        self.layout.addWidget(self.button_box)

    def result(self):
        return self.path.text()


class Worker(QtCore.QObject):
    finished = QtCore.pyqtSignal(tuple)
    startProgBar = QtCore.pyqtSignal(object)
    __slots__ = {'partial': 'partial function to run in the Worker',
                 'args': "arguments to emit alongside the partial function's output"}

    def __init__(self, partial, *args):
        self.partial = partial
        self.args = args
        super().__init__()

    def run(self):
        args = []
        try:
            result = self.partial()
        except Exception as e:
            self.finished.emit((e,))
            return

        if result is not None:
            args = self.args
        self.finished.emit((result, *args))


class AltTQDM(QtCore.QObject):
    barUpdate = QtCore.pyqtSignal(int)
    barFinished = QtCore.pyqtSignal()
    __slots__ = {'tqdm': 'underlying tqdm object/progress bar'}

    def __init__(self, iter_obj: Iterable = None, desc: str = '', unit: str = '',
                 bar_format: str = '', total: int = None):
        self.tqdm = tqdm(iter_obj, desc=desc, unit=unit, bar_format=bar_format, total=total)
        super().__init__()

    def __iter__(self):
        for i, item in enumerate(self.tqdm):
            self.barUpdate.emit(1)
            yield item
        self.barFinished.emit()

    def __enter__(self):
        return self

    def update(self, n: int = 1):
        self.barUpdate.emit(n)

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.barFinished.emit()


class AltParallel(QtCore.QObject):
    barUpdate = QtCore.pyqtSignal(int)
    barTotalUpdate = QtCore.pyqtSignal(int)
    barFinished = QtCore.pyqtSignal()
    __slots__ = {'parallel': 'underlying Parallel object',
                 'desc': 'progress bar description',
                 'total': 'progreess bar total',
                 'prev_total': 'previous bar total',
                 'prev_report': 'last sent update of bar progress'}

    def __init__(self, n_jobs: int = -1, total=None, desc: str = '', unit: str = 'it',
                 bar_format: str = '', *args, **kwargs):
        kwargs['backend'] = 'multiprocessing'
        self.parallel = Parallel(*args, n_jobs=n_jobs, **kwargs)
        self.desc = desc
        self.total = total
        self.prev_total = 0
        self.prev_report = 0
        super().__init__()
        self.barUpdate.emit(1)

    def __call__(self, *args, **kwargs):
        parallel_backend('multiprocessing')
        self.parallel.print_progress = functools.partial(self.print_progress)
        return self.parallel.__call__(*args, **kwargs)

    def update(self, n: int = 1):
        self.barUpdate.emit(n)
        print(f"{self.desc}: finished {self.parallel.n_completed_tasks} of {self.prev_total} tasks\r")

    def total_update(self, n: int):
        self.barTotalUpdate.emit(n)

    def print_progress(self):
        if self.total is None and self.prev_total < self.parallel.n_dispatched_tasks:
            self.prev_total = self.parallel.n_dispatched_tasks
            self.total_update(self.parallel.n_dispatched_tasks)

        if self.parallel.n_completed_tasks - self.prev_report > 0:
            self.update(self.parallel.n_completed_tasks - self.prev_report)
            self.prev_report = self.parallel.n_completed_tasks
            if self.prev_report == self.total:
                self.barFinished.emit()


class TextWithCopyButton(QtWidgets.QWidget):
    __slots__ = {'text': 'text',
                 'copy_button': 'copy button',
                 'copied_label': 'label indicating when copy button was pressed',
                 'layout': 'widget layout'}

    def __init__(self, text: str, parent=None):
        super().__init__(parent)
        self.text = text
        self.text_edit = QtWidgets.QTextBrowser(self)
        self.copy_button = QtWidgets.QPushButton('Copy to clipboard')
        self.copied_label = QtWidgets.QLabel()
        self.layout = QtWidgets.QVBoxLayout(self)
        self.init_ui()

    def init_ui(self):
        self.text_edit.setHtml(self.text)
        self.text_edit.setReadOnly(True)
        self.layout.addWidget(self.text_edit)
        self.copy_button.clicked.connect(self.copy_to_clipboard)
        self.layout.addWidget(self.copy_button)
        self.layout.addWidget(self.copied_label)

    def copy_to_clipboard(self):
        cb = QtWidgets.QApplication.clipboard()
        cb.clear(mode=cb.Clipboard)
        cb.setText(self.text_edit.toPlainText())
        self.copied_label.setText('Copied to clipboard')


class ToggleSwitchCore(QtWidgets.QPushButton):
    """
    Based upon a StackOverflow response by the user Heike:
    https://stackoverflow.com/questions/56806987/switch-button-in-pyqt
    """
    stateChanged = QtCore.pyqtSignal(bool)
    RADIUS = 11
    WIDTH = 42
    BORDER = 2
    __slots__ = {}

    def __init__(self, parent=None):
        super().__init__(parent)
        self.setCheckable(True)
        self.setMinimumWidth((self.WIDTH + self.BORDER) * 2)
        self.setMinimumHeight((self.RADIUS + self.BORDER) * 2)
        self.clicked.connect(self.state_changed)

    def state_changed(self):
        self.stateChanged.emit(self.isChecked())

    def setChecked(self, is_checked):
        super(ToggleSwitchCore, self).setChecked(is_checked)
        self.state_changed()

    def paintEvent(self, event):
        label = " True" if self.isChecked() else "False"
        bg_color = QtGui.QColor('#72e5bf') if self.isChecked() else QtGui.QColor('#e96e3a')

        radius = self.RADIUS
        width = self.WIDTH
        center = self.rect().center()

        painter = QtGui.QPainter(self)
        painter.setRenderHint(QtGui.QPainter.Antialiasing)
        painter.translate(center)
        painter.setBrush(QtGui.QColor("#cccccc"))

        pen = QtGui.QPen(QtGui.QColor("#222228"))
        pen.setWidth(self.BORDER)
        painter.setPen(pen)

        painter.drawRoundedRect(QtCore.QRect(-width, -radius, 2 * width, 2 * radius), radius, radius)
        painter.setBrush(QtGui.QBrush(bg_color))
        sw_rect = QtCore.QRect(-radius, -radius, width + radius, 2 * radius)
        if not self.isChecked():
            sw_rect.moveLeft(-width)
        painter.drawRoundedRect(sw_rect, radius, radius)
        painter.drawText(sw_rect, QtCore.Qt.AlignCenter, label)


class ToggleSwitch(QtWidgets.QWidget):
    IS_CHECK_BOX_LIKE = True
    __slots__ = {'switch': 'toggle switch core',
                 'layout': 'widget layout'}

    def __init__(self, parent=None):
        super().__init__(parent)
        self.switch = ToggleSwitchCore(self)
        self.layout = QtWidgets.QHBoxLayout(self)

        self.layout.addWidget(self.switch)
        self.layout.addStretch(1)

        self.setChecked = self.switch.setChecked
        self.stateChanged = self.switch.stateChanged
        self.isChecked = self.switch.isChecked


class ComboBoxOrOtherWidget(QtWidgets.QWidget):
    IS_COMBO_BOX_LIKE = True
    OTHER_TEXT = 'Other...'
    __slots__ = {'layout': 'layout',
                 'combo': 'Literals combo box',
                 'items': 'combo box items',
                 'other': 'other widget',
                 'default': 'default_value'}

    def __init__(self, items: Tuple[str, ...], other: QtWidgets.QWidget, default: str = None, parent=None):
        super().__init__(parent)
        self.layout = QtWidgets.QHBoxLayout(self)
        self.combo = QtWidgets.QComboBox(self)
        self.items = items
        self.other = other
        self.default = default

        self.currentIndexChanged = self.combo.currentIndexChanged
        self.init_ui()

    def clear(self):
        try:
            self.other.clear()
        except AttributeError:
            pass

    def init_ui(self):
        self.layout.addWidget(self.combo)
        self.layout.addWidget(self.other)
        self.combo.addItems(self.items)
        self.combo.addItem(self.OTHER_TEXT)
        self.currentIndexChanged.connect(self.check_other)
        if self.default is not None:
            self.combo.setCurrentText(self.default)
        else:
            self.combo.setCurrentText(self.OTHER_TEXT)
        self.check_other()

    def check_other(self):
        self.other.setEnabled(self.combo.currentText() == self.OTHER_TEXT)

    def currentText(self):
        if self.combo.currentText() == self.OTHER_TEXT:
            return get_val_from_widget(self.other)
        return self.combo.currentText()


class HelpButton(QtWidgets.QToolButton):
    __slots__ = {'param_name': 'name of the parameter',
                 'desc': 'description of the parameter'}

    def __init__(self, parent=None):
        super().__init__(parent)
        self.setIcon(self.style().standardIcon(QtWidgets.QStyle.SP_MessageBoxQuestion))
        self.param_name = ''
        self.desc = ''

    def _disconnect_help(self):
        try:
            self.clicked.disconnect()
        except TypeError:
            pass

    def connect_desc_help(self, desc: str):
        self._disconnect_help()
        self.clicked.connect(self._show_help_desc)
        self.desc = desc

    def connect_param_help(self, param_name: str, desc: str):
        self._disconnect_help()
        self.clicked.connect(self._show_help_param)
        self.param_name = param_name
        self.desc = desc

    def _show_help_param(self):
        QtWidgets.QToolTip.showText(QtGui.QCursor.pos(), f"<b>{self.param_name}:</b> <br>{self.desc}")

    def _show_help_desc(self):
        QtWidgets.QToolTip.showText(QtGui.QCursor.pos(), self.desc)


class ColorPicker(QtWidgets.QWidget):
    IS_LINE_EDIT_LIKE = True
    __slots__ = {'picker_window': 'color picker dialog',
                 'layout': 'widget layout',
                 'pick_button': 'button for opening dialog',
                 'color_line': 'text edit containing the color name/hex code'}

    def __init__(self, default: Union[str, None] = None, parent=None):
        super().__init__(parent)
        self.picker_window = QtWidgets.QColorDialog(self)
        self.layout = QtWidgets.QGridLayout(self)
        self.pick_button = QtWidgets.QPushButton("Pick color", self)
        self.color_line = QtWidgets.QLineEdit(self)
        self.init_ui()

        self.textChanged = self.color_line.textChanged
        self.textChanged.connect(self.update_color)

        if default is not None:
            self.color_line.setText(default)

    def init_ui(self):
        self.pick_button.clicked.connect(self.get_color)

        self.layout.addWidget(self.pick_button, 0, 0)
        self.layout.addWidget(self.color_line, 0, 2)

    def get_color(self):
        color = QtWidgets.QColorDialog.getColor()
        if color.isValid():
            self.color_line.setText(color.name())
            self.update_color()

    def update_color(self):
        try:
            color = self.text()
            text_color = matplotlib.colors.to_hex(
                [abs(i - j) for i, j in zip(matplotlib.colors.to_rgb('white'), matplotlib.colors.to_rgb(color))])
            self.color_line.setStyleSheet("QLineEdit {background : " + color + "; \ncolor : " + text_color + ";}")
        except ValueError:
            pass

    def setText(self, text: str):
        self.color_line.setText(text)
        self.update_color()

    def text(self):
        try:
            return matplotlib.colors.to_hex(self.color_line.text())
        except ValueError:
            return None


class MultipleChoiceList(QtWidgets.QWidget):
    __slots__ = {'layout': 'layout',
                 'items': 'list items',
                 'list': 'list widget',
                 'select_all_button': 'select all button',
                 'clear_all_button': 'clear all button',
                 'current_layout_row': 'current bottom-most row in the layout'}

    def __init__(self, items: Sequence, icons: Sequence = None, parent=None):
        super().__init__(parent)
        self.layout = QtWidgets.QGridLayout(self)
        self.setLayout(self.layout)

        self.items = []
        self.list_items = []
        self.list = QtWidgets.QListWidget(self)
        self.list.setSelectionMode(QtWidgets.QAbstractItemView.MultiSelection)
        self.layout.addWidget(self.list, 1, 1, 4, 1)

        self.select_all_button = QtWidgets.QPushButton('Select all', self)
        self.select_all_button.clicked.connect(self.select_all)
        self.layout.addWidget(self.select_all_button, 5, 1)

        self.clear_all_button = QtWidgets.QPushButton('Clear selection', self)
        self.clear_all_button.clicked.connect(self.clear_all)
        self.layout.addWidget(self.clear_all_button, 6, 1)

        self.current_layout_row = 6

        for row in range(2, 4):
            self.layout.setRowStretch(row, 2)

        self.itemSelectionChanged = self.list.itemSelectionChanged

        self.add_items(items, icons)

    def get_sorted_selection(self):
        items = self.list.selectedItems()
        return sorted(items, key=self.list.row)

    def get_sorted_selected_names(self):
        return [name for name in self.items if name in {item.text() for item in self.list.selectedItems()}]

    def get_sorted_names(self):
        return self.items

    def add_item(self, item, icon=None):
        self.items.append(item)
        list_item = QtWidgets.QListWidgetItem(item)
        if icon is not None:
            list_item.setIcon(icon)
        self.list_items.append(list_item)
        self.list.addItem(list_item)

    def add_items(self, items, icons=None):
        for i, item in enumerate(items):
            self.add_item(item, None if icons is None else icons[i])

    def select_all(self):
        for ind in range(self.list.count()):
            item = self.list.item(ind)
            if not item.isSelected():
                item.setSelected(True)

    def clear_all(self):
        for item in self.list.selectedItems():
            item.setSelected(False)


class MultiChoiceListWithDelete(MultipleChoiceList):
    itemDeleted = QtCore.pyqtSignal(int)

    def __init__(self, items: Sequence, icons: Sequence = None, parent=None, delete_text: str = 'delete'):
        super().__init__(items, icons, parent)
        self.delete_text = delete_text
        self.delete_button = QtWidgets.QPushButton(f'{delete_text.capitalize()} selected')
        self.delete_button.clicked.connect(self.delete_selected)
        self.delete_all_button = QtWidgets.QPushButton(f'{delete_text.capitalize()} all')
        self.delete_all_button.clicked.connect(self.delete_all)

        self.layout.addWidget(self.delete_button, self.current_layout_row + 1, 1)
        self.layout.addWidget(self.delete_all_button, self.current_layout_row + 2, 1)
        self.current_layout_row += 2

    def delete_selected(self):
        for item in self.list.selectedItems():
            row = self.list.row(item)
            self.list.takeItem(row)
            self.items.pop(row)
            self.list_items.pop(row)
            self.itemDeleted.emit(row)

    def delete_all(self):
        accepted = QtWidgets.QMessageBox.question(self, f"{self.delete_text.capitalize()} all items?",
                                                  f"Are you sure you want to {self.delete_text} all items?",
                                                  QtWidgets.QMessageBox.Yes | QtWidgets.QMessageBox.No)
        if accepted == QtWidgets.QMessageBox.Yes:
            for n_item in reversed(range(len(self.items))):
                self.itemDeleted.emit(n_item)
            self.delete_all_quietly()

    def delete_all_quietly(self):
        self.items = []
        self.list_items = []
        self.list.clear()


class MultiChoiceListWithReorder(MultipleChoiceList):
    itemOrderChanged = QtCore.pyqtSignal()
    __slots__ = {'top_button': 'send to top',
                 'up_button': 'send up',
                 'down_button': 'send down',
                 'bottm_button': 'send to bottom',
                 'reorder_label': 'label'}

    def __init__(self, items: Sequence, icons: Sequence = None, parent=None):
        super().__init__(items, icons, parent)
        self.top_button = QtWidgets.QPushButton('Top')
        self.top_button.clicked.connect(self.top_selected)

        self.up_button = QtWidgets.QPushButton('Up')
        self.up_button.clicked.connect(self.up_selected)

        self.down_button = QtWidgets.QPushButton('Down')
        self.down_button.clicked.connect(self.down_selected)

        self.bottom_button = QtWidgets.QPushButton('Bottom')
        self.bottom_button.clicked.connect(self.bottom_selected)

        self.reorder_label = QtWidgets.QLabel('<b>Reorder:</b>')
        self.reorder_label.setAlignment(QtCore.Qt.AlignmentFlag.AlignHCenter)
        self.layout.addWidget(self.reorder_label, self.current_layout_row + 1, 1)
        self.layout.addWidget(self.top_button, self.current_layout_row + 2, 1)
        self.layout.addWidget(self.up_button, self.current_layout_row + 3, 1)
        self.layout.addWidget(self.down_button, self.current_layout_row + 4, 1)
        self.layout.addWidget(self.bottom_button, self.current_layout_row + 5, 1)
        self.current_layout_row += 5

    def up_selected(self):
        for item in self.get_sorted_selection():
            row = self.list.row(item)
            if row > 0:
                self.items[row - 1], self.items[row] = self.items[row], self.items[row - 1]
                self.list_items[row - 1], self.list_items[row] = self.list_items[row], self.list_items[row - 1]
                current_item = self.list.takeItem(row)
                self.list.insertItem(row - 1, current_item)
                current_item.setSelected(True)
        self.itemOrderChanged.emit()

    def down_selected(self):
        for item in reversed(self.get_sorted_selection()):
            row = self.list.row(item)
            end = len(self.items) - 1
            if row < end:
                self.items[row + 1], self.items[row] = self.items[row], self.items[row + 1]
                self.list_items[row + 1], self.list_items[row] = self.list_items[row], self.list_items[row + 1]
                current_item = self.list.takeItem(row)
                self.list.insertItem(row + 1, current_item)
                current_item.setSelected(True)
        self.itemOrderChanged.emit()

    def top_selected(self):
        for i, item in enumerate(self.get_sorted_selection()):
            row = self.list.row(item)
            if row > i:
                self.items.insert(i, self.items.pop(row))
                self.list_items.insert(i, self.list_items.pop(row))
                current_item = self.list.takeItem(row)
                self.list.insertItem(i, current_item)
                current_item.setSelected(True)
        self.itemOrderChanged.emit()

    def bottom_selected(self):
        for i, item in enumerate(reversed(self.get_sorted_selection())):
            row = self.list.row(item)
            this_end = len(self.items) - 1 - i
            if row < this_end:
                self.items.insert(this_end, self.items.pop(row))
                self.list_items.insert(this_end, self.list_items.pop(row))

                current_item = self.list.takeItem(row)
                self.list.insertItem(this_end, current_item)
                current_item.setSelected(True)
        self.itemOrderChanged.emit()


class MultiChoiceListWithDeleteReorder(MultiChoiceListWithReorder, MultiChoiceListWithDelete):
    itemDeleted = QtCore.pyqtSignal(int)

    def __init__(self, items: Sequence, icons: Sequence = None, parent=None):
        super().__init__(items, icons, parent=parent)


class FileListWidgetItem(QtWidgets.QListWidgetItem):
    __slots__ = {'file_path': 'path of the file',
                 'display_name': 'display name of the file path'}

    def __init__(self, file_path, parent=None):
        self.file_path = file_path
        self.display_name = Path(file_path).stem
        super().__init__(self.display_name, parent)

    def filePath(self):
        return self.file_path


class OrderedFileList(MultiChoiceListWithDeleteReorder):
    __slots__ = {'add_files_button': 'add files button'}

    def __init__(self, parent=None):
        super().__init__([], parent=parent)
        self.add_files_button = QtWidgets.QPushButton('Add files...', self)
        self.add_files_button.clicked.connect(self.add_files)
        self.layout.addWidget(self.add_files_button, 0, 1)

    def add_files(self):
        filenames, _ = QtWidgets.QFileDialog.getOpenFileNames(self, "Load fastq files")
        if filenames:
            self.add_items(filenames)

    def add_item(self, item, icon=None):
        self.items.append(item)
        list_item = FileListWidgetItem(item)
        if icon is not None:
            list_item.setIcon(icon)
        self.list_items.append(list_item)
        self.list.addItem(list_item)

    def get_sorted_names(self):
        return [name for name in self.items]


class MandatoryComboBox(QtWidgets.QComboBox):
    __slots__ = {'default_choice': 'default choice'}

    def __init__(self, default_choice: str, parent=None):
        super().__init__(parent)
        self.default_choice = default_choice
        self.addItem(self.default_choice)
        self.currentTextChanged.connect(self.set_bg_color)
        self.set_bg_color()

    def clear(self) -> None:
        super().clear()
        self.addItem(self.default_choice)

    def set_bg_color(self):
        if self.is_legal():
            self.setStyleSheet("MandatoryComboBox{border: 1.5px solid #57C4AD;}")
        else:
            self.setStyleSheet("MandatoryComboBox{border: 1.5px solid #DB4325;}")

    def disable_bg_color(self):
        self.setStyleSheet("MandatoryComboBox{}")

    def is_legal(self):
        return self.currentText() != self.default_choice

    def setEnabled(self, to_enable: bool):
        self.setDisabled(not to_enable)

    def setDisabled(self, to_disable: bool):
        if to_disable:
            self.disable_bg_color()
        else:
            self.set_bg_color()
        super().setDisabled(to_disable)


class GeneSetComboBox(MandatoryComboBox):
    IS_SPIN_BOX_LIKE = True
    boxOpened = QtCore.pyqtSignal()
    __slots__ = {'available_objects': 'available gene sets and tables'}

    def __init__(self, parent=None):
        super().__init__('Choose gene set...', parent)
        self.available_objects = dict()

    def showPopup(self) -> None:
        self.boxOpened.emit()
        super().showPopup()

    def update_gene_sets(self, available_objects: Dict[str, Tuple[QtWidgets.QWidget, QtGui.QIcon]]):
        # if the gene sets didn't change, don't do anything
        if available_objects.keys() == self.available_objects.keys() and [item[0] for item in available_objects] == [
            item[0] for item in self.available_objects]:
            return
        # clear the previous list of gene sets, and add the current available gene sets to the list
        self.clear()
        self.available_objects = available_objects
        for obj_name in self.available_objects:
            self.addItem(self.available_objects[obj_name][1], obj_name)

    def value(self):
        set_name = self.currentText()
        if set_name == self.default_choice:
            return set()
        assert set_name in self.available_objects, f"Tab '{set_name}' does not exist!"
        return self.available_objects[set_name][0].obj()


class MinMaxDialog(QtWidgets.QDialog):
    __slots__ = {}

    def __init__(self, parent=None):
        super().__init__(parent)
        self.setWindowFlag(QtCore.Qt.WindowMinMaxButtonsHint)


class TrueFalseBoth(QtWidgets.QWidget):
    IS_MULTI_INPUT = True
    STYLESHEET = '''QPushButton::checked {background-color : green;
            color: white;
            border: 1px solid #32ba32;
            border-radius: 4px;}'''
    selectionChanged = QtCore.pyqtSignal()
    __slots__ = {'layout': 'layout',
                 'true_button': 'True button',
                 'false_button': 'False button'}

    def __init__(self, default=None, parent=None):
        super().__init__(parent)
        self.layout = QtWidgets.QHBoxLayout(self)
        self.true_button = QtWidgets.QPushButton("True")
        self.true_button.setCheckable(True)
        self.true_button.setMinimumSize(90, 30)
        self.true_button.setStyleSheet(self.STYLESHEET)

        self.false_button = QtWidgets.QPushButton("False")
        self.false_button.setCheckable(True)
        self.false_button.setMinimumSize(90, 30)
        self.false_button.setStyleSheet(self.STYLESHEET)

        self.true_button.toggled.connect(self.selectionChanged.emit)
        self.false_button.toggled.connect(self.selectionChanged.emit)

        self.layout.addWidget(self.true_button)
        self.layout.addWidget(self.false_button)
        self.layout.addStretch(1)

        self.setValue(default)

    def setValue(self, default):
        default = parsing.data_to_list(default)
        if True in default:
            self.true_button.click()
        if False in default:
            self.false_button.click()

    def value(self):
        checked = []
        if self.true_button.isChecked():
            checked.append(True)
        if self.false_button.isChecked():
            checked.append(False)
        return checked


class PathLineEdit(QtWidgets.QWidget):
    IS_LINE_EDIT_LIKE = True
    textChanged = QtCore.pyqtSignal(bool)
    __slots__ = {'file_path': 'line edit',
                 'open_button': 'open button',
                 'is_file': 'is the current text pointing to an exisintg file',
                 '_is_legal': 'is the current path legal',
                 'layout': 'layout'}

    def __init__(self, contents: str = 'auto', button_text: str = 'Load', is_file: bool = True, parent=None):
        super().__init__(parent)
        self.file_path = QtWidgets.QLineEdit('', self)
        self.open_button = QtWidgets.QPushButton(button_text, self)
        self.is_file = is_file
        self._is_legal = False

        self.layout = QtWidgets.QGridLayout(self)
        self.setLayout(self.layout)
        self.layout.addWidget(self.open_button, 1, 0)
        self.layout.addWidget(self.file_path, 1, 1)

        self.file_path.textChanged.connect(self._check_legality)
        if self.is_file:
            self.open_button.clicked.connect(self.choose_file)
        else:
            self.open_button.clicked.connect(self.choose_folder)

        if contents == 'auto':
            if self.is_file:
                contents = 'No file chosen'
            else:
                contents = 'No folder chosen'
        self.file_path.setText(contents)

    def clear(self):
        self.file_path.clear()

    @property
    def is_legal(self):
        return self._is_legal

    def _check_legality(self):
        current_path = self.file_path.text()
        if (self.is_file and validation.is_legal_file_path(current_path)) or (
            not self.is_file and validation.is_legal_dir_path(current_path)):
            self._is_legal = True
        else:
            self._is_legal = False
        self.set_file_path_bg_color()
        self.textChanged.emit(self.is_legal)

    def set_file_path_bg_color(self):
        if self.is_legal:
            self.file_path.setStyleSheet("QLineEdit{border: 1.5px solid #57C4AD;}")
        else:
            self.file_path.setStyleSheet("QLineEdit{border: 1.5px solid #DB4325;}")

    def disable_bg_color(self):
        self.file_path.setStyleSheet("QLineEdit{}")

    def setEnabled(self, to_enable: bool):
        self.setDisabled(not to_enable)

    def setDisabled(self, to_disable: bool):
        if to_disable:
            self.disable_bg_color()
        else:
            self.set_file_path_bg_color()
        super().setDisabled(to_disable)

    def choose_file(self):
        filename, _ = QtWidgets.QFileDialog.getOpenFileName(self, "Choose a file")
        if filename:
            self.file_path.setText(filename)

    def choose_folder(self):
        dirname = QtWidgets.QFileDialog.getExistingDirectory(self, "Choose a folder")
        if dirname:
            self.file_path.setText(dirname)

    def text(self):
        return self.file_path.text()

    def setText(self, text: str):
        return self.file_path.setText(text)


class StrIntLineEdit(QtWidgets.QLineEdit):
    IS_LINE_EDIT_LIKE = True

    def __init__(self, text: Union[str, int] = '', parent=None):
        super().__init__(str(text), parent)

    def text(self) -> Union[str, int]:
        val = super().text()
        if val.lstrip('-').isnumeric():
            return int(val)
        return val

    def setText(self, p_str):
        super().setText(str(p_str))


class RadioButtonBox(QtWidgets.QGroupBox):
    selectionChanged = QtCore.pyqtSignal()
    __slots__ = {'button_box': 'button group',
                 'radio_layout': 'radio button layout',
                 'radio_buttons': 'radio buttons'}

    def __init__(self, title: str, actions, is_flat=False, parent=None):
        super().__init__(title, parent)
        self.setFlat(is_flat)
        self.button_box = QtWidgets.QButtonGroup()
        self.buttonClicked = self.button_box.buttonClicked
        self.checkedButton = self.button_box.checkedButton
        self.radio_layout = QtWidgets.QGridLayout()
        self.radio_buttons = {}
        self.setLayout(self.radio_layout)
        self.add_items(actions)

    def add_items(self, actions: Iterable):
        for action in actions:
            if isinstance(action, str):
                self.add_item(action)
            elif isinstance(action, tuple):
                title = action[0]
                sub_actions = action[1]
                self.radio_layout.addWidget(QtWidgets.QLabel(title), self.radio_layout.count(), 0, 1, 10)
                for subaction in sub_actions:
                    self.add_item(subaction, indent=True)

    def add_item(self, action: str, indent: bool = False):
        self.radio_buttons[action] = QtWidgets.QRadioButton(action)
        self.button_box.addButton(self.radio_buttons[action])
        if indent:
            self.radio_layout.addWidget(self.radio_buttons[action], self.radio_layout.count(), 1, 1, 10)
        else:
            self.radio_layout.addWidget(self.radio_buttons[action], self.radio_layout.count(), 0, 1, 9)

    def set_selection(self, selection: Union[str, int]):
        if isinstance(selection, str):
            for button_name, button in self.radio_buttons.items():
                if button_name == selection:
                    button.click()
        else:
            for i, button in enumerate(self.radio_buttons.values()):
                if i == selection:
                    button.click()
        self.selectionChanged.emit()


class OptionalWidget(QtWidgets.QWidget):
    __slots__ = {'layout': 'layout',
                 'other': 'other widget',
                 'default': 'default value',
                 'checkbox': 'disable checkbox'}

    def __init__(self, other: QtWidgets.QWidget, default=EMPTY, parent=None):
        super().__init__(parent)
        self.layout = QtWidgets.QHBoxLayout(self)
        self.other = other
        self.default = default
        self.checkbox = QtWidgets.QCheckBox('Disable input?')
        self.toggled = self.checkbox.toggled

        self.init_ui()

    def clear(self):
        self.checkbox.setChecked(False)
        try:
            self.other.clear()
        except AttributeError:
            pass

    def init_ui(self):
        self.toggled.connect(self.other.setDisabled)
        self.layout.addWidget(self.checkbox)
        self.layout.addWidget(self.other)
        if self.default is None:
            self.checkbox.setChecked(True)

    def check_other(self):
        self.other.setDisabled(self.checkbox.isChecked())

    def value(self):
        if self.checkbox.isChecked():
            return None
        return get_val_from_widget(self.other)

    def setValue(self, val):
        if val is None:
            self.checkbox.setChecked(True)
        else:
            try:
                self.checkbox.setChecked(False)
                self.other.setValue(val)
            except AttributeError:
                raise ValueError(f'Unable to set default value of {type(self.other)} to {val}')


class ComparisonPicker(QtWidgets.QWidget):
    __slots__ = {'design_mat': 'design matrix',
                 'layout': 'layout',
                 'factor': 'factor combo box',
                 'numerator': 'numerator combo box',
                 'denominator': 'denominator combo box'}

    def __init__(self, design_mat: pd.DataFrame, parent=None):
        super().__init__(parent)
        self.design_mat = design_mat
        self.layout = QtWidgets.QHBoxLayout(self)
        self.factor = QtWidgets.QComboBox(self)
        self.numerator = QtWidgets.QComboBox(self)
        self.denominator = QtWidgets.QComboBox(self)

        self.init_ui()

    def init_ui(self):
        self.layout.addWidget(self.factor)
        self.layout.addWidget(self.numerator)
        self.layout.addWidget(self.denominator)

        self.factor.currentTextChanged.connect(self.update_combos)
        self.factor.addItems([str(item) for item in self.design_mat.columns])

    def update_combos(self):
        this_factor = self.factor.currentText()
        if this_factor in self.design_mat.columns:
            options = sorted({str(item) for item in self.design_mat[this_factor]})
        else:
            options = ["Select a factor..."]
        self.numerator.clear()
        self.denominator.clear()

        self.numerator.addItems(options)
        self.denominator.addItems(options)

    def value(self) -> Tuple[str, str, str]:
        return self.factor.currentText(), self.numerator.currentText(), self.denominator.currentText()


class ComparisonPickerGroup(QtWidgets.QWidget):
    __slots__ = {'design_mat': 'design matrix',
                 'widgets': 'widgets',
                 'layout': 'layouts',
                 'inputs': 'inputs',
                 'input_labels': 'labels for inputs'}

    def __init__(self, design_mat: pd.DataFrame, parent=None):
        super().__init__(parent)
        self.design_mat = design_mat
        self.widgets = {}
        self.layout = QtWidgets.QGridLayout(self)

        self.inputs = []
        self.input_labels = []

        self.init_ui()
        self.add_comparison_widget()

    def init_ui(self):
        self.widgets['add_widget'] = QtWidgets.QPushButton('Add comparison')
        self.widgets['add_widget'].clicked.connect(self.add_comparison_widget)
        self.layout.addWidget(self.widgets['add_widget'], 0, 0, 1, 4)

        self.widgets['remove_widget'] = QtWidgets.QPushButton('Remove comparison')
        self.widgets['remove_widget'].clicked.connect(self.remove_comparison_widget)
        self.layout.addWidget(self.widgets['remove_widget'], 1, 0, 1, 4)

        self.layout.addWidget(QtWidgets.QLabel('<b>Design factor</b>', self), 2, 1)
        self.layout.addWidget(QtWidgets.QLabel('<b>Numerator</b>', self), 2, 2)
        self.layout.addWidget(QtWidgets.QLabel('<b>Denominator</b>', self), 2, 3)

    @QtCore.pyqtSlot()
    def add_comparison_widget(self):
        self.inputs.append(ComparisonPicker(self.design_mat, self))
        self.layout.addWidget(self.inputs[-1], len(self.inputs) + 2, 1, 1, 3)

        self.input_labels.append(QtWidgets.QLabel(f'Comparison #{len(self.inputs)}:', self.inputs[-1]))
        self.layout.addWidget(self.input_labels[-1], len(self.inputs) + 2, 0)

        self.layout.setRowStretch(len(self.inputs) + 2, 1)
        self.layout.setRowStretch(len(self.inputs) + 3, 2)

    @QtCore.pyqtSlot()
    def remove_comparison_widget(self):
        if len(self.inputs) > 0:
            self.inputs.pop(-1).deleteLater()
            self.input_labels.pop(-1).deleteLater()
            self.layout.setRowStretch(len(self.inputs) + 3, 2)

    def get_comparison_values(self):
        values = []
        for widget in self.inputs:
            values.append(widget.value())
        return values


class QMultiInput(QtWidgets.QPushButton):
    IS_MULTI_INPUT = True
    CHILD_QWIDGET = None
    valueChanged = QtCore.pyqtSignal()
    __slots__ = {'label': 'multi widget label',
                 'dialog_widgets': 'dict of the dialog widgets',
                 'dialog_started': 'stores whether the dialog was already started >=1 times',
                 'dialog_layout': 'layout of the dialog window'}

    def __init__(self, label: str = '', text='Set input', parent=None):
        super().__init__(text, parent)
        self.label = label
        self.dialog_widgets = {}
        self.dialog_started: bool = False
        self.init_ui()
        self.dialog_layout = QtWidgets.QGridLayout(self.dialog_widgets['box'])
        self.clicked.connect(self.start_dialog)
        self.init_dialog_ui()

    def clear(self):
        while len(self.dialog_widgets['inputs']) > 0:
            self.remove_widget()

    @QtCore.pyqtSlot()
    def start_dialog(self):
        self.dialog_widgets['box'].exec()

    def init_ui(self):
        self.dialog_widgets['box'] = QtWidgets.QDialog(self.parent())

    def init_dialog_ui(self):
        self.dialog_widgets['box'].setWindowTitle(f"Input for parameter '{self.label}'")
        self.dialog_widgets['box'].resize(400, 175)
        self.dialog_widgets['inputs'] = []
        self.dialog_widgets['input_labels'] = []

        self.dialog_widgets['add_widget'] = QtWidgets.QPushButton('Add field')
        self.dialog_widgets['add_widget'].clicked.connect(self.add_widget)
        self.dialog_layout.addWidget(self.dialog_widgets['add_widget'], 0, 0, 1, 2)

        self.dialog_widgets['remove_widget'] = QtWidgets.QPushButton('Remove field')
        self.dialog_widgets['remove_widget'].clicked.connect(self.remove_widget)
        self.dialog_layout.addWidget(self.dialog_widgets['remove_widget'], 1, 0, 1, 2)

        self.dialog_widgets['done'] = QtWidgets.QPushButton('Done')
        self.dialog_widgets['done'].clicked.connect(self.dialog_widgets['box'].close)
        self.dialog_widgets['done'].clicked.connect(self.valueChanged.emit)
        self.dialog_layout.addWidget(self.dialog_widgets['done'], 2, 0, 1, 2)

        self.add_widget()

    @QtCore.pyqtSlot()
    def add_widget(self):
        self.dialog_widgets['inputs'].append(self.CHILD_QWIDGET())
        self.dialog_layout.addWidget(self.dialog_widgets['inputs'][-1], len(self.dialog_widgets['inputs']) + 2, 1)

        self.dialog_widgets['input_labels'].append(
            QtWidgets.QLabel(f'{self.label}:', self.dialog_widgets['inputs'][-1]))
        self.dialog_layout.addWidget(self.dialog_widgets['input_labels'][-1], len(self.dialog_widgets['inputs']) + 2, 0)

    @QtCore.pyqtSlot()
    def remove_widget(self):
        if len(self.dialog_widgets['inputs']) > 0:
            self.dialog_widgets['inputs'].pop(-1).deleteLater()
            self.dialog_widgets['input_labels'].pop(-1).deleteLater()

    def value(self):
        values = []
        for widget in self.dialog_widgets['inputs']:
            values.append(self.get_widget_value(widget))
        return values[0] if len(values) == 1 else values

    def get_widget_value(self, widget):
        raise NotImplementedError

    def set_widget_value(self, ind: int, val):
        raise NotImplementedError

    def setValue(self, value: Iterable):
        values = parsing.data_to_list(value)
        while len(self.dialog_widgets['inputs']) > 0:
            self.remove_widget()
        for item in values:
            if item is None:
                continue
            self.add_widget()
            self.set_widget_value(-1, item)


class MultiColorPicker(QMultiInput):
    CHILD_QWIDGET = ColorPicker

    def get_widget_value(self, widget: type(CHILD_QWIDGET)):
        return widget.text()

    def set_widget_value(self, ind: int, val):
        self.dialog_widgets['inputs'][ind].setText(val)


class QMultiSpinBox(QMultiInput):
    __slots__ = {'minimum': 'minimum value for spin boxes',
                 'maximum': 'maximum value for spinboxes'}
    CHILD_QWIDGET = QtWidgets.QSpinBox

    def __init__(self, label: str = '', text='Set input', parent=None, minimum=-2147483648, maximum=2147483647):
        self.minimum = minimum
        self.maximum = maximum
        super().__init__(label, text, parent)

    def add_widget(self):
        super().add_widget()
        self.dialog_widgets['inputs'][-1].setMinimum(self.minimum)
        self.dialog_widgets['inputs'][-1].setMaximum(self.maximum)

    def get_widget_value(self, widget: type(CHILD_QWIDGET)):
        return widget.value()

    def set_widget_value(self, ind: int, val):
        self.dialog_widgets['inputs'][ind].setValue(val)


class QMultiDoubleSpinBox(QMultiSpinBox):
    __slots__ = {'minimum': 'minimum value for spinboxes',
                 'maximum': 'maximum value for spinboxes',
                 'step_size': 'default step size of spinboxes'}
    CHILD_QWIDGET = QtWidgets.QDoubleSpinBox

    def __init__(self, label: str = '', text='Set input', parent=None, minimum=float("-inf"), maximum=float("inf"),
                 step_size: float = 0.05):
        self.minimum = minimum
        self.maximum = maximum
        self.step_size = step_size
        super().__init__(label, text, parent)

    def add_widget(self):
        super().add_widget()
        self.dialog_widgets['inputs'][-1].setMinimum(self.minimum)
        self.dialog_widgets['inputs'][-1].setMaximum(self.maximum)
        self.dialog_widgets['inputs'][-1].setSingleStep(self.step_size)


class QMultiLineEdit(QMultiInput):
    CHILD_QWIDGET = QtWidgets.QLineEdit

    def get_widget_value(self, widget: type(CHILD_QWIDGET)):
        return widget.text()

    def set_widget_value(self, ind: int, val):
        self.dialog_widgets['inputs'][ind].setText(str(val))

    def value(self):
        res = super().value()
        if res == '':
            return []
        return res


class QMultiPathLineEdit(QMultiLineEdit):
    CHILD_QWIDGET = PathLineEdit

    def __init__(self, is_file: bool = True, label: str = '', text='Set input', parent=None):
        self.is_file = is_file
        super().__init__(label, text, parent)

    def add_widget(self):
        self.dialog_widgets['inputs'].append(self.CHILD_QWIDGET(is_file=self.is_file))
        self.dialog_layout.addWidget(self.dialog_widgets['inputs'][-1], len(self.dialog_widgets['inputs']) + 2, 1)

        self.dialog_widgets['input_labels'].append(
            QtWidgets.QLabel(f'{self.label}:', self.dialog_widgets['inputs'][-1]))
        self.dialog_layout.addWidget(self.dialog_widgets['input_labels'][-1], len(self.dialog_widgets['inputs']) + 2, 0)


class TwoLayerMultiLineEdit(QMultiInput):
    CHILD_QWIDGET = QMultiLineEdit

    def get_widget_value(self, widget: type(CHILD_QWIDGET)):
        return widget.value()

    def set_widget_value(self, ind: int, val):
        self.dialog_widgets['inputs'][ind].setValue(val)


class QMultiStrIntLineEdit(QMultiLineEdit):
    CHILD_QWIDGET = StrIntLineEdit


class QMultiComboBox(QMultiInput):
    CHILD_QWIDGET = QtWidgets.QComboBox
    __slots__ = {'items': 'combo box items'}

    def __init__(self, label: str, text: str = 'Set Input', parent=None, items=()):
        self.items = items
        super().__init__(label, text, parent)

    def add_widget(self):
        super().add_widget()
        self.dialog_widgets['inputs'][-1].addItems(self.items)

    def get_widget_value(self, widget: type(CHILD_QWIDGET)):
        return widget.currentText()

    def set_widget_value(self, ind: int, val):
        self.dialog_widgets['inputs'][ind].setCurrentText(val)


class QMultiToggleSwitch(QMultiInput):
    CHILD_QWIDGET = ToggleSwitch

    def get_widget_value(self, widget: type(CHILD_QWIDGET)):
        return widget.isChecked()

    def set_widget_value(self, ind: int, val: bool):
        self.dialog_widgets['inputs'][ind].setChecked(val)


class ThreadStdOutStreamTextQueueReceiver(QtCore.QObject):
    queue_stdout_element_received_signal = QtCore.pyqtSignal(str)
    __slots__ = {'queue': 'queue'}

    def __init__(self, q: Queue, *args, **kwargs):
        QtCore.QObject.__init__(self, *args, **kwargs)
        self.queue = q

    @QtCore.pyqtSlot()
    def run(self):
        self.queue_stdout_element_received_signal.emit('Welcome to RNAlysis!\n')
        while True:
            text = self.queue.get()
            self.queue_stdout_element_received_signal.emit(text)


class StdOutTextEdit(QtWidgets.QTextEdit):
    __slots__ = {'carriage': 'carriage',
                 'prev_coord': 'previous coordinate'}

    def __init__(self, parent=None):
        super().__init__(parent)
        self.setReadOnly(True)
        self.setLineWidth(50)
        self.setMinimumWidth(500)
        self.document().setMaximumBlockCount(500)
        self.setFont(QtGui.QFont('Consolas', 11))
        self.carriage = False
        self.prev_coord = 0

    @QtCore.pyqtSlot(str)
    def append_text(self, text: str):

        self.moveCursor(QtGui.QTextCursor.End)
        if text == '\n':
            return
        text = text.replace("<", "&lt;").replace(">", "&gt;")

        if self.carriage:
            self.carriage = False
            diff = self.document().characterCount() - self.prev_coord
            cursor = self.textCursor()
            cursor.movePosition(QtGui.QTextCursor.PreviousCharacter, QtGui.QTextCursor.MoveAnchor, n=diff)
            cursor.movePosition(QtGui.QTextCursor.End, QtGui.QTextCursor.KeepAnchor)
            cursor.removeSelectedText()

        if text.endswith('\r'):
            self.carriage = True
            self.prev_coord = self.document().characterCount()
            text.rstrip('\r')

        if text.startswith('Warning: '):
            self.insertHtml(f'<div style="color:red;">{text}</div><br>')
        else:
            self.insertHtml(f'<div style="color:black;">{text}</div><br>')
        self.verticalScrollBar().setValue(self.verticalScrollBar().maximum())


class WriteStream(QtCore.QObject):
    message = QtCore.pyqtSignal(str)
    __slots__ = {'buffer': 'buffer'}

    def __init__(self, q: Queue, parent=None):
        super(WriteStream, self).__init__(parent)

        self.buffer = q

    def write(self, text):
        self.buffer.put(text)

    def flush(self):
        pass

    def close(self):
        pass


class NewParam:
    __slots__ = {'annotation': 'annotation',
                 'default': 'default'}

    def __init__(self, annotation, default=EMPTY):
        self.annotation = annotation
        self.default = default


def param_to_widget(param, name: str,
                    actions_to_connect: Union[Iterable[Callable], Callable] = tuple(), pipeline_mode: bool = False):
    column_annotations = {param_typing.GroupedColumns: TableColumnGroupPicker,
                          param_typing.ColumnNames: TableColumnPicker,
                          param_typing.ColumnName: TableSingleColumnPicker}
    actions_to_connect = parsing.data_to_tuple(actions_to_connect)

    if param.default == EMPTY:
        is_default = False
    else:
        is_default = True
    # if the annotation is a Union that contains None, turn the None into a checkbox,
    # and the rest of the annotation into an optional widget (AKA OptionalWidget widget)
    if get_origin(param.annotation) == Union and type(None) in get_args(param.annotation):
        args = get_args(param.annotation)
        none_ind = args.index(type(None))
        without_none = tuple(args[0:none_ind] + args[none_ind + 1:])
        if param.default is None:
            this_default = None
            other_default = EMPTY
        else:
            this_default = EMPTY
            other_default = param.default
        widget = OptionalWidget(
            param_to_widget(NewParam(Union[without_none], other_default), name, actions_to_connect, pipeline_mode),
            this_default)
        for action in actions_to_connect:
            widget.toggled.connect(action)

    # if the annotation is a Union that contains Literals, turn the Literals into a combobox,
    # and the rest of the annotation into an 'other...' option in the combo box (AKA ComboBoxOrOther widget)
    elif get_origin(param.annotation) == Union and Literal in [get_origin(ann) for ann in get_args(param.annotation)]:
        args = get_args(param.annotation)
        literal_ind = [get_origin(ann) for ann in args].index(Literal)
        literal = args[literal_ind]
        without_literal = tuple(args[0:literal_ind] + args[literal_ind + 1:])
        if param.default in get_args(literal):
            this_default = param.default
            other_default = EMPTY
        else:
            this_default = None
            other_default = param.default
        widget = ComboBoxOrOtherWidget(get_args(literal),
                                       param_to_widget(NewParam(Union[without_literal], other_default), name,
                                                       actions_to_connect, pipeline_mode), this_default)
        for action in actions_to_connect:
            widget.currentIndexChanged.connect(action)

    elif param.annotation == param_typing.Color:
        widget = ColorPicker(param.default if is_default else '')
        for action in actions_to_connect:
            widget.textChanged.connect(action)

    elif param.annotation == param_typing.ColorList:
        widget = MultiColorPicker(name)
        widget.setValue(param.default if is_default else '')
        for action in actions_to_connect:
            widget.valueChanged.connect(action)

    elif param.annotation in column_annotations:
        if pipeline_mode:
            param.annotation = param_typing.type_to_supertype(param.annotation)
            widget = param_to_widget(param, name, actions_to_connect, pipeline_mode)
        else:
            cls = column_annotations[param.annotation]
            widget = cls()
        for action in actions_to_connect:
            widget.valueChanged.connect(action)

    elif param.annotation == Sequence[str]:
        widget = GeneSetComboBox()
        for action in actions_to_connect:
            widget.currentTextChanged.connect(action)

    elif param.annotation == bool:
        widget = ToggleSwitch()
        default = param.default if is_default else False
        widget.setChecked(default)
        for action in actions_to_connect:
            widget.stateChanged.connect(action)

    elif param.annotation == int:
        widget = QtWidgets.QSpinBox()
        widget.setMinimum(-2147483648)
        widget.setMaximum(2147483647)
        default = param.default if is_default else 0
        widget.setValue(default)
        for action in actions_to_connect:
            widget.valueChanged.connect(action)

    elif param.annotation == param_typing.PositiveInt:
        widget = QtWidgets.QSpinBox()
        widget.setMinimum(1)
        widget.setMaximum(2147483647)
        default = param.default if is_default else 0
        widget.setValue(default)
        for action in actions_to_connect:
            widget.valueChanged.connect(action)

    elif param.annotation == param_typing.NonNegativeInt:
        widget = QtWidgets.QSpinBox()
        widget.setMinimum(0)
        widget.setMaximum(2147483647)
        default = param.default if is_default else 0
        widget.setValue(default)
        for action in actions_to_connect:
            widget.valueChanged.connect(action)

    elif param.annotation == param_typing.NegativeInt:
        widget = QtWidgets.QSpinBox()
        widget.setMinimum(-2147483648)
        widget.setMaximum(-1)
        default = param.default if is_default else 0
        widget.setValue(default)
        for action in actions_to_connect:
            widget.valueChanged.connect(action)


    elif param.annotation == float:
        widget = QtWidgets.QDoubleSpinBox()
        widget.setMinimum(float("-inf"))
        widget.setMaximum(float("inf"))
        widget.setSingleStep(0.05)
        default = param.default if is_default else 0.0
        widget.setValue(default)
        for action in actions_to_connect:
            widget.valueChanged.connect(action)

    elif param.annotation == param_typing.Fraction:
        widget = QtWidgets.QDoubleSpinBox()
        widget.setMinimum(0.0)
        widget.setMaximum(1.0)
        widget.setSingleStep(0.05)
        default = param.default if is_default else 0.0
        widget.setValue(default)
        for action in actions_to_connect:
            widget.valueChanged.connect(action)

    elif param.annotation == str:
        widget = QtWidgets.QLineEdit(param.default if is_default else '')
        for action in actions_to_connect:
            widget.textChanged.connect(action)

    elif get_origin(param.annotation) == Literal:
        widget = QtWidgets.QComboBox()
        widget.addItems(get_args(param.annotation))
        for action in actions_to_connect:
            widget.currentIndexChanged.connect(action)
        if is_default:
            widget.setCurrentText(str(param.default))

    elif param.annotation in (Union[str, List[str]], Union[str, Iterable[str]], List[str], Iterable[str]):
        widget = QMultiLineEdit(name)
        if is_default:
            widget.setValue(param.default)
        for action in actions_to_connect:
            widget.valueChanged.connect(action)

    elif param.annotation == List[Union[str, Path]]:
        is_file = 'folder' not in name
        widget = QMultiPathLineEdit(is_file, name)
        if is_default:
            widget.setValue(param.default)
        for action in actions_to_connect:
            widget.valueChanged.connect(action)

    elif param.annotation == List[Iterable[str]]:
        widget = TwoLayerMultiLineEdit(name)
        if is_default:
            widget.setValue(param.default)
        for action in actions_to_connect:
            widget.valueChanged.connect(action)

    elif param.annotation in (Union[float, List[float]], Union[float, Iterable[float]]):
        widget = QMultiDoubleSpinBox(name)
        if is_default:
            widget.setValue(param.default)
        for action in actions_to_connect:
            widget.valueChanged.connect(action)

    elif param.annotation in (Union[param_typing.Fraction, List[param_typing.Fraction]],
                              Union[param_typing.Fraction, Iterable[param_typing.Fraction]]):
        widget = QMultiDoubleSpinBox(name, minimum=0, maximum=1, step_size=0.05)
        if is_default:
            widget.setValue(param.default)
        for action in actions_to_connect:
            widget.valueChanged.connect(action)

    elif param.annotation in (Union[int, List[int]], Union[int, Iterable[int]], List[int], Iterable[int]):
        widget = QMultiSpinBox(name)
        if is_default:
            widget.setValue(param.default)
        for action in actions_to_connect:
            widget.valueChanged.connect(action)

    elif param.annotation in (Union[param_typing.PositiveInt, List[param_typing.PositiveInt]],
                              Union[param_typing.PositiveInt, Iterable[param_typing.PositiveInt]],
                              List[param_typing.PositiveInt], Iterable[param_typing.PositiveInt]):
        widget = QMultiSpinBox(name, minimum=1)
        if is_default:
            widget.setValue(param.default)
        for action in actions_to_connect:
            widget.valueChanged.connect(action)

    elif param.annotation in (Union[param_typing.NonNegativeInt, List[param_typing.NonNegativeInt]],
                              Union[param_typing.NonNegativeInt, Iterable[param_typing.NonNegativeInt]],
                              List[param_typing.NonNegativeInt], Iterable[param_typing.NonNegativeInt]):
        widget = QMultiSpinBox(name, minimum=0)
        if is_default:
            widget.setValue(param.default)
        for action in actions_to_connect:
            widget.valueChanged.connect(action)

    elif param.annotation in (Union[param_typing.NegativeInt, List[param_typing.NegativeInt]],
                              Union[param_typing.NegativeInt, Iterable[param_typing.NegativeInt]],
                              List[param_typing.NegativeInt], Iterable[param_typing.NegativeInt]):
        widget = QMultiSpinBox(name, maximum=-1)
        if is_default:
            widget.setValue(param.default)
        for action in actions_to_connect:
            widget.valueChanged.connect(action)

    elif param.annotation in (Union[bool, List[bool]], Union[bool, Iterable[bool]], List[bool], Iterable[bool]):
        widget = QMultiToggleSwitch(name)
        if is_default:
            widget.setValue(param.default)
        for action in actions_to_connect:
            widget.valueChanged.connect(action)

    elif param.annotation == Union[str, int]:
        widget = StrIntLineEdit(param.default if is_default else '')
        for action in actions_to_connect:
            widget.textChanged.connect(action)

    elif param.annotation in (Union[str, int, Iterable[str], Iterable[int]], Union[str, int, List[str], List[int]],
                              Union[List[int], List[str]]):
        widget = QMultiStrIntLineEdit(name)
        widget.setValue(param.default if is_default else '')
        for action in actions_to_connect:
            widget.valueChanged.connect(action)

    elif param.annotation in (Path, Union[str, Path], Union[None, str, Path]):
        is_file = 'folder' not in name
        widget = PathLineEdit(is_file=is_file)
        for action in actions_to_connect:
            widget.textChanged.connect(action)
        if is_default:
            widget.setText(str(param.default))

    elif get_origin(param.annotation) in (collections.abc.Iterable, list, tuple, set) and get_origin(
        get_args(param.annotation)[0]) == Literal:
        widget = QMultiComboBox(name, items=get_args(get_args(param.annotation)[0]))
        if is_default:
            widget.setValue(param.default)
        for action in actions_to_connect:
            widget.valueChanged.connect(action)

    elif param.annotation == Union[bool, Tuple[bool, bool]]:
        widget = TrueFalseBoth(param.default)
        for action in actions_to_connect:
            widget.selectionChanged.connect(action)

    else:
        widget = QtWidgets.QTextEdit()
        default = param.default if is_default else ''
        widget.setText(str(default))
        widget.setToolTip('Enter a Python expression here')
        widget.setToolTipDuration(0)
        for action in actions_to_connect:
            widget.textChanged.connect(action)
    return widget


def get_val_from_widget(widget):
    if isinstance(widget, (QtWidgets.QSpinBox, QtWidgets.QDoubleSpinBox)) or (hasattr(widget, 'IS_SPIN_BOX_LIKE')):
        val = widget.value()
    elif isinstance(widget, QtWidgets.QLineEdit) or (hasattr(widget, 'IS_LINE_EDIT_LIKE')):
        val = widget.text()
    elif isinstance(widget, QtWidgets.QCheckBox) or (hasattr(widget, 'IS_CHECK_BOX_LIKE')):
        val = widget.isChecked()
    elif isinstance(widget, QtWidgets.QTextEdit):
        val = ast.literal_eval(widget.toPlainText())
    elif isinstance(widget, QtWidgets.QComboBox) or hasattr(widget, 'IS_COMBO_BOX_LIKE'):
        val = widget.currentText()
    elif hasattr(widget, 'IS_MULTI_INPUT'):
        val = widget.value()
    else:
        try:
            val = widget.value()
        except AttributeError:
            raise TypeError(f"Invalid QtWidget type {type(widget)}.")
    return val


def clear_layout(layout, exceptions: set = frozenset()):
    while layout.count() > len(exceptions):
        child = layout.takeAt(0)
        if child.widget() and child.widget() not in exceptions:
            child.widget().deleteLater()


class ReactiveHeaderView(QtWidgets.QHeaderView):
    LOOKUP_DATABASES_PATH = Path(__file__).parent.parent.joinpath('data_files/lookup_databases.json')
    __slots__ = {'context_menu': 'context menu'}

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.setSectionsClickable(True)
        self.context_menu = None
        self.db_actions = None

    def contextMenu(self, value: str):
        self.context_menu = QtWidgets.QMenu(self)
        copy_action = QtWidgets.QAction(f'Copy "{value}"')
        copy_action.triggered.connect(functools.partial(QtWidgets.QApplication.clipboard().setText, value))
        self.context_menu.addAction(copy_action)

        with open(self.LOOKUP_DATABASES_PATH) as f:
            databases = json.load(f)

        self.db_actions = []
        for name in settings.get_databases_settings():
            action = QtWidgets.QAction(f'Search "{value}" on {name}')
            self.db_actions.append(action)
            open_url_partial = functools.partial(QtGui.QDesktopServices.openUrl,
                                                 QtCore.QUrl(f'{databases[name]}{value}'))
            action.triggered.connect(open_url_partial)
            self.context_menu.addAction(action)

        self.context_menu.exec(QtGui.QCursor.pos())

    def mousePressEvent(self, event: QtGui.QMouseEvent):
        if event.button() == QtCore.Qt.RightButton:
            point = event.pos()
            if point.isNull():
                return
            ind = self.logicalIndexAt(point)
            if ind == -1:
                return
            self.contextMenu(str(self.model().headerData(ind, self.orientation())))

        else:
            super().mousePressEvent(event)


class ReactiveListWidget(QtWidgets.QListWidget):
    LOOKUP_DATABASES_PATH = Path(__file__).parent.parent.joinpath('data_files/lookup_databases.json')
    __slots__ = {'context_menu': 'context menu'}

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.context_menu = None
        self.setUniformItemSizes(True)
        self.setSpacing(2)
        self.db_actions = None

    def contextMenu(self, value: str):
        self.context_menu = QtWidgets.QMenu(self)
        copy_action = QtWidgets.QAction(f'Copy "{value}"')
        copy_action.triggered.connect(functools.partial(QtWidgets.QApplication.clipboard().setText, value))
        self.context_menu.addAction(copy_action)

        with open(self.LOOKUP_DATABASES_PATH) as f:
            databases = json.load(f)

        self.db_actions = []
        for name in settings.get_databases_settings():
            action = QtWidgets.QAction(f'Search "{value}" on {name}')
            self.db_actions.append(action)
            open_url_partial = functools.partial(QtGui.QDesktopServices.openUrl,
                                                 QtCore.QUrl(f'{databases[name]}{value}'))
            action.triggered.connect(open_url_partial)
            self.context_menu.addAction(action)

        self.context_menu.exec(QtGui.QCursor.pos())

    def mousePressEvent(self, event: QtGui.QMouseEvent):
        if event.button() == QtCore.Qt.RightButton:
            point = event.pos()
            if point.isNull():
                return
            item = self.itemAt(point)
            if item == -1:
                return
            self.contextMenu(item.text())

        else:
            super().mousePressEvent(event)


class ReactiveTableView(QtWidgets.QTableView):
    __slots__ = {'context_menu': 'context menu'}

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.setHorizontalHeader(ReactiveHeaderView(QtCore.Qt.Orientation.Horizontal))
        self.setVerticalHeader(ReactiveHeaderView(QtCore.Qt.Orientation.Vertical))
        self.context_menu = None

    def contextMenu(self, value: str):
        self.context_menu = QtWidgets.QMenu(self)
        copy_action = QtWidgets.QAction(f'Copy "{value}"')
        copy_action.triggered.connect(functools.partial(QtWidgets.QApplication.clipboard().setText, value))
        self.context_menu.addAction(copy_action)

        self.context_menu.exec(QtGui.QCursor.pos())

    def mousePressEvent(self, event: QtGui.QMouseEvent):
        if event.button() == QtCore.Qt.RightButton:
            point = event.pos()
            if point.isNull():
                return
            ind = self.indexAt(point)
            if ind == -1:
                return
            self.contextMenu(str(self.model().data(ind)))

        else:
            super().mousePressEvent(event)
