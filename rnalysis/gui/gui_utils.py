import ast
import functools
import inspect
import traceback
import typing
import typing_extensions
import collections
from pathlib import Path
from queue import Queue
import time
import matplotlib
import pandas as pd
from PyQt5 import QtCore, QtWidgets, QtGui

from rnalysis import __version__
from rnalysis.utils import io, parsing, settings, validation
from joblib import Parallel


class PathInputDialog(QtWidgets.QDialog):
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


class ProgressSerialGui:
    def __init__(self, iter_obj: typing.Iterable, desc: str = '', unit: str = '', bar_format: str = '',
                 total: int = None, parent=None):
        self.iter_obj = iter_obj
        self.desc = desc
        self.parent = parent
        if total is not None:
            self.total = total
        else:
            try:
                self.total = len(iter_obj)
            except TypeError:
                self.total = 2
        self.dialog = QtWidgets.QProgressDialog(self.desc, "Cancel", 0, self.total, parent)
        self.dialog.setMinimumDuration(0)
        self.dialog.setWindowTitle(self.desc)
        self.dialog.setValue(0)
        self.dialog.setWindowModality(QtCore.Qt.WindowModal)

    def __iter__(self):
        start_time = time.time()
        for i, item in enumerate(self.iter_obj):
            elapsed_time = time.time() - start_time
            completed_items = i + 1
            remaining_time = (elapsed_time / completed_items) * (self.total - completed_items)
            self.dialog.setValue(completed_items)
            self.dialog.setLabelText(self.desc + '\n' + f"Elapsed time: {elapsed_time:.2f} seconds"
                                     + '\n' + f"Remaining time: {remaining_time:.2f} seconds")
            QtWidgets.QApplication.processEvents()
            yield item
        self.dialog.close()


class ProgressParallelGui(Parallel):
    def __init__(self, total=None, desc: str = '', unit: str = 'it', bar_format: str = '', parent=None, *args,
                 **kwargs):
        self.total = total
        self.desc = desc
        # self.dialog = QtWidgets.QProgressDialog(desc, "Cancel", 0, 2 if total is None else total, parent)
        # self.dialog.setMinimumDuration(0)
        # self.dialog.setWindowTitle(self.desc)
        # self.dialog.setValue(0)
        # self.dialog.setWindowModality(QtCore.Qt.WindowModal)
        super().__init__(verbose=100, *args, **kwargs)

    def __call__(self, *args, **kwargs):
        return Parallel.__call__(self, *args, **kwargs)

    def print_progress(self):
        super().print_progress()
        QtWidgets.QApplication.processEvents()


class ComboBoxOrOtherWidget(QtWidgets.QWidget):
    IS_COMBO_BOX_LIKE = True
    OTHER_TEXT = 'Other...'

    def __init__(self, items: typing.List[str], other: QtWidgets.QWidget, default: str = None, parent=None):
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
    def __init__(self, parent=None):
        super().__init__(parent)
        self.setIcon(self.style().standardIcon(QtWidgets.QStyle.SP_MessageBoxQuestion))
        self.param_name = ''
        self.desc = ''

    def connect_help(self, param_name: str, desc: str):
        self.clicked.connect(self.show_help)
        self.param_name = param_name
        self.desc = desc

    def show_help(self):
        QtWidgets.QToolTip.showText(QtGui.QCursor.pos(), f"<b>{self.param_name}:</b> <br>{self.desc}")


class ColorPicker(QtWidgets.QWidget):
    IS_LINE_EDIT_LIKE = True

    def __init__(self, default: typing.Union[str, None] = None, parent=None):
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
    def __init__(self, items: typing.Sequence, icons: typing.Sequence = None, parent=None):
        super().__init__(parent)
        self.layout = QtWidgets.QGridLayout(self)
        self.setLayout(self.layout)

        self.items = items
        self.list_items = []
        self.list = QtWidgets.QListWidget(self)
        for i, item in enumerate(self.items):
            list_item = QtWidgets.QListWidgetItem(item)
            if icons is not None:
                list_item.setIcon(icons[i])
            self.list_items.append(list_item)
            self.list.addItem(list_item)
        self.list.setSelectionMode(QtWidgets.QAbstractItemView.MultiSelection)
        self.layout.addWidget(self.list, 1, 1, 4, 1)

        self.select_all_button = QtWidgets.QPushButton('Select all', self)
        self.select_all_button.clicked.connect(self.select_all)
        self.layout.addWidget(self.select_all_button, 5, 1)

        self.clear_all_button = QtWidgets.QPushButton('Clear all', self)
        self.clear_all_button.clicked.connect(self.clear_all)
        self.layout.addWidget(self.clear_all_button, 6, 1)

        for row in range(2, 4):
            self.layout.setRowStretch(row, 2)

        self.itemSelectionChanged = self.list.itemSelectionChanged
        self.selectedItems = self.list.selectedItems

    def select_all(self):
        for ind in range(self.list.count()):
            item = self.list.item(ind)
            if not item.isSelected():
                item.setSelected(True)

    def clear_all(self):
        for item in self.list.selectedItems():
            item.setSelected(False)


class MandatoryComboBox(QtWidgets.QComboBox):
    def __init__(self, default_choice: str, parent=None):
        super().__init__(parent)
        self.default_choice = default_choice
        self.addItem(self.default_choice)
        self.currentTextChanged.connect(self.set_bg_color)

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

    def setDisabled(self, to_disable: bool):
        if to_disable:
            self.disable_bg_color()
        else:
            self.set_bg_color()
        super().setDisabled(to_disable)


class MinMaxDialog(QtWidgets.QDialog):
    def __init__(self, parent=None):
        super().__init__(parent)
        self.setWindowFlag(QtCore.Qt.WindowMinMaxButtonsHint)


class PathLineEdit(QtWidgets.QWidget):
    IS_LINE_EDIT_LIKE = True
    textChanged = QtCore.pyqtSignal(bool)

    def __init__(self, contents: str = 'No file chosen', button_text: str = 'Load', parent=None):
        super().__init__(parent)
        self.file_path = QtWidgets.QLineEdit('', self)
        self.open_button = QtWidgets.QPushButton(button_text, self)
        self._is_legal = False

        self.layout = QtWidgets.QGridLayout(self)
        self.setLayout(self.layout)
        self.layout.addWidget(self.open_button, 1, 0)
        self.layout.addWidget(self.file_path, 1, 1)

        self.file_path.textChanged.connect(self._check_legality)
        self.open_button.clicked.connect(self.choose_file)

        self.file_path.setText(contents)

    def clear(self):
        self.file_path.clear()

    @property
    def is_legal(self):
        return self._is_legal

    def _check_legality(self):
        current_path = self.file_path.text()
        if validation.is_legal_file_path(current_path):
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

    def choose_file(self):
        filename, _ = QtWidgets.QFileDialog.getOpenFileName(self, "Choose a file", str(Path.home()),
                                                            "All Files (*)")
        if filename:
            self.file_path.setText(filename)

    def text(self):
        return self.file_path.text()

    def setText(self, text: str):
        return self.file_path.setText(text)


class StrIntLineEdit(QtWidgets.QLineEdit):
    IS_LINE_EDIT_LIKE = True

    def text(self) -> typing.Union[str, int]:
        val = super().text()
        if val.lstrip('-').isnumeric():
            return int(val)
        return val


class RadioButtonBox(QtWidgets.QGroupBox):
    selectionChanged = QtCore.pyqtSignal()

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

    def add_items(self, actions: typing.Iterable):
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

    def set_selection(self, selection: typing.Union[str, int]):
        if isinstance(selection, str):
            for button_name, button in self.radio_buttons.items():
                button.setChecked(button_name == selection)
        else:
            for i, button in enumerate(self.radio_buttons.values()):
                button.setChecked(i == selection)
        self.selectionChanged.emit()


class SpinBoxWithDisable(QtWidgets.QSpinBox):
    def changeEvent(self, e):
        if e.type() == QtCore.QEvent.EnabledChange:
            self.lineEdit().setVisible(self.isEnabled())
        return super().changeEvent(e)


class DoubleSpinBoxWithDisable(QtWidgets.QDoubleSpinBox):
    def changeEvent(self, e):
        if e.type() == QtCore.QEvent.EnabledChange:
            self.lineEdit().setVisible(self.isEnabled())
        return super().changeEvent(e)


class OptionalLineEdit(QtWidgets.QWidget):
    IS_LINE_EDIT_LIKE = True

    def __init__(self, parent=None):
        super().__init__(parent)
        self.layout = QtWidgets.QGridLayout(self)
        self.line = QtWidgets.QLineEdit()
        self.checkbox = QtWidgets.QCheckBox('Disable input?')
        self.checkbox.toggled.connect(self.line.setDisabled)

        self.layout.addWidget(self.checkbox, 0, 0)
        self.layout.addWidget(self.line, 0, 1)

        self.textChanged = self.line.textChanged

    def clear(self):
        self.checkbox.setChecked(False)
        self.line.clear()

    def setText(self, val):
        if val is None:
            self.checkbox.setChecked(True)
        else:
            self.line.setText(val)

    def text(self):
        if self.checkbox.isChecked():
            return None
        return self.line.text()


class OptionalSpinBox(QtWidgets.QWidget):
    IS_SPIN_BOX_LIKE = True

    def __init__(self, parent=None):
        super().__init__(parent)
        self.layout = QtWidgets.QGridLayout()
        self.spinbox = SpinBoxWithDisable(self)
        self.checkbox = QtWidgets.QCheckBox('Disable input?')
        self.checkbox.toggled.connect(self.spinbox.setDisabled)

        self.layout.addWidget(self.checkbox, 0, 0)
        self.layout.addWidget(self.spinbox, 0, 1)
        self.setLayout(self.layout)

        self.valueChanged = self.spinbox.valueChanged

    def clear(self):
        self.checkbox.setChecked(False)
        self.spinbox.clear()

    def setValue(self, val):
        if val is None:
            self.checkbox.setChecked(True)
        else:
            self.spinbox.setValue(val)

    def value(self):
        if self.checkbox.isChecked():
            return None
        return self.spinbox.value()

    def setMinimum(self, min_val: int):
        self.spinbox.setMinimum(min_val)

    def setMaximum(self, max_val: int):
        self.spinbox.setMaximum(max_val)


class OptionalDoubleSpinBox(OptionalSpinBox):
    IS_SPIN_BOX_LIKE = True

    def __init__(self, parent=None):
        super().__init__(parent)
        self.checkbox.disconnect()
        self.layout.removeWidget(self.spinbox)
        self.spinbox = DoubleSpinBoxWithDisable(self)
        self.checkbox.toggled.connect(self.spinbox.setDisabled)
        self.valueChanged = self.spinbox.valueChanged
        self.layout.addWidget(self.spinbox, 0, 1)

    def setSingleStep(self, step: float):
        self.spinbox.setSingleStep(step)


class QMultiInput(QtWidgets.QPushButton):
    IS_MULTI_INPUT = True
    CHILD_QWIDGET = None
    valueChanged = QtCore.pyqtSignal()

    def __init__(self, label: str, text='Set input', parent=None):
        super().__init__(text, parent)
        self.label = label
        self.dialog_widgets = {}
        self.dialog_started: bool = False
        self.init_ui()
        self.dialog_layout = QtWidgets.QGridLayout(self.dialog_widgets['box'])
        self.clicked.connect(self.start_dialog)
        self.init_dialog_ui()

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
        self.dialog_widgets['done'].clicked.connect(self.value_changed)
        self.dialog_layout.addWidget(self.dialog_widgets['done'], 2, 0, 1, 2)

        self.add_widget()

    def value_changed(self):
        self.valueChanged.emit()

    @QtCore.pyqtSlot()
    def add_widget(self):
        self.dialog_widgets['inputs'].append(self.CHILD_QWIDGET())
        self.dialog_layout.addWidget(self.dialog_widgets['inputs'][-1], len(self.dialog_widgets['inputs']) + 2, 1)

        self.dialog_widgets['input_labels'].append(
            QtWidgets.QLabel(f'{self.label}:', self.dialog_widgets['inputs'][-1]))
        self.dialog_layout.addWidget(self.dialog_widgets['input_labels'][-1], len(self.dialog_widgets['inputs']) + 2, 0)

    @QtCore.pyqtSlot()
    def remove_widget(self):
        self.dialog_widgets['inputs'].pop(-1).deleteLater()
        self.dialog_widgets['input_labels'].pop(-1).deleteLater()

    def get_values(self):
        values = []
        for widget in self.dialog_widgets['inputs']:
            values.append(self.get_widget_value(widget))
        return values[0] if len(values) == 1 else values

    def get_widget_value(self, widget):
        raise NotImplementedError

    def set_widget_value(self, ind: int, val):
        raise NotImplementedError

    def set_defaults(self, defaults: typing.Iterable):
        defaults = parsing.data_to_list(defaults)
        while len(self.dialog_widgets['inputs']) > 0:
            self.remove_widget()
        for item in defaults:
            self.add_widget()
            self.set_widget_value(-1, item)


class QMultiInputNamedDict(QMultiInput):
    def __init__(self, label: str, text='Set input', parent=None):
        super().__init__(label, text, parent)
        self.dialog_widgets['inputs_keys'] = []

    @QtCore.pyqtSlot()
    def add_widget(self):
        self.dialog_widgets['inputs'].append(self.CHILD_QWIDGET())
        self.dialog_widgets['inputs_keys'].append(QtWidgets.QLineEdit())
        self.dialog_layout.addWidget(self.dialog_widgets['inputs'][-1], len(self.dialog_widgets['inputs']) + 2, 2)
        self.dialog_layout.addWidget(self.dialog_widgets['inputs_keys'][-1],
                                     len(self.dialog_widgets['inputs_keys']) + 2, 1)

        self.dialog_widgets['input_labels'].append(
            QtWidgets.QLabel(f'{self.label}:', self.dialog_widgets['inputs'][-1]))
        self.dialog_layout.addWidget(self.dialog_widgets['input_labels'][-1], len(self.dialog_widgets['inputs']) + 2, 0)


class MultiColorPicker(QMultiInput):
    CHILD_QWIDGET = ColorPicker

    def get_widget_value(self, widget: type(CHILD_QWIDGET)):
        return widget.text()

    def set_widget_value(self, ind: int, val):
        self.dialog_widgets['inputs'][ind].setText(val)


class QMultiSpinBox(QMultiInput):
    CHILD_QWIDGET = QtWidgets.QSpinBox

    def add_widget(self):
        super().add_widget()
        self.dialog_widgets['inputs'][-1].setMinimum(-2147483648)
        self.dialog_widgets['inputs'][-1].setMaximum(2147483647)

    def get_widget_value(self, widget: type(CHILD_QWIDGET)):
        return widget.value()

    def set_widget_value(self, ind: int, val):
        self.dialog_widgets['inputs'][ind].setValue(val)


class QMultiDoubleSpinBox(QMultiSpinBox):
    CHILD_QWIDGET = QtWidgets.QDoubleSpinBox

    def add_widget(self):
        super(QMultiInput).add_widget()
        self.dialog_widgets['inputs'][-1].setMinimum(float("-inf"))
        self.dialog_widgets['inputs'][-1].setMaximum(float("inf"))


class QMultiLineEdit(QMultiInput):
    CHILD_QWIDGET = QtWidgets.QLineEdit

    def get_widget_value(self, widget: type(CHILD_QWIDGET)):
        return widget.text()

    def set_widget_value(self, ind: int, val):
        self.dialog_widgets['inputs'][ind].setText(str(val))


class QMultiStrIntLineEdit(QMultiLineEdit):
    CHILD_QWIDGET = StrIntLineEdit


class QMultiComboBox(QMultiInput):
    CHILD_QWIDGET = QtWidgets.QComboBox

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


class QMultiBoolComboBox(QMultiComboBox):
    def __init__(self, label: str, text: str = 'Set Input', parent=None):
        super().__init__(label, text, parent, items=['True', 'False'])

    def get_widget_value(self, widget: QtWidgets.QComboBox):
        return ast.literal_eval(widget.currentText())

    def set_widget_value(self, ind: int, val):
        self.dialog_widgets['inputs'][ind].setCurrentText(str(val))


class DataFrameModel(QtCore.QAbstractTableModel):
    """
    Basde upon:
    https://stackoverflow.com/a/44605011
    """
    DtypeRole = QtCore.Qt.UserRole + 1000
    ValueRole = QtCore.Qt.UserRole + 1001

    def __init__(self, df=pd.DataFrame(), parent=None):
        super(DataFrameModel, self).__init__(parent)
        self._dataframe = df

    def setDataFrame(self, dataframe):
        self.beginResetModel()
        self._dataframe = dataframe.copy()
        self.endResetModel()

    def dataFrame(self):
        return self._dataframe

    dataFrame = QtCore.pyqtProperty(pd.DataFrame, fget=dataFrame, fset=setDataFrame)

    @QtCore.pyqtSlot(int, QtCore.Qt.Orientation, result=str)
    def headerData(self, section: int, orientation: QtCore.Qt.Orientation, role: int = QtCore.Qt.DisplayRole):
        if role == QtCore.Qt.DisplayRole:
            if orientation == QtCore.Qt.Horizontal:
                return self._dataframe.columns[section]
            else:
                return str(self._dataframe.index[section])
        return QtCore.QVariant()

    def rowCount(self, parent=QtCore.QModelIndex()):
        if parent.isValid():
            return 0
        return len(self._dataframe.index)

    def columnCount(self, parent=QtCore.QModelIndex()):
        if parent.isValid():
            return 0
        return self._dataframe.columns.size

    def data(self, index, role=QtCore.Qt.DisplayRole):
        if not index.isValid() or not (0 <= index.row() < self.rowCount()
                                       and 0 <= index.column() < self.columnCount()):
            return QtCore.QVariant()
        row = self._dataframe.index[index.row()]
        col = self._dataframe.columns[index.column()]
        dt = self._dataframe[col].dtype

        try:
            val = self._dataframe.loc[row, col]
        except IndexError:
            print(row, col, self._dataframe.shape)
        if role == QtCore.Qt.DisplayRole:
            return str(val)
        elif role == DataFrameModel.ValueRole:
            return val
        if role == DataFrameModel.DtypeRole:
            return dt
        return QtCore.QVariant()

    def roleNames(self):
        roles = {
            QtCore.Qt.DisplayRole: b'display',
            DataFrameModel.DtypeRole: b'dtype',
            DataFrameModel.ValueRole: b'value'
        }
        return roles


class DataFramePreviewModel(DataFrameModel):
    def __init__(self, df=pd.DataFrame(), parent=None):
        df = df.iloc[0:2, 0:3]
        df['...'] = '...'
        df.loc['...', :] = '...'
        super().__init__(df, parent)


class DataView(QtWidgets.QWidget):
    def __init__(self, data, name: str, parent=None):
        super().__init__(parent)
        self.data = data
        self.name = name
        self.label_font = QtGui.QFont()
        self.layout = QtWidgets.QVBoxLayout(self)

    def init_ui(self):
        self.setGeometry(1000, 200, 1000, 500)
        self.setWindowTitle(f"View of '{self.name}'")
        self.label_font.setPointSize(14)
        self.label.setFont(self.label_font)
        self.layout.addWidget(self.label)
        self.layout.addWidget(self.save_button)
        self.layout.addWidget(self.data_view)
        self.save_button.clicked.connect(self.save)
        # self.table.setFlags(self.table.flags() & ~QtCore.Qt.ItemIsEditable)
        self.setLayout(self.layout)

    def closeEvent(self, a0: QtGui.QCloseEvent) -> None:
        self.deleteLater()

    def save(self):
        raise NotImplementedError


class GeneSetView(DataView):
    def __init__(self, data: set, name: str, parent=None):
        super().__init__(data, name, parent)
        self.label = QtWidgets.QLabel(f"Gene set '{name}': {len(self.data)} features")

        self.data_view = QtWidgets.QListWidget()
        self.save_button = QtWidgets.QPushButton('Save gene set', self)

        self.init_ui()

    def init_ui(self):
        super().init_ui()
        self.data_view.addItems([str(item) for item in self.data])

    def save(self):
        default_name = str(self.name) + '.txt'
        filename, _ = QtWidgets.QFileDialog.getSaveFileName(self, "Save gene set",
                                                            str(Path.home().joinpath(default_name)),
                                                            "Text document (*.txt);;"
                                                            "All Files (*)")
        if filename:
            save_gene_set(self.data, filename)
            print(f"Successfully saved at {io.get_datetime()} under {filename}")


class DataFrameView(DataView):
    def __init__(self, data: pd.DataFrame, name: str, parent=None):
        super().__init__(data, name, parent)
        self.label = QtWidgets.QLabel(f"Table '{name}': {self.data.shape[0]} rows, {self.data.shape[1]} columns")

        self.data_view = QtWidgets.QTableView()
        self.save_button = QtWidgets.QPushButton('Save table', self)

        self.init_ui()

    def init_ui(self):
        super().init_ui()
        self.data_view.setModel(DataFrameModel(self.data))

    def save(self):
        default_name = str(self.name) + '.csv'
        filename, _ = QtWidgets.QFileDialog.getSaveFileName(self, "Save table",
                                                            str(Path.home().joinpath(default_name)),
                                                            "Comma-Separated Values (*.csv);;"
                                                            "Tab-Separated Values (*.tsv);;"
                                                            "All Files (*)")
        if filename:
            io.save_csv(self.data, filename)
            print(f"Successfully saved at {io.get_datetime()} under {filename}")


class ErrorMessage(QtWidgets.QDialog):
    def __init__(self, exc_type, exc_value, exc_tb, parent=None):
        super().__init__(parent)
        self.exception = exc_type, exc_value, exc_tb
        self.layout = QtWidgets.QVBoxLayout(self)
        self.widgets = {}
        self.init_ui()

    def init_ui(self):
        self.setWindowTitle("Error")
        self.widgets['error_label'] = QtWidgets.QLabel('<i>RNAlysis</i> has encountered the following error:')
        self.layout.addWidget(self.widgets['error_label'])
        self.setWindowIcon(self.style().standardIcon(QtWidgets.QStyle.SP_MessageBoxCritical))

        tb = "\n".join(traceback.format_exception(*self.exception))
        self.widgets['error_text'] = QtWidgets.QPlainTextEdit(tb)
        self.widgets['error_text'].setReadOnly(True)
        self.layout.addWidget(self.widgets['error_text'])

        self.widgets['ok_button'] = QtWidgets.QPushButton('OK')
        self.widgets['ok_button'].clicked.connect(self.close)
        self.layout.addWidget(self.widgets['ok_button'])

        self.widgets['copy_button'] = QtWidgets.QPushButton('Copy to clipboard')
        self.widgets['copy_button'].clicked.connect(self.copy_to_clipboard)
        self.layout.addWidget(self.widgets['copy_button'])

        self.widgets['copied_label'] = QtWidgets.QLabel()
        self.layout.addWidget(self.widgets['copied_label'])

    def copy_to_clipboard(self):
        cb = QtWidgets.QApplication.clipboard()
        cb.clear(mode=cb.Clipboard)
        cb.setText("".join(traceback.format_exception(*self.exception)), mode=cb.Clipboard)
        self.widgets['copied_label'].setText('Copied to clipboard')


class AboutWindow(QtWidgets.QMessageBox):
    def __init__(self, parent=None):
        super().__init__(parent)
        self.setIcon(QtWidgets.QMessageBox.Information)
        text = f"""<p><b><i>RNAlysis</i> version {__version__}</b>
                </p>
                <br>
                <img src="../../docs/source/logo.png" width="500"/>
                <p>
                Development lead: Guy Teichman (<a href="mailto:guyteichman@gmail.com">guyteichman@gmail.com</a>)
                </p>
                <p>
                Contributors: Or Ganon, Netta Dunsky, Shachar Shani
                </p>"""
        self.setText(text)
        self.setWindowTitle("About RNAlysis")
        self.setStandardButtons(QtWidgets.QMessageBox.Ok)
        self.buttonClicked.connect(self.close)


class ThreadStdOutStreamTextQueueReceiver(QtCore.QObject):
    queue_stdout_element_received_signal = QtCore.pyqtSignal(str)

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
    def __init__(self, parent):
        super(StdOutTextEdit, self).__init__()
        self.setParent(parent)
        self.setReadOnly(True)
        self.setLineWidth(50)
        self.setMinimumWidth(500)
        self.setFont(QtGui.QFont('Consolas', 11))

    @QtCore.pyqtSlot(str)
    def append_text(self, text: str):
        self.moveCursor(QtGui.QTextCursor.End)
        if text == '\n':
            return
        if text.startswith('Warning: '):
            self.insertHtml(f'<div style="color:red;">{text}</div><br>')
        else:
            self.insertHtml(f'<div style="color:black;">{text}</div><br>')
        self.verticalScrollBar().setValue(self.verticalScrollBar().maximum())

        QtWidgets.QApplication.processEvents()


class WriteStream(QtCore.QObject):
    message = QtCore.pyqtSignal(str)

    def __init__(self, q: Queue, parent=None):
        super(WriteStream, self).__init__(parent)

        self.queue = q

    def write(self, text):
        self.queue.put(text)

    def flush(self):
        pass


class SettingsWindow(MinMaxDialog):
    THEMES = {'Light': 'light', 'Dark': 'dark'}
    FONT_SIZES = [str(i) for i in [6, 7, 8, 9, 10, 11, 12, 14, 16, 18, 20, 22, 24, 28, 32, 36, 48, 72]]

    def __init__(self, parent=None):
        super().__init__(parent)
        self.settings_changed: bool = False
        self.layout = QtWidgets.QVBoxLayout(self)
        self.setWindowTitle("Settings")

        self.appearance_group = QtWidgets.QGroupBox('Appearance settings')
        self.appearance_grid = QtWidgets.QGridLayout(self.appearance_group)
        self.appearance_widgets = {}

        self.tables_group = QtWidgets.QGroupBox('User-defined tables')
        self.tables_grid = QtWidgets.QGridLayout(self.tables_group)
        self.tables_widgets = {}

        self.button_box = QtWidgets.QDialogButtonBox(
            QtWidgets.QDialogButtonBox.Ok | QtWidgets.QDialogButtonBox.Cancel |
            QtWidgets.QDialogButtonBox.Apply | QtWidgets.QDialogButtonBox.RestoreDefaults)

        self.layout.addWidget(self.appearance_group)
        self.layout.addWidget(self.tables_group)
        self.layout.addWidget(self.button_box)

        self.init_appearance_ui()
        self.init_tables_ui()
        self.init_buttons()
        self.set_choices()

        self.settings_changed = False

    def _trigger_settings_changed(self):
        self.settings_changed = True

    def set_choices(self):
        current_font, current_font_size, current_theme = settings.get_gui_settings()
        current_theme = {val: key for key, val in self.THEMES.items()}[current_theme]

        self.appearance_widgets['app_theme'].setCurrentText(current_theme)
        self.appearance_widgets['app_font'].setCurrentText(current_font)
        self.appearance_widgets['app_font_size'].setCurrentText(str(current_font_size))

        attr_ref_path = settings.get_attr_ref_path('predefined') if settings.is_setting_in_file(
            settings.__attr_file_key__) else 'No file chosen'
        biotype_ref_path = settings.get_biotype_ref_path('predefined') if settings.is_setting_in_file(
            settings.__biotype_file_key__) else 'No file chosen'

        self.tables_widgets['attr_ref_path'].setText(attr_ref_path)
        self.tables_widgets['biotype_ref_path'].setText(biotype_ref_path)

    def init_appearance_ui(self):
        self.appearance_widgets['app_theme'] = QtWidgets.QComboBox(self.appearance_group)
        self.appearance_widgets['app_theme'].addItems(self.THEMES.keys())

        self.appearance_widgets['app_font'] = QtWidgets.QComboBox(self.appearance_group)
        self.appearance_widgets['app_font'].addItems(list(QtGui.QFontDatabase().families()))
        self.appearance_widgets['app_font'].setEditable(True)
        self.appearance_widgets['app_font'].completer().setCompletionMode(QtWidgets.QCompleter.PopupCompletion)
        self.appearance_widgets['app_font'].setInsertPolicy(QtWidgets.QComboBox.NoInsert)

        self.appearance_widgets['app_font_size'] = QtWidgets.QComboBox(self.appearance_group)
        self.appearance_widgets['app_font_size'].addItems(self.FONT_SIZES)

        for widget_name in ['app_theme', 'app_font', 'app_font_size']:
            self.appearance_widgets[widget_name].currentIndexChanged.connect(self._trigger_settings_changed)

        self.appearance_grid.addWidget(QtWidgets.QLabel('Theme'), 0, 0)
        self.appearance_grid.addWidget(self.appearance_widgets['app_theme'], 0, 1)
        self.appearance_grid.addWidget(QtWidgets.QLabel('Application Font'), 1, 0)
        self.appearance_grid.addWidget(self.appearance_widgets['app_font'], 1, 1)
        self.appearance_grid.addWidget(QtWidgets.QLabel('Application Font Size'), 2, 0)
        self.appearance_grid.addWidget(self.appearance_widgets['app_font_size'], 2, 1)

    def save_settings(self):
        if self.settings_changed:
            self.settings_changed = False
            font = self.appearance_widgets['app_font'].currentText()
            font_size = int(self.appearance_widgets['app_font_size'].currentText())
            theme = self.THEMES[self.appearance_widgets['app_theme'].currentText()]
            settings.set_gui_settings(font, font_size, theme)

            attr_ref_path = self.tables_widgets['attr_ref_path'].text() if self.tables_widgets[
                'attr_ref_path'].is_legal else ''
            biotype_ref_path = self.tables_widgets['biotype_ref_path'].text() if self.tables_widgets[
                'biotype_ref_path'].is_legal else ''

            settings.set_table_settings(attr_ref_path, biotype_ref_path)

            self.parent().update_style_sheet()
            print('Settings saved successfully')

    def reset_settings(self):
        self.settings_changed = False
        settings.reset_settings()
        print("Settings reset successfully")
        self.parent().update_style_sheet()
        self.set_choices()

    def save_and_exit(self):
        self.save_settings()
        self.close()

    def init_tables_ui(self):
        self.tables_widgets['attr_ref_path'] = PathLineEdit()
        self.tables_widgets['biotype_ref_path'] = PathLineEdit()

        for widget in self.tables_widgets.values():
            widget.textChanged.connect(self._trigger_settings_changed)

        self.tables_grid.addWidget(self.tables_widgets['attr_ref_path'], 0, 1)
        self.tables_grid.addWidget(self.tables_widgets['biotype_ref_path'], 1, 1)
        self.tables_grid.addWidget(QtWidgets.QLabel('Attribute Reference Table path:'), 0, 0)
        self.tables_grid.addWidget(QtWidgets.QLabel('Biotype Reference Table path:'), 1, 0)

    def init_buttons(self):
        self.button_box.accepted.connect(self.save_and_exit)
        self.button_box.rejected.connect(self.close)
        self.button_box.clicked.connect(self.handle_button_click)

    def handle_button_click(self, button):
        role = self.button_box.buttonRole(button)
        if role == QtWidgets.QDialogButtonBox.ApplyRole:
            self.save_settings()
        elif role == QtWidgets.QDialogButtonBox.ResetRole:
            self.reset_settings()

    def closeEvent(self, event):
        to_exit = True
        if self.settings_changed:
            quit_msg = "Are you sure you want to close settings without saving?"

            reply = QtWidgets.QMessageBox.question(self, 'Close settings without saving?',
                                                   quit_msg, QtWidgets.QMessageBox.No, QtWidgets.QMessageBox.Yes)
            to_exit = reply == QtWidgets.QMessageBox.Yes

        if to_exit:
            event.accept()
        else:
            event.ignore()


class NewParam:
    def __init__(self, annotation, default=inspect._empty):
        self.annotation = annotation
        self.default = default


def param_to_widget(param, name: str,
                    actions_to_connect: typing.Union[typing.Iterable[typing.Callable], typing.Callable] = tuple()):
    actions_to_connect = parsing.data_to_tuple(actions_to_connect)

    if param.default == inspect._empty:
        is_default = False
    else:
        is_default = True

    if 'color' in name:
        if param.annotation == str:
            widget = ColorPicker(param.default if is_default else '')
            for action in actions_to_connect:
                widget.textChanged.connect(action)

        elif param.annotation in (typing.Tuple[str], typing.Iterable[str], typing.List[str]):
            widget = MultiColorPicker(name)
            widget.set_defaults(param.default if is_default else '')
            for action in actions_to_connect:
                widget.valueChanged.connect(action)

    elif param.annotation == bool:
        widget = QtWidgets.QCheckBox(text=name)
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
    elif param.annotation == float:
        widget = QtWidgets.QDoubleSpinBox()
        widget.setMinimum(float("-inf"))
        widget.setMaximum(float("inf"))
        widget.setSingleStep(0.05)
        default = param.default if is_default else 0.0
        widget.setValue(default)
        for action in actions_to_connect:
            widget.valueChanged.connect(action)
    elif param.annotation == str:
        widget = QtWidgets.QLineEdit(param.default if is_default else '')
        for action in actions_to_connect:
            widget.textChanged.connect(action)
    elif typing_extensions.get_origin(param.annotation) == typing.Union and typing_extensions.Literal in [
        typing_extensions.get_origin(ann) for ann in typing_extensions.get_args(param.annotation)]:
        args = typing_extensions.get_args(param.annotation)
        literal_ind = [typing_extensions.get_origin(ann) for ann in args].index(
            typing_extensions.Literal)
        literal = args[literal_ind]
        without_literal = tuple(args[0:literal_ind] + args[literal_ind + 1:])
        if param.default in typing_extensions.get_args(literal):
            this_default = param.default
            other_default = inspect._empty
        else:
            this_default = None
            other_default = param.default
        widget = ComboBoxOrOtherWidget(typing_extensions.get_args(literal),
                                       param_to_widget(NewParam(typing.Union[without_literal], other_default), name,
                                                       actions_to_connect), this_default)
        for action in actions_to_connect:
            widget.currentIndexChanged.connect(action)
    elif typing_extensions.get_origin(param.annotation) == typing.Literal:
        widget = QtWidgets.QComboBox()
        widget.addItems(typing_extensions.get_args(param.annotation))
        for action in actions_to_connect:
            widget.currentIndexChanged.connect(action)
    elif param.annotation == typing.Union[int, None]:
        widget = OptionalSpinBox()
        widget.setMinimum(-2147483648)
        widget.setMaximum(2147483647)
        default = param.default if is_default else 0
        widget.setValue(default)
        for action in actions_to_connect:
            widget.valueChanged.connect(action)
    elif param.annotation == typing.Union[float, None]:
        widget = OptionalDoubleSpinBox()
        widget.setMinimum(float("-inf"))
        widget.setMaximum(float("inf"))
        widget.setSingleStep(0.05)
        default = param.default if is_default else 0.0
        widget.setValue(default)
        for action in actions_to_connect:
            widget.valueChanged.connect(action)
    elif param.annotation in (typing.Union[str, typing.List[str]], typing.Union[str, typing.Iterable[str]]):
        widget = QMultiLineEdit(name)
        if is_default:
            widget.set_defaults(param.default)
        for action in actions_to_connect:
            widget.valueChanged.connect(action)
    elif param.annotation in (typing.Union[float, typing.List[float]],
                              typing.Union[float, typing.Iterable[float]]):
        widget = QMultiDoubleSpinBox(name)
        if is_default:
            widget.set_defaults(param.default)
        for action in actions_to_connect:
            widget.valueChanged.connect(action)
    elif param.annotation in (typing.Union[int, typing.List[int]], typing.Union[int, typing.Iterable[int]]):
        widget = QMultiSpinBox(name)
        if is_default:
            widget.set_defaults(param.default)
        for action in actions_to_connect:
            widget.valueChanged.connect(action)
    elif param.annotation in (typing.Union[bool, typing.List[bool]], typing.Union[bool, typing.Iterable[bool]]):
        widget = QMultiBoolComboBox(name)
        if is_default:
            widget.set_defaults(param.default)
        for action in actions_to_connect:
            widget.valueChanged.connect(action)
    elif param.annotation == typing.Union[str, int]:
        widget = StrIntLineEdit(param.default if is_default else '')
        for action in actions_to_connect:
            widget.textChanged.connect(action)
    elif param.annotation in (typing.Union[str, int, typing.Iterable[str], typing.Iterable[int]],
                              typing.Union[str, int, typing.List[str], typing.List[int]]):
        widget = QMultiStrIntLineEdit(name)
        widget.set_defaults(param.default if is_default else '')
        for action in actions_to_connect:
            widget.valueChanged.connect(action)
    elif param.annotation in (Path, typing.Union[str, Path], typing.Union[None, str, Path]):
        widget = PathLineEdit()
        for action in actions_to_connect:
            widget.textChanged.connect(action)
    elif param.annotation == typing.Union[str, None]:
        widget = OptionalLineEdit()
        for action in actions_to_connect:
            widget.textChanged.connect(action)
    elif typing_extensions.get_origin(param.annotation) in (
        collections.abc.Iterable, typing.List, typing.Tuple, typing.Set) and typing_extensions.get_origin(
        typing_extensions.get_args(param.annotation)[0]) == typing_extensions.Literal:
        widget = QMultiComboBox(name, items=typing_extensions.get_args(typing_extensions.get_args(param.annotation)[0]))
        if is_default:
            widget.set_defaults(param.default)
        for action in actions_to_connect:
            widget.valueChanged.connect(action)
    # elif param.annotation == typing.Dict[str, typing.List[str]]:
    #     pass
    # elif param.annotation == typing.Dict[str, typing.List[int]]:
    #     pass
    else:
        widget = QtWidgets.QTextEdit()
        default = param.default if is_default else ''
        widget.setText(str(default))
        for action in actions_to_connect:
            widget.textChanged.connect(action)
    return widget


def get_val_from_widget(widget):
    if isinstance(widget, (QtWidgets.QSpinBox, QtWidgets.QDoubleSpinBox)) or (hasattr(widget, 'IS_SPIN_BOX_LIKE')):
        val = widget.value()
    elif isinstance(widget, QtWidgets.QLineEdit) or (hasattr(widget, 'IS_LINE_EDIT_LIKE')):
        val = widget.text()
    elif isinstance(widget, QtWidgets.QCheckBox):
        val = widget.isChecked()
    elif isinstance(widget, QtWidgets.QTextEdit):
        val = ast.literal_eval(widget.toPlainText())
    elif isinstance(widget, QtWidgets.QComboBox) or hasattr(widget, 'IS_COMBO_BOX_LIKE'):
        val = widget.currentText()
    elif hasattr(widget, 'IS_MULTI_INPUT'):
        val = widget.get_values()
    else:
        raise TypeError(f"Invalid QtWidget type {type(widget)}.")
    return val


def save_gene_set(gene_set: set, path):
    with open(path, 'w') as f:
        f.writelines(
            [f"{item}\n" if (i + 1) < len(gene_set) else f"{item}" for i, item in enumerate(gene_set)])


def clear_layout(layout, exceptions: set = frozenset()):
    while layout.count() > len(exceptions):
        child = layout.takeAt(0)
        if child.widget() and child.widget() not in exceptions:
            child.widget().deleteLater()
