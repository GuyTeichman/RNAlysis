import functools
import itertools
import json
import time
import traceback
import warnings
from pathlib import Path
from queue import Queue
from typing import Callable, Union, Tuple

import pandas as pd
import yaml
from PyQt5 import QtCore, QtWidgets, QtGui

from rnalysis import __version__
from rnalysis.gui import gui_style, gui_widgets
from rnalysis.utils import settings, io, generic


class CheckableFileSystemModel(QtWidgets.QFileSystemModel):
    checkStateChanged = QtCore.pyqtSignal(str, bool)
    finishedDataChange = QtCore.pyqtSignal()

    def __init__(self):
        super().__init__()
        self.checkStates = {}
        self.rowsInserted.connect(self.checkAdded)
        self.rowsRemoved.connect(self.checkParent)
        self.rowsAboutToBeRemoved.connect(self.checkRemoved)

    def checkState(self, index):
        return self.checkStates.get(self.filePath(index), QtCore.Qt.Unchecked)

    def setCheckState(self, index, state, emitStateChange=True):
        path = self.filePath(index)
        if self.checkStates.get(path) == state:
            return
        self.checkStates[path] = state
        if emitStateChange:
            self.checkStateChanged.emit(path, bool(state))

    def checkAdded(self, parent, first, last):
        # if a file/directory is added, ensure it follows the parent state as long
        # as the parent is already tracked; note that this happens also when
        # expanding a directory that has not been previously loaded
        if not parent.isValid():
            return
        if self.filePath(parent) in self.checkStates:
            state = self.checkState(parent)
            for row in range(first, last + 1):
                index = self.index(row, 0, parent)
                path = self.filePath(index)
                if path not in self.checkStates:
                    self.checkStates[path] = state
        self.checkParent(parent)

    def checkRemoved(self, parent, first, last):
        # remove items from the internal dictionary when a file is deleted;
        # note that this *has* to happen *before* the model actually updates,
        # that's the reason this function is connected to rowsAboutToBeRemoved
        for row in range(first, last + 1):
            path = self.filePath(self.index(row, 0, parent))
            if path in self.checkStates:
                self.checkStates.pop(path)

    def checkParent(self, parent):
        # verify the state of the parent according to the children states
        if not parent.isValid():
            self.finishedDataChange.emit()
            return
        childStates = [self.checkState(self.index(r, 0, parent)) for r in range(self.rowCount(parent))]
        newState = QtCore.Qt.Checked if all(childStates) else QtCore.Qt.Unchecked
        oldState = self.checkState(parent)
        if newState != oldState:
            self.setCheckState(parent, newState)
            self.dataChanged.emit(parent, parent)
        self.checkParent(parent.parent())

    def flags(self, index):
        return super().flags(index) | QtCore.Qt.ItemIsUserCheckable

    def data(self, index, role=QtCore.Qt.DisplayRole):
        if role == QtCore.Qt.CheckStateRole and index.column() == 0:
            return self.checkState(index)
        return super().data(index, role)

    def setData(self, index, value, role, checkParent=True, emitStateChange=True):
        if role == QtCore.Qt.CheckStateRole and index.column() == 0:
            self.setCheckState(index, value, emitStateChange)
            for row in range(self.rowCount(index)):
                # set the data for the children, but do not emit the state change,
                # and don't check the parent state (to avoid recursion)
                self.setData(index.child(row, 0), value, QtCore.Qt.CheckStateRole,
                             checkParent=False, emitStateChange=False)
            self.dataChanged.emit(index, index)
            if checkParent:
                self.checkParent(index.parent())
            return True

        return super().setData(index, value, role)


class FilterProxy(QtCore.QSortFilterProxyModel):
    '''
    Based on StackOverflow answer by user ekhumoro:
    https://stackoverflow.com/questions/72587813/how-to-filter-by-no-extension-in-qfilesystemmodel
    '''

    def __init__(self, disables=False, parent=None):
        super().__init__(parent)
        self._disables = bool(disables)

    def filterAcceptsRow(self, row, parent):
        index = self.sourceModel().index(row, 0, parent)
        if not self._disables:
            return self.matchIndex(index)
        return index.isValid()

    def matchIndex(self, index):
        return (super().filterAcceptsRow(index.row(), index.parent()))

    def flags(self, index):
        flags = super().flags(index)
        if (self._disables and
            not self.matchIndex(self.mapToSource(index))):
            flags &= ~QtCore.Qt.ItemIsEnabled
        return flags


class MultiFileSelectionDialog(gui_widgets.MinMaxDialog):
    '''
    Based on a Stack Overflow answer by user 'musicamante':
    https://stackoverflow.com/questions/63309406/qfilesystemmodel-with-checkboxes
    '''

    def __init__(self):
        super().__init__()
        self.layout = QtWidgets.QGridLayout(self)
        self.tree_mycomputer = QtWidgets.QTreeView()
        self.tree_home = QtWidgets.QTreeView()
        self.open_button = QtWidgets.QPushButton('Open')
        self.cancel_button = QtWidgets.QPushButton('Cancel')
        self.logger = QtWidgets.QPlainTextEdit()
        self.init_models()
        self.init_ui()

    def init_models(self):
        model_mycomputer = CheckableFileSystemModel()
        model_mycomputer.setRootPath('')
        model_home = CheckableFileSystemModel()
        model_home.setRootPath('.')
        proxy = FilterProxy(False, self)
        proxy.setFilterRegularExpression(r'^(?![.])(?!.*[-_.]$).+')
        proxy.setSourceModel(model_home)
        self.tree_mycomputer.setModel(model_mycomputer)
        self.tree_mycomputer.setRootIndex(model_mycomputer.index(model_mycomputer.myComputer()))
        self.tree_home.setModel(proxy)
        self.tree_home.setRootIndex(proxy.mapFromSource(
            model_home.index(QtCore.QStandardPaths.standardLocations(QtCore.QStandardPaths.HomeLocation)[0])))
        model_mycomputer.finishedDataChange.connect(self.update_log)
        model_home.finishedDataChange.connect(self.update_log)

    def init_ui(self):
        self.setWindowTitle('Choose files:')
        self.resize(1000, 750)

        self.tree_mycomputer.setSortingEnabled(True)
        self.tree_home.setSortingEnabled(True)
        self.tree_mycomputer.header().setSectionResizeMode(QtWidgets.QHeaderView.Stretch)
        self.tree_home.header().setSectionResizeMode(QtWidgets.QHeaderView.Stretch)

        self.layout.addWidget(self.tree_mycomputer, 0, 0, 4, 8)
        self.layout.addWidget(self.tree_home, 4, 0, 4, 8)
        self.layout.setRowStretch(3, 2)
        self.layout.setRowStretch(7, 2)

        self.layout.addWidget(self.open_button, 8, 7)
        self.layout.addWidget(self.cancel_button, 9, 7)

        self.layout.addWidget(self.logger, 8, 0, 2, 7)
        self.logger.setReadOnly(True)

        self.open_button.clicked.connect(self.accept)
        self.cancel_button.clicked.connect(self.reject)

        self.update_log()

    def update_log(self):
        self.logger.setPlainText("\n".join(self.result()))
        self.logger.verticalScrollBar().setValue(
            self.logger.verticalScrollBar().maximum())

    def result(self):
        files_to_open = []
        queue = Queue()
        for pth, state in itertools.chain(self.tree_mycomputer.model().checkStates.items(),
                                          self.tree_home.model().sourceModel().checkStates.items()):
            if state:
                queue.put(pth)

        while not queue.empty():
            this_path = Path(queue.get())
            if this_path.is_file():
                files_to_open.append(str(this_path))
            else:
                try:
                    for item in this_path.iterdir():
                        queue.put(item)
                except PermissionError:
                    warnings.warn(f'Cannot access items under {this_path} - permission denied. ')
        return files_to_open


class DataFrameModel(QtCore.QAbstractTableModel):
    """
    Basde upon:
    https://stackoverflow.com/a/44605011
    """
    DtypeRole = QtCore.Qt.UserRole + 1000
    ValueRole = QtCore.Qt.UserRole + 1001

    def __init__(self, df=pd.DataFrame(), parent=None):
        super().__init__(parent)
        if isinstance(df, pd.Series):
            df = df.to_frame()
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
        shape = df.shape
        if len(shape) == 1:
            shape = (shape[0], 1)

        n_rows = min(2, shape[0])
        n_cols = min(3, shape[1])
        if isinstance(df, pd.DataFrame):
            df = df.iloc[:n_rows, :n_cols]
        elif isinstance(df, pd.Series):
            df = df.iloc[:n_rows]

        if n_rows < shape[0]:
            if isinstance(df, pd.DataFrame):
                df.loc['...', :] = '...'
            elif isinstance(df, pd.Series):
                df.loc['...'] = '...'
        if n_cols < shape[1]:
            df['...'] = '...'
        super().__init__(df, parent)


class DataView(gui_widgets.MinMaxDialog):
    def __init__(self, data, name: str, parent=None):
        super().__init__(parent)
        self.data = data
        self.name = name
        self.label_font = QtGui.QFont()
        self.layout = QtWidgets.QVBoxLayout(self)

    def init_ui(self):
        self.setGeometry(1000, 200, 1000, 500)
        self.setWindowTitle(f"View of '{self.name}'")
        self.label_font.setPointSize(12)
        self.label.setFont(self.label_font)
        self.layout.addWidget(self.label)
        self.layout.addWidget(self.save_button)
        self.layout.addWidget(self.data_view)
        self.save_button.clicked.connect(self.save)
        self.setLayout(self.layout)
        self.setStyleSheet(gui_style.get_stylesheet())

    def closeEvent(self, a0: QtGui.QCloseEvent) -> None:
        self.deleteLater()

    def save(self):
        raise NotImplementedError


class GeneSetView(DataView):
    def __init__(self, data: set, name: str, parent=None):
        super().__init__(data, name, parent)
        self.label = QtWidgets.QLabel(f"Gene set '{name}': {len(self.data)} features")

        self.data_view = gui_widgets.ReactiveListWidget()
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
            io.save_gene_set(self.data, filename)
            print(f"Successfully saved at {io.get_datetime()} under {filename}")


class DataFrameView(DataView):
    def __init__(self, data: pd.DataFrame, name: str, parent=None):
        super().__init__(data, name, parent)
        shape = self.data.shape
        if len(shape) == 1:
            shape = (shape[0], 1)
        self.label = QtWidgets.QLabel(f"Table '{name}': {shape[0]} rows, {shape[1]} columns")

        self.data_view = gui_widgets.ReactiveTableView()
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
        self.setWindowIcon(self.style().standardIcon(QtWidgets.QStyle.SP_MessageBoxCritical))

        self.widgets['error_label'] = QtWidgets.QLabel('<i>RNAlysis</i> has encountered the following error:')
        self.layout.addWidget(self.widgets['error_label'])

        self.widgets['error_summary'] = QtWidgets.QLabel(f'<b>{";".join(self.exception[1].args)}</b>')
        self.widgets['error_summary'].setTextInteractionFlags(QtCore.Qt.TextInteractionFlag.TextSelectableByMouse)
        self.widgets['error_summary'].setWordWrap(True)
        self.layout.addWidget(self.widgets['error_summary'])

        self.layout.addSpacing(3)
        self.widgets['full_text_label'] = QtWidgets.QLabel('Full error report:')
        self.layout.addWidget(self.widgets['full_text_label'])

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
        img_path = str(Path(__file__).parent.joinpath('splash_transparent.png'))
        text = f"""<br>
                <p align="center"><b><i><b>RNAlysis</i>/<b> version {__version__}</b>
                </p>
                <br>
                <img src="{img_path}" width="500"/>
                <p>
                Development lead: Guy Teichman (<a href="mailto:guyteichman@gmail.com">guyteichman@gmail.com</a>)
                </p>
                <p>
                Contributors: Dror Cohen, Or Ganon, Netta Dunsky, Shachar Shani
                </p>"""
        self.setText(text)
        self.setWindowTitle("About RNAlysis")
        self.setStandardButtons(QtWidgets.QMessageBox.Ok)
        self.buttonClicked.connect(self.close)


class SettingsWindow(gui_widgets.MinMaxDialog):
    LOOKUP_DATABASES_PATH = Path(__file__).parent.parent.joinpath('data_files/lookup_databases.json')
    THEMES = gui_style.get_stylesheet_names()
    FONT_SIZES = [str(i) for i in [6, 7, 8, 9, 10, 11, 12, 14, 16, 18, 20, 22, 24, 28, 32, 36, 48, 72]]
    REPORT_GEN_OPTIONS = {'Ask me every time': None, 'Always': True, 'Never': False}
    styleSheetUpdated = QtCore.pyqtSignal()

    def __init__(self, parent=None):
        super().__init__(parent)
        self.settings_changed: bool = False
        self.layout = QtWidgets.QVBoxLayout(self)
        self.setWindowTitle("Settings")

        self.appearance_group = QtWidgets.QGroupBox('User Interface settings')
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

    def exec(self):
        self.set_choices()
        super().exec()

    def set_choices(self):
        current_font, current_font_size, current_theme, current_dbs, current_show_tutorial, current_report_gen = \
            settings.get_gui_settings()
        current_theme = {val: key for key, val in self.THEMES.items()}[current_theme]

        self.appearance_widgets['app_theme'].setCurrentText(current_theme)
        self.appearance_widgets['app_font'].setCurrentText(current_font)
        self.appearance_widgets['app_font_size'].setCurrentText(str(current_font_size))
        self.appearance_widgets['show_tutorial'].setChecked(current_show_tutorial)
        self.appearance_widgets['report_gen'].setCurrentText(
            {val: key for key, val in self.REPORT_GEN_OPTIONS.items()}[current_report_gen])

        for i in range(self.appearance_widgets['databases'].count()):
            item = self.appearance_widgets['databases'].item(i)
            if item.text() in current_dbs:
                item.setCheckState(QtCore.Qt.Checked)

        attr_ref_path = settings.get_attr_ref_path('predefined') if settings.is_setting_in_file(
            settings.__attr_file_key__) else 'No file chosen'
        biotype_ref_path = settings.get_biotype_ref_path('predefined') if settings.is_setting_in_file(
            settings.__biotype_file_key__) else 'No file chosen'

        self.tables_widgets['attr_ref_path'].setText(attr_ref_path)
        self.tables_widgets['biotype_ref_path'].setText(biotype_ref_path)

    def init_appearance_ui(self):
        self.appearance_widgets['app_theme'] = QtWidgets.QComboBox(self.appearance_group)
        self.appearance_widgets['app_theme'].addItems(self.THEMES.keys())

        self.appearance_widgets['app_font'] = QtWidgets.QFontComboBox(self.appearance_group)
        self.appearance_widgets['app_font'].setFontFilters(QtWidgets.QFontComboBox.ScalableFonts)
        self.appearance_widgets['app_font'].setEditable(True)
        self.appearance_widgets['app_font'].completer().setCompletionMode(QtWidgets.QCompleter.PopupCompletion)
        self.appearance_widgets['app_font'].setInsertPolicy(QtWidgets.QComboBox.NoInsert)

        self.appearance_widgets['app_font_size'] = QtWidgets.QComboBox(self.appearance_group)
        self.appearance_widgets['app_font_size'].addItems(self.FONT_SIZES)

        self.appearance_widgets['databases'] = QtWidgets.QListWidget(self)
        with open(self.LOOKUP_DATABASES_PATH) as f:
            for key in json.load(f).keys():
                item = QtWidgets.QListWidgetItem(key)
                item.setCheckState(QtCore.Qt.Unchecked)
                self.appearance_widgets['databases'].addItem(item)

        self.appearance_widgets['show_tutorial'] = QtWidgets.QCheckBox("Show tutorial page on startup")
        self.appearance_widgets['report_gen'] = QtWidgets.QComboBox()
        self.appearance_widgets['report_gen'].addItems(self.REPORT_GEN_OPTIONS.keys())
        self.appearance_widgets['report_gen'].setInsertPolicy(QtWidgets.QComboBox.NoInsert)

        for widget_name in ['app_theme', 'app_font', 'app_font_size', 'report_gen']:
            self.appearance_widgets[widget_name].currentIndexChanged.connect(self._trigger_settings_changed)
        self.appearance_widgets['show_tutorial'].stateChanged.connect(self._trigger_settings_changed)

        self.appearance_grid.addWidget(QtWidgets.QLabel('Theme'), 0, 0)
        self.appearance_grid.addWidget(self.appearance_widgets['app_theme'], 0, 1)
        self.appearance_grid.addWidget(QtWidgets.QLabel('Application Font'), 1, 0)
        self.appearance_grid.addWidget(self.appearance_widgets['app_font'], 1, 1)
        self.appearance_grid.addWidget(QtWidgets.QLabel('Application Font Size'), 2, 0)
        self.appearance_grid.addWidget(self.appearance_widgets['app_font_size'], 2, 1)
        self.appearance_grid.addWidget(QtWidgets.QLabel('Right-click databases'), 3, 0)
        self.appearance_grid.addWidget(self.appearance_widgets['databases'], 3, 1)
        self.appearance_grid.addWidget(self.appearance_widgets['show_tutorial'], 4, 0, 1, 2)
        self.appearance_grid.addWidget(QtWidgets.QLabel('Turn on report generation automatically'), 5, 0)
        self.appearance_grid.addWidget(self.appearance_widgets['report_gen'], 5, 1)

    def save_settings(self):
        if self.settings_changed:
            self.settings_changed = False
            font = self.appearance_widgets['app_font'].currentText()
            font_size = int(self.appearance_widgets['app_font_size'].currentText())
            theme = self.THEMES[self.appearance_widgets['app_theme'].currentText()]
            show_tutorial = self.appearance_widgets['show_tutorial'].isChecked()
            prompt_report_gen = self.REPORT_GEN_OPTIONS[self.appearance_widgets['report_gen'].currentText()]
            dbs = []
            for i in range(self.appearance_widgets['databases'].count()):
                item = self.appearance_widgets['databases'].item(i)
                if item.checkState():
                    dbs.append(item.text())

            settings.set_gui_settings(font, font_size, theme, dbs, show_tutorial, prompt_report_gen)

            attr_ref_path = self.tables_widgets['attr_ref_path'].text() if self.tables_widgets[
                'attr_ref_path'].is_legal else ''
            biotype_ref_path = self.tables_widgets['biotype_ref_path'].text() if self.tables_widgets[
                'biotype_ref_path'].is_legal else ''

            settings.set_table_settings(attr_ref_path, biotype_ref_path)

            self.styleSheetUpdated.emit()
            print('Settings saved successfully')

    def reset_settings(self):
        self.settings_changed = False
        settings.reset_settings()
        self.styleSheetUpdated.emit()
        self.set_choices()
        print("Settings reset successfully")

    def save_and_exit(self):
        self.save_settings()
        self.close()

    def init_tables_ui(self):
        self.tables_widgets['attr_ref_path'] = gui_widgets.PathLineEdit()
        self.tables_widgets['biotype_ref_path'] = gui_widgets.PathLineEdit()

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


class HowToCiteWindow(gui_widgets.MinMaxDialog):
    CITATION_RNALYSIS = """
    Teichman, G., Cohen, D., Ganon, O., Dunsky, N., Shani, S., Gingold, H., and Rechavi, O. (2022).
    RNAlysis: analyze your RNA sequencing data without writing a single line of code. BioRxiv 2022.11.25.517851.
    <br>
    <a href=https://doi.org/10.1101/2022.11.25.517851>doi.org/10.1101/2022.11.25.517851</a>
    """
    CITATION_CUTADAPT = """
    Martin, M. (2011). Cutadapt removes adapter sequences from high-throughput sequencing reads.
    EMBnet.journal, 17(1), pp. 10-12.
    <br>
    <a href=https://doi.org/10.14806/ej.17.1.200>doi.org/10.14806/ej.17.1.200</a>
    """
    CITATION_KALLISTO = """
    Bray, N., Pimentel, H., Melsted, P. et al.
    Near-optimal probabilistic RNA-seq quantification.
    Nat Biotechnol 34, 525–527 (2016).
    <br>
    <a href=https://doi.org/10.1038/nbt.3519>doi.org/10.1038/nbt.3519</a>
    """
    CITATION_DESEQ2 = """
    Love MI, Huber W, Anders S (2014).
    “Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2.”
    Genome Biology, 15, 550.
    <br>
    <a href=https://doi.org/10.1186/s13059-014-0550-8>doi.org/10.1186/s13059-014-0550-8</a>
    """
    CITATION_HDBSCAN = """
    L. McInnes, J. Healy, S. Astels, hdbscan:
    Hierarchical density based clustering In:
    Journal of Open Source Software, The Open Journal, volume 2, number 11. 2017
    <br>
    <a href=https://doi.org/10.1371/journal.pcbi.0030039>doi.org/10.1371/journal.pcbi.0030039</a>"""
    CITATION_XLMHG = """
    <p>
    Eden, E., Lipson, D., Yogev, S., and Yakhini, Z. (2007).
     Discovering Motifs in Ranked Lists of DNA Sequences. PLOS Comput. Biol. 3, e39.
    <br>
    <a href=https://doi.org/10.1371/journal.pcbi.0030039>doi.org/10.1371/journal.pcbi.0030039</a>
    </p>
    <p>
    Wagner, F. (2017). The XL-mHG test for gene set enrichment. ArXiv.
    <br>
    <a href=https://doi.org/10.48550/arXiv.1507.07905>doi.org/10.48550/arXiv.1507.07905</a>
    </p>"""
    CITATION_FILE_PATH = Path(__file__).parent.parent.joinpath('data_files/tool_citations.json')

    def __init__(self, parent=None):
        super().__init__(parent)
        img_path = str(Path(__file__).parent.joinpath('logo_small.png'))
        text = f"""<p align="center"><b><i>RNAlysis</i> version {__version__}</b>
                </p>
                <br>
                <img src="{img_path}" width="250"/>"""
        self.label = QtWidgets.QLabel(text)

        self.citation_labels = []
        self.citations = []
        with open(self.CITATION_FILE_PATH, encoding='utf-8') as f:
            citation_dict = json.load(f)
        for data in citation_dict.values():
            txt = f"If you use {data['name']} in your research, please cite:"
            self.citation_labels.append(QtWidgets.QLabel(txt))
            self.citations.append(gui_widgets.TextWithCopyButton(data['citation']))

        self.ok_button = QtWidgets.QPushButton('Close')

        self.main_layout = QtWidgets.QVBoxLayout(self)
        self.scroll = QtWidgets.QScrollArea()
        self.scroll_widget = QtWidgets.QWidget(self.scroll)
        self.layout = QtWidgets.QVBoxLayout(self.scroll_widget)
        self.init_ui()

    def init_ui(self):
        self.scroll.setHorizontalScrollBarPolicy(QtCore.Qt.ScrollBarAlwaysOff)
        self.scroll.setWidgetResizable(True)
        self.scroll.setWidget(self.scroll_widget)
        self.layout.setSizeConstraint(QtWidgets.QLayout.SetMinAndMaxSize)

        self.main_layout.addWidget(self.scroll)

        self.setWindowTitle("How to cite RNAlysis")
        self.layout.addWidget(self.label)

        for label, citation in zip(self.citation_labels, self.citations):
            self.layout.addWidget(label)
            self.layout.addWidget(citation)

        self.ok_button.clicked.connect(self.close)
        self.layout.addWidget(self.ok_button)


def splash_screen():
    img_path = str(Path(__file__).parent.joinpath('splash.png'))
    splash_pixmap = QtGui.QPixmap(img_path)
    splash_font = QtGui.QFont('Calibri', 16)
    splash = QtWidgets.QSplashScreen(splash_pixmap)
    splash.setFont(splash_font)
    splash.showMessage(f"<i>RNAlysis</i> version {__version__}", QtCore.Qt.AlignBottom | QtCore.Qt.AlignHCenter)
    splash.show()
    return splash


class FuncExternalWindow(gui_widgets.MinMaxDialog):
    IGNORED_WIDGETS = {'help_link'}
    paramsAccepted = QtCore.pyqtSignal(list, dict, object)
    geneSetsRequested = QtCore.pyqtSignal(object)
    __slots__ = {'func_name': 'name of the function to be applied',
                 'func': 'function to be applied',
                 'signature': 'signature of the function',
                 'desc': 'description of the function',
                 'param_desc': 'description of the function parameters',
                 'help_link': 'link to the documentation page of the function',
                 'excluded_params': 'parameters to be excluded from the window',
                 'scroll': 'scroll area',
                 'scroll_widget': 'widget containing the scroll area',
                 'scroll_layout': 'layout for the scroll widget',
                 'param_group': 'widget group for the parameter widgets',
                 'param_grid': 'layout for the parameter widgets',
                 'param_widgets': 'parameter widgets',
                 'start_button': 'start button',
                 'import_button': 'import button',
                 'export_button': 'export button',
                 'close_button': 'close button',
                 'args': 'function args'}

    def __init__(self, func_name: str, func: Callable, help_link: Union[None, str], excluded_params: set, threaded=True,
                 parent=None):
        super().__init__(parent)
        self.func_name = func_name
        self.func = func
        self.signature = generic.get_method_signature(self.func)
        self.desc, self.param_desc = io.get_method_docstring(self.func)
        self.help_link = help_link
        self.excluded_params = excluded_params.copy()
        self.threaded = threaded
        self.args = []

        self.widgets = {}
        self.main_layout = QtWidgets.QVBoxLayout(self)

        self.scroll = QtWidgets.QScrollArea()
        self.scroll_widget = QtWidgets.QWidget(self.scroll)
        self.scroll_layout = QtWidgets.QGridLayout(self.scroll_widget)

        self.param_group = QtWidgets.QGroupBox(f"1. Set {func_name} parameters")
        self.param_grid = QtWidgets.QGridLayout(self.param_group)
        self.param_widgets = {}

        self.start_button = QtWidgets.QPushButton(f'Start {self.func_name}')
        self.import_button = QtWidgets.QPushButton('Import parameters')
        self.export_button = QtWidgets.QPushButton('Export parameters')
        self.close_button = QtWidgets.QPushButton('Close')

    def init_ui(self):
        self.scroll.setHorizontalScrollBarPolicy(QtCore.Qt.ScrollBarAlwaysOff)
        self.scroll.setWidgetResizable(True)
        self.scroll.setWidget(self.scroll_widget)
        self.scroll_layout.setSizeConstraint(QtWidgets.QLayout.SetMinAndMaxSize)

        self.main_layout.addWidget(self.scroll)

        if self.help_link is not None:
            self.param_widgets['help_link'] = QtWidgets.QLabel(
                text=f'<a href="{self.help_link}">Open documentation for <b>{self.func_name}</b></a>')
            self.param_widgets['help_link'].setOpenExternalLinks(True)
            self.main_layout.addWidget(self.param_widgets['help_link'])

        self.scroll_layout.addWidget(self.param_group, 0, 0)
        self.scroll_layout.addWidget(self.start_button, 1, 0, 1, 2)
        self.scroll_layout.addWidget(self.import_button, 2, 0, 1, 2)
        self.scroll_layout.addWidget(self.export_button, 3, 0, 1, 2)
        self.scroll_layout.addWidget(self.close_button, 4, 0, 1, 2)

        if self.threaded:
            self.start_button.clicked.connect(self.run_function_threaded)
        else:
            self.start_button.clicked.connect(self.run_function_in_main_loop)

        self.import_button.clicked.connect(self.import_parameters)
        self.export_button.clicked.connect(self.export_parameters)
        self.close_button.clicked.connect(self.close)

        self.init_param_ui()

    def init_param_ui(self):
        i = 0

        for name, param in self.signature.items():
            if name in self.excluded_params:
                continue
            this_desc = self.param_desc.get(name, '')
            self.param_widgets[name] = gui_widgets.param_to_widget(param, name)
            self.connect_widget(self.param_widgets[name])

            label = QtWidgets.QLabel(f'{name}:', self.param_widgets[name])
            label.setToolTip(this_desc)
            help_button = gui_widgets.HelpButton()
            self.param_grid.addWidget(help_button, i, 2)
            self.param_grid.addWidget(label, i, 0)
            self.param_grid.addWidget(self.param_widgets[name], i, 1)
            help_button.connect_param_help(name, this_desc)
            i += 1

        self.param_grid.setRowStretch(i, 1)

    def connect_widget(self, widget: QtWidgets.QWidget):
        if isinstance(widget, (gui_widgets.ComboBoxOrOtherWidget, gui_widgets.OptionalWidget)):
            self.connect_widget(widget.other)
        elif isinstance(widget, gui_widgets.GeneSetComboBox):
            widget.boxOpened.connect(functools.partial(self.geneSetsRequested.emit, widget))

    def get_analysis_kwargs(self):
        kwargs = {}
        for param_name, widget in self.param_widgets.items():
            if param_name in self.IGNORED_WIDGETS:
                continue
            val = gui_widgets.get_val_from_widget(widget)
            kwargs[param_name] = val
        return kwargs

    def get_analysis_args(self):
        args = self.args
        return args

    def import_parameters(self):
        filename, _ = QtWidgets.QFileDialog.getOpenFileName(self, "Choose a parameter file",
                                                            filter="YAML file (*.yaml)")
        if not filename:
            return
        with open(filename) as f:
            params = yaml.safe_load(f)

        for key in ['name', 'args', 'kwargs']:
            assert key in params, f"Invalid parameter file: key '{key}' missing."
        assert params.get('name') == self.func_name, f"Parameter file for function '{params.get('name')}' " \
                                                     f"does not match this window's function: '{self.func_name}'"
        args = params['args']
        kwargs = params['kwargs']

        self.args = args

        for key, val in kwargs.items():
            widget = self.param_widgets.get(key)
            if widget is not None:
                gui_widgets.set_widget_value(widget, val)

        return params

    def export_parameters(self):
        args = self.get_analysis_args()
        kwargs = self.get_analysis_kwargs()
        parameter_dict = {'name': self.func_name, 'args': args, 'kwargs': kwargs}

        default_name = f'parameters {self.func_name}.yaml'
        filename, _ = QtWidgets.QFileDialog.getSaveFileName(self, "Export parameter file",
                                                            str(Path.home().joinpath(default_name)),
                                                            "YAML file (*.yaml)")
        if filename:
            with open(filename, 'w') as f:
                yaml.safe_dump(parameter_dict, f)
                print(f"Successfully saved at {io.get_datetime()} under {filename}")

    def run_function_threaded(self):
        args = self.get_analysis_args()
        kwargs = self.get_analysis_kwargs()
        try:
            self.paramsAccepted.emit(args, kwargs, self.showNormal)
            self.showMinimized()
        except Exception as e:
            self.showNormal()
            raise e

    def run_function_in_main_loop(self):
        args = self.get_analysis_args()
        kwargs = self.get_analysis_kwargs()
        self.func(*args, **kwargs)


class PairedFuncExternalWindow(FuncExternalWindow):
    __slots__ = {'pairs_group': 'widget group for picking file pairs',
                 'pairs_grid': 'layout for widget group',
                 'pairs_widgets': 'widgets for picking file pairs'}
    EXCLUDED_PARAMS = {'r1_files', 'r2_files'}

    def __init__(self, func_name: str, func: Callable, help_link: Union[str, None], excluded_params: set, parent=None):
        super().__init__(func_name, func, help_link, excluded_params, parent=parent)

        self.pairs_group = QtWidgets.QGroupBox("2. Choose FASTQ file pairs")
        self.pairs_grid = QtWidgets.QGridLayout(self.pairs_group)
        self.pairs_widgets = {}

    def init_ui(self):
        self.scroll_layout.addWidget(self.pairs_group, 0, 1)
        self.setMinimumSize(1250, 650)
        super().init_ui()
        self.init_pairs_ui()

    def init_pairs_ui(self):
        self.pairs_widgets['r1_files'] = gui_widgets.OrderedFileList(self)
        self.pairs_widgets['r2_files'] = gui_widgets.OrderedFileList(self)

        self.pairs_grid.addWidget(self.pairs_widgets['r1_files'], 1, 0)
        self.pairs_grid.addWidget(self.pairs_widgets['r2_files'], 1, 1)
        self.pairs_grid.addWidget(QtWidgets.QLabel("<b>R1 files:</b>"), 0, 0)
        self.pairs_grid.addWidget(QtWidgets.QLabel("<b>R2 files:</b>"), 0, 1)

    def get_analysis_kwargs(self):
        kwargs = super().get_analysis_kwargs()
        kwargs['r1_files'] = self.pairs_widgets['r1_files'].get_sorted_names()
        kwargs['r2_files'] = self.pairs_widgets['r2_files'].get_sorted_names()
        return kwargs

    def import_parameters(self):
        params = super().import_parameters()
        for key in ['r1_files', 'r2_files']:
            assert key in params['kwargs'], f"cannot find value for key '{key}'"
            value = params['kwargs'][key]
            gui_widgets.set_widget_value(self.pairs_widgets[key], value)


class StatusBar(QtWidgets.QStatusBar):
    taskQueueRequested = QtCore.pyqtSignal()

    def __init__(self, parent=None):
        super().__init__(parent)
        self.progbar_desc = ''
        self.progbar_total = 0
        self.progbar_start_time = time.time()
        self.progbar_completed_items = 0
        self._is_running = False

        self.n_tasks_button = QtWidgets.QPushButton(self)
        self.desc_label = QtWidgets.QLabel(self)
        self.progress_bar = QtWidgets.QProgressBar(self)
        self.elapsed_label = QtWidgets.QLabel(self)
        self.remaining_label = QtWidgets.QLabel(self)
        self.update_timer = QtCore.QTimer(self)

        self.init_ui()

    def init_ui(self):
        self.n_tasks_button.clicked.connect(self.taskQueueRequested.emit)

        self.addWidget(self.n_tasks_button)
        self.addWidget(self.desc_label)
        self.addWidget(self.progress_bar)
        self.addWidget(self.elapsed_label)
        self.addWidget(self.remaining_label)

        self.update_n_tasks(0)
        self.reset_progress()

        self.update_timer.timeout.connect(self.update_time)
        self.update_timer.start(2000)

    def update_n_tasks(self, n_tasks: int):
        if n_tasks <= 0:
            self._is_running = False
            self.n_tasks_button.setVisible(False)
        else:
            self.n_tasks_button.setVisible(True)
            self.n_tasks_button.setText(f'{n_tasks} tasks running... ')

    def update_desc(self, desc: str):
        self.desc_label.setText(f'{desc}:')
        self.desc_label.setVisible(True)

    def update_time(self):
        if not self._is_running:
            return

        elapsed_time = time.time() - self.progbar_start_time
        if self.progbar_completed_items == 0:
            remaining_time = elapsed_time * self.progbar_total
        else:
            remaining_time = (elapsed_time / self.progbar_completed_items) * abs(
                self.progbar_total - self.progbar_completed_items)
        self.elapsed_label.setText(f"{generic.format_time(elapsed_time)} elapsed ")
        self.remaining_label.setText(f"{generic.format_time(remaining_time)} remaining ")
        self.elapsed_label.setVisible(True)
        self.remaining_label.setVisible(True)

    def reset_progress(self):
        self._is_running = False
        self.desc_label.setVisible(False)
        self.progress_bar.setVisible(False)
        self.elapsed_label.setVisible(False)
        self.remaining_label.setVisible(False)

    def start_progress(self, total: int, description: str):
        self._is_running = True
        self.progbar_start_time = time.time()
        self.progbar_completed_items = 0

        self.update_bar_total(total)
        self.progress_bar.setValue(0)
        self.update_desc(description)
        self.progress_bar.setVisible(True)
        self.move_progress_bar(0)

    def update_bar_total(self, total: int):
        self.progbar_total = total
        self.progress_bar.setRange(0, total)

    def move_progress_bar(self, value: int):
        self.progbar_completed_items += value
        self.update_time()
        self.progress_bar.setValue(self.progbar_completed_items)


class TaskQueueWindow(gui_widgets.MinMaxDialog):
    cancelRequested = QtCore.pyqtSignal(int, str)

    def __init__(self, parent=None):
        super().__init__(parent)
        self.tasks = []
        self.main_layout = QtWidgets.QVBoxLayout(self)
        self.list = gui_widgets.MultiChoiceListWithDelete(self.tasks, parent=self, delete_text='cancel')

        self.main_layout.addWidget(self.list)
        self.setWindowTitle('Task queue')

        self.list.itemDeleted.connect(self.request_cancel)

    def update_tasks(self, tasks: list):
        if tasks == self.tasks:
            return
        self.tasks = tasks
        self.list.delete_all_quietly()
        self.list.add_items(self.tasks)

    @QtCore.pyqtSlot(int)
    def request_cancel(self, index: int):
        name = self.tasks.pop(index)
        self.cancelRequested.emit(index, name)


class ApplyTablePipelineWindow(gui_widgets.MinMaxDialog):
    __slots__ = {'available_objects': 'available objects',
                 'layout': 'layout',
                 'label': 'main label of the window',
                 'button_box': 'button box for accept/cancel buttons',
                 'list': 'multiple choice list for choosing objects to apply to'}

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
        self.setWindowTitle('Choose tables to apply Pipeline to')
        self.button_box.accepted.connect(self.accept)
        self.button_box.rejected.connect(self.reject)
        self.layout.addWidget(self.label)
        self.layout.addWidget(self.list)
        self.layout.addWidget(self.button_box, QtCore.Qt.AlignCenter)

    def result(self):
        return [item.text() for item in self.list.get_sorted_selection()]


class ReportGenerationMessageBox(QtWidgets.QMessageBox):
    def __init__(self, parent=None):
        super().__init__(parent=parent)
        self.setWindowTitle("Report Generation")
        self.setText("Do you want to enable report generation for this session?\n"
                     "(this will slow down the program slightly)")

        self.yes_button = self.addButton(QtWidgets.QPushButton("Yes"), QtWidgets.QMessageBox.YesRole)
        self.no_button = self.addButton(QtWidgets.QPushButton("No"), QtWidgets.QMessageBox.NoRole)

        self.checkbox = QtWidgets.QCheckBox("Don't ask me again")
        self.setCheckBox(self.checkbox)

    def exec(self) -> Tuple[bool, bool]:
        super().exec()
        result = self.clickedButton()
        choice = self.buttonRole(result) == QtWidgets.QMessageBox.ButtonRole.YesRole
        checkbox_checked = self.checkbox.isChecked()
        return choice, checkbox_checked
