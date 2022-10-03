import itertools
import traceback
from pathlib import Path
import warnings
from queue import Queue

import pandas as pd
from PyQt5 import QtCore, QtWidgets, QtGui

from rnalysis import __version__
from rnalysis.gui import gui_style, gui_widgets
from rnalysis.utils import settings, io


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
            io.save_gene_set(self.data, filename)
            print(f"Successfully saved at {io.get_datetime()} under {filename}")


class DataFrameView(DataView):
    def __init__(self, data: pd.DataFrame, name: str, parent=None):
        super().__init__(data, name, parent)
        shape = self.data.shape
        if len(shape) == 1:
            shape = (shape[0], 1)
        self.label = QtWidgets.QLabel(f"Table '{name}': {shape[0]} rows, {shape[1]} columns")

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


class SettingsWindow(gui_widgets.MinMaxDialog):
    THEMES = gui_style.get_stylesheet_names()
    FONT_SIZES = [str(i) for i in [6, 7, 8, 9, 10, 11, 12, 14, 16, 18, 20, 22, 24, 28, 32, 36, 48, 72]]
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
        current_font, current_font_size, current_theme, current_show_tutorial = settings.get_gui_settings()
        current_theme = {val: key for key, val in self.THEMES.items()}[current_theme]

        self.appearance_widgets['app_theme'].setCurrentText(current_theme)
        self.appearance_widgets['app_font'].setCurrentText(current_font)
        self.appearance_widgets['app_font_size'].setCurrentText(str(current_font_size))
        self.appearance_widgets['show_tutorial'].setChecked(current_show_tutorial)

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
        # self.appearance_widgets['app_font'].addItems(list(QtGui.QFontDatabase().families()))
        self.appearance_widgets['app_font'].setEditable(True)
        self.appearance_widgets['app_font'].completer().setCompletionMode(QtWidgets.QCompleter.PopupCompletion)
        self.appearance_widgets['app_font'].setInsertPolicy(QtWidgets.QComboBox.NoInsert)

        self.appearance_widgets['app_font_size'] = QtWidgets.QComboBox(self.appearance_group)
        self.appearance_widgets['app_font_size'].addItems(self.FONT_SIZES)

        self.appearance_widgets['show_tutorial'] = QtWidgets.QCheckBox("Show tutorial page on startup")

        for widget_name in ['app_theme', 'app_font', 'app_font_size']:
            self.appearance_widgets[widget_name].currentIndexChanged.connect(self._trigger_settings_changed)
        self.appearance_widgets['show_tutorial'].stateChanged.connect(self._trigger_settings_changed)

        self.appearance_grid.addWidget(QtWidgets.QLabel('Theme'), 0, 0)
        self.appearance_grid.addWidget(self.appearance_widgets['app_theme'], 0, 1)
        self.appearance_grid.addWidget(QtWidgets.QLabel('Application Font'), 1, 0)
        self.appearance_grid.addWidget(self.appearance_widgets['app_font'], 1, 1)
        self.appearance_grid.addWidget(QtWidgets.QLabel('Application Font Size'), 2, 0)
        self.appearance_grid.addWidget(self.appearance_widgets['app_font_size'], 2, 1)
        self.appearance_grid.addWidget(self.appearance_widgets['show_tutorial'], 3, 0, 1, 2)

    def save_settings(self):
        if self.settings_changed:
            self.settings_changed = False
            font = self.appearance_widgets['app_font'].currentText()
            font_size = int(self.appearance_widgets['app_font_size'].currentText())
            theme = self.THEMES[self.appearance_widgets['app_theme'].currentText()]
            show_tutorial = self.appearance_widgets['show_tutorial'].isChecked()
            settings.set_gui_settings(font, font_size, theme, show_tutorial)

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
    CITATION = f"Teichman, G. (2021) RNAlysis: RNA Sequencing analysis software (Python package version {__version__})."

    def __init__(self, parent=None):
        super().__init__(parent)
        text = f"""<p><b><i>RNAlysis</i> version {__version__}</b>
                </p>
                <br>
                <img src="../../docs/source/logo.png" width="500"/>"""
        self.label = QtWidgets.QLabel(text)
        self.text_edit = QtWidgets.QTextEdit(self.CITATION)
        self.ok_button = QtWidgets.QPushButton('OK')
        self.copy_button = QtWidgets.QPushButton('Copy to clipboard')
        self.copied_label = QtWidgets.QLabel()

        self.layout = QtWidgets.QVBoxLayout(self)
        self.init_ui()

    def init_ui(self):
        self.setWindowTitle("How to cite RNAlysis")

        self.label.setWordWrap(True)
        self.layout.addWidget(self.label)
        self.text_edit.setReadOnly(True)
        self.layout.addWidget(self.text_edit)

        self.ok_button.clicked.connect(self.close)
        self.layout.addWidget(self.ok_button)

        self.copy_button.clicked.connect(self.copy_to_clipboard)
        self.layout.addWidget(self.copy_button)
        self.layout.addWidget(self.copied_label)

    def copy_to_clipboard(self):
        cb = QtWidgets.QApplication.clipboard()
        cb.clear(mode=cb.Clipboard)
        cb.setText(self.text_edit.toPlainText())
        self.copied_label.setText('Copied to clipboard')


def splash_screen():
    img_path = str(Path(__file__).parent.joinpath('splash.png'))
    splash_pixmap = QtGui.QPixmap(img_path)
    splash_font = QtGui.QFont('Calibri', 16)
    splash = QtWidgets.QSplashScreen(splash_pixmap)
    splash.setFont(splash_font)
    splash.showMessage(f"<i>RNAlysis</i> version {__version__}", QtCore.Qt.AlignBottom | QtCore.Qt.AlignHCenter)
    splash.show()
    return splash
