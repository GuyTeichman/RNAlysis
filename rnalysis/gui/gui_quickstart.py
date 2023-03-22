from pathlib import Path

from PyQt5 import QtCore, QtWidgets, QtGui

from rnalysis.utils import settings, io


class QuickStartWizard(QtWidgets.QWizard):
    TITLES = (
        "Load a table",
        "Examine your table",
        "Filter your table",
        "Undo the operations you applied to your table",
        "Apply your operations 'in-place' or apply on a new table",
        "Save the changes you made to your table",
        "Work on multiple tables at the same time",
        "Different types of tables offer different ways to filter and analyze your data",
        "Create and save graphs",
        "Quickly look-up genes in your database of choice",
        "Sort your tabs and change their icons",
        "Restore tabs you closed",
        "Import lists of genes as Gene Sets",
        "Visualize the intersections between your tables and gene sets",
        "Apply set operations to your tables and gene sets",
        "Perform enrichment analysis on your tables and gene sets",
        "Create Pipelines to streamline your data analysis",
        "Apply Pipelines to one or more of your tables",
        "Export and share Pipelines to make your analysis more reproducible",
        "Interface with other bioinformatic tools")

    CONTENTS = (
        "Choose a file from your computer, and click 'start' to load it into <i>RNAlysis</i>. ",

        "You will now be able to see an overview of your data, including the table's name, type, and dimensions. "
        "Click the 'View full table' button to see your table in its entirety. ",

        "Choose a filtering function, set your desired parameters, and click 'Apply' to filter your table. "
        "The changes you make will not affect your original file until you save them. ",

        "At any moment, you can use the 'Command history' window "
        "to undo or redo an operation you applied to your table. ",

        "Instead of applying operations 'in-place', "
        "you can choose to apply the operation to a copy of your table in a new tab. "
        "The table in the original tab won't be modified. ",

        "To save result of your filtering operations, "
        "click the 'Save table' button and choose where to save the modified table. ",

        "You can work on multiple tables at the same time by opening a new tab and loading another table. ",

        "When loading a table, you can specify its type. "
        "Different types of tables support different types of functions: for example, "
        "count matrices support clustering analysis. ",

        "Some functions can generate graphs of your data. "
        "You can resize those graphs, and save them to your computer in multiple file formats. ",

        "Easily get information about your genes with a right-click. "
        "Select from a range of databases such as NCBI Genes, UniProtKB and others, "
        "which can be configured from the settings menu. ",

        "To help organize your workspace, you can sort tabs by "
        "right-clicking a tab and choosing a sorting method. "
        "You can also change specific tabs' colors, to help you differentiate them. ",

        "If you accidentally closed one of your tabs - don't worry! "
        "You can restore closed tabs through the 'Edit' menu. ",

        "In addition to tables, <i>RNAlysis</i> can also import lists of genes as Gene Sets. "
        "We will soon review what we can do with those gene sets. ",

        "In the 'Visualize Gene Sets' window you can create Venn diagrams and UpSet plots "
        "that will display the various intersections between your tables and gene sets. ",

        "In the 'Set Operations' window you can extract specific subsets from your data. "
        "Either use predefined set operations, or click on specific subsets in the preview pane to select them. ",

        "In the 'Enrichment Analysis' window, you can perform various types of enrichment analysis "
        "on the tables and gene sets you filtered. ",

        "You can group multiple operations in a specific order and with specific parameters into a Pipeline. "
        "Just add those functions to the Pipeline in the order you choose. ",

        "You can apply a Pipeline to a group of tables through the 'Pipelines' menu. "
        "Using Pipelines to analyze multiple datasets can make your workflow faster and less error-prone. ",

        "Pipelines you export can be imported from any computer, and can be shared with others "
        "to help make your analysis easier to understand and more reproducible.",

        "<i>RNAlysis</i> offers a graphic interface for many bioinformatic tools. "
        "Analyze your data at any stage - adapter trimming, alignment, feature counting, or differential expression. ",
    )

    VIDEO_FILES = (
        'load_table.webp',
        'view_table.webp',
        'filter_table.webp',
        'undo_actions.webp',
        'apply_inplace.webp',
        'save_table.webp',
        'new_tab.webp',
        'table_types.webp',
        'generate_graphs.webp',
        'quick_search.webp',
        'sort_tabs.webp',
        'restore_tabs.webp',
        'import_gene_sets.webp',
        'visualize_gene_sets.webp',
        'set_operations.webp',
        'enrichment_analysis.webp',
        'create_pipeline.webp',
        'apply_pipeline.webp',
        'export_pipeline.webp',
        'external_windows.webp',
    )

    def __init__(self, parent=None):
        super().__init__(parent)
        self.addPage(StartPage(self))

        for title, content, filename in zip(self.TITLES, self.CONTENTS, self.VIDEO_FILES):
            page = QuickStartPage(title, content, filename, self)
            self.addPage(page)

        self.addPage(EndPage(self))

        self.prev_id = self.currentId()

        self.currentIdChanged.connect(self.play_tutorial)

        self.setWizardStyle(QtWidgets.QWizard.ModernStyle)
        self.setWindowTitle('Welcome to RNAlysis!')
        self.setPixmap(self.LogoPixmap,
                       QtGui.QPixmap(Path.cwd().parent.parent.joinpath('docs/source/favicon.ico').as_posix()))
        self.setField('dont_show_again', not settings.get_show_tutorial_settings())

    def play_tutorial(self, ind: int):
        page = self.page(ind)
        if isinstance(page, QuickStartPage):
            self.currentPage().start_movie()
        try:
            prev_page = self.page(self.prev_id)
            if isinstance(prev_page, QuickStartPage):
                prev_page.stop_movie()
        except IndexError:
            pass

        self.prev_id = ind

    def show(self):
        QtWidgets.QApplication.processEvents()
        super().show()
        self.next()
        self.back()

        rect = self.geometry().getRect()
        self.setGeometry(rect[0], 100, *rect[2:])

    def accept(self):
        self.save_settings()
        super().accept()

    def reject(self):
        self.save_settings()
        super().reject()

    def save_settings(self):
        show_tutorial = not self.field('dont_show_again')
        settings.set_show_tutorial_settings(show_tutorial)

    def closeEvent(self, *args, **kwargs):
        self.save_settings()
        super().closeEvent(*args, **kwargs)


class EndPage(QtWidgets.QWizardPage):
    USER_GUIDE_URL = "https://guyteichman.github.io/RNAlysis/build/user_guide.html"
    DOCUMENTATION_URL = "://guyteichman.github.io/RNAlysis"
    TEXT = f'If you want to learn more about the various features <i>RNAlysis</i> has to offer, ' \
           f'you can read more about them in the <a href="{USER_GUIDE_URL}"><b>user guide </b></a>' \
           f'or the complete <a href="{DOCUMENTATION_URL}"><b>documentation</b></a>. Good luck!'

    def __init__(self, parent=None):
        super().__init__(parent)
        self.setTitle("<b>You are now ready to use <i>RNAlysis</i>!</b>")
        self.layout = QtWidgets.QVBoxLayout(self)
        self.label = QtWidgets.QLabel(self.TEXT)
        self.label.setWordWrap(True)
        self.layout.addWidget(self.label)
        self.setFinalPage(True)


class StartPage(QtWidgets.QWizardPage):
    def __init__(self, parent=None):
        super().__init__(parent)
        self.setTitle('<b><i>Welcome to RNAlysis</i></b>')
        self.setSubTitle('This quick-start guide will lead you through the basic usage of <i>RNAlysis</i>. ')
        self.layout = QtWidgets.QGridLayout(self)
        self.dont_show_again = QtWidgets.QCheckBox("Do not show this window again")
        self.dont_show_again.stateChanged.connect(self.setFinalPage)
        self.dont_show_again.stateChanged.connect(self.completeChanged.emit)
        self.registerField('dont_show_again', self.dont_show_again)
        img_path = str(Path(__file__).parent.joinpath('logo_small.png'))
        text = f'<img src="{img_path}" class="center" style="vertial-align:middle"/><br><br>'

        self.image = QtWidgets.QLabel(text)
        self.layout.addWidget(self.image)
        self.layout.addWidget(self.dont_show_again)


class QuickStartMovie(QtWidgets.QWidget):
    clicked = QtCore.pyqtSignal()

    def __init__(self, video_path: Path, parent=None):
        super().__init__(parent)
        full_video_path = str(video_path.absolute())
        self.layout = QtWidgets.QGridLayout(self)
        self.label = QtWidgets.QLabel(self)
        self.video = QtGui.QMovie(full_video_path)
        self.video.setCacheMode(self.video.CacheAll)
        self.label.setMovie(self.video)

        self.play_button = QtWidgets.QToolButton()
        self.play_button.clicked.connect(self.change_play_state)

        self.stop_button = QtWidgets.QToolButton()
        self.stop_button.setIcon(self.style().standardIcon(QtWidgets.QStyle.SP_MediaStop))
        self.stop_button.clicked.connect(self.stop)

        self.speed_button = QtWidgets.QToolButton()
        self.speed_button.setIcon(self.style().standardIcon(QtWidgets.QStyle.SP_MediaSeekForward))
        self.speed_button.setCheckable(True)
        self.speed_button.clicked.connect(self.change_speed)

        self.position_slider = QtWidgets.QSlider(QtCore.Qt.Horizontal)
        self.position_slider.setRange(0, self.video.frameCount())
        self.position_slider.valueChanged.connect(self.set_frame)
        self.position_slider.setTracking(False)

        self.video.frameChanged.connect(self.update_slider)

        self.press_pos = None
        self.paused = False
        self.clicked.connect(self.change_play_state)
        self.label.setStyleSheet("""
                QLabel {border-style: outset;
                border-width: 3px;
                border-color: black;}""")

        pixmap = QtGui.QPixmap(full_video_path)
        size = pixmap.size()

        size.scale(750, 750, QtCore.Qt.KeepAspectRatio)
        self.video.setScaledSize(size)

        self.layout.addWidget(self.label, 0, 0, 4, 8)
        self.layout.addWidget(self.play_button, 4, 0)
        self.layout.addWidget(self.stop_button, 4, 1)
        self.layout.addWidget(self.speed_button, 4, 2)
        self.layout.addWidget(self.position_slider, 4, 3, 1, 5)

        self.update_play_button()

    def update_slider(self, frame: int):
        if self.position_slider.isSliderDown():
            return
        self.position_slider.setSliderPosition(frame)

    def change_play_state(self):
        self.paused = not self.paused
        self.video.setPaused(self.paused)
        self.update_play_button()

    def update_play_button(self):
        if self.paused:
            self.play_button.setIcon(self.style().standardIcon(QtWidgets.QStyle.SP_MediaPlay))
        else:
            self.play_button.setIcon(self.style().standardIcon(QtWidgets.QStyle.SP_MediaPause))
        QtWidgets.QApplication.processEvents()

    def mousePressEvent(self, event):
        if event.button() == QtCore.Qt.LeftButton:
            self.press_pos = event.pos()

    def mouseReleaseEvent(self, event):
        # ensure that the left button was pressed *and* released within the
        # geometry of the widget; if so, emit the signal;
        if self.press_pos is not None and event.button() == QtCore.Qt.LeftButton and event.pos() in self.rect():
            self.clicked.emit()
        self.press_pos = None

    def pause(self):
        self.video.setPaused(True)
        self.paused = True
        self.update_play_button()

    def stop(self):
        self.video.stop()
        self.video.start()
        self.video.stop()
        self.paused = True
        self.update_play_button()

    def start(self):
        self.video.start()
        self.paused = False
        self.update_play_button()

    def restart(self):
        self.video.stop()
        self.video.start()
        self.paused = False
        self.update_play_button()
        self.speed_button.setChecked(False)

    def set_frame(self, desired_frame: int):
        if desired_frame == self.video.currentFrameNumber():
            return
        else:
            self.video.jumpToFrame(desired_frame)

    def change_speed(self):
        if self.speed_button.isChecked():
            self.video.setSpeed(400)
        else:
            self.video.setSpeed(100)


class QuickStartPage(QtWidgets.QWizardPage):
    def __init__(self, title: str, content: str, video_name: str, parent=None):
        super().__init__(parent)
        self.setTitle('<b>' + title + '</b>')
        self.setSubTitle(content)
        self.layout = QtWidgets.QVBoxLayout(self)
        self.video_path = Path.joinpath(Path(__file__).parent, io.get_tutorial_videos_dir().joinpath(video_name))
        self.label = QtWidgets.QLabel()
        self.movie = QuickStartMovie(self.video_path)
        self.layout.addWidget(self.movie)

    def start_movie(self):
        self.movie.restart()

    def stop_movie(self):
        self.movie.pause()
