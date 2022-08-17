import pytest

from rnalysis.gui.gui_tutorial import *

LEFT_CLICK = QtCore.Qt.LeftButton
RIGHT_CLICK = QtCore.Qt.RightButton


def widget_setup(qtbot, widget_class, *args, **kwargs):
    widget = widget_class(*args, **kwargs)
    widget.show()
    qtbot.add_widget(widget)
    return qtbot, widget


def test_TutorialMovie_init(qtbot):
    qtbot, window = widget_setup(qtbot, TutorialMovie, 'tests/test_files/test_video.webp')
    assert window.video.isValid()


def test_TutorialMovie_set_frame(qtbot):
    qtbot, window = widget_setup(qtbot, TutorialMovie, 'tests/test_files/test_video.webp')
    window.video.start()
    window.pause()
    window.set_frame(17)
    assert window.video.state() == QtGui.QMovie.Paused
    assert window.video.currentFrameNumber() == 17
    window.set_frame(19)
    assert window.video.state() == QtGui.QMovie.Paused
    assert window.video.currentFrameNumber() == 19
    window.set_frame(19)
    assert window.video.state() == QtGui.QMovie.Paused
    assert window.video.currentFrameNumber() == 19
    window.set_frame(1)
    assert window.video.state() == QtGui.QMovie.Paused
    assert window.video.currentFrameNumber() == 1


def test_TutorialMovie_start(qtbot):
    qtbot, window = widget_setup(qtbot, TutorialMovie, 'tests/test_files/test_video.webp')
    assert window.video.state() == QtGui.QMovie.NotRunning
    window.start()
    assert window.video.state() == QtGui.QMovie.Running


def test_TutorialMovie_restart(qtbot):
    qtbot, window = widget_setup(qtbot, TutorialMovie, 'tests/test_files/test_video.webp')
    window.restart()
    assert window.video.currentFrameNumber() == 0
    assert window.video.state() == QtGui.QMovie.Running


def test_TutorialMovie_stop(qtbot):
    qtbot, window = widget_setup(qtbot, TutorialMovie, 'tests/test_files/test_video.webp')
    window.video.start()
    window.stop()
    assert window.video.state() == QtGui.QMovie.NotRunning


def test_TutorialMovie_pause(qtbot):
    qtbot, window = widget_setup(qtbot, TutorialMovie, 'tests/test_files/test_video.webp')
    window.video.start()
    window.pause()
    assert window.video.state() == QtGui.QMovie.Paused


def test_TutorialMovie_stop_button(qtbot):
    qtbot, window = widget_setup(qtbot, TutorialMovie, 'tests/test_files/test_video.webp')
    qtbot.mouseClick(window.stop_button, LEFT_CLICK)
    assert window.video.state() == QtGui.QMovie.NotRunning


def test_TutorialMovie_pause_button(qtbot):
    qtbot, window = widget_setup(qtbot, TutorialMovie, 'tests/test_files/test_video.webp')
    window.video.start()
    qtbot.mouseClick(window.play_button, LEFT_CLICK)
    assert window.video.state() == QtGui.QMovie.Paused
    qtbot.mouseClick(window.play_button, LEFT_CLICK)
    assert window.video.state() == QtGui.QMovie.Running


def test_TutorialMovie_pause_click_video(qtbot):
    qtbot, window = widget_setup(qtbot, TutorialMovie, 'tests/test_files/test_video.webp')
    window.video.start()
    qtbot.mouseClick(window, LEFT_CLICK)
    assert window.video.state() == QtGui.QMovie.Paused
    qtbot.mouseClick(window, LEFT_CLICK)
    assert window.video.state() == QtGui.QMovie.Running


def test_TutorialMovie_speed_button(qtbot):
    qtbot, window = widget_setup(qtbot, TutorialMovie, 'tests/test_files/test_video.webp')
    window.video.start()
    default_speed = window.video.speed()
    qtbot.mouseClick(window.speed_button, LEFT_CLICK)
    assert window.video.speed() > default_speed
    qtbot.mouseClick(window.speed_button, LEFT_CLICK)
    assert window.video.speed() == default_speed


def test_WelcomeWizard_init(qtbot):
    qtbot, window = widget_setup(qtbot, WelcomeWizard)
    for i in range(len(window.TITLES) + 1):
        qtbot.mouseClick(window.button(QtWidgets.QWizard.NextButton), LEFT_CLICK)
    qtbot.mouseClick(window.button(QtWidgets.QWizard.FinishButton), LEFT_CLICK)


@pytest.mark.parametrize('n_next', [0, 1, len(WelcomeWizard.TITLES)])
def test_WelcomeWizard_cancel(qtbot, n_next):
    qtbot, window = widget_setup(qtbot, WelcomeWizard)
    for i in range(n_next):
        qtbot.mouseClick(window.button(QtWidgets.QWizard.NextButton), LEFT_CLICK)
    qtbot.mouseClick(window.button(QtWidgets.QWizard.CancelButton), LEFT_CLICK)


def test_WelcomeWizard_do_not_show_again_finish(qtbot, monkeypatch):
    saved = []

    def save_func(show_tutorial):
        assert not show_tutorial
        saved.append(True)

    monkeypatch.setattr(settings, 'set_show_tutorial_settings', save_func)
    monkeypatch.setattr(settings, 'get_show_tutorial_settings', lambda: True)
    qtbot, window = widget_setup(qtbot, WelcomeWizard)
    assert not window.currentPage().dont_show_again.isChecked()
    window.currentPage().dont_show_again.setChecked(True)
    qtbot.mouseClick(window.button(QtWidgets.QWizard.FinishButton), LEFT_CLICK)

    assert saved[0]


@pytest.mark.parametrize("show", [True, False])
def test_WelcomeWizard_do_not_show_again_cancel(qtbot, monkeypatch, show):
    saved = []

    def save_func(show_tutorial):
        assert show_tutorial != show
        saved.append(True)

    monkeypatch.setattr(settings, 'set_show_tutorial_settings', save_func)
    monkeypatch.setattr(settings, 'get_show_tutorial_settings', lambda: show)
    qtbot, window = widget_setup(qtbot, WelcomeWizard)
    assert window.currentPage().dont_show_again.isChecked() != show
    window.currentPage().dont_show_again.setChecked(show)
    qtbot.mouseClick(window.button(QtWidgets.QWizard.CancelButton), LEFT_CLICK)
    assert saved[0]
