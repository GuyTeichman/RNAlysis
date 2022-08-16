from rnalysis.gui.gui_tutorial import *

LEFT_CLICK = QtCore.Qt.LeftButton
RIGHT_CLICK = QtCore.Qt.RightButton


def widget_setup(qtbot, widget_class, *args, **kwargs):
    widget = widget_class(*args, **kwargs)
    widget.show()
    qtbot.add_widget(widget)
    return qtbot, widget


def test_TutorialMovie_init(qtbot):
    assert False


def test_TutorialMovie_buttons(qtbot):
    assert False


def test_WelcomeWizard_init(qtbot):
    qtbot, window = widget_setup(qtbot, WelcomeWizard)
    for i in range(len(window.TITLES) + 1):
        qtbot.mouseClick(window.button(QtWidgets.QWizard.NextButton), LEFT_CLICK)
    qtbot.mouseClick(window.button(QtWidgets.QWizard.FinishButton), LEFT_CLICK)


def test_WelcomeWizard_do_not_show_again_saved(qtbot):
    assert False
