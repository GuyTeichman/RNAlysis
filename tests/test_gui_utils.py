import pytest
from pytestqt import qtbot
from PyQt5 import QtCore, QtWidgets
import sys
from rnalysis.gui.gui_utils import *

LEFT_CLICK = QtCore.Qt.LeftButton
RIGHT_CLICK = QtCore.Qt.RightButton


def test_AboutWindow(qtbot, monkeypatch):
    exit_calls = []

    def mock_close(*args, **kwargs):
        exit_calls.append(1)

    monkeypatch.setattr(AboutWindow, 'close', mock_close)

    window = AboutWindow()
    window.show()
    qtbot.add_widget(window)
    qtbot.mouseClick(window.buttons()[0], LEFT_CLICK)
    assert exit_calls == [1]


@pytest.mark.parametrize("item,expected", [
    ('hello', 'hello'),
    ('', ''),
    ('55', 55),
    ('high5', 'high5'),
    ('-37', -37)
])
def test_StrIntLineEdit(qtbot, item, expected):
    widget = StrIntLineEdit()
    widget.show()
    qtbot.add_widget(widget)
    qtbot.keyClicks(widget, item)
    assert widget.text() == expected


@pytest.mark.parametrize("item,expected", [
    ('hello', 'hello'),
    ('57', '57'),
    ('', '')
])
def test_OptionalLineEdit(qtbot, item, expected):
    widget = OptionalLineEdit()
    widget.show()
    qtbot.add_widget(widget)
    qtbot.keyClicks(widget.line, item)

    assert widget.line.isEnabled()
    assert widget.text() == item

    qtbot.mouseClick(widget.checkbox, LEFT_CLICK)

    assert widget.text() is None
    assert not widget.line.isEnabled()

    qtbot.mouseClick(widget.checkbox, LEFT_CLICK)

    assert widget.line.isEnabled()
    assert widget.text() == item


@pytest.mark.parametrize("item,expected", [
    ('hello', 'hello'),
    (None, None),
    ('57', '57'),
    ('', '')
])
def test_OptionalLineEdit_setText(qtbot, item, expected):
    widget = OptionalLineEdit()
    widget.show()
    qtbot.add_widget(widget)

    widget.setText(item)
    assert widget.text() == expected

    if item is None:
        assert not widget.line.isEnabled()
    else:
        assert widget.line.isEnabled()


@pytest.mark.parametrize("item,expected", [
    ('57', 57),
    ('-32', -32),
    ('0', 0),
])
def test_OptionalSpinBox(qtbot, item, expected):
    widget = OptionalSpinBox()
    widget.setMinimum(-100)
    widget.setMaximum(100)
    widget.show()
    qtbot.add_widget(widget)
    widget.spinbox.clear()
    qtbot.keyClicks(widget.spinbox, item)

    path = qtbot.screenshot(widget)
    assert widget.spinbox.isEnabled()
    assert widget.value() == expected

    qtbot.mouseClick(widget.checkbox, LEFT_CLICK)

    assert widget.value() is None
    assert not widget.spinbox.isEnabled()

    qtbot.mouseClick(widget.checkbox, LEFT_CLICK)

    assert widget.spinbox.isEnabled()
    assert widget.value() == expected


@pytest.mark.parametrize("item,expected", [
    (57, 57),
    (-32, -32),
    (0, 0),
    (None, None)
])
def test_OptionalSpinBox_setValue(qtbot, item, expected):
    widget = OptionalSpinBox()
    widget.setMinimum(-100)
    widget.setMaximum(100)
    widget.show()
    qtbot.add_widget(widget)

    widget.setValue(item)
    assert widget.value() == expected

    if item is None:
        assert not widget.spinbox.isEnabled()
    else:
        assert widget.spinbox.isEnabled()


@pytest.mark.parametrize("item,expected", [
    ('57', 57),
    ('-32', -32),
    ('0', 0),
    ('3.14', 3.14),
    ('-0.75', -0.75)
])
def test_OptionalDoubleSpinBox(qtbot, item, expected):
    widget = OptionalDoubleSpinBox()
    widget.setMinimum(-100)
    widget.setMaximum(100)
    widget.setSingleStep(0.01)
    widget.show()
    qtbot.add_widget(widget)
    widget.spinbox.clear()
    qtbot.keyClicks(widget.spinbox, item)

    path = qtbot.screenshot(widget)
    assert widget.spinbox.isEnabled()
    assert widget.value() == expected

    qtbot.mouseClick(widget.checkbox, LEFT_CLICK)

    assert widget.value() is None
    assert not widget.spinbox.isEnabled()

    qtbot.mouseClick(widget.checkbox, LEFT_CLICK)

    assert widget.spinbox.isEnabled()
    assert widget.value() == expected


@pytest.mark.parametrize("item,expected", [
    (57, 57),
    (-32, -32),
    (0, 0),
    (3.14, 3.14),
    (-0.75, -0.75),
    (None, None)
])
def test_OptionalDoubleSpinBox_setValue(qtbot, item, expected):
    widget = OptionalDoubleSpinBox()
    widget.setMinimum(-100)
    widget.setMaximum(100)
    widget.show()
    qtbot.add_widget(widget)

    widget.setValue(item)
    assert widget.value() == expected

    if item is None:
        assert not widget.spinbox.isEnabled()
    else:
        assert widget.spinbox.isEnabled()
