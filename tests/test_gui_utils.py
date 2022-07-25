import pytest
from pytestqt import qtbot
from PyQt5 import QtCore, QtWidgets, QtGui
import sys
from rnalysis.gui.gui_utils import *

LEFT_CLICK = QtCore.Qt.LeftButton
RIGHT_CLICK = QtCore.Qt.RightButton


def widget_setup(qtbot, widget_class, *args, **kwargs):
    widget = widget_class(*args, **kwargs)
    widget.show()
    qtbot.add_widget(widget)
    return qtbot, widget


def test_AboutWindow(qtbot, monkeypatch):
    exit_calls = []

    def mock_close(*args, **kwargs):
        exit_calls.append(1)

    monkeypatch.setattr(AboutWindow, 'close', mock_close)

    qtbot, window = widget_setup(qtbot, AboutWindow)
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
    qtbot, widget = widget_setup(qtbot, StrIntLineEdit)

    qtbot.keyClicks(widget, item)
    assert widget.text() == expected


@pytest.mark.parametrize("item,expected", [
    ('hello', 'hello'),
    ('57', '57'),
    ('', '')
])
def test_OptionalLineEdit(qtbot, item, expected):
    qtbot, widget = widget_setup(qtbot, OptionalLineEdit)
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
    qtbot, widget = widget_setup(qtbot, OptionalLineEdit)
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
    qtbot, widget = widget_setup(qtbot, OptionalSpinBox)
    widget.setMinimum(-100)
    widget.setMaximum(100)
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
    qtbot, widget = widget_setup(qtbot, OptionalSpinBox)
    widget.setMinimum(-100)
    widget.setMaximum(100)

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
    qtbot, widget = widget_setup(qtbot, OptionalDoubleSpinBox)
    widget.setMinimum(-100)
    widget.setMaximum(100)
    widget.setSingleStep(0.01)
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
    qtbot, widget = widget_setup(qtbot, OptionalDoubleSpinBox)
    widget.setMinimum(-100)
    widget.setMaximum(100)

    widget.setValue(item)
    assert widget.value() == expected

    if item is None:
        assert not widget.spinbox.isEnabled()
    else:
        assert widget.spinbox.isEnabled()


@pytest.mark.parametrize("item,expected", [
    ('black', '#000000'),
    ('#123456', '#123456'),
    ('r', '#FF0000')
])
def test_ColorPicker_written_colors(qtbot, item, expected):
    qtbot, widget = widget_setup(qtbot, ColorPicker)
    widget.color_line.clear()
    qtbot.keyClicks(widget.color_line, item)
    assert widget.text().lower() == expected.lower()


def test_ColorPicker_validColor(qtbot, monkeypatch):
    color = "#ccab56"

    def mock_get_color():
        return QtGui.QColor(color)

    monkeypatch.setattr(QtWidgets.QColorDialog, 'getColor', mock_get_color)
    qtbot, widget = widget_setup(qtbot, ColorPicker)
    widget.color_line.clear()
    qtbot.keyClicks(widget.color_line, 'black')

    qtbot.mouseClick(widget.pick_button, LEFT_CLICK)

    assert widget.text().lower() == color.lower()


def test_ColorPicker_set_default(qtbot, monkeypatch):
    prev_color = '#12345f'

    class MockColor(QtGui.QColor):
        def isValid(self):
            return False

    def mock_get_color():
        return MockColor('#fffffff')

    monkeypatch.setattr(QtWidgets.QColorDialog, 'getColor', mock_get_color)
    qtbot, widget = widget_setup(qtbot, ColorPicker)
    widget.color_line.clear()
    qtbot.keyClicks(widget.color_line, prev_color)

    qtbot.mouseClick(widget.pick_button, LEFT_CLICK)

    assert widget.text().lower() == prev_color.lower()


def test_MandatoryComboBox_is_legal(qtbot):
    qtbot, widget = widget_setup(qtbot, MandatoryComboBox, 'default_choice')
    assert not widget.is_legal()

    widget.addItems(['item1', 'item2'])
    qtbot.keyClicks(widget, 'item1')

    assert widget.currentText() == 'item1'
    assert widget.is_legal()


def test_MandatoryComboBox_clear(qtbot):
    default = 'default_chocie'
    qtbot, widget = widget_setup(qtbot, MandatoryComboBox, default)
    widget.addItems(['a', 'b', 'c'])
    widget.clear()
    assert widget.count() == 1
    assert widget.itemText(0) == default



