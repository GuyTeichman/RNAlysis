import pytest
from pytestqt import qtbot
from PyQt5 import QtCore, QtWidgets
import sys
from rnalysis.gui.gui_utils import *


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
    print(type(widget), isinstance(widget, QtWidgets.QWidget))
    qtbot.add_widget(widget)
    qtbot.keyClicks(widget, item)
    assert widget.text() == expected
