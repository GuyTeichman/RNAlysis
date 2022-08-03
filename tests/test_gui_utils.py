import os

import pytest
from rnalysis.gui.gui_utils import *
from utils.io import save_gene_set

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


def test_clear_layout(qtbot):
    qtbot, widget = widget_setup(qtbot, QtWidgets.QWidget)
    layout = QtWidgets.QGridLayout(widget)
    layout.addWidget(QtWidgets.QSpinBox(), 0, 0)
    layout.addWidget(QtWidgets.QLineEdit(), 1, 2)
    layout.addWidget(QtWidgets.QLabel("test"), 3, 3)

    clear_layout(layout)

    assert layout.count() == 0


@pytest.mark.parametrize("gene_set,expected_split", [
    ({1, 2, 3}, ['1', '2', '3']),
    ({'geneA', 'geneB', 'geneC', 'geneD'}, ["geneA", "geneB", "geneC", "geneD"])
])
def test_save_gene_set(gene_set, expected_split):
    pth = 'tests/test_files/tmp_saved_gene_set.txt'
    try:
        save_gene_set(gene_set, pth)
        with open(pth) as f:
            split = f.read().split('\n')
        assert sorted(split) == sorted(expected_split)
    finally:
        try:
            os.unlink(pth)
        except FileNotFoundError:
            pass


class FilledComboBox(QtWidgets.QComboBox):
    def __init__(self, parent=None):
        super().__init__(parent)
        self.addItems(['test1', 'test2', '12'])

    def clear(self):
        pass


@pytest.mark.parametrize("widget_class,keyboard_interact,expected_val", [
    (QtWidgets.QCheckBox, False, True),
    (QtWidgets.QLineEdit, True, "12"),
    (QtWidgets.QSpinBox, True, 12),
    (QtWidgets.QTextEdit, True, 12),
    (FilledComboBox, True, '12')

])
def test_get_val_from_widget_native_types(qtbot, widget_class, keyboard_interact, expected_val):
    qtbot, widget = widget_setup(qtbot, widget_class)
    if keyboard_interact:
        widget.clear()
        qtbot.keyClicks(widget, "12")
    else:
        qtbot.mouseClick(widget, LEFT_CLICK)
    assert get_val_from_widget(widget) == expected_val


@pytest.mark.parametrize("widget_class,keyboard_interact,attr,expected_val,kwargs", [
    (OptionalSpinBox, True, 'spinbox', 12, {}),
    (OptionalDoubleSpinBox, True, 'spinbox', 12, {}),
    (PathLineEdit, True, 'file_path', '12', {}),
    (StrIntLineEdit, True, None, 12, {}),
    (OptionalLineEdit, True, 'line', '12', {}),
    (ComboBoxOrOtherWidget, True, 'combo', '12',
     {'items': ['opt1', 'opt2', '12'], 'default': 'opt1'}),
    (ComboBoxOrOtherWidget, True, 'other', '12',
     {'items': ['opt1', 'opt2', 'opt3'], 'default': None}),
])
def test_get_val_from_widget_nonnative_types(qtbot, widget_class, keyboard_interact, attr, expected_val, kwargs):
    if widget_class == ComboBoxOrOtherWidget:
        kwargs['other'] = QtWidgets.QLineEdit()
        kwargs['other'].setText('12')
    qtbot, widget = widget_setup(qtbot, widget_class, **kwargs)
    interact_with = widget if attr is None else getattr(widget, attr)
    if keyboard_interact:
        widget.clear()
        qtbot.keyClicks(interact_with, "12")
    else:
        qtbot.mouseClick(interact_with, LEFT_CLICK)
    assert get_val_from_widget(widget) == expected_val


@pytest.mark.parametrize("widget_class,default,excepted_val_empty,expected_val,kwargs", [
    (QMultiSpinBox, [0, 2, 3], 0, [0, 2, 3], {}),
    (QMultiDoubleSpinBox, [0.1, 3.2, 5], 0.0, [0.1, 3.2, 5], {}),
    (QMultiLineEdit, ['', 'text', 'other text'], '', ['', 'text', 'other text'], {}),
    (QMultiStrIntLineEdit, ['3', '-7', 'text', 'othertext', 'param5'], '', [3, -7, 'text', 'othertext', 'param5'], {}),
    (QMultiBoolComboBox, [True, True, False, True], True, [True, True, False, True], {}),
    (MultiColorPicker, ['r', 'black', '#0000ff'], None, ['#ff0000', '#000000', '#0000ff'], {}),
    (QMultiComboBox, ['option3', 'option2', 'option2'], 'option1', ['option3', 'option2', 'option2'],
     {'items': ['option1', 'option2', 'option3']})
])
def test_get_val_from_widget_multiinput_types(qtbot, widget_class, default, excepted_val_empty, expected_val, kwargs):
    qtbot, widget = widget_setup(qtbot, widget_class, label='label', **kwargs)
    assert get_val_from_widget(widget) == excepted_val_empty

    widget.set_defaults(default)
    assert get_val_from_widget(widget) == expected_val


@pytest.mark.parametrize("widget_class", (QtWidgets.QWidget, QtWidgets.QDial))
def test_get_val_from_widget_bad_widget(qtbot, widget_class):
    qtbot, widget = widget_setup(qtbot, widget_class)
    with pytest.raises(TypeError):
        get_val_from_widget(widget)


def test_ErrorMessage_message(qtbot):
    err_text = 'my error text'
    try:
        raise ValueError(err_text)
    except ValueError as e:
        err_tb = e.__traceback__
        err_value = e
    qtbot, dialog = widget_setup(qtbot, ErrorMessage, ValueError, err_value, err_tb)
    assert 'ValueError' in dialog.widgets['error_text'].toPlainText()
    assert err_text in dialog.widgets['error_text'].toPlainText()


def test_ErrorMessage_close(qtbot, monkeypatch):
    closed = []

    def mock_close(*args, **kwargs):
        closed.append(1)

    monkeypatch.setattr(ErrorMessage, 'close', mock_close)
    err_text = 'my error text'
    try:
        raise ValueError(err_text)
    except ValueError as e:
        err_tb = e.__traceback__
        err_value = e
    qtbot, dialog = widget_setup(qtbot, ErrorMessage, ValueError, err_value, err_tb)
    qtbot.mouseClick(dialog.widgets['ok_button'], LEFT_CLICK)

    assert closed == [1]


def test_ErrorMessage_copy_to_clipboard(qtbot, monkeypatch):
    err_text = 'my error text'
    try:
        raise ValueError(err_text)
    except ValueError as e:
        err_tb = e.__traceback__
        err_value = e
    qtbot, dialog = widget_setup(qtbot, ErrorMessage, ValueError, err_value, err_tb)
    qtbot.mouseClick(dialog.widgets['copy_button'], LEFT_CLICK)

    assert 'ValueError' in QtWidgets.QApplication.clipboard().text()
    assert err_text in QtWidgets.QApplication.clipboard().text()


def test_HelpButton(qtbot, monkeypatch):
    param_name = 'myparam'
    desc = 'mydesc'
    help_shown = []

    def mock_show_text(pos, text):
        assert param_name in text
        assert desc in text
        assert pos == QtGui.QCursor.pos()
        help_shown.append(1)

    monkeypatch.setattr(QtWidgets.QToolTip, 'showText', mock_show_text)

    qtbot, widget = widget_setup(qtbot, HelpButton)
    widget.connect_param_help(param_name, desc)

    qtbot.mouseClick(widget, LEFT_CLICK)

    assert help_shown == [1]


def test_PathLineEdit_is_legal(qtbot):
    qtbot, widget = widget_setup(qtbot, PathLineEdit)

    assert not widget.is_legal

    widget.clear()
    assert not widget.is_legal

    qtbot.keyClicks(widget.file_path, 'path/that/doesnt/exist.png')
    assert not widget.is_legal

    widget.clear()
    qtbot.keyClicks(widget.file_path, 'tests/test_files/test_deseq.csv')
    assert widget.is_legal


def test_PathLineEdit_text(qtbot):
    qtbot, widget = widget_setup(qtbot, PathLineEdit)
    widget.setText('test123')
    assert widget.text() == 'test123'

    widget.clear()
    qtbot.keyClicks(widget.file_path, 'test456')
    assert widget.text() == 'test456'


def test_PathLineEdit_choose_file(qtbot, monkeypatch):
    pth = 'path/to/a/file'

    def mock_choose_file(*args, **kwargs):
        return pth, None

    monkeypatch.setattr(QtWidgets.QFileDialog, 'getOpenFileName', mock_choose_file)
    qtbot, widget = widget_setup(qtbot, PathLineEdit)
    qtbot.mouseClick(widget.open_button, LEFT_CLICK)

    assert widget.text() == pth


def test_PathLineEdit_choose_file_not_chosen(qtbot, monkeypatch):
    pth = 'path/to/a/file'

    def mock_choose_file(*args, **kwargs):
        return False, None

    monkeypatch.setattr(QtWidgets.QFileDialog, 'getOpenFileName', mock_choose_file)
    qtbot, widget = widget_setup(qtbot, PathLineEdit)
    widget.setText(pth)
    qtbot.mouseClick(widget.open_button, LEFT_CLICK)

    assert widget.text() == pth


def test_MultipleChoiceList_select_all(qtbot):
    items = ['item1', 'item2', 'item3']
    qtbot, widget = widget_setup(qtbot, MultipleChoiceList, items)
    qtbot.mouseClick(widget.select_all_button, LEFT_CLICK)

    selected_items = [item.text() for item in widget.selectedItems()]
    assert selected_items == items


def test_MultipleChoiceList_clear_selection(qtbot):
    items = ['item1', 'item2', 'item3']
    qtbot, widget = widget_setup(qtbot, MultipleChoiceList, items)
    for item in widget.list_items:
        item.setSelected(True)

    assert len(widget.selectedItems()) == len(items)

    qtbot.mouseClick(widget.clear_all_button, LEFT_CLICK)

    selected_items = widget.selectedItems()
    assert len(selected_items) == 0


def test_MultipleChoiceList_emit(qtbot):
    selection_changed = []

    def selection_changed_slot(*args, **kwargs):
        selection_changed.append(1)

    items = ['item1', 'item2', 'item3']
    qtbot, widget = widget_setup(qtbot, MultipleChoiceList, items)
    widget.itemSelectionChanged.connect(selection_changed_slot)

    assert len(selection_changed) == 0

    qtbot.mouseClick(widget.select_all_button, LEFT_CLICK)
    assert selection_changed == [1] * len(items)
    qtbot.mouseClick(widget.clear_all_button, LEFT_CLICK)
    assert selection_changed == [1] * len(items) * 2
    widget.list_items[0].setSelected(True)
    assert selection_changed == [1] * len(items) * 2 + [1]


def test_MultipleChoiceList_icons(qtbot):
    pixmap = QtGui.QPixmap(32, 32)
    pixmap.fill(QtCore.Qt.transparent)
    icon = QtGui.QIcon(pixmap)

    pixmap2 = QtGui.QPixmap(34, 34)
    pixmap2.fill(QtCore.Qt.white)
    icon2 = QtGui.QIcon(pixmap2)
    items = ['item1', 'item2', 'item3']
    icons = [icon, icon2, icon2]
    qtbot, widget = widget_setup(qtbot, MultipleChoiceList, items, icons)

    actual_icons = [item.icon() for item in widget.list_items]

    assert [icon.pixmap(32, 32).toImage() for icon in icons] == [icon.pixmap(32, 32).toImage() for icon in actual_icons]


def test_RadioButtonBox_set_selection(qtbot):
    actions = ['action1', 'action2', 'action3']
    qtbot, widget = widget_setup(qtbot, RadioButtonBox, 'title', actions)
    widget.set_selection('action2')
    assert widget.checkedButton().text() == 'action2'


def test_RadioButtonBox_add_buttons(qtbot):
    assert False


def test_RadioButtonBox_add_buttons_indented(qtbot):
    assert False


def test_RadioButtonBox_emit(qtbot):
    assert False


def test_QMultiInput_dialog(qtbot):
    assert False
