import pytest
import joblib
from rnalysis.gui.gui_widgets import *
import time
from typing_extensions import Literal
from typing import Set
from rnalysis.utils import io

LEFT_CLICK = QtCore.Qt.LeftButton
RIGHT_CLICK = QtCore.Qt.RightButton


def widget_setup(qtbot, widget_class, *args, **kwargs):
    widget = widget_class(*args, **kwargs)
    widget.show()
    qtbot.add_widget(widget)
    return qtbot, widget


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


def test_MandatoryComboBox_disable(qtbot):
    default = 'default_chocie'
    qtbot, widget = widget_setup(qtbot, MandatoryComboBox, default)
    widget.addItems(['a', 'b', 'c'])
    widget.setDisabled(True)
    assert not widget.isEnabled()

    widget.setEnabled(True)
    assert widget.isEnabled()

    widget.setEnabled(False)
    assert not widget.isEnabled()


def test_clear_layout(qtbot):
    qtbot, widget = widget_setup(qtbot, QtWidgets.QWidget)
    layout = QtWidgets.QGridLayout(widget)
    layout.addWidget(QtWidgets.QSpinBox(), 0, 0)
    layout.addWidget(QtWidgets.QLineEdit(), 1, 2)
    layout.addWidget(QtWidgets.QLabel("test"), 3, 3)

    clear_layout(layout)

    assert layout.count() == 0


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
    (ToggleSwitch, False, 'switch', True, {}),
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
    (QMultiLineEdit, ['', 'text', 'other text'], [], ['', 'text', 'other text'], {}),
    (QMultiStrIntLineEdit, ['3', '-7', 'text', 'othertext', 'param5'], [], [3, -7, 'text', 'othertext', 'param5'], {}),
    (QMultiBoolComboBox, [True, True, False, True], True, [True, True, False, True], {}),
    (MultiColorPicker, ['r', 'black', '#0000ff'], None, ['#ff0000', '#000000', '#0000ff'], {}),
    (QMultiComboBox, ['option3', 'option2', 'option2'], 'option1', ['option3', 'option2', 'option2'],
     {'items': ['option1', 'option2', 'option3']}),
    (TrueFalseBoth, [True, False], [], [True, False], None)
])
def test_get_val_from_widget_multiinput_types(qtbot, widget_class, default, excepted_val_empty, expected_val, kwargs):
    if kwargs is None:
        qtbot, widget = widget_setup(qtbot, widget_class)
    else:
        qtbot, widget = widget_setup(qtbot, widget_class, label='label', **kwargs)
    assert get_val_from_widget(widget) == excepted_val_empty

    widget.set_defaults(default)
    assert get_val_from_widget(widget) == expected_val


@pytest.mark.parametrize("widget_class", (QtWidgets.QWidget, QtWidgets.QDateTimeEdit))
def test_get_val_from_widget_bad_widget(qtbot, widget_class):
    qtbot, widget = widget_setup(qtbot, widget_class)
    with pytest.raises(TypeError):
        get_val_from_widget(widget)


def test_HelpButton_param_help(qtbot, monkeypatch):
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


def test_HelpButton_desc_help(qtbot, monkeypatch):
    desc = 'mydesc'
    help_shown = []

    def mock_show_text(pos, text):
        assert text == desc
        assert pos == QtGui.QCursor.pos()
        help_shown.append(1)

    monkeypatch.setattr(QtWidgets.QToolTip, 'showText', mock_show_text)

    qtbot, widget = widget_setup(qtbot, HelpButton)
    widget.connect_desc_help(desc)

    qtbot.mouseClick(widget, LEFT_CLICK)

    assert help_shown == [1]


def test_ComparisonPicker_init(qtbot):
    _ = widget_setup(qtbot, ComparisonPicker, io.load_csv('tests/test_files/test_design_matrix.csv', 0))


def test_ComparisonPicker_get_value(qtbot):
    qtbot, widget = widget_setup(qtbot, ComparisonPicker, io.load_csv('tests/test_files/test_design_matrix.csv', 0))
    assert widget.get_value() == ('condition', 'cond1', 'cond1')

    widget.factor.setCurrentText('replicate')
    widget.numerator.setCurrentText('rep3')
    widget.denominator.setCurrentText('rep2')

    assert widget.get_value() == ('replicate', 'rep3', 'rep2')


def test_ComparisonPickerGroup_init(qtbot):
    _ = widget_setup(qtbot, ComparisonPickerGroup, io.load_csv('tests/test_files/test_design_matrix.csv', 0))


def test_ComparisonPickerGroup_add_widget(qtbot):
    qtbot, widget = widget_setup(qtbot, ComparisonPickerGroup,
                                 io.load_csv('tests/test_files/test_design_matrix.csv', 0))
    widget.add_comparison_widget()
    assert len(widget.inputs) == 2
    widget.add_comparison_widget()
    assert len(widget.inputs) == 3


def test_ComparisonPickerGroup_remove_widget(qtbot):
    qtbot, widget = widget_setup(qtbot, ComparisonPickerGroup,
                                 io.load_csv('tests/test_files/test_design_matrix.csv', 0))
    widget.remove_comparison_widget()
    assert len(widget.inputs) == 0
    widget.remove_comparison_widget()
    assert len(widget.inputs) == 0


def test_ComparisonPickerGroup_get_comparison_values(qtbot):
    qtbot, widget = widget_setup(qtbot, ComparisonPickerGroup,
                                 io.load_csv('tests/test_files/test_design_matrix.csv', 0))
    widget.add_comparison_widget()

    widget.inputs[0].factor.setCurrentText('replicate')
    widget.inputs[0].numerator.setCurrentText('rep3')
    widget.inputs[0].denominator.setCurrentText('rep2')
    assert widget.get_comparison_values() == [('replicate', 'rep3', 'rep2'), ('condition', 'cond1', 'cond1')]


def test_PathLineEdit_disable(qtbot):
    qtbot, widget = widget_setup(qtbot, PathLineEdit)

    widget.setDisabled(True)
    assert not widget.isEnabled()

    widget.setEnabled(True)
    assert widget.isEnabled()

    widget.setEnabled(False)
    assert not widget.isEnabled()


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


@pytest.mark.parametrize('path,is_legal_expected', [
    ('path/that/doesnt/exist', False),
    ('tests/test_files/test_deseq.csv', False),
    ('tests/test_files', True),
    ('path/that/doesnt/exist.png', False)])
def test_PathLineEdit_is_legal_dirmode(qtbot, path, is_legal_expected):
    qtbot, widget = widget_setup(qtbot, PathLineEdit, is_file=False)

    assert not widget.is_legal

    widget.clear()
    qtbot.keyClicks(widget.file_path, path)
    assert widget.is_legal == is_legal_expected


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


def test_PathLineEdit_choose_folder(qtbot, monkeypatch):
    pth = 'path/to/a/folder'

    def mock_choose_folder(*args, **kwargs):
        return pth

    monkeypatch.setattr(QtWidgets.QFileDialog, 'getExistingDirectory', mock_choose_folder)
    qtbot, widget = widget_setup(qtbot, PathLineEdit, is_file=False)
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

    selected_items = [item.text() for item in widget.get_sorted_selection()]
    assert selected_items == items


def test_MultipleChoiceList_clear_selection(qtbot):
    items = ['item1', 'item2', 'item3']
    qtbot, widget = widget_setup(qtbot, MultipleChoiceList, items)
    for item in widget.list_items:
        item.setSelected(True)

    assert len(widget.get_sorted_selection()) == len(items)

    qtbot.mouseClick(widget.clear_all_button, LEFT_CLICK)

    selected_items = widget.get_sorted_selection()
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


def test_TableColumnPicker_select_all(qtbot):
    cols = ['a', 'b', 'c', 'd', 'e']
    qtbot, widget = widget_setup(qtbot, TableColumnPicker)
    widget.add_columns(cols)

    changed = []
    widget.valueChanged.connect(functools.partial(changed.append, True))

    qtbot.mouseClick(widget.clear_button, LEFT_CLICK)
    qtbot.mouseClick(widget.select_all_button, LEFT_CLICK)
    qtbot.mouseClick(widget.done_button, LEFT_CLICK)
    assert widget.get_values() == cols
    assert len(changed) == 1


def test_TableColumnPicker_clear_selection(qtbot):
    cols = ['a', 'b', 'c', 'd', 'e']
    qtbot, widget = widget_setup(qtbot, TableColumnPicker)
    widget.add_columns(cols)

    changed = []
    widget.valueChanged.connect(functools.partial(changed.append, True))

    qtbot.mouseClick(widget.clear_button, LEFT_CLICK)
    qtbot.mouseClick(widget.done_button, LEFT_CLICK)
    assert widget.get_values() == []
    assert len(changed) == 1


@pytest.mark.parametrize('selections', [['e'], ['b', 'e'], ['a', 'b', 'c']])
def test_TableColumnPicker_custom_selection(qtbot, selections):
    cols = ['a', 'b', 'c', 'd', 'e']
    qtbot, widget = widget_setup(qtbot, TableColumnPicker)
    widget.add_columns(cols)

    changed = []
    widget.valueChanged.connect(functools.partial(changed.append, True))

    qtbot.mouseClick(widget.clear_button, LEFT_CLICK)
    for selection in selections:
        ind = cols.index(selection)
        qtbot.mouseClick(widget.column_checks[ind].switch, LEFT_CLICK)
    qtbot.mouseClick(widget.done_button, LEFT_CLICK)
    assert widget.get_values() == selections
    assert len(changed) == 1


@pytest.mark.parametrize('selections', [['e'], ['b', 'e'], ['a', 'b', 'c']])
def test_TableSingleColumnPicker_custom_selection(qtbot, selections):
    cols = ['a', 'b', 'c', 'd', 'e']
    qtbot, widget = widget_setup(qtbot, TableSingleColumnPicker)
    widget.add_columns(cols)

    changed = []
    widget.valueChanged.connect(functools.partial(changed.append, True))

    for selection in selections:
        ind = cols.index(selection)
        qtbot.mouseClick(widget.column_checks[ind].switch, LEFT_CLICK)
    qtbot.mouseClick(widget.done_button, LEFT_CLICK)
    assert widget.get_values() == selections[-1]
    assert len(changed) == 1


def test_TableColumnGroupPicker_select_all(qtbot):
    cols = ['a', 'b', 'c', 'd', 'e']
    truth = [[item] for item in cols]
    qtbot, widget = widget_setup(qtbot, TableColumnGroupPicker)
    widget.add_columns(cols)

    changed = []
    widget.valueChanged.connect(functools.partial(changed.append, True))

    qtbot.mouseClick(widget.clear_button, LEFT_CLICK)
    qtbot.mouseClick(widget.select_all_button, LEFT_CLICK)
    qtbot.mouseClick(widget.done_button, LEFT_CLICK)
    assert widget.get_values() == truth
    assert len(changed) == 1


def test_TableColumnGroupPicker_clear_selection(qtbot):
    cols = ['a', 'b', 'c', 'd', 'e']
    qtbot, widget = widget_setup(qtbot, TableColumnGroupPicker)
    widget.add_columns(cols)

    changed = []
    widget.valueChanged.connect(functools.partial(changed.append, True))

    qtbot.mouseClick(widget.clear_button, LEFT_CLICK)
    qtbot.mouseClick(widget.done_button, LEFT_CLICK)
    assert widget.get_values() == []
    assert len(changed) == 1


@pytest.mark.parametrize('selections', [['e'], ['b', 'e'], ['a', 'b', 'c']])
def test_TableColumnGroupPicker_custom_selection(qtbot, selections):
    cols = ['a', 'b', 'c', 'd', 'e']
    qtbot, widget = widget_setup(qtbot, TableColumnGroupPicker)
    widget.add_columns(cols)
    truth = [[item] for item in selections]

    changed = []
    widget.valueChanged.connect(functools.partial(changed.append, True))

    qtbot.mouseClick(widget.clear_button, LEFT_CLICK)
    for selection in selections:
        ind = cols.index(selection)
        qtbot.mouseClick(widget.column_checks[ind].switch, LEFT_CLICK)
    qtbot.mouseClick(widget.done_button, LEFT_CLICK)
    assert widget.get_values() == truth
    assert len(changed) == 1


@pytest.mark.parametrize('selections', [[['d', 'e']], [['a'], ['b'], ['c', 'd', 'e']], [['a', 'b'], ['c']]])
def test_TableColumnGroupPicker_custom_selection_grouped(qtbot, selections):
    cols = ['a', 'b', 'c', 'd', 'e']
    qtbot, widget = widget_setup(qtbot, TableColumnGroupPicker)
    widget.add_columns(cols)

    changed = []
    widget.valueChanged.connect(functools.partial(changed.append, True))

    qtbot.mouseClick(widget.clear_button, LEFT_CLICK)
    for i, grp in enumerate(selections):
        for selection in grp:
            ind = cols.index(selection)
            qtbot.mouseClick(widget.column_checks[ind].switch, LEFT_CLICK)
            qtbot.keyClicks(widget.column_combos[ind], str(i + 1))

    qtbot.mouseClick(widget.done_button, LEFT_CLICK)
    assert widget.get_values() == selections
    assert len(changed) == 1

    widget.reset()

    assert widget.get_values() == [[item] for item in cols]


def test_PathInputDialog(qtbot):
    pth = str(Path('tests/test_files/test_deseq.csv').absolute())
    qtbot, widget = widget_setup(qtbot, PathInputDialog)
    widget.path.clear()
    qtbot.keyClicks(widget.path.file_path, pth)
    qtbot.mouseClick(widget.button_box.buttons()[0], LEFT_CLICK)
    assert widget.result() == pth


def test_ToggleSwitch(qtbot):
    qtbot, widget = widget_setup(qtbot, ToggleSwitch)
    assert not widget.isChecked()
    qtbot.mouseClick(widget.switch, LEFT_CLICK)
    assert widget.isChecked()


def test_ToggleSwitchCore(qtbot):
    qtbot, widget = widget_setup(qtbot, ToggleSwitchCore)
    assert not widget.isChecked()
    qtbot.mouseClick(widget, LEFT_CLICK)
    assert widget.isChecked()

    widget.paintEvent(QtGui.QPaintEvent(QtCore.QRect(0, 0, 1, 1)))


def test_MultiChoiceListWithDelete_delete_all(qtbot, monkeypatch):
    monkeypatch.setattr(QtWidgets.QMessageBox, "question", lambda *args: QtWidgets.QMessageBox.Yes)

    items = ['item1', 'item2', 'item3']
    qtbot, widget = widget_setup(qtbot, MultiChoiceListWithDelete, items)
    qtbot.mouseClick(widget.delete_all_button, LEFT_CLICK)

    selected_items = [item.text() for item in widget.get_sorted_selection()]
    assert selected_items == []
    assert widget.items == []


def test_MultiChoiceListWithDelete_delete(qtbot):
    items = ['item1', 'item2', 'item3']
    qtbot, widget = widget_setup(qtbot, MultiChoiceListWithDelete, items)

    qtbot.mouseClick(widget.delete_button, LEFT_CLICK)

    assert widget.items == items

    qtbot.mouseClick(widget.select_all_button, LEFT_CLICK)
    qtbot.mouseClick(widget.delete_button, LEFT_CLICK)

    selected_items = [item.text() for item in widget.get_sorted_selection()]
    assert selected_items == []
    assert widget.items == []


def test_MultiChoiceListWithDelete_select_all(qtbot):
    items = ['item1', 'item2', 'item3']
    qtbot, widget = widget_setup(qtbot, MultiChoiceListWithDelete, items)
    qtbot.mouseClick(widget.select_all_button, LEFT_CLICK)

    selected_items = [item.text() for item in widget.get_sorted_selection()]
    assert selected_items == items


def test_MultiChoiceListWithDelete_clear_selection(qtbot):
    items = ['item1', 'item2', 'item3']
    qtbot, widget = widget_setup(qtbot, MultiChoiceListWithDelete, items)
    for item in widget.list_items:
        item.setSelected(True)

    assert len(widget.get_sorted_selection()) == len(items)

    qtbot.mouseClick(widget.clear_all_button, LEFT_CLICK)

    selected_items = widget.get_sorted_selection()
    assert len(selected_items) == 0


def test_MultiChoiceListWithDelete_emit(qtbot):
    selection_changed = []

    def selection_changed_slot(*args, **kwargs):
        selection_changed.append(1)

    items = ['item1', 'item2', 'item3']
    qtbot, widget = widget_setup(qtbot, MultiChoiceListWithDelete, items)
    widget.itemSelectionChanged.connect(selection_changed_slot)

    assert len(selection_changed) == 0

    qtbot.mouseClick(widget.select_all_button, LEFT_CLICK)
    assert selection_changed == [1] * len(items)
    qtbot.mouseClick(widget.clear_all_button, LEFT_CLICK)
    assert selection_changed == [1] * len(items) * 2
    widget.list_items[0].setSelected(True)
    assert selection_changed == [1] * len(items) * 2 + [1]


def test_MultiChoiceListWithReorder_up(qtbot):
    items = ['item1', 'item2', 'item3']
    qtbot, widget = widget_setup(qtbot, MultiChoiceListWithReorder, items)

    assert widget.items == items

    widget.list_items[-1].setSelected(True)
    qtbot.mouseClick(widget.up_button, LEFT_CLICK)

    selected_items = [item.text() for item in widget.get_sorted_selection()]
    assert selected_items == ['item3']
    assert widget.items == ['item1', 'item3', 'item2']

    for i in range(3):
        qtbot.mouseClick(widget.up_button, LEFT_CLICK)
        selected_items = [item.text() for item in widget.get_sorted_selection()]
        assert selected_items == ['item3']
        assert widget.items == ['item3', 'item1', 'item2']


def test_MultiChoiceListWithReorder_up_multi(qtbot):
    items = ['item1', 'item2', 'item3', 'item4']
    qtbot, widget = widget_setup(qtbot, MultiChoiceListWithReorder, items)

    assert widget.items == items

    widget.list_items[1].setSelected(True)
    widget.list_items[-1].setSelected(True)
    qtbot.mouseClick(widget.up_button, LEFT_CLICK)

    selected_items = [item.text() for item in widget.get_sorted_selection()]
    assert selected_items == ['item2', 'item4']
    assert widget.items == ['item2', 'item1', 'item4', 'item3']

    qtbot.mouseClick(widget.up_button, LEFT_CLICK)
    selected_items = [item.text() for item in widget.get_sorted_selection()]
    assert selected_items == ['item2', 'item4']
    assert widget.items == ['item2', 'item4', 'item1', 'item3']


def test_MultiChoiceListWithReorder_top(qtbot):
    items = ['item1', 'item2', 'item3', 'item4']
    qtbot, widget = widget_setup(qtbot, MultiChoiceListWithReorder, items)

    assert widget.items == items

    widget.list_items[1].setSelected(True)
    widget.list_items[-1].setSelected(True)
    qtbot.mouseClick(widget.top_button, LEFT_CLICK)

    selected_items = [item.text() for item in widget.get_sorted_selection()]
    assert selected_items == ['item2', 'item4']
    assert widget.items == ['item2', 'item4', 'item1', 'item3']


def test_MultiChoiceListWithReorder_down(qtbot):
    items = ['item1', 'item2', 'item3']
    qtbot, widget = widget_setup(qtbot, MultiChoiceListWithReorder, items)

    assert widget.items == items

    widget.list_items[0].setSelected(True)
    qtbot.mouseClick(widget.down_button, LEFT_CLICK)

    selected_items = [item.text() for item in widget.get_sorted_selection()]
    assert selected_items == ['item1']
    assert widget.items == ['item2', 'item1', 'item3']

    for i in range(3):
        qtbot.mouseClick(widget.down_button, LEFT_CLICK)
        selected_items = [item.text() for item in widget.get_sorted_selection()]
        assert selected_items == ['item1']
        assert widget.items == ['item2', 'item3', 'item1']


def test_MultiChoiceListWithReorder_down_multi(qtbot):
    items = ['item1', 'item2', 'item3', 'item4']
    qtbot, widget = widget_setup(qtbot, MultiChoiceListWithReorder, items)

    assert widget.items == items

    widget.list_items[0].setSelected(True)
    widget.list_items[2].setSelected(True)
    qtbot.mouseClick(widget.down_button, LEFT_CLICK)

    selected_items = [item.text() for item in widget.get_sorted_selection()]
    assert selected_items == ['item1', 'item3']
    assert widget.items == ['item2', 'item1', 'item4', 'item3']

    qtbot.mouseClick(widget.down_button, LEFT_CLICK)
    selected_items = [item.text() for item in widget.get_sorted_selection()]
    assert selected_items == ['item1', 'item3']
    assert widget.items == ['item2', 'item4', 'item1', 'item3']


def test_MultiChoiceListWithReorder_bottom(qtbot):
    items = ['item1', 'item2', 'item3', 'item4']
    qtbot, widget = widget_setup(qtbot, MultiChoiceListWithReorder, items)

    assert widget.items == items

    widget.list_items[0].setSelected(True)
    widget.list_items[2].setSelected(True)
    qtbot.mouseClick(widget.bottom_button, LEFT_CLICK)

    selected_items = [item.text() for item in widget.get_sorted_selection()]
    assert selected_items == ['item1', 'item3']
    assert widget.items == ['item2', 'item4', 'item1', 'item3']


def test_OrderedFileList_init(qtbot):
    _ = widget_setup(qtbot, OrderedFileList)


def test_OrderedFileList_add_item(qtbot):
    pth = Path('tests/test_files/test_deseq.csv')
    qtbot, widget = widget_setup(qtbot, OrderedFileList)
    widget.add_item(pth)
    assert widget.items == [pth]
    assert len(widget.list_items) == 1
    assert isinstance(widget.list_items[0], FileListWidgetItem)
    assert widget.list_items[0].file_path == pth
    assert widget.list_items[0].display_name == pth.stem
    assert widget.list_items[0].filePath() == pth
    assert widget.list_items[0].text() == pth.stem
    assert widget.list.count() == 1


@pytest.mark.parametrize("filenames", [
    ['tests/test_files/test_deseq.csv'],
    [],
    ['tests/test_files/test_deseq.csv', 'tests/test_files/counted.tsv'],
    ['tests/test_files/a.fa', 'tests/test_files/c.fa', 'tests/test_files/b.fa', 'aa.fa']
])
def test_OrderedFileList_add_files(qtbot, monkeypatch, filenames):
    def mock_get_file_names(*args, **kwargs):
        return filenames, None

    monkeypatch.setattr(QtWidgets.QFileDialog, 'getOpenFileNames', mock_get_file_names)
    qtbot, widget = widget_setup(qtbot, OrderedFileList)
    widget.add_files()

    assert widget.list.count() == len(filenames)
    assert widget.items == filenames


@pytest.mark.parametrize("filenames", [
    ['tests/test_files/test_deseq.csv'],
    [],
    ['tests/test_files/test_deseq.csv', 'tests/test_files/counted.tsv'],
    ['tests/test_files/a.fa', 'tests/test_files/c.fa', 'tests/test_files/b.fa', 'aa.fa']

])
def test_OrderedFileList_get_sorted_selected_names(qtbot, monkeypatch, filenames):
    def mock_get_file_names(*args, **kwargs):
        return filenames, None

    monkeypatch.setattr(QtWidgets.QFileDialog, 'getOpenFileNames', mock_get_file_names)

    qtbot, widget = widget_setup(qtbot, OrderedFileList)
    widget.add_files()
    res = widget.get_sorted_names()

    assert res == filenames


def test_TrueFalseBoth(qtbot):
    selection_changed = []
    qtbot, widget = widget_setup(qtbot, TrueFalseBoth)
    widget.selectionChanged.connect(functools.partial(selection_changed.append, True))
    assert widget.get_values() == []
    assert selection_changed == []

    qtbot.mouseClick(widget.true_button, LEFT_CLICK)
    assert widget.get_values() == [True]
    assert selection_changed == [True]

    qtbot.mouseClick(widget.false_button, LEFT_CLICK)
    assert widget.get_values() == [True, False]
    assert selection_changed == [True, True]

    qtbot.mouseClick(widget.true_button, LEFT_CLICK)
    assert widget.get_values() == [False]
    assert selection_changed == [True, True, True]


def test_StdOutTextEdit_print(qtbot):
    qtbot, widget = widget_setup(qtbot, StdOutTextEdit)
    widget.append_text('hello world')
    assert widget.toPlainText() == 'hello world' + '\n'
    widget.append_text('<tag> </tag>')
    assert widget.toPlainText() == 'hello world' + '\n' + '<tag> </tag>' + '\n'


def test_StdOutTextEdit_warnings(qtbot):
    qtbot, widget = widget_setup(qtbot, StdOutTextEdit)
    widget.append_text('Warning: warning text goes here')
    assert widget.toPlainText() == 'Warning: warning text goes here' + '\n'


def test_StdOutTextEdit_carriage(qtbot):
    line1, line2, line3 = 'first line', 'second line\r', 'last line'
    qtbot, widget = widget_setup(qtbot, StdOutTextEdit)
    widget.append_text(line1)
    widget.append_text(line2)
    assert widget.toPlainText() == line1 + '\n' + 'second line ' + '\n'
    widget.append_text(line3)
    assert widget.toPlainText() == line1 + '\n' + line3 + '\n'


def test_StdOutTextEdit_empty_line(qtbot):
    qtbot, widget = widget_setup(qtbot, StdOutTextEdit)
    for i in range(10):
        widget.append_text('\n')
    assert widget.toPlainText() == ''


def test_NewParam():
    _ = NewParam('annotation')
    _ = NewParam('annotation', 'default')


def test_RadioButtonBox_set_selection(qtbot):
    actions = ['action1', 'action2', 'action3']
    qtbot, widget = widget_setup(qtbot, RadioButtonBox, 'title', actions)
    widget.set_selection('action2')
    assert widget.checkedButton().text() == 'action2'


def test_RadioButtonBox_add_buttons(qtbot):
    actions = ['action1', 'action2', 'action3']
    qtbot, widget = widget_setup(qtbot, RadioButtonBox, 'title', [])
    widget.add_items(actions)
    assert len(widget.radio_buttons) == len(actions)


def test_RadioButtonBox_add_buttons_indented(qtbot):
    actions = [('set1', ('action1', 'action2')), ('set2', ('action3',))]
    qtbot, widget = widget_setup(qtbot, RadioButtonBox, 'title', [])
    widget.add_items(actions)
    assert len(widget.radio_buttons) == 3


def test_RadioButtonBox_emit(qtbot):
    actions = ['action1', 'action2', 'action3']
    qtbot, widget = widget_setup(qtbot, RadioButtonBox, 'title', actions)
    with qtbot.waitSignal(widget.selectionChanged, timeout=1000) as blocker:
        widget.set_selection('action2')
    assert widget.checkedButton().text() == 'action2'

    with qtbot.waitSignal(widget.buttonClicked, timeout=1000) as blocker:
        qtbot.mouseClick(widget.radio_buttons['action1'], LEFT_CLICK)
    assert widget.checkedButton().text() == 'action1'
    assert blocker.args[0] == widget.radio_buttons['action1']


@pytest.mark.parametrize('param_type,default,expected_widget', [
    (str, 'default', QtWidgets.QLineEdit),
    (int, 15, QtWidgets.QSpinBox),
    (float, 3.14, QtWidgets.QDoubleSpinBox),
    (Literal['a1', 'b1', 'c1'], 'b1', QtWidgets.QComboBox),
    (Union[str, float, None, bool], True, QtWidgets.QTextEdit),
    (Union[str, float, None, bool], 17, QtWidgets.QTextEdit)
])
def test_param_to_widget_native_types(qtbot, param_type, default, expected_widget):
    _run_param_to_widget(qtbot, param_type, default, 'param_name', expected_widget)


@pytest.mark.parametrize('param_type,default,name,expected_widget', [
    (bool, True, 'status', ToggleSwitch),
    (param_typing.Color, '#000000', 'linecolor', ColorPicker),
    (param_typing.ColorList, ['#000000', '#aabbcc'], 'colors', MultiColorPicker),
    (Union[str, None], None, 'name', OptionalLineEdit),
    (Union[str, None], 'text', 'name', OptionalLineEdit),
    (Union[int, None], None, 'name', OptionalSpinBox),
    (Union[int, None], 5, 'name', OptionalSpinBox),
    (Union[float, None], None, 'name', OptionalDoubleSpinBox),
    (Union[float, None], -0.5, 'name', OptionalDoubleSpinBox),
    (Union[str, List[str]], ['a', 'b'], 'name', QMultiLineEdit),
    (Union[str, Iterable[str]], None, 'name', QMultiLineEdit),
    (Union[int, List[int]], [2, 7], 'name', QMultiSpinBox),
    (Union[int, Iterable[int]], None, 'name', QMultiSpinBox),
    (Union[float, List[float]], [2.5, -0.72], 'name', QMultiDoubleSpinBox),
    (Union[float, Iterable[float]], None, 'name', QMultiDoubleSpinBox),
    (Union[bool, List[bool]], [True, False, True], 'name', QMultiBoolComboBox),
    (Union[bool, Iterable[bool]], None, 'name', QMultiBoolComboBox),
    (List[Literal['a', 'b', 'c']], ['a', 'b', 'a'], 'name', QMultiComboBox),
    (Set[Literal['a', 'b', 'c']], 'c', 'name', QMultiComboBox),
    (Iterable[Literal['a', 'b', 'c']], None, 'name', QMultiComboBox),
    (Union[str, int], 5, 'name', StrIntLineEdit),
    (Union[str, int], 'text', 'name', StrIntLineEdit),
    (Union[str, Path], str(Path('tests/test_files/test_deseq.csv').absolute()), 'name', PathLineEdit),
    (Union[bool, Tuple[bool, bool]], True, 'name', TrueFalseBoth),
    (Union[bool, Tuple[bool, bool]], [True, False], 'name', TrueFalseBoth),
    (Union[str, int, List[str], List[int]], [3, 5, -2], 'name', QMultiStrIntLineEdit),
    (Union[str, int, Iterable[str], Iterable[int]], ['a', 'b', 'c'], 'name', QMultiStrIntLineEdit)

])
def test_param_to_widget_nonnative_types(qtbot, param_type, default, name, expected_widget):
    _run_param_to_widget(qtbot, param_type, default, name, expected_widget)


@pytest.mark.parametrize('param_type,name,expected_widget,expected_widget_pipeline', [
    (param_typing.GroupedColumns, 'samples', TableColumnGroupPicker, TwoLayerMultiLineEdit),
    (param_typing.GroupedColumns, 'sample_grouping', TableColumnGroupPicker, TwoLayerMultiLineEdit),
    (param_typing.ColumnNames, 'sample_names', TableColumnPicker, QMultiLineEdit),
    (param_typing.ColumnNames, 'sample1', TableColumnPicker, QMultiLineEdit),
    (param_typing.ColumnNames, 'sample2', TableColumnPicker, QMultiLineEdit),
    (param_typing.ColumnNames, 'columns', TableColumnPicker, QMultiLineEdit),
    (param_typing.ColumnName, 'by', TableSingleColumnPicker, QtWidgets.QLineEdit),
    (param_typing.ColumnName, 'column', TableSingleColumnPicker, QtWidgets.QLineEdit),

])
def test_param_to_widget_pipeline_mode_types(qtbot, param_type, name, expected_widget, expected_widget_pipeline):
    _run_param_to_widget(qtbot, param_type, 'default', name, expected_widget, expected_widget_pipeline,
                         test_pipeline_mode=True)


@pytest.mark.parametrize('param_type,default,literal_default,expected_sub_widget', [
    (Union[Literal['all'], str], 'text', 'all', QtWidgets.QLineEdit),
    (Union[Literal['any'], Union[str, None]], None, 'any', OptionalLineEdit),
    (Union[Literal['any', 'none'], Union[int, List[int]]], [15, 16], 'any', QMultiSpinBox),

])
def test_param_to_widget_with_literals(qtbot, param_type, default, literal_default, expected_sub_widget):
    expected_widget = ComboBoxOrOtherWidget
    _run_param_to_widget(qtbot, param_type, default, 'name', expected_widget)
    _run_param_to_widget(qtbot, param_type, literal_default, 'name', expected_widget)

    param = NewParam(param_type)
    widget = param_to_widget(param, 'name')
    widget.show()
    qtbot.add_widget(widget)
    assert type(widget.other) == expected_sub_widget


@QtCore.pyqtSlot(name='my_signal')
def _func_to_connect():
    return


def _run_param_to_widget(qtbot, param_type, default, name, expected_widget, expected_widget_pipeline=None,
                         test_pipeline_mode: bool = False):
    param = NewParam(param_type)
    plain_widget = param_to_widget(param, name)
    plain_widget.show()
    qtbot.add_widget(plain_widget)
    assert type(plain_widget) == expected_widget

    param_with_default = NewParam(param_type, default)
    widget = param_to_widget(param_with_default, name, _func_to_connect)
    widget.show()
    qtbot.add_widget(widget)
    assert type(widget) == expected_widget

    if test_pipeline_mode:
        widget = param_to_widget(param, name, pipeline_mode=True)
        widget.show()
        qtbot.add_widget(widget)
        assert type(widget) == expected_widget_pipeline
    else:
        val = get_val_from_widget(widget)
        if isinstance(val, list) and len(val) <= 1:
            if default is None:
                assert len(val) == 0
            else:
                assert (val[0] == default)
        else:
            assert val == default


def test_worker(qtbot):
    def func(a, b, c):
        time.sleep(0.5)
        return a * b * c, a + b + c

    partial = functools.partial(func, 5, 6, 7)
    try:
        worker = Worker(partial, 'other input')

        with qtbot.waitSignal(worker.finished, timeout=4000) as blocker:
            worker.run()
        assert blocker.args == [((5 * 6 * 7, 5 + 6 + 7), 'other input')]
    finally:
        try:
            worker.deleteLater()
        except Exception:
            pass


def test_AltTqdm_init():
    assert list(range(10)) == list(AltTQDM(range(10)))


def test_AltTqdm_iter_signals(qtbot):
    bar_updates = []
    iterator = AltTQDM(range(5))
    iterator.barUpdate.connect(bar_updates.append)
    with qtbot.waitSignal(iterator.barFinished) as blocker:
        for i in iterator:
            assert len(bar_updates) == i + 1
            assert bar_updates[-1] == 1
    assert bar_updates == [1, 1, 1, 1, 1]


def test_AltTqdm_enter_signals():
    bar_updates = []
    bar_finished = []
    with AltTQDM(total=10) as pbar:
        pbar.barUpdate.connect(bar_updates.append)
        pbar.barFinished.connect(functools.partial(bar_finished.append, True))
        for i in range(10):
            if i % 2 == 1:
                pbar.update(2)
                assert len(bar_updates) == (i + 1) // 2
                assert bar_updates[-1] == 2
    assert bar_updates == [2, 2, 2, 2, 2]
    assert bar_finished == [True]


def _power(a, b):
    return a ** b


def test_AltParallel_init():
    res = AltParallel(n_jobs=1, desc='description')(
        joblib.delayed(_power)(a, b) for a, b in zip(range(5), [2, 2, 2, 2, 2]))
    assert res == [0, 1, 4, 9, 16]


def test_AltParallel_signals(qtbot):
    bar_updates = []
    parallel = AltParallel(n_jobs=1, desc='description', total=5)
    parallel.barUpdate.connect(bar_updates.append)
    with qtbot.waitSignal(parallel.barFinished) as blocker:
        res = parallel(joblib.delayed(_power)(a, b) for a, b in zip(range(5), [2, 2, 2, 2, 2]))
    assert bar_updates == [1, 1, 1, 1, 1]
    assert res == [0, 1, 4, 9, 16]


def test_TextWithCopyButton_copy_to_clipboard(qtbot, monkeypatch):
    rich_text = '''<p>
    Auther, F. (2001). The name of the paper. Journal.
    <br>
    <a href=www.google.com>doi.org/10.12345/journalName.496351</a>
    </p>'''
    truth = 'Auther, F. (2001). The name of the paper. Journal. \ndoi.org/10.12345/journalName.496351 '
    qtbot, dialog = widget_setup(qtbot, TextWithCopyButton, rich_text)
    qtbot.mouseClick(dialog.copy_button, LEFT_CLICK)

    cb = QtWidgets.QApplication.clipboard().text()
    assert cb == truth
