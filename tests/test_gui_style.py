import builtins

import pytest

from rnalysis.gui import gui_style
from rnalysis.utils import settings


def test_get_stylesheet_names(monkeypatch):
    sheets = {'first': None, 'second': 'val', 'third': 3}
    truth = {'First': 'first', 'Second': 'second', 'Third': 'third'}
    monkeypatch.setattr(gui_style, 'STYLESHEETS', sheets)
    assert gui_style.get_stylesheet_names() == truth


@pytest.mark.parametrize(
    'fontsize,fontname,truth',
    [
        (
            18,
            'Arial',
            """QWidget {
font: Arial;
font-size: 18pt;
}

MainWindow{
font: Arial;
font-size: 36pt;
}

QToolButton {
font: bold;
font-size: 27pt;
}""",
        ),
        (
            7,
            'Times New Roman',
            """QWidget {
font: Times New Roman;
font-size: 7pt;
}

MainWindow{
font: Times New Roman;
font-size: 14pt;
}

QToolButton {
font: bold;
font-size: 10pt;
}""",
        ),
        (
            50,
            'David',
            """QWidget {
font: David;
font-size: 50pt;
}

MainWindow{
font: David;
font-size: 100pt;
}

QToolButton {
font: bold;
font-size: 75pt;
}""",
        ),
    ],
)
def test_get_parametric_stylesheet(monkeypatch, fontsize, fontname, truth):
    param_stylesheet = """QWidget {
font: $FONTPLACEHOLDER;
font-size: $FONTSIZEPLACEHOLDER*1;
}

MainWindow{
font: $FONTPLACEHOLDER;
font-size: $FONTSIZEPLACEHOLDER*2;
}

QToolButton {
font: bold;
font-size: $FONTSIZEPLACEHOLDER*1.5;
}"""

    class MockIO:
        def __init__(self):
            pass

        def read(self, n=1):
            return param_stylesheet

        def __enter__(self):
            return self

        def __exit__(self, exc_type, exc_val, exc_tb):
            pass

    def mock_open(pth, mode='r'):
        return MockIO()

    monkeypatch.setattr(builtins, 'open', mock_open)
    stylesheet = gui_style.get_parametric_stylesheet(fontsize, fontname)
    assert stylesheet == truth


def test_primary_button_uses_material_blue():
    # The single main "go" action of each screen is a filled primary button; its fill is
    # Material Design Blue 700 with white text (AA-safe), replacing the previous warm gold.
    # Reads the real parametric_style.qss (open is not mocked here).
    stylesheet = gui_style.get_parametric_stylesheet(12, 'Arial').lower()
    assert 'qpushbutton[class="primary"]' in stylesheet
    assert '#1976d2' in stylesheet  # Blue 700 base fill
    assert '#1565c0' in stylesheet  # Blue 800 hover fill
    assert '#0d47a1' in stylesheet  # Blue 900 pressed fill
    assert 'color: #ffffff' in stylesheet  # white label text on the blue fill
    # the old gold accent must be gone from the primary-button styling
    assert '#f2b134' not in stylesheet


def test_get_stylesheet(monkeypatch):
    param_stylesheet = '1\n2\n3\n\n'
    monkeypatch.setattr(gui_style, 'get_parametric_stylesheet', lambda size, name: param_stylesheet)
    monkeypatch.setattr(
        settings, 'get_gui_settings', lambda: ('Arial', 18, 'light', ['UniProtKB', 'PomBase'], True, None)
    )

    stylesheet = gui_style.get_stylesheet()

    assert stylesheet.startswith(param_stylesheet + '\n')
    assert len(stylesheet) > len(param_stylesheet + '\n')
