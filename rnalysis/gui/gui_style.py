import re
import qdarkstyle
from rnalysis.utils import settings
from pathlib import Path

FONTPLACEHOLDER = "$FONTPLACEHOLDER"
FONTSIZEPLACEHOLDER = "$FONTSIZEPLACEHOLDER"
STYLESHEETS = {'base': None, 'light': qdarkstyle.LightPalette, 'dark': qdarkstyle.DarkPalette}
PARAMETRIC_STYLESHEET_PATH = 'styles/parametric_style.qss'


def get_stylesheet_names():
    return {sheet.capitalize(): sheet for sheet in STYLESHEETS.keys()}


def get_parametric_stylesheet(font_base_size: int, font_name: str):
    assert isinstance(font_base_size, int) and font_base_size >= 1
    assert isinstance(font_name, str)
    with open(Path.joinpath(Path(__file__).parent, PARAMETRIC_STYLESHEET_PATH)) as f:
        style_text = f.read()

    style_text = style_text.replace(FONTPLACEHOLDER, font_name)

    fontsize_pattern = f"\\{FONTSIZEPLACEHOLDER}\\*[0-9]*\\.?[0-9]+"
    for fontsize_match in reversed(list(re.finditer(fontsize_pattern, style_text))):
        multiplier = float(fontsize_match.group(0)[len(FONTSIZEPLACEHOLDER) + 1:])
        start, end = fontsize_match.span()
        style_text = style_text[:start] + f"{int(font_base_size * multiplier)}pt" + style_text[end:]
    style_text = style_text.replace(FONTSIZEPLACEHOLDER, str(font_base_size))
    return style_text


def get_stylesheet():
    font_name, font_base_size, stylesheet_name, _ = settings.get_gui_settings()
    palette = STYLESHEETS[stylesheet_name]
    param_stylesheet = get_parametric_stylesheet(font_base_size, font_name)
    if palette is None:
        other_stylesheet = """
QWidget::item:selected { background-color: #7FABFF; } QWidget::item:hover:!selected { background-color: #53A7EF; }
QSplitter::handle { background-color: #C9CDD0; border: 0px solid #FAFAFA; spacing: 0px; padding: 1px; margin: 0px; }
QSplitter::handle:hover { background-color: #788D9C; }
"""
    else:
        other_stylesheet = qdarkstyle.load_stylesheet(qt_api='pyqt5', palette=palette)
    return param_stylesheet + '\n' + other_stylesheet
