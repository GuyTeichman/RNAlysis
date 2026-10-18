"""
Dependency-free docstring parsing shared by the GUI and the developer tooling.

This module must only import from the Python standard library: the ``check-gui-docstrings`` pre-commit hook
(``packaging/check_gui_docstrings.py``) loads it by file path, in an environment where none of RNAlysis's
dependencies (polars, numpy, PyQt6, ...) are installed and where importing the ``rnalysis`` package is impossible.
``rnalysis.utils.parsing.parse_docstring`` delegates here, so the GUI and the hook share one parser.
"""

import re
from typing import Dict, Tuple

# a parameter's description is any text (non-greedy) up to the next reST field line (':type x:', ':return:',
# ':Examples:', etc.), the end of its paragraph, or the end of the docstring.
# Field lines may be indented by stray whitespace (a common docstring typo), so that is tolerated too.
_PARAM_PATTERN = re.compile(r'^[ \t]*:param ([a-zA-Z_0-9]+):([\s\S]*?)(?=^[ \t]*:.*:|\n\n|\Z)', re.MULTILINE)


def parse_docstring(docstring: str) -> Tuple[str, Dict[str, str]]:
    """
    Parse a given docstring (str) to retrieve the description text, as well as a dictionary of parameter descriptions.

    :param docstring: the docstring to be parsed
    :type docstring: str
    :return: a string matching the description, and a dictionary containing the parameter descriptions
    """
    docstring = re.sub(' +', ' ', docstring)
    desc, _, params_str = docstring.partition('\n\n')
    params = {match.group(1): match.group(2).replace('. ', '. \n') for match in _PARAM_PATTERN.finditer(params_str)}
    return desc, params
