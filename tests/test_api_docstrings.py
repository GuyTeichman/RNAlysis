import re
from pathlib import Path

import pytest

RNALYSIS = Path(__file__).resolve().parent.parent / 'rnalysis'
PUBLIC_API_MODULES = ('filtering.py', 'enrichment.py', 'general.py')

# The public-API docstrings become the auto-generated API reference and the GUI help. Since the v4.0
# migration off pandas, their :Examples: outputs must be shown in polars repr, not pandas (issue #146,
# the docstring analog of the guide fix #110). These markers are unambiguous pandas-repr artifacts.
PANDAS_REPR_FOOTER = re.compile(r'\[\d+ rows x \d+ columns\]')

# A fixture path under tests/ that is not the real tests/test_files/ location is a broken example that
# raises FileNotFoundError if a user copies it. All bundled fixtures live under tests/test_files/.
BROKEN_FIXTURE_PATH = re.compile(r'tests/(?!test_files/)[\w./-]+\.csv')


@pytest.mark.parametrize('module', PUBLIC_API_MODULES)
def test_public_api_docstrings_use_polars_repr(module):
    text = (RNALYSIS / module).read_text(encoding='utf-8')
    assert '<BLANKLINE>' not in text, f"{module} docstrings still contain a pandas '<BLANKLINE>' repr marker"
    assert not PANDAS_REPR_FOOTER.search(text), (
        f"{module} docstrings still contain a pandas '[N rows x M columns]' repr footer; "
        f'regenerate the :Examples: outputs in polars repr (issue #146)'
    )


@pytest.mark.parametrize('module', PUBLIC_API_MODULES)
def test_public_api_docstring_example_paths_exist(module):
    text = (RNALYSIS / module).read_text(encoding='utf-8')
    broken = sorted(set(BROKEN_FIXTURE_PATH.findall(text)))
    assert not broken, (
        f'{module} docstring examples reference nonexistent fixture paths {broken} (use tests/test_files/)'
    )
