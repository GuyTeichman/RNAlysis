import re
from pathlib import Path

import pytest

DOCS_SOURCE = Path(__file__).resolve().parent.parent / 'docs' / 'source'

# Since the v4.0 migration off pandas, RNAlysis returns polars objects (see issue #100). The
# hand-written user guides must not describe function outputs or underlying data structures as
# pandas DataFrames/Series, or users will expect a pandas API that no longer exists.
STALE_PANDAS_PHRASES = ('pandas DataFrame', 'pandas Series')

# The example REPL outputs in the guides must be shown in polars repr, not pandas repr (issue #110).
# These markers are unambiguous artifacts of pandas' DataFrame/Series repr and never appear in polars'.
PANDAS_REPR_FOOTER = re.compile(r'\[\d+ rows x \d+ columns\]')


@pytest.mark.parametrize('guide', ['user_guide.rst', 'user_guide_gui.rst'])
def test_user_guides_describe_polars_output_not_pandas(guide):
    text = (DOCS_SOURCE / guide).read_text(encoding='utf-8')
    found = [phrase for phrase in STALE_PANDAS_PHRASES if phrase in text]
    assert not found, f'{guide} still describes outputs as {found}; RNAlysis returns polars objects since v4.0'


@pytest.mark.parametrize('guide', ['user_guide.rst', 'user_guide_gui.rst'])
def test_user_guide_example_outputs_use_polars_repr(guide):
    text = (DOCS_SOURCE / guide).read_text(encoding='utf-8')
    assert '<BLANKLINE>' not in text, f"{guide} still contains a pandas '<BLANKLINE>' repr marker (issue #110)"
    assert not PANDAS_REPR_FOOTER.search(text), (
        f"{guide} still contains a pandas '[N rows x M columns]' repr footer; "
        f'regenerate the example outputs in polars repr (issue #110)'
    )
