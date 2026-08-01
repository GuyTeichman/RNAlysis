from pathlib import Path

import pytest

DOCS_SOURCE = Path(__file__).resolve().parent.parent / 'docs' / 'source'

# Since the v4.0 migration off pandas, RNAlysis returns polars objects (see issue #100). The
# hand-written user guides must not describe function outputs or underlying data structures as
# pandas DataFrames/Series, or users will expect a pandas API that no longer exists.
STALE_PANDAS_PHRASES = ('pandas DataFrame', 'pandas Series')


@pytest.mark.parametrize('guide', ['user_guide.rst', 'user_guide_gui.rst'])
def test_user_guides_describe_polars_output_not_pandas(guide):
    text = (DOCS_SOURCE / guide).read_text(encoding='utf-8')
    found = [phrase for phrase in STALE_PANDAS_PHRASES if phrase in text]
    assert not found, f"{guide} still describes outputs as {found}; RNAlysis returns polars objects since v4.0"
