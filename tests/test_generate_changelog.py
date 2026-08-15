"""Unit tests for packaging/generate_changelog.py (HISTORY.rst -> latest_changelog.md).

The script lives in ``packaging/`` (which is *not* a package and would shadow the PyPI ``packaging``
distribution if imported by name -- see ``tests/test_contributors.py`` for the same pattern), so it
is loaded directly from its file path.
"""
import importlib.util
from pathlib import Path

import pytest

_MOD_PATH = Path(__file__).resolve().parent.parent / 'packaging' / 'generate_changelog.py'
_spec = importlib.util.spec_from_file_location('rnalysis_generate_changelog_script', _MOD_PATH)
generate_changelog = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(generate_changelog)

DATED_HISTORY = """\
=======
History
=======

1.2.3 (2026-01-01)
-------------------

Fixed
******
* Fixed something new.


1.2.2 (2025-12-01)
-------------------

Fixed
******
* Fixed something older.
"""

UNDATED_TOP_HISTORY = """\
=======
History
=======

1.2.3 (unreleased)
-------------------

Fixed
******
* Fixed something new.


1.2.2 (2025-12-01)
-------------------

Fixed
******
* Fixed something older.
"""


def test_get_change_log_for_latest_returns_newest_dated_section(tmp_path):
    history = tmp_path / 'HISTORY.rst'
    history.write_text(DATED_HISTORY)
    changelog = generate_changelog.get_change_log_for('latest', history_path=history)
    assert '1.2.3' in changelog
    assert 'Fixed something new.' in changelog
    assert '1.2.2' not in changelog


def test_get_change_log_for_latest_raises_when_top_section_is_undated(tmp_path):
    # an undated top section must fail loudly instead of silently generating
    # the previous version's changelog (which the release workflow would then auto-commit)
    history = tmp_path / 'HISTORY.rst'
    history.write_text(UNDATED_TOP_HISTORY)
    with pytest.raises(ValueError, match='unreleased/undated'):
        generate_changelog.get_change_log_for('latest', history_path=history)


def test_get_change_log_for_explicit_version_still_works_with_undated_top(tmp_path):
    # requesting a specific released version is unaffected by an undated top section
    history = tmp_path / 'HISTORY.rst'
    history.write_text(UNDATED_TOP_HISTORY)
    changelog = generate_changelog.get_change_log_for('1.2.2', history_path=history)
    assert 'Fixed something older.' in changelog
