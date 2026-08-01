"""Unit tests for packaging/contributors.py (pure text helpers).

The script lives in ``packaging/`` (which is *not* a package and would shadow the PyPI ``packaging``
distribution if imported by name), so it is loaded directly from its file path.
"""
import importlib.util
from pathlib import Path

import pytest

_MOD_PATH = Path(__file__).resolve().parent.parent / "packaging" / "contributors.py"
_spec = importlib.util.spec_from_file_location("rnalysis_contributors_script", _MOD_PATH)
contributors = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(contributors)


README_SAMPLE = """\
Contributors
*************

.. contributors-list-start

* Dror Cohen
* `Mintxoklet <https://github.com/Mintxoklet>`_
* `Bipin Kumar <https://github.com/kbipinkumar>`_

.. contributors-list-end

----
"""


def test_extract_existing_handles_finds_only_linked_handles():
    assert contributors.extract_existing_handles(README_SAMPLE) == {"mintxoklet", "kbipinkumar"}


def test_merge_appends_only_new_logins_case_insensitively():
    new_text, added = contributors.merge_contributor_bullets(
        README_SAMPLE, ["MINTXOKLET", "newperson", "kbipinkumar", "another"])
    assert added == ["newperson", "another"]
    assert "* `newperson <https://github.com/newperson>`_" in new_text
    assert "* `another <https://github.com/another>`_" in new_text
    assert new_text.count("kbipinkumar") == 1        # existing entry not duplicated
    assert "* Dror Cohen" in new_text                # hand-written name preserved
    # new bullets land inside the region, before the end marker
    assert new_text.index("newperson") < new_text.index(contributors.CONTRIB_END)


def test_merge_is_idempotent():
    once, added1 = contributors.merge_contributor_bullets(README_SAMPLE, ["newperson"])
    twice, added2 = contributors.merge_contributor_bullets(once, ["newperson"])
    assert added1 == ["newperson"]
    assert added2 == []
    assert once == twice


def test_merge_dedups_within_input():
    _, added = contributors.merge_contributor_bullets(README_SAMPLE, ["dupe", "dupe", "Dupe"])
    assert added == ["dupe"]


def test_merge_no_new_returns_original_unchanged():
    text, added = contributors.merge_contributor_bullets(README_SAMPLE, ["Mintxoklet"])
    assert added == []
    assert text == README_SAMPLE


def test_merge_raises_without_markers():
    with pytest.raises(ValueError):
        contributors.merge_contributor_bullets("no markers here", ["x"])


def test_format_thanks_sorts_dedups_and_formats():
    line = contributors.format_thanks(["Zed", "amy", "amy", "Bob"])
    assert line == ("Thanks to `@amy <https://github.com/amy>`_, "
                    "`@Bob <https://github.com/Bob>`_, "
                    "`@Zed <https://github.com/Zed>`_ for contributing to this release! \U0001f389")


def test_format_thanks_empty():
    assert contributors.format_thanks([]) == ""
    assert contributors.format_thanks(["", None]) == ""
