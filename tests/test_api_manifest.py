"""
Rename tripwire for RNAlysis' public API.

Pipeline YAML files store function *names*, and the GUI is generated from the public API by
reflection - so renaming or removing a public function silently breaks every Pipeline a user
exported before the rename, and silently moves or removes a GUI button. Neither is visible in a
code diff.

This module pins the public API surface in a checked-in manifest
(``tests/test_files/public_api_manifest.yaml``) and fails if any pinned name disappears.
Adding new names never fails the test - it only emits a warning asking for the manifest to be
updated. To regenerate the manifest after a deliberate change, run::

    python -m tests.test_api_manifest
"""
import types
import warnings
from pathlib import Path

import pytest
import yaml

from rnalysis import enrichment, fastq, filtering

MANIFEST_PATH = Path(__file__).parent.joinpath('test_files/public_api_manifest.yaml')

#: classes whose public attributes are Pipeline-able and/or reflected into the GUI
MANIFEST_CLASSES = {
    'filtering.Filter': filtering.Filter,
    'filtering.CountFilter': filtering.CountFilter,
    'filtering.DESeqFilter': filtering.DESeqFilter,
    'filtering.FoldChangeFilter': filtering.FoldChangeFilter,
}
#: modules whose public functions are Pipeline-able and/or reflected into the GUI
MANIFEST_MODULES = {'fastq': fastq, 'enrichment': enrichment}

MANIFEST_HEADER = """\
# Public API manifest - a rename tripwire, not a source of truth.
#
# Every name here is part of RNAlysis' public API: it may appear inside a Pipeline YAML file a
# user exported years ago, and it drives what the GUI shows (see the reflection contract in
# CLAUDE.md). tests/test_api_manifest.py fails if any of these names disappears.
#
# Regenerate after a deliberate change with:  python -m tests.test_api_manifest
"""

_RENAME_INSTRUCTIONS = """
A public API name that this manifest pins is gone. If this was a deliberate rename or removal,
it breaks every Pipeline YAML that references the old name (hard rule 4), so before updating the
manifest you must:
  1. Add a back-compat alias so old Pipeline YAMLs referencing the old name keep importing
     (rnalysis/utils/generic.py: GenericPipeline._resolve_function).
  2. Update the pinned Pipeline/session fixtures under tests/test_files/ that use the old name.
  3. Add an entry to HISTORY.rst describing the rename and the alias.
Then regenerate this manifest with:  python -m tests.test_api_manifest
"""


def collect_public_api() -> dict:
    """
    Collect the current public API surface of the manifested classes and modules.

    :return: a dictionary mapping each namespace name to its sorted list of public names
    :rtype: dict
    """
    api = {}
    for namespace_name, cls in MANIFEST_CLASSES.items():
        api[namespace_name] = sorted(name for name in dir(cls)
                                     if not name.startswith('_') and not hasattr(object, name))
    for namespace_name, module in MANIFEST_MODULES.items():
        api[namespace_name] = sorted(name for name in dir(module)
                                     if not name.startswith('_')
                                     and isinstance(getattr(module, name), types.FunctionType)
                                     and getattr(module, name).__module__ == module.__name__)
    return api


def load_manifest() -> dict:
    """
    Load the checked-in public API manifest.

    :return: a dictionary mapping each namespace name to its list of pinned public names
    :rtype: dict
    """
    with open(MANIFEST_PATH) as f:
        return yaml.safe_load(f)


def write_manifest():  # pragma: no cover
    """Regenerate the checked-in public API manifest from the live API."""
    with open(MANIFEST_PATH, 'w') as f:
        f.write(MANIFEST_HEADER)
        yaml.safe_dump(collect_public_api(), f, default_flow_style=False)
    print(f"Wrote public API manifest to '{MANIFEST_PATH}'")


def test_public_api_manifest_has_no_renames_or_removals():
    manifest = load_manifest()
    live = collect_public_api()

    removed = {}
    for namespace_name, pinned_names in manifest.items():
        if namespace_name not in live:
            removed[namespace_name] = sorted(pinned_names)
            continue
        missing = sorted(set(pinned_names) - set(live[namespace_name]))
        if missing:
            removed[namespace_name] = missing

    added = {namespace_name: sorted(set(names) - set(manifest.get(namespace_name, [])))
             for namespace_name, names in live.items()}
    added = {namespace_name: names for namespace_name, names in added.items() if names}

    if added and not removed:
        warnings.warn(f"New public API names are missing from '{MANIFEST_PATH.name}': {added}. "
                      f"Please regenerate the manifest with 'python -m tests.test_api_manifest'.")

    assert not removed, (f"{_RENAME_INSTRUCTIONS}\n"
                         f"Public API names that disappeared: {removed}\n"
                         f"Public API names that were added (a likely rename target): {added or 'none'}")


def test_public_api_manifest_covers_every_namespace():
    manifest = load_manifest()
    assert set(manifest.keys()) == set(collect_public_api().keys()), \
        ("The public API manifest does not cover the same namespaces as the live API. "
         "Regenerate it with 'python -m tests.test_api_manifest'.")


# --- meta-tests: the tripwire itself must actually trip -----------------------------------------


def test_manifest_tripwire_fails_on_a_rename(monkeypatch):
    manifest = load_manifest()
    manifest['filtering.Filter'] = manifest['filtering.Filter'] + ['filter_by_a_former_name']
    monkeypatch.setitem(globals(), 'load_manifest', lambda: manifest)

    with pytest.raises(AssertionError) as err:
        test_public_api_manifest_has_no_renames_or_removals()
    message = str(err.value)
    assert 'filter_by_a_former_name' in message
    assert 'alias' in message and 'fixtures' in message and 'HISTORY.rst' in message


def test_manifest_tripwire_ignores_added_names(monkeypatch):
    manifest = load_manifest()
    manifest['filtering.Filter'] = [name for name in manifest['filtering.Filter'] if name != 'describe']
    monkeypatch.setitem(globals(), 'load_manifest', lambda: manifest)

    with pytest.warns(UserWarning, match='describe'):
        test_public_api_manifest_has_no_renames_or_removals()


if __name__ == '__main__':  # pragma: no cover
    write_manifest()
