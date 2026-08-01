"""Guards for the packaging metadata.

The optional-dependencies file (``requirements_extra.txt``) must stay a *valid* pip requirements
file: Dependabot's dependency-graph parser reads every ``requirements*.txt`` in the repo, and a
custom grammar (the old ``pkg>=x: tag`` format) makes that job fail with an ``InvalidRequirement``
error. These tests pin both properties: the file parses as pip requirements, and ``setup.py``'s
``get_extra_requires`` still maps every feature tag to the right package(s).
"""
import importlib.util
from pathlib import Path

import pytest
from packaging.requirements import InvalidRequirement, Requirement

REPO_ROOT = Path(__file__).resolve().parent.parent
EXTRAS_FILE = REPO_ROOT / 'requirements_extra.txt'

# feature tag -> set of package names it must pull in (versions deliberately omitted so bumping a
# pin here doesn't break this test). Derived from the intended feature<->package design, not by
# re-running the parser, so it can actually disagree with the code.
EXPECTED_FEATURES = {
    'hdbscan': {'hdbscan'},
    'single-set': {'xlmhglite'},
    'randomization': {'numba'},
    'fastq': {'cutadapt'},
    'reports': {'pyvis', 'networkx'},
    'picardtools': {'install-jdk'},
}


def _load_get_extra_requires():
    """Import ``get_extra_requires`` from setup.py without triggering its ``setup()`` call."""
    spec = importlib.util.spec_from_file_location('rnalysis_setup', REPO_ROOT / 'setup.py')
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module.get_extra_requires


def _requirement_lines():
    for raw in EXTRAS_FILE.read_text(encoding='utf-8').splitlines():
        requirement = raw.partition('#')[0].strip()
        if requirement:
            yield raw, requirement


def test_requirements_extra_is_valid_pip_syntax():
    """Every non-comment line must parse as a PEP 508 requirement (what pip/Dependabot do)."""
    for raw, requirement in _requirement_lines():
        try:
            Requirement(requirement)
        except InvalidRequirement as err:  # pragma: no cover - failure path is the point
            pytest.fail(f'{requirement!r} (from line {raw!r}) is not valid pip syntax: {err}')


def test_get_extra_requires_maps_features_to_packages():
    """Each feature tag resolves to its intended package(s)."""
    get_extra_requires = _load_get_extra_requires()
    extras = {tag: {Requirement(dep).name for dep in deps}
              for tag, deps in get_extra_requires(str(EXTRAS_FILE)).items()}

    for tag, expected_packages in EXPECTED_FEATURES.items():
        assert extras.get(tag) == expected_packages, f'feature {tag!r} maps to the wrong package(s)'


def test_each_package_is_also_its_own_extra():
    """Back-compat: `pip install RNAlysis[<package>]` must keep working for every optional dep."""
    get_extra_requires = _load_get_extra_requires()
    extras = get_extra_requires(str(EXTRAS_FILE))
    for _, requirement in _requirement_lines():
        name = Requirement(requirement).name
        assert name in extras, f'package {name!r} is not exposed as its own extra'


def test_all_extra_collects_every_optional_dependency():
    get_extra_requires = _load_get_extra_requires()
    extras = get_extra_requires(str(EXTRAS_FILE))
    all_packages = {Requirement(dep).name for dep in extras['all']}
    file_packages = {Requirement(req).name for _, req in _requirement_lines()}
    assert all_packages == file_packages
