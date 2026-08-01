import os
import sys

import pytest

# Force headless offscreen rendering for Qt on macOS CI environments
if sys.platform == 'darwin':
    os.environ["QT_QPA_PLATFORM"] = "offscreen"


def pytest_configure(config):
    os.environ["QT_DEBUG_PLUGINS"] = "1"


# --- test tiers (unit / integration / e2e) ---------------------------------------------------
# The `integration` tier is split into two sub-tiers so CI can time (and, later, schedule) them
# separately:
#   * integration_net   — bound to external WEB services; platform-independent.
#   * integration_tools — bound to external CLIs / R; exercise OS-specific binaries.
# Both sub-tiers also carry the umbrella `integration` marker, so `-m integration` still selects all
# of them (back-compat with the old single tier).
_INTEGRATION_NET_MODULES = frozenset({
    'test_io',                       # web APIs (UniProt/Ensembl/PANTHER/PhylomeDB/OrthoInspector/KEGG/GO)
})
_INTEGRATION_TOOLS_MODULES = frozenset({
    'test_fastq',                    # kallisto / bowtie2 / cutadapt CLIs
    'test_differential_expression',  # R (DESeq2 / limma-voom)
    'test_feature_counting',         # R (Rsubread featureCounts)
    'test_installs',                 # installs R packages
})
# Leaf tiers — an explicit one of these on a test/class/module wins over auto-assignment. The
# umbrella `integration` is intentionally NOT listed here: a test tagged only `integration` still
# gets sub-classified into net/tools so it lands in exactly one CI step.
_TIER_MARKERS = ('unit', 'integration_net', 'integration_tools', 'e2e')


def _tier_markers_for(module: str, fixturenames, has_availability_skip: bool) -> tuple:
    """Return the marker name(s) a test should carry, from its module name and traits.

    Pure function (no pytest state) so it can be unit-tested directly. Integration tests get two
    markers: the umbrella `integration` plus their leaf sub-tier. Module membership wins over the
    availability-skip heuristic (a tools test with an 'available' skipif stays in tools).
    """
    if module.startswith('test_gui') or 'qtbot' in fixturenames:
        return ('e2e',)
    if module in _INTEGRATION_TOOLS_MODULES:
        return ('integration', 'integration_tools')
    if module in _INTEGRATION_NET_MODULES or has_availability_skip:
        return ('integration', 'integration_net')
    return ('unit',)


def _has_availability_skip(item) -> bool:
    """True if the test is gated by a service-availability skipif (its reason mentions 'available').

    This catches the live-network tests scattered through otherwise-mocked modules
    (e.g. the Ensembl/UniProt tests in test_enrichment.py) without per-test marking.
    """
    for marker in item.iter_markers(name='skipif'):
        if 'available' in str(marker.kwargs.get('reason', '')).lower():
            return True
    return False


def pytest_collection_modifyitems(config, items):
    """Auto-assign every test to a tier (unit / integration_net / integration_tools / e2e).

    An explicit leaf-tier marker (@pytest.mark.{unit,integration_net,integration_tools,e2e}) on the
    test/class/module always wins; this only fills in a tier for tests that don't declare one.
    Integration tests additionally carry the umbrella `integration` marker. See pytest.ini.
    """
    for item in items:
        if any(item.get_closest_marker(m) for m in _TIER_MARKERS):
            continue  # respect an explicit leaf tier

        module = item.module.__name__.rsplit('.', 1)[-1] if item.module is not None else ''
        markers = _tier_markers_for(module, getattr(item, 'fixturenames', ()), _has_availability_skip(item))
        for marker in markers:
            item.add_marker(getattr(pytest.mark, marker))
