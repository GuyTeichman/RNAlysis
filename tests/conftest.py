import os
import sys

import pytest

# Force headless offscreen rendering for Qt on macOS CI environments
if sys.platform == 'darwin':
    os.environ["QT_QPA_PLATFORM"] = "offscreen"


def pytest_configure(config):
    os.environ["QT_DEBUG_PLUGINS"] = "1"


# --- test tiers (unit / integration / e2e) ---------------------------------------------------
# Modules whose tests are dominated by external services (R, external CLI tools, or web APIs that
# aren't mocked). Tests in these modules default to the `integration` tier.
_INTEGRATION_MODULES = frozenset({
    'test_io',                       # web APIs (UniProt/Ensembl/PANTHER/PhylomeDB/OrthoInspector/KEGG/GO)
    'test_fastq',                    # kallisto / bowtie2 / cutadapt CLIs
    'test_differential_expression',  # R (DESeq2 / limma-voom)
    'test_feature_counting',         # R (Rsubread featureCounts)
    'test_installs',                 # installs R packages
})
_TIER_MARKERS = ('unit', 'integration', 'e2e')


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
    """Auto-assign every test to exactly one tier (unit / integration / e2e).

    An explicit @pytest.mark.{unit,integration,e2e} on the test/class/module always wins; this only
    fills in a tier for tests that don't declare one. See pytest.ini for the tier definitions.
    """
    for item in items:
        if any(item.get_closest_marker(m) for m in _TIER_MARKERS):
            continue  # respect an explicit tier

        module = item.module.__name__.rsplit('.', 1)[-1] if item.module is not None else ''
        if module.startswith('test_gui') or 'qtbot' in getattr(item, 'fixturenames', ()):
            tier = 'e2e'
        elif module in _INTEGRATION_MODULES or _has_availability_skip(item):
            tier = 'integration'
        else:
            tier = 'unit'
        item.add_marker(getattr(pytest.mark, tier))
