"""Unit tests for the CI test-tier auto-classification in conftest.py (`_tier_markers_for`)."""

import pytest

from tests.conftest import _tier_markers_for


@pytest.mark.parametrize(
    'module,fixturenames,availability_skip,expected',
    [
        # GUI modules and any qtbot-using test -> e2e
        ('test_gui', (), False, ('e2e',)),
        ('test_gui_windows', (), False, ('e2e',)),
        ('test_filtering', ('qtbot',), False, ('e2e',)),
        # web-API module -> integration_net (with umbrella)
        ('test_io', (), False, ('integration', 'integration_net')),
        # a live-network test in an otherwise-mocked module (availability skipif) -> net
        ('test_enrichment', (), True, ('integration', 'integration_net')),
        # CLI / R modules -> integration_tools (with umbrella)
        ('test_fastq', (), False, ('integration', 'integration_tools')),
        ('test_differential_expression', (), False, ('integration', 'integration_tools')),
        ('test_feature_counting', (), False, ('integration', 'integration_tools')),
        ('test_installs', (), False, ('integration', 'integration_tools')),
        # module membership wins over the availability-skip heuristic (tools stays tools)
        ('test_fastq', (), True, ('integration', 'integration_tools')),
        # everything else -> unit
        ('test_filtering', (), False, ('unit',)),
        ('test_enrichment', (), False, ('unit',)),
        ('', (), False, ('unit',)),
    ],
)
def test_tier_markers_for(module, fixturenames, availability_skip, expected):
    assert _tier_markers_for(module, fixturenames, availability_skip) == expected
