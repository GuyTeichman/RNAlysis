from rnalysis.utils.param_typing import DEFAULT_ORGANISMS


def test_default_organisms_are_spelled_correctly():
    # These names are used as live taxonomy-search queries (io.map_taxon_id); a misspelling degrades
    # the match. Guard the historically-misspelled entry in particular.
    assert 'Arabidopsis thaliana' in DEFAULT_ORGANISMS
    assert 'Arabodopsis thaliana' not in DEFAULT_ORGANISMS
