import pytest
import requests

from rnalysis.utils import param_typing
from rnalysis.utils.param_typing import DEFAULT_ORGANISMS


def test_default_organisms_are_spelled_correctly():
    # These names are used as live taxonomy-search queries (io.map_taxon_id); a misspelling degrades
    # the match. Guard the historically-misspelled entry in particular.
    assert 'Arabidopsis thaliana' in DEFAULT_ORGANISMS
    assert 'Arabodopsis thaliana' not in DEFAULT_ORGANISMS


# The get_*_taxons / get_gene_id_types helpers run at import time (they populate Literal[...] type
# annotations on Filter/FeatureSet methods), so any uncaught error crashes `import rnalysis.filtering`
# for every user and every test collection. They must degrade to an empty tuple + a warning instead.
# See the flaky CI collection error where a PantherDB 200-with-garbage body raised JSONDecodeError.
@pytest.mark.parametrize('getter_name,io_func_name', [
    ('get_gene_id_types', 'get_legal_gene_id_types'),
    ('get_panther_taxons', 'get_legal_panther_taxons'),
    ('get_ensembl_taxons', 'get_legal_ensembl_taxons'),
    ('get_phylomedb_taxons', 'get_legal_phylomedb_taxons'),
])
@pytest.mark.parametrize('error', [
    ValueError('Expecting value: line 1 column 1 (char 0)'),  # json.loads on an empty/garbage 200 body
    KeyError('search'),                                       # 200 but an unexpected JSON shape
    requests.exceptions.ConnectionError('no network'),
    requests.exceptions.HTTPError('500 server error'),
    requests.exceptions.ReadTimeout('timed out'),            # the request-timeout path added here
])
def test_legal_value_getter_degrades_on_bad_response(monkeypatch, getter_name, io_func_name, error):
    getter = getattr(param_typing, getter_name)

    def _raise(*args, **kwargs):
        raise error

    monkeypatch.setattr(param_typing.io, io_func_name, _raise)
    getter.cache_clear()  # these are lru_cached; force re-evaluation through the mock
    try:
        with pytest.warns(UserWarning):
            result = getter()
        assert result == tuple()
    finally:
        getter.cache_clear()
