"""Unit tests for the live-service availability probes in ``tests/__init__.py``.

These probes gate the ``@pytest.mark.skipif(not <SERVICE>_AVAILABLE, ...)`` network tests. A probe
that reports a *healthy* service as unavailable silently skips real coverage on every CI run, so the
probes themselves are worth pinning down. All tests here mock the HTTP layer -- no real network -- so
they belong to the ``unit`` tier (the module name keeps them out of ``integration_net``).
"""
import pytest
import requests

import tests as tests_pkg


# A healthy PantherDB answers an ortholog/matchortho query with a non-empty `search.mapping.mapped`.
_MAPPED_PAYLOAD = {'search': {'mapping': {'mapped': [{'target_gene': 'HUMAN|HGNC=1|UniProtKB=P12345'}]}}}


class _FakeResp:
    def __init__(self, status_code, json_payload=None, raise_json=False):
        self.status_code = status_code
        self._json_payload = json_payload
        self._raise_json = raise_json

    def json(self):
        if self._raise_json:
            raise requests.exceptions.JSONDecodeError('Expecting value', '', 0)
        return self._json_payload


@pytest.mark.unit
class TestIsPantherdbAvailable:
    """`is_pantherdb_available` must report a healthy PantherDB as up -- and, crucially, a PantherDB
    whose ortholog endpoint is degraded as *down*, so the live ``TestPantherOrthologMapper`` suite
    skips instead of failing.

    Two regressions are guarded here. (1) The original probe GET-requested the bare ``.../ortholog/``
    parent, which 404s even on a healthy service, so it returned False on *every* run and the suite
    was permanently skipped. (2) Its replacement POSTed the param-free ``supportedgenomes`` endpoint,
    which PantherDB can answer 200 while its ``ortholog/matchortho`` endpoint returns empty/invalid
    200 bodies -- so the suite ran and its ``assert len(...) > 0`` checks failed instead of skipping.
    The probe now exercises the exact ``ortholog/matchortho`` capability the tests need and treats an
    empty/malformed mapping as unavailable.
    """

    @staticmethod
    def _healthy_pantherdb(url, *args, **kwargs):
        # A healthy service answers the ortholog/matchortho query with a non-empty mapping.
        if 'matchortho' in url:
            return _FakeResp(200, json_payload=_MAPPED_PAYLOAD)
        return _FakeResp(404)

    def test_reports_up_when_service_healthy(self, monkeypatch):
        monkeypatch.setattr(tests_pkg.requests, 'get', self._healthy_pantherdb)
        monkeypatch.setattr(tests_pkg.requests, 'post', self._healthy_pantherdb)
        assert tests_pkg.is_pantherdb_available() is True

    def test_probes_ortholog_matchortho_endpoint(self, monkeypatch):
        calls = []

        def record(url, *args, **kwargs):
            calls.append(('post', url, kwargs.get('timeout')))
            return _FakeResp(200, json_payload=_MAPPED_PAYLOAD)

        def forbid_get(url, *args, **kwargs):
            raise AssertionError(f'probe must POST the ortholog endpoint, not GET (got {url})')

        monkeypatch.setattr(tests_pkg.requests, 'post', record)
        monkeypatch.setattr(tests_pkg.requests, 'get', forbid_get)

        assert tests_pkg.is_pantherdb_available() is True
        assert len(calls) == 1
        _, url, timeout = calls[0]
        assert 'ortholog/matchortho' in url
        assert timeout is not None, 'probe must set a request timeout so a hung service cannot block collection'

    def test_reports_down_on_empty_mapping(self, monkeypatch):
        # HTTP 200 but no orthologs mapped -- a degraded ortholog service the suite must skip, not fail.
        payload = {'search': {'mapping': {'mapped': []}}}
        monkeypatch.setattr(tests_pkg.requests, 'post', lambda *a, **k: _FakeResp(200, json_payload=payload))
        assert tests_pkg.is_pantherdb_available() is False

    def test_reports_down_on_malformed_payload(self, monkeypatch):
        monkeypatch.setattr(tests_pkg.requests, 'post', lambda *a, **k: _FakeResp(200, json_payload={'search': {}}))
        assert tests_pkg.is_pantherdb_available() is False

    def test_reports_down_on_empty_body(self, monkeypatch):
        monkeypatch.setattr(tests_pkg.requests, 'post', lambda *a, **k: _FakeResp(200, raise_json=True))
        assert tests_pkg.is_pantherdb_available() is False

    def test_reports_down_on_connection_error(self, monkeypatch):
        def boom(*args, **kwargs):
            raise requests.exceptions.ConnectionError('service unreachable')

        monkeypatch.setattr(tests_pkg.requests, 'post', boom)
        monkeypatch.setattr(tests_pkg.requests, 'get', boom)
        assert tests_pkg.is_pantherdb_available() is False

    @pytest.mark.parametrize('status_code', [404, 500])
    def test_reports_down_on_4xx_5xx(self, monkeypatch, status_code):
        monkeypatch.setattr(tests_pkg.requests, 'post', lambda *a, **k: _FakeResp(status_code))
        monkeypatch.setattr(tests_pkg.requests, 'get', lambda *a, **k: _FakeResp(status_code))
        assert tests_pkg.is_pantherdb_available() is False
