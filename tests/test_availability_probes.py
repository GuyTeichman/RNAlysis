"""Unit tests for the live-service availability probes in ``tests/__init__.py``.

These probes gate the ``@pytest.mark.skipif(not <SERVICE>_AVAILABLE, ...)`` network tests. A probe
that reports a *healthy* service as unavailable silently skips real coverage on every CI run, so the
probes themselves are worth pinning down. All tests here mock the HTTP layer -- no real network -- so
they belong to the ``unit`` tier (the module name keeps them out of ``integration_net``).
"""
import pytest
import requests

import tests as tests_pkg


class _FakeResp:
    def __init__(self, status_code):
        self.status_code = status_code


@pytest.mark.unit
class TestIsPantherdbAvailable:
    """`is_pantherdb_available` must report a healthy PantherDB as up.

    Regression guard for the probe that GET-requested the bare ``.../ortholog/`` parent endpoint,
    which returns HTTP 404 even when PantherDB is perfectly healthy -- so the probe returned False
    on *every* run and the whole ``TestPantherOrthologMapper`` suite was permanently skipped.
    """

    @staticmethod
    def _healthy_pantherdb(url, *args, **kwargs):
        # Mirror the real service: the param-free supportedgenomes endpoint answers 200 with JSON,
        # while the bare .../ortholog/ parent (no query params) 404s even when the service is up.
        if 'supportedgenomes' in url:
            return _FakeResp(200)
        return _FakeResp(404)

    def test_reports_up_when_service_healthy(self, monkeypatch):
        monkeypatch.setattr(tests_pkg.requests, 'get', self._healthy_pantherdb)
        monkeypatch.setattr(tests_pkg.requests, 'post', self._healthy_pantherdb)
        assert tests_pkg.is_pantherdb_available() is True

    def test_probes_supportedgenomes_endpoint(self, monkeypatch):
        calls = []

        def record(url, *args, **kwargs):
            calls.append(('post', url, kwargs.get('timeout')))
            return _FakeResp(200)

        def forbid_get(url, *args, **kwargs):
            raise AssertionError(f'probe must not GET the 404-ing ortholog parent (got {url})')

        monkeypatch.setattr(tests_pkg.requests, 'post', record)
        monkeypatch.setattr(tests_pkg.requests, 'get', forbid_get)

        assert tests_pkg.is_pantherdb_available() is True
        assert len(calls) == 1
        _, url, timeout = calls[0]
        assert 'supportedgenomes' in url
        assert timeout is not None, 'probe must set a request timeout so a hung service cannot block collection'

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
