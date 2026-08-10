"""Unit tests for packaging/snapshot_api_payload.py (pure logic; no live network).

The script lives in ``packaging/`` (which is *not* a Python package and would shadow the PyPI
``packaging`` distribution if imported by name -- see ``tests/test_contributors.py`` for the same
pattern), so it is loaded directly from its file path.

Every test here mocks the HTTP layer (``requests_mock``, the convention already used in
``tests/test_io.py``) or exercises pure functions directly. Nothing in this module makes a real
network call, so it runs in the fast ``unit`` tier (auto-assigned by tests/conftest.py, since this
module is not one of the modules bound to a network/CLI/GUI tier).
"""
import importlib.util
from pathlib import Path

import pytest
import requests
import requests_mock

_MOD_PATH = Path(__file__).resolve().parent.parent / 'packaging' / 'snapshot_api_payload.py'
_spec = importlib.util.spec_from_file_location('rnalysis_snapshot_api_payload', _MOD_PATH)
snap = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(snap)


# ---------------------------------------------------------------------------
# parse_key_value_pairs (--param)
# ---------------------------------------------------------------------------

def test_parse_key_value_pairs_basic():
    assert snap.parse_key_value_pairs(['a=1', 'b=2'], '--param') == {'a': '1', 'b': '2'}


def test_parse_key_value_pairs_empty_list():
    assert snap.parse_key_value_pairs([], '--param') == {}


def test_parse_key_value_pairs_value_may_contain_equals():
    assert snap.parse_key_value_pairs(['filter=a=b'], '--param') == {'filter': 'a=b'}


def test_parse_key_value_pairs_rejects_missing_equals():
    with pytest.raises(ValueError, match='--param'):
        snap.parse_key_value_pairs(['not-a-pair'], '--param')


def test_parse_key_value_pairs_rejects_empty_key():
    with pytest.raises(ValueError):
        snap.parse_key_value_pairs(['=value'], '--param')


# ---------------------------------------------------------------------------
# parse_headers (--header)
# ---------------------------------------------------------------------------

def test_parse_headers_basic():
    assert snap.parse_headers(['Accept: application/json']) == {'Accept': 'application/json'}


def test_parse_headers_splits_on_first_colon_only():
    # header values may themselves contain a colon (e.g. a bearer token) -- must survive intact
    result = snap.parse_headers(['Authorization: Bearer abc:def:ghi'])
    assert result == {'Authorization': 'Bearer abc:def:ghi'}


def test_parse_headers_strips_whitespace():
    assert snap.parse_headers(['X-Test  :   value  ']) == {'X-Test': 'value'}


def test_parse_headers_rejects_missing_colon():
    with pytest.raises(ValueError, match='--header'):
        snap.parse_headers(['no-colon-here'])


def test_parse_headers_rejects_empty_key():
    with pytest.raises(ValueError):
        snap.parse_headers([': value'])


# ---------------------------------------------------------------------------
# guess_extension / derive_filename
# ---------------------------------------------------------------------------

@pytest.mark.parametrize('content_type,expected', [
    ('application/json', '.json'),
    ('application/json; charset=utf-8', '.json'),
    ('text/xml', '.xml'),
    ('application/xml', '.xml'),
    ('text/csv', '.csv'),
    ('text/html; charset=utf-8', '.html'),
    ('text/plain', '.txt'),
    (None, '.txt'),
    ('application/octet-stream', '.txt'),
])
def test_guess_extension(content_type, expected):
    assert snap.guess_extension(content_type) == expected


def test_derive_filename_sanitizes_url_into_safe_characters():
    name = snap.derive_filename('https://rest.uniprot.org/idmapping/status/123?fields=x')
    assert name == 'rest_uniprot_org_idmapping_status_123.txt'


def test_derive_filename_uses_content_type_extension():
    name = snap.derive_filename('https://api.example.org/search', content_type='application/json')
    assert name.endswith('.json')


def test_derive_filename_truncates_long_urls():
    long_path = 'a' * 300
    name = snap.derive_filename(f'https://example.org/{long_path}')
    stem = name.rsplit('.', 1)[0]
    assert len(stem) <= snap._MAX_STEM_LENGTH


def test_derive_filename_never_empty_stem():
    # a URL with nothing but punctuation in host+path should still produce a usable name
    name = snap.derive_filename('https://.../')
    assert name and name != '.txt'


# ---------------------------------------------------------------------------
# resolve_output_path
# ---------------------------------------------------------------------------

def test_resolve_output_path_defaults_to_derived_name_under_default_dir(tmp_path):
    default_dir = tmp_path / 'test_files'
    path = snap.resolve_output_path('https://rest.uniprot.org/idmapping/status/123', None,
                                     default_dir=default_dir)
    assert path == default_dir / 'rest_uniprot_org_idmapping_status_123.txt'


def test_resolve_output_path_bare_filename_lands_under_default_dir(tmp_path):
    default_dir = tmp_path / 'test_files'
    path = snap.resolve_output_path('https://example.org/x', 'my_fixture.json', default_dir=default_dir)
    assert path == default_dir / 'my_fixture.json'


def test_resolve_output_path_with_directory_is_used_as_is(tmp_path):
    default_dir = tmp_path / 'test_files'
    explicit = tmp_path / 'elsewhere' / 'fixture.json'
    path = snap.resolve_output_path('https://example.org/x', str(explicit), default_dir=default_dir)
    assert path == explicit


# ---------------------------------------------------------------------------
# write_payload
# ---------------------------------------------------------------------------

def test_write_payload_creates_parent_dirs_and_writes_bytes(tmp_path):
    out = tmp_path / 'nested' / 'dir' / 'payload.json'
    snap.write_payload(out, b'{"hello": "world"}')
    assert out.read_bytes() == b'{"hello": "world"}'


def test_write_payload_overwrites_existing_file(tmp_path):
    out = tmp_path / 'payload.txt'
    out.write_bytes(b'old content')
    snap.write_payload(out, b'new content')
    assert out.read_bytes() == b'new content'


# ---------------------------------------------------------------------------
# main() -- HTTP mocked via requests_mock, never a live call
# ---------------------------------------------------------------------------

def test_main_writes_fetched_payload_to_derived_path(tmp_path, monkeypatch):
    monkeypatch.setattr(snap, 'DEFAULT_OUTPUT_DIR', tmp_path)
    with requests_mock.Mocker() as m:
        m.get('https://rest.uniprot.org/idmapping/status/abc123', content=b'{"jobStatus": "FINISHED"}',
              headers={'Content-Type': 'application/json'})
        exit_code = snap.main(['https://rest.uniprot.org/idmapping/status/abc123'])

    assert exit_code == 0
    out_file = tmp_path / 'rest_uniprot_org_idmapping_status_abc123.json'
    assert out_file.read_bytes() == b'{"jobStatus": "FINISHED"}'


def test_main_respects_explicit_out_param_and_header(tmp_path, monkeypatch):
    monkeypatch.setattr(snap, 'DEFAULT_OUTPUT_DIR', tmp_path)
    out_path = tmp_path / 'custom' / 'fixture.txt'
    captured = {}

    with requests_mock.Mocker() as m:
        m.get('https://api.example.org/search', content=b'payload-bytes')
        exit_code = snap.main([
            'https://api.example.org/search',
            '--param', 'query=BRCA1',
            '--header', 'Accept: application/json',
            '--out', str(out_path),
        ])
        captured['qs'] = m.last_request.qs
        captured['headers'] = m.last_request.headers

    assert exit_code == 0
    assert out_path.read_bytes() == b'payload-bytes'
    assert captured['qs'] == {'query': ['brca1']}  # requests_mock lower-cases query values in qs
    assert captured['headers']['Accept'] == 'application/json'


def test_main_returns_nonzero_and_does_not_write_on_http_error(tmp_path, monkeypatch):
    monkeypatch.setattr(snap, 'DEFAULT_OUTPUT_DIR', tmp_path)
    with requests_mock.Mocker() as m:
        m.get('https://api.example.org/missing', status_code=404, content=b'not found')
        exit_code = snap.main(['https://api.example.org/missing'])

    assert exit_code != 0
    assert list(tmp_path.iterdir()) == []


def test_main_rejects_malformed_param(tmp_path, monkeypatch, capsys):
    monkeypatch.setattr(snap, 'DEFAULT_OUTPUT_DIR', tmp_path)
    with pytest.raises(SystemExit):
        snap.main(['https://api.example.org/x', '--param', 'not-a-pair'])


def test_build_arg_parser_requires_url():
    parser = snap.build_arg_parser()
    with pytest.raises(SystemExit):
        parser.parse_args([])


def test_build_arg_parser_parses_repeated_flags():
    parser = snap.build_arg_parser()
    args = parser.parse_args(['https://x.test/y', '--param', 'a=1', '--param', 'b=2',
                              '--header', 'Accept: */*'])
    assert args.url == 'https://x.test/y'
    assert args.param == ['a=1', 'b=2']
    assert args.header == ['Accept: */*']
