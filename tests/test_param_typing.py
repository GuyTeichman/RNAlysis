import json
import os
import subprocess
import sys
from pathlib import Path

import pytest

from rnalysis.utils import param_typing
from rnalysis.utils.param_typing import DEFAULT_ORGANISMS

REPO_ROOT = Path(__file__).resolve().parent.parent

# param_typing getter -> the vocabulary key it reads out of the packaged snapshot
# (rnalysis/data_files/api_vocabularies.json, written by packaging/generate_api_vocabularies.py).
GETTERS = {
    'get_gene_id_types': 'gene_id_types',
    'get_panther_taxons': 'panther_taxons',
    'get_phylomedb_taxons': 'phylomedb_taxons',
    'get_ensembl_taxons': 'ensembl_taxons',
}


def clear_vocabulary_caches():
    """Drop every lru_cache in the snapshot-reading path, so the next call re-reads the file."""
    param_typing._load_api_vocabularies.cache_clear()
    for getter_name in GETTERS:
        getattr(param_typing, getter_name).cache_clear()


@pytest.fixture
def snapshot_path(tmp_path, monkeypatch):
    """Point param_typing at a throwaway snapshot file, and leave its caches clean either way."""
    path = tmp_path / 'api_vocabularies.json'
    monkeypatch.setattr(param_typing, 'API_VOCABULARIES_PATH', path)
    clear_vocabulary_caches()
    yield path
    clear_vocabulary_caches()


def write_snapshot(path, vocabularies):
    payload = {'_comment': 'test fixture', 'generated_at': '2026-01-02T03:04:05+00:00',
               'rnalysis_version': '9.9.9',
               'vocabularies': {key: {'source': 'test', 'generated_at': '2026-01-02T03:04:05+00:00',
                                      'values': list(values)} for key, values in vocabularies.items()}}
    path.write_text(json.dumps(payload), encoding='utf-8')


def test_default_organisms_are_spelled_correctly():
    # These names are used as live taxonomy-search queries (io.map_taxon_id); a misspelling degrades
    # the match. Guard the historically-misspelled entry in particular.
    assert 'Arabidopsis thaliana' in DEFAULT_ORGANISMS
    assert 'Arabodopsis thaliana' not in DEFAULT_ORGANISMS


@pytest.mark.parametrize('getter_name,key', list(GETTERS.items()))
def test_getter_returns_the_snapshot_values(snapshot_path, getter_name, key):
    write_snapshot(snapshot_path, {key: ['Zebra finch', 'Homo sapiens']})

    assert getattr(param_typing, getter_name)() == ('Zebra finch', 'Homo sapiens')


@pytest.mark.parametrize('getter_name,key', list(GETTERS.items()))
def test_getter_preserves_snapshot_order(snapshot_path, getter_name, key):
    # UniProt returns gene-ID types in a deliberate priority order that becomes the order of the
    # GUI dropdown; the snapshot preserves it and so must the getter.
    values = ['UniProtKB AC/ID', 'UniParc', 'Ensembl', 'ArrayExpress']
    write_snapshot(snapshot_path, {key: values})

    assert getattr(param_typing, getter_name)() == tuple(values)


# The getters populate Literal[...] type annotations on Filter/FeatureSet methods, so they run while
# `import rnalysis.filtering` executes the class bodies: an uncaught error here crashes the import for
# every user and every test-collection run. They must degrade to an empty tuple + a warning instead.
# Since the data now comes from a packaged file rather than the network, every one of these failures
# means a broken installation, and the warning must say so (telling the user to check their internet
# connection would be worse than useless).
@pytest.mark.parametrize('getter_name,key', list(GETTERS.items()))
def test_getter_degrades_when_the_snapshot_is_missing(snapshot_path, getter_name, key):
    assert not snapshot_path.exists()

    with pytest.warns(UserWarning, match='install'):
        assert getattr(param_typing, getter_name)() == tuple()


@pytest.mark.parametrize('getter_name,key', list(GETTERS.items()))
def test_getter_degrades_on_a_corrupt_snapshot(snapshot_path, getter_name, key):
    snapshot_path.write_text('{"vocabularies": {not json', encoding='utf-8')

    with pytest.warns(UserWarning, match='install'):
        assert getattr(param_typing, getter_name)() == tuple()


@pytest.mark.parametrize('getter_name,key', list(GETTERS.items()))
@pytest.mark.parametrize('payload', [
    {},                                                        # no 'vocabularies' at all
    {'vocabularies': {}},                                      # the vocabulary is missing
    {'vocabularies': {'KEY': {'source': 'test'}}},             # the entry has no 'values'
    {'vocabularies': {'KEY': {'values': 'Homo sapiens'}}},     # a string instead of a list
    {'vocabularies': {'KEY': {'values': [{'name': 'Homo sapiens'}]}}},  # not a list of strings
    {'vocabularies': 'nope'},
])
def test_getter_degrades_on_an_unexpected_snapshot_shape(snapshot_path, getter_name, key, payload):
    # A Literal[...] annotation can only be built from hashable values, so a wrong-shaped snapshot
    # must not reach the annotation - it would crash the import instead of emptying one dropdown.
    payload = json.loads(json.dumps(payload).replace('"KEY"', f'"{key}"'))
    snapshot_path.write_text(json.dumps(payload), encoding='utf-8')

    with pytest.warns(UserWarning, match='install'):
        assert getattr(param_typing, getter_name)() == tuple()


@pytest.mark.parametrize('getter_name,key', list(GETTERS.items()))
def test_getter_does_not_touch_the_network(monkeypatch, getter_name, key):
    """The whole point of the snapshot: no live service is contacted to build a type annotation."""
    from rnalysis.utils import io

    def _fail(*args, **kwargs):
        raise AssertionError('a param_typing getter performed network I/O')

    for io_func_name in ('get_legal_gene_id_types', 'get_legal_panther_taxons',
                         'get_legal_phylomedb_taxons', 'get_legal_ensembl_taxons'):
        monkeypatch.setattr(io, io_func_name, _fail)
    clear_vocabulary_caches()
    try:
        values = getattr(param_typing, getter_name)()
        assert values == tuple(load_packaged_snapshot()['vocabularies'][key]['values'])
    finally:
        clear_vocabulary_caches()


def load_packaged_snapshot():
    return json.loads(param_typing.API_VOCABULARIES_PATH.read_text(encoding='utf-8'))


def test_packaged_snapshot_is_shipped_and_well_formed():
    assert param_typing.API_VOCABULARIES_PATH.is_file(), \
        'the packaged vocabulary snapshot is missing - regenerate it with ' \
        'packaging/generate_api_vocabularies.py'
    payload = load_packaged_snapshot()

    assert 'generate_api_vocabularies.py' in payload['_comment']
    assert payload['generated_at']
    assert payload['rnalysis_version']
    assert set(payload['vocabularies']) == set(GETTERS.values())
    for key, vocabulary in payload['vocabularies'].items():
        assert vocabulary['source'], f'{key} does not name its source service'
        assert all(isinstance(value, str) for value in vocabulary['values'])


@pytest.mark.parametrize('key', list(GETTERS.values()))
def test_packaged_snapshot_is_populated(key):
    # A regeneration that ran while a service was down (or against a parser that no longer matches
    # the service's format, as in issue #263) would empty a GUI dropdown for a whole release - fail
    # loudly here rather than shipping it.
    assert len(load_packaged_snapshot()['vocabularies'][key]['values']) > 10


IMPORT_WITHOUT_NETWORK_SCRIPT = '''
import socket
import sys
import time


class NetworkAccessDuringImport(RuntimeError):
    pass


def blocked(*args, **kwargs):
    raise NetworkAccessDuringImport('rnalysis.filtering tried to use the network at import time')


# Block the points where a socket actually reaches the network, rather than replacing socket.socket
# itself - the class is subclassed by the standard library (ssl.SSLSocket), so it has to stay a class.
socket.socket.connect = blocked
socket.socket.connect_ex = blocked
socket.create_connection = blocked
socket.getaddrinfo = blocked

start = time.perf_counter()
from rnalysis import filtering  # noqa: F401
from rnalysis.utils import param_typing

elapsed = time.perf_counter() - start
counts = {name: len(getattr(param_typing, name)())
          for name in ('get_gene_id_types', 'get_panther_taxons', 'get_ensembl_taxons')}
print('IMPORT_OK', elapsed, counts, file=sys.stderr)
'''


def test_import_filtering_works_without_network(tmp_path):
    """`import rnalysis.filtering` must not need a socket - it is pure local work now.

    Run in a subprocess because the check is about *import* time, and rnalysis is already imported
    in this process. Blocking is done inside the child (no firewall needed, nothing shared with the
    parent), so this is safe to run in CI.
    """
    script = tmp_path / 'import_without_network.py'
    script.write_text(IMPORT_WITHOUT_NETWORK_SCRIPT, encoding='utf-8')

    result = subprocess.run([sys.executable, str(script)], capture_output=True, text=True, timeout=300,
                            cwd=str(REPO_ROOT), env={**os.environ, 'PYTHONPATH': str(REPO_ROOT)})

    assert result.returncode == 0, f'importing rnalysis.filtering without network failed:\n{result.stderr}'
    assert 'IMPORT_OK' in result.stderr
    # the vocabularies came from the packaged snapshot, not from a degraded empty fallback
    assert "'get_panther_taxons': 0" not in result.stderr
