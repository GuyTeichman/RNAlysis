"""Unit tests for packaging/generate_api_vocabularies.py (everything except the live fetches).

The script lives in ``packaging/`` (which is *not* a package and would shadow the PyPI ``packaging``
distribution if imported by name), so it is loaded directly from its file path -- same trick as
tests/test_contributors.py.

The script itself is only ever run by hand at release time, but the snapshot it writes ships in
every install and is the *only* source of the GUI's taxon/gene-ID-type dropdowns (see
rnalysis/utils/param_typing.py), so its assembly logic is worth pinning: metadata that marks the
file auto-generated, fetch order preserved (UniProt's gene-ID-type order is meaningful), retries,
and the carry-over that keeps a temporarily-dead service from silently emptying a dropdown.
"""
import importlib.util
import json
from datetime import datetime
from pathlib import Path

import pytest

_MOD_PATH = Path(__file__).resolve().parent.parent / 'packaging' / 'generate_api_vocabularies.py'
_spec = importlib.util.spec_from_file_location('rnalysis_generate_api_vocabularies_script', _MOD_PATH)
generate_api_vocabularies = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(generate_api_vocabularies)

gav = generate_api_vocabularies


def _spec_for(fetch, key='gene_id_types', source='Example service (https://example.org)'):
    return {key: gav.VocabularySpec(source=source, fetch=fetch)}


def test_build_snapshot_marks_the_file_auto_generated():
    payload = gav.build_snapshot({}, version='9.9.9', generated_at='2026-01-02T03:04:05+00:00')

    assert 'generate_api_vocabularies.py' in payload['_comment']
    assert 'edit' in payload['_comment'].lower()
    assert payload['rnalysis_version'] == '9.9.9'
    # ISO-8601, parseable by anything that reads this file later
    assert datetime.fromisoformat(payload['generated_at'])


def test_build_snapshot_records_source_and_values_per_vocabulary():
    entries = {'panther_taxons': gav.VocabularyEntry(source='PantherDB (https://example.org)',
                                                     values=('Homo sapiens', 'Mus musculus'),
                                                     generated_at='2026-01-02T03:04:05+00:00')}
    payload = gav.build_snapshot(entries, version='9.9.9', generated_at='2026-01-02T03:04:05+00:00')

    vocab = payload['vocabularies']['panther_taxons']
    assert vocab['source'] == 'PantherDB (https://example.org)'
    assert vocab['values'] == ['Homo sapiens', 'Mus musculus']
    assert vocab['generated_at'] == '2026-01-02T03:04:05+00:00'


def test_build_snapshot_preserves_fetch_order():
    # UniProt returns gene-ID types in a deliberate priority order (io.get_legal_gene_id_types'
    # GROUP_PRIORITIZATION) that drives the order of the GUI dropdown -- never re-sort it.
    values = ('UniProtKB AC/ID', 'UniParc', 'Ensembl', 'ArrayExpress')
    entries = {'gene_id_types': gav.VocabularyEntry(source='UniProtKB', values=values,
                                                    generated_at='2026-01-02T03:04:05+00:00')}
    payload = gav.build_snapshot(entries, version='9.9.9', generated_at='2026-01-02T03:04:05+00:00')

    assert payload['vocabularies']['gene_id_types']['values'] == list(values)


def test_collect_vocabularies_fetches_each_vocabulary():
    entries, failures = gav.collect_vocabularies(_spec_for(lambda: ('a', 'b')), retries=1)

    assert failures == []
    assert entries['gene_id_types'].values == ('a', 'b')
    assert entries['gene_id_types'].source == 'Example service (https://example.org)'


def test_collect_vocabularies_retries_a_failing_fetch():
    calls = []

    def flaky():
        calls.append(1)
        if len(calls) < 3:
            raise ConnectionError('service hiccup')
        return ('a',)

    entries, failures = gav.collect_vocabularies(_spec_for(flaky), retries=3, retry_delay=0)

    assert len(calls) == 3
    assert failures == []
    assert entries['gene_id_types'].values == ('a',)


def test_collect_vocabularies_reports_a_service_that_stays_down():
    def always_fails():
        raise ConnectionError('service is down')

    entries, failures = gav.collect_vocabularies(_spec_for(always_fails), retries=2, retry_delay=0)

    assert entries == {}
    assert len(failures) == 1
    assert failures[0].key == 'gene_id_types'
    assert 'service is down' in failures[0].reason


def test_collect_vocabularies_keeps_going_after_one_service_fails():
    specs = {'gene_id_types': gav.VocabularySpec(source='A', fetch=lambda: ('a',)),
             'panther_taxons': gav.VocabularySpec(source='B', fetch=_raise_connection_error),
             'ensembl_taxons': gav.VocabularySpec(source='C', fetch=lambda: ('c',))}

    entries, failures = gav.collect_vocabularies(specs, retries=1, retry_delay=0)

    assert set(entries) == {'gene_id_types', 'ensembl_taxons'}
    assert [failure.key for failure in failures] == ['panther_taxons']


def _raise_connection_error():
    raise ConnectionError('service is down')


def test_collect_vocabularies_warns_about_an_empty_but_successful_fetch(capsys):
    # A service that answers politely with nothing is not an exception, but it would ship an empty
    # dropdown -- say so at release time instead of burying it in the JSON.
    entries, failures = gav.collect_vocabularies(_spec_for(tuple), retries=1, retry_delay=0)

    assert failures == []
    assert entries['gene_id_types'].values == ()
    assert 'gene_id_types' in capsys.readouterr().err


def test_carry_over_previous_keeps_values_of_a_dead_service():
    # A service being down at release time must not silently ship an empty dropdown -- keep the
    # previous release's values and mark them stale, so the maintainer sees what happened.
    previous = {'vocabularies': {'panther_taxons': {'source': 'PantherDB', 'values': ['Homo sapiens'],
                                                    'generated_at': '2025-01-01T00:00:00+00:00'}}}
    entries = {'gene_id_types': gav.VocabularyEntry(source='UniProtKB', values=('a',),
                                                    generated_at='2026-01-02T03:04:05+00:00')}

    merged, carried = gav.carry_over_previous(entries, previous, ['panther_taxons'])

    assert carried == ['panther_taxons']
    assert merged['panther_taxons'].values == ('Homo sapiens',)
    assert merged['panther_taxons'].generated_at == '2025-01-01T00:00:00+00:00'
    assert merged['panther_taxons'].stale is True
    assert merged['gene_id_types'].stale is False


def test_carry_over_previous_ignores_a_vocabulary_missing_from_the_previous_file():
    merged, carried = gav.carry_over_previous({}, {'vocabularies': {}}, ['panther_taxons'])

    assert merged == {}
    assert carried == []


def test_build_snapshot_flags_a_stale_entry():
    entries = {'panther_taxons': gav.VocabularyEntry(source='PantherDB', values=('Homo sapiens',),
                                                     generated_at='2025-01-01T00:00:00+00:00', stale=True)}
    payload = gav.build_snapshot(entries, version='9.9.9', generated_at='2026-01-02T03:04:05+00:00')

    assert payload['vocabularies']['panther_taxons']['stale'] is True


def test_write_snapshot_round_trips_non_ascii_values(tmp_path):
    # Species names carry non-ASCII characters; the file is written and read as UTF-8 everywhere
    # (Windows' default cp1252 would mangle them).
    entries = {'phylomedb_taxons': gav.VocabularyEntry(source='PhylomeDB', values=('Ceratitis capitata "Bréal"',),
                                                       generated_at='2026-01-02T03:04:05+00:00')}
    payload = gav.build_snapshot(entries, version='9.9.9', generated_at='2026-01-02T03:04:05+00:00')
    out_path = tmp_path / 'api_vocabularies.json'

    gav.write_snapshot(payload, out_path)

    assert json.loads(out_path.read_text(encoding='utf-8')) == payload
    assert out_path.read_text(encoding='utf-8').endswith('\n')


def test_read_previous_snapshot_returns_empty_for_a_missing_or_corrupt_file(tmp_path):
    assert gav.read_previous_snapshot(tmp_path / 'does_not_exist.json') == {}

    corrupt = tmp_path / 'corrupt.json'
    corrupt.write_text('{not json', encoding='utf-8')
    assert gav.read_previous_snapshot(corrupt) == {}


def test_default_specs_cover_every_getter_in_param_typing():
    # The snapshot is the only source for these four getters; a vocabulary missing here would
    # ship as an empty dropdown.
    assert set(gav.default_specs()) == {'gene_id_types', 'panther_taxons', 'phylomedb_taxons', 'ensembl_taxons'}


def test_main_writes_a_snapshot_and_reports_failures(tmp_path, monkeypatch, capsys):
    out_path = tmp_path / 'api_vocabularies.json'
    specs = {'gene_id_types': gav.VocabularySpec(source='A', fetch=lambda: ('a',)),
             'panther_taxons': gav.VocabularySpec(source='B', fetch=_raise_connection_error)}
    monkeypatch.setattr(gav, 'default_specs', lambda: specs)

    exit_code = gav.main(['--out', str(out_path), '--retries', '1', '--retry-delay', '0'])

    payload = json.loads(out_path.read_text(encoding='utf-8'))
    assert payload['vocabularies']['gene_id_types']['values'] == ['a']
    assert 'panther_taxons' not in payload['vocabularies']
    assert exit_code == 1  # one service stayed down
    assert 'panther_taxons' in capsys.readouterr().err


@pytest.mark.parametrize('key', ['gene_id_types', 'panther_taxons', 'phylomedb_taxons', 'ensembl_taxons'])
def test_default_specs_name_their_source(key):
    spec = gav.default_specs()[key]
    assert spec.source
    assert callable(spec.fetch)
