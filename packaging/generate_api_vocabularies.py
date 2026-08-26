#!/usr/bin/env python
"""Regenerate the packaged snapshot of the remote vocabularies used in RNAlysis' type annotations.

Four vocabularies -- UniProtKB's gene-ID types, and PantherDB's / Ensembl's / PhylomeDB's legal
taxons -- appear inside ``Literal[...]`` type annotations on public ``Filter``/``FeatureSet``
methods, which is how the reflection-generated GUI knows what to put in those dropdowns (see
CLAUDE.md: "type annotations are load-bearing UI"). Annotations are evaluated when the class body
runs, so fetching those lists live would put four web requests on the *import* path of
``rnalysis.filtering`` -- slow, and empty dropdowns on an offline machine (a cluster, a locked-down
lab PC). Instead they are snapshotted here, at release time, into

    rnalysis/data_files/api_vocabularies.json

which ``rnalysis/utils/param_typing.py`` reads (no network, milliseconds). The legal values are
therefore pinned per RNAlysis version: identical on every machine running that version, refreshed
with each release. Users are never blocked by a slightly stale list -- every one of these
parameters is annotated ``Union[str, Literal[...]]``, which the GUI renders as a combo box with a
free-text "other..." field.

This is a *release-time* helper, like ``packaging/checksum_videos.py``: run it by hand from a
checkout, commit the regenerated JSON, and only then tag the release.

Usage
-----
    python packaging/generate_api_vocabularies.py                    # regenerate in place
    python packaging/generate_api_vocabularies.py --out /tmp/x.json  # write elsewhere (dry run)
    python packaging/generate_api_vocabularies.py --retries 5        # flakier network

Notes
-----
* This makes REAL, live requests to UniProtKB, PantherDB, Ensembl and PhylomeDB (FTP) -- it must
  be run from a machine with network access, and it reuses ``rnalysis.utils.io``'s own
  ``get_legal_*`` functions so the snapshot is exactly what the live path would have produced.
* Values keep the order the service returned them in. UniProt's gene-ID-type order in particular
  is deliberate (``io.get_legal_gene_id_types``' ``GROUP_PRIORITIZATION``) and shows up as the
  order of the GUI dropdown -- do not re-sort it.
* If a service stays down after the retries, the previous snapshot's values for it are carried
  over and marked ``"stale": true`` (a dead service must not silently empty a dropdown), the
  failure is reported on stderr, and the exit code is the number of failed services. Decide
  consciously whether to ship a stale entry or wait for the service to come back.
"""

import argparse
import dataclasses
import json
import sys
import time
from datetime import datetime, timezone
from pathlib import Path
from typing import Callable, Dict, List, Sequence, Tuple

REPO_ROOT = Path(__file__).resolve().parent.parent
# Append (never insert) the repo root: this script lives in packaging/, and putting the repo root
# *first* would let `import packaging` resolve to that directory instead of the PyPI distribution.
if str(REPO_ROOT) not in sys.path:
    sys.path.append(str(REPO_ROOT))

import rnalysis  # noqa: E402  (needs the sys.path bootstrap above when RNAlysis isn't pip-installed)
from rnalysis.utils import io, parsing  # noqa: E402

DEFAULT_OUTPUT_PATH = REPO_ROOT / 'rnalysis' / 'data_files' / 'api_vocabularies.json'
DEFAULT_RETRIES = 3
DEFAULT_RETRY_DELAY = 5.0

GENERATED_COMMENT = (
    'AUTO-GENERATED at release time by packaging/generate_api_vocabularies.py - do not edit by hand. '
    'This is a snapshot of the remote vocabularies (UniProtKB gene-ID types, PantherDB/Ensembl/PhylomeDB '
    'taxons) that RNAlysis bakes into Literal[...] type annotations, so that importing rnalysis performs '
    'no network I/O and the GUI dropdowns are populated offline. Regenerate it with that script, '
    'as part of cutting a release.'
)


@dataclasses.dataclass(frozen=True)
class VocabularySpec:
    """How to fetch one vocabulary: a human-readable source, and the (live, networked) fetcher."""

    source: str
    fetch: Callable[[], Sequence[str]]


@dataclasses.dataclass(frozen=True)
class VocabularyEntry:
    """One fetched (or carried-over) vocabulary, as it will be written to the snapshot."""

    source: str
    values: Tuple[str, ...]
    generated_at: str
    stale: bool = False


@dataclasses.dataclass(frozen=True)
class VocabularyFailure:
    """A vocabulary whose service stayed down for every attempt."""

    key: str
    reason: str


def fetch_gene_id_types() -> Tuple[str, ...]:
    """UniProtKB's ID-mapping fields, in UniProt's priority order (see io.get_legal_gene_id_types)."""
    return parsing.data_to_tuple(io.get_legal_gene_id_types()[0].keys())


def fetch_panther_taxons() -> Tuple[str, ...]:
    return parsing.data_to_tuple(io.get_legal_panther_taxons())


def fetch_phylomedb_taxons() -> Tuple[str, ...]:
    return parsing.data_to_tuple(io.get_legal_phylomedb_taxons())


def fetch_ensembl_taxons() -> Tuple[str, ...]:
    return parsing.data_to_tuple(io.get_legal_ensembl_taxons())


def default_specs() -> Dict[str, VocabularySpec]:
    """The four vocabularies ``rnalysis.utils.param_typing`` reads out of the snapshot.

    The keys are the snapshot's vocabulary names and must stay in sync with param_typing's getters
    (``get_gene_id_types``, ``get_panther_taxons``, ``get_phylomedb_taxons``, ``get_ensembl_taxons``).
    """
    return {
        'gene_id_types': VocabularySpec(
            source='UniProtKB ID-mapping fields (https://rest.uniprot.org/configure/idmapping/fields)',
            fetch=fetch_gene_id_types,
        ),
        'panther_taxons': VocabularySpec(
            source='PantherDB supported genomes (https://www.pantherdb.org/services/oai/pantherdb/supportedgenomes)',
            fetch=fetch_panther_taxons,
        ),
        'phylomedb_taxons': VocabularySpec(
            source='PhylomeDB MetaPhOrs species list (ftp://phylomedb.org/metaphors/latest/species.txt.gz)',
            fetch=fetch_phylomedb_taxons,
        ),
        'ensembl_taxons': VocabularySpec(
            source='Ensembl REST info/species (https://rest.ensembl.org/info/species)', fetch=fetch_ensembl_taxons
        ),
    }


def utcnow_iso() -> str:
    return datetime.now(timezone.utc).isoformat(timespec='seconds')


def collect_vocabularies(
    specs: Dict[str, VocabularySpec],
    retries: int = DEFAULT_RETRIES,
    retry_delay: float = DEFAULT_RETRY_DELAY,
    generated_at: str = None,
) -> Tuple[Dict[str, VocabularyEntry], List[VocabularyFailure]]:
    """Fetch every vocabulary, retrying each one, and keep going when one service stays down.

    Returns the successfully-fetched entries and one ``VocabularyFailure`` per service that failed
    every attempt -- the caller decides what to do with a partial result (see ``carry_over_previous``).
    """
    generated_at = generated_at or utcnow_iso()
    entries: Dict[str, VocabularyEntry] = {}
    failures: List[VocabularyFailure] = []

    for key, spec in specs.items():
        last_error = None
        for attempt in range(1, max(retries, 1) + 1):
            try:
                values = parsing.data_to_tuple(spec.fetch())
            except Exception as error:  # any failure of a flaky external service is retryable here
                last_error = error
                print(f'[warn] {key}: attempt {attempt} failed ({error!r})', file=sys.stderr)
                if attempt < max(retries, 1) and retry_delay:
                    time.sleep(retry_delay)
                continue
            entries[key] = VocabularyEntry(source=spec.source, values=values, generated_at=generated_at)
            print(f'[ok] {key}: {len(values)} values from {spec.source}')
            if not values:
                # Not an error (the service answered), but it would ship an empty dropdown -- either
                # the service has changed shape, or the io.get_legal_* parser no longer matches it.
                print(
                    f'[warn] {key}: the fetch succeeded but returned 0 values - its GUI dropdown '
                    f'will be empty. Check io.get_legal_* against the current response format '
                    f'before shipping this snapshot.',
                    file=sys.stderr,
                )
            break
        else:
            failures.append(VocabularyFailure(key=key, reason=repr(last_error)))

    return entries, failures


def read_previous_snapshot(path: Path) -> dict:
    """Read an existing snapshot, or return ``{}`` if it is missing or unreadable."""
    try:
        with open(path, encoding='utf-8') as snapshot_file:
            payload = json.load(snapshot_file)
    except (OSError, ValueError):
        return {}
    return payload if isinstance(payload, dict) else {}


def carry_over_previous(
    entries: Dict[str, VocabularyEntry], previous: dict, failed_keys: Sequence[str]
) -> Tuple[Dict[str, VocabularyEntry], List[str]]:
    """Fill failed vocabularies in from the previous snapshot, marked ``stale``.

    A service being down for an hour must not ship a release whose taxon dropdown is empty; keeping
    the previous release's values (flagged, so the maintainer can see it) is strictly better.
    """
    merged = dict(entries)
    carried = []
    previous_vocabularies = previous.get('vocabularies', {}) if isinstance(previous, dict) else {}

    for key in failed_keys:
        previous_entry = previous_vocabularies.get(key)
        if not isinstance(previous_entry, dict) or not previous_entry.get('values'):
            continue
        merged[key] = VocabularyEntry(
            source=previous_entry.get('source', ''),
            values=parsing.data_to_tuple(previous_entry['values']),
            generated_at=previous_entry.get('generated_at', ''),
            stale=True,
        )
        carried.append(key)

    return merged, carried


def build_snapshot(entries: Dict[str, VocabularyEntry], version: str = None, generated_at: str = None) -> dict:
    """Assemble the JSON payload: the vocabularies plus the metadata that dates and sources them."""
    payload = {
        '_comment': GENERATED_COMMENT,
        'generated_at': generated_at or utcnow_iso(),
        'rnalysis_version': version or rnalysis.__version__,
        'vocabularies': {},
    }
    for key, entry in entries.items():
        vocabulary = {'source': entry.source, 'generated_at': entry.generated_at, 'values': list(entry.values)}
        if entry.stale:
            vocabulary['stale'] = True
        payload['vocabularies'][key] = vocabulary
    return payload


def write_snapshot(payload: dict, path: Path) -> None:
    """Write the snapshot as UTF-8 JSON (species names are not all ASCII)."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, 'w', encoding='utf-8') as snapshot_file:
        json.dump(payload, snapshot_file, indent=2, ensure_ascii=False)
        snapshot_file.write('\n')


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description='Regenerate rnalysis/data_files/api_vocabularies.json from the live services '
        '(UniProtKB, PantherDB, Ensembl, PhylomeDB). Run at release time, then commit the result.',
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        '--out',
        default=str(DEFAULT_OUTPUT_PATH),
        metavar='PATH',
        help=f'output file path (default: {DEFAULT_OUTPUT_PATH})',
    )
    parser.add_argument(
        '--retries',
        type=int,
        default=DEFAULT_RETRIES,
        help=f'attempts per service before giving up on it (default: {DEFAULT_RETRIES})',
    )
    parser.add_argument(
        '--retry-delay',
        type=float,
        default=DEFAULT_RETRY_DELAY,
        help=f'seconds to wait between attempts (default: {DEFAULT_RETRY_DELAY})',
    )
    return parser


def main(argv=None) -> int:
    args = build_arg_parser().parse_args(argv)
    out_path = Path(args.out)

    generated_at = utcnow_iso()
    entries, failures = collect_vocabularies(
        default_specs(), retries=args.retries, retry_delay=args.retry_delay, generated_at=generated_at
    )
    if failures:
        entries, carried = carry_over_previous(
            entries, read_previous_snapshot(out_path), [failure.key for failure in failures]
        )
        for failure in failures:
            carried_note = (
                " (kept the previous snapshot's values, marked stale)"
                if failure.key in carried
                else ' (NOT in the snapshot - its dropdown will be empty)'
            )
            print(f'[FAIL] {failure.key}: service stayed down: {failure.reason}{carried_note}', file=sys.stderr)

    write_snapshot(build_snapshot(entries, generated_at=generated_at), out_path)
    print(f'[ok] wrote {out_path} ({len(entries)} vocabularies, RNAlysis {rnalysis.__version__})')
    if failures:
        print('Review the failures above before committing this snapshot.', file=sys.stderr)
    return len(failures)


if __name__ == '__main__':
    raise SystemExit(main())
