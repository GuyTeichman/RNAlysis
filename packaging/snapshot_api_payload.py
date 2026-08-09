#!/usr/bin/env python
"""Fetch a URL and save its raw response body to disk as a test fixture.

RNAlysis leans on external web services (UniProt, Ensembl, PANTHER, PhylomeDB, OrthoInspector,
KEGG, GO, ...) that change formats, rename fields, and go down without notice -- see CLAUDE.md:
"External web APIs are the most fragile part of the codebase." Whenever you fix or harden code
that talks to one of these services, the accompanying test must NOT depend on the live network --
mock it. But an honest mock needs a REAL, current payload to mock against.

This script is that capture step: it performs a plain GET, writes the raw response body to disk
byte-for-byte (no parsing, no reformatting), and prints where it landed so you can wire it into a
test fixture (conventionally under ``tests/test_files/``) and mock the request against it -- e.g.
with ``requests_mock``, the convention already used throughout ``tests/test_io.py``.

This is a provider-neutral helper (usable by any agent, by hand, or from CI): it has no RNAlysis
import beyond the standard library and ``requests`` (already a core dependency).

Usage
-----
    # basic capture -- filename derived from the URL, saved under tests/test_files/
    python packaging/snapshot_api_payload.py https://rest.uniprot.org/idmapping/status/abc123

    # with query parameters and custom headers
    python packaging/snapshot_api_payload.py https://api.example.org/search \\
        --param query=BRCA1 --param format=json \\
        --header "Accept: application/json" \\
        --out tests/test_files/example_search_response.json

    # a bare filename (no path separators) still lands under tests/test_files/
    python packaging/snapshot_api_payload.py https://example.org/data --out my_fixture.txt

Notes
-----
* GET only, on purpose -- this captures read responses for test fixtures, not a general-purpose
  HTTP client. Add other verbs only if a real workflow needs them.
* When --out is omitted, the filename is derived from the URL's host + path, sanitized to safe
  filesystem characters, and given an extension guessed from the response's Content-Type header
  (falls back to .txt if the type is missing or unrecognized).
* This makes a REAL, live HTTP request over the network -- do not run it in an environment
  without network access (e.g. a sandboxed CI job). Capture the payload once, locally, then
  commit the fixture file and reference it from a mocked test.
* A non-2xx response is treated as a failure (nothing is written) -- if you specifically want to
  capture an error payload (e.g. PANTHER's empty-200 quirk, or a 404 body), that's a case this
  script does not special-case; inspect the response and adapt as needed.
"""
from __future__ import annotations

import argparse
import re
import sys
import urllib.parse
from pathlib import Path
from typing import Dict, Optional, Sequence

import requests

REPO_ROOT = Path(__file__).resolve().parent.parent
DEFAULT_OUTPUT_DIR = REPO_ROOT / 'tests' / 'test_files'

_EXTENSION_BY_CONTENT_TYPE = {
    'application/json': '.json',
    'text/json': '.json',
    'application/xml': '.xml',
    'text/xml': '.xml',
    'text/csv': '.csv',
    'text/html': '.html',
    'text/plain': '.txt',
}
_DEFAULT_EXTENSION = '.txt'
_MAX_STEM_LENGTH = 100


def parse_key_value_pairs(pairs: Sequence[str], flag_name: str) -> Dict[str, str]:
    """Parse repeated ``key=value`` CLI arguments (used for --param) into a dict.

    Splits on the FIRST ``=`` only, so a value that itself contains ``=`` survives intact.
    Raises ``ValueError`` on malformed input; a caller building a CLI should catch it and hand
    the message to ``argparse.ArgumentParser.error`` (kept as a plain function, not tied to
    argparse, so it stays directly unit-testable).
    """
    result: Dict[str, str] = {}
    for pair in pairs:
        if '=' not in pair:
            raise ValueError(f"{flag_name} value {pair!r} is not in key=value form")
        key, _, value = pair.partition('=')
        if not key:
            raise ValueError(f"{flag_name} value {pair!r} has an empty key")
        result[key] = value
    return result


def parse_headers(pairs: Sequence[str]) -> Dict[str, str]:
    """Parse repeated ``Key: value`` CLI arguments (used for --header) into a dict.

    Splits on the FIRST colon only, so a header value that itself contains a colon (e.g. an
    ``Authorization: Bearer abc:def`` token) survives intact.
    """
    result: Dict[str, str] = {}
    for pair in pairs:
        if ':' not in pair:
            raise ValueError(f"--header value {pair!r} is not in 'Key: value' form")
        key, _, value = pair.partition(':')
        key = key.strip()
        value = value.strip()
        if not key:
            raise ValueError(f"--header value {pair!r} has an empty key")
        result[key] = value
    return result


def guess_extension(content_type: Optional[str]) -> str:
    """Map a Content-Type header value to a file extension, defaulting to ``.txt``."""
    if not content_type:
        return _DEFAULT_EXTENSION
    media_type = content_type.split(';', 1)[0].strip().lower()
    return _EXTENSION_BY_CONTENT_TYPE.get(media_type, _DEFAULT_EXTENSION)


def derive_filename(url: str, content_type: Optional[str] = None) -> str:
    """Derive a descriptive, filesystem-safe fixture filename from a URL and response type.

    e.g. ``https://rest.uniprot.org/idmapping/status/123?fields=x`` (JSON content type) ->
    ``rest_uniprot_org_idmapping_status_123.json``. The query string is deliberately dropped --
    it's usually noise (job IDs, tokens) rather than something worth encoding in a filename; use
    ``--out`` for an exact name when that matters.
    """
    parsed = urllib.parse.urlsplit(url)
    raw_stem = f'{parsed.netloc}{parsed.path}'
    stem = re.sub(r'[^A-Za-z0-9]+', '_', raw_stem).strip('_') or 'payload'
    if len(stem) > _MAX_STEM_LENGTH:
        stem = stem[:_MAX_STEM_LENGTH].rstrip('_')
    return f'{stem}{guess_extension(content_type)}'


def resolve_output_path(url: str, out: Optional[str], content_type: Optional[str] = None,
                         default_dir: Path = DEFAULT_OUTPUT_DIR) -> Path:
    """Decide where the payload should be written.

    - ``out`` is ``None``: ``<default_dir>/<filename derived from the URL>``.
    - ``out`` is a bare filename (no directory component): ``<default_dir>/<that filename>``, so
      a maintainer can pass ``--out my_fixture.json`` without typing the full
      ``tests/test_files/`` path.
    - ``out`` contains a directory component (relative or absolute): used as-is (relative paths
      resolve against the current working directory, matching normal CLI behavior).
    """
    if out is None:
        return default_dir / derive_filename(url, content_type)
    out_path = Path(out)
    if out_path.parent == Path('.') and not out_path.is_absolute():
        return default_dir / out_path.name
    return out_path


def write_payload(path: Path, content: bytes) -> None:
    """Write raw bytes to ``path``, creating parent directories as needed."""
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_bytes(content)


def fetch(url: str, params: Dict[str, str], headers: Dict[str, str], timeout: float) -> requests.Response:
    """Perform the live GET request. The only function in this module that touches the network."""
    response = requests.get(url, params=params, headers=headers, timeout=timeout)
    response.raise_for_status()
    return response


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description='Fetch a URL and save its raw response body to disk as a test fixture '
                     '(e.g. to mock an external web API in a test without hitting the live service).',
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument('url', help='URL to GET, e.g. https://rest.uniprot.org/idmapping/status/abc123')
    parser.add_argument('--param', action='append', default=[], metavar='KEY=VALUE',
                         help='query parameter to add to the request; may be repeated')
    parser.add_argument('--header', action='append', default=[], metavar='"Key: value"',
                         help='HTTP header to send; may be repeated')
    parser.add_argument('--out', default=None, metavar='PATH',
                         help='output file path. A bare filename (no directory) lands under '
                              'tests/test_files/; a path with a directory is used as-is. '
                              'Default: a filename derived from the URL.')
    parser.add_argument('--timeout', type=float, default=30.0,
                         help='request timeout in seconds (default: 30)')
    return parser


def main(argv=None) -> int:
    parser = build_arg_parser()
    args = parser.parse_args(argv)

    try:
        params = parse_key_value_pairs(args.param, '--param')
        headers = parse_headers(args.header)
    except ValueError as exc:
        parser.error(str(exc))  # prints usage + message and raises SystemExit(2)

    try:
        response = fetch(args.url, params, headers, args.timeout)
    except requests.RequestException as exc:
        print(f'[FAIL] could not fetch {args.url}: {exc}', file=sys.stderr)
        return 1

    content_type = response.headers.get('Content-Type')
    out_path = resolve_output_path(args.url, args.out, content_type, DEFAULT_OUTPUT_DIR)
    write_payload(out_path, response.content)

    print(f'[ok] {args.url} -> {out_path} ({len(response.content)} bytes, status {response.status_code})')
    print('Remember: commit this fixture and mock the request against it (e.g. requests_mock) -- '
          'do not let a test depend on the live service.')
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
