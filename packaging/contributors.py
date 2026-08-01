#!/usr/bin/env python3
"""Maintain RNAlysis contributor acknowledgments.

Two capabilities, both dependency-free (stdlib only) so the CI workflow needs no ``pip install``:

* ``--update-readme``        append any GitHub contributor missing from the contributor roster in
                             README.rst. Never removes anyone; hand-added names without a GitHub
                             handle (e.g. collaborators who never opened a PR) are preserved.
* ``--emit-release-thanks``  print a reStructuredText "Thanks to ..." block crediting everyone whose
                             pull request was merged since the previous release, for pasting into the
                             new section of HISTORY.rst.

Both are best-effort: the pure text helpers are unit-tested (tests/test_contributors.py); the GitHub
API calls degrade gracefully. Run ``python packaging/contributors.py --help`` for usage.
"""
from __future__ import annotations

import argparse
import json
import os
import sys
import re
import urllib.error
import urllib.request
from pathlib import Path
from typing import Iterable, List, Tuple
from urllib.parse import quote

GITHUB_API = "https://api.github.com"
DEFAULT_REPO = "GuyTeichman/RNAlysis"
README_PATH = Path(__file__).resolve().parent.parent / "README.rst"

# reStructuredText comment markers that bracket the auto-maintained contributor bullet list.
CONTRIB_START = ".. contributors-list-start"
CONTRIB_END = ".. contributors-list-end"

# Logins to never auto-add (the maintainer is already credited as Development Lead).
EXCLUDE_LOGINS = frozenset({"guyteichman"})

_HANDLE_RE = re.compile(r"github\.com/([A-Za-z0-9](?:[A-Za-z0-9-]*[A-Za-z0-9])?)", re.IGNORECASE)


# --------------------------------------------------------------------------- pure text helpers

def extract_existing_handles(region: str) -> set:
    """Return the lower-cased GitHub handles already linked in a chunk of text."""
    return {m.lower() for m in _HANDLE_RE.findall(region)}


def _contributor_bullet(login: str) -> str:
    return f"* `{login} <https://github.com/{login}>`_"


def merge_contributor_bullets(readme_text: str, logins: Iterable[str]) -> Tuple[str, List[str]]:
    """Append GitHub *logins* not already present to README's contributor region.

    Returns ``(new_text, added_logins)``. Never removes or reorders existing entries; existing
    hand-written names (with or without a handle) are preserved. De-duplicates against both the
    existing handles and repeats within *logins*, case-insensitively. Idempotent. Raises
    ``ValueError`` if the region markers are missing.
    """
    try:
        start = readme_text.index(CONTRIB_START)
        end = readme_text.index(CONTRIB_END, start)
    except ValueError as e:
        raise ValueError(
            f"README is missing the {CONTRIB_START!r} / {CONTRIB_END!r} markers") from e

    existing = extract_existing_handles(readme_text[start:end])

    added: List[str] = []
    seen = set(existing)
    for login in logins:
        key = login.lower()
        if key in seen:
            continue
        seen.add(key)
        added.append(login)

    if not added:
        return readme_text, []

    head = readme_text[:end].rstrip("\n")
    new_bullets = "\n".join(_contributor_bullet(login) for login in added)
    new_text = f"{head}\n{new_bullets}\n\n{readme_text[end:]}"
    return new_text, added


def format_thanks(logins: Iterable[str]) -> str:
    """Return an rst 'Thanks to ...' line crediting *logins* (sorted, de-duplicated); '' if none."""
    uniq = sorted({login for login in logins if login}, key=str.lower)
    if not uniq:
        return ""
    links = ", ".join(f"`@{login} <https://github.com/{login}>`_" for login in uniq)
    return f"Thanks to {links} for contributing to this release! \U0001f389"


# --------------------------------------------------------------------------- GitHub API (best-effort)

def _api_get(url: str, token: str | None = None):
    req = urllib.request.Request(url, headers={
        "Accept": "application/vnd.github+json",
        "User-Agent": "rnalysis-contributors-script",
        "X-GitHub-Api-Version": "2022-11-28",
    })
    if token:
        req.add_header("Authorization", f"Bearer {token}")
    with urllib.request.urlopen(req, timeout=30) as resp:
        payload = json.loads(resp.read().decode("utf-8"))
        link = resp.headers.get("Link", "")
    return payload, link


def _next_link(link_header: str) -> str | None:
    for part in link_header.split(","):
        section = part.split(";")
        if len(section) >= 2 and 'rel="next"' in section[1]:
            return section[0].strip().strip("<>")
    return None


def fetch_all_contributor_logins(repo: str = DEFAULT_REPO, token: str | None = None) -> List[str]:
    """Return GitHub contributor logins (most contributions first), excluding bots."""
    url = f"{GITHUB_API}/repos/{repo}/contributors?per_page=100&anon=0"
    logins: List[str] = []
    while url:
        data, link = _api_get(url, token)
        for c in data:
            login = c.get("login")
            if login and c.get("type") != "Bot" and not login.endswith("[bot]"):
                logins.append(login)
        url = _next_link(link)
    return logins


def get_previous_release_date(repo: str = DEFAULT_REPO, token: str | None = None) -> str | None:
    try:
        data, _ = _api_get(f"{GITHUB_API}/repos/{repo}/releases/latest", token)
    except urllib.error.HTTPError:
        return None
    return data.get("published_at")


def fetch_merged_pr_authors_since(since_iso: str, repo: str = DEFAULT_REPO,
                                  token: str | None = None) -> List[str]:
    """Return logins of authors whose PRs merged after *since_iso* (excludes bots)."""
    q = quote(f"repo:{repo} is:pr is:merged merged:>{since_iso}")
    url = f"{GITHUB_API}/search/issues?q={q}&per_page=100"
    authors: List[str] = []
    while url:
        data, link = _api_get(url, token)
        for pr in (data.get("items", []) if isinstance(data, dict) else []):
            user = (pr.get("user") or {}).get("login")
            if user and not user.endswith("[bot]"):
                authors.append(user)
        url = _next_link(link)
    return authors


# --------------------------------------------------------------------------- CLI

def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--repo", default=os.environ.get("GITHUB_REPOSITORY", DEFAULT_REPO))
    parser.add_argument("--readme", default=str(README_PATH))
    parser.add_argument("--since", help="ISO date for --emit-release-thanks (default: last release date)")
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument("--update-readme", action="store_true",
                      help="append missing GitHub contributors to README.rst")
    mode.add_argument("--emit-release-thanks", action="store_true",
                      help="print an rst 'Thanks to ...' block for PRs merged since the last release")
    args = parser.parse_args(argv)

    token = os.environ.get("GITHUB_TOKEN") or os.environ.get("GH_TOKEN")

    if args.update_readme:
        logins = [login for login in fetch_all_contributor_logins(args.repo, token)
                  if login.lower() not in EXCLUDE_LOGINS]
        readme = Path(args.readme)
        new_text, added = merge_contributor_bullets(readme.read_text(encoding="utf-8"), logins)
        if added:
            readme.write_text(new_text, encoding="utf-8")
            print("Added contributors: " + ", ".join(added))
        else:
            print("No new contributors to add.")
        return 0

    since = args.since or get_previous_release_date(args.repo, token)
    if not since:
        print("Could not determine the previous release date; pass --since YYYY-MM-DD.", file=sys.stderr)
        return 1
    thanks = format_thanks(fetch_merged_pr_authors_since(since, args.repo, token))
    if thanks:
        print(thanks)
    else:
        print(f"No merged pull requests found since {since}.", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
