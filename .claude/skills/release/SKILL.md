---
name: release
description: >-
  Ordered checklist for cutting an RNAlysis release: dating the HISTORY.rst entry, merging
  `development` into `master`, running `bumpversion`, regenerating and committing the Sphinx
  docs (required — `docs/build/` is tracked in git and is what the app's in-GUI help links
  point at), tagging to trigger the standalone-build workflow, and the (manual) PyPI publish.
  Use when Guy — the maintainer — says "cut a release", "release a new version", "ship 4.x.0",
  "bump the version", or asks to merge `development` into `master`. This is a maintainer-only,
  infrequent, multi-file, partly-irreversible procedure (tag push and PyPI upload can't be
  cleanly undone) — codifying the exact order is what prevents a missed step, such as a version
  string left stale in one of five files, stale/orphaned docs pages shipped as the release's
  in-app help, a tag pushed before `HISTORY.rst` has a real date, or a PyPI upload from the
  wrong commit. Not for routine feature/fix PRs, which follow the normal branch-off-`development`
  flow in `.claude/workflows.md` — this skill covers only the release act itself, and several of
  its steps require the maintainer's own admin credentials, not just any contributor's.
---

# Cutting an RNAlysis release

This is a checklist, not a script — every command it lists already exists in the repo
(`bumpversion`, `packaging/*.py`, `gh`, the Pyinstaller workflow); this skill only fixes the
*order* and the details that are easy to get wrong. **Do not run this speculatively.** Confirm
with Guy before the irreversible steps (pushing the tag, `twine upload`) if there's any doubt.

Read `.claude/workflows.md`'s "Release (maintainer)" and "Branching & releases" sections first —
this skill is the expanded version of those.

## A load-bearing fact: three narrow, deliberate exceptions to "never commit to master/development"

CLAUDE.md hard rule 2 says never commit directly to `master` or `development`. The *release
procedure itself* has always broken that rule in exactly three mechanical, tool-generated spots —
verified both from `CONTRIBUTING.rst`'s "Deploying" section and from git history (e.g. `c16b014e`
"Bump version: 4.1.2 → 4.2.0" sits directly in `master`'s linear history, not behind a PR merge):

1. **The `bumpversion` commit is pushed straight to `master`** (Step 2).
2. **The docs regeneration commit is pushed straight to `master`** (Step 3 — new; docs are a
   release artifact in this repo, not a side project, so it belongs in the same direct-push
   lineage as the version bump it depends on).
3. **The Pyinstaller workflow's checksum/changelog commit is pushed straight to `master`** (via
   the `VID_CHECKSUM_PAT` secret; see `.github/workflows/pyinstaller.yml`).

This only works because `master`'s branch protection requires 1 PR approval **but has
`enforce_admins: false`** (confirmed via `gh api repos/GuyTeichman/RNAlysis/branches/master/protection`)
— i.e. it's bypassable by a repo admin (Guy), not by an arbitrary contributor or bot token. Treat
these steps as **maintainer-only** — don't attempt them with anything less than Guy's own
push access, and don't generalize "the release process pushes to master" into license to commit
there for anything else. `development` has **no** branch protection at all (`404` from the same
API check) — the "never commit directly" rule there is a convention, not an enforced gate, which
makes discipline about it even more important.

---

## Preconditions

1. **Everything committed, `development` up to date and green.**
   ```bash
   git fetch origin
   git log origin/development -1
   gh run list --branch development --workflow "Build CI" --limit 5
   ```
   All three OSes × three Python versions should be green on the latest run. One caveat from the
   CI layout: the network-dependent test tier runs on **only** ubuntu/Python 3.13, so a single red
   network job is a 1-of-1 signal for that tier, not a 1-of-9 flake — look at *why* it's red
   (a dead external API is a known, tolerated failure mode here; a real regression is not).

2. **`HISTORY.rst` has a real, dated entry for the version you're releasing.** At the top of the
   file the newest section currently looks like:
   ```rst
   4.3.0 (unreleased)
   -------------------
   ```
   This **must** become a real ISO date (`4.3.0 (2026-08-10)`) before the release reaches
   `master`. `packaging/generate_changelog.py` extracts the changelog with the regex
   `\(\d{4}-\d{2}-\d{2}\)` — it will **not** match `(unreleased)`, so skipping this step means the
   GitHub Release body and the in-app "what's new" changelog end up wrong or empty. Handy
   coincidence: `unreleased` and `YYYY-MM-DD` are both 10 characters, so the `---` underline
   usually doesn't need re-sizing — but it's an RST section header, so double-check the underline
   is still at least as long as the header text.
   - If there's no unreleased section at the top (nothing merged since the last release), stop —
     there's nothing to release.
   - Make this edit via a normal small PR into `development` (red/green not applicable — it's a
     doc change, but it still goes through the standard branch-off-`development`,
     PR-into-`development` flow; don't hand-edit `master`).

3. **Nothing uncommitted** locally (`git status` clean) on whatever checkout you're releasing
   from.

---

## Step 1 — Merge `development` into `master`

A release *is* `development` reaching `master` (`master` only receives `development` at a version
release — CLAUDE.md rule 2). Do it the same way any change reaches a protected branch: a PR, not a
raw push (precedent: PR #52, "Merge pull request #52 from GuyTeichman/development").

```bash
gh pr create --base master --head development --title "Release <version>" --body-file <tmpfile>
```

Wait for CI on that PR, then merge it **as a merge commit** (not squash/rebase) so `master`'s
history stays a real superset of `development`'s — squashing here would flatten every feature
commit that shipped in this release into one. `master` requires 1 approving review; get it
reviewed like any other PR reaching a protected branch. After it merges:

```bash
git fetch origin
git switch master && git pull
```

## Step 2 — Bump the version

```bash
bumpversion patch    # or minor / major, per semver judgement on what shipped
```

Per `.bumpversion.cfg` (`commit = True`, `tag = False`), this **edits and commits** the version
string across exactly these five files, and does **not** tag:

- `rnalysis/__init__.py` (`__version__ = "..."`)
- `setup.py` (`version="..."`)
- `RNAlysis.spec` (`name='RNAlysis-...'`, the PyInstaller build's output folder name)
- `.github/workflows/pyinstaller.yml` (the version string is baked into the build matrix —
  `CMD_BUILD` zip names and `OUT_FILE_NAME` for both macOS and Windows targets; `bumpversion`'s
  plain string replacement catches every occurrence in the file, so nothing extra to edit by hand)
- `docs/source/conf.py` (`release = '...'`)

This is the direct-push-to-`master` exception described above — run it on `master`, after the
merge in Step 1, then `git push` (per CONTRIBUTING.rst's own documented "Deploying" section). Do
**not** push yet if you're about to do Step 3 next anyway — one combined push is fine.

## Step 3 — Regenerate and commit the docs (required)

**Not optional, not "if you have time."** `docs/build/` (1354 files as of this writing) is
**tracked in git** and is exactly what's live at `https://guyteichman.github.io/RNAlysis/` —
confirmed via `gh api repos/GuyTeichman/RNAlysis/pages`: GitHub Pages serves directly from the
`/docs` folder on `master`, `build_type: "legacy"`, no build step of its own. It's also what the
app's in-GUI help links point at. A stale `docs/build/` after a release means broken in-app help
until the next person happens to rebuild it — and it must happen **after** Step 2, because
`bumpversion` just edited `docs/source/conf.py`'s `release` string, which every generated page
embeds, and because `HISTORY.rst` just got its real release date.

**The commands documented elsewhere for this are wrong for this repo layout — verified live,
don't trust them:**
- CLAUDE.md's Commands section says `make -C docs html` — there is **no `docs/Makefile`**
  (`docs/` only has `source/`, `build/`, `make.bat`, `index.html`), so that command fails outright.
- The root `Makefile`'s `docs` target runs `sphinx-apidoc --seperate -o docs/ rnalysis` (note the
  `--seperate` typo and the wrong output directory) followed by `$(MAKE) -C docs html` and opens
  `docs/_build/html/index.html` — three ways stale: no `docs/Makefile` to invoke, and the real
  committed output directory is `docs/build/`, not `docs/_build/`.
- Despite the repo map/gotchas calling the per-function `.rst` files "sphinx-apidoc generated",
  **no `sphinx-apidoc` step is actually used or needed.** `docs/source/conf.py` sets
  `extensions = [..., 'sphinx.ext.autosummary', ...]` and `autosummary_generate = True`;
  `docs/source/index.rst` has a top-level `.. autosummary:: :toctree:` over the four public
  modules, and the custom templates `docs/source/_templates/autosummary/{module,class}.rst`
  recursively emit a stub `.rst` per class and per public method by introspecting the **live**
  code. A single `sphinx-build` invocation regenerates the stubs *and* builds the HTML in one
  pass — verified locally: running it live regenerated 4 stale class-level `.rst` files and added
  8 missing per-method stub files for `annotate_from_gtf`/`filter_by_gtf_attribute` (merged by
  PR #201), which are **currently missing from `development`'s committed docs right now** — live,
  concrete proof this step is easy to skip and silently drifts stale otherwise.

**The verified, actually-working commands** (run from the repo root, in an environment with
`pip install -e .[all]` + `requirements_dev.txt` — i.e. `Sphinx`/`sphinx_rtd_theme` installed and
the real `rnalysis` package importable, since `autodoc`/`autosummary` introspect live code):

```bash
# Clean slate first: autosummary_generate only adds/overwrites stub files, it never deletes
# one for a method/class that was renamed or removed since the last regen — and a full HTML
# rebuild only overwrites pages for sources that still exist, so a stale page for a removed
# source lingers in docs/build/ (a checked-in, publicly served directory) unless removed by
# hand. Deleting both before rebuilding avoids shipping orphaned pages as "current" docs.
rm -f docs/source/rnalysis.*.rst   # auto-generated API stubs only; narrative pages are untouched
rm -rf docs/build

# Regenerates the stubs (autosummary_generate=True) AND builds the HTML in this one pass.
python -m sphinx -b html docs/source docs/build
```

Verified working (Sphinx 7.3.7, this repo, `development` HEAD at the time of writing): completes
with `build succeeded` and a couple hundred pre-existing warnings (docstring/xref nits — normal
background noise for this project, not new). Watch for `build succeeded` and treat actual Sphinx
**errors** (not warnings) as blocking; a large new spike in warning count vs. the prior release is
worth a look but isn't automatically fatal.

Then commit the result — this is the second direct-push-to-`master` exception described above:

```bash
git add docs/
git commit -m "Regenerate docs for <version>"
git push origin master
```

## Step 4 — Tag and push, to trigger the standalone build

`bumpversion` deliberately does not tag (`tag = False`) — tag manually, matching the existing `V*`
convention (`V4.2.0`, `V4.1.2`, ...; this is also what `.github/workflows/pyinstaller.yml`'s
`on.push.tags: ['V*']` trigger listens for):

```bash
git tag V<version>
git push origin master
git push origin V<version>       # (or: git push --tags)
```

Pushing the tag kicks off the **Pyinstaller Release** workflow (`.github/workflows/pyinstaller.yml`),
which — with no further action needed from you — will:

1. Run `packaging/checksum_videos.py` (MD5-checksums every `rnalysis/gui/videos/*.webp` quick-start
   video into `rnalysis/gui/videos/checksums/*.txt` — this is how the app detects a stale/corrupt
   bundled video) and `packaging/generate_changelog.py` (extracts the just-dated `HISTORY.rst`
   section for this version, converts RST→MyST, writes `rnalysis/data_files/latest_changelog.md` —
   this is both the GitHub Release body and the in-app changelog).
2. **Auto-commit** those two generated artifacts straight to `master` (the third deliberate
   direct-push exception, via the `VID_CHECKSUM_PAT` secret).
3. Create the GitHub Release (tag name, title "Stable release `<tag>`", body = the generated
   changelog).
4. Build standalone packages on `macos-latest` and `windows-latest` via `pyinstaller RNAlysis.spec`,
   zip them (`RNAlysis-<version>_macos_M1.zip`, `RNAlysis-<version>_windows.zip`), and attach them
   to the Release as assets.

Watch it: `gh run watch --exit-status $(gh run list --workflow pyinstaller.yml --limit 1 --json databaseId -q '.[0].databaseId')`.

## Step 5 — Publish to PyPI (manual — not automated)

`CONTRIBUTING.rst`'s "Deploying" section says GitHub Actions "will build and publish... the PyPI
package", but **no workflow in this repo does that** — `pyinstaller.yml` only builds the
standalone zips, and `twine` is a `requirements_dev.txt` dependency with no CI caller. In practice
this is a manual step from a checkout of the just-tagged `master`:

```bash
make dist       # python setup.py sdist bdist_wheel  (per the Makefile's `dist`/`release` targets)
make release    # == make dist && twine upload dist/*   (needs PyPI credentials configured for twine)
```
(`make release` runs `dist` as a prerequisite, so calling it alone is enough.) Treat
`CONTRIBUTING.rst`'s claim of an automated PyPI publish as stale documentation, not as something
to go implement as part of a routine release — flag it to Guy if he wants that gap closed
separately; that's a CI/tooling change of its own, out of scope for the act of releasing.

## Step 6 — Sync `master` back into `development`

The release commits (`bumpversion`, the docs regen, the checksum/changelog auto-commit) exist
only on `master` until they're folded back — otherwise the next feature branch cut from
`development` starts from a stale version string and stale docs. Git history shows this done as a
plain merge (`Merge remote-tracking branch 'origin/master'`), not a PR — a pure content sync with
no new code:

```bash
git switch development && git pull
git merge origin/master
git push origin development
```
(`development` has no branch protection, so this is a normal push — but keep it to *only* this
sync; don't use the lack of protection as license to push unrelated work there directly.)

---

## Post-release verification (definition of done)

- [ ] `HISTORY.rst`'s new-version header carries a real date, on both `master` and `development`.
- [ ] The version string matches across all five `bumpversion`-managed files (spot-check
      `rnalysis/__init__.py`) on `master`.
- [ ] **Docs regenerated and committed; `docs/build/` reflects the new version string** (Step 3 —
      spot-check e.g. `docs/build/index.html` or any page's footer/sidebar version, and confirm
      no leftover stale pages for anything renamed/removed this release). Live at
      `https://guyteichman.github.io/RNAlysis/` once pushed — no separate deploy step.
- [ ] The `V<version>` tag exists on `origin` and points at the right commit.
- [ ] The GitHub Release exists, its body is the correct changelog section (not empty, not a
      stale/previous version), and both `..._macos_M1.zip` and `..._windows.zip` assets are
      attached and non-trivial in size.
- [ ] `rnalysis/gui/videos/checksums/*.txt` and `rnalysis/data_files/latest_changelog.md` were
      auto-committed to `master` by the Pyinstaller workflow.
- [ ] PyPI shows the new version (`https://pypi.org/project/RNAlysis/`), if Step 5 was run.
- [ ] `master` has been merged back into `development` (Step 6).
- [ ] Nothing else broke: a quick `pip install RNAlysis==<version>` in a scratch venv, or
      launching a fresh standalone build, is worth doing before telling anyone it's out.

If any of these is missing, the release isn't done yet — go back to the corresponding step rather
than leaving it for "later."
