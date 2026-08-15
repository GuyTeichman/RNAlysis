---
name: release
description: >-
  Ordered checklist for cutting an RNAlysis release: dating the HISTORY.rst entry, merging
  `development` into `master`, running `bumpversion`, regenerating and committing the packaged
  API-vocabulary snapshot and the Sphinx docs (both required — `docs/build/` is tracked in git and
  is what the app's in-GUI help links point at, and the vocabulary snapshot is what populates the
  GUI's taxon/gene-ID dropdowns for that version), tagging to trigger the standalone-build
  workflow, and the (manual) PyPI publish.
  Use when Guy — the maintainer — says "cut a release", "release a new version", "ship 4.x.0",
  "bump the version", or asks to merge `development` into `master`. This is a maintainer-only,
  infrequent, multi-file, partly-irreversible procedure (tag push and PyPI upload can't be
  cleanly undone) — codifying the exact order is what prevents a missed step, such as a version
  string left stale in one of five files, stale/orphaned docs pages shipped as the release's
  in-app help, a tag pushed before `HISTORY.rst` has a real date, or a PyPI upload from the
  wrong commit. Also covers a separate **beta lane** — a private freeze smoke-test that builds the
  standalone app through CI (no GitHub Release, no version bump, no commits, no HISTORY date) to
  verify the PyInstaller freeze still works and hand you an executable to test by hand. Use for the
  stable procedure when Guy — the maintainer — says "cut a release", "release a new version",
  "ship 4.x.0", "bump the version", or asks to merge `development` into `master`; use for the beta
  lane when he says "run a beta", "beta build", "freeze smoke-test", or "test the
  PyInstaller/standalone build". Not for routine feature/fix PRs, which follow the normal
  branch-off-`development` flow in `.claude/workflows.md` — this skill covers only the release act
  itself, and several of its steps require the maintainer's own admin credentials, not just any
  contributor's.
---

# Cutting an RNAlysis release

This is a checklist, not a script — every command it lists already exists in the repo
(`bumpversion`, `packaging/*.py`, `gh`, the Pyinstaller workflow); this skill only fixes the
*order* and the details that are easy to get wrong. **Do not run this speculatively.** Confirm
with Guy before the irreversible steps (pushing the tag, `twine upload`) if there's any doubt.

Read `.claude/workflows.md`'s "Release (maintainer)" and "Branching & releases" sections first —
this skill is the expanded version of those.

---

## Two lanes: pick the right one

This skill documents **two distinct procedures**. Don't confuse them:

- **Stable lane** — Preconditions → Steps 1–7 below. The real thing: `development` reaches
  `master`, the version is bumped, the vocabulary snapshot and docs are regenerated, a tag is
  pushed, standalone builds are published as a GitHub Release, and the PyPI package is uploaded.
  Has **irreversible** steps. This is what "cut a release" means.
- **Beta lane** — the [Beta lane](#beta-lane--freeze-smoke-test-no-release) section directly below.
  A throwaway **freeze smoke-test**: its only job is to prove `pyinstaller RNAlysis.spec` still
  produces a launchable macOS + Windows app and to hand you those executables to click through by
  hand. It **publishes nothing and commits nothing.**

**The one rule that binds the two lanes together — and the reason this section exists at all:**
the `HISTORY.rst` date is stamped **only** in the stable lane (Preconditions #2). A beta must
**never** turn `X.Y.Z (unreleased)` into a dated header, bump the version, or push to
`master`/`development`. The 4.3.0 cycle is the cautionary tale: a beta was tested informally, the
HISTORY header acquired a date while the version was in fact unreleased, and "is 4.3.0 out?" became
ambiguous — to the maintainer, and to any other agent reading the repo. The whole point of the beta
lane is that someone reading the repo mid-beta sees the **exact same** "still unreleased" state as
before. If a beta leaves a dated HISTORY, a bumped version, a new commit, or a tag behind, it has
**failed its contract**, not succeeded.

---

## Beta lane — freeze smoke-test (no release)

**When to use it:** you want to confirm the PyInstaller freeze still works — it breaks in ways the
normal test suite can't catch (missing bundled data files, hidden imports, Qt plugins, one-file vs
one-dir issues) — and to manually exercise the frozen app, *without* cutting a release.

**Mechanism:** `.github/workflows/pyinstaller.yml`'s `workflow_dispatch` trigger has a `beta`
boolean input (**default `true`**). When `beta=true` it **skips the `createrelease` job entirely**
— no commit to `master`, no GitHub Release — and runs only the cross-platform `build` job,
uploading the frozen apps as **run artifacts** you download from the Actions run page. It is the
*exact same* `build` job a real release uses (same PyInstaller spec, same dependency install), so a
green beta is a faithful prediction of the release freeze — the whole reason it's wired into the
same workflow rather than a separate, driftable copy.

```bash
# Dispatch a beta build of whatever is on `development` (or any branch/ref you want to
# freeze-test). `beta` defaults to true, so this alone is enough:
gh workflow run pyinstaller.yml --ref development

# ...or be explicit:
gh workflow run pyinstaller.yml --ref development -f beta=true

# Watch it:
gh run watch --exit-status $(gh run list --workflow pyinstaller.yml --limit 1 --json databaseId -q '.[0].databaseId')
```

When it finishes, open the run's **Artifacts** section, download `RNAlysis-BETA-macos-M1-<sha>` /
`RNAlysis-BETA-windows-<sha>`, unzip, and launch. The build is unsigned, so on macOS clear the
quarantine bit first: `xattr -dr com.apple.quarantine RNAlysis.app` (same as any unsigned local
build).

**What a beta deliberately does NOT do — its isolation contract:**

- ❌ **No `bumpversion`, no `.bumpversion.cfg` edit.** The frozen app self-reports whatever version
  is currently committed (e.g. `4.3.0`); the beta identity lives only in the downloaded artifact's
  name, not inside the app. (Making the *app itself* read `X.Y.Zb1` was considered and dropped —
  see issue #283 — because it requires a committed bump, which is exactly the residue this contract
  forbids.)
- ❌ **No `HISTORY.rst` date stamp.** `X.Y.Z (unreleased)` stays `(unreleased)`. The date is a
  stable-lane act only.
- ❌ **No commit to `master` or `development`, and no new tag.** Nothing is pushed anywhere.
- ❌ **No API-vocabulary regen (Step 3), no docs regen (Step 4).** Those are versioned release
  artifacts; a beta freezes the branch exactly as it stands.
- ❌ **No PyPI upload, no GitHub Release** (not even a pre-release).

**Verify the guardrail held** (worth doing once, after your first beta, to trust the wiring):

```bash
git fetch origin
git log origin/master -1        # unchanged — no auto-commit landed on master
gh release list --limit 3       # no new release/pre-release created for this build
```

**Cleanup:** none. Run artifacts expire on their own (14-day retention, set on the upload step);
there is no branch, tag, commit, or release to delete. That "nothing to clean up" property *is* the
design goal.

---

## A load-bearing fact: four narrow, deliberate exceptions to "never commit to master/development"

CLAUDE.md hard rule 2 says never commit directly to `master` or `development`. The *release
procedure itself* has always broken that rule in a few mechanical, tool-generated spots —
verified both from `CONTRIBUTING.rst`'s "Deploying" section and from git history (e.g. `c16b014e`
"Bump version: 4.1.2 → 4.2.0" sits directly in `master`'s linear history, not behind a PR merge):

1. **The `bumpversion` commit is pushed straight to `master`** (Step 2).
2. **The API-vocabulary snapshot commit is pushed straight to `master`** (Step 3 — same reasoning
   as the docs below: a generated, version-stamped release artifact that depends on the bump).
3. **The docs regeneration commit is pushed straight to `master`** (Step 4 — new; docs are a
   release artifact in this repo, not a side project, so it belongs in the same direct-push
   lineage as the version bump it depends on).
4. **The Pyinstaller workflow's checksum/changelog commit is pushed straight to `master`** (via
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

## Step 3 — Regenerate and commit the API vocabulary snapshot (required)

`rnalysis/data_files/api_vocabularies.json` is the packaged snapshot of the remote vocabularies —
UniProtKB's gene-ID types, and PantherDB's / Ensembl's / PhylomeDB's legal taxons — that
`rnalysis/utils/param_typing.py` bakes into `Literal[...]` annotations, and therefore into the GUI's
dropdowns. It is a **versioned release artifact**, like the quick-start video checksums: the legal
values are pinned per RNAlysis version, so they only ever refresh when a release regenerates them.
Skip this step and the release ships the previous release's taxon lists (a new Ensembl species won't
be in the dropdown — users can still type it into the combo box's free-text field, so this is a
staleness bug, not a blocker).

```bash
python packaging/generate_api_vocabularies.py
git diff --stat rnalysis/data_files/api_vocabularies.json
git commit -am "Regenerate the API vocabulary snapshot for <version>"
```

- Run it **after** Step 2: the file records `rnalysis_version`, which must be the version being
  released.
- It makes **live** requests (UniProtKB REST, PantherDB, Ensembl REST, PhylomeDB FTP), so run it
  from a machine with network access, and **read its output**. It retries each service, then:
  a service that stays down keeps its previous values, marked `"stale": true` (a dead service must
  never silently empty a dropdown) — decide consciously whether to ship that or wait for the
  service to come back. A service that answers but returns *zero* values is reported too, and means
  the matching `io.get_legal_*` parser no longer matches the service's response format — a real bug
  to fix before releasing, not something to ship.
- Sanity-check the diff: entry counts should move by a handful, not collapse. As of 4.3.0 the
  snapshot holds 73 gene-ID types, 144 PantherDB taxons and 356 Ensembl taxons (PhylomeDB is
  empty — a known bug in `io.get_legal_phylomedb_taxons`, not a failed fetch).
- Same direct-push-to-`master` lineage as Steps 2 and 4; one combined push at the end is fine.

## Step 4 — Regenerate and commit the docs (required)

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

## Step 5 — Tag and push, to trigger the standalone build

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

## Step 6 — Publish to PyPI (manual — not automated)

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

## Step 7 — Sync `master` back into `development`

The release commits (`bumpversion`, the vocabulary snapshot, the docs regen, the
checksum/changelog auto-commit) exist
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
- [ ] **`rnalysis/data_files/api_vocabularies.json` regenerated and committed** (Step 3 — its
      `rnalysis_version` is the version just released, `generated_at` is today, no entry is marked
      `"stale": true` unless you consciously chose to ship one, and no entry collapsed to an empty
      list that wasn't empty before).
- [ ] **Docs regenerated and committed; `docs/build/` reflects the new version string** (Step 4 —
      spot-check e.g. `docs/build/index.html` or any page's footer/sidebar version, and confirm
      no leftover stale pages for anything renamed/removed this release). Live at
      `https://guyteichman.github.io/RNAlysis/` once pushed — no separate deploy step.
- [ ] The `V<version>` tag exists on `origin` and points at the right commit.
- [ ] The GitHub Release exists, its body is the correct changelog section (not empty, not a
      stale/previous version), and both `..._macos_M1.zip` and `..._windows.zip` assets are
      attached and non-trivial in size.
- [ ] `rnalysis/gui/videos/checksums/*.txt` and `rnalysis/data_files/latest_changelog.md` were
      auto-committed to `master` by the Pyinstaller workflow.
- [ ] PyPI shows the new version (`https://pypi.org/project/RNAlysis/`), if Step 6 was run.
- [ ] `master` has been merged back into `development` (Step 7).
- [ ] Nothing else broke: a quick `pip install RNAlysis==<version>` in a scratch venv, or
      launching a fresh standalone build, is worth doing before telling anyone it's out.

If any of these is missing, the release isn't done yet — go back to the corresponding step rather
than leaving it for "later."
