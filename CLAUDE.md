# CLAUDE.md — working in RNAlysis

Operational guide for Claude Code (and any AI/agent) working in this repository.
Read this first. Deeper material lives in [`.claude/context.md`](.claude/context.md),
[`.claude/workflows.md`](.claude/workflows.md), and [`.claude/design-philosophy.md`](.claude/design-philosophy.md).

> **Not using Claude Code?** This guide applies to you too. [`AGENTS.md`](AGENTS.md) is the
> provider-neutral entry point for any AI agent or automated contributor — it distills the
> non-negotiable rules and points back here for the full detail. `CLAUDE.md` (this file) is the
> canonical, complete version; keep the two in sync when you change the rules.

---

## What RNAlysis is

RNAlysis is desktop software that lets biologists analyze RNA-sequencing data —
from raw FASTQ through filtering, normalization, differential expression, clustering,
and enrichment — **without writing a single line of code**. It ships both a PyQt6
graphical app and a programmatic Python API, and it targets scientists with **zero
programming experience**. Correctness and reproducibility are non-negotiable
(see [design-philosophy](.claude/design-philosophy.md)).

---

## The one thing to understand first: the API↔GUI reflection contract

**The programmatic API is the source of truth. The GUI is generated from it by reflection.**
Internalize this before changing anything, because it explains most of the conventions.

When you add a public method to a `Filter`/`FeatureSet` subclass (or a public function in
`fastq`/`enrichment`), the GUI **automatically discovers and renders it** — you do not write
Qt code for it. Three things drive that rendering:

1. **Discovery** — `gui.py::get_all_actions()` calls `dir()` on the object and keeps every
   public (non-`_`) callable not listed in that widget's `EXCLUDED_FUNCS`. Methods are sorted
   into GUI tabs by **name heuristics**: contains `normalize` → *Normalize*; contains `filter`
   or `split` → *Filter*; membership in `SUMMARY_FUNCS`/`CLUSTERING_FUNCS`/`GENERAL_FUNCS` →
   those tabs; otherwise → *Visualize*.

2. **Parameter widgets** — `gui_widgets.py::param_to_widget()` maps each parameter's **type
   annotation** to a specific Qt widget. **Type annotations are load-bearing UI.** Examples:
   - `bool` → toggle switch · `int` → spin box · `PositiveInt`/`NonNegativeInt`/`NegativeInt` → range-bounded spin box
   - `Fraction` → 0–1 slider/spin · `Color` → color picker · `ColorList` → multi-color picker · `ColorMap` → colormap combo
   - `ColumnName`/`ColumnNames`/`GroupedColumns` → table-column pickers
   - `Literal['a','b']` → combo box · `Union[..., Literal[...]]` → combo box with an "other…" field
   - `Union[T, None]` (i.e. `Optional[T]`) → a checkbox that enables/disables a `T` widget
   - `Sequence[str]` → gene-set combo box
   These custom types live in [`rnalysis/utils/param_typing.py`](rnalysis/utils/param_typing.py).
   **Choosing the right annotation is how you design the GUI.**

3. **Labels & help text** — the `@readable_name('Human readable label')` decorator
   (`utils/generic.py`) sets the button/window title; the reStructuredText `:param x: ...`
   docstring lines become the per-parameter help/tooltips. A public method **without**
   `@readable_name` still appears (using its function name), so hide internals with a leading
   underscore or an `EXCLUDED_FUNCS`/`EXCLUDED_PARAMS` entry.

**Consequence:** a rename, a signature change, or a changed/loosened type annotation on a
public API method silently changes the GUI. Treat any such change as *risky* (see rules below).

---

## Repository map

```
rnalysis/
  __init__.py            # version, settings keys, FROZEN_ENV flag
  general.py             # misc top-level helpers (parsing IDs, saving tables, settings paths)
  filtering.py           # PUBLIC API — Filter, CountFilter, DESeqFilter, FoldChangeFilter, Pipeline
  enrichment.py          # PUBLIC API — FeatureSet, RankedSet (GO/KEGG/user enrichment, set ops)
  fastq.py               # PUBLIC API — adapter trimming, alignment, counting pipelines (wrap CLI tools)
  utils/                 # IMPLEMENTATION (not user-facing)
    io.py                # web services (async aiohttp + rate-limit + cache), subprocess, R runner
    enrichment_runner.py # enrichment statistics engines
    clustering.py        # clustering algorithms
    differential_expression.py  # builds/executes R scripts for DESeq2 & limma-voom
    feature_counting.py  # featureCounts (Rsubread) bridge
    genome_annotation.py # GTF/GFF parsing & gene-length logic
    ontology.py          # GO/KEGG DAG handling
    param_typing.py      # custom annotation types that drive GUI widgets  ← load-bearing
    generic.py           # readable_name, GenericPipeline, parallel helpers, misc
    parsing.py, validation.py, settings.py, installs.py
  gui/                   # PyQt6 app — reflects the API; you rarely add per-function code here
    gui.py               # main window, action discovery, worker dispatch
    gui_widgets.py       # param_to_widget() and all custom widgets
    gui_windows.py       # FuncExternalWindow (auto-built function dialogs)
    gui_report.py        # interactive analysis report (networkx + pyvis)
    gui_graphics.py, gui_quickstart.py, gui_style.py, main.py
  data_files/
    r_templates/         # parametrized R scripts (DESeq2, limma, Rsubread, install, sessioninfo)
    report_templates/, report_misc/
tests/                   # pytest + pytest-qt; one test_*.py per module; fixtures in test_files/
docs/source/             # Sphinx sources; per-function .rst are sphinx-apidoc GENERATED (don't hand-edit)
packaging/               # changelog + video-checksum helpers for releases
```

Class hierarchy (public API): `Filter` → `CountFilter`, `DESeqFilter`, `FoldChangeFilter`;
`FeatureSet` → `RankedSet`; `GenericPipeline` → `filtering.Pipeline`, `fastq.*Pipeline`.

---

## Hard rules (do not violate without explicit sign-off)

1. **TDD, always.** Red → green → refactor. Write a failing test first, make it pass, then
   clean up. New code must be covered. Use the `tdd` skill. Run the *relevant* test module(s)
   before calling anything done, and say plainly what you could **not** run (R / network / Qt paths).
2. **Branch model: feature branch → `development` → `master`.** `development` is the long-lived
   integration branch; `master` is the stable/released default branch and only receives
   `development` at a version release. **Never commit directly to `master` or `development`** —
   cut a feature/fix branch off `development` and open a PR **into `development`**.
3. **Plan-first for risky changes.** If a change touches a public API signature, the API→GUI
   reflection contract, or a serialized format (Pipeline YAML / session file / exported
   parameters), propose a short plan and wait for approval. Small localized fixes: just do it
   (TDD) and show the diff.
4. **Back-compat of serialized artifacts is a hard invariant.** Old Pipeline YAMLs, saved
   sessions, and exported parameter files **must keep loading** in new versions. Breaking this
   is a bug — flag it loudly.
5. **Results must not change arbitrarily between versions.** A given analysis with given
   parameters should reproduce. If a change alters numerical output, it must be intentional,
   justified, and called out in `HISTORY.rst`.
6. **Every new function needs reasonable defaults.** Assume the user is a non-programmer,
   non-bioinformatician. Defaults should produce a sensible, correct analysis out of the box.
7. **Cross-platform or it doesn't ship.** Must work on Windows, macOS, and Linux (clusters),
   both run-from-source **and** frozen (PyInstaller). See multiprocessing gotcha below.
8. **At PR time, an independent agent with clean context reviews the diff** before it's
   considered done (use the `code-review` skill or spawn a fresh subagent — not your own biased read).

---

## Commands

Environment: Windows, PowerShell primary (a Bash tool is also available). Python 3.10–3.15.

```bash
# Install for development
pip install -e .[all]
pip install -r requirements_dev.txt      # pytest, pytest-qt, etc.

# Run the app
rnalysis-gui                              # installed entry point
python rnalysis_app.py                    # from a checkout
python -c "from rnalysis import gui; gui.run_gui()"

# Tests (full suite is long and needs R, kallisto, bowtie2, network, a Qt display)
pytest tests/test_filtering.py            # a single module (preferred while iterating)
pytest tests/test_filtering.py -k <name>  # a subset
coverage run --source=rnalysis/ -m pytest tests/ && coverage report -m   # what CI does

# Docs (per-function .rst are generated — see workflows.md)
make -C docs html

# Release (maintainer)
bumpversion patch        # or minor / major — updates version across files per .bumpversion.cfg
```

CI (`.github/workflows/build_ci.yml`) runs on every PR across
{ubuntu, windows, macos} × Python {3.12, 3.13, 3.14}, installing R, kallisto, and bowtie2.

---

## Gotchas / tread carefully

- **External web APIs are the most fragile part of the codebase.** UniProt, Ensembl, PANTHER,
  PhylomeDB, OrthoInspector, KEGG, and GO change formats, rename fields, and go down. Assume
  breakage: keep the retry (`tenacity`) + rate-limit (`aiolimiter`) + cache scaffolding in
  `io.py`, degrade gracefully, and never let one dead service crash the app. (The current open
  PR is fixing ≥2 such services.)
- **Multiprocessing is fragile because of the frozen-vs-source split.** `FROZEN_ENV`
  (PyInstaller) restricts `PARALLEL_BACKENDS` to `multiprocessing`/`sequential`; from source it
  also allows `loky`/`threading`. The existing choices are deliberate — change them only with
  care and test **both** frozen and source paths in your head before touching them.
- **GUI/qtbot tests are flaky in CI in ways that don't reproduce locally.** This is a known,
  unresolved problem. Don't assume a green local run means green CI; don't "fix" a flaky GUI
  test by weakening its assertions without understanding it.
- **R interface = permission pain.** On fresh machines / CI, installing or updating R packages
  (DESeq2, limma, Rsubread) often fails on write permissions to the R library. Error paths in
  `installs.py` already surface this — preserve those messages.
- **Docstrings are load-bearing** (they become GUI help). Keep the `:param name: ...` /
  `:type name: ...` reST format when editing public functions.
- **`docs/source/rnalysis.*.rst` are generated** by `sphinx-apidoc`. Don't hand-edit them.

---

## Style

- 4-space indent, LF, UTF-8, trailing whitespace trimmed, final newline (`.editorconfig`).
- No enforced formatter (no `black`/`ruff`/`pre-commit` config). Match the style of the file
  you're editing. `flake8` is referenced but unconfigured — keep lines reasonable, imports tidy.
- Prefer **Polars** for tabular work (the codebase migrated off Pandas in 4.0). Lazy evaluation
  is currently underused and welcome where it helps (see roadmap in context.md).
- Give public API parameters **precise `param_typing` annotations** — that's how you design the
  GUI. Add `@readable_name(...)` and a reST docstring with a runnable `:Examples:` block.

---

## Relevant skills

`tdd` (mandatory workflow), `diagnosing-bugs` (hard bugs / CI-only flakiness),
`code-review` (the clean-context PR review), `research` (nailing down an external API's real
current behavior before coding against it).
