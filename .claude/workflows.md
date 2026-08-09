# workflows.md — typical recipes for working in RNAlysis

Step-by-step playbooks. Assumes the rules in [`../CLAUDE.md`](../CLAUDE.md) (TDD, plan-first for
risky changes, branch-not-`master`, back-compat invariant, clean-context PR review).

Every workflow ends the same way → **[Finishing a change](#finishing-a-change-pr--review)**.

---

## Add or change a public API function (the common case)

Because the GUI is reflected from the API, adding one well-annotated method gives you a working
GUI dialog for free. Do it TDD.

1. **Decide it's not risky, or plan first.** New method with a new name = usually fine to just
   build. Changing an existing public signature / renaming / loosening a type annotation =
   *risky* (it mutates the GUI and possibly serialized params) → propose a short plan, wait for OK.
2. **Red.** Add a failing test in the matching `tests/test_*.py` using fixtures from
   `tests/test_files/`. Assert on real expected values (correctness is paramount).
3. **Green.** Implement the method on the right class (`Filter`/`CountFilter`/`DESeqFilter`/
   `FoldChangeFilter`/`FeatureSet`/`RankedSet`) or module. Requirements:
   - `@readable_name('Human-readable label')` from `utils.generic`.
   - **Precise type annotations** using `utils.param_typing` types — this designs the GUI widget
     (see the reflection table in `CLAUDE.md`). `Fraction`, `PositiveInt`, `Color`, `ColorMap`,
     `ColumnName(s)`, `Literal[...]`, `Optional[...]`, etc.
   - A reST docstring: one-line summary, `:param x:` / `:type x:` for **every** parameter
     (these become GUI help), and a runnable `:Examples:` block.
   - **Reasonable defaults** for every parameter — a non-expert should get a sensible result.
   - Follow the tab heuristic if you want a specific GUI tab: name it `*normalize*` / `*filter*`
     / `*split*`, or add it to the relevant `SUMMARY_FUNCS`/`CLUSTERING_FUNCS`/`GENERAL_FUNCS`
     set in `gui.py`.
   - To hide a helper from the GUI: prefix with `_`, or add to the widget's `EXCLUDED_FUNCS`.
4. **Verify GUI reflection** if the type is unusual: check `param_to_widget()` already handles
   the annotation. A brand-new annotation type means teaching `param_to_widget` a new branch
   (and testing it in `test_gui_widgets.py`).
5. **Refactor**, run the module's tests, update `HISTORY.rst`, and add the feature to
   `README.rst` if it's user-facing (per `CONTRIBUTING.rst`).
6. **Screenshot the dialog** (if the change is visible in the GUI — new/renamed/retyped
   parameter, new/edited `@readable_name`, edited `:param:` help). Rebuild the reflected dialog
   from the API and attach it to the PR — see [Attach GUI screenshots to a PR](#attach-gui-screenshots-to-a-pr).

## Fix a bug

1. Use the `diagnosing-bugs` skill for anything non-obvious. Reproduce first.
2. **Red:** write a test that fails *because of the bug*.
3. **Green:** fix it. **Refactor:** clean up. Run the relevant module's tests.
4. If the fix changes any numerical/analysis output, that's a reproducibility event — make it
   intentional and note it in `HISTORY.rst`.

## Change code that touches an external web API (UniProt/Ensembl/PANTHER/KEGG/GO/…)

Highest-risk area (see context.md).

1. **Confirm current real behavior** of the service before coding — hit it live or use the
   `research` skill. Don't trust stale code or memory; these change often.
2. Keep the `io.py` scaffolding: `tenacity` retries, `aiolimiter` rate limiting, and the
   response caches. Degrade gracefully; a dead service must not crash the app.
3. Tests here can't depend on the live service — mock the responses. Capture a real current
   response payload as a fixture in `tests/test_files/`.

## Add/modify differential expression or feature counting (R bridge)

1. R logic lives in **`rnalysis/data_files/r_templates/*.R`** as parametrized templates;
   `utils/differential_expression.py` and `utils/feature_counting.py` fill and run them via
   `io.run_r_script`. Prefer editing the template + the Python that parametrizes it over
   inlining R strings.
2. Preserve the permission-aware error messages in `installs.py` (`install_deseq2`, etc.).
3. Keep exported DE parameters loadable by older sessions (back-compat invariant).

## Add a GUI-only element (rare)

Most functions need no GUI code. Genuinely GUI-only work (a new window, a new widget type, a
report feature) goes in `gui/` and is tested in the matching `tests/test_gui_*.py`. Beware the
known qtbot CI flakiness — write deterministic, well-scoped assertions.

## Attach GUI screenshots to a PR

**Required for any visible GUI change** (rule 9 in `CLAUDE.md`). Because the GUI is reflected from
the API, a code diff doesn't *show* the visual change — so attach a screenshot of the affected
dialog. The dialog can be rebuilt straight from the API, headlessly (no clicking through the app):

```bash
# in the dev env (PyQt6 installed); one or more dotted targets:
python packaging/capture_gui_dialog.py filtering.CountFilter.filter_low_reads --out .gui_shots
python packaging/capture_gui_dialog.py --help     # options: --exclude, --show-all, --offscreen
```

- **When:** a new/renamed/removed parameter, a changed type annotation, a new/edited
  `@readable_name`, or edited `:param:` help text, on a public `Filter`/`FeatureSet`/`fastq`/
  `enrichment` function. For a rename/retype, capture **before** (from `development`) **and after**.
  **Skip** for behind-the-scenes changes with no visible dialog effect.
- **Where:** publish PNGs to the **`assets` branch** under `pr-<N>/` (namespaced so parallel PRs
  never collide), then link the raw URLs from the PR. Do the publish in a throwaway git worktree so
  it never disturbs your working tree.
- **Fidelity/limits:** faithful for the parameter form; does **not** capture window title bars or
  display-level flows (loaded data, result plots, hover/menus) — those stay a human task.

The `gui-screenshots` skill automates the full capture → publish → attach flow, including the
assets-branch worktree recipe and the collision-safe push. On a headless Linux CI runner, wrap the
capture in `xvfb-run` so fonts resolve (Qt's bare `offscreen` platform renders text as tofu).

## Regenerate documentation

- Per-function `docs/source/rnalysis.*.rst` are **generated by `sphinx-apidoc`** — never
  hand-edit. Regenerate with `sphinx-apidoc` (see `Makefile`'s `docs` target) and build with
  `make -C docs html`. Narrative pages (`quick_start.rst`, `faq.rst`, `installation.rst`, etc.)
  are hand-written and safe to edit.

## Release (maintainer)

1. Ensure everything's committed, tests pass in CI, and `HISTORY.rst` has an entry.
2. `bumpversion patch|minor|major` — updates the version string across `rnalysis/__init__.py`,
   `setup.py`, `RNAlysis.spec`, the PyInstaller workflow, and `docs/source/conf.py` (per
   `.bumpversion.cfg`; it commits, does not tag).
3. Standalone builds come from `.github/workflows/pyinstaller.yml`; `packaging/` holds the
   changelog + quick-start-video checksum helpers.

---

## Branching & releases

- **`development`** = long-lived integration branch; **`master`** = stable/released default branch.
- Work happens on **feature/fix branches cut from `development`**, merged back into `development`
  via PR. `master` receives `development` only at a **version release** (then `bumpversion` +
  the PyInstaller build workflow produce the release).
- So: never commit directly to `master` or `development`; branch off `development` and target
  your PR at `development`.

## Finishing a change (PR + review)

1. Confirm you're on a **feature branch cut from `development`** (never `master`/`development`
   directly), and that the PR will **target `development`**.
2. Relevant tests green locally; state clearly what you couldn't run (R / network / Qt).
3. `HISTORY.rst` updated; `README.rst` updated if user-facing.
4. **Visible GUI change? Attach dialog screenshot(s)** (rule 9) — see
   [Attach GUI screenshots to a PR](#attach-gui-screenshots-to-a-pr). Note anything you couldn't
   render (e.g. a display-level flow) and why.
5. Open the PR (`gh`). PR body ends with the Claude Code attribution footer.
6. **Trigger an independent clean-context review of the diff** — the `code-review` skill or a
   freshly spawned subagent that hasn't been part of writing the change. Address its findings
   before considering the work done.
7. Commit messages end with the repo's required co-author + session trailers (the harness
   supplies these).
