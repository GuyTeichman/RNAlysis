# AGENTS.md — working in RNAlysis with an AI agent

Provider-neutral entry point for **any** AI coding agent or automated contributor
(Claude, Cursor, Aider, GitHub Copilot, Codex, Windsurf, Gemini, custom bots, …) working in
this repository. It follows the [agents.md](https://agents.md) convention.

**The full working guide is [`CLAUDE.md`](CLAUDE.md).** Despite the name it is written for
*all* agents, not only Claude Code — it holds the complete repository map, the API↔GUI
reflection contract, the command list, the style rules, and the known gotchas. Read it. Deeper
material lives in [`.claude/context.md`](.claude/context.md),
[`.claude/workflows.md`](.claude/workflows.md), and
[`.claude/design-philosophy.md`](.claude/design-philosophy.md). Human contributors should also
read [`CONTRIBUTING.rst`](CONTRIBUTING.rst).

What follows is the short version: the one architectural idea you must understand before
changing anything, and the non-negotiable rules. If you read nothing else here, read this — but
prefer reading `CLAUDE.md` in full.

---

## What RNAlysis is

Desktop software that lets biologists analyze RNA-sequencing data — from raw FASTQ through
filtering, normalization, differential expression, clustering, and enrichment — **without
writing any code**. It ships a PyQt6 graphical app and a programmatic Python API, and targets
scientists with **zero programming experience**. Correctness and reproducibility are
non-negotiable.

---

## The one idea to understand first: the API↔GUI reflection contract

**The programmatic Python API is the source of truth. The GUI is generated from it by
reflection.** When you add a public method to a `Filter`/`FeatureSet` subclass (or a public
function in `fastq`/`enrichment`), the GUI **auto-discovers and renders it** — you do not write
Qt code for it. What drives that rendering:

- **Type annotations are load-bearing UI.** Each parameter's annotation maps to a specific Qt
  widget (e.g. `bool` → toggle, `Fraction` → 0–1 slider, `Color` → color picker,
  `Literal[...]` → combo box, `Optional[T]` → a checkbox that enables a `T` widget). The custom
  types live in `rnalysis/utils/param_typing.py`. **Choosing the right annotation is how you
  design the GUI.**
- **Docstrings are load-bearing help text.** The reST `:param x: ...` lines become the
  per-parameter tooltips; keep that format.
- **`@readable_name('…')`** sets the human-readable button/window label.

**Consequence:** renaming a public method, changing its signature, or loosening a type
annotation silently changes the GUI (and possibly serialized parameters). Treat any such change
as *risky*. The full contract — discovery heuristics, the complete annotation→widget table, and
how to hide internals — is in [`CLAUDE.md`](CLAUDE.md).

---

## Non-negotiable rules

Hard invariants. Do not violate them without explicit sign-off from a maintainer.

1. **Test-driven development, always.** Write a failing test first, make it pass, then refactor.
   New code must be covered. Run the *relevant* `tests/test_*.py` module(s) before calling
   anything done, and state plainly what you could **not** run (paths needing R, network,
   kallisto/bowtie2, or a Qt display).
2. **Branch model: feature branch → `development` → `master`.** `development` is the long-lived
   integration branch; `master` is stable/released and only receives `development` at a version
   release. **Never commit directly to `master` or `development`.** Cut a feature/fix branch off
   `development` and open your PR **into `development`**.
3. **Plan first for risky changes.** If a change touches a public API signature, the API→GUI
   reflection contract, or a serialized format (Pipeline YAML / session file / exported
   parameters), propose a short plan and wait for approval before coding. Small localized fixes:
   just do it (TDD) and show the diff.
4. **Back-compat of serialized artifacts is a hard invariant.** Old Pipeline YAMLs, saved
   sessions, and exported parameter files **must keep loading** in new versions. Breaking this
   is a bug — flag it loudly.
5. **Results must reproduce across versions.** A given analysis with given parameters must
   produce the same output. If a change alters numerical output, it must be intentional,
   justified, and called out in `HISTORY.rst`.
6. **Every new function needs reasonable defaults.** Assume a non-programmer, non-bioinformatician
   user; defaults should produce a sensible, correct analysis out of the box.
7. **Cross-platform or it doesn't ship.** Must work on Windows, macOS, and Linux, both
   run-from-source **and** frozen (PyInstaller). Mind the multiprocessing/frozen split described
   in `CLAUDE.md`.
8. **Independent review before done.** At PR time, have an agent or reviewer with clean context
   — not the one that wrote the change — review the diff, and address the findings.
9. **Visible GUI changes ship with screenshots.** Because the GUI is reflected from the API, a
   change to a public `Filter`/`FeatureSet`/`fastq`/`enrichment` function that is **visible in a
   dialog** (a new/renamed/removed parameter, a changed type annotation, a new/edited
   `@readable_name`, or edited `:param:` help text) must include a screenshot of the affected
   dialog(s) on the PR — before/after for a renamed or retyped parameter. Generate them from the
   API headlessly with `packaging/capture_gui_dialog.py` (run `--help` for usage); publish the PNGs
   however your workflow links artifacts. Behind-the-scenes changes with no visible dialog effect
   are exempt.

---

## Fast start

```bash
pip install -e .[all]                      # install for development
pip install -r requirements_dev.txt        # pytest, pytest-qt, etc.

rnalysis-gui                               # run the app (installed entry point)

pytest tests/test_filtering.py             # run one module while iterating
pytest tests/test_filtering.py -k <name>   # run a subset of a module
```

The full command list, the CI matrix, the repository map, the style rules, and the
external-API / multiprocessing / Qt-flakiness gotchas are all in [`CLAUDE.md`](CLAUDE.md). A few
things worth stating up front: don't hand-edit generated files (`docs/source/rnalysis.*.rst`
come from `sphinx-apidoc`), and prefer **Polars** over Pandas for tabular work (the codebase
migrated off Pandas in 4.0).

---

## A note on tool-specific instructions

`CLAUDE.md` refers to Claude Code "skills" (`tdd`, `diagnosing-bugs`, `code-review`,
`research`, `safe-optimization`). Those are convenience shortcuts for one particular tool — but
the *workflows they encode* apply to every agent: disciplined TDD, methodical debugging, a
clean-context review of each diff, verifying an external API's real current behavior before
coding against it, and profiling-then-proving a performance change is safe (rule 5;
`packaging/bench_equal.py` is the provider-neutral equality/benchmark engine). Use whatever
tooling you have to achieve the same outcomes.
`research`, `gui-screenshots`). Those are convenience shortcuts for one particular tool — but the
*workflows they encode* apply to every agent: disciplined TDD, methodical debugging, a
clean-context review of each diff, verifying an external API's real current behavior before coding
against it, and attaching a dialog screenshot to any visible GUI change (rule 9;
`packaging/capture_gui_dialog.py` is the provider-neutral engine). Use whatever tooling you have to
achieve the same outcomes.
