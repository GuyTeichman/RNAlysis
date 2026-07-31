# design-philosophy.md — the "why" behind RNAlysis

The principles that should guide every change. When a decision is unclear, resolve it in favor
of these. Companion to [`../CLAUDE.md`](../CLAUDE.md), [`context.md`](context.md), and
[`workflows.md`](workflows.md).

---

## Who we build for

**The target user is a biologist with zero programming experience**, often not a
bioinformatician either. That single fact drives most of what follows: the GUI exists so these
users never touch code; defaults exist so they get correct results without understanding every
knob; transparency exists so they can *report* exactly what their analysis did.

---

## Foundational principles

### 1. Correctness and reliability are paramount
This software is used to make **biological discoveries**. Wrong output is worse than no output.

- **Don't ship a feature you're not sure works.** Confidence comes from tests (TDD) and from
  understanding the method, not from "it ran without error."
- **Results must not change arbitrarily between versions.** A given analysis with given
  parameters should reproduce across releases. Any intentional change to numerical output must
  be justified and documented in `HISTORY.rst`.

### 2. Reproducibility is a promise, not a feature
Pipelines, saved sessions, exported parameters, and the interactive analysis report all exist so
that an analysis can be re-run, audited, and reported exactly.

- **Backward-compatible loading is a hard invariant:** old Pipeline YAMLs, sessions, and
  exported parameter files must keep loading in new versions.
- The auto-generated report captures the full analysis path for transparency and publication.

### 3. Reasonable defaults for everything
Every function parameter should have a default that yields a sensible, correct analysis for a
non-expert. Advanced users can override; beginners shouldn't have to. This is the intended
standard for all past and future functions.

### 4. Universality / cross-system compatibility
RNAlysis should "just work" on the machines scientists actually use.

- Supported everywhere: **Windows, macOS, and Linux** (the last often headless **clusters**).
- Two run modes that must both work: **from source** (PyPI) and **frozen** (PyInstaller
  standalone). The frozen environment is why multiprocessing and packaging are handled carefully
  (`FROZEN_ENV`, restricted `PARALLEL_BACKENDS`).
- A user should ideally never have to ask "is this feature supported on my OS?"

### 5. Graphical *and* programmatic parity
The **majority of features must work both in the GUI and in Python scripts.** The API is the
source of truth and the GUI is generated from it — so parity is largely structural, but keep it
that way: don't build GUI-only capabilities that scripts can't reach, and don't add API behavior
that the reflection layer can't express.

### 6. Customization *with* clarity → transparency
Give users a high degree of control over parameters — and document every one — for two reasons:
advanced users can tune analyses, and, more importantly, **fully specified + documented
parameters make analyses transparent and accurately reportable**. A user should be able to say
exactly how GO annotations were fetched, filtered, and tested, which statistic was used, and
what the defaults were.

---

## How the principles show up in the code (so you don't fight them)

- **Type-annotations-as-UI.** Precise `param_typing` annotations are how the GUI is designed and
  how parameter transparency is enforced. Loose annotations = a worse GUI and murkier reporting.
- **Docstrings-as-help.** `:param:` docstrings surface directly to users; they're part of the
  product, not just developer notes.
- **`@readable_name`** turns API methods into human-facing GUI labels — the same operation has
  one name for scripts and a friendly name for the GUI.
- **Deliberate parallelism/packaging choices** exist to satisfy universality across frozen and
  source environments; improve them without regressing either.

---

## Decision heuristics

- Torn between a clever refactor and a boring, obviously-correct one near analysis output? Pick
  correctness and reproducibility.
- Tempted to change a default? Only if the new default is *more* reasonable for a non-expert and
  doesn't silently change existing users' results.
- Considering an OS- or environment-specific shortcut? Don't — unless it's guarded and the other
  environments still work.
- Adding a parameter? Give it a sensible default, a precise type annotation, and a docstring
  line. All three.
