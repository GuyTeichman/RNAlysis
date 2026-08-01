# context.md — RNAlysis project context

Undiscoverable-from-the-code context: who this is for, what state the project is in, what's
fragile, and where it's headed. Pairs with [`../CLAUDE.md`](../CLAUDE.md),
[`workflows.md`](workflows.md), and [`design-philosophy.md`](design-philosophy.md).

---

## Who & what

- **Author/maintainer:** Guy Teichman (guyteichman@gmail.com). Open-source, MIT-licensed,
  volunteer-driven, with several named contributors (see `README.rst`).
- **Users:** biologists doing real RNA-seq analysis for **biological discovery** — largely
  **non-programmers** and often **non-bioinformaticians**. Many run it on lab Windows/Mac
  machines; power users run pieces on **Linux clusters**.
- **Citation:** Teichman et al. (2023), *BMC Biology* 21:74. A published scientific tool —
  wrong results are worse than missing features.
- **Distribution:** PyPI (`pip install RNAlysis[all]`, launch `rnalysis-gui`) and standalone
  PyInstaller apps (`RNAlysis.exe` / `RNAlysis.dmg`) for users who can't install Python.
- **Current version:** 4.2.0 (see `HISTORY.rst`). Recent arc: 4.0 = Pandas→Polars + report
  overhaul; 4.1 = Qt5→Qt6 + NumPy 2; 4.2 = Py 3.14 support + report perf + external-API fixes.

## How the pieces fit (quick recall)

Three public API modules — `filtering` (tables), `enrichment` (gene sets), `fastq` (reads) —
built on a `utils/` implementation layer, with a PyQt6 `gui/` that **reflects the API** rather
than reimplementing it. External muscle: **R** (DESeq2, limma-voom, Rsubread featureCounts) via
templated scripts, **CLI tools** (kallisto, bowtie2, ShortStack, CutAdapt, Picard/samtools) via
subprocess, and **web databases** (UniProt, Ensembl, PANTHER, PhylomeDB, OrthoInspector, KEGG,
GO) via async HTTP. Cross-cutting features: **Pipelines** (record funcs+params → YAML, replay),
**sessions** (save/restore whole analysis state), and **interactive reports** (a networkx/pyvis
graph of the whole analysis path, for reproducibility).

---

## Fragile spots — tread carefully

1. **External web APIs (most fragile).** They change response formats, rename fields, and go
   down without notice. The **current open PR** is fixing changes in the API/outputs of **≥2**
   external services. When coding against one, verify its *current real* behavior (use the
   `research` skill / a quick live probe) rather than trusting old code or memory. Keep the
   retry/rate-limit/cache scaffolding; degrade gracefully; never let a dead service crash the app.

2. **Multiprocessing across the frozen-vs-source split.** The app runs both from a terminal and
   as a frozen PyInstaller executable, and multiprocessing behaves differently between them.
   `FROZEN_ENV` gates `PARALLEL_BACKENDS` (frozen → `multiprocessing`/`sequential`; source →
   also `loky`/`threading`). The current design reflects hard-won choices — there's room to
   improve, but **don't break it**; reason about both environments before touching parallelism.

3. **Tests, especially GUI/qtbot tests.** `pytest-qt`/`qtbot` produces **CI-only crashes that
   never reproduce locally**; the root cause isn't pinned down yet. The suite is also **very
   long**. Known desired improvement: split into **unit / integration / e2e** via pytest markers
   and run sequentially with **fail-fast**, and/or make slow tests faster. Until then, a green
   local run is not proof of green CI.

4. **R interface / environment setup.** On test machines and fresh installs, installing or
   updating R packages frequently fails on **write permissions** to the R library folder. A
   recognized sore spot to harden. `installs.py` already emits actionable messages here.

## Stable spots

- **The GUI framework** and the **Polars migration** are solid. Don't churn them without reason.
- Caveat: Polars usage is *functional but dated* — the migration was Guy's first Polars project
  years ago; both his skill and Polars itself have advanced since. Notably **lazy evaluation is
  barely used** and could improve performance/UX. Modernizing is welcome (see roadmap).

---

## Near-term roadmap (what Guy most wants help with)

Priority themes, not a committed schedule:

1. **External-API resilience** — detect/handle upstream changes and outages; strengthen
   retries, caching, and graceful degradation. (Overlaps the current open PR.)
2. **Test-suite restructuring** — unit/integration/e2e markers, fail-fast ordering, and
   stabilizing + speeding up the slow/flaky GUI tests.
3. **R setup robustness** — make install/update/permissions painless on fresh machines and CI.
4. **Polars modernization** — introduce lazy evaluation and newer idioms for performance/UX.
5. **Biological file parsing (GTF/GFF/etc.)** — parse these **directly with Polars**: robust to
   real-world format variation, faster, and better integrated into RNAlysis. **Guy has existing
   code elsewhere to seed this** — ask for it before starting. Ties themes 1 and 4 together.

---

## Working agreement (how Guy wants me to operate)

- **TDD** (red-green-refactor); relevant tests run; new code covered.
- **Plan-first** when a change is risky (public API signature, the API→GUI reflection contract,
  or a serialized format). Otherwise implement and show the diff.
- **Serialized-format back-compat and result reproducibility are hard invariants.**
- **Branch model:** feature branch off `development` → PR into `development` → `master` at
  release. Never commit directly to `master` or `development`.
- **A clean-context independent agent reviews the diff at PR time.**
- Ship only features we're confident are correct. Provide reasonable defaults. Keep it working
  on Windows, macOS, and Linux, both from source and frozen.
