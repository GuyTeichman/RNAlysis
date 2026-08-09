---
name: safe-optimization
description: >-
  Make RNAlysis code faster while proving the output does not change. Use before touching any
  code for performance reasons -- a slow function, a profiling request, "speed this up", "this
  is too slow on large datasets", vectorizing a loop, parallelizing across CPU cores, adding/
  changing a numba @jit, or swapping an algorithm/data structure for a faster one. Also use when
  asked to "benchmark", "profile", or write a perf PR. Skip for changes that are not primarily
  about speed (a correctness fix, a new feature) even if they happen to also be faster.
---

# Safe performance optimization

Correctness and reproducibility are non-negotiable in RNAlysis (`CLAUDE.md` rule #5, `AGENTS.md`
rule #5): *"Results must reproduce across versions. A given analysis with given parameters must
produce the same output."* A performance change that is fast but wrong, or fast but silently
different, is a regression dressed up as an improvement. This skill is the procedure for making
code faster **and** proving it did not change what the code computes.

The engine for step 4 (proof) is a plain, provider-neutral script: **`packaging/bench_equal.py`**.
This skill adds the *when* and the *procedure* around it.

---

## When to use / skip

Use this skill for any change whose primary motivation is speed: a slow function, a profiling
request, vectorizing a loop, parallelizing across CPU cores, adding/changing a `numba` `@jit`,
swapping an algorithm or data structure, introducing Polars lazy evaluation to a hot path, etc.

Skip it for changes that are not primarily about speed -- a correctness fix, a new feature, a
refactor for readability -- even if they incidentally run faster. Those still go through the
normal `tdd` workflow; retrofitting this skill's ceremony onto them is not the point. (If a
correctness fix happens to also touch a real hotspot, it is fine to fold in a benchmark, but the
red-green test for the *bug* comes first.)

---

## Step 1 -- Profile first. Do not guess.

Optimizing the wrong thing wastes effort and adds risk to the codebase for nothing. This
codebase's own history has punished guessing more than once:

- In `pca`/clustering, **Box-Cox** (`rnalysis/utils/generic.py::box_cox`), not PCA itself, was
  the dominant per-gene cost -- the obvious suspect (PCA) was not the real one.
- In GO/KEGG enrichment, **set-intersection work was not a hotspot at all**; the real cost was
  the stats-test computation (recomputed per ontology term) and an `elim`-propagation deep-copy
  repeated on every run (`rnalysis/utils/enrichment_runner.py`) -- see `HISTORY.rst`'s entry on
  memoizing the hypergeometric p-value and dropping the `elim` deep-copy.

Profile with whatever tool fits the question -- `cProfile`/`pstats` for a first cut ("which
function"), `line_profiler` for "which line", `py-spy`/`scalene` for sampling a running process
without instrumentation overhead. A couple of representative inputs at realistic scale (not a
toy 3-row table) is what surfaces real hotspots; a tiny input mostly measures constant overhead.
Only once profiling data names a specific function/line should you move to step 2.

## Step 2 -- Establish a baseline

Before changing anything, capture the *current* code's output on one or more representative
inputs. "Representative" means inputs that exercise the real shape of the data this function
sees in production: realistic row/column counts, some `NaN`/null values if the real data has
them, edge cases (empty input, a single row, all-identical values) if those are plausible.

In practice this is either:
- a fixture already in `tests/test_files/`, loaded the way the matching `tests/test_*.py`
  module loads it, or
- synthetic data built to match the real shape (`np.random.default_rng(<fixed seed>)` /
  `pl.DataFrame(...)`) when no fixture is representative enough.

You do not need to persist the baseline anywhere durable -- `bench_equal.compare()` (step 4)
calls the *current* code itself as the baseline and compares it against your optimized version
in the same run, so "capturing" it is really just picking the inputs and keeping a reference to
the pre-optimization function (e.g. via `git stash`/a second import/a renamed copy) long enough
to run the comparison.

## Step 3 -- Optimize

Make the change the profiling data justified. Common levers in this codebase: vectorize a
Python loop with NumPy/Polars, push work into a Polars lazy pipeline instead of eager, memoize
a value recomputed across many iterations, parallelize across CPU cores (mind the
frozen-vs-source gotcha below), or add a `numba` `@jit` to a numeric inner loop (mind the RNG
gotcha below).

## Step 4 -- Prove bit-identical output

This is the step that turns "I made it faster" into "I made it faster, safely." Use
`packaging/bench_equal.py`'s `assert_equal`/`compare` to compare the old and new implementations
on the **same** inputs from step 2:

```python
import importlib.util
import sys
from pathlib import Path

spec = importlib.util.spec_from_file_location('bench_equal', Path('packaging/bench_equal.py'))
bench_equal = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = bench_equal
spec.loader.exec_module(bench_equal)

result = bench_equal.compare(old_box_cox, new_box_cox, args=(representative_data,))
print(f'{result.speedup:.1f}x faster, bit-identical')
```

or, for several representative inputs at once, write a small throwaway "bench spec" file and run
it from the CLI (`python packaging/bench_equal.py my_bench_spec.py --repeats 5`) -- see
`bench_equal.py`'s module docstring for the spec-file format (`BASELINE`, `CANDIDATE`, `INPUTS`).

**Exact equality is the default and should stay the default.** `assert_equal`/`compare` compare
NumPy arrays and Polars `DataFrame`/`Series` element-for-element (shape and dtype must match
too), with `NaN`/null treated as equal to itself -- the same standard the test suite already
holds analysis output to.

Only pass `rtol`/`atol` when the optimization *itself* legitimately reorders floating-point
operations (a parallel reduction, a different BLAS/LAPACK code path, `float32` vs `float64`
accumulation) such that bit-identical output is not achievable even though the computation is
equivalent. `HISTORY.rst` already has a real precedent for this: the parallelized Box-Cox
transform is documented as *"Results are unchanged up to floating-point precision (verified
against the test suite's reference outputs)"* -- not "unchanged", because parallel reduction
order isn't associative in floating point. If you need a tolerance, that is itself the
reproducibility event from hard invariant #5: **it must be intentional, justified, and called
out in `HISTORY.rst`** (state the tolerance and why exact equality isn't achievable), not a
quiet flag you flip to make a red assertion go green.

If `compare()`/`assert_equal` raises `AssertionError`: the optimization is not safe yet. Either
the optimization has a bug, or it is a genuine behavior change that needs its own justification
and sign-off (per the plan-first rule for risky changes) -- not a benchmark to report.

## Step 5 -- Benchmark and record the speedup

Once equality holds, report the number. `BenchmarkResult.speedup` (from `compare()`) gives it to
you directly, using best-of-N wall-clock time on each side (robust to one noisy call). Quote it
in the PR description and, for a user-visible perf change, in `HISTORY.rst` -- concretely, e.g.:
*"the hypergeometric test and `elim` propagation benefit the most (up to ~20x and ~2-9x faster
respectively on large ontologies)"*. Prefer a range across the representative inputs from step 2
over a single cherry-picked best case.

## Step 6 -- Known gotchas to check

**`numba` RNG must be seeded *inside* the jitted function.** `numba` compiles its own internal
RNG state that is separate from NumPy's global RNG -- seeding NumPy's global RNG from ordinary
Python code before calling into `@jit(nopython=True)` code is a **no-op** for any `np.random.*`
call made *inside* that jitted function. This codebase has a live example of the trap:
`rnalysis/utils/enrichment_runner.py`'s `PermutationTest.run()` calls `np.random.seed(...)` in
plain Python, then calls the jitted `_calc_permutation_pval` (decorated
`@generic.numba.jit(nopython=True)`), which itself calls `np.random.choice(...)` -- the seed set
in `run()` does not reach the random draws inside `_calc_permutation_pval`. This was the root
cause of a previously flaky "reproducibility" test. If you jit a function that needs
reproducible randomness, seed it with `np.random.seed(...)` (or pass the seed in and call it)
**from inside** the jitted function itself, and prove it with a test that calls the jitted
function twice with the same seed and asserts equal output.

**The frozen-vs-source multiprocessing split.** The app runs both from a source checkout and as
a frozen PyInstaller executable (`RNAlysis.exe`/`.dmg`), and they do not support the same
parallel backends. `rnalysis/__init__.py` sets `FROZEN_ENV = getattr(sys, 'frozen', False) and
hasattr(sys, '_MEIPASS')`; `rnalysis/utils/param_typing.py` gates on it directly:

```python
PARALLEL_BACKENDS = ('multiprocessing', 'sequential') if FROZEN_ENV else (
    'multiprocessing', 'loky', 'threading', 'sequential')
```

`loky`/`threading` are only available from source -- a frozen build cannot use them. If your
optimization parallelizes something, either accept a `parallel_backend` parameter typed
`Literal[PARALLEL_BACKENDS]` (as the existing `filtering`/`enrichment` functions do) so the GUI
exposes only the legal choices per environment, or, if the backend is chosen internally rather
than user-facing, follow the pattern in `rnalysis/utils/generic.py::box_cox_parallel_backend()`
(picks `'multiprocessing'` when `FROZEN_ENV`, `'loky'` otherwise) rather than hardcoding a
backend that breaks one of the two shipping forms. Reason through *both* environments before
touching parallelism -- you cannot test the frozen build's behavior by running from source.

---

## Definition of done

A performance change is not done until, in the PR:

1. **Profiling evidence** names the actual hotspot (not an assumption) -- a one-line summary or
   a `cProfile`/`line_profiler` excerpt is enough.
2. **The baseline and candidate were compared on representative inputs** via
   `bench_equal.compare()`/`assert_equal` -- exact by default; any `rtol`/`atol` used is
   justified in the PR text.
3. **A measured speedup number** is reported (best-of-N, from `BenchmarkResult.speedup` or the
   CLI's printed ratio).
4. If output changed at all (even "only" up to floating-point precision), **`HISTORY.rst` has an
   entry** stating that plainly -- per hard invariant #5, this is never a silent side effect.
5. The relevant `tests/test_*.py` module(s) still pass, per the normal `tdd`/finishing-a-change
   workflow in `.claude/workflows.md`.
6. If the optimization touches parallelism, the frozen-vs-source split (Step 6) was reasoned
   through, not ignored.
