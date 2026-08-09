#!/usr/bin/env python
"""Benchmark two implementations of a function and ASSERT their outputs are equal.

This is the harness half of the ``safe-optimization`` skill (``.claude/skills/safe-optimization/
SKILL.md``): "faster" only counts once you can prove the output did not change. The core idea,
independent of any particular optimization:

    1. capture a baseline (the current/"old" implementation) on representative inputs,
    2. run the candidate ("new"/optimized) implementation on the *same* inputs,
    3. assert the two outputs are equal -- exact by default, with an opt-in float tolerance,
    4. report how much faster the candidate is.

Equality is exact (bit-for-bit on the values, NaN/null treated as matching itself) unless you
pass ``rtol``/``atol``, in which case numeric comparisons use them (still exact for
non-numeric data -- a tolerance never loosens a string/bool/dtype/shape mismatch). Per the
project's hard invariant that results must reproduce across versions (``CLAUDE.md`` rule #5),
exact is the right default: reach for a tolerance only when the optimization *itself*
legitimately reorders floating-point operations (parallel reductions, a different BLAS path,
etc.), and when you do, that is a reproducibility event worth a line in ``HISTORY.rst``, not a
detail to wave away.

Provider-neutral: no dependency on Claude Code, RNAlysis, or any AI tooling. It understands
plain Python scalars/lists/tuples/dicts, NumPy arrays, and Polars Series/DataFrames -- the data
shapes this codebase's hot paths actually pass around.

Library usage (typical -- import from a throwaway benchmarking script)
------------------------------------------------------------------------
    from bench_equal import compare  # loaded by path, see the skill for the snippet

    result = compare(old_box_cox, new_box_cox, args=(data,), rtol=1e-10, atol=1e-12)
    print(f"{result.speedup:.1f}x faster, bit-identical up to the given tolerance")

``compare()`` raises ``AssertionError`` if the outputs differ -- that IS the proof step; let it
raise (test-fail) if the optimization is not actually safe.

CLI usage (batch over several representative inputs, from a small "bench spec" file)
-----------------------------------------------------------------------------------------
    python packaging/bench_equal.py my_bench_spec.py --repeats 5

where ``my_bench_spec.py`` is a plain Python file defining:

    BASELINE  = <callable>          # required -- the current/"old" implementation
    CANDIDATE = <callable>          # required -- the optimized/"new" implementation
    INPUTS    = [...]               # optional -- see below (default: one no-arg call)
    # -- or, if inputs are expensive to build --
    def make_inputs(): ...          # optional -- called once, must return the same shape as INPUTS

Each item of ``INPUTS`` (or the list returned by ``make_inputs()``) is either:
  * a ``list``/``tuple`` of positional arguments, e.g. ``(data, 5)``; or
  * a ``dict`` with optional ``'args'`` (tuple) and ``'kwargs'`` (dict) keys, for functions that
    need keyword arguments, e.g. ``{'args': (data,), 'kwargs': {'axis': 0}}``.

Exit code is the number of inputs that failed to match (0 == every input was safe).
"""
from __future__ import annotations

import argparse
import importlib.util
import math
import statistics
import sys
import time
from collections.abc import Mapping
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Callable, Optional, Sequence

import numpy as np
import polars as pl

_NUMBER_TYPES = (bool, int, float, np.bool_, np.integer, np.floating)


# ============================================================================================
# Timing
# ============================================================================================

@dataclass(frozen=True)
class Timing:
    """Wall-clock durations (seconds) collected over one or more repeats of a single call."""
    times: tuple

    @property
    def best(self) -> float:
        """The minimum observed duration -- the standard microbenchmark statistic (what
        ``timeit.repeat(...)`` + ``min()`` reports), because it is the call least distorted by
        scheduler/GC noise. Prefer this over the mean when reporting a speedup."""
        return min(self.times)

    @property
    def mean(self) -> float:
        return statistics.mean(self.times)

    @property
    def median(self) -> float:
        return statistics.median(self.times)


def time_call(func: Callable, args: Sequence = (), kwargs: Optional[Mapping] = None, *,
              repeats: int = 5, warmup: int = 1):
    """Call ``func(*args, **kwargs)`` ``repeats`` times, timing each call.

    ``warmup`` untimed calls run first (JIT warm-up for numba, import caching, first-call
    allocator overhead, etc. -- excluded from the timed statistics on purpose).

    Returns ``(result, timing)``: ``result`` is the return value of the *last* timed call, and
    ``timing`` is a :class:`Timing` holding every individual timed call's duration. This assumes
    ``func`` is deterministic/pure for the given input -- exactly the premise the rest of this
    harness (and the ``safe-optimization`` skill) is built on; a function whose output legitimately
    varies call-to-call (uncontrolled randomness, wall-clock timestamps, ...) is not a fit for
    this comparison and should be made deterministic (e.g. an explicit seed) before benchmarking.
    """
    if repeats < 1:
        raise ValueError('repeats must be >= 1')
    kwargs = dict(kwargs or {})
    for _ in range(max(warmup, 0)):
        func(*args, **kwargs)
    times = []
    result = None
    for _ in range(repeats):
        start = time.perf_counter()
        result = func(*args, **kwargs)
        times.append(time.perf_counter() - start)
    return result, Timing(tuple(times))


def speedup(baseline: Timing, candidate: Timing) -> float:
    """How many times faster ``candidate`` is than ``baseline``, using each side's best-of-N time.

    Returns ``float('inf')`` if the candidate's best time is zero (immeasurably fast) rather than
    raising a ``ZeroDivisionError`` -- a benchmark script should be able to print that result, not
    crash on it.
    """
    if candidate.best <= 0:
        return math.inf
    return baseline.best / candidate.best


# ============================================================================================
# Equality
# ============================================================================================

def assert_equal(actual: Any, expected: Any, *, rtol: float = 0.0, atol: float = 0.0,
                  _path: str = 'root') -> None:
    """Recursively assert ``actual == expected``, raising ``AssertionError`` with a precise
    location on the first mismatch.

    Exact by default (``rtol=atol=0.0``). Handles, recursively:
      * NumPy arrays -- shape and dtype must match; values compared with ``rtol``/``atol``
        (NaNs treated as equal to each other).
      * Polars ``DataFrame``/``Series`` -- columns/length/dtype must match; numeric columns
        compared with ``rtol``/``atol`` (nulls/NaNs treated as equal to each other), non-numeric
        columns compared exactly regardless of tolerance.
      * ``dict``/``Mapping`` -- same keys, then each value compared recursively.
      * ``list``/``tuple`` -- same length, then each element compared recursively (order matters).
      * numbers (``int``/``float``/NumPy scalars, including ``bool``) -- compared with
        ``rtol``/``atol``; NaN equals NaN.
      * anything else -- plain ``==``.

    A dtype/shape/column mismatch is always a failure, even with a tolerance set: a tolerance
    loosens *numeric* comparison, it never papers over a structural difference (that is usually
    the more interesting bug). If an optimization legitimately changes a dtype (e.g. an upcast),
    cast explicitly before comparing and call it out in ``HISTORY.rst``.
    """
    if isinstance(actual, np.ndarray) or isinstance(expected, np.ndarray):
        _assert_ndarray_equal(np.asarray(actual), np.asarray(expected), rtol, atol, _path)
        return

    if isinstance(actual, pl.DataFrame) or isinstance(expected, pl.DataFrame):
        _assert_dataframe_equal(actual, expected, rtol, atol, _path)
        return

    if isinstance(actual, pl.Series) or isinstance(expected, pl.Series):
        actual_s = actual if isinstance(actual, pl.Series) else pl.Series(actual)
        expected_s = expected if isinstance(expected, pl.Series) else pl.Series(expected)
        _assert_series_equal(actual_s, expected_s, rtol, atol, _path)
        return

    if isinstance(actual, Mapping) and isinstance(expected, Mapping):
        if actual.keys() != expected.keys():
            raise AssertionError(f'dict keys differ at {_path}: {sorted(actual.keys())!r} != '
                                  f'{sorted(expected.keys())!r}')
        for key in actual:
            assert_equal(actual[key], expected[key], rtol=rtol, atol=atol, _path=f'{_path}[{key!r}]')
        return

    if isinstance(actual, (list, tuple)) and isinstance(expected, (list, tuple)):
        if len(actual) != len(expected):
            raise AssertionError(f'length differs at {_path}: {len(actual)} != {len(expected)}')
        for i, (a_item, e_item) in enumerate(zip(actual, expected)):
            assert_equal(a_item, e_item, rtol=rtol, atol=atol, _path=f'{_path}[{i}]')
        return

    if isinstance(actual, _NUMBER_TYPES) and isinstance(expected, _NUMBER_TYPES):
        _assert_numbers_equal(float(actual), float(expected), rtol, atol, _path)
        return

    if actual != expected:
        raise AssertionError(f'values differ at {_path}: {actual!r} != {expected!r}')


def equal(actual: Any, expected: Any, *, rtol: float = 0.0, atol: float = 0.0) -> bool:
    """Non-raising wrapper around :func:`assert_equal` -- returns ``True``/``False``."""
    try:
        assert_equal(actual, expected, rtol=rtol, atol=atol)
        return True
    except AssertionError:
        return False


def _assert_numbers_equal(a: float, b: float, rtol: float, atol: float, path: str) -> None:
    if math.isnan(a) and math.isnan(b):
        return
    if not math.isclose(a, b, rel_tol=rtol, abs_tol=atol):
        raise AssertionError(f'numbers differ at {path} (rtol={rtol}, atol={atol}): {a!r} != {b!r}')


def _assert_ndarray_equal(a: np.ndarray, b: np.ndarray, rtol: float, atol: float, path: str) -> None:
    if a.shape != b.shape:
        raise AssertionError(f'array shape differs at {path}: {a.shape} != {b.shape}')
    if a.dtype != b.dtype:
        raise AssertionError(f'array dtype differs at {path}: {a.dtype} != {b.dtype} '
                              f'(cast explicitly if this is expected)')
    try:
        matches = np.allclose(a, b, rtol=rtol, atol=atol, equal_nan=True)
    except TypeError:
        # non-numeric dtype (object/str/...): equal_nan is meaningless there, fall back to exact
        matches = bool(np.array_equal(a, b))
    if not matches:
        max_diff = _safe_max_abs_diff(a, b)
        raise AssertionError(f'array values differ at {path} beyond tolerance '
                              f'(rtol={rtol}, atol={atol}); max abs diff={max_diff}')


def _safe_max_abs_diff(a: np.ndarray, b: np.ndarray):
    try:
        return np.nanmax(np.abs(a.astype(np.float64) - b.astype(np.float64)))
    except (TypeError, ValueError):
        return 'n/a (non-numeric dtype)'


def _assert_dataframe_equal(actual: Any, expected: Any, rtol: float, atol: float, path: str) -> None:
    if not isinstance(actual, pl.DataFrame) or not isinstance(expected, pl.DataFrame):
        raise AssertionError(f'type mismatch at {path}: {type(actual).__name__} vs {type(expected).__name__}')
    if actual.columns != expected.columns:
        raise AssertionError(f'columns differ at {path}: {actual.columns!r} != {expected.columns!r}')
    if actual.height != expected.height:
        raise AssertionError(f'row count differs at {path}: {actual.height} != {expected.height}')
    for col in actual.columns:
        _assert_series_equal(actual[col], expected[col], rtol, atol, f'{path}.{col}')


def _assert_series_equal(a: pl.Series, b: pl.Series, rtol: float, atol: float, path: str) -> None:
    if a.len() != b.len():
        raise AssertionError(f'series length differs at {path}: {a.len()} != {b.len()}')
    if a.dtype != b.dtype:
        raise AssertionError(f'series dtype differs at {path}: {a.dtype} != {b.dtype} '
                              f'(cast explicitly if this is expected)')
    if (rtol or atol) and a.dtype.is_numeric():
        if not np.allclose(a.to_numpy(), b.to_numpy(), rtol=rtol, atol=atol, equal_nan=True):
            raise AssertionError(f'series values differ at {path} beyond tolerance (rtol={rtol}, atol={atol})')
        return
    if not a.equals(b, check_dtypes=True, check_names=False, null_equal=True):
        raise AssertionError(f'series are not exactly equal at {path}')


# ============================================================================================
# compare(): the timing + equality combinator
# ============================================================================================

@dataclass(frozen=True)
class BenchmarkResult:
    result: Any
    baseline_timing: Timing
    candidate_timing: Timing

    @property
    def speedup(self) -> float:
        return speedup(self.baseline_timing, self.candidate_timing)


def compare(baseline: Callable, candidate: Callable, args: Sequence = (),
            kwargs: Optional[Mapping] = None, *, repeats: int = 5, warmup: int = 1,
            rtol: float = 0.0, atol: float = 0.0) -> BenchmarkResult:
    """Time both implementations on the same input and ASSERT their outputs are equal.

    Raises ``AssertionError`` (via :func:`assert_equal`) if the candidate's output does not
    match the baseline's -- that assertion failing means the optimization is not (yet) safe to
    ship, full stop. On success, returns a :class:`BenchmarkResult` with both timings and the
    (now known-equal) result, so a caller can report the speedup.
    """
    baseline_result, baseline_timing = time_call(baseline, args, kwargs, repeats=repeats, warmup=warmup)
    candidate_result, candidate_timing = time_call(candidate, args, kwargs, repeats=repeats, warmup=warmup)
    assert_equal(candidate_result, baseline_result, rtol=rtol, atol=atol, _path='result')
    return BenchmarkResult(baseline_result, baseline_timing, candidate_timing)


# ============================================================================================
# CLI
# ============================================================================================

def _load_spec(path: Path):
    spec = importlib.util.spec_from_file_location(path.stem, path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module  # let the spec file use dataclasses/relative lookups safely
    spec.loader.exec_module(module)
    return module


def _resolve_raw_inputs(module) -> list:
    if hasattr(module, 'INPUTS'):
        return list(module.INPUTS)
    if hasattr(module, 'make_inputs'):
        return list(module.make_inputs())
    return [()]  # default: call with no arguments, once


def _normalize_input(item) -> tuple:
    """Turn one ``INPUTS`` entry into ``(args, kwargs)``.

    A ``list``/``tuple`` is positional args. A ``dict`` may carry ``'args'``/``'kwargs'``. Any
    other single value is treated as the sole positional argument.
    """
    if isinstance(item, dict):
        return tuple(item.get('args', ())), dict(item.get('kwargs', {}))
    if isinstance(item, (list, tuple)):
        return tuple(item), {}
    return (item,), {}


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Benchmark two implementations of a function and assert their outputs are "
                    "equal, reporting the speedup. See the 'safe-optimization' skill for the "
                    "full safe-performance-optimization workflow this supports.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    parser.add_argument('spec', type=Path,
                        help="path to a bench-spec .py file defining BASELINE, CANDIDATE, and "
                             "INPUTS (or make_inputs()) -- see this script's module docstring")
    parser.add_argument('--repeats', type=int, default=5,
                        help="timed calls per input, per implementation (default: 5)")
    parser.add_argument('--warmup', type=int, default=1,
                        help="untimed warm-up calls before timing (default: 1)")
    parser.add_argument('--rtol', type=float, default=0.0,
                        help="relative tolerance for numeric comparisons (default: 0.0, i.e. exact)")
    parser.add_argument('--atol', type=float, default=0.0,
                        help="absolute tolerance for numeric comparisons (default: 0.0, i.e. exact)")
    return parser


def main(argv=None) -> int:
    args = build_parser().parse_args(argv)

    module = _load_spec(args.spec)
    baseline = getattr(module, 'BASELINE', None)
    candidate = getattr(module, 'CANDIDATE', None)
    if baseline is None or candidate is None:
        print(f'[FAIL] spec {args.spec} must define both BASELINE and CANDIDATE callables', file=sys.stderr)
        return 1

    inputs = [_normalize_input(item) for item in _resolve_raw_inputs(module)]

    failures = 0
    for i, (call_args, call_kwargs) in enumerate(inputs):
        try:
            result = compare(baseline, candidate, call_args, call_kwargs, repeats=args.repeats,
                              warmup=args.warmup, rtol=args.rtol, atol=args.atol)
            print(f'[ok]   input {i}: outputs equal, {result.speedup:.2f}x faster '
                  f'(baseline {result.baseline_timing.best * 1000:.3f} ms -> '
                  f'candidate {result.candidate_timing.best * 1000:.3f} ms)')
        except AssertionError as exc:
            failures += 1
            print(f'[FAIL] input {i}: outputs differ -- NOT safe to ship: {exc}', file=sys.stderr)

    if failures:
        print(f'\n{failures} of {len(inputs)} input(s) produced a different output.', file=sys.stderr)
    return failures


if __name__ == '__main__':
    raise SystemExit(main())
