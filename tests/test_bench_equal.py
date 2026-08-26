"""Tests for ``packaging/bench_equal.py``, the safe-optimization benchmark/equality harness.

Network-free and Qt-free by design (matches the ``safe-optimization`` skill's pure-logic
scope): timing wrapper, the numpy/polars-aware equality comparator, the speedup calculation,
and the CLI's bench-spec loading.

``packaging/`` is not an importable package (it has no ``__init__.py``, and its name collides
with the ``packaging`` PyPI distribution used elsewhere in this suite -- see
``tests/test_packaging.py``), so the module under test is loaded directly from its file path,
the same way ``test_packaging.py`` loads ``setup.py``.
"""

import importlib.util
import math
import sys
import time
from pathlib import Path

import numpy as np
import polars as pl
import pytest

REPO_ROOT = Path(__file__).resolve().parent.parent
BENCH_EQUAL_PATH = REPO_ROOT / 'packaging' / 'bench_equal.py'


def _load_bench_equal():
    spec = importlib.util.spec_from_file_location('rnalysis_bench_equal', BENCH_EQUAL_PATH)
    module = importlib.util.module_from_spec(spec)
    # Register in sys.modules *before* exec: bench_equal.py's dataclasses use postponed
    # (string) annotations, and dataclasses resolves those via sys.modules[cls.__module__] --
    # without this the resolution finds nothing and raises AttributeError on class creation.
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


@pytest.fixture(scope='module')
def bench_equal():
    return _load_bench_equal()


# ----------------------------------------------------------------------------------------
# time_call / Timing
# ----------------------------------------------------------------------------------------


class TestTimeCall:
    def test_returns_the_function_result(self, bench_equal):
        result, timing = bench_equal.time_call(lambda x, y: x + y, args=(2, 3), repeats=2, warmup=0)
        assert result == 5
        assert len(timing.times) == 2

    def test_passes_through_kwargs(self, bench_equal):
        result, _ = bench_equal.time_call(lambda x, y=0: x + y, args=(2,), kwargs={'y': 10}, repeats=1, warmup=0)
        assert result == 12

    def test_records_one_duration_per_repeat(self, bench_equal):
        calls = []

        def func():
            calls.append(1)
            return len(calls)

        _, timing = bench_equal.time_call(func, repeats=4, warmup=2)
        assert len(calls) == 6  # 2 warmup + 4 timed
        assert len(timing.times) == 4

    def test_durations_are_nonnegative(self, bench_equal):
        _, timing = bench_equal.time_call(lambda: sum(range(1000)), repeats=3, warmup=0)
        assert all(t >= 0 for t in timing.times)

    def test_timing_measures_a_slower_call_as_slower(self, bench_equal):
        # not a strict timing-precision assertion (flaky on loaded CI boxes) -- just checks
        # the wrapper actually times real elapsed wall-clock time, not e.g. always returning 0.
        _, fast = bench_equal.time_call(lambda: None, repeats=1, warmup=0)
        _, slow = bench_equal.time_call(lambda: time.sleep(0.02), repeats=1, warmup=0)
        assert slow.best > fast.best

    def test_rejects_zero_repeats(self, bench_equal):
        with pytest.raises(ValueError):
            bench_equal.time_call(lambda: None, repeats=0)

    def test_timing_stats(self, bench_equal):
        Timing = bench_equal.Timing
        timing = Timing((0.1, 0.2, 0.3))
        assert timing.best == pytest.approx(0.1)
        assert timing.median == pytest.approx(0.2)
        assert timing.mean == pytest.approx(0.2)


# ----------------------------------------------------------------------------------------
# speedup
# ----------------------------------------------------------------------------------------


class TestSpeedup:
    def test_candidate_twice_as_fast(self, bench_equal):
        baseline = bench_equal.Timing((0.20,))
        candidate = bench_equal.Timing((0.10,))
        assert bench_equal.speedup(baseline, candidate) == pytest.approx(2.0)

    def test_candidate_slower_gives_ratio_below_one(self, bench_equal):
        baseline = bench_equal.Timing((0.10,))
        candidate = bench_equal.Timing((0.20,))
        assert bench_equal.speedup(baseline, candidate) == pytest.approx(0.5)

    def test_zero_candidate_time_does_not_raise(self, bench_equal):
        baseline = bench_equal.Timing((0.10,))
        candidate = bench_equal.Timing((0.0,))
        assert math.isinf(bench_equal.speedup(baseline, candidate))

    def test_uses_best_of_n_not_mean(self, bench_equal):
        # one slow outlier call should not drag down a candidate that's otherwise fast --
        # best-of-N is the standard microbenchmark convention (what timeit.repeat + min does)
        # because it is robust to scheduler/GC noise inflating individual calls.
        baseline = bench_equal.Timing((0.10, 0.10, 0.10))
        candidate = bench_equal.Timing((0.01, 0.01, 0.50))
        assert bench_equal.speedup(baseline, candidate) == pytest.approx(10.0)


# ----------------------------------------------------------------------------------------
# assert_equal / equal -- plain Python values
# ----------------------------------------------------------------------------------------


class TestAssertEqualScalars:
    def test_equal_ints_pass(self, bench_equal):
        bench_equal.assert_equal(3, 3)

    def test_unequal_ints_fail(self, bench_equal):
        with pytest.raises(AssertionError):
            bench_equal.assert_equal(3, 4)

    def test_equal_floats_pass(self, bench_equal):
        bench_equal.assert_equal(1.5, 1.5)

    def test_unequal_floats_fail_when_exact(self, bench_equal):
        with pytest.raises(AssertionError):
            bench_equal.assert_equal(1.5, 1.5 + 1e-6)

    def test_close_floats_pass_with_tolerance(self, bench_equal):
        bench_equal.assert_equal(1.5, 1.5 + 1e-9, rtol=1e-6, atol=1e-6)

    def test_close_floats_fail_when_beyond_tolerance(self, bench_equal):
        with pytest.raises(AssertionError):
            bench_equal.assert_equal(1.5, 1.6, rtol=1e-6, atol=1e-6)

    def test_nan_equals_nan(self, bench_equal):
        bench_equal.assert_equal(float('nan'), float('nan'))

    def test_strings_equal_and_unequal(self, bench_equal):
        bench_equal.assert_equal('abc', 'abc')
        with pytest.raises(AssertionError):
            bench_equal.assert_equal('abc', 'abd')

    def test_nested_list_and_dict_structures(self, bench_equal):
        a = {'x': [1, 2, {'y': 3.0}], 'z': (1, 2)}
        b = {'x': [1, 2, {'y': 3.0}], 'z': (1, 2)}
        bench_equal.assert_equal(a, b)

    def test_nested_structure_mismatch_is_detected(self, bench_equal):
        a = {'x': [1, 2, {'y': 3.0}]}
        b = {'x': [1, 2, {'y': 3.1}]}
        with pytest.raises(AssertionError):
            bench_equal.assert_equal(a, b)

    def test_dict_key_mismatch_fails(self, bench_equal):
        with pytest.raises(AssertionError):
            bench_equal.assert_equal({'a': 1}, {'b': 1})

    def test_sequence_length_mismatch_fails(self, bench_equal):
        with pytest.raises(AssertionError):
            bench_equal.assert_equal([1, 2, 3], [1, 2])


class TestEqualWrapper:
    def test_returns_true_without_raising(self, bench_equal):
        assert bench_equal.equal(3, 3) is True

    def test_returns_false_without_raising(self, bench_equal):
        assert bench_equal.equal(3, 4) is False


# ----------------------------------------------------------------------------------------
# assert_equal -- numpy arrays
# ----------------------------------------------------------------------------------------


class TestAssertEqualNumpy:
    def test_identical_arrays_pass(self, bench_equal):
        a = np.array([1.0, 2.0, 3.0])
        bench_equal.assert_equal(a, a.copy())

    def test_different_values_fail_exact(self, bench_equal):
        a = np.array([1.0, 2.0, 3.0])
        b = np.array([1.0, 2.0, 3.000001])
        with pytest.raises(AssertionError):
            bench_equal.assert_equal(a, b)

    def test_tiny_difference_passes_with_tolerance(self, bench_equal):
        a = np.array([1.0, 2.0, 3.0])
        b = a + 1e-12
        bench_equal.assert_equal(a, b, rtol=1e-8, atol=1e-8)

    def test_tiny_difference_still_fails_when_exact(self, bench_equal):
        a = np.array([1.0, 2.0, 3.0])
        b = a + 1e-12
        with pytest.raises(AssertionError):
            bench_equal.assert_equal(a, b)

    def test_shape_mismatch_fails(self, bench_equal):
        with pytest.raises(AssertionError):
            bench_equal.assert_equal(np.zeros((2, 3)), np.zeros((3, 2)))

    def test_dtype_mismatch_fails(self, bench_equal):
        with pytest.raises(AssertionError):
            bench_equal.assert_equal(np.array([1, 2, 3], dtype=np.int32), np.array([1, 2, 3], dtype=np.int64))

    def test_nan_in_same_position_is_equal(self, bench_equal):
        a = np.array([1.0, np.nan, 3.0])
        b = np.array([1.0, np.nan, 3.0])
        bench_equal.assert_equal(a, b)

    def test_nan_vs_number_is_not_equal(self, bench_equal):
        a = np.array([1.0, np.nan, 3.0])
        b = np.array([1.0, 2.0, 3.0])
        with pytest.raises(AssertionError):
            bench_equal.assert_equal(a, b)

    def test_integer_arrays_exact(self, bench_equal):
        bench_equal.assert_equal(np.array([1, 2, 3]), np.array([1, 2, 3]))
        with pytest.raises(AssertionError):
            bench_equal.assert_equal(np.array([1, 2, 3]), np.array([1, 2, 4]))


# ----------------------------------------------------------------------------------------
# assert_equal -- polars DataFrames / Series
# ----------------------------------------------------------------------------------------


class TestAssertEqualPolars:
    def test_identical_dataframes_pass(self, bench_equal):
        df = pl.DataFrame({'a': [1, 2, 3], 'b': [1.0, 2.0, 3.0]})
        bench_equal.assert_equal(df, df.clone())

    def test_different_values_fail_exact(self, bench_equal):
        df1 = pl.DataFrame({'a': [1, 2, 3]})
        df2 = pl.DataFrame({'a': [1, 2, 4]})
        with pytest.raises(AssertionError):
            bench_equal.assert_equal(df1, df2)

    def test_float_tolerance_passes_close_values(self, bench_equal):
        df1 = pl.DataFrame({'a': [1.0, 2.0, 3.0]})
        df2 = pl.DataFrame({'a': [1.0 + 1e-12, 2.0, 3.0]})
        bench_equal.assert_equal(df1, df2, rtol=1e-8, atol=1e-8)

    def test_column_name_mismatch_fails(self, bench_equal):
        df1 = pl.DataFrame({'a': [1, 2]})
        df2 = pl.DataFrame({'b': [1, 2]})
        with pytest.raises(AssertionError):
            bench_equal.assert_equal(df1, df2)

    def test_row_count_mismatch_fails(self, bench_equal):
        df1 = pl.DataFrame({'a': [1, 2, 3]})
        df2 = pl.DataFrame({'a': [1, 2]})
        with pytest.raises(AssertionError):
            bench_equal.assert_equal(df1, df2)

    def test_dtype_mismatch_fails(self, bench_equal):
        df1 = pl.DataFrame({'a': [1, 2, 3]}, schema={'a': pl.Int32})
        df2 = pl.DataFrame({'a': [1, 2, 3]}, schema={'a': pl.Int64})
        with pytest.raises(AssertionError):
            bench_equal.assert_equal(df1, df2)

    def test_null_in_same_position_is_equal(self, bench_equal):
        df1 = pl.DataFrame({'a': [1, None, 3]})
        df2 = pl.DataFrame({'a': [1, None, 3]})
        bench_equal.assert_equal(df1, df2)

    def test_series_equal_and_unequal(self, bench_equal):
        s1 = pl.Series('a', [1.0, 2.0, 3.0])
        s2 = pl.Series('a', [1.0, 2.0, 3.0])
        bench_equal.assert_equal(s1, s2)
        s3 = pl.Series('a', [1.0, 2.0, 3.1])
        with pytest.raises(AssertionError):
            bench_equal.assert_equal(s1, s3)


# ----------------------------------------------------------------------------------------
# compare() -- the timing + equality combinator
# ----------------------------------------------------------------------------------------


class TestCompare:
    def test_equal_implementations_report_speedup(self, bench_equal):
        def slow(n):
            total = 0
            for i in range(n):
                total += i
            return total

        def fast(n):
            return n * (n - 1) // 2

        result = bench_equal.compare(slow, fast, args=(20000,), repeats=2, warmup=1)
        assert result.result == slow(20000)
        assert result.speedup > 0

    def test_diverging_implementations_raise(self, bench_equal):
        def baseline(n):
            return n + 1

        def buggy_candidate(n):
            return n + 2

        with pytest.raises(AssertionError):
            bench_equal.compare(baseline, buggy_candidate, args=(5,), repeats=1, warmup=0)

    def test_tolerance_is_forwarded(self, bench_equal):
        def baseline():
            return 1.0

        def candidate():
            return 1.0 + 1e-10

        # exact fails
        with pytest.raises(AssertionError):
            bench_equal.compare(baseline, candidate, repeats=1, warmup=0)
        # tolerant passes
        bench_equal.compare(baseline, candidate, repeats=1, warmup=0, rtol=1e-6, atol=1e-6)


# ----------------------------------------------------------------------------------------
# CLI: bench-spec loading and main()
# ----------------------------------------------------------------------------------------


class TestCli:
    def _write_spec(self, tmp_path, body: str) -> Path:
        spec_path = tmp_path / 'spec.py'
        spec_path.write_text(body, encoding='utf-8')
        return spec_path

    def test_main_returns_zero_on_matching_implementations(self, bench_equal, tmp_path, capsys):
        spec = self._write_spec(
            tmp_path,
            """
def _baseline(n):
    return sum(range(n))

def _candidate(n):
    return n * (n - 1) // 2

BASELINE = _baseline
CANDIDATE = _candidate
INPUTS = [(100,), (1000,)]
""",
        )
        exit_code = bench_equal.main([str(spec), '--repeats', '1', '--warmup', '0'])
        assert exit_code == 0
        out = capsys.readouterr().out
        assert '[ok]' in out

    def test_main_returns_nonzero_on_mismatch(self, bench_equal, tmp_path):
        spec = self._write_spec(
            tmp_path,
            """
BASELINE = lambda n: n
CANDIDATE = lambda n: n + 1
INPUTS = [(1,), (2,)]
""",
        )
        exit_code = bench_equal.main([str(spec), '--repeats', '1', '--warmup', '0'])
        assert exit_code == 2  # both inputs mismatch

    def test_main_supports_make_inputs_callable(self, bench_equal, tmp_path):
        spec = self._write_spec(
            tmp_path,
            """
BASELINE = lambda n: n * 2
CANDIDATE = lambda n: n + n

def make_inputs():
    return [(3,), (4,)]
""",
        )
        exit_code = bench_equal.main([str(spec), '--repeats', '1', '--warmup', '0'])
        assert exit_code == 0

    def test_main_supports_dict_style_inputs_with_kwargs(self, bench_equal, tmp_path):
        spec = self._write_spec(
            tmp_path,
            """
def _baseline(n, offset=0):
    return n + offset

BASELINE = _baseline
CANDIDATE = _baseline
INPUTS = [{'args': (5,), 'kwargs': {'offset': 2}}]
""",
        )
        exit_code = bench_equal.main([str(spec), '--repeats', '1', '--warmup', '0'])
        assert exit_code == 0

    def test_main_defaults_to_single_no_arg_call_without_inputs(self, bench_equal, tmp_path):
        spec = self._write_spec(
            tmp_path,
            """
BASELINE = lambda: 42
CANDIDATE = lambda: 42
""",
        )
        exit_code = bench_equal.main([str(spec), '--repeats', '1', '--warmup', '0'])
        assert exit_code == 0

    def test_main_reports_missing_baseline_or_candidate(self, bench_equal, tmp_path):
        spec = self._write_spec(tmp_path, 'CANDIDATE = lambda: 1\n')
        exit_code = bench_equal.main([str(spec), '--repeats', '1'])
        assert exit_code != 0

    def test_help_does_not_raise(self, bench_equal, capsys):
        with pytest.raises(SystemExit) as exc_info:
            bench_equal.main(['--help'])
        assert exc_info.value.code == 0
        out = capsys.readouterr().out
        assert 'usage' in out.lower()

    def test_rtol_atol_cli_flags_are_forwarded(self, bench_equal, tmp_path):
        spec = self._write_spec(
            tmp_path,
            """
BASELINE = lambda: 1.0
CANDIDATE = lambda: 1.0 + 1e-10
""",
        )
        # fails without tolerance flags
        assert bench_equal.main([str(spec), '--repeats', '1', '--warmup', '0']) != 0
        # passes with tolerance flags
        assert bench_equal.main([str(spec), '--repeats', '1', '--warmup', '0', '--rtol', '1e-6', '--atol', '1e-6']) == 0
