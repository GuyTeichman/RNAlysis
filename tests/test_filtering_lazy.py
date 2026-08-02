"""Reproducibility harness for the lazy-evaluation work in ``filtering.Filter`` (epic #66, issue #139).

Converting eager operations to fused lazy query plans must not change results (hard invariant #5). CSV
reference files can't prove that -- serializing to text loses float precision, which is why the older
tests fall back to ``np.allclose``. This harness instead compares the production (lazy) path against an
independent *in-memory* eager oracle and asserts they are byte-for-byte identical, so any float or dtype
drift introduced by a lazy rewrite is caught immediately.
"""
import polars as pl

from rnalysis.filtering import CountFilter


def assert_bit_identical(result: pl.DataFrame, reference: pl.DataFrame):
    """Assert two frames are the same table down to schema and null placement.

    Stricter than a value-only / ``np.allclose`` check: identical column names, order and dtypes
    (schema), identical shape, and identical values including nulls. This is the equivalence contract a
    lazy rewrite must satisfy to preserve reproducibility.
    """
    assert result.schema == reference.schema, f"schema differs:\n{result.schema}\n!=\n{reference.schema}"
    assert result.shape == reference.shape, f"shape differs: {result.shape} != {reference.shape}"
    assert result.equals(reference), "values differ between the lazy result and the eager oracle"


def _eager_fold_change(count_filter: CountFilter, numerator, denominator) -> pl.DataFrame:
    """Independent eager oracle for :meth:`CountFilter.fold_change` -- the pre-lazy computation, kept
    here so the production (lazy) path can be proven bit-identical to a straightforward eager one."""
    df = count_filter.df
    return df.select(pl.first()).with_columns(
        ((df.select(numerator).mean_horizontal() + 1) / (
            df.select(denominator).mean_horizontal() + 1)).alias('Fold Change'))


def test_fold_change_lazy_matches_eager():
    h = CountFilter('tests/test_files/counted_fold_change.csv')
    numerator = ['cond1_rep1', 'cond1_rep2']
    denominator = ['cond2_rep1', 'cond2_rep2']
    expected = _eager_fold_change(h, numerator, denominator)
    result = h.fold_change(numerator, denominator).df
    assert_bit_identical(result, expected)


def _eager_filter_low_reads(cf: CountFilter, threshold, n_samples) -> pl.DataFrame:
    """Eager oracle for CountFilter.filter_low_reads: two scans of self.df (build the count mask on a
    selected copy, then filter the full frame) -- what the method did before being fused into one pass."""
    df = cf.df
    mask = df.select(pl.col(cf._numeric_columns)).with_columns(
        pl.col(cf._numeric_columns) >= threshold).sum_horizontal() >= n_samples
    return df.filter(mask)


def _eager_split_by_reads(cf: CountFilter, threshold):
    df = cf.df
    high = df.filter(df.select(pl.col(cf._numeric_columns)).max_horizontal() >= threshold)
    low = df.filter(df.select(pl.col(cf._numeric_columns)).max_horizontal() < threshold)
    return high, low


def _eager_filter_by_row_sum(cf: CountFilter, threshold) -> pl.DataFrame:
    df = cf.df
    return df.filter(df.select(pl.col(cf._numeric_columns)).sum_horizontal() >= threshold)


def test_filter_low_reads_lazy_matches_eager():
    cf = CountFilter('tests/test_files/counted.csv')
    expected = _eager_filter_low_reads(cf, 5, 2)
    result = cf.filter_low_reads(5, n_samples=2, inplace=False).df
    assert_bit_identical(result, expected)


def test_split_by_reads_lazy_matches_eager():
    cf = CountFilter('tests/test_files/counted.csv')
    exp_high, exp_low = _eager_split_by_reads(cf, 5)
    high, low = cf.split_by_reads(5)
    assert_bit_identical(high.df, exp_high)
    assert_bit_identical(low.df, exp_low)


def test_filter_by_row_sum_lazy_matches_eager():
    cf = CountFilter('tests/test_files/counted.csv')
    expected = _eager_filter_by_row_sum(cf, 5)
    result = cf.filter_by_row_sum(5, inplace=False).df
    assert_bit_identical(result, expected)


def _eager_norm_apply(cf: CountFilter, scaling_factors: pl.DataFrame) -> pl.DataFrame:
    """Old N-pass ``_norm_scaling_factors`` (single-row scaling-factors case): divides each numeric
    column by its scalar factor via one eager ``self.df.select`` per column -- the pre-fusion logic."""
    numeric = cf._numeric_columns
    out = pl.DataFrame().lazy()
    for column in cf.df.columns:
        if column in numeric:
            out = out.with_columns(cf.df.select(pl.col(column).truediv(scaling_factors[column])))
        else:
            out = out.with_columns(cf.df[column].alias(column))
    return out.collect()


def _eager_normalize_to_rpm(cf: CountFilter) -> pl.DataFrame:
    sf = cf.df.select(pl.col(cf._numeric_columns)).sum() / (10 ** 6)
    return _eager_norm_apply(cf, sf)


def _eager_normalize_to_quantile(cf: CountFilter, quantile: float = 0.75) -> pl.DataFrame:
    data = cf.df.select(pl.col(cf._numeric_columns))
    expressed = data.filter(data.sum_horizontal() != 0)
    quantiles = expressed.quantile(quantile, interpolation='linear')
    sf = quantiles / quantiles.mean_horizontal()
    return _eager_norm_apply(cf, sf)


def test_normalize_to_rpm_lazy_matches_eager():
    cf = CountFilter('tests/test_files/counted.csv')
    expected = _eager_normalize_to_rpm(cf)
    result = cf.normalize_to_rpm(inplace=False).df
    assert_bit_identical(result, expected)


def test_normalize_to_quantile_lazy_matches_eager():
    cf = CountFilter('tests/test_files/counted.csv')
    expected = _eager_normalize_to_quantile(cf, 0.75)
    result = cf.normalize_to_quantile(0.75, inplace=False).df
    assert_bit_identical(result, expected)
