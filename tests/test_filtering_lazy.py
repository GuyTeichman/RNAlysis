"""Reproducibility harness for the lazy-evaluation work in ``filtering.Filter`` (epic #66, issue #139).

Converting eager operations to fused lazy query plans must not change results (hard invariant #5). CSV
reference files can't prove that -- serializing to text loses float precision, which is why the older
tests fall back to ``np.allclose``. This harness instead compares the production (lazy) path against an
independent *in-memory* eager oracle and asserts they are byte-for-byte identical, so any float or dtype
drift introduced by a lazy rewrite is caught immediately.
"""
import numpy as np
import polars as pl
import polars.selectors as cs

from rnalysis.filtering import CountFilter, Filter


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


def _eager_opposite(original_df: pl.DataFrame, kept_df: pl.DataFrame) -> pl.DataFrame:
    """Old ``_inplace`` opposite logic: the rows of self.df whose index is NOT among the kept rows,
    computed with the (now-deprecated) ``is_in`` anti-filter that preserves self.df's row order."""
    return original_df.filter(~pl.first().is_in(kept_df.select(pl.first()).to_series()))


def test_inplace_opposite_matches_eager_row_sum():
    cf = CountFilter('tests/test_files/counted.csv')
    kept = cf.filter_by_row_sum(5, opposite=False, inplace=False).df
    expected = _eager_opposite(cf.df, kept)
    result = cf.filter_by_row_sum(5, opposite=True, inplace=False).df
    assert_bit_identical(result, expected)


def test_inplace_opposite_matches_eager_percentile():
    f = Filter('tests/test_files/test_deseq.csv')
    kept = f.filter_percentile(0.75, 'log2FoldChange', opposite=False, inplace=False).df
    expected = _eager_opposite(f.df, kept)
    result = f.filter_percentile(0.75, 'log2FoldChange', opposite=True, inplace=False).df
    assert_bit_identical(result, expected)


def test_inplace_opposite_matches_eager_null_index():
    # edge case: a null in the index column. The old ~is_in anti-filter dropped null-index rows from the
    # opposite; the anti-join must match that (a naive anti-join would keep them) to stay bit-identical.
    df = pl.DataFrame({'': ['g1', None, 'g3', 'g4'], 'val': [1.0, 2.0, 3.0, 4.0]})
    f = Filter.from_dataframe(df, 'nulltest')
    kept = f.filter_percentile(0.75, 'val', opposite=False, inplace=False).df
    expected = _eager_opposite(f.df, kept)
    result = f.filter_percentile(0.75, 'val', opposite=True, inplace=False).df
    assert_bit_identical(result, expected)


def _per_gene_scaling_factors(cf: CountFilter) -> pl.DataFrame:
    """A per-gene scaling-factors table (index + one distinct value per gene per numeric column) -- the
    ``shape[0] > 1`` case that normalize_to_rpkm/tpm feed to _norm_scaling_factors."""
    numeric = cf._numeric_columns
    return cf.df.select(pl.first()).with_columns(
        [(pl.int_range(1, cf.df.height + 1) + j).cast(pl.Float64).alias(col) for j, col in enumerate(numeric)])


def _eager_norm_per_gene(cf: CountFilter, scaling_factors: pl.DataFrame) -> pl.DataFrame:
    """Old per-gene ``_norm_scaling_factors`` (one left-join per numeric column) -- the pre-fusion logic."""
    numeric = cf._numeric_columns
    out = pl.DataFrame().lazy()
    for column in cf.df.columns:
        if column in numeric:
            merged = cf.df.select(cs.first() | cs.by_name(column)).join(
                scaling_factors.select(cs.first() | cs.by_name(column)), left_on=cf.df.columns[0],
                right_on=scaling_factors.columns[0], how='left')
            merged_div = merged.with_columns((pl.nth(-2).truediv(pl.nth(-1))).alias('div'))
            out = out.with_columns((merged_div.select(pl.col('div').alias(column))))
        else:
            out = out.with_columns(cf.df[column].alias(column))
    return out.collect()


def test_norm_scaling_factors_per_gene_lazy_matches_eager():
    cf = CountFilter('tests/test_files/counted.csv')
    sf = _per_gene_scaling_factors(cf)
    expected = _eager_norm_per_gene(cf, sf)
    result = cf._norm_scaling_factors(sf)
    assert_bit_identical(result, expected)


def _eager_avg_mean(cf: CountFilter, grouping, names) -> pl.DataFrame:
    """Old ``_avg_subsamples`` mean path: one eager ``self.df.select().mean_horizontal()`` per group
    (single columns copied through) -- the pre-fusion logic."""
    out = cf.df.select(pl.first())
    for group, new_name in zip(grouping, names):
        if isinstance(group, str):
            out = out.with_columns(cf.df[group].alias(new_name))
        else:
            out = out.with_columns(cf.df.select(pl.col(group)).mean_horizontal().alias(new_name))
    return out


def test_avg_subsamples_mean_lazy_matches_eager():
    cf = CountFilter('tests/test_files/counted.csv')
    grouping = [['cond1', 'cond2'], 'cond3', ['cond4']]  # multi-col group, single str, single-item list
    names = ['groupA', 'cond3_solo', 'groupB']
    expected = _eager_avg_mean(cf, grouping, names)
    result = cf._avg_subsamples(grouping, 'mean', names)
    assert_bit_identical(result, expected)


def test_avg_subsamples_median_correct():
    # regression for a pre-existing crash: average_replicate_samples(function='median') used the
    # removed DataFrame.median_horizontal. It now computes a correct row-wise median; verify against numpy.
    cf = CountFilter('tests/test_files/counted.csv')
    grouping = [['cond1', 'cond2'], ['cond3', 'cond4']]
    names = ['g1', 'g2']
    result = cf._avg_subsamples(grouping, 'median', names)
    expected = cf.df.select(pl.first())
    for group, name in zip(grouping, names):
        med = [float(np.median(row)) for row in cf.df.select(group).to_numpy()]
        expected = expected.with_columns(pl.Series(name, med))
    assert_bit_identical(result, expected)
