import json
import shutil
from unittest.mock import Mock, patch

import matplotlib
import pytest
import yaml

from rnalysis import __version__
from rnalysis.exceptions import InvalidTypeError, InvalidValueError, RNAlysisInputError
from rnalysis.filtering import *
from tests import __attr_ref__, __biotype_ref__

matplotlib.use('Agg')


@pytest.fixture(scope='session')
def _basic_filter():
    return Filter('tests/test_files/counted.csv')


@pytest.fixture
def basic_filter(_basic_filter):
    return _basic_filter.__copy__()


@pytest.fixture(scope='session')
def _basic_countfilter():
    return CountFilter('tests/test_files/counted.csv')


@pytest.fixture
def basic_countfilter(_basic_countfilter):
    return _basic_countfilter.__copy__()


@pytest.fixture(scope='session')
def _basic_deseqfilter():
    return DESeqFilter('tests/test_files/test_deseq.csv')


@pytest.fixture
def basic_deseqfilter(_basic_deseqfilter):
    return _basic_deseqfilter.__copy__()


@pytest.fixture(scope='session')
def _big_countfilter():
    return CountFilter('tests/test_files/big_counted.csv')


@pytest.fixture
def big_countfilter(_big_countfilter):
    return _big_countfilter.__copy__()


@pytest.fixture(scope='session')
def _clustering_countfilter(_big_countfilter):
    return _big_countfilter.normalize_to_rpm(inplace=False).filter_low_reads(105, 3, inplace=False)


@pytest.fixture
def clustering_countfilter(_clustering_countfilter):
    return _clustering_countfilter.__copy__()


def unlink_tree(dir):
    for item in Path(dir).iterdir():
        if 'gitignore' in item.name:
            continue
        if item.is_file():
            item.unlink()
        else:
            shutil.rmtree(item)


def test_filter_api():
    f = Filter('tests/test_files/uncounted.csv')
    assert f.__str__() == "Generic table named 'uncounted.csv'"
    assert f.__repr__() == "Filter('tests/test_files/uncounted.csv')"


def test_countfilter_api(basic_countfilter):
    assert basic_countfilter.__str__() == "Count matrix named 'counted.csv'"
    assert basic_countfilter.__repr__() == "CountFilter('tests/test_files/counted.csv')"


def test_deseqfilter_api(basic_deseqfilter):
    assert basic_deseqfilter.__str__() == "Differential expression table named 'test_deseq.csv'"
    assert basic_deseqfilter.__repr__() == "DESeqFilter('tests/test_files/test_deseq.csv')"


def test_foldchangefilter_api():
    fc = FoldChangeFilter("tests/test_files/fc_1.csv", 'a', 'b')
    assert fc.__str__() == "FoldChangeFilter (numerator: 'a', denominator: 'b') of file fc_1.csv"
    assert fc.__repr__() == "FoldChangeFilter('tests/test_files/fc_1.csv', 'a', 'b')"


def test_filter_contains(basic_countfilter, basic_deseqfilter, basic_filter):
    objs = [basic_deseqfilter, basic_countfilter, basic_filter,
            FoldChangeFilter('tests/test_files/fc_1.csv', 'num', 'denom')]
    neither = ['WBGene' + str(i) for i in range(10)]
    for obj in objs:
        for ind in obj.df.select(pl.first()):
            assert ind in obj
        for false_ind in neither:
            assert false_ind not in obj


@pytest.mark.parametrize('filter_obj',
                         [Filter('tests/test_files/test_deseq.csv'), CountFilter('tests/test_files/counted.csv'),
                          DESeqFilter('tests/test_files/counted.csv'),
                          FoldChangeFilter('tests/test_files/fc_1.csv', 'num', 'denom')])
def test_filter_len(filter_obj):
    assert len(filter_obj) == filter_obj.df.shape[0]


def test_filter_inplace():
    d = DESeqFilter('tests/test_files/test_deseq_no_nans.csv')
    d_copy = DESeqFilter('tests/test_files/test_deseq_no_nans.csv')
    truth = io.load_table('tests/test_files/counted.csv')
    d_inplace_false = d._inplace(truth, opposite=False, inplace=False, suffix='suffix')
    assert d_inplace_false.df.equals(truth)
    assert d.df.equals(d_copy.df)
    d._inplace(truth, opposite=False, inplace=True, suffix='other_suffix')
    assert d.df.equals(truth)


def test_head(basic_filter, basic_deseqfilter):
    df = io.load_table('tests/test_files/test_deseq.csv')
    assert np.all(df.head(7) == basic_deseqfilter.head(7))
    assert np.all(df.head(1) == basic_deseqfilter.head(1))

    df2 = io.load_table('tests/test_files/counted.csv')
    assert np.all(df2.head() == basic_filter.head())
    assert np.all(df2.head(1000) == basic_filter.head(1000))


def test_tail(basic_filter, basic_deseqfilter):
    df = io.load_table('tests/test_files/test_deseq.csv')
    assert np.all(df.tail(7) == basic_deseqfilter.tail(7))
    assert np.all(df.tail(1) == basic_deseqfilter.tail(1))

    df2 = io.load_table('tests/test_files/counted.csv')
    assert np.all(df2.tail() == basic_filter.tail())
    assert np.all(df2.tail(1000) == basic_filter.tail(1000))


def test_describe(basic_countfilter, basic_deseqfilter):
    fc_df = io.load_table('tests/test_files/fc_1.csv', 0, squeeze=True)
    count_df = basic_countfilter.df.clone()
    deseq_df = basic_deseqfilter.df.clone()

    fc = FoldChangeFilter('tests/test_files/fc_1.csv', 'a', 'b')
    default_percentiles = (0.01, 0.25, 0.5, 0.75, 0.99)
    deciles = [i / 10 for i in range(1, 10)]
    for filter_obj, df in zip([basic_countfilter, basic_deseqfilter, fc], [count_df, deseq_df, fc_df]):
        assert filter_obj.describe().equals(df.describe(percentiles=default_percentiles))
        assert filter_obj.describe(deciles).equals(df.describe(percentiles=deciles))


def test_print_features_api(basic_countfilter):
    basic_countfilter.print_features()


def test_from_string(monkeypatch):
    monkeypatch.setattr('builtins.input', lambda x: 'first-word\nsecond word\nthird. word\n')
    assert Filter._from_string("msg") == ["first-word", "second word", "third. word"]
    monkeypatch.setattr('builtins.input', lambda x: 'first-word,second word,third. word')
    assert Filter._from_string("msg", delimiter=',') == ["first-word", "second word", "third. word"]


@pytest.mark.parametrize("map_to,map_from,remove_unmapped_genes,,expected", [
    ('UniProtKB AC/ID', 'WormBase', False, 'tests/test_files/counted_translated_with_unmapped.csv'),
    ('UniProtKB', 'auto', True, 'tests/test_files/counted_translated_remove_unmapped.csv'),
])
def test_filter_translate_gene_ids(map_to, map_from, remove_unmapped_genes, expected, monkeypatch, basic_countfilter):
    def mock_map_gene_ids(self, ids):
        if self.map_from == 'WormBase':
            return io.GeneIDDict(
                {'WBGene00007063': 'A0A0K3AWR5', 'WBGene00007064': 'A0A2X0T1Z3', 'WBGene00007067': 'D3NQA2',
                 'WBGene00077503': 'H2L2B5', 'WBGene00007071': 'Q17405', 'WBGene00014997': 'Q7JNR0',
                 'WBGene00043988': 'A4F2Z7', 'WBGene00043989': 'G5EFZ2', 'WBGene00007075': 'G5EDW3',
                 'WBGene00007076': 'G5EFZ2'})
        return io.GeneIDDict({})

    monkeypatch.setattr(io.GeneIDTranslator, 'run', mock_map_gene_ids)
    truth = io.load_table(expected)

    res = basic_countfilter.translate_gene_ids(map_to, map_from, remove_unmapped_genes, inplace=False)
    assert res.df.sort(pl.first()).equals(truth.sort(pl.first()))
    basic_countfilter.translate_gene_ids(map_to, map_from, remove_unmapped_genes, inplace=True)
    assert basic_countfilter.df.sort(pl.first()).equals(truth.sort(pl.first()))


def test_countfilter_normalize_to_rpm_htseqcount(basic_countfilter):
    truth = io.load_table(r"tests/test_files/test_norm_reads_rpm_htseqcount.csv")
    not_inplace, factors = basic_countfilter.normalize_to_rpm_htseqcount("tests/test_files/uncounted.csv",
                                                                         inplace=False,
                                                                         return_scaling_factors=True)
    assert np.isclose(truth.drop(cs.first()), not_inplace.df.drop(cs.first())).all()

    norm_df = basic_countfilter._norm_scaling_factors(factors)
    assert np.isclose(truth.drop(cs.first()), norm_df.drop(cs.first())).all()

    basic_countfilter.normalize_to_rpm_htseqcount("tests/test_files/uncounted.csv")
    assert np.isclose(truth.drop(cs.first()), basic_countfilter.df.drop(cs.first())).all()


def test_countfilter_normalize_to_rpm(basic_countfilter):
    truth = io.load_table(r"tests/test_files/test_norm_to_rpm.csv")
    not_inplace, factors = basic_countfilter.normalize_to_rpm(inplace=False, return_scaling_factors=True)
    assert np.isclose(truth.drop(cs.first()), not_inplace.df.drop(cs.first())).all()

    norm_df = basic_countfilter._norm_scaling_factors(factors)
    assert np.isclose(truth.drop(cs.first()), norm_df.drop(cs.first())).all()

    basic_countfilter.normalize_to_rpm()
    assert np.isclose(truth.drop(cs.first()), basic_countfilter.df.drop(cs.first())).all()


@pytest.mark.parametrize('gtf_path,feature_type,method',
                         [
                             ('tests/test_files/test_gtf_wormbase.gtf', 'gene', 'merged_exons'),
                             ('tests/test_files/test_gff3_wormbase.gff3', 'gene', 'geometric_mean'),
                             ('tests/test_files/test_gtf_wormbase.gtf', 'transcript', 'mean'),
                             ('tests/test_files/test_gff3_wormbase.gff3', 'transcript', 'max')
                         ])
def test_countfilter_normalize_to_rpkm(monkeypatch, gtf_path, feature_type, method, basic_countfilter):
    def mock_get_feature_lengths(this_gtf_path, this_feature_type, this_method):
        assert this_gtf_path == gtf_path
        assert this_feature_type == feature_type
        assert this_method == method
        with open('tests/test_files/feature_lengths.json') as f:
            return json.load(f)

    monkeypatch.setattr(genome_annotation, 'get_genomic_feature_lengths', mock_get_feature_lengths)
    truth = io.load_table(r"tests/test_files/test_norm_to_rpkm.csv")
    not_inplace, factors = basic_countfilter.normalize_to_rpkm(gtf_path, feature_type, method, inplace=False,
                                                               return_scaling_factors=True)
    assert np.isclose(truth.drop(cs.first()), not_inplace.df.drop(cs.first())).all()

    norm_df = basic_countfilter._norm_scaling_factors(factors)
    assert np.isclose(truth.drop(cs.first()), norm_df.drop(cs.first())).all()

    basic_countfilter.normalize_to_rpkm(gtf_path, feature_type, method, )
    assert np.isclose(truth.drop(cs.first()), basic_countfilter.df.drop(cs.first())).all()


@pytest.mark.parametrize('gtf_path,feature_type,method',
                         [
                             ('tests/test_files/test_gtf_wormbase.gtf', 'gene', 'merged_exons'),
                             ('tests/test_files/test_gff3_wormbase.gff3', 'gene', 'geometric_mean'),
                             ('tests/test_files/test_gtf_wormbase.gtf', 'transcript', 'mean'),
                             ('tests/test_files/test_gff3_wormbase.gff3', 'transcript', 'max')
                         ])
def test_countfilter_normalize_to_tpm(monkeypatch, gtf_path, feature_type, method, basic_countfilter):
    def mock_get_feature_lengths(this_gtf_path, this_feature_type, this_method):
        assert this_gtf_path == gtf_path
        assert this_feature_type == feature_type
        assert this_method == method
        with open('tests/test_files/feature_lengths.json') as f:
            return json.load(f)

    monkeypatch.setattr(genome_annotation, 'get_genomic_feature_lengths', mock_get_feature_lengths)
    truth = io.load_table(r"tests/test_files/test_norm_to_tpm.csv")
    not_inplace, factors = basic_countfilter.normalize_to_tpm(gtf_path, feature_type, method, inplace=False,
                                                              return_scaling_factors=True)
    assert np.isclose(truth.drop(cs.first()), not_inplace.df.drop(cs.first())).all()

    norm_df = basic_countfilter._norm_scaling_factors(factors)
    assert np.isclose(truth.drop(cs.first()), norm_df.drop(cs.first())).all()

    basic_countfilter.normalize_to_tpm(gtf_path, feature_type, method)
    assert np.isclose(truth.drop(cs.first()), basic_countfilter.df.drop(cs.first())).all()


def test_countfilter_normalize_rle(basic_countfilter):
    truth = io.load_table(r"tests/test_files/test_norm_rle.csv")
    not_inplace, factors = basic_countfilter.normalize_rle(inplace=False, return_scaling_factors=True)
    assert np.allclose(truth.drop(cs.first()), not_inplace.df.drop(cs.first()))

    norm_df = basic_countfilter._norm_scaling_factors(factors)
    assert np.allclose(truth.drop(cs.first()), norm_df.drop(cs.first()))
    basic_countfilter.normalize_rle()
    assert np.allclose(truth.drop(cs.first()), basic_countfilter.df.drop(cs.first()))


def test_countfilter_normalize_tmm(basic_countfilter):
    truth = io.load_table(r"tests/test_files/test_norm_tmm_truth.csv")
    not_inplace, factors = basic_countfilter.normalize_tmm(ref_column='cond1', inplace=False,
                                                           return_scaling_factors=True)

    print(truth)
    print(not_inplace.df)
    assert np.allclose(truth.drop(cs.first()), not_inplace.df.drop(cs.first()))

    norm_df = basic_countfilter._norm_scaling_factors(factors)
    assert np.allclose(truth.drop(cs.first()), norm_df.drop(cs.first()))

    basic_countfilter.normalize_tmm(ref_column='cond1')
    assert np.allclose(truth.drop(cs.first()), basic_countfilter.df.drop(cs.first()))


def test_countfilter_normalize_median_of_ratios(basic_countfilter):
    truth = io.load_table(r"tests/test_files/test_norm_mrn.csv")
    not_inplace, factors = basic_countfilter.normalize_median_of_ratios([['cond1', 'cond2'], ['cond3', 'cond4']],
                                                                        inplace=False,
                                                                        return_scaling_factors=True)
    assert np.isclose(truth.drop(cs.first()), not_inplace.df.drop(cs.first())).all()

    norm_df = basic_countfilter._norm_scaling_factors(factors)
    assert np.isclose(truth.drop(cs.first()), norm_df.drop(cs.first())).all()

    basic_countfilter.normalize_median_of_ratios([['cond1', 'cond2'], ['cond3', 'cond4']])
    assert np.isclose(truth.drop(cs.first()), basic_countfilter.df.drop(cs.first())).all()


@pytest.mark.parametrize('quantile,truth_path', [
    (0.75, "tests/test_files/test_norm_quantile_75.csv"),
    (0.32, "tests/test_files/test_norm_quantile_32.csv"),
])
def test_countfilter_normalize_to_quantile(quantile, truth_path, basic_countfilter):
    truth = io.load_table(truth_path)
    not_inplace, factors = basic_countfilter.normalize_to_quantile(quantile, inplace=False, return_scaling_factors=True)
    assert np.isclose(truth.drop(cs.first()), not_inplace.df.drop(cs.first())).all()

    norm_df = basic_countfilter._norm_scaling_factors(factors)
    assert np.isclose(truth.drop(cs.first()), norm_df.drop(cs.first())).all()

    basic_countfilter.normalize_to_quantile(quantile)
    assert np.isclose(truth.drop(cs.first()), basic_countfilter.df.drop(cs.first())).all()


class TestFilterConcatenate:
    def setup_method(self):
        # Create Filter objects for testing
        df1 = pl.DataFrame({'index': ['a', 'b', 'c'], 'A': [1, 2, 3], 'B': [4, 5, 6]})
        df2 = pl.DataFrame({'index': ['d', 'e', 'f'], 'A': [7, 8, 9], 'B': [10, 11, 12]})
        self.filter1 = Filter.from_dataframe(df1, 'filter1.csv')
        self.filter2 = Filter.from_dataframe(df2, 'filter2.csv')

    def test_concatenate_with_valid_filter(self):
        # Call the concatenate method
        result = self.filter1.concatenate(self.filter2)

        # Assert that the resulting Filter object has the expected data and filename
        expected_df = pl.DataFrame(
            {'index': ['a', 'b', 'c', 'd', 'e', 'f'], 'A': [1, 2, 3, 7, 8, 9], 'B': [4, 5, 6, 10, 11, 12]})
        assert result.df.equals(expected_df)
        assert result.fname == Path('filter1_filter2.csv')

    def test_concatenate_with_invalid_filter(self):
        # Create a non-Filter object
        non_filter_obj = "not a Filter object"

        # Assert that concatenating with a non-Filter object raises an AssertionError
        with pytest.raises(InvalidTypeError):
            self.filter1.concatenate(non_filter_obj)

    def test_concatenate_with_different_columns(self):
        # Create a Filter object with different columns
        df3 = pl.DataFrame({'index': ['d', 'e', 'f'], 'C': [7, 8, 9], 'D': [10, 11, 12]})
        filter3 = Filter.from_dataframe(df3, 'filter3.csv')

        # Assert that concatenating Filter objects with different columns raises an AssertionError
        with pytest.raises(InvalidValueError):
            self.filter1.concatenate(filter3)

    def test_concatenate_with_overlapping_indices(self):
        # Create a Filter object with overlapping indices
        df4 = pl.DataFrame({'index': ['c', 'd', 'e'], 'A': [7, 8, 9], 'B': [10, 11, 12]})
        filter4 = Filter.from_dataframe(df4, 'filter4.csv')

        # Assert that concatenating Filter objects with overlapping indices raises an AssertionError
        with pytest.raises(InvalidValueError):
            self.filter1.concatenate(filter4)


def test_countfilter_norm_reads_with_scaling_factors(basic_countfilter):
    truth = io.load_table(r"tests/test_files/test_norm_scaling_factors.csv")
    factors = io.load_table("tests/test_files/scaling_factors.csv")
    h_norm = basic_countfilter.normalize_with_scaling_factors("tests/test_files/scaling_factors.csv", inplace=False)
    assert np.isclose(truth.drop(cs.first()), h_norm.df.drop(cs.first())).all()
    basic_countfilter.normalize_with_scaling_factors(factors)
    assert np.isclose(truth.drop(cs.first()), basic_countfilter.df.drop(cs.first())).all()


@pytest.mark.parametrize('input_path,threshold,n_samples,opposite,truth_path', [
    ("tests/test_files/counted_low_rpm.csv", 5, 1, False, "tests/test_files/counted_low_rpm_truth.csv"),
    ("tests/test_files/counted_low_rpm.csv", 6, 2, False, "tests/test_files/counted_low_rpm_2samples_truth.csv"),
    ("tests/test_files/counted.csv", 60, 1, True, "tests/test_files/counted_below60_rpm.csv"),
])
def test_filter_low_reads(input_path, threshold, n_samples, opposite, truth_path):
    truth = io.load_table(truth_path)
    h = CountFilter(input_path)
    h.filter_low_reads(threshold, n_samples, opposite=opposite)
    assert np.isclose(truth.sort(pl.first()).drop(cs.first()), h.df.sort(pl.first()).drop(cs.first())).all()


@pytest.mark.parametrize('interactive', [True, False])
@pytest.mark.parametrize('show_cursor', [True, False])
@pytest.mark.parametrize('title,alpha,log2fc_threshold', [
    ('auto', 0.05, None),
    ('title', 0.1, 0),
    ('title', 0.001, 1)])
def test_deseqfilter_volcano_plot_api(interactive, show_cursor, title, alpha, log2fc_threshold, basic_deseqfilter):
    basic_deseqfilter.volcano_plot(alpha, log2fc_threshold, title, interactive=interactive, show_cursor=show_cursor)
    plt.close('all')


@pytest.mark.parametrize('ref_column,columns,split_plots', [
    ('auto', 'all', False),
    ('cond2', ['cond1', 'cond2'], True),
    ('cond1', 'cond2', False)
])
def test_countfilter_ma_plot_api(ref_column, columns, split_plots, basic_countfilter):
    basic_countfilter.ma_plot(ref_column, columns, split_plots)
    plt.close('all')


@pytest.mark.parametrize('args,kwargs', [
    (tuple(), dict(log2=False)),
    ((['cond1', 'cond3'],), dict(log2=True)),
])
def test_countfilter_pairplot_api(basic_countfilter, args, kwargs):
    basic_countfilter.pairplot(*args, **kwargs)
    plt.close('all')


def test_countfilter_pairplot_shows_spearman_of_numeric_columns(basic_countfilter):
    # Regression: pairplot annotated each cell with spearmanr of sample_df indexed by seaborn's
    # numeric-only axis positions, but sample_df also contains the (non-numeric) gene-index column
    # at position 0 -- so the first row/column's correlation was computed on the index column
    # (garbage under old numpy, a hard crash under numpy 2.x). Each cell must show the Spearman
    # correlation of *its own* pair of sample columns. We identify each cell's pair from its axis
    # labels (which seaborn sets to the plotted column names) -- an independent path from the fix.
    import re
    from scipy.stats import spearmanr
    fig = basic_countfilter.pairplot(show_corr=True, log2=False)
    try:
        checked = 0
        for ax in fig.axes:
            annotations = [t.get_text() for t in ax.texts if 'Spearman' in t.get_text()]
            if not annotations:
                continue
            x_col, y_col = ax.get_xlabel(), ax.get_ylabel()
            assert x_col and y_col and x_col != y_col  # a real off-diagonal sample pair, not the index
            expected = f"{spearmanr(basic_countfilter.df[:, [x_col, y_col]])[0]:.2f}"
            shown = re.search(r'Spearman.*?=\s*(-?\d+\.\d+)', annotations[0]).group(1)
            assert shown == expected, f"cell ({y_col}, {x_col}): shown {shown} != expected {expected}"
            checked += 1
        assert checked >= 1  # at least one off-diagonal cell was actually verified
    finally:
        plt.close('all')


@pytest.mark.parametrize('args,kwargs,xfail', [
    (tuple(), dict(), False),
    ((['cond1', 'cond2'],), dict(metric='euclidean', linkage='Ward'), False),
    ((['cond1', 'cond2'],), dict(metric='euclidean', linkage='single'), False),
    (tuple(), dict(linkage='invalid'), True),
    (tuple(), dict(metric='invalid'), True),
    (tuple(), dict(linkage=5), True),

])
def test_countfilter_clustergram_api(basic_countfilter, args, kwargs, xfail):
    try:
        if xfail:
            with pytest.raises(RNAlysisInputError):
                basic_countfilter.clustergram(*args, **kwargs)
        else:
            basic_countfilter.clustergram(*args, **kwargs)
    finally:
        plt.close('all')


def test_countfilter_box_plot_api(basic_countfilter):
    basic_countfilter.enhanced_box_plot(ylabel='A different label')
    basic_countfilter.enhanced_box_plot(samples=['cond1', 'cond3'], scatter=True)
    plt.close('all')


@pytest.mark.parametrize('highlight',
                         [None, {'WBGene00007063', 'WBGene00007064'}, DESeqFilter('tests/test_files/test_deseq.csv')])
@pytest.mark.parametrize('interactive', [True, False])
@pytest.mark.parametrize('show_cursor', [True, False])
@pytest.mark.parametrize('s1,s2,xlabel,ylabel,title', [
    ('cond1', 'cond2', 'auto', 'auto', 'auto'),
    ('cond3', ['cond2', 'cond1', 'cond4'], 'x', 'y', 'title')])
def test_countfilter_scatter_sample_vs_sample_api(basic_countfilter, s1, s2, xlabel, ylabel, title, interactive,
                                                  show_cursor, highlight):
    basic_countfilter.scatter_sample_vs_sample(s1, s2, xlabel, ylabel, title, interactive=interactive,
                                               show_cursor=show_cursor, highlight=highlight)
    plt.close('all')


@pytest.mark.parametrize('kwargs,xfail', [({}, False),
                                          (dict(samples=['cond1', ['cond2', 'cond3']], n_components=2, labels=False,
                                                power_transform=False), False),
                                          (dict(n_components=2.0), True),
                                          (dict(n_components=1), True)])
def test_countfilter_pca_api(kwargs, xfail, basic_countfilter):
    basic_countfilter.filter_low_reads(1)
    try:
        if xfail:
            with pytest.raises(RNAlysisInputError):
                basic_countfilter.pca(**kwargs)

        else:
            basic_countfilter.pca()
    finally:
        plt.close('all')


def test_countfilter_enhanced_box_plot_api(basic_countfilter):
    basic_countfilter.box_plot(notch=True, ylabel='A different label')
    basic_countfilter.box_plot(samples=['cond1', 'cond3'], scatter=True)
    plt.close('all')


def test_countfilter_violin_plot_api(basic_countfilter):
    basic_countfilter.violin_plot(ylabel='A different label')
    basic_countfilter.violin_plot(samples=['cond1', 'cond4'])
    plt.close('all')


def _filter_biotype_tester(filter_obj, truth_protein_coding, truth_pirna, from_table: bool = True, **kwargs):
    if from_table:
        protein_coding = filter_obj.filter_biotype_from_ref_table(ref=__biotype_ref__, inplace=False, **kwargs)
        pirna = filter_obj.filter_biotype_from_ref_table('piRNA', ref=__biotype_ref__, inplace=False, **kwargs)
    else:
        gtf_path = 'tests/test_files/test_gtf_for_biotypes.gtf'
        protein_coding = filter_obj.filter_biotype_from_gtf(gtf_path, inplace=False, **kwargs)
        pirna = filter_obj.filter_biotype_from_gtf(gtf_path, 'piRNA', inplace=False, **kwargs)
    assert np.all(truth_protein_coding.sort(pl.first()) == protein_coding.df.sort(pl.first()))
    assert np.all(truth_pirna.sort(pl.first()) == pirna.df.sort(pl.first()))


def test_countfilter_filter_biotype_from_ref_table():
    truth_protein_coding = io.load_table('tests/test_files/counted_biotype_protein_coding.csv')
    truth_pirna = io.load_table('tests/test_files/counted_biotype_piRNA.csv')
    h = CountFilter("tests/test_files/counted_biotype.csv")
    _filter_biotype_tester(h, truth_protein_coding=truth_protein_coding, truth_pirna=truth_pirna)


def test_countfilter_filter_biotype_from_ref_table_opposite():
    truth_no_pc = io.load_table(r'tests/test_files/counted_biotype_no_protein_coding.csv')
    truth_no_pirna = io.load_table(r'tests/test_files/counted_biotype_no_piRNA.csv')
    h = CountFilter("tests/test_files/counted_biotype.csv")
    _filter_biotype_tester(h, truth_protein_coding=truth_no_pc, truth_pirna=truth_no_pirna, opposite=True)


def test_countfilter_filter_biotype_from_gtf():
    truth_protein_coding = io.load_table('tests/test_files/counted_biotype_protein_coding.csv')
    truth_pirna = io.load_table('tests/test_files/counted_biotype_piRNA.csv')
    h = CountFilter("tests/test_files/counted_biotype.csv")
    _filter_biotype_tester(h, truth_protein_coding=truth_protein_coding, truth_pirna=truth_pirna, from_table=False)


def test_countfilter_filter_biotype_from_gtf_opposite():
    truth_no_pc = io.load_table(r'tests/test_files/counted_biotype_no_protein_coding.csv')
    truth_no_pirna = io.load_table(r'tests/test_files/counted_biotype_no_piRNA.csv')
    h = CountFilter("tests/test_files/counted_biotype.csv")
    _filter_biotype_tester(h, truth_protein_coding=truth_no_pc, truth_pirna=truth_no_pirna, from_table=False,
                           opposite=True)


# ----------------------------------------------------------------------------------------------------------------------
# #154: generic filter/annotate by any GTF/GFF attribute -- column-9 attributes (biotype, ...) AND fixed columns
# (chromosome/source/strand). counted_attributes.csv has GENE1-4 (present in the GTF/GFF) plus GENE5 (absent).
# ----------------------------------------------------------------------------------------------------------------------
_ATTR_GTF = 'tests/test_files/test_gtf_attributes.gtf'
_ATTR_GFF3 = 'tests/test_files/test_gff3_attributes.gff3'
_ATTR_COUNTS = 'tests/test_files/counted_attributes.csv'


def test_filter_by_gtf_attribute_reproduces_biotype_special_case():
    # filter_by_gtf_attribute on the biotype attribute must reproduce the dedicated filter_biotype_from_gtf.
    gtf = 'tests/test_files/test_gtf_for_biotypes.gtf'
    counts = CountFilter('tests/test_files/counted_biotype.csv')
    generic = counts.filter_by_gtf_attribute(gtf, 'gene_biotype', 'protein_coding', inplace=False)
    special = counts.filter_biotype_from_gtf(gtf, 'protein_coding', inplace=False)
    assert np.all(generic.df.sort(pl.first()) == special.df.sort(pl.first()))


@pytest.mark.parametrize('path', [_ATTR_GTF, _ATTR_GFF3])
@pytest.mark.parametrize('attribute,value,truth', [
    ('chromosome', 'chr2', {'GENE2', 'GENE4'}),
    ('chromosome', ['chr1', 'chrX'], {'GENE1', 'GENE3'}),
    ('strand', '+', {'GENE1', 'GENE3'}),
    ('source', 'havana', {'GENE2', 'GENE4'}),
    ('gene_biotype', 'protein_coding', {'GENE1', 'GENE3', 'GENE4'}),
])
def test_filter_by_gtf_attribute(path, attribute, value, truth):
    counts = CountFilter(_ATTR_COUNTS)
    res = counts.filter_by_gtf_attribute(path, attribute, value, inplace=False)
    assert parsing.data_to_set(res.df.select(pl.first())) == truth


def test_filter_by_gtf_attribute_opposite_includes_unmapped():
    # opposite returns every table feature NOT kept, including GENE5 which is absent from the annotation file.
    counts = CountFilter(_ATTR_COUNTS)
    res = counts.filter_by_gtf_attribute(_ATTR_GTF, 'chromosome', 'chr2', opposite=True, inplace=False)
    assert parsing.data_to_set(res.df.select(pl.first())) == {'GENE1', 'GENE3', 'GENE5'}


def _annotation_mapping(filter_obj, column_name):
    return dict(zip(filter_obj.df.select(pl.first()).to_series().to_list(), filter_obj.df[column_name].to_list()))


@pytest.mark.parametrize('path', [_ATTR_GTF, _ATTR_GFF3])
def test_annotate_from_gtf_biotype(path):
    counts = CountFilter(_ATTR_COUNTS)
    res = counts.annotate_from_gtf(path, 'gene_biotype', inplace=False)
    assert 'gene_biotype' in res.df.columns
    # numeric columns of the original table are preserved
    assert {'cond1', 'cond2'}.issubset(set(res.df.columns))
    assert _annotation_mapping(res, 'gene_biotype') == {
        'GENE1': 'protein_coding', 'GENE2': 'lincRNA', 'GENE3': 'protein_coding',
        'GENE4': 'protein_coding', 'GENE5': None}  # GENE5 is absent from the annotation -> null


def test_annotate_from_gtf_fixed_column_custom_name():
    counts = CountFilter(_ATTR_COUNTS)
    res = counts.annotate_from_gtf(_ATTR_GTF, 'chromosome', column_name='chrom', inplace=False)
    assert 'chrom' in res.df.columns
    assert _annotation_mapping(res, 'chrom') == {
        'GENE1': 'chr1', 'GENE2': 'chr2', 'GENE3': 'chrX', 'GENE4': 'chr2', 'GENE5': None}


def test_annotate_from_gtf_inplace_default_column_name():
    counts = CountFilter(_ATTR_COUNTS)
    assert counts.annotate_from_gtf(_ATTR_GFF3, 'chromosome', inplace=True) is None
    assert _annotation_mapping(counts, 'chromosome') == {
        'GENE1': 'chr1', 'GENE2': 'chr2', 'GENE3': 'chrX', 'GENE4': 'chr2', 'GENE5': None}


def test_annotate_from_gtf_overwrites_existing_column():
    counts = CountFilter(_ATTR_COUNTS)
    res = counts.annotate_from_gtf(_ATTR_GTF, 'chromosome', column_name='cond2', inplace=False)
    assert res.df.columns.count('cond2') == 1  # overwritten, not duplicated
    assert _annotation_mapping(res, 'cond2') == {
        'GENE1': 'chr1', 'GENE2': 'chr2', 'GENE3': 'chrX', 'GENE4': 'chr2', 'GENE5': None}


def test_annotate_from_gtf_rejects_feature_id_column_name():
    # naming the annotation column after the table's feature-id (index) column must not silently drop the IDs
    counts = Filter.from_dataframe(pl.DataFrame({'gene': ['GENE1', 'GENE2'], 'cond1': [1.0, 2.0]}), 'annot_test')
    with pytest.raises(InvalidValueError):
        counts.annotate_from_gtf(_ATTR_GTF, 'chromosome', column_name='gene')


def test_filter_by_attribute(basic_deseqfilter):
    truth = io.load_table('tests/test_files/test_deseq_filter_by_attr1.csv').sort(pl.first())
    d_notinplace = basic_deseqfilter.filter_by_attribute('attribute1', ref=__attr_ref__, inplace=False)
    assert np.all(truth == d_notinplace.df.sort(pl.first()))
    basic_deseqfilter.filter_by_attribute('attribute1', ref=__attr_ref__)
    assert np.all(truth == basic_deseqfilter.df.sort(pl.first()))


def test_filter_by_attribute_from_string(monkeypatch, basic_deseqfilter):
    monkeypatch.setattr('builtins.input', lambda x: 'attribute1\nattribute2\n')
    union_truth = io.load_table('tests/test_files/counted_filter_by_bigtable_union_truth.csv')
    h = CountFilter('tests/test_files/counted_filter_by_bigtable.csv')
    assert np.all(union_truth.sort(pl.first()) == h.filter_by_attribute(mode='union',
                                                                        ref=__attr_ref__,
                                                                        inplace=False).df.sort(pl.first()))

    monkeypatch.setattr('builtins.input', lambda x: 'attribute1\nattribute2')
    assert np.all(union_truth.sort(pl.first()) == h.filter_by_attribute(mode='union',
                                                                        ref=__attr_ref__,
                                                                        inplace=False).df.sort(pl.first()))

    monkeypatch.setattr('builtins.input', lambda x: 'attribute1')
    deseq_truth = io.load_table('tests/test_files/test_deseq_filter_by_attr1.csv')
    assert np.all(
        deseq_truth.sort(pl.first()) == basic_deseqfilter.filter_by_attribute(ref=__attr_ref__, inplace=False).df.sort(
            pl.first()))


def test_filter_by_attribute_union():
    union_truth = io.load_table('tests/test_files/counted_filter_by_bigtable_union_truth.csv')
    h = CountFilter('tests/test_files/counted_filter_by_bigtable.csv')
    union = h.filter_by_attribute(['attribute1', 'attribute2'], mode='union',
                                  ref=__attr_ref__, inplace=False)
    assert np.all(union.df.sort(pl.first()) == union_truth.sort(pl.first()))


def test_filter_by_attribute_intersection():
    intersection_truth = io.load_table(r'tests/test_files/counted_filter_by_bigtable_intersect_truth.csv')
    h = CountFilter('tests/test_files/counted_filter_by_bigtable.csv')
    intersection = h.filter_by_attribute(['attribute1', 'attribute2'], mode='intersection',
                                         ref=__attr_ref__,
                                         inplace=False)
    assert np.all(intersection.df.sort(pl.first()) == intersection_truth.sort(pl.first()))


def test_filter_by_attribute_invalid_mode():
    h = CountFilter('tests/test_files/counted_filter_by_bigtable.csv')
    with pytest.raises(InvalidValueError):
        h.filter_by_attribute(['attribute1', 'attribute2'], mode='difference',
                              ref=__attr_ref__)


def test_split_by_attribute():
    h = CountFilter('tests/test_files/counted_filter_by_bigtable.csv')
    attrs = ['attribute2', 'attribute3', 'attribute4', 'attribute1']
    newobjs = h.split_by_attribute(attrs, ref=__attr_ref__)
    assert len(newobjs) == len(attrs)
    for i, attr in enumerate(attrs):
        assert np.all(
            newobjs[i].df.sort(pl.first()) == h.filter_by_attribute(attr, ref=__attr_ref__,
                                                                    inplace=False).df.sort(pl.first()))


def test_split_by_attribute_multiple(basic_deseqfilter):
    attrs = ['attribute2', 'attribute3', 'attribute4', 'attribute1']
    newobjs = basic_deseqfilter.split_by_attribute(attrs, ref=__attr_ref__)
    assert len(newobjs) == len(attrs)
    for i, attr in enumerate(attrs):
        assert np.all(
            newobjs[i].df.sort(pl.first()) == basic_deseqfilter.filter_by_attribute(attr, ref=__attr_ref__,
                                                                                    inplace=False).df.sort(pl.first()))


def test_split_by_attribute_only_one_attribute(basic_deseqfilter):
    newobj = basic_deseqfilter.split_by_attribute(['attribute1'], ref=__attr_ref__)
    assert len(newobj) == 1
    assert np.all(
        newobj[0].df.sort(pl.first()) == basic_deseqfilter.filter_by_attribute('attribute1', ref=__attr_ref__,
                                                                               inplace=False).df.sort(pl.first()))
    with pytest.raises(InvalidTypeError):
        basic_deseqfilter.split_by_attribute('attribute1', ref=__attr_ref__)


def test_split_by_attribute_faulty_attributes(basic_filter):
    with pytest.raises(InvalidTypeError):
        basic_filter.split_by_attribute(['attribute1', ['attribute2', 'attribute3']],
                                        ref=__attr_ref__)
    with pytest.raises(InvalidTypeError):
        basic_filter.split_by_attribute(['attribute1', 2], ref=__attr_ref__)


def test_deseq_filter_significant():
    deseqfilter = DESeqFilter('tests/test_files/test_deseq_sig.csv')
    truth = io.load_table("tests/test_files/test_deseq_sig_truth.csv").sort(pl.first())
    result = deseqfilter.filter_significant(alpha=0.05, inplace=False)
    assert result.df.sort(pl.first()).equals(truth)


def test_deseq_filter_significant_opposite():
    deseqfilter = DESeqFilter('tests/test_files/test_deseq_sig.csv')
    truth = io.load_table(r'tests/test_files/test_deseq_not_sig_truth.csv').sort(pl.first())
    deseqfilter.filter_significant(alpha=0.05, opposite=True)
    assert deseqfilter.df.sort(pl.first()).equals(truth)


def test_filter_top_n_ascending_number(basic_deseqfilter):
    truth = io.load_table("tests/test_files/test_deseq_top10.csv").sort(pl.first())
    basic_deseqfilter.filter_top_n('padj', 10)
    assert np.isclose(truth.drop(cs.first()), basic_deseqfilter.df.sort(pl.first()).drop(cs.first())).all()


def test_filter_top_n_ascending_text():
    deseqfilter = DESeqFilter("tests/test_files/test_deseq_textcol.csv")
    truth = io.load_table("tests/test_files/test_deseq_top10_text_ascend.csv").sort(pl.first())
    deseqfilter.filter_top_n('textcol', 10, True)
    assert deseqfilter.df.sort(pl.first()).equals(truth)


def test_filter_top_n_multiple_columns():
    deseqfilter = DESeqFilter("tests/test_files/test_deseq_textcol.csv")
    truth = io.load_table("tests/test_files/test_deseq_textcol_top15_text_basemean.csv").sort(pl.first())
    deseqfilter.filter_top_n(['textcol', 'baseMean'], 15, True)
    assert deseqfilter.df.sort(pl.first()).equals(truth)


def test_filter_top_n_descending_number(basic_deseqfilter):
    truth = io.load_table("tests/test_files/test_deseq_bottom7.csv").sort(pl.first())
    basic_deseqfilter.filter_top_n('log2FoldChange', 7, False)
    assert np.isclose(truth.drop(cs.first()), basic_deseqfilter.df.sort(pl.first()).drop(cs.first())).all()


def test_filter_top_n_descending_text():
    deseqfilter = DESeqFilter("tests/test_files/test_deseq_textcol.csv")
    truth = io.load_table("tests/test_files/test_deseq_bottom10_text_descend.csv").sort(pl.first())
    deseqfilter.filter_top_n('textcol', 10, False)
    assert np.all(truth == deseqfilter.df.sort(pl.first()))


def test_filter_top_n_nonexisting_column(basic_deseqfilter):
    colname = 'somecol'
    with pytest.raises(InvalidValueError):
        basic_deseqfilter.filter_top_n(colname, 5)
    with pytest.raises(InvalidValueError):
        basic_deseqfilter.filter_top_n([basic_deseqfilter.df.columns[0], colname])
    assert colname not in basic_deseqfilter.df.columns


def test_deseq_filter_abs_log2_fold_change():
    truth = io.load_table("tests/test_files/test_deseq_fc_4_truth.csv").sort(pl.first())
    d = DESeqFilter("tests/test_files/test_deseq_fc.csv")
    fc4 = d.filter_abs_log2_fold_change(4, inplace=False)
    assert np.all(fc4.df.sort(pl.first()) == truth)


def test_deseq_filter_fold_change_direction():
    pos_truth = io.load_table("tests/test_files/test_deseq_fc_pos_truth.csv")
    neg_truth = io.load_table("tests/test_files/test_deseq_fc_neg_truth.csv")
    d = DESeqFilter("tests/test_files/test_deseq_fc.csv")
    pos = d.filter_fold_change_direction('pos', inplace=False)
    neg = d.filter_fold_change_direction('neg', inplace=False)
    assert np.all(pos.df == pos_truth)
    assert np.all(neg.df == neg_truth)


def test_deseq_split_fold_change():
    d = DESeqFilter("tests/test_files/test_deseq_fc.csv")
    pos_truth = io.load_table("tests/test_files/test_deseq_fc_pos_truth.csv")
    neg_truth = io.load_table("tests/test_files/test_deseq_fc_neg_truth.csv")
    d = DESeqFilter("tests/test_files/test_deseq_fc.csv")
    pos, neg = d.split_fold_change_direction()
    assert np.all(pos.df == pos_truth)
    assert np.all(neg.df == neg_truth)


def test_intersection():
    intersection_truth = {'WBGene00021375', 'WBGene00044258', 'WBGene00219304', 'WBGene00194708', 'WBGene00018199',
                          'WBGene00019174', 'WBGene00021019', 'WBGene00013816', 'WBGene00045366', 'WBGene00219307',
                          'WBGene00045410', 'WBGene00010100', 'WBGene00077437', 'WBGene00007674', 'WBGene00023036',
                          'WBGene00012648', 'WBGene00022486'}
    set1 = DESeqFilter('tests/test_files/test_deseq_set_ops_1.csv')
    set2 = DESeqFilter('tests/test_files/test_deseq_set_ops_2.csv')

    assert set1.intersection(set2, inplace=False) == intersection_truth


@pytest.mark.parametrize("this_set,other_sets,majority_threshold,truth",
                         [(Filter("tests/test_files/test_deseq.csv"),
                           [{'WBGene00000001', 'WBGene00000002', 'WBGene00000003'},
                            {'WBGene00000002', 'WBGene00000004'}], 2 / 3,
                           {'WBGene00000002', 'WBGene00000003', 'WBGene00000004'}),
                          (DESeqFilter('tests/test_files/test_deseq_set_ops_1.csv'),
                           [DESeqFilter('tests/test_files/test_deseq_set_ops_2.csv'),
                            {'WBGene00008447', 'WBGene00021018'}], 2 / 3, {'WBGene00008447',
                                                                           'WBGene00013816',
                                                                           'WBGene00018199',
                                                                           'WBGene00019174',
                                                                           'WBGene00021018',
                                                                           'WBGene00021019',
                                                                           'WBGene00021375',
                                                                           'WBGene00044258',
                                                                           'WBGene00045366',
                                                                           'WBGene00045410',
                                                                           'WBGene00194708',
                                                                           'WBGene00219304',
                                                                           'WBGene00219307',
                                                                           'WBGene00010100',
                                                                           'WBGene00077437',
                                                                           'WBGene00023036',
                                                                           'WBGene00012648',
                                                                           'WBGene00022486',
                                                                           'WBGene00007674'}),
                          (DESeqFilter('tests/test_files/test_deseq_set_ops_1.csv'),
                           [DESeqFilter('tests/test_files/test_deseq_set_ops_2.csv'),
                            {'WBGene00008447', 'WBGene00021018'}], 1, set())])
def test_majority_vote_intersection(this_set, other_sets, majority_threshold, truth):
    result = this_set.majority_vote_intersection(*other_sets, majority_threshold=majority_threshold)
    assert result == truth


def test_union():
    intersection_truth = {'WBGene00021375', 'WBGene00044258', 'WBGene00219304', 'WBGene00194708', 'WBGene00018199',
                          'WBGene00019174', 'WBGene00021019', 'WBGene00013816', 'WBGene00045366', 'WBGene00219307',
                          'WBGene00045410', 'WBGene00010100', 'WBGene00077437', 'WBGene00007674', 'WBGene00023036',
                          'WBGene00012648', 'WBGene00022486'}
    set2_unique = {'WBGene00018193', 'WBGene00021589', 'WBGene00001118', 'WBGene00010755', 'WBGene00020407',
                   'WBGene00044799', 'WBGene00021654', 'WBGene00012919', 'WBGene00021605'}
    set1_unique = {'WBGene00008447', 'WBGene00021018', 'WBGene00012452', 'WBGene00010507', 'WBGene00022730',
                   'WBGene00012961', 'WBGene00022438', 'WBGene00016635', 'WBGene00044478'}

    set1 = DESeqFilter('tests/test_files/test_deseq_set_ops_1.csv')
    set2 = DESeqFilter('tests/test_files/test_deseq_set_ops_2.csv')
    union_truth = intersection_truth.union(set1_unique.union(set2_unique))
    assert set1.union(set2) == union_truth


def test_difference():
    set2_unique = {'WBGene00018193', 'WBGene00021589', 'WBGene00001118', 'WBGene00010755', 'WBGene00020407',
                   'WBGene00044799', 'WBGene00021654', 'WBGene00012919', 'WBGene00021605'}
    set1_unique = {'WBGene00008447', 'WBGene00021018', 'WBGene00012452', 'WBGene00010507', 'WBGene00022730',
                   'WBGene00012961', 'WBGene00022438', 'WBGene00016635', 'WBGene00044478'}

    set1 = DESeqFilter('tests/test_files/test_deseq_set_ops_1.csv')
    set2 = DESeqFilter('tests/test_files/test_deseq_set_ops_2.csv')

    assert set1.difference(set2, inplace=False) == set1_unique
    assert set2.difference(set1, inplace=False) == set2_unique


def test_symmetric_difference():
    set2_unique = {'WBGene00018193', 'WBGene00021589', 'WBGene00001118', 'WBGene00010755', 'WBGene00020407',
                   'WBGene00044799', 'WBGene00021654', 'WBGene00012919', 'WBGene00021605'}
    set1_unique = {'WBGene00008447', 'WBGene00021018', 'WBGene00012452', 'WBGene00010507', 'WBGene00022730',
                   'WBGene00012961', 'WBGene00022438', 'WBGene00016635', 'WBGene00044478'}

    set1 = DESeqFilter('tests/test_files/test_deseq_set_ops_1.csv')
    set2 = DESeqFilter('tests/test_files/test_deseq_set_ops_2.csv')

    assert set1.symmetric_difference(set2) == set2.symmetric_difference(set1)
    assert set1.symmetric_difference(set2) == set1_unique.union(set2_unique)


def test_set_ops_symmetric_difference_more_than_two_objects():
    set1 = DESeqFilter('tests/test_files/test_deseq_set_ops_1.csv')
    set2 = DESeqFilter('tests/test_files/test_deseq_set_ops_2.csv')
    set3 = {1, 2, 3}
    with pytest.raises(TypeError):
        _ = set1._set_ops([set2, set3], 'set', set.symmetric_difference)


def test_deseq_feature_set():
    truth = {'WBGene00008447', 'WBGene00021018', 'WBGene00012452', 'WBGene00010507', 'WBGene00022730',
             'WBGene00012648',
             'WBGene00012961', 'WBGene00022438', 'WBGene00016635', 'WBGene00044478', 'WBGene00021375',
             'WBGene00044258', 'WBGene00219304', 'WBGene00194708', 'WBGene00018199', 'WBGene00022486',
             'WBGene00019174', 'WBGene00021019', 'WBGene00013816', 'WBGene00045366', 'WBGene00219307',
             'WBGene00045410', 'WBGene00010100', 'WBGene00077437', 'WBGene00007674', 'WBGene00023036'}
    d = DESeqFilter('tests/test_files/test_deseq_set_ops_1.csv')
    assert d.index_set == truth


def test_deseq_feature_string():
    truth = {'WBGene00008447', 'WBGene00021018', 'WBGene00012452', 'WBGene00010507', 'WBGene00022730',
             'WBGene00012648',
             'WBGene00012961', 'WBGene00022438', 'WBGene00016635', 'WBGene00044478', 'WBGene00021375',
             'WBGene00044258', 'WBGene00219304', 'WBGene00194708', 'WBGene00018199', 'WBGene00022486',
             'WBGene00019174', 'WBGene00021019', 'WBGene00013816', 'WBGene00045366', 'WBGene00219307',
             'WBGene00045410', 'WBGene00010100', 'WBGene00077437', 'WBGene00007674', 'WBGene00023036'}
    d = DESeqFilter('tests/test_files/test_deseq_set_ops_1.csv')
    assert set(d.index_string.split("\n")) == truth


def test_set_ops_multiple_variable_types():
    set2_unique = {'WBGene00018193', 'WBGene00021589', 'WBGene00001118', 'WBGene00010755', 'WBGene00020407',
                   'WBGene00044799', 'WBGene00021654', 'WBGene00012919', 'WBGene00021605'}
    set1_unique = {'WBGene00008447', 'WBGene00021018', 'WBGene00012452', 'WBGene00010507', 'WBGene00022730',
                   'WBGene00012961', 'WBGene00022438', 'WBGene00016635', 'WBGene00044478'}

    set1 = CountFilter('tests/test_files/test_deseq_set_ops_1.csv')
    set2 = DESeqFilter('tests/test_files/test_deseq_set_ops_2.csv')

    assert set1.symmetric_difference(set2) == set2.symmetric_difference(set1)
    assert set1.symmetric_difference(set2) == set1_unique.union(set2_unique)


def test_countfilter_rpm_negative_threshold(basic_countfilter):
    with pytest.raises(InvalidValueError):
        basic_countfilter.filter_low_reads(threshold=-3)


def test_countfilter_threshold_invalid(basic_countfilter):
    with pytest.raises(InvalidTypeError):
        basic_countfilter.filter_low_reads("5")


def test_countfilter_split_by_reads(basic_countfilter):
    high_truth = io.load_table(r"tests/test_files/counted_above60_rpm.csv")
    low_truth = io.load_table(r"tests/test_files/counted_below60_rpm.csv")
    high, low = basic_countfilter.split_by_reads(threshold=60)
    assert high.df.sort(pl.first()).equals(high_truth.sort(pl.first()))
    assert low.df.sort(pl.first()).equals(low_truth.sort(pl.first()))


def test_filter_percentile():
    truth = io.load_table(r'tests/test_files/test_deseq_percentile_0.25.csv').sort(pl.first())
    h = DESeqFilter(r'tests/test_files/test_deseq_percentile.csv')
    h.filter_percentile(0.25, 'padj', inplace=True)
    assert truth.equals(h.df.sort(pl.first()))
    h.filter_percentile(1, 'baseMean')
    assert truth.equals(h.df.sort(pl.first()))
    h.filter_percentile(0, 'padj')
    assert len(h) == 1
    assert (h.df['padj'] == truth['padj'].min()).all()


def test_filter_percentile_bad_input():
    h = DESeqFilter(r'tests/test_files/test_deseq_percentile.csv')
    with pytest.raises(InvalidValueError):
        h.filter_percentile(-0.2, 'pvalue')
    with pytest.raises(InvalidValueError):
        h.filter_percentile(1.1, 'baseMean')
    with pytest.raises(InvalidTypeError):
        h.filter_percentile('0.5', 'log2FoldChange')


def test_split_by_percentile():
    truth_below = io.load_table(r'tests/test_files/test_deseq_percentile_0.25.csv').sort(pl.first())
    truth_above = io.load_table(r'tests/test_files/test_deseq_percentile_0.75.csv').sort(pl.first())
    h = DESeqFilter(r'tests/test_files/test_deseq_percentile.csv')
    below, above = h.split_by_percentile(0.25, 'padj')
    assert np.all(truth_below == below.df.sort(pl.first()))
    assert np.all(truth_above == above.df.sort(pl.first()))


def test_countfilter_filter_biotype_from_ref_table_multiple():
    truth = io.load_table('tests/test_files/counted_biotype_piRNA_protein_coding.csv').sort(pl.first())
    h = CountFilter("tests/test_files/counted_biotype.csv")
    both = h.filter_biotype_from_ref_table(['protein_coding', 'piRNA'], ref=__biotype_ref__, inplace=False)
    assert both.df.sort(pl.first()).equals(truth)


def test_countfilter_filter_biotype_from_ref_table_multiple_opposite():
    truth = io.load_table('tests/test_files/counted_biotype_piRNA_protein_coding_opposite.csv').sort(pl.first())
    h = CountFilter("tests/test_files/counted_biotype.csv")
    neither = h.filter_biotype_from_ref_table(['protein_coding', 'piRNA'], ref=__biotype_ref__,
                                              inplace=False, opposite=True)
    assert np.all(truth == neither.df.sort(pl.first()))


def test_deseq_filter_biotype_from_ref_table():
    truth_protein_coding = io.load_table('tests/test_files/test_deseq_biotype_protein_coding.csv')
    truth_pirna = io.load_table('tests/test_files/test_deseq_biotype_piRNA.csv')
    d = DESeqFilter("tests/test_files/test_deseq_biotype.csv")
    _filter_biotype_tester(d, truth_protein_coding=truth_protein_coding, truth_pirna=truth_pirna)


def test_deseq_filter_biotype_from_ref_table_opposite():
    truth_no_pirna = io.load_table(r'tests/test_files/test_deseq_biotype_piRNA_opposite.csv').sort(pl.first())
    d = DESeqFilter("tests/test_files/test_deseq_biotype.csv")
    d.filter_biotype_from_ref_table('piRNA', ref=__biotype_ref__, opposite=True, inplace=True)
    assert np.all(d.df.sort(pl.first()) == truth_no_pirna)


def test_deseq_filter_biotype_from_ref_table_multiple():
    truth = io.load_table('tests/test_files/test_deseq_biotype_piRNA_protein_coding.csv').sort(pl.first())
    d = DESeqFilter("tests/test_files/test_deseq_biotype.csv")
    both = d.filter_biotype_from_ref_table(['protein_coding', 'piRNA'], ref=__biotype_ref__,
                                           inplace=False)
    assert np.all(truth == both.df.sort(pl.first()))


def test_deseq_filter_biotype_from_ref_table_multiple_opposite():
    truth = io.load_table('tests/test_files/test_deseq_biotype_piRNA_protein_coding_opposite.csv').sort(pl.first())
    d = DESeqFilter("tests/test_files/test_deseq_biotype.csv")
    neither = d.filter_biotype_from_ref_table(['protein_coding', 'piRNA'], ref=__biotype_ref__,
                                              inplace=False,
                                              opposite=True)
    assert np.all(truth == neither.df.sort(pl.first()))


def test_deseqfilter_union_multiple():
    intersection_truth = {'WBGene00021375', 'WBGene00044258', 'WBGene00219304', 'WBGene00194708', 'WBGene00018199',
                          'WBGene00019174', 'WBGene00021019', 'WBGene00013816', 'WBGene00045366', 'WBGene00219307',
                          'WBGene00045410', 'WBGene00010100', 'WBGene00077437', 'WBGene00007674', 'WBGene00023036',
                          'WBGene00012648', 'WBGene00022486'}
    set2_unique = {'WBGene00018193', 'WBGene00021589', 'WBGene00001118', 'WBGene00010755', 'WBGene00020407',
                   'WBGene00044799', 'WBGene00021654', 'WBGene00012919', 'WBGene00021605'}
    set1_unique = {'WBGene00008447', 'WBGene00021018', 'WBGene00012452', 'WBGene00010507', 'WBGene00022730',
                   'WBGene00012961', 'WBGene00022438', 'WBGene00016635', 'WBGene00044478'}
    set3_unique = {'WBGene44444444', 'WBGene99999999', 'WBGene98765432'}

    set1 = DESeqFilter('tests/test_files/test_deseq_set_ops_1.csv')
    set2 = DESeqFilter('tests/test_files/test_deseq_set_ops_2.csv')
    set3 = {'WBGene00077437', 'WBGene00007674', 'WBGene00023036', 'WBGene00012648', 'WBGene44444444',
            'WBGene99999999',
            'WBGene98765432'}
    union_truth = intersection_truth.union(set1_unique, set2_unique, set3_unique)
    assert set1.union(set2, set3) == union_truth


def test_deseqfilter_intersection_multiple():
    intersection_truth = {'WBGene00077437', 'WBGene00007674', 'WBGene00023036',
                          'WBGene00012648', 'WBGene00022486'}
    set1 = DESeqFilter('tests/test_files/test_deseq_set_ops_1.csv')
    set2 = DESeqFilter('tests/test_files/test_deseq_set_ops_2.csv')
    set3 = {'WBGene00077437', 'WBGene00007674', 'WBGene00023036', 'WBGene00012648', 'WBGene00022486',
            'WBGene99999999',
            'WBGene98765432'}

    assert set1.intersection(set2, set3, inplace=False) == intersection_truth


def test_deseqfilter_difference_multiple():
    set2_unique = {'WBGene00021589', 'WBGene00001118', 'WBGene00010755', 'WBGene00020407',
                   'WBGene00044799', 'WBGene00021654', 'WBGene00012919', 'WBGene00021605'}
    set1_unique = {'WBGene00021018', 'WBGene00012452', 'WBGene00010507', 'WBGene00022730',
                   'WBGene00012961', 'WBGene00022438', 'WBGene00016635', 'WBGene00044478'}

    set1 = DESeqFilter('tests/test_files/test_deseq_set_ops_1.csv')
    set2 = DESeqFilter('tests/test_files/test_deseq_set_ops_2.csv')
    set3 = {'WBGene00018193', 'WBGene00008447', 'WBGene12345678'}

    assert set1.difference(set2, set3, inplace=False) == set1_unique
    assert set2.difference(set3, set1, inplace=False) == set2_unique


def test_intersection_inplace():
    set1_truth = io.load_table('tests/test_files/test_deseq_set_ops_1_inplace_intersection.csv').sort(pl.first())
    set2_truth = io.load_table('tests/test_files/test_deseq_set_ops_2_inplace_intersection.csv').sort(pl.first())
    set1 = DESeqFilter('tests/test_files/test_deseq_set_ops_1.csv')
    set2 = DESeqFilter('tests/test_files/test_deseq_set_ops_2.csv')
    set1_int = set1.__copy__()
    set2_int = set2.__copy__()
    set1_int.intersection(set2, inplace=True)
    set2_int.intersection(set1, inplace=True)
    assert np.all(set1_truth == set1_int.df.sort(pl.first()))
    assert np.all(set2_truth == set2_int.df.sort(pl.first()))


def test_difference_inplace():
    set1_truth = io.load_table('tests/test_files/test_deseq_set_ops_1_inplace_difference.csv').sort(pl.first())
    set2_truth = io.load_table('tests/test_files/test_deseq_set_ops_2_inplace_difference.csv').sort(pl.first())
    set1 = DESeqFilter('tests/test_files/test_deseq_set_ops_1.csv')
    set2 = DESeqFilter('tests/test_files/test_deseq_set_ops_2.csv')
    set1_diff = set1.__copy__()
    set2_diff = set2.__copy__()
    set1_diff.difference(set2, inplace=True)
    set2_diff.difference(set1, inplace=True)
    assert np.all(set1_truth == set1_diff.df.sort(pl.first()))
    assert np.all(set2_truth == set2_diff.df.sort(pl.first()))


def test_countfilter_fold_change():
    truth_num_name = f"Mean of {['cond1_rep1', 'cond1_rep2']}"
    truth_denom_name = f"Mean of {['cond2_rep1', 'cond2_rep2']}"
    truth = io.load_table(r'tests/test_files/counted_fold_change_truth.csv')
    h = CountFilter(r'tests/test_files/counted_fold_change.csv')
    fc = h.fold_change(['cond1_rep1', 'cond1_rep2'], ['cond2_rep1', 'cond2_rep2'])
    assert truth_num_name == fc.numerator
    assert truth_denom_name == fc.denominator
    print(fc.df == truth)
    assert np.allclose(fc.df.select(pl.last()), truth.select(pl.last()))


def test_fcfilter_filter_abs_fc():
    truth = io.load_table('tests/test_files/fcfilter_abs_fold_change_truth.csv').sort(pl.first())
    f = FoldChangeFilter('tests/test_files/counted_fold_change_truth.csv', 'numer', 'denom')
    f.filter_abs_log2_fold_change(1)
    assert np.all(np.squeeze(f.df.sort(pl.first())) == np.squeeze(truth))


def test_fcfilter_fold_change_direction():
    truth_pos = io.load_table('tests/test_files/fc_1_pos_fc.csv', 0, squeeze=True)
    truth_neg = io.load_table('tests/test_files/fc_1_neg_fc.csv', 0, squeeze=True)
    fc = FoldChangeFilter('tests/test_files/fc_1.csv', 'name', 'name')
    pos = fc.filter_fold_change_direction('pos', inplace=False)
    neg = fc.filter_fold_change_direction('neg', inplace=False)
    assert truth_pos.equals(pos.df)
    assert truth_neg.equals(neg.df)


def test_fcfilter_split_fold_change_direction():
    truth_pos = io.load_table('tests/test_files/fc_1_pos_fc.csv', 0, squeeze=True)
    truth_neg = io.load_table('tests/test_files/fc_1_neg_fc.csv', 0, squeeze=True)
    fc = FoldChangeFilter('tests/test_files/fc_1.csv', 'name', 'name')
    pos, neg = fc.split_fold_change_direction()
    assert truth_pos.equals(pos.df)
    assert truth_neg.equals(neg.df)


def test_fcfilter_filter_fold_change_direction_bad_input():
    fc = FoldChangeFilter('tests/test_files/fc_1.csv', 'name', 'name')
    with pytest.raises(ValueError):
        fc.filter_fold_change_direction('bad_input')


@pytest.mark.parametrize('args', [
    ('log2FoldChange', '|x|>', 4),
    ('log2FoldChange', 'abs_gt', 4),
    ('log2FoldChange', 'abs greater than', 4)
])
def test_number_filters_absgt(basic_deseqfilter, args):
    truth = io.load_table(r'tests/test_files/test_deseq_absgt.csv').sort(pl.first())
    result = basic_deseqfilter.number_filters(*args, inplace=False)
    assert np.all(np.squeeze(truth) == np.squeeze(result.df.sort(pl.first())))


@pytest.mark.parametrize('args', [
    ('baseMean', 'greater tHAn', 1000),
    ('baseMean', 'GT', 1000),
    ('baseMean', '>', 1000)
])
def test_number_filters_gt(basic_deseqfilter, args):
    truth = io.load_table(r'tests/test_files/test_deseq_gt.csv').sort(pl.first())
    result = basic_deseqfilter.number_filters(*args, inplace=False)
    assert np.all(np.squeeze(truth) == np.squeeze(result.df.sort(pl.first())))


@pytest.mark.parametrize('args', [
    ('lfcSE', 'lt', 0.2),
    ('lfcSE', '<', 0.2),
    ('lfcSE', 'Lesser than', 0.2)
])
def test_number_filters_lt(basic_deseqfilter, args):
    truth = io.load_table(r'tests/test_files/test_deseq_lt.csv').sort(pl.first())
    result = basic_deseqfilter.number_filters(*args, inplace=False)
    assert np.all(np.squeeze(truth) == np.squeeze(result.df.sort(pl.first())))


@pytest.mark.parametrize('args', [
    ('cond2', 'eq', 0),
    ('cond2', '=', 0),
    ('cond2', 'Equals', 0)
])
def test_number_filters_eq(basic_countfilter, args):
    truth = io.load_table(r'tests/test_files/counted_eq.csv').sort(pl.first())
    result = basic_countfilter.number_filters(*args, inplace=False)
    assert np.all(np.squeeze(truth) == np.squeeze(result.df.sort(pl.first())))


@pytest.mark.parametrize('args', [
    ('Cond2', 'lt', 5),
    ('cond2', 'contains', 6),
    ('cond2', 'equals', '55')
])
def test_number_filters_invalid_input(basic_countfilter, args):
    with pytest.raises(RNAlysisInputError):
        basic_countfilter.number_filters(*args)


def test_text_filters_eq():
    truth = io.load_table('tests/test_files/text_filters_eq.csv').sort(pl.first())
    d = CountFilter('tests/test_files/text_filters.csv')
    filt_1 = d.text_filters('class', 'eQ', 'B', inplace=False)
    filt_2 = d.text_filters('class', '=', 'B', inplace=False)
    filt_3 = d.text_filters('class', 'Equals', 'B', inplace=False)
    assert np.all(filt_1.df.sort(pl.first()) == filt_2.df.sort(pl.first()))
    assert np.all(filt_2.df.sort(pl.first()) == filt_3.df.sort(pl.first()))
    assert np.all(np.squeeze(truth) == np.squeeze(filt_1.df.sort(pl.first())))


def test_text_filters_ct():
    truth = io.load_table('tests/test_files/text_filters_ct.csv').sort(pl.first())
    d = CountFilter('tests/test_files/text_filters.csv')
    filt_1 = d.text_filters('name', 'ct', 'C3.', inplace=False)
    filt_2 = d.text_filters('name', 'IN', 'C3.', inplace=False)
    filt_3 = d.text_filters('name', 'contaiNs', 'C3.', inplace=False)
    assert np.all(filt_1.df.sort(pl.first()) == filt_2.df.sort(pl.first()))
    assert np.all(filt_2.df.sort(pl.first()) == filt_3.df.sort(pl.first()))
    assert np.all(np.squeeze(truth) == np.squeeze(filt_1.df.sort(pl.first())))


def test_text_filters_sw():
    truth = io.load_table('tests/test_files/text_filters_sw.csv').sort(pl.first())
    d = CountFilter('tests/test_files/text_filters.csv')
    filt_1 = d.text_filters('name', 'sw', '2R', inplace=False)
    filt_2 = d.text_filters('name', 'Starts With', '2R', inplace=False)
    assert np.all(filt_1.df.sort(pl.first()) == filt_2.df.sort(pl.first()))
    assert np.all(np.squeeze(truth) == np.squeeze(filt_1.df.sort(pl.first())))


def test_text_filters_ew():
    truth = io.load_table('tests/test_files/text_filters_ew.csv').sort(pl.first())
    d = CountFilter('tests/test_files/text_filters.csv')
    filt_1 = d.text_filters('name', 'ew', '3', inplace=False)
    filt_2 = d.text_filters('name', 'ends With', '3', inplace=False)
    assert np.all(filt_1.df.sort(pl.first()) == filt_2.df.sort(pl.first()))
    assert np.all(np.squeeze(truth) == np.squeeze(filt_1.df.sort(pl.first())))


@pytest.mark.parametrize('args', [
    ('Cond2', 'contains', '5'),
    ('cond2', 'lt', '6'),
    ('cond2', 'equals', 55)
])
def test_text_filters_invalid_input(basic_countfilter, args):
    with pytest.raises(RNAlysisInputError):
        basic_countfilter.text_filters(*args)


class TestCountFilterFromFolder:
    def test_count_filter_from_folder(self):
        counted_fname = '__allexpr_temporary_testfile.csv'
        uncounted_fname = '__allfeature_temporary_testfile.csv'

        truth_all_expr = io.load_table('tests/test_files/test_count_from_folder_all_expr.csv').sort(pl.first())
        truth_all_feature = io.load_table('tests/test_files/test_count_from_folder_all_feature.csv').sort(pl.first())
        counts = CountFilter.from_folder_htseqcount('tests/test_files/test_count_from_folder', norm_to_rpm=False,
                                                    save_csv=True,
                                                    counted_fname=counted_fname, uncounted_fname=uncounted_fname)

        try:
            assert counts.df.sort(pl.first()).equals(truth_all_expr)

            all_feature = io.load_table(f'tests/test_files/test_count_from_folder/{uncounted_fname}').sort(pl.first())
            assert all_feature.equals(truth_all_feature)

        finally:
            os.remove('tests/test_files/test_count_from_folder/__allexpr_temporary_testfile.csv')
            os.remove('tests/test_files/test_count_from_folder/__allfeature_temporary_testfile.csv')

    def test_count_filter_from_folder_save_without_suffix(self):
        counted_fname = '__allexpr_temporary_testfile.csv'
        uncounted_fname = '__allfeature_temporary_testfile.csv'
        try:
            _ = CountFilter.from_folder_htseqcount('tests/test_files/test_count_from_folder', norm_to_rpm=False,
                                                   save_csv=True,
                                                   counted_fname=counted_fname, uncounted_fname=uncounted_fname)
        finally:
            os.remove(f'tests/test_files/test_count_from_folder/{counted_fname}')
            os.remove(f'tests/test_files/test_count_from_folder/{uncounted_fname}')

    def test_count_filter_from_folder_norm(self):
        truth_norm = io.load_table('tests/test_files/test_count_from_folder_norm.csv').drop(cs.first())
        counts_norm = CountFilter.from_folder_htseqcount('tests/test_files/test_count_from_folder',
                                                         norm_to_rpm=True, save_csv=False)
        assert np.all(np.isclose(counts_norm.df.drop(cs.first()), truth_norm, atol=0, rtol=0.0001))


def test_biotypes_from_ref_table():
    truth = io.load_table('tests/test_files/biotypes_truth.csv').sort(pl.first())
    c = CountFilter('tests/test_files/counted_biotype.csv')
    df = c.biotypes_from_ref_table(ref=__biotype_ref__).sort(pl.first())
    print(truth)
    print(df)
    assert df.equals(truth)


def test_biotypes_from_ref_table_long_form():
    truth = pl.read_csv('tests/test_files/biotypes_long_format_truth.csv').sort(pl.first())
    c = CountFilter('tests/test_files/counted_biotype.csv')
    df = c.biotypes_from_ref_table(long_format=True, ref=__biotype_ref__).sort(pl.first())
    print(df)
    print(truth)
    assert np.allclose(df.drop(cs.first()), truth.drop(cs.first()), equal_nan=True)


def test_biotypes_from_gtf():
    truth = io.load_table('tests/test_files/biotypes_truth.csv').sort(pl.first())
    c = CountFilter('tests/test_files/counted_biotype.csv')
    df = c.biotypes_from_gtf('tests/test_files/test_gtf_for_biotypes.gtf').sort(pl.first())
    assert df.equals(truth)


def test_biotypes_from_gtf_long_form():
    truth = pl.read_csv('tests/test_files/biotypes_long_format_truth.csv').sort(pl.first())
    c = CountFilter('tests/test_files/counted_biotype.csv')
    df = c.biotypes_from_gtf('tests/test_files/test_gtf_for_biotypes.gtf', long_format=True).sort(pl.first())
    assert np.allclose(df.drop(cs.first()), truth.drop(cs.first()), equal_nan=True)


def test_filter_by_row_sum(basic_countfilter):
    truth = io.load_table('tests/test_files/test_filter_row_sum.csv').sort(pl.first())
    basic_countfilter.filter_by_row_sum(29)
    assert np.all(basic_countfilter.df.sort(pl.first()) == truth)


def test_sort_inplace(basic_countfilter):
    truth = io.load_table('tests/test_files/counted_sorted.csv')
    basic_countfilter.sort(by='cond3', ascending=True, inplace=True)
    assert basic_countfilter.df.equals(truth)


def test_sort_not_inplacebasic_countfilter(basic_countfilter):
    truth = io.load_table('tests/test_files/counted_sorted.csv')
    c_copy = io.load_table('tests/test_files/counted.csv')
    c_sorted = basic_countfilter.sort(by='cond3', ascending=True, inplace=False)
    assert c_sorted.df.equals(truth)
    assert basic_countfilter.df.equals(c_copy)


def test_sort_by_multiple_columns(basic_countfilter):
    truth = io.load_table('tests/test_files/counted_sorted_multiple_truth.csv')
    basic_countfilter.sort(by=['cond3', 'cond4', 'cond1', 'cond2'], ascending=[True, False, True, False], inplace=True)
    assert np.all(truth == basic_countfilter.df)


def test_sort_with_na_first():
    truth_first = io.load_table('tests/test_files/test_deseq_with_nan_sorted_nanfirst_truth.csv')
    truth_last = io.load_table('tests/test_files/test_deseq_with_nan_sorted_nanlast_truth.csv')
    c = CountFilter('tests/test_files/test_deseq_with_nan.csv')
    c.sort(by='padj', ascending=True, na_position='first', inplace=True)
    assert truth_first.equals(c.df)
    c.sort(by='padj', ascending=True, na_position='last', inplace=True)
    assert truth_last.equals(c.df)


def test_sort_descending(basic_countfilter):
    basic_countfilter.sort(by='cond3', ascending=False, inplace=True)
    assert basic_countfilter.df['cond3'].to_pandas().is_monotonic_decreasing


def test_filter_missing_values():
    truth = io.load_table('tests/test_files/test_deseq_with_nan_all_removed.csv')
    f = Filter('tests/test_files/test_deseq_with_nan.csv')
    f.filter_missing_values()
    assert np.all(f.df.sort(pl.first()) == truth.sort(pl.first()))


def test_filter_missing_values_foldchangefilter():
    truth = io.load_table('tests/test_files/fc_1_nan_removed.csv', 0, squeeze=True)
    f = FoldChangeFilter('tests/test_files/fc_1_nan.csv', 'num', 'denom')
    res_all = f.filter_missing_values(inplace=False)
    assert truth.equals(res_all.df)
    res_foldchange = f.filter_missing_values(inplace=False)
    assert truth.equals(res_foldchange.df)


def test_filter_missing_values_one_columns():
    truth = io.load_table('tests/test_files/test_deseq_with_nan_basemean_removed.csv')
    f = Filter('tests/test_files/test_deseq_with_nan.csv')
    f.filter_missing_values('baseMean')
    assert truth.equals(f.df)


def test_filter_missing_values_multiple_columns():
    truth = io.load_table('tests/test_files/test_deseq_with_nan_basemean_pvalue_removed.csv')
    f = Filter('tests/test_files/test_deseq_with_nan.csv')
    f.filter_missing_values(['baseMean', 'pvalue'])
    assert truth.equals(f.df)


def test_filter_missing_values_invalid_type():
    f = Filter('tests/test_files/test_deseq_with_nan.csv')
    with pytest.raises(TypeError):
        f.filter_missing_values(columns={'baseMean': True, 'log2FolgChange': False})


def test_filter_missing_values_nonexistent_column():
    f = Filter('tests/test_files/test_deseq_with_nan.csv')
    with pytest.raises(InvalidValueError):
        f.filter_missing_values('pval')
    with pytest.raises(InvalidValueError):
        f.filter_missing_values(['padj', 'pval'])


def test_pipeline_api():
    pl = Pipeline()
    Pipeline('countfilter')
    Pipeline(DESeqFilter)
    pl = Pipeline(filter_type='FoldChangeFilter')
    assert pl.__len__() == 0


def test_pipeline_repr():
    pl = Pipeline('countfilter')
    assert repr(pl) == "Pipeline('CountFilter')"
    pl.add_function('sort')
    pl.add_function(CountFilter.filter_biotype_from_ref_table, 'protein_coding', opposite=True)
    assert repr(pl) == "Pipeline('CountFilter'): CountFilter.sort()-->" \
                       "CountFilter.filter_biotype_from_ref_table('protein_coding', opposite=True)"


def test_pipeline_str():
    pl = Pipeline('countfilter')
    assert str(pl) == "Pipeline for Count matrix"
    pl.add_function('sort')
    pl.add_function(CountFilter.filter_biotype_from_ref_table, 'protein_coding', opposite=True)
    assert str(pl) == "Pipeline for Count matrix:\n" \
                      "\tSort table rows: ()\n" \
                      "\tFilter by feature biotype (based on a reference table): ('protein_coding', opposite=True)"


def test_pipeline_add_function():
    pl = Pipeline()
    pl.add_function(DESeqFilter.filter_biotype_from_ref_table, biotype='protein_coding')
    assert len(pl.functions) == 1 and len(pl.params) == 1 and len(pl) == 1
    assert pl.functions[0] == DESeqFilter.filter_biotype_from_ref_table
    assert pl.params[0] == ((), {'biotype': 'protein_coding'})

    pl = Pipeline()
    pl.add_function('filter_biotype_from_ref_table', 'piRNA')
    assert len(pl.functions) == 1 and len(pl.params) == 1 and len(pl) == 1
    assert pl.functions[0] == Filter.filter_biotype_from_ref_table
    assert pl.params[0] == (('piRNA',), {})

    pl_deseq = Pipeline('DEseqFilter')
    pl_deseq.add_function(Filter.number_filters, 'log2FoldChange', operator='>', value=5)
    assert len(pl_deseq.functions) == 1 and len(pl_deseq.params) == 1 and len(pl_deseq) == 1
    assert pl_deseq.functions[0] == DESeqFilter.number_filters
    assert pl_deseq.params[0] == (('log2FoldChange',), {'operator': '>', 'value': 5})


def test_pipeline_add_multiple_functions():
    pl_deseq = Pipeline('DEseqFilter')
    pl_deseq.add_function(Filter.number_filters, 'log2FoldChange', operator='>', value=5)
    pl_deseq.add_function(DESeqFilter.filter_significant)
    pl_deseq.add_function('sort', by='log2FoldChange')

    assert len(pl_deseq.functions) == 3 and len(pl_deseq.params) == 3 and len(pl_deseq) == 3
    assert pl_deseq.functions == [DESeqFilter.number_filters, DESeqFilter.filter_significant, DESeqFilter.sort]
    assert pl_deseq.params == [(('log2FoldChange',), {'operator': '>', 'value': 5}), ((), {}),
                               ((), {'by': 'log2FoldChange'})]


def test_pipeline_remove_last_function():
    pl = Pipeline()
    pl.add_function(DESeqFilter.filter_biotype_from_ref_table, biotype='protein_coding',
                    ref=__biotype_ref__)
    pl.remove_last_function()
    assert len(pl.functions) == 0 and len(pl.params) == 0 and len(pl) == 0


def test_pipeline_remove_last_from_empty_pipeline():
    pl = Pipeline()
    with pytest.raises(InvalidValueError):
        pl.remove_last_function()


def test_pipeline_apply_empty_pipeline(basic_deseqfilter):
    pl = Pipeline()
    with pytest.raises(InvalidValueError):
        pl.apply_to(basic_deseqfilter)


def test_pipeline_apply_to(basic_countfilter, basic_deseqfilter):
    pl = Pipeline('deseqfilter')
    pl.add_function('filter_significant', 10 ** -70, opposite=True)
    deseq_truth = basic_deseqfilter.__copy__()
    deseq_truth.filter_significant(10 ** -70, opposite=True)
    deseq_pipelined = pl.apply_to(basic_deseqfilter, inplace=False)
    pl.apply_to(basic_deseqfilter)
    basic_deseqfilter.sort('log2FoldChange')
    deseq_truth.sort('log2FoldChange')
    deseq_pipelined.sort('log2FoldChange')
    assert np.all(basic_deseqfilter.df == deseq_truth.df)
    assert np.all(deseq_pipelined.df == deseq_truth.df)

    pl2 = Pipeline('countfilter')
    pl2.add_function(Filter.filter_biotype_from_ref_table, biotype='protein_coding',
                     ref=__biotype_ref__)
    cnt_truth = basic_countfilter.__copy__()
    cnt_truth.filter_biotype_from_ref_table('protein_coding', ref=__biotype_ref__)
    cnt_pipelined = pl2.apply_to(basic_countfilter, inplace=False)
    pl2.apply_to(basic_countfilter, inplace=True)
    basic_countfilter.sort(basic_countfilter.columns[0])
    cnt_truth.sort(basic_countfilter.columns[0])
    cnt_pipelined.sort(basic_countfilter.columns[0])
    assert np.all(basic_countfilter.df == cnt_truth.df)
    assert np.all(cnt_pipelined.df == cnt_truth.df)


def test_pipeline_apply_to_with_multiple_functions():
    d = DESeqFilter('tests/test_files/test_deseq_with_nan.csv')
    d_copy = d.__copy__()
    p = Pipeline('deseqfilter')
    p.add_function('filter_missing_values')
    p.add_function(Filter.filter_biotype_from_ref_table, biotype='protein_coding', ref=__biotype_ref__)
    p.add_function('number_filters', 'log2FoldChange', 'gt', 0.75)
    p.add_function('sort', 'baseMean', ascending=False)
    p.add_function('biotypes_from_ref_table')
    p.remove_last_function()

    d_pipelined = p.apply_to(d_copy, inplace=False)
    p.apply_to(d_copy)
    d.filter_missing_values()
    d.filter_biotype_from_ref_table('protein_coding', __biotype_ref__)
    d.number_filters('log2FoldChange', 'gt', 0.75)
    d.sort('baseMean', ascending=False)
    assert d.df.equals(d_pipelined.df)
    assert d.df.equals(d_copy.df)


def test_pipeline_apply_to_invalid_object(basic_countfilter):
    pl = Pipeline('deseqfilter')
    pl.add_function(DESeqFilter.filter_significant, alpha=10 ** -70)
    with pytest.raises(InvalidTypeError):
        pl.apply_to(basic_countfilter)


def test_pipeline_init_invalid_filter_type():
    with pytest.raises(RNAlysisInputError):
        Pipeline(filter_type='otherFilter')

    class otherFilter:
        def __init__(self):
            self.value = 'value'
            self.othervalue = 'othervalue'

    with pytest.raises(RNAlysisInputError):
        Pipeline(filter_type=otherFilter)
    with pytest.raises(RNAlysisInputError):
        Pipeline(filter_type=max)
    with pytest.raises(RNAlysisInputError):
        Pipeline(filter_type=5)


def test_pipeline_add_function_out_of_module():
    pl = Pipeline()
    with pytest.raises(InvalidValueError):
        pl.add_function(len)

    with pytest.raises(InvalidValueError):
        pl.add_function(np.sort, arg='val')


def test_pipeline_add_function_invalid_type():
    pl = Pipeline()
    with pytest.raises(InvalidValueError):
        pl.add_function('string', arg='val')


def test_pipeline_add_function_mismatch_filter_type():
    pl_deseq = Pipeline('DESeqFilter')
    pl_deseq.add_function(CountFilter.filter_biotype_from_ref_table, biotype='protein_coding',
                          ref=__biotype_ref__)
    with pytest.raises(InvalidValueError):
        pl_deseq.add_function(CountFilter.filter_low_reads, threshold=5)


def _get_pipeline_with_plot(inplace: bool):
    p = Pipeline('DESeqFilter')
    p.add_function('filter_missing_values')
    p.add_function(Filter.filter_top_n, 'baseMean', n=15)
    p.add_function(DESeqFilter.volcano_plot, alpha=0.001)
    p.add_function('biotypes_from_ref_table', ref=__biotype_ref__)
    p.add_function('filter_top_n', 'log2FoldChange', 6)
    d = DESeqFilter('tests/test_files/test_deseq_with_nan.csv')
    res_truth = {}
    d_copy = d.__copy__()
    d.filter_missing_values()
    d.filter_top_n('baseMean', n=15)
    res_truth['volcano_plot_1'] = d.volcano_plot(alpha=0.001)
    res_truth['biotypes_from_ref_table_1'] = d.biotypes_from_ref_table(ref=__biotype_ref__)
    d.filter_top_n('log2FoldChange', 6)
    if inplace:
        res = p.apply_to(d_copy, inplace=True)
        assert d.df.equals(d_copy.df)
    else:
        d_pipelined, res = p.apply_to(d_copy, inplace=False)
        assert d.df.equals(d_pipelined.df)
    assert res.keys() == res_truth.keys()
    assert res['biotypes_from_ref_table_1'].equals(res_truth['biotypes_from_ref_table_1'])
    assert type(res['volcano_plot_1']) == type(res_truth['volcano_plot_1'])


def test_pipeline_apply_to_with_plot_inplace():
    _get_pipeline_with_plot(inplace=True)


def test_pipeline_apply_to_with_plot_not_inplace():
    _get_pipeline_with_plot(inplace=False)


def test_pipeline_apply_to_with_split_function():
    # TODO: prettify me!
    pl_d = Pipeline('DESeqFilter')
    pl_d.add_function(DESeqFilter.filter_missing_values)
    pl_d.add_function(DESeqFilter.split_fold_change_direction)
    pl_d.add_function(DESeqFilter.filter_top_n, by='padj', n=3)
    pl_d.add_function('sort', by='baseMean')
    d = DESeqFilter('tests/test_files/test_deseq_with_nan.csv')
    d_pipeline_res = pl_d.apply_to(d, inplace=False)
    d_res = d.filter_missing_values(inplace=False)
    d_res = d_res.split_fold_change_direction()
    for i in d_res:
        i.filter_top_n(by='padj', n=3)
        i.sort(by='baseMean')
    for i, j in zip(d_res, d_pipeline_res):
        assert i == j

    pl_c = Pipeline('CountFilter')
    pl_c.add_function(CountFilter.filter_top_n, by='cond2', n=2, opposite=True)
    pl_c.add_function(CountFilter.split_hdbscan, min_cluster_size=3, return_probabilities=True)
    pl_c.add_function(CountFilter.filter_top_n, by='cond1', n=5)
    c = CountFilter('tests/test_files/counted.csv')
    c_pipeline_res, c_pipeline_dict = pl_c.apply_to(c, inplace=False)
    c_res = c.filter_top_n(by='cond2', n=2, opposite=True, inplace=False)
    c_res, prob = c_res.split_hdbscan(min_cluster_size=3, return_probabilities=True)
    for i in c_res:
        i.filter_top_n(by='cond1', n=5)
    for i, j in zip(c_res, c_pipeline_res):
        assert i == j
    assert np.all(c_pipeline_dict['split_hdbscan_1'] == prob)

    pl_c.remove_last_function()
    pl_c.remove_last_function()
    pl_c.add_function('split_kmeans', n_clusters=[2, 3, 4], random_seed=42)
    pl_c.add_function(CountFilter.filter_top_n, by='cond1', n=5)
    c = CountFilter('tests/test_files/counted.csv')
    c_pipeline_res = pl_c.apply_to(c, inplace=False)
    c_res = c.filter_top_n(by='cond2', n=2, opposite=True, inplace=False)
    c_res = c_res.split_kmeans(n_clusters=[2, 3, 4], random_seed=42)
    for tpl in c_res:
        for this_obj in tpl:
            this_obj.filter_top_n(by='cond1', n=5)
    for i, j in zip(c_res, c_pipeline_res):
        assert i == j


def test_pipeline_apply_to_with_split_function_inplace_raise_error(basic_deseqfilter):
    pl = Pipeline('DESeqFilter')
    pl.add_function('filter_missing_values')
    pl.add_function(DESeqFilter.filter_significant, 0.05, opposite=True)
    pl.add_function(DESeqFilter.split_fold_change_direction)
    with pytest.raises(InvalidValueError):
        pl.apply_to(basic_deseqfilter, inplace=True)


def test_pipeline_apply_to_multiple_splits(basic_countfilter):
    pl_c = Pipeline('CountFilter')
    pl_c.add_function(CountFilter.filter_top_n, by='cond2', n=2, opposite=True)
    pl_c.add_function(CountFilter.split_hdbscan, min_cluster_size=3, return_probabilities=True)
    pl_c.add_function(CountFilter.split_kmedoids, n_clusters=2, random_seed=42)
    pl_c.add_function(CountFilter.split_by_reads, 15)

    c_pipeline_res, c_pipeline_dict = pl_c.apply_to(basic_countfilter, inplace=False)
    c_res = basic_countfilter.filter_top_n(by='cond2', n=2, opposite=True, inplace=False)
    c_res, prob = c_res.split_hdbscan(min_cluster_size=3, return_probabilities=True)
    c_res_cont = []
    for this_obj in c_res:
        c_res_cont.append(this_obj.split_kmedoids(n_clusters=2, random_seed=42))
    c_res_cont_fin = []
    for tpl in c_res_cont:
        c_res_cont_fin.append([])
        for this_obj in tpl:
            c_res_cont_fin[-1].append(this_obj.split_by_reads(15))
    assert len(c_pipeline_res) == len(c_res_cont_fin)
    for i, j in zip(c_res_cont_fin, c_pipeline_res):
        assert parsing.data_to_tuple(i) == j
    assert np.all(c_pipeline_dict['split_hdbscan_1'] == prob)


def test_pipeline_apply_to_filter_normalize_split_plot(big_countfilter):
    # Note: n=17500 (instead of a smaller value) is deliberately chosen so that filter_top_n
    # leaves only ~2500 of the ~20000 rows before the clustering/plotting steps. This is a
    # self-consistency check (Pipeline.apply_to vs. the same calls made manually) with no
    # comparison to a committed truth file, so shrinking the data flowing through the slow
    # clustering steps (split_hdbscan, clustergram, split_kmedoids) is safe as long as it still
    # produces multiple non-trivial clusters, which it does (verified to leave >=2 clusters of
    # size >= 7, so split_kmedoids(n_clusters=[2, 7]) remains valid on each).
    scaling_factor_path = 'tests/test_files/big_scaling_factors.csv'
    pl_c = Pipeline('CountFilter')
    pl_c.add_function(CountFilter.normalize_with_scaling_factors, scaling_factor_path)
    pl_c.add_function('biotypes_from_ref_table', ref=__biotype_ref__)
    pl_c.add_function(CountFilter.filter_top_n, by='cond3rep1', n=17500, opposite=True)
    pl_c.add_function(CountFilter.split_hdbscan, min_cluster_size=40, return_probabilities=True)
    pl_c.add_function(CountFilter.filter_low_reads, threshold=10)
    pl_c.add_function(CountFilter.clustergram)
    pl_c.add_function(CountFilter.split_kmedoids, n_clusters=[2, 7], random_seed=42, n_init=1, max_iter=100)
    pl_c.add_function(CountFilter.sort, by='cond2rep3')
    pl_c.add_function(CountFilter.biotypes_from_ref_table, 'long', __biotype_ref__)

    c = big_countfilter
    c_pipeline_res, c_pipeline_dict = pl_c.apply_to(c, inplace=False)
    c_dict = dict()
    c.normalize_with_scaling_factors(scaling_factor_path)
    c_dict['biotypes_from_ref_table_1'] = c.biotypes_from_ref_table(ref=__biotype_ref__)
    c_res = c.filter_top_n(by='cond3rep1', n=17500, opposite=True, inplace=False)
    c_res, prob = c_res.split_hdbscan(min_cluster_size=40, return_probabilities=True)
    c_dict['split_hdbscan_1'] = prob
    for obj in c_res:
        obj.filter_low_reads(threshold=10)
    clustergrams = []
    for i, obj in enumerate(c_res):
        clustergrams.append(obj.clustergram())
    c_dict['clustergram_1'] = tuple(clustergrams)
    c_res_cont = []
    for obj in c_res:
        k_out = obj.split_kmedoids(n_clusters=[2, 7], n_init=1, max_iter=100, random_seed=42)
        c_res_cont.append(k_out)
    for tpl_outer in c_res_cont:
        for tpl_inner in tpl_outer:
            for obj in tpl_inner:
                obj.sort(by='cond2rep3')
    biotypes = []
    for tpl_outer in c_res_cont:
        for tpl_inner in tpl_outer:
            for obj in tpl_inner:
                biotypes.append(obj.biotypes_from_ref_table('long', __biotype_ref__))
    c_dict['biotypes_from_ref_table_2'] = tuple(biotypes)

    assert len(c_res_cont) == len(c_pipeline_res)
    for i, j in zip(c_res_cont, c_pipeline_res):
        assert parsing.data_to_tuple(i) == j
    assert np.all(c_dict.keys() == c_pipeline_dict.keys())
    assert c_dict['biotypes_from_ref_table_1'].equals(c_pipeline_dict['biotypes_from_ref_table_1'])
    assert len(c_dict['clustergram_1']) == len(c_pipeline_dict['clustergram_1'])
    for i, j in zip(c_dict['clustergram_1'], c_pipeline_dict['clustergram_1']):
        assert type(i) == type(j)
    assert len(c_dict['biotypes_from_ref_table_2']) == len(c_pipeline_dict['biotypes_from_ref_table_2'])
    for i, j in zip(c_dict['biotypes_from_ref_table_2'], c_pipeline_dict['biotypes_from_ref_table_2']):
        assert i.equals(j)


def test_gap_statistic_api(clustering_countfilter):
    # The gap statistic evaluates every candidate k from 2 up to max_n_clusters_estimate (default
    # 'auto' -> 20 here), fitting KMeans on 10 reference datasets per k. The assertion below only
    # checks the return type, not a specific k or a truth file, so capping max_n_clusters_estimate
    # at 5 (still >1, so the multi-k search code path is still exercised) cuts the number of KMeans
    # fits substantially while testing the same behavior.
    res = clustering_countfilter.split_kmeans(n_clusters='gap', max_n_clusters_estimate=5)
    assert isinstance(res, tuple)


def test_silhouette_api(clustering_countfilter):
    res = clustering_countfilter.split_kmeans(n_clusters='silhouette')
    assert isinstance(res, tuple)


def test_split_kmeans_api(clustering_countfilter):
    c = clustering_countfilter
    res = c.split_kmeans(n_clusters=4, n_init=2, max_iter=50)
    assert isinstance(res, tuple)
    assert len(res) == 4
    _test_correct_clustering_split(c, res)
    res2 = c.split_kmeans(n_clusters=[2, 3, 7], n_init=1, max_iter=10, random_seed=42, plot_style='std_bar')
    assert isinstance(res2, tuple)
    assert np.all([isinstance(i, tuple) for i in res2])
    [_test_correct_clustering_split(c, i) for i in res2]


def test_split_hierarchical_api(clustering_countfilter):
    c = clustering_countfilter
    res = c.split_hierarchical(4)
    assert isinstance(res, tuple)
    assert len(res) == 4
    _test_correct_clustering_split(c, res)
    res2 = c.split_hierarchical([2, 3, 7], metric='Euclidean', linkage='AVERAGE', plot_style='std_bar')
    assert isinstance(res2, tuple)
    assert np.all([isinstance(i, tuple) for i in res2])
    [_test_correct_clustering_split(c, i) for i in res2]

    res3 = c.split_hierarchical(None, distance_threshold=10)
    assert isinstance(res3, tuple)
    _test_correct_clustering_split(c, res3)

    with pytest.raises(InvalidValueError):
        c.split_hierarchical(5, metric='badinput')
    with pytest.raises(InvalidValueError):
        c.split_hierarchical(4, linkage='badinput')


def test_split_hdbscan_api(clustering_countfilter):
    c = clustering_countfilter
    res = c.split_hdbscan(100)
    _test_correct_clustering_split(c, res, True)
    res2 = c.split_hdbscan(4, 5, 'manhattan', 0.2, 'leaf', plot_style='std_area', return_probabilities=True)
    assert isinstance(res2, list)
    assert isinstance(res2[0], tuple)
    assert isinstance(res2[1], np.ndarray)
    with pytest.raises(InvalidValueError):
        c.split_hdbscan(min_cluster_size=1)
    with pytest.raises(InvalidValueError):
        c.split_hdbscan(c.shape[0] + 1)


def test_split_clicom_api(clustering_countfilter):
    # split_clicom's cluster-similarity-matrix step scales roughly with the square of the number
    # of features, and this test only checks structural correctness of the split (no comparison to
    # a committed truth file), so we shrink the working data to the first 300 (of 1345) features to
    # keep the same code path (all 3 clustering methods, multi-value parameter combos, and the
    # power_transform=(False, True) branch) while running much faster. 300 rows is still comfortably
    # more than the largest n_clusters/min_cluster_size used below.
    c = CountFilter.from_dataframe(clustering_countfilter.df.head(300), clustering_countfilter.fname,
                                   clustering_countfilter.is_normalized, suppress_warnings=True)
    res = c.split_clicom({'method': 'hdbscan', 'min_cluster_size': [75, 40], 'metric': ['ys1', 'spearman']},
                         {'method': 'hierarchical', 'n_clusters': [7], 'metric': ['euclidean', 'jackknife'],
                          'linkage': ['average', 'ward']},
                         {'method': 'kmedoids', 'n_clusters': [7, 16], 'metric': 'spearman'},
                         power_transform=(False, True), evidence_threshold=0.5, min_cluster_size=5)
    _test_correct_clustering_split(c, res, True)


def test_split_kmedoids_api(clustering_countfilter):
    c = clustering_countfilter
    res = c.split_kmedoids(n_clusters=4, n_init=3)
    assert isinstance(res, tuple)
    assert len(res) == 4
    _test_correct_clustering_split(c, res)
    res2 = c.split_kmedoids(n_clusters=[2, 3, 7], n_init=1, max_iter=10, random_seed=42, plot_style='std_bar')
    assert isinstance(res2, tuple)
    assert np.all([isinstance(i, tuple) for i in res2])
    [_test_correct_clustering_split(c, i) for i in res2]


def _test_correct_clustering_split(counts, res, missing_indices: bool = False):
    assert isinstance(res, tuple)
    for obj in res:
        assert len(counts.intersection(obj)) == obj.shape[0]
        # all of the indices in the split object are in the original too
        for obj2 in res:
            if obj != obj2:
                assert len(obj.intersection(obj2)) == 0
                # the clusters don't overlap with each other at all
    if not missing_indices:
        assert len(res[0].union(*res[1:])) == counts.shape[0]
        # if all values are clustered, make sure that all clusters sum up to the original object


def test_fc_randomization():
    truth = io.load_table('tests/test_files/fc_randomization_truth.csv')
    fc1 = FoldChangeFilter("tests/test_files/fc_1.csv", 'a', 'b')
    fc2 = FoldChangeFilter("tests/test_files/fc_2.csv", "c", "d")
    res = fc1.randomization_test(fc2, random_seed=0)
    assert truth['significant'].equals(res['significant'])
    assert np.allclose(truth.drop(cs.first()), res.drop(cs.first()))


def test_filter_save_csv():
    expected_path = Path('tests/test_files/test_deseq_with_nan_removemissingvals.csv')
    d = DESeqFilter('tests/test_files/test_deseq_with_nan.csv')
    d.filter_missing_values()
    try:
        d.save_csv()
        assert expected_path.exists()
        d_loaded = DESeqFilter(expected_path)
        assert d_loaded.df.equals(d.df)
    finally:
        expected_path.unlink(missing_ok=True)


def test_filter_save_csv_after_operations():
    expected_path = Path('tests/test_files/d_significant.csv')
    expected_path_2 = Path('tests/test_files/d_significant_with_suffix.csv')
    d = DESeqFilter('tests/test_files/test_deseq_with_nan.csv')
    d.filter_missing_values()
    d_sig = d.filter_significant(inplace=False)
    try:
        d_sig.save_csv(alt_filename='d_significant')
        d_sig.save_csv(alt_filename='d_significant_with_suffix.csv')
        assert expected_path.exists()
        assert expected_path_2.exists()
        assert DESeqFilter(expected_path).df.equals(d_sig.df)
        assert DESeqFilter(expected_path_2).df.equals(d_sig.df)
    finally:
        expected_path.unlink(missing_ok=True)
        expected_path_2.unlink(missing_ok=True)


def test_assert_padj_col():
    d = DESeqFilter('tests/test_files/test_deseq.csv')
    d._assert_padj_col()
    d = DESeqFilter('tests/test_files/test_deseq.csv', padj_col='baseMean')
    d._assert_padj_col()
    d = DESeqFilter('tests/test_files/test_deseq.csv', padj_col='not there')
    with pytest.raises(KeyError):
        d._assert_padj_col()


def test_assert_log2fc_col():
    d = DESeqFilter('tests/test_files/test_deseq.csv')
    d._assert_log2fc_col()
    d = DESeqFilter('tests/test_files/test_deseq.csv', log2fc_col='baseMean')
    d._assert_log2fc_col()
    d = DESeqFilter('tests/test_files/test_deseq.csv', log2fc_col='not there')
    with pytest.raises(KeyError):
        d._assert_log2fc_col()


@pytest.mark.parametrize('input_set,return_type,expected_output',
                         [({'gene1', 'gene3', 'gene2'}, 'set', {'gene1', 'gene3', 'gene2'}),
                          ({'gene1', 'gene3', 'gene2'}, 'str', 'gene1\ngene2\ngene3')])
def test_return_type(input_set, return_type, expected_output):
    assert Filter._return_type(input_set, return_type) == expected_output


def test_set_ops_wrong_type(basic_filter):
    with pytest.raises(TypeError):
        basic_filter._set_ops([{1, 2, 3}, 'string'], 'set')


@pytest.mark.parametrize('sample_list,truth_path',
                         [([['cond1', 'cond2'], ('cond3', 'cond4')], 'tests/test_files/counted_averaged_1.csv'),
                          ([['cond1'], ['cond2', 'cond3', 'cond4']], 'tests/test_files/counted_averaged_2.csv'),
                          (['cond1', ['cond2', 'cond3', 'cond4']], 'tests/test_files/counted_averaged_2.csv'),
                          (['cond1', 'cond2', 'cond3', 'cond4'], 'tests/test_files/counted.csv')])
def test_avg_subsamples(sample_list, truth_path, basic_countfilter):
    truth = io.load_table(truth_path)
    res = basic_countfilter._avg_subsamples(sample_list, new_column_names='auto')

    assert np.all(res.columns == truth.columns)
    assert np.all(res.select(pl.first()).equals(truth.select(pl.first())))
    assert np.isclose(res.drop(cs.first()), truth.drop(cs.first()), atol=0, rtol=0.001).all()


@pytest.mark.parametrize('function', ['mean', 'median', 'geometric_mean'])
def test_avg_subsamples_functions(function, basic_countfilter):
    sample_list = [['cond1', 'cond2'], ['cond3', 'cond4']]
    res = basic_countfilter._avg_subsamples(sample_list, function=function, new_column_names='auto')

    src = basic_countfilter.df
    numeric = res.drop(cs.first()).to_numpy()
    for i, group in enumerate(sample_list):
        group_vals = src.select(group).to_numpy()
        if function == 'mean':
            expected = group_vals.mean(axis=1)
        elif function == 'median':
            expected = np.median(group_vals, axis=1)
        else:
            with np.errstate(divide='ignore'):  # rows containing a zero -> geometric mean 0
                expected = np.exp(np.log(group_vals).mean(axis=1))  # geometric mean, row-wise
        assert np.allclose(numeric[:, i], expected, rtol=1e-6)


@pytest.mark.parametrize('input_file,expected_triplicates',
                         [('tests/test_files/counted.csv', [['cond1', 'cond2', 'cond3'], ['cond4']]),
                          ('tests/test_files/counted_6cols.csv',
                           [['cond1', 'cond2', 'cond3'], ['cond4', 'cond5', 'cond6']])])
def test_triplicates(input_file, expected_triplicates):
    counted = CountFilter(input_file)
    assert counted.triplicates == expected_triplicates


def _log2_plus1(df: pl.DataFrame):
    return np.log2(df + 1)


def _log10_plus1(df: pl.DataFrame):
    return np.log10(df + 1)


def _box_cox_plus1(df: pl.DataFrame):
    return PowerTransformer(method='box-cox').fit_transform(df + 1)


def _multiply_by_3_reduce_2(df: pl.DataFrame):
    return (df * 3) - 2


@pytest.mark.parametrize('filter_obj,columns,function,kwargs,truth_path', [
    (Filter('tests/test_files/counted_6cols.csv'), 'all', 'log2', {}, 'counted_6cols_log2.csv'),
    (DESeqFilter('tests/test_files/test_deseq.csv'), 'baseMean', 'log2', {}, 'test_deseq_log2.csv'),
    (CountFilter('tests/test_files/counted.csv'), ['cond1', 'cond4'], 'log2', {}, 'counted_log2.csv'),
    (Filter('tests/test_files/counted_6cols.csv'), 'all', 'log10', {}, 'counted_6cols_log10.csv'),
    (DESeqFilter('tests/test_files/test_deseq.csv'), 'baseMean', 'log10', {}, 'test_deseq_log10.csv'),
    (CountFilter('tests/test_files/counted.csv'), ['cond1', 'cond4'], 'log10', {}, 'counted_log10.csv'),
    (Filter('tests/test_files/counted_6cols.csv'), 'all', 'box-cox', {}, 'counted_6cols_boxcox.csv'),
    (DESeqFilter('tests/test_files/test_deseq.csv'), 'baseMean', 'box-cox', {}, 'test_deseq_boxcox.csv'),
    (CountFilter('tests/test_files/counted.csv'), ['cond1', 'cond4'], 'box-cox', {}, 'counted_boxcox.csv'),
    (Filter('tests/test_files/counted_6cols.csv'), 'all',
     lambda x, mult, red: (x * mult) - red, {'mult': 3, 'red': 2}, 'counted_6cols_multiply3_reduce2.csv'),
    (DESeqFilter('tests/test_files/test_deseq.csv'), 'baseMean',
     lambda x, mult, red: (x * mult) - red, {'mult': 3, 'red': 2}, 'test_deseq_multiply3_reduce2.csv'),
    (CountFilter('tests/test_files/counted.csv'), ['cond1', 'cond4'],
     lambda x, mult, red: (x * mult) - red, {'mult': 3, 'red': 2}, 'counted_multiply3_reduce2.csv')])
def test_transform(filter_obj, columns, function, kwargs, truth_path):
    truth = io.load_table(Path('tests/test_files/').joinpath(truth_path))
    filter_obj.transform(function, columns, **kwargs)
    print(filter_obj.df)
    print(truth)
    assert np.allclose(filter_obj.df.drop(cs.first()), truth.drop(cs.first()))


def test_import_pipeline():
    truth = Pipeline('countfilter')
    truth.add_function(CountFilter.biotypes_from_ref_table)
    truth.add_function('number_filters', 5, by='col_name')

    imported = Pipeline.import_pipeline('tests/test_files/test_pipeline.yaml')

    assert truth == imported


def test_export_pipeline():
    with open('tests/test_files/test_pipeline.yaml') as f:
        truth = yaml.safe_load(f)

    p = Pipeline('countfilter')
    p.add_function(CountFilter.biotypes_from_ref_table)
    p.add_function('number_filters', 5, by='col_name')
    outname = 'tests/test_files/temp_exported_pipeline.yaml'
    p.export_pipeline(outname)

    try:
        with open(outname) as f:
            exported = yaml.safe_load(f)
        assert set(exported.keys()) == set(truth.keys())
        for key in exported.keys():
            if key == 'metadata':
                assert set(exported[key].keys()) == set(truth[key].keys())
                for internal_key in exported[key]:
                    if internal_key == 'export_time':
                        continue
                    elif internal_key == 'rnalysis_version':
                        assert exported[key][internal_key] == __version__
            else:
                assert exported[key] == truth[key]
    finally:
        os.remove(outname)


# --- Pipeline YAML serialization robustness -----------------------------------------------------
# Pipeline YAMLs are a serialized artifact users keep for years (hard rule 4), so exporting must
# never crash on a value a user can plausibly pass, importing must never blame the user's typo on
# Python internals, and a Pipeline that cannot be loaded must say which RNAlysis version wrote it.


def test_export_pipeline_sanitizes_numpy_scalar_params():
    p = Pipeline('countfilter')
    p.add_function('filter_low_reads', threshold=np.float64(5.0))

    exported = yaml.safe_load(p.export_pipeline(None))

    threshold = exported['params'][0][1]['threshold']
    assert threshold == 5.0
    assert type(threshold) is float
    assert Pipeline.import_pipeline(p.export_pipeline(None)) == p


def test_export_pipeline_sanitizes_arrays_paths_and_nested_containers():
    p = Pipeline('countfilter')
    p.add_function('number_filters', np.int64(5), by=Path('some/dir/table.csv'),
                   nested={'values': np.array([1.5, 2.5]), 'flags': [np.bool_(True), (np.int32(3),)]})

    exported = yaml.safe_load(p.export_pipeline(None))

    args, kwargs = exported['params'][0]
    assert args == [5] and type(args[0]) is int
    assert kwargs['by'] == 'some/dir/table.csv'
    assert kwargs['nested'] == {'values': [1.5, 2.5], 'flags': [True, [3]]}


def test_export_pipeline_unrepresentable_param_raises_typed_error():
    class Unrepresentable:
        pass

    p = Pipeline('countfilter')
    p.add_function('filter_low_reads', threshold=Unrepresentable())

    with pytest.raises(RNAlysisInputError) as err:
        p.export_pipeline(None)
    assert 'filter_low_reads' in str(err.value)
    assert 'threshold' in str(err.value)


def test_export_pipeline_unrepresentable_positional_param_raises_typed_error():
    class Unrepresentable:
        pass

    p = Pipeline('countfilter')
    p.add_function('number_filters', Unrepresentable())

    with pytest.raises(RNAlysisInputError) as err:
        p.export_pipeline(None)
    assert 'number_filters' in str(err.value)


def test_import_pipeline_nonexistent_path_raises_file_not_found():
    missing = 'tests/test_files/no_such_pipeline_file.yaml'
    with pytest.raises(FileNotFoundError) as err:
        Pipeline.import_pipeline(missing)
    assert 'no_such_pipeline_file.yaml' in str(err.value)


def test_import_pipeline_nonexistent_path_object_raises_file_not_found():
    missing = Path('tests/test_files/no_such_pipeline_file.yaml')
    with pytest.raises(FileNotFoundError) as err:
        Pipeline.import_pipeline(missing)
    assert 'no_such_pipeline_file.yaml' in str(err.value)


def test_import_pipeline_from_yaml_string():
    truth = Pipeline('countfilter')
    truth.add_function(CountFilter.biotypes_from_ref_table)
    truth.add_function('number_filters', 5, by='col_name')

    with open('tests/test_files/test_pipeline.yaml') as f:
        content = f.read()

    assert Pipeline.import_pipeline(content) == truth


def test_import_pipeline_without_metadata_still_loads():
    truth = Pipeline('countfilter')
    truth.add_function(CountFilter.biotypes_from_ref_table)
    truth.add_function('number_filters', 5, by='col_name')

    assert Pipeline.import_pipeline('tests/test_files/test_pipeline_no_metadata.yaml') == truth


def test_import_pipeline_unknown_function_reports_exported_version():
    content = ("filter_type: countfilter\n"
               "metadata:\n"
               "   export_time: 2022/12/01, 16:47:40\n"
               "   rnalysis_version: 3.2.2\n"
               "functions:\n"
               "- a_function_that_never_existed\n"
               "params:\n"
               "- - []\n"
               "  - {}\n")

    with pytest.raises(RNAlysisInputError) as err:
        Pipeline.import_pipeline(content)
    message = str(err.value)
    assert '3.2.2' in message
    assert 'a_function_that_never_existed' in message
    assert __version__ in message


def test_import_pipeline_unknown_function_without_metadata_reports_unknown_version():
    content = ("filter_type: countfilter\n"
               "functions:\n"
               "- a_function_that_never_existed\n"
               "params:\n"
               "- - []\n"
               "  - {}\n")

    with pytest.raises(RNAlysisInputError) as err:
        Pipeline.import_pipeline(content)
    assert 'a_function_that_never_existed' in str(err.value)
    assert __version__ in str(err.value)


def test_import_pipeline_private_function_name_rejected():
    content = ("filter_type: countfilter\n"
               "functions:\n"
               "- _init_from_dict\n"
               "params:\n"
               "- - []\n"
               "  - {}\n")

    with pytest.raises(RNAlysisInputError):
        Pipeline.import_pipeline(content)


@pytest.mark.parametrize('content', [
    "metadata:\n   rnalysis_version: 9.9.9\nfilter_type: countfilter\nfunctions:\n- describe\n",
    "metadata:\n   rnalysis_version: 9.9.9\nfilter_type: countfilter\nparams:\n- - []\n  - {}\n",
    "metadata:\n   rnalysis_version: 9.9.9\nfunctions:\n- describe\nparams:\n- - []\n  - {}\n",
    "metadata:\n   rnalysis_version: 9.9.9\nfilter_type: not_a_filter_type\nfunctions: []\nparams: []\n",
    "metadata:\n   rnalysis_version: 9.9.9\nfilter_type: countfilter\nfunctions:\n- describe\n- describe\n"
    "params:\n- - []\n  - {}\n",
])
def test_import_pipeline_malformed_dict_reports_exported_version(content):
    with pytest.raises(RNAlysisInputError) as err:
        Pipeline.import_pipeline(content)
    assert '9.9.9' in str(err.value)


def test_import_pipeline_yaml_string_that_is_not_a_pipeline():
    with pytest.raises(RNAlysisInputError):
        Pipeline.import_pipeline('a: 1\nb: 2\n')


def test_import_pipeline_unparseable_yaml_string_raises_typed_error():
    with pytest.raises(RNAlysisInputError) as err:
        Pipeline.import_pipeline('functions: [\nparams: }\n')
    assert 'YAML' in str(err.value)


def test_import_pipeline_unparseable_yaml_file_raises_typed_error(tmp_path):
    broken = tmp_path.joinpath('broken_pipeline.yaml')
    broken.write_text('functions: [\nparams: }\n')

    with pytest.raises(RNAlysisInputError) as err:
        Pipeline.import_pipeline(broken)
    assert 'YAML' in str(err.value)


def test_import_pipeline_filter_type_is_case_insensitive():
    # Pipeline('CountFilter') is accepted, so a hand-written Pipeline YAML must be too
    truth = Pipeline('countfilter')
    truth.add_function('describe')

    imported = Pipeline.import_pipeline('filter_type: CountFilter\nfunctions:\n- describe\nparams:\n- - []\n  - {}\n')

    assert imported == truth


def test_export_pipeline_unhashable_set_item_raises_typed_error():
    p = Pipeline('countfilter')
    p.add_function('filter_low_reads', threshold=5, opposite={(1, 2)})

    with pytest.raises(RNAlysisInputError) as err:
        p.export_pipeline(None)
    assert 'filter_low_reads' in str(err.value)
    assert 'opposite' in str(err.value)


def test_pipeline_tuple_params_round_trip_as_lists():
    # documented, deliberate limitation: YAML has no tuple type, so nested tuples come back as
    # lists. Only the top-level positional-argument tuple is restored.
    p = Pipeline('countfilter')
    p.add_function('filter_low_reads', 5, opposite=(True,))

    imported = Pipeline.import_pipeline(p.export_pipeline(None))

    assert imported.params[0][0] == (5,)
    assert imported.params[0][1] == {'opposite': [True]}


@pytest.mark.parametrize('components,gene_fraction,truth_paths', [
    (1, 0.32, ['tests/test_files/counted_pc1_0.32_top.csv', 'tests/test_files/counted_pc1_0.32_bottom.csv'])
])
def test_split_by_principal_components(components, gene_fraction, truth_paths, basic_countfilter):
    truth = [CountFilter(pth) for pth in truth_paths]
    basic_countfilter.filter_low_reads(1)
    res = basic_countfilter.split_by_principal_components(components, gene_fraction)
    assert len(res) == len(truth)
    for i in range(len(truth)):
        assert res[i].df.sort(pl.first()).equals(truth[i].df.sort(pl.first()))


@pytest.mark.parametrize('ids,mode,truth_path', [
    ('kegg_id2', 'union', 'tests/test_files/counted_filter_by_kegg_truth_1.csv'),
    (['kegg_id2'], 'intersection', 'tests/test_files/counted_filter_by_kegg_truth_1.csv'),
    (['kegg_id3', 'kegg_id1'], 'union', 'tests/test_files/counted_filter_by_kegg_truth_2.csv'),
    (['kegg_id1', 'kegg_id3'], 'intersection', 'tests/test_files/counted_filter_by_kegg_truth_3.csv'),
    (['kegg_id1', 'kegg_id2', 'kegg_id3'], 'union', 'tests/test_files/counted_filter_by_kegg_truth_4.csv'),
    (['kegg_id1', 'kegg_id2', 'kegg_id3'], 'intersection', 'tests/test_files/counted_filter_by_kegg_truth_5.csv'),

])
def test_filter_by_kegg_annotations(monkeypatch, ids, mode, truth_path, basic_filter):
    truth = io.load_table(truth_path)

    def annotation_iter(self):
        annotations = [
            ['kegg_id1', None, {'WBGene00007063', 'WBGene00007064', 'WBGene00043988', 'other1', 'other2'}],
            ['kegg_id2', None, {'WBGene00007063', 'WBGene00044951', 'WBGene00007067', 'other3'}],
            ['kegg_id3', None, {'WBGene00007063', 'WBGene00007066', 'WBGene00043988', 'other1'}]]
        for annotation in annotations:
            if annotation[0] in ids:
                yield annotation

    monkeypatch.setattr(io.KEGGAnnotationIterator, 'get_kegg_organism_code', lambda *args, **kwargs: 'cel')
    monkeypatch.setattr(io.KEGGAnnotationIterator, 'get_pathway_annotations', annotation_iter)
    monkeypatch.setattr(io, 'get_taxon_and_id_type', lambda *args, **kwargs: ((6239, 'elegans'), 'KEGG'))

    res = basic_filter.filter_by_kegg_annotations(ids, mode, gene_id_type='WormBase', inplace=False)

    try:
        assert res.df.sort(pl.first()).equals(truth.sort(pl.first()))
    except Exception as e:
        raise e


@pytest.mark.parametrize('ids,mode,truth_path', [
    ('go_id2', 'union', 'tests/test_files/counted_filter_by_kegg_truth_1.csv'),
    (['go_id2'], 'intersection', 'tests/test_files/counted_filter_by_kegg_truth_1.csv'),
    (['go_id3', 'go_id1'], 'union', 'tests/test_files/counted_filter_by_kegg_truth_2.csv'),
    (['go_id1', 'go_id3'], 'intersection', 'tests/test_files/counted_filter_by_kegg_truth_3.csv'),
    (['go_id1', 'go_id2', 'go_id3'], 'union', 'tests/test_files/counted_filter_by_kegg_truth_4.csv'),
    (['go_id1', 'go_id2', 'go_id3'], 'intersection', 'tests/test_files/counted_filter_by_kegg_truth_5.csv'),

])
def test_filter_by_go_annotations(monkeypatch, ids, mode, truth_path, basic_filter):
    truth = io.load_table(truth_path)

    class MockGOTerm:
        def __init__(self, go_id: str):
            self.id = go_id
            self.namespace = 'biological_process'

    class MockDAGTree:
        def __init__(self):
            pass

        def __getitem__(self, item):
            if item in ['go_id1', 'go_id2', 'go_id3']:
                return MockGOTerm(item)

        def __contains__(self, item):
            if item in ['go_id1', 'go_id2', 'go_id3']:
                return True
            return False

        def upper_induced_graph_iter(self, item):
            return []

    def annotation_iter(self):
        annotations = [
            ['go_id1', {'WBGene00007063', 'WBGene00007064', 'WBGene00043988', 'other1', 'other2'}],
            ['go_id2', {'WBGene00007063', 'WBGene00044951', 'WBGene00007067', 'other3'}],
            ['go_id3', {'WBGene00007063', 'WBGene00007066', 'WBGene00043988', 'other1'}]]
        for annotation in annotations:
            for gene_id in annotation[1]:
                annotation_dict = dict(annotation_class=annotation[0], bioentity_internal_id=gene_id, source='WB')
                yield annotation_dict

    monkeypatch.setattr(io.GOlrAnnotationIterator, '_get_n_annotations', lambda *args, **kwargs: 13)
    monkeypatch.setattr(io.GOlrAnnotationIterator, '_annotation_generator_func', annotation_iter)
    monkeypatch.setattr(ontology, 'fetch_go_basic', lambda: MockDAGTree())

    res = basic_filter.filter_by_go_annotations(ids, mode, gene_id_type='WormBase', organism=6239, inplace=False)
    assert res.df.sort(pl.first()).equals(truth.sort(pl.first()))


@pytest.mark.parametrize("pth,cols,truth_pth", [
    ('tests/test_files/counted.csv', 'cond2', 'tests/test_files/counted_drop_columns_truth_1.csv'),
    ('tests/test_files/counted.csv', ['cond4', 'cond2'], 'tests/test_files/counted_drop_columns_truth_2.csv'),
])
def test_drop_columns(pth, cols, truth_pth):
    obj = Filter(pth)
    truth = Filter(truth_pth)
    obj.drop_columns(cols)
    assert obj.df.sort(pl.first()).equals(truth.df.sort(pl.first()))


@pytest.mark.parametrize('comparisons,expected_paths,script_path', [
    ([('replicate', 'rep2', 'rep3')], ['tests/test_files/deseq2_tests/case1/DESeq2_replicate_rep2_vs_rep3_truth.csv'],
     'tests/test_files/deseq2_tests/case1/expected_deseq_script_1.R'),
    ([('condition', 'cond2', 'cond1'), ('condition', 'cond3', 'cond2')],
     ['tests/test_files/deseq2_tests/case2/DESeq2_condition_cond2_vs_cond1_truth.csv',
      'tests/test_files/deseq2_tests/case2/DESeq2_condition_cond3_vs_cond2_truth.csv'],
     'tests/test_files/deseq2_tests/case2/expected_deseq_script_2.R')
])
def test_differential_expression_deseq2(big_countfilter, monkeypatch, comparisons, expected_paths, script_path):
    outdir = 'tests/test_files/deseq2_tests/outdir'
    sample_table_path = 'tests/test_files/test_design_matrix.csv'
    truth = parsing.data_to_tuple([DESeqFilter(file) for file in expected_paths])
    c = big_countfilter

    def mock_run_analysis(self):
        assert self.r_installation_folder == 'auto'
        assert self.comparisons == comparisons
        assert CountFilter(self.data_path).df.equals(c.df.sort(pl.first()))
        assert io.load_table(self.design_mat_path).equals(io.load_table(sample_table_path))

        return Path(script_path).parent

    monkeypatch.setattr(differential_expression.DESeqRunner, 'run', mock_run_analysis)
    try:
        res = c.differential_expression_deseq2(sample_table_path, comparisons, output_folder=outdir)
        assert sorted(res, key=lambda filter_obj: filter_obj.fname.name) == sorted(truth, key=lambda
            filter_obj: filter_obj.fname.name)
        with open(Path(outdir).joinpath(Path(script_path).name)) as outfile, open(script_path) as truthfile:
            assert outfile.read() == truthfile.read()
    finally:
        unlink_tree(outdir)


@pytest.mark.parametrize('random_effect', [None, 'replicate'])
@pytest.mark.parametrize('comparisons,expected_paths,script_path', [
    ([('replicate', 'rep2', 'rep3')], ['tests/test_files/limma_tests/case1/LimmaVoom_replicate_rep2_vs_rep3_truth.csv'],
     'tests/test_files/limma_tests/case1/expected_limma_script_1.R'),
    ([('condition', 'cond2', 'cond1'), ('condition', 'cond3', 'cond2')],
     ['tests/test_files/limma_tests/case2/LimmaVoom_condition_cond2_vs_cond1_truth.csv',
      'tests/test_files/limma_tests/case2/LimmaVoom_condition_cond3_vs_cond2_truth.csv'],
     'tests/test_files/limma_tests/case2/expected_limma_script_2.R')
])
def test_differential_expression_limma(big_countfilter, monkeypatch, comparisons, expected_paths, script_path,
                                       random_effect):
    outdir = 'tests/test_files/limma_tests/outdir'
    sample_table_path = 'tests/test_files/test_design_matrix.csv'
    truth = parsing.data_to_tuple(
        [DESeqFilter(file, log2fc_col='logFC', padj_col='adj.P.Val') for file in expected_paths])
    c = big_countfilter

    def mock_run_analysis(self):
        assert self.r_installation_folder == 'auto'
        assert self.comparisons == comparisons
        assert CountFilter(self.data_path).df.equals(c.df)
        assert io.load_table(self.design_mat_path).equals(io.load_table(sample_table_path))
        assert self.random_effect == random_effect

        return Path(script_path).parent

    monkeypatch.setattr(differential_expression.LimmaVoomRunner, 'run', mock_run_analysis)
    try:
        res = c.differential_expression_limma_voom(sample_table_path, comparisons, output_folder=outdir,
                                                   random_effect=random_effect)
        assert sorted(res, key=lambda filter_obj: filter_obj.fname.name) == sorted(truth, key=lambda
            filter_obj: filter_obj.fname.name)
        with open(Path(outdir).joinpath(Path(script_path).name)) as outfile, open(script_path) as truthfile:
            assert outfile.read() == truthfile.read()
    finally:
        unlink_tree(outdir)


@pytest.mark.parametrize('keep,exp_path', [
    ('first', 'tests/test_files/counted_duplicates_first.csv'),
    ('last', 'tests/test_files/counted_duplicates_last.csv'),
    ('neither', 'tests/test_files/counted_duplicates_neither.csv')
])
def test_filter_duplicate_ids(keep, exp_path):
    f = Filter('tests/test_files/counted_duplicates.csv')
    truth = io.load_table(exp_path).sort(pl.first())
    res = f.filter_duplicate_ids(keep, inplace=False)
    assert res.df.sort(pl.first()).equals(truth)
    f.filter_duplicate_ids(keep)
    assert f.df.sort(pl.first()).equals(truth)


def test_filter_by_row_name():
    names = ['WBGene00044951', 'WBGene00014997', 'WBGene00007069']
    f = Filter('tests/test_files/counted_duplicates.csv')
    truth = io.load_table('tests/test_files/counted_drop_names.csv')
    res = f.filter_by_row_name(names, inplace=False)
    assert res.df.equals(truth)
    f.filter_by_row_name(names)
    assert f.df.equals(truth)


class TestOrthologDictTableGenerator:
    def test_create_one2many_table_ortholog(self):
        # Test when the mode is 'ortholog'
        mapping_dict = {
            'gene1': ['ortholog1', 'ortholog2'],
            'gene2': ['ortholog3'],
        }
        expected_table = pl.DataFrame([
            ('gene1', 'ortholog1'),
            ('gene1', 'ortholog2'),
            ('gene2', 'ortholog3'),
        ], schema=['gene', 'ortholog'])

        result_table = Filter._create_one2many_table(io.OrthologDict(mapping_dict))
        assert result_table.equals(expected_table)

    def test_create_one2many_table_paralog(self):
        # Test when the mode is 'paralog'
        mapping_dict = {
            'gene1': ['paralog1', 'paralog2'],
            'gene2': ['paralog3'],
        }
        expected_table = pl.DataFrame([
            ('gene1', 'paralog1'),
            ('gene1', 'paralog2'),
            ('gene2', 'paralog3'),
        ], schema=['gene', 'paralog'])

        result_table = Filter._create_one2many_table(io.OrthologDict(mapping_dict), mode='paralog')
        assert result_table.equals(expected_table)

    def test_create_one2many_table_empty(self):
        # Test when the mapping_dict is empty
        mapping_dict = {}
        expected_table = pl.DataFrame(schema=['gene', 'ortholog'])

        result_table = Filter._create_one2many_table(io.OrthologDict(mapping_dict))
        assert result_table.equals(expected_table)

    def test_create_one2many_table_no_mapping(self):
        # Test when no mapping_dict is provided (None)
        expected_table = pl.DataFrame(schema=['gene', 'ortholog'])

        result_table = Filter._create_one2many_table(io.OrthologDict())
        assert result_table.equals(expected_table)


@patch('rnalysis.utils.io.PantherOrthologMapper')
def test_find_paralogs_panther(mock_mapper):
    df_truth = pl.DataFrame({'gene': ['gene1', 'gene1', 'gene2'], 'paralog': ['paralog1', 'paralog2', 'paralog3']})
    mock_mapper_instance = Mock()
    mock_mapper.return_value = mock_mapper_instance
    mock_mapper_instance.get_paralogs.return_value = io.OrthologDict(
        {'gene1': ['paralog1', 'paralog2'], 'gene2': ['paralog3']})

    filter_obj = Filter('tests/test_files/test_map_orthologs.csv')
    result = filter_obj.find_paralogs_panther(organism='Homo sapiens', gene_id_type='UniProtKB')

    assert result.equals(df_truth)
    mock_mapper.assert_called_once_with(9606, 9606, 'UniProtKB')
    mock_mapper_instance.get_paralogs.assert_called_once_with(parsing.data_to_tuple(filter_obj.index_set))


@pytest.mark.parametrize('filter_percent_identity', [True, False])
@patch('rnalysis.utils.io.EnsemblOrthologMapper')
def test_find_paralogs_ensembl(mock_mapper, filter_percent_identity):
    df_truth = pl.DataFrame({'gene': ['gene1', 'gene1', 'gene2'], 'paralog': ['paralog1', 'paralog2', 'paralog3']})
    mock_mapper_instance = Mock()
    mock_mapper.return_value = mock_mapper_instance
    mock_mapper_instance.get_paralogs.return_value = io.OrthologDict(
        {'gene1': ['paralog1', 'paralog2'], 'gene2': ['paralog3']})

    filter_obj = Filter('tests/test_files/test_map_orthologs.csv')
    result = filter_obj.find_paralogs_ensembl(organism='Homo sapiens', gene_id_type='UniProtKB',
                                              filter_percent_identity=filter_percent_identity)

    assert result.equals(df_truth)
    mock_mapper.assert_called_once_with(9606, 9606, 'UniProtKB')
    mock_mapper_instance.get_paralogs.assert_called_once_with(parsing.data_to_tuple(filter_obj.index_set),
                                                              filter_percent_identity)


@pytest.mark.parametrize('remove_unmapped', [True, False])
@pytest.mark.parametrize('non_unique_mode,filter_least_diverged', [
    ('first', True),
    ('last', False),
])
@patch('rnalysis.utils.io.PantherOrthologMapper')
def test_map_orthologs_panther(mock_mapper, non_unique_mode, filter_least_diverged, remove_unmapped):
    if remove_unmapped:
        filter_obj_path = 'tests/test_files/test_map_orthologs_remove_unmapped_truth.csv'
    else:
        filter_obj_path = 'tests/test_files/test_map_orthologs_truth.csv'
    filter_obj_truth = Filter(filter_obj_path)
    one2many_truth = pl.DataFrame(
        {'gene': ['gene1', 'gene1', 'gene2'], 'ortholog': ['ortholog1', 'ortholog2', 'ortholog3']})
    mock_mapper_instance = Mock()
    mock_mapper.return_value = mock_mapper_instance
    mock_mapper_instance.get_orthologs.return_value = (
        io.OrthologDict({'gene1': 'ortholog2', 'gene2': 'ortholog3'})
        , io.OrthologDict({'gene1': ['ortholog1', 'ortholog2'], 'gene2': ['ortholog3']}))

    filter_obj = Filter('tests/test_files/test_map_orthologs.csv')
    one2one, one2many = filter_obj.map_orthologs_panther('Homo sapiens', 'caenorhabditis elegans',
                                                         gene_id_type='UniProtKB', inplace=False,
                                                         non_unique_mode=non_unique_mode,
                                                         filter_least_diverged=filter_least_diverged,
                                                         remove_unmapped_genes=remove_unmapped)

    assert one2many.equals(one2many_truth)
    assert one2one == filter_obj_truth
    mock_mapper.assert_called_once_with(9606, 6239, 'UniProtKB')
    mock_mapper_instance.get_orthologs.assert_called_once_with(
        parsing.data_to_tuple(filter_obj.index_set), non_unique_mode, filter_least_diverged)


@pytest.mark.parametrize('remove_unmapped', [True, False])
@pytest.mark.parametrize('non_unique_mode,filter_percent_identity', [
    ('first', True),
    ('last', False),
])
@patch('rnalysis.utils.io.EnsemblOrthologMapper')
def test_map_orthologs_ensembl(mock_mapper, non_unique_mode, filter_percent_identity, remove_unmapped):
    if remove_unmapped:
        filter_obj_path = 'tests/test_files/test_map_orthologs_remove_unmapped_truth.csv'
    else:
        filter_obj_path = 'tests/test_files/test_map_orthologs_truth.csv'
    filter_obj_truth = Filter(filter_obj_path)
    one2many_truth = pl.DataFrame(
        {'gene': ['gene1', 'gene1', 'gene2'], 'ortholog': ['ortholog1', 'ortholog2', 'ortholog3']})
    mock_mapper_instance = Mock()
    mock_mapper.return_value = mock_mapper_instance
    mock_mapper_instance.get_orthologs.return_value = (
        io.OrthologDict({'gene1': 'ortholog2', 'gene2': 'ortholog3'})
        , io.OrthologDict({'gene1': ['ortholog1', 'ortholog2'], 'gene2': ['ortholog3']}))

    filter_obj = Filter('tests/test_files/test_map_orthologs.csv')
    one2one, one2many = filter_obj.map_orthologs_ensembl('Homo sapiens', 'caenorhabditis elegans',
                                                         gene_id_type='UniProtKB', inplace=False,
                                                         non_unique_mode=non_unique_mode,
                                                         filter_percent_identity=filter_percent_identity,
                                                         remove_unmapped_genes=remove_unmapped)

    assert one2many.equals(one2many_truth)
    assert one2one == filter_obj_truth
    mock_mapper.assert_called_once_with(9606, 6239, 'UniProtKB')
    mock_mapper_instance.get_orthologs.assert_called_once_with(
        parsing.data_to_tuple(filter_obj.index_set), non_unique_mode, filter_percent_identity)


@pytest.mark.parametrize('remove_unmapped', [True, False])
@pytest.mark.parametrize('non_unique_mode,filter_consistency_score,threshold', [
    ('first', True, 0.5),
    ('last', False, 0.93),
])
@patch('rnalysis.utils.io.PhylomeDBOrthologMapper')
def test_map_orthologs_phylomedb(mock_mapper, non_unique_mode, filter_consistency_score, remove_unmapped, threshold):
    if remove_unmapped:
        filter_obj_path = 'tests/test_files/test_map_orthologs_remove_unmapped_truth.csv'
    else:
        filter_obj_path = 'tests/test_files/test_map_orthologs_truth.csv'
    filter_obj_truth = Filter(filter_obj_path)
    one2many_truth = pl.DataFrame(
        {'gene': ['gene1', 'gene1', 'gene2'], 'ortholog': ['ortholog1', 'ortholog2', 'ortholog3']})
    mock_mapper_instance = Mock()
    mock_mapper.return_value = mock_mapper_instance
    mock_mapper_instance.get_orthologs.return_value = (
        io.OrthologDict({'gene1': 'ortholog2', 'gene2': 'ortholog3'})
        , io.OrthologDict({'gene1': ['ortholog1', 'ortholog2'], 'gene2': ['ortholog3']}))

    filter_obj = Filter('tests/test_files/test_map_orthologs.csv')
    one2one, one2many = filter_obj.map_orthologs_phylomedb('Homo sapiens', 'caenorhabditis elegans',
                                                           gene_id_type='UniProtKB', inplace=False,
                                                           non_unique_mode=non_unique_mode,
                                                           filter_consistency_score=filter_consistency_score,
                                                           consistency_score_threshold=threshold,
                                                           remove_unmapped_genes=remove_unmapped)

    assert one2many.equals(one2many_truth)
    assert one2one == filter_obj_truth
    mock_mapper.assert_called_once_with(9606, 6239, 'UniProtKB')
    mock_mapper_instance.get_orthologs.assert_called_once_with(
        parsing.data_to_tuple(filter_obj.index_set), non_unique_mode, threshold, filter_consistency_score)


@pytest.mark.parametrize('remove_unmapped', [True, False])
@pytest.mark.parametrize('non_unique_mode', ['first', 'last'])
@patch('rnalysis.utils.io.OrthoInspectorOrthologMapper')
def test_map_orthologs_orthoinspector(mock_mapper, non_unique_mode, remove_unmapped):
    if remove_unmapped:
        filter_obj_path = 'tests/test_files/test_map_orthologs_remove_unmapped_truth.csv'
    else:
        filter_obj_path = 'tests/test_files/test_map_orthologs_truth.csv'
    filter_obj_truth = Filter(filter_obj_path)
    one2many_truth = pl.DataFrame(
        {'gene': ['gene1', 'gene1', 'gene2'], 'ortholog': ['ortholog1', 'ortholog2', 'ortholog3']})
    mock_mapper_instance = Mock()
    mock_mapper.return_value = mock_mapper_instance
    mock_mapper_instance.get_orthologs.return_value = (
        io.OrthologDict({'gene1': 'ortholog2', 'gene2': 'ortholog3'})
        , io.OrthologDict({'gene1': ['ortholog1', 'ortholog2'], 'gene2': ['ortholog3']}))

    filter_obj = Filter('tests/test_files/test_map_orthologs.csv')
    one2one, one2many = filter_obj.map_orthologs_orthoinspector('Homo sapiens', 'caenorhabditis elegans',
                                                                gene_id_type='UniProtKB', inplace=False,
                                                                non_unique_mode=non_unique_mode,
                                                                remove_unmapped_genes=remove_unmapped)

    assert one2many.equals(one2many_truth)
    assert one2one == filter_obj_truth
    mock_mapper.assert_called_once_with(9606, 6239, 'UniProtKB')
    mock_mapper_instance.get_orthologs.assert_called_once_with(
        parsing.data_to_tuple(filter_obj.index_set), non_unique_mode)


@pytest.mark.parametrize('component,ascending,power_transform', [
    (1, True, True),
    (2, False, False),
    (3, True, True)
])
def test_sort_by_principal_component(clustering_countfilter, component, ascending, power_transform):
    # sort_by_principal_component(power_transform=True) runs a separate Box-Cox optimization per gene
    # (the data is transposed before the power transform), so its cost scales ~linearly with the
    # number of genes. This test only checks structural/self-consistency properties (isinstance,
    # inequality/equality between calls, and raised AssertionErrors) rather than a committed truth
    # file, so shrinking the working data to the first 300 (of 1345) genes keeps the same code path
    # while cutting the power_transform=True cases from ~25-30s each to a few seconds.
    c = CountFilter.from_dataframe(clustering_countfilter.df.head(300), clustering_countfilter.fname,
                                   clustering_countfilter.is_normalized, suppress_warnings=True)
    c_sorted = c.sort_by_principal_component(component, ascending=ascending, power_transform=power_transform,
                                             inplace=False)
    assert isinstance(c_sorted, CountFilter)
    assert not c_sorted.df.equals(c.df)

    c.sort_by_principal_component(component, ascending=ascending, power_transform=power_transform)
    assert c.df.equals(c_sorted.df)

    with pytest.raises(InvalidValueError):
        c.sort_by_principal_component(0)
    with pytest.raises(InvalidValueError):
        c.sort_by_principal_component(c.shape[0] + 1)


@pytest.mark.parametrize('adjusted_pvals,bin_size,title', [
    (False, 0.05, 'auto'),
    (True, 0.1, 'Custom Title')
])
def test_pval_histogram(basic_deseqfilter, adjusted_pvals, bin_size, title, monkeypatch):
    monkeypatch.setattr(plt, 'show', lambda: None)
    fig = basic_deseqfilter.pval_histogram(adjusted_pvals=adjusted_pvals, bin_size=bin_size, title=title)
    assert isinstance(fig, plt.Figure)
    plt.close(fig)


@pytest.mark.parametrize('features,split_plots',
                         [('WBGene00007063', False), (['WBGene00007063', 'WBGene00007064'], True),
                          (['WBGene00007063', 'WBGene00007064'], False)])
@pytest.mark.parametrize('samples', [[['cond1', 'cond2', 'cond4'], ['cond3']], [['cond1', 'cond2']], 'all'])
@pytest.mark.parametrize('avg_function,spread_function,count_unit', [
    ('mean', 'sem', 'Normalized reads'),
    ('median', 'std', 'TPM'),
    ('geometric_mean', 'range', 'RPKM'),
    ('geometric_mean', 'gstd', 'RPKM'),
    ('geometric_mean', 'gsem', 'RPKM'),
    ('median', 'iqr', 'RPKM')])
def test_plot_expression(features, samples, avg_function, spread_function, count_unit, split_plots, monkeypatch,
                         basic_countfilter):
    monkeypatch.setattr(plt, 'show', lambda: None)
    figs = basic_countfilter.plot_expression(features, samples, avg_function, spread_function, count_unit=count_unit,
                                             split_plots=split_plots)
    assert isinstance(figs, list)
    assert all(isinstance(fig, plt.Figure) for fig in figs)
    for fig in figs:
        plt.close(fig)


@pytest.mark.parametrize('x_logscale', [True, False])
@pytest.mark.parametrize('y_logscale', [True, False])
@pytest.mark.parametrize('column,bins,x_label', [
    ('baseMean', 100, 'auto'),
    ('padj', 16, 'adj. p-value')])
def test_histogram(column, bins, x_label, x_logscale, y_logscale, basic_deseqfilter):
    fig = basic_deseqfilter.histogram(column, bins, x_label, x_logscale, y_logscale)
    assert isinstance(fig, plt.Figure)
    plt.close(fig)


# ---------------------------------------------------------------------------
# infer_table_type: auto-detection of the most likely table type on load.
# Returns one of the keys used by the GUI's FILTER_OBJ_TYPES combo.
# Detection must be CONSERVATIVE (fall back to 'Other table' when unsure) and
# must NEVER raise on malformed input.
# ---------------------------------------------------------------------------
COUNT_KEY = 'Count matrix'
DE_KEY = 'Differential expression'
FC_KEY = 'Fold change'
OTHER_KEY = 'Other table'


@pytest.mark.parametrize('fname,expected', [
    # --- real count matrices (raw counts / normalized expression) ---
    ('tests/test_files/counted.csv', COUNT_KEY),
    ('tests/test_files/counted.tsv', COUNT_KEY),
    ('tests/test_files/counted_6cols.csv', COUNT_KEY),
    ('tests/test_files/big_counted.csv', COUNT_KEY),
    ('tests/test_files/elegans_developmental_stages.tsv', COUNT_KEY),
    # count matrix whose *filename* misleadingly says "fold_change" but content is counts:
    ('tests/test_files/counted_fold_change.csv', COUNT_KEY),
    # --- real differential-expression tables (DESeq2 + limma-voom) ---
    ('tests/test_files/test_deseq.csv', DE_KEY),
    ('tests/test_files/sample_deseq.csv', DE_KEY),
    ('tests/test_files/test_deseq_with_nan.csv', DE_KEY),
    ('tests/test_files/deseq2_tests/case1/DESeq2_replicate_rep2_vs_rep3_truth.csv', DE_KEY),
    ('tests/test_files/limma_tests/case1/LimmaVoom_replicate_rep2_vs_rep3_truth.csv', DE_KEY),
    ('tests/test_files/limma_tests/case2/LimmaVoom_condition_cond2_vs_cond1_truth.csv', DE_KEY),
    # --- real fold-change tables (single fold-change column) ---
    ('tests/test_files/fc_1.csv', FC_KEY),
    ('tests/test_files/fc_2.csv', FC_KEY),
    # --- generic / ambiguous tables -> Other ---
    ('tests/test_files/biotype_ref_table_for_tests.csv', OTHER_KEY),  # single string column
    ('tests/test_files/attribute_reference_table.csv', OTHER_KEY),  # single all-null attribute column
    ('tests/test_files/attr_ref_table_for_tests.csv', OTHER_KEY),  # numeric attributes, has negatives
    ('tests/test_files/attr_ref_table_for_non_categorical.csv', OTHER_KEY),  # numeric attributes, has negatives
])
def test_infer_table_type_fixtures(fname, expected):
    assert infer_table_type(fname) == expected


@pytest.mark.parametrize('fname', [
    'tests/test_files/test_deseq.csv',
    'tests/test_files/sample_deseq.csv',
    'tests/test_files/deseq2_tests/case1/DESeq2_replicate_rep2_vs_rep3_truth.csv',
    'tests/test_files/limma_tests/case1/LimmaVoom_replicate_rep2_vs_rep3_truth.csv',
])
def test_infer_table_type_de_never_detected_as_count(fname):
    # A differential-expression table must NEVER be preselected as a count matrix.
    assert infer_table_type(fname) != COUNT_KEY


@pytest.mark.parametrize('fname', [
    'tests/test_files/counted.csv',
    'tests/test_files/counted_6cols.csv',
    'tests/test_files/big_counted.csv',
])
def test_infer_table_type_count_never_detected_as_de(fname):
    # A count matrix must NEVER be preselected as a differential-expression table.
    assert infer_table_type(fname) != DE_KEY


def test_infer_table_type_accepts_polars_frame():
    df = pl.DataFrame({'gene': ['a', 'b', 'c'], 'cond1': [5, 10, 0], 'cond2': [3, 7, 22]})
    assert infer_table_type(df=df) == COUNT_KEY


def test_infer_table_type_frame_takes_precedence_over_fname():
    de = pl.DataFrame({'gene': ['a', 'b'], 'log2FoldChange': [-1.0, 2.0], 'padj': [0.01, 0.2]})
    # a count-matrix path is given, but the explicit frame is a DE table -> DE wins
    assert infer_table_type(fname='tests/test_files/counted.csv', df=de) == DE_KEY


@pytest.mark.parametrize('columns,expected', [
    # DESeq2-style
    (['log2FoldChange', 'padj'], DE_KEY),
    (['baseMean', 'log2FoldChange', 'lfcSE', 'stat', 'pvalue', 'padj'], DE_KEY),
    # limma-voom style
    (['logFC', 'AveExpr', 't', 'P.Value', 'adj.P.Val', 'B'], DE_KEY),
    # spaced/variant fold-change header still recognized alongside a p-value column
    (['log2 fold change', 'pvalue'], DE_KEY),
    # a lone log2fc column without a p-value column is NOT a DE table
    (['log2FoldChange'], FC_KEY),  # single fold-change-like column -> fold change
    # a lone p-value column without a fold-change column -> Other
    (['padj'], OTHER_KEY),
])
def test_infer_table_type_de_and_fc_by_columns(columns, expected):
    data = {'gene': ['g1', 'g2', 'g3']}
    for i, col in enumerate(columns):
        data[col] = [0.5 * (i + 1), -0.5 * (i + 1), 0.1 * (i + 1)]
    df = pl.DataFrame(data)
    assert infer_table_type(df=df) == expected


def test_infer_table_type_single_fold_change_column():
    df = pl.DataFrame({'gene': ['a', 'b', 'c'], 'Fold Change': [1.4, 2.6, 0.6]})
    assert infer_table_type(df=df) == FC_KEY


def test_infer_table_type_single_non_fold_numeric_column_is_other():
    # a single numeric column that isn't fold-change-like is ambiguous -> Other
    df = pl.DataFrame({'gene': ['a', 'b', 'c'], 'gene_length': [1200, 3400, 900]})
    assert infer_table_type(df=df) == OTHER_KEY


def test_infer_table_type_negative_numeric_matrix_is_other():
    # multi-column numeric with negatives (e.g. attribute table / transformed data) -> Other
    df = pl.DataFrame({'gene': ['a', 'b'], 'attr1': [-5, 3], 'attr2': [7, -2]})
    assert infer_table_type(df=df) == OTHER_KEY


def test_infer_table_type_binary_matrix_is_other():
    # a purely presence/absence {0,1} numeric matrix is not count-like -> Other
    df = pl.DataFrame({'gene': ['a', 'b', 'c'], 'attr1': [0, 1, 1], 'attr2': [1, 0, 1]})
    assert infer_table_type(df=df) == OTHER_KEY


def test_infer_table_type_non_numeric_junk_is_other():
    df = pl.DataFrame({'gene': ['a', 'b'], 'x': ['foo', 'bar'], 'y': ['baz', 'qux']})
    assert infer_table_type(df=df) == OTHER_KEY


def test_infer_table_type_tiny_count_matrix():
    df = pl.DataFrame({'gene': ['a'], 'cond1': [500], 'cond2': [300]})
    assert infer_table_type(df=df) == COUNT_KEY


def test_infer_table_type_all_null_numeric_column_is_other():
    df = pl.DataFrame({'gene': ['a', 'b'], 'x': [None, None]}, schema={'gene': pl.Utf8, 'x': pl.Float64})
    assert infer_table_type(df=df) == OTHER_KEY


def test_infer_table_type_empty_frame_is_other():
    assert infer_table_type(df=pl.DataFrame()) == OTHER_KEY


def test_infer_table_type_index_only_frame_is_other():
    df = pl.DataFrame({'gene': ['a', 'b', 'c']})
    assert infer_table_type(df=df) == OTHER_KEY


def test_infer_table_type_never_raises_on_bad_input():
    # None inputs, nonexistent files, and garbage must all degrade to 'Other table'.
    assert infer_table_type() == OTHER_KEY
    assert infer_table_type(None, None) == OTHER_KEY
    assert infer_table_type('tests/test_files/this_file_does_not_exist_12345.csv') == OTHER_KEY


def test_infer_table_type_never_raises_on_malformed_file(tmp_path):
    junk = tmp_path / 'junk.csv'
    junk.write_bytes(b'\x00\x01\x02not,a,real\x00table\xff\xfe')
    assert infer_table_type(str(junk)) == OTHER_KEY


def test_infer_table_type_fname_uses_lightweight_sampled_read(monkeypatch):
    # The fname path must do a lightweight, row-limited read (header + small sample) rather than a full
    # load of the whole file on the UI thread.
    from rnalysis import filtering as filtering_mod
    calls = []
    real_load = filtering_mod.io.load_table

    def spy(filename, *args, **kwargs):
        calls.append(kwargs)
        return real_load(filename, *args, **kwargs)

    monkeypatch.setattr(filtering_mod.io, 'load_table', spy)
    assert infer_table_type('tests/test_files/counted.csv') == COUNT_KEY
    assert len(calls) == 1
    assert filtering_mod._DETECT_SAMPLE_ROWS is not None
    assert calls[0].get('nrows') == filtering_mod._DETECT_SAMPLE_ROWS


def test_infer_table_type_sampling_limitation(tmp_path):
    # Documented limitation: the Count-matrix value checks only see a sampled window of rows. A table that
    # is non-negative within the sample but turns negative later is (acceptably) pre-selected as a count
    # matrix -- still a user-overridable default, not a hard classification.
    from rnalysis import filtering as filtering_mod
    n = filtering_mod._DETECT_SAMPLE_ROWS
    genes = [f'g{i}' for i in range(n + 50)]
    cond1 = [5.0] * (n + 50)
    cond2 = [3.0] * (n + 50)
    cond2[n + 10] = -7.0  # a negative value that appears only AFTER the sampled window
    full = pl.DataFrame({'gene': genes, 'cond1': cond1, 'cond2': cond2})
    path = tmp_path / 'late_negative.csv'
    full.write_csv(path)
    # reading only the first n rows sees a non-negative matrix -> Count matrix
    assert infer_table_type(str(path)) == COUNT_KEY
    # the full in-memory frame sees the late negative -> Other (proving the above is a sampling effect)
    assert infer_table_type(df=full) == OTHER_KEY
