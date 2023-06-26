import json
import shutil

import matplotlib
import pytest
import yaml

from rnalysis import __version__
from rnalysis.filtering import *
from tests import __attr_ref__, __biotype_ref__

matplotlib.use('Agg')


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


def test_countfilter_api():
    h = CountFilter('tests/test_files/counted.csv')
    assert h.__str__() == "Count matrix named 'counted.csv'"
    assert h.__repr__() == "CountFilter('tests/test_files/counted.csv')"


def test_deseqfilter_api():
    d = DESeqFilter('tests/test_files/test_deseq.csv')
    assert d.__str__() == "Differential expression table named 'test_deseq.csv'"
    assert d.__repr__() == "DESeqFilter('tests/test_files/test_deseq.csv')"


def test_foldchangefilter_api():
    fc = FoldChangeFilter("tests/test_files/fc_1.csv", 'a', 'b')
    assert fc.__str__() == "FoldChangeFilter (numerator: 'a', denominator: 'b') of file fc_1.csv"
    assert fc.__repr__() == "FoldChangeFilter('tests/test_files/fc_1.csv', 'a', 'b')"


def test_filter_contains():
    objs = [Filter('tests/test_files/test_deseq.csv'), CountFilter('tests/test_files/counted.csv'),
            DESeqFilter('tests/test_files/counted.csv'), FoldChangeFilter('tests/test_files/fc_1.csv', 'num', 'denom')]
    neither = ['WBGene' + str(i) for i in range(10)]
    for obj in objs:
        for ind in obj.df.index:
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
    assert np.all(d_inplace_false.df == truth)
    assert np.all(d.df == d_copy.df)
    d._inplace(truth, opposite=False, inplace=True, suffix='other_suffix')
    assert np.all(d.df == truth)


def test_head():
    df = io.load_table('tests/test_files/test_deseq.csv', 0)
    d = DESeqFilter('tests/test_files/test_deseq.csv')
    assert np.all(df.head(7) == d.head(7))
    assert np.all(df.head(1) == d.head(1))

    df2 = io.load_table('tests/test_files/counted.csv', 0)
    f = Filter('tests/test_files/counted.csv')
    assert np.all(df2.head() == f.head())
    assert np.all(df2.head(1000) == f.head(1000))


def test_tail():
    df = io.load_table('tests/test_files/test_deseq.csv', 0)
    d = DESeqFilter('tests/test_files/test_deseq.csv')
    assert np.all(df.tail(7) == d.tail(7))
    assert np.all(df.tail(1) == d.tail(1))

    df2 = io.load_table('tests/test_files/counted.csv', 0)
    f = Filter('tests/test_files/counted.csv')
    assert np.all(df2.tail() == f.tail())
    assert np.all(df2.tail(1000) == f.tail(1000))


def test_describe():
    fc_df = io.load_table('tests/test_files/fc_1.csv', 0, squeeze=True)
    count_df = io.load_table('tests/test_files/counted.csv', 0)
    deseq_df = io.load_table('tests/test_files/test_deseq.csv', 0)

    fc = FoldChangeFilter('tests/test_files/fc_1.csv', 'a', 'b')
    count = CountFilter('tests/test_files/counted.csv')
    deseq = DESeqFilter('tests/test_files/test_deseq.csv')
    default_percentiles = (0.01, 0.25, 0.5, 0.75, 0.99)
    deciles = [i / 10 for i in range(1, 10)]
    for filter_obj, df in zip([count, deseq, fc], [count_df, deseq_df, fc_df]):
        assert np.all(filter_obj.describe() == df.describe(percentiles=default_percentiles))
        assert np.all(filter_obj.describe(deciles) == df.describe(percentiles=deciles))


def test_print_features_api():
    count = CountFilter('tests/test_files/counted.csv')
    count.print_features()


def test_from_string(monkeypatch):
    monkeypatch.setattr('builtins.input', lambda x: 'first-word\nsecond word\nthird. word\n')
    assert Filter._from_string("msg") == ["first-word", "second word", "third. word"]
    monkeypatch.setattr('builtins.input', lambda x: 'first-word,second word,third. word')
    assert Filter._from_string("msg", delimiter=',') == ["first-word", "second word", "third. word"]


@pytest.mark.parametrize("map_to,map_from,remove_unmapped_genes,,expected", [
    ('UniProtKB AC/ID', 'WormBase', False, 'tests/test_files/counted_translated_with_unmapped.csv'),
    ('UniProtKB', 'auto', True, 'tests/test_files/counted_translated_remove_unmapped.csv'),
])
def test_filter_translate_gene_ids(map_to, map_from, remove_unmapped_genes, expected, monkeypatch):
    def mock_map_gene_ids(ids, trans_from, trans_to='UniProtKB AC', verbose=True):
        if trans_from == 'WormBase':
            return io.GeneIDDict(
                {'WBGene00007063': 'A0A0K3AWR5', 'WBGene00007064': 'A0A2X0T1Z3', 'WBGene00007067': 'D3NQA2',
                 'WBGene00077503': 'H2L2B5', 'WBGene00007071': 'Q17405', 'WBGene00014997': 'Q7JNR0',
                 'WBGene00043988': 'A4F2Z7', 'WBGene00043989': 'G5EFZ2', 'WBGene00007075': 'G5EDW3',
                 'WBGene00007076': 'G5EFZ2'})
        return io.GeneIDDict({})

    monkeypatch.setattr(io, 'map_gene_ids', mock_map_gene_ids)
    truth = io.load_table(expected, index_col=0)
    f = Filter('tests/test_files/counted.csv')

    res = f.translate_gene_ids(map_to, map_from, remove_unmapped_genes, inplace=False)
    assert res.df.sort_index().equals(truth.sort_index())
    f.translate_gene_ids(map_to, map_from, remove_unmapped_genes, inplace=True)
    assert f.df.sort_index().equals(truth.sort_index())


def test_countfilter_normalize_to_rpm_htseqcount():
    truth = io.load_table(r"tests/test_files/test_norm_reads_rpm_htseqcount.csv", 0)
    h = CountFilter("tests/test_files/counted.csv")
    not_inplace, factors = h.normalize_to_rpm_htseqcount("tests/test_files/uncounted.csv", inplace=False,
                                                         return_scaling_factors=True)
    assert np.isclose(truth, not_inplace.df).all()

    norm_df = h._norm_scaling_factors(factors)
    assert np.isclose(truth, norm_df).all()

    h.normalize_to_rpm_htseqcount("tests/test_files/uncounted.csv")
    assert np.isclose(truth, h.df).all()


def test_countfilter_normalize_to_rpm():
    truth = io.load_table(r"tests/test_files/test_norm_to_rpm.csv", 0)
    h = CountFilter("tests/test_files/counted.csv")
    not_inplace, factors = h.normalize_to_rpm(inplace=False, return_scaling_factors=True)
    assert np.isclose(truth, not_inplace.df).all()

    norm_df = h._norm_scaling_factors(factors)
    assert np.isclose(truth, norm_df).all()

    h.normalize_to_rpm()
    assert np.isclose(truth, h.df).all()


@pytest.mark.parametrize('gtf_path,feature_type,method',
                         [
                             ('tests/test_files/test_gtf_wormbase.gtf', 'gene', 'merged_exons'),
                             ('tests/test_files/test_gff3_wormbase.gff3', 'gene', 'geometric_mean'),
                             ('tests/test_files/test_gtf_wormbase.gtf', 'transcript', 'mean'),
                             ('tests/test_files/test_gff3_wormbase.gff3', 'transcript', 'max')
                         ])
def test_countfilter_normalize_to_rpkm(monkeypatch, gtf_path, feature_type, method):
    def mock_get_feature_lengths(this_gtf_path, this_feature_type, this_method):
        assert this_gtf_path == gtf_path
        assert this_feature_type == feature_type
        assert this_method == method
        with open('tests/test_files/feature_lengths.json') as f:
            return json.load(f)

    monkeypatch.setattr(genome_annotation, 'get_genomic_feature_lengths', mock_get_feature_lengths)
    truth = io.load_table(r"tests/test_files/test_norm_to_rpkm.csv", 0)
    h = CountFilter("tests/test_files/counted.csv")
    not_inplace, factors = h.normalize_to_rpkm(gtf_path, feature_type, method, inplace=False,
                                               return_scaling_factors=True)
    assert np.isclose(truth, not_inplace.df).all()

    norm_df = h._norm_scaling_factors(factors)
    assert np.isclose(truth, norm_df).all()

    h.normalize_to_rpkm(gtf_path, feature_type, method, )
    assert np.isclose(truth, h.df).all()


@pytest.mark.parametrize('gtf_path,feature_type,method',
                         [
                             ('tests/test_files/test_gtf_wormbase.gtf', 'gene', 'merged_exons'),
                             ('tests/test_files/test_gff3_wormbase.gff3', 'gene', 'geometric_mean'),
                             ('tests/test_files/test_gtf_wormbase.gtf', 'transcript', 'mean'),
                             ('tests/test_files/test_gff3_wormbase.gff3', 'transcript', 'max')
                         ])
def test_countfilter_normalize_to_tpm(monkeypatch, gtf_path, feature_type, method):
    def mock_get_feature_lengths(this_gtf_path, this_feature_type, this_method):
        assert this_gtf_path == gtf_path
        assert this_feature_type == feature_type
        assert this_method == method
        with open('tests/test_files/feature_lengths.json') as f:
            return json.load(f)

    monkeypatch.setattr(genome_annotation, 'get_genomic_feature_lengths', mock_get_feature_lengths)
    truth = io.load_table(r"tests/test_files/test_norm_to_tpm.csv", 0)
    h = CountFilter("tests/test_files/counted.csv")
    not_inplace, factors = h.normalize_to_tpm(gtf_path, feature_type, method, inplace=False,
                                              return_scaling_factors=True)
    assert np.isclose(truth, not_inplace.df).all()

    norm_df = h._norm_scaling_factors(factors)
    assert np.isclose(truth, norm_df).all()

    h.normalize_to_tpm(gtf_path, feature_type, method)
    assert np.isclose(truth, h.df).all()


def test_countfilter_normalize_rle():
    truth = io.load_table(r"tests/test_files/test_norm_rle.csv", 0)
    h = CountFilter("tests/test_files/counted.csv")
    not_inplace, factors = h.normalize_rle(inplace=False, return_scaling_factors=True)
    assert np.isclose(truth, not_inplace.df).all()

    norm_df = h._norm_scaling_factors(factors)
    assert np.isclose(truth, norm_df).all()
    h.normalize_rle()
    assert np.isclose(truth, h.df).all()


def test_countfilter_normalize_tmm():
    truth = io.load_table(r"tests/test_files/test_norm_tmm.csv", 0)
    h = CountFilter("tests/test_files/counted.csv")
    not_inplace, factors = h.normalize_tmm(ref_column='cond1', inplace=False, return_scaling_factors=True)
    assert np.isclose(truth, not_inplace.df).all()

    norm_df = h._norm_scaling_factors(factors)
    assert np.isclose(truth, norm_df).all()

    h.normalize_tmm(ref_column='cond1')
    assert np.isclose(truth, h.df).all()


def test_countfilter_normalize_median_of_ratios():
    truth = io.load_table(r"tests/test_files/test_norm_mrn.csv", 0)
    h = CountFilter("tests/test_files/counted.csv")
    not_inplace, factors = h.normalize_median_of_ratios([['cond1', 'cond2'], ['cond3', 'cond4']], inplace=False,
                                                        return_scaling_factors=True)
    assert np.isclose(truth, not_inplace.df).all()

    norm_df = h._norm_scaling_factors(factors)
    assert np.isclose(truth, norm_df).all()

    h.normalize_median_of_ratios([['cond1', 'cond2'], ['cond3', 'cond4']])
    assert np.isclose(truth, h.df).all()


@pytest.mark.parametrize('quantile,truth_path', [
    (0.75, "tests/test_files/test_norm_quantile_75.csv"),
    (0.32, "tests/test_files/test_norm_quantile_32.csv"),
])
def test_countfilter_normalize_to_quantile(quantile, truth_path):
    truth = io.load_table(truth_path, 0)
    h = CountFilter("tests/test_files/counted.csv")
    not_inplace, factors = h.normalize_to_quantile(quantile, inplace=False, return_scaling_factors=True)
    assert np.isclose(truth, not_inplace.df).all()

    norm_df = h._norm_scaling_factors(factors)
    assert np.isclose(truth, norm_df).all()

    h.normalize_to_quantile(quantile)
    assert np.isclose(truth, h.df).all()


def test_countfilter_norm_reads_with_scaling_factors():
    truth = io.load_table(r"tests/test_files/test_norm_scaling_factors.csv", 0)
    h = CountFilter("tests/test_files/counted.csv")
    factors = io.load_table("tests/test_files/scaling_factors.csv")
    h_norm = h.normalize_with_scaling_factors("tests/test_files/scaling_factors.csv", inplace=False)
    h.normalize_with_scaling_factors(factors)
    assert np.isclose(truth, h.df).all()
    assert h_norm.df.equals(h.df)


def test_filter_low_reads():
    truth = io.load_table("tests/test_files/counted_low_rpm_truth.csv", 0)
    h = CountFilter("tests/test_files/counted_low_rpm.csv")
    h.filter_low_reads(threshold=5)
    assert np.isclose(truth, h.df).all()


def test_filter_low_reads_reverse():
    h = CountFilter("tests/test_files/counted.csv")
    low_truth = io.load_table(r"tests/test_files/counted_below60_rpm.csv", 0)
    h.filter_low_reads(threshold=60, opposite=True)
    h.df.sort_index(inplace=True)
    low_truth.sort_index(inplace=True)
    print(h.shape)
    print(low_truth.shape)
    print(h.df)
    print(low_truth)

    assert np.all(h.df == low_truth)


@pytest.mark.parametrize('interactive', [True, False])
@pytest.mark.parametrize('show_cursor', [True, False])
@pytest.mark.parametrize('title,alpha,log2fc_threshold', [
    ('auto', 0.05, None),
    ('title', 0.1, 0),
    ('title', 0.001, 1)])
def test_deseqfilter_volcano_plot_api(interactive, show_cursor, title, alpha, log2fc_threshold):
    d = DESeqFilter("tests/test_files/test_deseq.csv")
    d.volcano_plot(alpha, log2fc_threshold, title, interactive=interactive, show_cursor=show_cursor)
    plt.close('all')


@pytest.mark.parametrize('ref_column,columns,split_plots', [
    ('auto', 'all', False),
    ('cond2', ['cond1', 'cond2'], True),
    ('cond1', 'cond2', False)
])
def test_countfilter_ma_plot_api(ref_column, columns, split_plots):
    c = CountFilter("tests/test_files/counted.csv")
    c.ma_plot(ref_column, columns, split_plots)
    plt.close('all')


def test_countfilter_pairplot_api():
    c = CountFilter("tests/test_files/counted.csv")
    c.pairplot(log2=False)
    c.pairplot(['cond1', 'cond3'], log2=True)
    plt.close('all')


def test_countfilter_clustergram_api():
    c = CountFilter("tests/test_files/counted.csv")
    try:
        c.clustergram()
        c.clustergram(c.columns[0:2], metric='euclidean', linkage='Ward')
        c.clustergram(c.columns[0:2], metric='euclidean', linkage='single')
        with pytest.raises(AssertionError):
            c.clustergram(linkage='invalid')
        with pytest.raises(AssertionError):
            c.clustergram(metric='invalid')
        with pytest.raises(AssertionError):
            c.clustergram(linkage=5)
    finally:
        plt.close('all')


def test_countfilter_box_plot_api():
    c = CountFilter("tests/test_files/counted.csv")
    c.enhanced_box_plot(ylabel='A different label')
    c.enhanced_box_plot(samples=['cond1', 'cond3'], scatter=True)
    plt.close('all')


def test_countfilter_plot_expression_api():
    c = CountFilter("tests/test_files/counted.csv")
    c.plot_expression(['WBGene00007063'], 'all')
    c.plot_expression('WBGene00007063', [[0, 1], ['cond3', 'cond4']])
    c.plot_expression(['WBGene00007064', 'WBGene00044951', 'WBGene00043988', 'WBGene00007066'], ['cond1'])
    c.plot_expression(['WBGene00007064', 'WBGene00044951'], [['cond1'], [1]])
    plt.close('all')


@pytest.mark.parametrize('highlight',
                         [None, {'WBGene00007063', 'WBGene00007064'}, DESeqFilter('tests/test_files/test_deseq.csv')])
@pytest.mark.parametrize('interactive', [True, False])
@pytest.mark.parametrize('show_cursor', [True, False])
@pytest.mark.parametrize('s1,s2,xlabel,ylabel,title', [
    ('cond1', 'cond2', 'auto', 'auto', 'auto'),
    ('cond3', ['cond2', 'cond1', 'cond4'], 'x', 'y', 'title')])
def test_countfilter_scatter_sample_vs_sample_api(s1, s2, xlabel, ylabel, title, interactive, show_cursor, highlight):
    c = CountFilter("tests/test_files/counted.csv")
    c.scatter_sample_vs_sample(s1, s2, xlabel, ylabel, title, interactive=interactive, show_cursor=show_cursor,
                               highlight=highlight)
    plt.close('all')


def test_countfilter_pca_api():
    c = CountFilter("tests/test_files/counted.csv")
    c.filter_low_reads(1)
    try:
        _ = c.pca()
        _ = c.pca(samples=['cond1', ['cond2', 'cond3']], n_components=2, labels=False,
                  power_transform=True)
        with pytest.raises(AssertionError):
            _ = c.pca(n_components=2.0)
        with pytest.raises(AssertionError):
            _ = c.pca(n_components=1)
    finally:
        plt.close('all')


def test_countfilter_enhanced_box_plot_api():
    c = CountFilter("tests/test_files/counted.csv")
    c.box_plot(notch=True, ylabel='A different label')
    c.box_plot(samples=['cond1', 'cond3'], scatter=True)
    plt.close('all')


def test_countfilter_violin_plot_api():
    c = CountFilter("tests/test_files/counted.csv")
    c.violin_plot(ylabel='A different label')
    c.violin_plot(samples=['cond1', 'cond4'])
    plt.close('all')


def _filter_biotype_from_table_tester(filter_obj, truth_protein_coding, truth_pirna):
    protein_coding = filter_obj.filter_biotype_from_ref_table(ref=__biotype_ref__, inplace=False)
    pirna = filter_obj.filter_biotype_from_ref_table('piRNA', ref=__biotype_ref__, inplace=False)
    pirna.df.sort_index(inplace=True)
    protein_coding.df.sort_index(inplace=True)
    truth_protein_coding.sort_index(inplace=True)
    truth_pirna.sort_index(inplace=True)
    assert np.all(truth_protein_coding == protein_coding.df)
    assert np.all(truth_pirna == pirna.df)


def test_htcount_filter_biotype_from_ref_table():
    truth_protein_coding = io.load_table('tests/test_files/counted_biotype_protein_coding.csv', 0)
    truth_pirna = io.load_table('tests/test_files/counted_biotype_piRNA.csv', 0)
    h = CountFilter("tests/test_files/counted_biotype.csv")
    _filter_biotype_from_table_tester(h, truth_protein_coding=truth_protein_coding, truth_pirna=truth_pirna)


def test_htcount_filter_biotype_from_ref_table_opposite():
    truth_no_pirna = io.load_table(r'tests/test_files/counted_biotype_no_piRNA.csv', 0)
    h = CountFilter("tests/test_files/counted_biotype.csv")
    h.filter_biotype_from_ref_table('piRNA', ref=__biotype_ref__, opposite=True, inplace=True)
    h.df.sort_index(inplace=True)
    truth_no_pirna.sort_index(inplace=True)
    assert np.all(h.df == truth_no_pirna)


def test_filter_by_attribute():
    truth = io.load_table('tests/test_files/test_deseq_filter_by_attr1.csv', 0)
    d = DESeqFilter('tests/test_files/test_deseq.csv')
    d_notinplace = d.filter_by_attribute('attribute1', ref=__attr_ref__, inplace=False)
    d.filter_by_attribute('attribute1', ref=__attr_ref__)
    truth.sort_index(inplace=True)
    d.df.sort_index(inplace=True)
    d_notinplace.df.sort_index(inplace=True)
    assert np.all(truth == d.df)
    assert np.all(truth == d_notinplace.df)


def test_filter_by_attribute_from_string(monkeypatch):
    monkeypatch.setattr('builtins.input', lambda x: 'attribute1\nattribute2\n')
    union_truth = io.load_table('tests/test_files/counted_filter_by_bigtable_union_truth.csv', 0)
    h = CountFilter('tests/test_files/counted_filter_by_bigtable.csv')
    assert np.all(union_truth.sort_index() == h.filter_by_attribute(mode='union',
                                                                    ref=__attr_ref__,
                                                                    inplace=False).df.sort_index())

    monkeypatch.setattr('builtins.input', lambda x: 'attribute1\nattribute2')
    assert np.all(union_truth.sort_index() == h.filter_by_attribute(mode='union',
                                                                    ref=__attr_ref__,
                                                                    inplace=False).df.sort_index())

    monkeypatch.setattr('builtins.input', lambda x: 'attribute1')
    deseq_truth = io.load_table('tests/test_files/test_deseq_filter_by_attr1.csv', 0)
    d = DESeqFilter('tests/test_files/test_deseq.csv')
    assert np.all(
        deseq_truth.sort_index() == d.filter_by_attribute(ref=__attr_ref__, inplace=False).df.sort_index())


def test_filter_by_attribute_union():
    union_truth = io.load_table('tests/test_files/counted_filter_by_bigtable_union_truth.csv', 0)
    h = CountFilter('tests/test_files/counted_filter_by_bigtable.csv')
    union = h.filter_by_attribute(['attribute1', 'attribute2'], mode='union',
                                  ref=__attr_ref__, inplace=False)
    assert np.all(union.df.sort_index() == union_truth.sort_index())


def test_filter_by_attribute_intersection():
    intersection_truth = io.load_table(r'tests/test_files/counted_filter_by_bigtable_intersect_truth.csv', 0)
    h = CountFilter('tests/test_files/counted_filter_by_bigtable.csv')
    intersection = h.filter_by_attribute(['attribute1', 'attribute2'], mode='intersection',
                                         ref=__attr_ref__,
                                         inplace=False)
    intersection.df.sort_index(inplace=True)
    intersection_truth.sort_index(inplace=True)
    assert np.all(intersection.df == intersection_truth)


def test_filter_by_attribute_invalid_mode():
    h = CountFilter('tests/test_files/counted_filter_by_bigtable.csv')
    with pytest.raises(AssertionError):
        h.filter_by_attribute(['attribute1', 'attribute2'], mode='difference',
                              ref=__attr_ref__)


def test_split_by_attribute():
    h = CountFilter('tests/test_files/counted_filter_by_bigtable.csv')
    attrs = ['attribute2', 'attribute3', 'attribute4', 'attribute1']
    newobjs = h.split_by_attribute(attrs, ref=__attr_ref__)
    assert len(newobjs) == len(attrs)
    for i, attr in enumerate(attrs):
        assert np.all(
            newobjs[i].df.sort_index() == h.filter_by_attribute(attr, ref=__attr_ref__,
                                                                inplace=False).df.sort_index())


def test_split_by_attribute_multiple():
    f = Filter('tests/test_files/test_deseq.csv')
    attrs = ['attribute2', 'attribute3', 'attribute4', 'attribute1']
    newobjs = f.split_by_attribute(attrs, ref=__attr_ref__)
    assert len(newobjs) == len(attrs)
    for i, attr in enumerate(attrs):
        assert np.all(
            newobjs[i].df.sort_index() == f.filter_by_attribute(attr, ref=__attr_ref__,
                                                                inplace=False).df.sort_index())


def test_split_by_attribute_only_one_attribute():
    f = Filter('tests/test_files/test_deseq.csv')
    newobj = f.split_by_attribute(['attribute1'], ref=__attr_ref__)
    assert len(newobj) == 1
    assert np.all(
        newobj[0].df.sort_index() == f.filter_by_attribute('attribute1', ref=__attr_ref__,
                                                           inplace=False).df.sort_index())
    with pytest.raises(AssertionError):
        f.split_by_attribute('attribute1', ref=__attr_ref__)


def test_split_by_attribute_faulty_attributes():
    f = Filter('tests/test_files/test_deseq.csv')
    with pytest.raises(AssertionError):
        f.split_by_attribute(['attribute1', ['attribute2', 'attribute3']],
                             ref=__attr_ref__)
    with pytest.raises(AssertionError):
        f.split_by_attribute(['attribute1', 2], ref=__attr_ref__)


def test_deseq_filter_significant():
    truth = io.load_table("tests/test_files/test_deseq_sig_truth.csv", 0)
    d = DESeqFilter("tests/test_files/test_deseq_sig.csv")
    d.filter_significant(alpha=0.05)
    assert np.all(d.df == truth)


def test_deseq_filter_significant_opposite():
    truth = io.load_table(r'tests/test_files/test_deseq_not_sig_truth.csv', 0).sort_index()
    d = DESeqFilter("tests/test_files/test_deseq_sig.csv")
    d.filter_significant(alpha=0.05, opposite=True)
    d.df.sort_index(inplace=True)
    assert d.df.equals(truth)


def test_filter_top_n_ascending_number():
    truth = io.load_table("tests/test_files/test_deseq_top10.csv", 0).sort_index()
    d = DESeqFilter("tests/test_files/test_deseq.csv")
    d.filter_top_n('padj', 10)
    d.df.sort_index(inplace=True)
    assert np.isclose(truth, d.df).all()


def test_filter_top_n_ascending_text():
    truth = io.load_table("tests/test_files/test_deseq_top10_text_ascend.csv", 0).sort_index()
    d = DESeqFilter("tests/test_files/test_deseq_textcol.csv")
    print(d.sort('textcol', inplace=False).df)
    d.filter_top_n('textcol', 10, True)
    d.df.sort_index(inplace=True)
    assert d.df.equals(truth)


def test_filter_top_n_multiple_columns():
    truth = io.load_table("tests/test_files/test_deseq_textcol_top15_text_basemean.csv", 0).sort_index()
    d = DESeqFilter("tests/test_files/test_deseq_textcol.csv")
    d.filter_top_n(['textcol', 'baseMean'], 15, True)
    d.df.sort_index(inplace=True)
    assert d.df.equals(truth)


def test_filter_top_n_descending_number():
    truth = io.load_table("tests/test_files/test_deseq_bottom7.csv", 0).sort_index()
    d = DESeqFilter("tests/test_files/test_deseq.csv")
    d.filter_top_n('log2FoldChange', 7, False)
    d.df.sort_index(inplace=True)
    assert np.isclose(truth, d.df).all()


def test_filter_top_n_descending_text():
    truth = io.load_table("tests/test_files/test_deseq_bottom10_text_descend.csv", 0).sort_index()
    d = DESeqFilter("tests/test_files/test_deseq_textcol.csv")
    d.filter_top_n('textcol', 10, False)
    d.df.sort_index(inplace=True)
    assert np.all(truth == d.df)


def test_filter_top_n_nonexisting_column():
    d = DESeqFilter("tests/test_files/test_deseq.csv")
    colname = 'somecol'
    with pytest.raises(AssertionError):
        d.filter_top_n(colname, 5)
        d.filter_top_n([d.df.columns[0], colname])
    assert colname not in d.df.columns


def test_deseq_filter_abs_log2_fold_change():
    truth = io.load_table("tests/test_files/test_deseq_fc_4_truth.csv", 0)
    d = DESeqFilter("tests/test_files/test_deseq_fc.csv")
    fc4 = d.filter_abs_log2_fold_change(4, inplace=False)
    fc4.df.sort_index(inplace=True)
    truth.sort_index(inplace=True)
    assert np.all(fc4.df == truth)


def test_deseq_filter_fold_change_direction():
    pos_truth = io.load_table("tests/test_files/test_deseq_fc_pos_truth.csv", 0)
    neg_truth = io.load_table("tests/test_files/test_deseq_fc_neg_truth.csv", 0)
    d = DESeqFilter("tests/test_files/test_deseq_fc.csv")
    pos = d.filter_fold_change_direction('pos', inplace=False)
    neg = d.filter_fold_change_direction('neg', inplace=False)
    assert np.all(pos.df == pos_truth)
    assert np.all(neg.df == neg_truth)


def test_deseq_split_fold_change():
    d = DESeqFilter("tests/test_files/test_deseq_fc.csv")
    pos_truth = io.load_table("tests/test_files/test_deseq_fc_pos_truth.csv", 0)
    neg_truth = io.load_table("tests/test_files/test_deseq_fc_neg_truth.csv", 0)
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


def test_htcount_rpm_negative_threshold():
    h = CountFilter("tests/test_files/counted.csv")
    with pytest.raises(AssertionError):
        h.filter_low_reads(threshold=-3)


def test_htcount_threshold_invalid():
    h = CountFilter("tests/test_files/counted.csv")
    with pytest.raises(AssertionError):
        h.filter_low_reads("5")


def test_htcount_split_by_reads():
    h = CountFilter("tests/test_files/counted.csv")
    high_truth = io.load_table(r"tests/test_files/counted_above60_rpm.csv", 0)
    low_truth = io.load_table(r"tests/test_files/counted_below60_rpm.csv", 0)
    high, low = h.split_by_reads(threshold=60)
    assert np.all(high.df == high_truth)
    assert np.all(low.df == low_truth)


def test_filter_percentile():
    truth = io.load_table(r'tests/test_files/test_deseq_percentile_0.25.csv', 0)
    h = DESeqFilter(r'tests/test_files/test_deseq_percentile.csv')
    h.filter_percentile(0.25, 'padj', inplace=True)
    h.df.sort_index(inplace=True)
    truth.sort_index(inplace=True)
    assert truth.equals(h.df)
    h.filter_percentile(1, 'baseMean')
    assert truth.equals(h.df)
    h.filter_percentile(0, 'padj')
    assert len(h) == 1
    assert h.df['padj'].values == truth['padj'].min()


def test_filter_percentile_bad_input():
    h = DESeqFilter(r'tests/test_files/test_deseq_percentile.csv')
    with pytest.raises(AssertionError):
        h.filter_percentile(-0.2, 'pvalue')
    with pytest.raises(AssertionError):
        h.filter_percentile(1.1, 'baseMean')
    with pytest.raises(AssertionError):
        h.filter_percentile('0.5', 'log2FoldChange')


def test_split_by_percentile():
    truth_below = io.load_table(r'tests/test_files/test_deseq_percentile_0.25.csv', 0)
    truth_above = io.load_table(r'tests/test_files/test_deseq_percentile_0.75.csv', 0)
    h = DESeqFilter(r'tests/test_files/test_deseq_percentile.csv')
    below, above = h.split_by_percentile(0.25, 'padj')
    for i in [truth_below, truth_above, below.df, above.df]:
        i.sort_index(inplace=True)
    assert np.all(truth_below == below.df)
    assert np.all(truth_above == above.df)


def test_htcount_filter_biotype_from_ref_table_multiple():
    truth = io.load_table('tests/test_files/counted_biotype_piRNA_protein_coding.csv', 0)
    h = CountFilter("tests/test_files/counted_biotype.csv")
    both = h.filter_biotype_from_ref_table(['protein_coding', 'piRNA'], ref=__biotype_ref__,
                                           inplace=False)
    both.df.sort_index(inplace=True)
    truth.sort_index(inplace=True)
    assert np.all(truth == both.df)


def test_htcount_filter_biotype_from_ref_table_multiple_opposite():
    truth = io.load_table('tests/test_files/counted_biotype_piRNA_protein_coding_opposite.csv', 0)
    h = CountFilter("tests/test_files/counted_biotype.csv")
    neither = h.filter_biotype_from_ref_table(['protein_coding', 'piRNA'], ref=__biotype_ref__,
                                              inplace=False,
                                              opposite=True)
    neither.df.sort_index(inplace=True)
    truth.sort_index(inplace=True)
    assert np.all(truth == neither.df)


def test_deseq_filter_biotype_from_ref_table():
    truth_protein_coding = io.load_table('tests/test_files/test_deseq_biotype_protein_coding.csv', 0)
    truth_pirna = io.load_table('tests/test_files/test_deseq_biotype_piRNA.csv', 0)
    d = DESeqFilter("tests/test_files/test_deseq_biotype.csv")
    _filter_biotype_from_table_tester(d, truth_protein_coding=truth_protein_coding, truth_pirna=truth_pirna)


def test_deseq_filter_biotype_from_ref_table_opposite():
    truth_no_pirna = io.load_table(r'tests/test_files/test_deseq_biotype_piRNA_opposite.csv', 0)
    d = DESeqFilter("tests/test_files/test_deseq_biotype.csv")
    d.filter_biotype_from_ref_table('piRNA', ref=__biotype_ref__, opposite=True, inplace=True)
    d.df.sort_index(inplace=True)
    truth_no_pirna.sort_index(inplace=True)
    assert np.all(d.df == truth_no_pirna)


def test_deseq_filter_biotype_from_ref_table_multiple():
    truth = io.load_table('tests/test_files/test_deseq_biotype_piRNA_protein_coding.csv', 0)
    d = DESeqFilter("tests/test_files/test_deseq_biotype.csv")
    both = d.filter_biotype_from_ref_table(['protein_coding', 'piRNA'], ref=__biotype_ref__,
                                           inplace=False)
    both.df.sort_index(inplace=True)
    truth.sort_index(inplace=True)
    assert np.all(truth == both.df)


def test_deseq_filter_biotype_from_ref_table_multiple_opposite():
    truth = io.load_table('tests/test_files/test_deseq_biotype_piRNA_protein_coding_opposite.csv', 0)
    d = DESeqFilter("tests/test_files/test_deseq_biotype.csv")
    neither = d.filter_biotype_from_ref_table(['protein_coding', 'piRNA'], ref=__biotype_ref__,
                                              inplace=False,
                                              opposite=True)
    neither.df.sort_index(inplace=True)
    truth.sort_index(inplace=True)
    assert np.all(truth == neither.df)


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
    set1_truth = io.load_table('tests/test_files/test_deseq_set_ops_1_inplace_intersection.csv', 0)
    set2_truth = io.load_table('tests/test_files/test_deseq_set_ops_2_inplace_intersection.csv', 0)
    set1 = DESeqFilter('tests/test_files/test_deseq_set_ops_1.csv')
    set2 = DESeqFilter('tests/test_files/test_deseq_set_ops_2.csv')
    set1_int = set1.__copy__()
    set2_int = set2.__copy__()
    set1_int.intersection(set2, inplace=True)
    set2_int.intersection(set1, inplace=True)
    set1_int.df.sort_index(inplace=True)
    set2_int.df.sort_index(inplace=True)
    set1_truth.sort_index(inplace=True)
    set2_truth.sort_index(inplace=True)

    assert np.all(set1_truth == set1_int.df)
    assert np.all(set2_truth == set2_int.df)


def test_difference_inplace():
    set1_truth = io.load_table('tests/test_files/test_deseq_set_ops_1_inplace_difference.csv', 0)
    set2_truth = io.load_table('tests/test_files/test_deseq_set_ops_2_inplace_difference.csv', 0)
    set1 = DESeqFilter('tests/test_files/test_deseq_set_ops_1.csv')
    set2 = DESeqFilter('tests/test_files/test_deseq_set_ops_2.csv')
    set1_diff = set1.__copy__()
    set2_diff = set2.__copy__()
    set1_diff.difference(set2, inplace=True)
    set2_diff.difference(set1, inplace=True)
    set1_diff.df.sort_index(inplace=True)
    set2_diff.df.sort_index(inplace=True)
    set1_truth.sort_index(inplace=True)
    set2_truth.sort_index(inplace=True)

    assert np.all(set1_truth == set1_diff.df)
    assert np.all(set2_truth == set2_diff.df)


def test_htcount_fold_change():
    truth_num_name = f"Mean of {['cond1_rep1', 'cond1_rep2']}"
    truth_denom_name = f"Mean of {['cond2_rep1', 'cond2_rep2']}"
    truth = io.load_table(r'tests/test_files/counted_fold_change_truth.csv', 0)
    truth = truth.squeeze()
    h = CountFilter(r'tests/test_files/counted_fold_change.csv')
    fc = h.fold_change(['cond1_rep1', 'cond1_rep2'], ['cond2_rep1', 'cond2_rep2'])
    assert truth_num_name == fc.numerator
    assert truth_denom_name == fc.denominator
    assert np.all(np.isclose(fc.df, truth))


def test_fcfilter_filter_abs_fc():
    truth = io.load_table('tests/test_files/fcfilter_abs_fold_change_truth.csv', 0)
    truth = truth.squeeze()
    truth.sort_index(inplace=True)
    f = FoldChangeFilter('tests/test_files/counted_fold_change_truth.csv', 'numer', 'denom')
    f.filter_abs_log2_fold_change(1)
    f.df.sort_index(inplace=True)
    print(f.df.values)
    print(truth.values)
    assert np.all(np.squeeze(f.df.values) == np.squeeze(truth.values))


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


def test_number_filters_gt():
    truth = io.load_table(r'tests/test_files/test_deseq_gt.csv', 0)
    d = DESeqFilter(r'tests/test_files/test_deseq.csv')
    filt_1 = d.number_filters('baseMean', '>', 1000, inplace=False)
    filt_2 = d.number_filters('baseMean', 'GT', 1000, inplace=False)
    filt_3 = d.number_filters('baseMean', 'greater tHAn', 1000, inplace=False)
    filt_1.df.sort_index(inplace=True)
    filt_2.df.sort_index(inplace=True)
    filt_3.df.sort_index(inplace=True)
    truth.sort_index(inplace=True)
    assert np.all(filt_1.df == filt_2.df)
    assert np.all(filt_2.df == filt_3.df)
    assert np.all(np.squeeze(truth) == np.squeeze(filt_1.df))


def test_number_filters_lt():
    truth = io.load_table(r'tests/test_files/test_deseq_lt.csv', 0)
    d = DESeqFilter(r'tests/test_files/test_deseq.csv')
    filt_1 = d.number_filters('lfcSE', 'Lesser than', 0.2, inplace=False)
    filt_2 = d.number_filters('lfcSE', 'lt', 0.2, inplace=False)
    filt_3 = d.number_filters('lfcSE', '<', 0.2, inplace=False)
    filt_1.df.sort_index(inplace=True)
    filt_2.df.sort_index(inplace=True)
    filt_3.df.sort_index(inplace=True)
    truth.sort_index(inplace=True)
    assert np.all(filt_1.df == filt_2.df)
    assert np.all(filt_2.df == filt_3.df)
    assert np.all(np.squeeze(truth) == np.squeeze(filt_1.df))


def test_number_filters_eq():
    truth = io.load_table(r'tests/test_files/counted_eq.csv', 0)
    d = CountFilter(r'tests/test_files/counted.csv')
    filt_1 = d.number_filters('cond2', 'eQ', 0, inplace=False)
    filt_2 = d.number_filters('cond2', '=', 0, inplace=False)
    filt_3 = d.number_filters('cond2', 'Equals', 0, inplace=False)
    filt_1.df.sort_index(inplace=True)
    filt_2.df.sort_index(inplace=True)
    filt_3.df.sort_index(inplace=True)
    truth.sort_index(inplace=True)
    assert np.all(filt_1.df == filt_2.df)
    assert np.all(filt_2.df == filt_3.df)
    assert np.all(np.squeeze(truth) == np.squeeze(filt_1.df))


def test_number_filters_invalid_input():
    d = CountFilter(r'tests/test_files/counted.csv')
    with pytest.raises(AssertionError):
        d.number_filters('Cond2', 'lt', 5)
    with pytest.raises(AssertionError):
        d.number_filters('cond2', 'contains', 6)
    with pytest.raises(AssertionError):
        d.number_filters('cond2', 'equals', '55')


def test_text_filters_eq():
    truth = io.load_table('tests/test_files/text_filters_eq.csv', 0)
    d = CountFilter('tests/test_files/text_filters.csv')
    filt_1 = d.text_filters('class', 'eQ', 'B', inplace=False)
    filt_2 = d.text_filters('class', '=', 'B', inplace=False)
    filt_3 = d.text_filters('class', 'Equals', 'B', inplace=False)
    filt_1.df.sort_index(inplace=True)
    filt_2.df.sort_index(inplace=True)
    filt_3.df.sort_index(inplace=True)
    truth.sort_index(inplace=True)
    assert np.all(filt_1.df == filt_2.df)
    assert np.all(filt_2.df == filt_3.df)
    assert np.all(np.squeeze(truth) == np.squeeze(filt_1.df))


def test_text_filters_ct():
    truth = io.load_table('tests/test_files/text_filters_ct.csv', 0)
    d = CountFilter('tests/test_files/text_filters.csv')
    filt_1 = d.text_filters('name', 'ct', 'C3.', inplace=False)
    filt_2 = d.text_filters('name', 'IN', 'C3.', inplace=False)
    filt_3 = d.text_filters('name', 'contaiNs', 'C3.', inplace=False)
    filt_1.df.sort_index(inplace=True)
    filt_2.df.sort_index(inplace=True)
    filt_3.df.sort_index(inplace=True)
    truth.sort_index(inplace=True)
    assert np.all(filt_1.df == filt_2.df)
    assert np.all(filt_2.df == filt_3.df)
    assert np.all(np.squeeze(truth) == np.squeeze(filt_1.df))


def test_text_filters_sw():
    truth = io.load_table('tests/test_files/text_filters_sw.csv', 0)
    d = CountFilter('tests/test_files/text_filters.csv')
    filt_1 = d.text_filters('name', 'sw', '2R', inplace=False)
    filt_2 = d.text_filters('name', 'Starts With', '2R', inplace=False)
    filt_1.df.sort_index(inplace=True)
    filt_2.df.sort_index(inplace=True)
    truth.sort_index(inplace=True)
    print(filt_1.df)
    assert np.all(filt_1.df == filt_2.df)
    assert np.all(np.squeeze(truth) == np.squeeze(filt_1.df))


def test_text_filters_ew():
    truth = io.load_table('tests/test_files/text_filters_ew.csv', 0)
    d = CountFilter('tests/test_files/text_filters.csv')
    filt_1 = d.text_filters('name', 'ew', '3', inplace=False)
    filt_2 = d.text_filters('name', 'ends With', '3', inplace=False)
    filt_1.df.sort_index(inplace=True)
    filt_2.df.sort_index(inplace=True)
    truth.sort_index(inplace=True)
    print(filt_1.df)
    assert np.all(filt_1.df == filt_2.df)
    assert np.all(np.squeeze(truth) == np.squeeze(filt_1.df))


def test_text_filters_invalid_input():
    d = CountFilter(r'tests/test_files/counted.csv')
    with pytest.raises(AssertionError):
        d.text_filters('Cond2', 'contains', '5')
    with pytest.raises(AssertionError):
        d.text_filters('cond2', 'lt', '6')
    with pytest.raises(AssertionError):
        d.text_filters('cond2', 'equals', 55)


def test_count_filter_from_folder():
    counted_fname = '__allexpr_temporary_testfile.csv'
    uncounted_fname = '__allfeature_temporary_testfile.csv'

    truth_all_expr = io.load_table('tests/test_files/test_count_from_folder_all_expr.csv', 0).sort_index()
    truth_all_feature = io.load_table('tests/test_files/test_count_from_folder_all_feature.csv', 0).sort_index()
    counts = CountFilter.from_folder_htseqcount('tests/test_files/test_count_from_folder', norm_to_rpm=False,
                                                save_csv=True,
                                                counted_fname=counted_fname, uncounted_fname=uncounted_fname)

    try:
        print('counts:')
        print(counts.df.sort_index())
        print('truth:')
        print(truth_all_expr)
        assert np.all(np.isclose(counts.df.sort_index(), truth_all_expr, atol=0, rtol=0.0001))

        all_feature = io.load_table(f'tests/test_files/test_count_from_folder/{uncounted_fname}', 0).sort_index()
        assert all_feature.equals(truth_all_feature)

    finally:
        os.remove('tests/test_files/test_count_from_folder/__allexpr_temporary_testfile.csv')
        os.remove('tests/test_files/test_count_from_folder/__allfeature_temporary_testfile.csv')


def test_count_filter_from_folder_save_without_suffix():
    counted_fname = '__allexpr_temporary_testfile.csv'
    uncounted_fname = '__allfeature_temporary_testfile.csv'
    try:
        _ = CountFilter.from_folder_htseqcount('tests/test_files/test_count_from_folder', norm_to_rpm=False,
                                               save_csv=True,
                                               counted_fname=counted_fname, uncounted_fname=uncounted_fname)
    finally:
        os.remove(f'tests/test_files/test_count_from_folder/{counted_fname}')
        os.remove(f'tests/test_files/test_count_from_folder/{uncounted_fname}')


def test_count_filter_from_folder_norm():
    truth_norm = io.load_table('tests/test_files/test_count_from_folder_norm.csv', 0)
    counts_norm = CountFilter.from_folder_htseqcount('tests/test_files/test_count_from_folder', norm_to_rpm=True,
                                                     save_csv=False)
    assert np.all(np.isclose(counts_norm.df, truth_norm, atol=0, rtol=0.0001))


def test_biotypes_from_ref_table():
    truth = io.load_table('tests/test_files/biotypes_truth.csv', 0).sort_index()
    c = CountFilter('tests/test_files/counted_biotype.csv')
    df = c.biotypes_from_ref_table(ref=__biotype_ref__).sort_index()
    print('\n')
    print(df, truth)
    assert df.equals(truth)


def test_biotypes_from_ref_table_long_form():
    truth = pd.read_csv('tests/test_files/biotypes_long_format_truth.csv', index_col=0, header=[0, 1]).sort_index()
    truth.columns = ['_'.join(col) for col in truth.columns.values]
    c = CountFilter('tests/test_files/counted_biotype.csv')
    df = c.biotypes_from_ref_table(long_format=True, ref=__biotype_ref__).sort_index()
    assert np.isclose(df, truth, equal_nan=True).all()


def test_filter_by_row_sum():
    truth = io.load_table('tests/test_files/test_filter_row_sum.csv', 0)
    h = CountFilter('tests/test_files/counted.csv')
    h.filter_by_row_sum(29)
    h.df.sort_index(inplace=True)
    truth.sort_index(inplace=True)
    assert np.all(h.df == truth)


def test_sort_inplace():
    c = CountFilter('tests/test_files/counted.csv')
    c.sort(by='cond3', ascending=True, inplace=True)
    assert c.df['cond3'].is_monotonic_increasing


def test_sort_not_inplace():
    c = CountFilter('tests/test_files/counted.csv')
    c_copy = io.load_table('tests/test_files/counted.csv', 0)
    c_sorted = c.sort(by='cond3', ascending=True, inplace=False)
    assert c_sorted.df['cond3'].is_monotonic_increasing
    assert np.all(c.df == c_copy)


def test_sort_by_multiple_columns():
    truth = io.load_table('tests/test_files/counted_sorted_multiple_truth.csv', 0)
    c = CountFilter('tests/test_files/counted.csv')
    c.sort(by=['cond3', 'cond4', 'cond1', 'cond2'], ascending=[True, False, True, False], inplace=True)
    assert np.all(truth == c.df)


def test_sort_with_na_first():
    truth_first = io.load_table('tests/test_files/test_deseq_with_nan_sorted_nanfirst_truth.csv', 0)
    truth_last = io.load_table('tests/test_files/test_deseq_with_nan_sorted_nanlast_truth.csv', 0)
    c = CountFilter('tests/test_files/test_deseq_with_nan.csv')
    c.sort(by='padj', ascending=True, na_position='first', inplace=True)
    assert truth_first.equals(c.df)
    c.sort(by='padj', ascending=True, na_position='last', inplace=True)
    assert truth_last.equals(c.df)


def test_sort_descending():
    c = CountFilter('tests/test_files/counted.csv')
    c.sort(by='cond3', ascending=False, inplace=True)
    assert c.df['cond3'].is_monotonic_decreasing


def test_filter_missing_values():
    truth = io.load_table('tests/test_files/test_deseq_with_nan_all_removed.csv', 0)
    f = Filter('tests/test_files/test_deseq_with_nan.csv')
    f.filter_missing_values()
    assert np.all(f.df.sort_index() == truth.sort_index())


def test_filter_missing_values_foldchangefilter():
    truth = io.load_table('tests/test_files/fc_1_nan_removed.csv', 0, squeeze=True)
    f = FoldChangeFilter('tests/test_files/fc_1_nan.csv', 'num', 'denom')
    res_all = f.filter_missing_values(inplace=False)
    assert truth.equals(res_all.df)
    res_foldchange = f.filter_missing_values(inplace=False)
    assert truth.equals(res_foldchange.df)


def test_filter_missing_values_one_columns():
    truth = io.load_table('tests/test_files/test_deseq_with_nan_basemean_removed.csv', 0)
    f = Filter('tests/test_files/test_deseq_with_nan.csv')
    f.filter_missing_values('baseMean')
    print(f.df.sort_index())
    print(truth.sort_index())
    print(f.df.sort_index() == truth.sort_index())
    assert truth.equals(f.df)


def test_filter_missing_values_multiple_columns():
    truth = io.load_table('tests/test_files/test_deseq_with_nan_basemean_pvalue_removed.csv', 0)
    f = Filter('tests/test_files/test_deseq_with_nan.csv')
    f.filter_missing_values(['baseMean', 'pvalue'])
    print(f.df.sort_index())
    print(truth.sort_index())
    print(f.df.sort_index() == truth.sort_index())
    assert truth.equals(f.df)


def test_filter_missing_values_invalid_type():
    f = Filter('tests/test_files/test_deseq_with_nan.csv')
    with pytest.raises(TypeError):
        f.filter_missing_values(columns={'baseMean': True, 'log2FolgChange': False})


def test_filter_missing_values_nonexistent_column():
    f = Filter('tests/test_files/test_deseq_with_nan.csv')
    with pytest.raises(AssertionError):
        f.filter_missing_values('pval')
    with pytest.raises(AssertionError):
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
    with pytest.raises(AssertionError):
        pl.remove_last_function()


def test_pipeline_apply_empty_pipeline():
    pl = Pipeline()
    d = DESeqFilter('tests/test_files/test_deseq.csv')
    with pytest.raises(AssertionError):
        pl.apply_to(d)


def test_pipeline_apply_to():
    pl = Pipeline('deseqfilter')
    pl.add_function('filter_significant', 10 ** -70, opposite=True)
    deseq = DESeqFilter('tests/test_files/test_deseq.csv')
    deseq_truth = deseq.__copy__()
    deseq_truth.filter_significant(10 ** -70, opposite=True)
    deseq_pipelined = pl.apply_to(deseq, inplace=False)
    pl.apply_to(deseq)
    deseq.sort('log2FoldChange')
    deseq_truth.sort('log2FoldChange')
    deseq_pipelined.sort('log2FoldChange')
    assert np.all(deseq.df == deseq_truth.df)
    assert np.all(deseq_pipelined.df == deseq_truth.df)

    pl2 = Pipeline('countfilter')
    pl2.add_function(Filter.filter_biotype_from_ref_table, biotype='protein_coding',
                     ref=__biotype_ref__)
    cnt = CountFilter('tests/test_files/counted.csv')
    cnt_truth = cnt.__copy__()
    cnt_truth.filter_biotype_from_ref_table('protein_coding', ref=__biotype_ref__)
    cnt_pipelined = pl2.apply_to(cnt, inplace=False)
    pl2.apply_to(cnt, inplace=True)
    cnt.sort(cnt.columns[0])
    cnt_truth.sort(cnt.columns[0])
    cnt_pipelined.sort(cnt.columns[0])
    assert np.all(cnt.df == cnt_truth.df)
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


def test_pipeline_apply_to_invalid_object():
    pl = Pipeline('deseqfilter')
    pl.add_function(DESeqFilter.filter_significant, alpha=10 ** -70)
    cnt = io.load_table('tests/test_files/counted.csv', 0)
    with pytest.raises(AssertionError):
        pl.apply_to(cnt)


def test_pipeline_init_invalid_filter_type():
    with pytest.raises(AssertionError):
        Pipeline(filter_type='otherFilter')

    class otherFilter:
        def __init__(self):
            self.value = 'value'
            self.othervalue = 'othervalue'

    with pytest.raises(AssertionError):
        Pipeline(filter_type=otherFilter)
    with pytest.raises(AssertionError):
        Pipeline(filter_type=max)
    with pytest.raises(AssertionError):
        Pipeline(filter_type=5)


def test_pipeline_add_function_out_of_module():
    pl = Pipeline()
    with pytest.raises(AssertionError):
        pl.add_function(len)

    with pytest.raises(AssertionError):
        pl.add_function(np.sort, arg='val')


def test_pipeline_add_function_invalid_type():
    pl = Pipeline()
    with pytest.raises(AssertionError):
        pl.add_function('string', arg='val')


def test_pipeline_add_function_mismatch_filter_type():
    pl_deseq = Pipeline('DESeqFilter')
    pl_deseq.add_function(CountFilter.filter_biotype_from_ref_table, biotype='protein_coding',
                          ref=__biotype_ref__)
    with pytest.raises(AssertionError):
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
    print(res)
    print(res_truth)
    assert res.keys() == res_truth.keys()
    assert res['biotypes_from_ref_table_1'].equals(res_truth['biotypes_from_ref_table_1'])
    assert type(res['volcano_plot_1']) == type(res_truth['volcano_plot_1'])


def test_pipeline_apply_to_with_plot_inplace():
    _get_pipeline_with_plot(inplace=True)


def test_pipeline_apply_to_with_plot_not_inplace():
    _get_pipeline_with_plot(inplace=False)


def test_pipeline_apply_to_with_split_function():
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


def test_pipeline_apply_to_with_split_function_inplace_raise_error():
    pl = Pipeline('DESeqFilter')
    pl.add_function('filter_missing_values')
    pl.add_function(DESeqFilter.filter_significant, 0.05, opposite=True)
    pl.add_function(DESeqFilter.split_fold_change_direction)
    with pytest.raises(AssertionError):
        pl.apply_to(DESeqFilter('tests/test_files/test_deseq.csv'), inplace=True)


def test_pipeline_apply_to_multiple_splits():
    pl_c = Pipeline('CountFilter')
    pl_c.add_function(CountFilter.filter_top_n, by='cond2', n=2, opposite=True)
    pl_c.add_function(CountFilter.split_hdbscan, min_cluster_size=3, return_probabilities=True)
    pl_c.add_function(CountFilter.split_kmedoids, n_clusters=2, random_seed=42)
    pl_c.add_function(CountFilter.split_by_reads, 15)

    c = CountFilter('tests/test_files/counted.csv')
    c_pipeline_res, c_pipeline_dict = pl_c.apply_to(c, inplace=False)
    c_res = c.filter_top_n(by='cond2', n=2, opposite=True, inplace=False)
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


def test_pipeline_apply_to_filter_normalize_split_plot():
    scaling_factor_path = 'tests/test_files/big_scaling_factors.csv'
    pl_c = Pipeline('CountFilter')
    pl_c.add_function(CountFilter.normalize_with_scaling_factors, scaling_factor_path)
    pl_c.add_function('biotypes_from_ref_table', ref=__biotype_ref__)
    pl_c.add_function(CountFilter.filter_top_n, by='cond3rep1', n=1500, opposite=True)
    pl_c.add_function(CountFilter.split_hdbscan, min_cluster_size=40, return_probabilities=True)
    pl_c.add_function(CountFilter.filter_low_reads, threshold=10)
    pl_c.add_function(CountFilter.clustergram)
    pl_c.add_function(CountFilter.split_kmedoids, n_clusters=[2, 7], random_seed=42, n_init=1, max_iter=100)
    pl_c.add_function(CountFilter.sort, by='cond2rep3')
    pl_c.add_function(CountFilter.biotypes_from_ref_table, 'long', __biotype_ref__)

    c = CountFilter('tests/test_files/big_counted.csv')
    c_pipeline_res, c_pipeline_dict = pl_c.apply_to(c, inplace=False)
    c_dict = dict()
    c.normalize_with_scaling_factors(scaling_factor_path)
    c_dict['biotypes_from_ref_table_1'] = c.biotypes_from_ref_table(ref=__biotype_ref__)
    c_res = c.filter_top_n(by='cond3rep1', n=1500, opposite=True, inplace=False)
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


def test_gap_statistic_api():
    c = CountFilter('tests/test_files/big_counted.csv')
    res = c.split_kmeans(n_clusters='gap')
    assert isinstance(res, tuple)


def test_silhouette_api():
    c = CountFilter('tests/test_files/big_counted.csv')
    res = c.split_kmeans(n_clusters='silhouette')
    assert isinstance(res, tuple)


def test_split_kmeans_api():
    c = CountFilter('tests/test_files/big_counted.csv')
    res = c.split_kmeans(n_clusters=4, n_init=2, max_iter=50)
    assert isinstance(res, tuple)
    assert len(res) == 4
    _test_correct_clustering_split(c, res)
    res2 = c.split_kmeans(n_clusters=[2, 3, 7], n_init=1, max_iter=10, random_seed=42, plot_style='std_bar')
    assert isinstance(res2, tuple)
    assert np.all([isinstance(i, tuple) for i in res2])
    [_test_correct_clustering_split(c, i) for i in res2]


def test_split_hierarchical_api():
    c = CountFilter('tests/test_files/big_counted.csv')
    c.filter_top_n(by='cond1rep1', n=2000)
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

    with pytest.raises(AssertionError):
        c.split_hierarchical(5, metric='badinput')
    with pytest.raises(AssertionError):
        c.split_hierarchical(4, linkage='badinput')


def test_split_hdbscan_api():
    c = CountFilter('tests/test_files/big_counted.csv').filter_top_n('cond1rep1', n=2000, inplace=False)
    res = c.split_hdbscan(100)
    _test_correct_clustering_split(c, res, True)
    res2 = c.split_hdbscan(4, 5, 'manhattan', 0.2, 'leaf', plot_style='std_area', return_probabilities=True)
    assert isinstance(res2, list)
    assert isinstance(res2[0], tuple)
    assert isinstance(res2[1], np.ndarray)
    with pytest.raises(AssertionError):
        c.split_hdbscan(min_cluster_size=1)
    with pytest.raises(AssertionError):
        c.split_hdbscan(c.shape[0] + 1)


def test_split_clicom_api():
    c = CountFilter('tests/test_files/big_counted.csv')
    c.filter_top_n(by='cond1rep1', n=300)
    res = c.split_clicom({'method': 'hdbscan', 'min_cluster_size': [20, 75, 40], 'metric': ['ys1', 'yr1', 'spearman']},
                         {'method': 'hierarchical', 'n_clusters': [7, 12], 'metric': ['euclidean', 'jackknife', 'yr1'],
                          'linkage': ['average', 'ward']},
                         {'method': 'kmedoids', 'n_clusters': [7, 16], 'metric': 'spearman'},
                         power_transform=(False, True), evidence_threshold=0.5, min_cluster_size=5)
    _test_correct_clustering_split(c, res, True)


def test_split_kmedoids_api():
    c = CountFilter('tests/test_files/big_counted.csv')
    c.filter_top_n(by='cond1rep1', n=2000)
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
    assert np.all(truth['significant'] == res['significant'])
    assert np.isclose(truth.iloc[:, :-1], res.iloc[:, :-1]).all()


def test_filter_save_csv():
    d = DESeqFilter('tests/test_files/test_deseq_with_nan.csv')
    d.filter_missing_values()
    d.save_csv()
    pth = Path('tests/test_files/test_deseq_with_nan_removemissingvals.csv')
    assert pth.exists()
    d_loaded = DESeqFilter(pth)
    assert np.isclose(d_loaded.df, d.df).all()
    pth.unlink()
    d_sig = d.filter_significant(inplace=False)
    d_sig.save_csv(alt_filename='d_significant')
    d_sig.save_csv(alt_filename='d_significant_with_suffix.csv')
    pth_sig = Path('tests/test_files/d_significant.csv')
    pth_sig_suffix = Path('tests/test_files/d_significant_with_suffix.csv')
    assert pth_sig.exists()
    assert pth_sig_suffix.exists()
    assert np.isclose(DESeqFilter(pth_sig).df, d_sig.df).all()
    assert np.isclose(DESeqFilter(pth_sig_suffix).df, d_sig.df).all()
    pth_sig.unlink()
    pth_sig_suffix.unlink()


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


def test_set_ops_wrong_type():
    with pytest.raises(TypeError):
        Filter('tests/test_files/counted.csv')._set_ops([{1, 2, 3}, 'string'], 'set')


@pytest.mark.parametrize('sample_list,truth_path',
                         [([['cond1', 'cond2'], ('cond3', 'cond4')], 'tests/test_files/counted_averaged_1.csv'),
                          ([['cond1'], ['cond2', 'cond3', 'cond4']], 'tests/test_files/counted_averaged_2.csv'),
                          (['cond1', ['cond2', 'cond3', 'cond4']], 'tests/test_files/counted_averaged_2.csv'),
                          (['cond1', 'cond2', 'cond3', 'cond4'], 'tests/test_files/counted.csv')])
def test_avg_subsamples(sample_list, truth_path):
    counts = CountFilter('tests/test_files/counted.csv')
    truth = io.load_table(truth_path, 0)
    res = counts._avg_subsamples(sample_list)

    assert np.all(res.columns == truth.columns)
    assert np.all(res.index == truth.index)
    assert np.isclose(res, truth, atol=0, rtol=0.001).all()


@pytest.mark.parametrize('input_file,expected_triplicates',
                         [('tests/test_files/counted.csv', [['cond1', 'cond2', 'cond3'], ['cond4']]),
                          ('tests/test_files/counted_6cols.csv',
                           [['cond1', 'cond2', 'cond3'], ['cond4', 'cond5', 'cond6']])])
def test_triplicates(input_file, expected_triplicates):
    counted = CountFilter(input_file)
    assert counted.triplicates == expected_triplicates


def _log2_plus1(df: pd.DataFrame):
    return np.log2(df + 1)


def _log10_plus1(df: pd.DataFrame):
    return np.log10(df + 1)


def _box_cox_plus1(df: pd.DataFrame):
    return PowerTransformer(method='box-cox').fit_transform(df + 1)


def _multiply_by_3_reduce_2(df: pd.DataFrame):
    return (df * 3) - 2


@pytest.mark.parametrize('filter_obj,columns',
                         [(Filter('tests/test_files/counted_6cols.csv'), 'all'),
                          (DESeqFilter('tests/test_files/test_deseq.csv'), 'baseMean'),
                          (CountFilter('tests/test_files/counted.csv'), ['cond1', 'cond4'])])
@pytest.mark.parametrize('function,kwargs,matching_function',
                         [('log2', {}, _log2_plus1),
                          ('log10', {}, _log10_plus1),
                          (_log2_plus1, {}, _log2_plus1),
                          ('box-cox', {}, _box_cox_plus1),
                          (lambda x, mult, red: (x * mult) - red, {'mult': 3, 'red': 2}, _multiply_by_3_reduce_2)])
def test_transform(filter_obj, columns, function, kwargs, matching_function):
    truth = filter_obj.__copy__()
    if columns == 'all':
        truth.df = matching_function(truth.df)
    else:
        if len(truth.df[columns].shape) == 1:
            truth.df[columns] = matching_function(truth.df[columns].to_frame())
        else:
            truth.df[columns] = matching_function(truth.df[columns])
    filter_obj.transform(function, columns, **kwargs)

    filter_obj.df.equals(truth.df)


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


@pytest.mark.parametrize('components,gene_fraction,truth_paths', [
    (1, 0.32, ['tests/test_files/counted_pc1_0.32_top.csv', 'tests/test_files/counted_pc1_0.32_bottom.csv'])
])
def test_split_by_principal_components(components, gene_fraction, truth_paths):
    truth = [CountFilter(pth) for pth in truth_paths]
    c = CountFilter('tests/test_files/counted.csv')
    c.filter_low_reads(1)
    res = c.split_by_principal_components(components, gene_fraction)
    assert len(res) == len(truth)
    for i in range(len(truth)):
        assert res[i].df.sort_index().equals(truth[i].df.sort_index())


@pytest.mark.parametrize('ids,mode,truth_path', [
    ('kegg_id2', 'union', 'tests/test_files/counted_filter_by_kegg_truth_1.csv'),
    (['kegg_id2'], 'intersection', 'tests/test_files/counted_filter_by_kegg_truth_1.csv'),
    (['kegg_id3', 'kegg_id1'], 'union', 'tests/test_files/counted_filter_by_kegg_truth_2.csv'),
    (['kegg_id1', 'kegg_id3'], 'intersection', 'tests/test_files/counted_filter_by_kegg_truth_3.csv'),
    (['kegg_id1', 'kegg_id2', 'kegg_id3'], 'union', 'tests/test_files/counted_filter_by_kegg_truth_4.csv'),
    (['kegg_id1', 'kegg_id2', 'kegg_id3'], 'intersection', 'tests/test_files/counted_filter_by_kegg_truth_5.csv'),

])
def test_filter_by_kegg_annotations(monkeypatch, ids, mode, truth_path):
    truth = io.load_table(truth_path, index_col=0)

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

    f = Filter('tests/test_files/counted.csv')
    res = f.filter_by_kegg_annotations(ids, mode, gene_id_type='KEGG', inplace=False)

    try:
        assert res.df.sort_index().equals(truth.sort_index())
    except Exception as e:
        print(res.df)
        print(truth)
        raise e


@pytest.mark.parametrize('ids,mode,truth_path', [
    ('go_id2', 'union', 'tests/test_files/counted_filter_by_kegg_truth_1.csv'),
    (['go_id2'], 'intersection', 'tests/test_files/counted_filter_by_kegg_truth_1.csv'),
    (['go_id3', 'go_id1'], 'union', 'tests/test_files/counted_filter_by_kegg_truth_2.csv'),
    (['go_id1', 'go_id3'], 'intersection', 'tests/test_files/counted_filter_by_kegg_truth_3.csv'),
    (['go_id1', 'go_id2', 'go_id3'], 'union', 'tests/test_files/counted_filter_by_kegg_truth_4.csv'),
    (['go_id1', 'go_id2', 'go_id3'], 'intersection', 'tests/test_files/counted_filter_by_kegg_truth_5.csv'),

])
def test_filter_by_go_annotations(monkeypatch, ids, mode, truth_path):
    truth = io.load_table(truth_path, index_col=0)

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

    f = Filter('tests/test_files/counted.csv')
    res = f.filter_by_go_annotations(ids, mode, gene_id_type='WormBase', inplace=False)

    try:
        assert res.df.sort_index().equals(truth.sort_index())
    except Exception as e:
        print(res.df)
        print(truth)
        raise e


@pytest.mark.parametrize("pth,cols,truth_pth", [
    ('tests/test_files/counted.csv', 'cond2', 'tests/test_files/counted_drop_columns_truth_1.csv'),
    ('tests/test_files/counted.csv', ['cond4', 'cond2'], 'tests/test_files/counted_drop_columns_truth_2.csv'),
])
def test_drop_columns(pth, cols, truth_pth):
    obj = Filter(pth)
    truth = Filter(truth_pth)
    obj.drop_columns(cols)
    assert obj.df.sort_index().equals(truth.df.sort_index())


@pytest.mark.parametrize('comparisons,expected_paths,script_path', [
    ([('replicate', 'rep2', 'rep3')], ['tests/test_files/deseq2_tests/case1/DESeq2_replicate_rep2_vs_rep3_truth.csv'],
     'tests/test_files/deseq2_tests/case1/expected_deseq_script_1.R'),
    ([('condition', 'cond2', 'cond1'), ('condition', 'cond3', 'cond2')],
     ['tests/test_files/deseq2_tests/case2/DESeq2_condition_cond2_vs_cond1_truth.csv',
      'tests/test_files/deseq2_tests/case2/DESeq2_condition_cond3_vs_cond2_truth.csv'],
     'tests/test_files/deseq2_tests/case2/expected_deseq_script_2.R')
])
def test_differential_expression_deseq2(monkeypatch, comparisons, expected_paths, script_path):
    outdir = 'tests/test_files/deseq2_tests/outdir'
    sample_table_path = 'tests/test_files/test_design_matrix.csv'
    truth = parsing.data_to_tuple([DESeqFilter(file) for file in expected_paths])
    c = CountFilter('tests/test_files/big_counted.csv')

    def mock_run_analysis(data_path, design_mat_path, comps, r_installation_folder):
        assert r_installation_folder == 'auto'
        assert comps == comparisons
        assert CountFilter(data_path) == c
        assert io.load_table(design_mat_path, 0).equals(io.load_table(sample_table_path, 0))

        return Path(script_path).parent

    monkeypatch.setattr(differential_expression, 'run_deseq2_analysis', mock_run_analysis)
    try:
        res = c.differential_expression_deseq2(sample_table_path, comparisons, output_folder=outdir)
        assert sorted(res, key=lambda filter_obj: filter_obj.fname.name) == sorted(truth, key=lambda
            filter_obj: filter_obj.fname.name)
        with open(Path(outdir).joinpath(Path(script_path).name)) as outfile, open(script_path) as truthfile:
            assert outfile.read() == truthfile.read()
    finally:
        unlink_tree(outdir)


@pytest.mark.parametrize('comparisons,expected_paths,script_path', [
    ([('replicate', 'rep2', 'rep3')], ['tests/test_files/limma_tests/case1/LimmaVoom_replicate_rep2_vs_rep3_truth.csv'],
     'tests/test_files/limma_tests/case1/expected_limma_script_1.R'),
    ([('condition', 'cond2', 'cond1'), ('condition', 'cond3', 'cond2')],
     ['tests/test_files/limma_tests/case2/LimmaVoom_condition_cond2_vs_cond1_truth.csv',
      'tests/test_files/limma_tests/case2/LimmaVoom_condition_cond3_vs_cond2_truth.csv'],
     'tests/test_files/limma_tests/case2/expected_limma_script_2.R')
])
def test_differential_expression_limma(monkeypatch, comparisons, expected_paths, script_path):
    outdir = 'tests/test_files/limma_tests/outdir'
    sample_table_path = 'tests/test_files/test_design_matrix.csv'
    truth = parsing.data_to_tuple(
        [DESeqFilter(file, log2fc_col='logFC', padj_col='adj.P.Val') for file in expected_paths])
    c = CountFilter('tests/test_files/big_counted.csv')

    def mock_run_analysis(data_path, design_mat_path, comps, r_installation_folder):
        assert r_installation_folder == 'auto'
        assert comps == comparisons
        assert CountFilter(data_path) == c
        assert io.load_table(design_mat_path, 0).equals(io.load_table(sample_table_path, 0))

        return Path(script_path).parent

    monkeypatch.setattr(differential_expression, 'run_limma_analysis', mock_run_analysis)
    try:
        res = c.differential_expression_limma_voom(sample_table_path, comparisons, output_folder=outdir)
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
    truth = io.load_table(exp_path, 0)
    res = f.filter_duplicate_ids(keep, inplace=False)
    assert res.df.equals(truth)
    f.filter_duplicate_ids(keep)
    assert f.df.equals(truth)


def test_filter_by_row_name():
    names = ['WBGene00044951', 'WBGene00014997', 'WBGene00007069']
    f = Filter('tests/test_files/counted_duplicates.csv')
    truth = io.load_table('tests/test_files/counted_drop_names.csv', 0)
    res = f.filter_by_row_name(names, inplace=False)
    assert res.df.equals(truth)
    f.filter_by_row_name(names)
    assert f.df.equals(truth)
