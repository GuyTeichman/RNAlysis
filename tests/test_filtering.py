import pytest
import numpy as np
import pandas as pd
from pathlib import Path
from rnalysis import utils
from sklearn_extra.cluster import KMedoids
import matplotlib
from rnalysis.filtering import *
from rnalysis.filtering import _KMedoidsIter
import os
from tests import __attr_ref__, __biotype_ref__

matplotlib.use('Agg')


def test_filter_api():
    f = Filter('test_files/uncounted.csv')
    assert f.__repr__() == "Filter of file uncounted.csv"


def test_countfilter_api():
    h = CountFilter('test_files/counted.csv')
    assert h.__repr__() == "CountFilter of file counted.csv"


def test_deseqfilter_api():
    d = DESeqFilter('test_files/test_deseq.csv')
    assert d.__repr__() == "DESeqFilter of file test_deseq.csv"


def test_foldchangefilter_api():
    fc = FoldChangeFilter("test_files/fc_1.csv", 'a', 'b')
    assert fc.__repr__() == "FoldChangeFilter (numerator: 'a', denominator: 'b') of file fc_1.csv"


def test_filter_contains():
    objs = [Filter('test_files/test_deseq.csv'), CountFilter('test_files/counted.csv'),
            DESeqFilter('test_files/counted.csv'), FoldChangeFilter('test_files/fc_1.csv', 'num', 'denom')]
    neither = ['WBGene' + str(i) for i in range(10)]
    for obj in objs:
        for ind in obj.df.index:
            assert ind in obj
        for false_ind in neither:
            assert false_ind not in obj


def test_filter_len():
    objs = [Filter('test_files/test_deseq.csv'), CountFilter('test_files/counted.csv'),
            DESeqFilter('test_files/counted.csv'), FoldChangeFilter('test_files/fc_1.csv', 'num', 'denom')]
    for obj in objs:
        assert len(obj) == obj.df.shape[0]


def test_filter_inplace():
    d = DESeqFilter('test_files/test_deseq_no_nans.csv')
    d_copy = DESeqFilter('test_files/test_deseq_no_nans.csv')
    truth = utils.load_csv('test_files/counted.csv')
    d_inplace_false = d._inplace(truth, opposite=False, inplace=False, suffix='suffix')
    assert np.all(d_inplace_false.df == truth)
    assert np.all(d.df == d_copy.df)
    d._inplace(truth, opposite=False, inplace=True, suffix='other_suffix')
    assert np.all(d.df == truth)


def test_head():
    df = utils.load_csv('test_files/test_deseq.csv', 0)
    d = DESeqFilter('test_files/test_deseq.csv')
    assert np.all(df.head(7) == d.head(7))
    assert np.all(df.head(1) == d.head(1))

    df2 = utils.load_csv('test_files/counted.csv', 0)
    f = Filter('test_files/counted.csv')
    assert np.all(df2.head() == f.head())
    assert np.all(df2.head(1000) == f.head(1000))


def test_tail():
    df = utils.load_csv('test_files/test_deseq.csv', 0)
    d = DESeqFilter('test_files/test_deseq.csv')
    assert np.all(df.tail(7) == d.tail(7))
    assert np.all(df.tail(1) == d.tail(1))

    df2 = utils.load_csv('test_files/counted.csv', 0)
    f = Filter('test_files/counted.csv')
    assert np.all(df2.tail() == f.tail())
    assert np.all(df2.tail(1000) == f.tail(1000))


def test_describe():
    fc_df = utils.load_csv('test_files/fc_1.csv', 0, squeeze=True)
    count_df = utils.load_csv('test_files/counted.csv', 0)
    deseq_df = utils.load_csv('test_files/test_deseq.csv', 0)

    fc = FoldChangeFilter('test_files/fc_1.csv', 'a', 'b')
    count = CountFilter('test_files/counted.csv')
    deseq = DESeqFilter('test_files/test_deseq.csv')
    default_percentiles = (0.01, 0.25, 0.5, 0.75, 0.99)
    deciles = [i / 10 for i in range(1, 10)]
    for filter_obj, df in zip([count, deseq, fc], [count_df, deseq_df, fc_df]):
        assert np.all(filter_obj.describe() == df.describe(percentiles=default_percentiles))
        assert np.all(filter_obj.describe(deciles) == df.describe(percentiles=deciles))


def test_print_features_api():
    count = CountFilter('test_files/counted.csv')
    count.print_features()


def test_from_string(monkeypatch):
    monkeypatch.setattr('builtins.input', lambda x: 'first-word\nsecond word\nthird. word\n')
    assert Filter._from_string("msg") == ["first-word", "second word", "third. word"]
    monkeypatch.setattr('builtins.input', lambda x: 'first-word,second word,third. word')
    assert Filter._from_string("msg", delimiter=',') == ["first-word", "second word", "third. word"]


def test_color_gen():
    gen = Filter._color_gen()
    preset_colors = ['tab:blue', 'tab:red', 'tab:green', 'tab:orange', 'tab:purple', 'tab:brown', 'tab:pink',
                     'tab:gray', 'tab:olive', 'tab:cyan', 'gold', 'maroon', 'mediumslateblue', 'fuchsia',
                     'mediumblue', 'black', 'lawngreen']
    for i in range(150):
        color = next(gen)
        assert color in preset_colors or (isinstance(color, np.ndarray) and len(color) == 3 and
                                          np.max(color) <= 1 and np.min(color) >= 0)


def test_countfilter_normalize_to_rpm():
    truth = utils.load_csv(r"test_files/test_norm_reads_rpm.csv", 0)
    h = CountFilter("test_files/counted.csv")
    not_inplace = h.normalize_to_rpm("test_files/uncounted.csv", inplace=False)
    assert np.isclose(truth, not_inplace.df).all()
    h.normalize_to_rpm("test_files/uncounted.csv")
    assert np.isclose(truth, h.df).all()


def test_countfilter_norm_reads_with_scaling_factors():
    truth = utils.load_csv(r"test_files/test_norm_scaling_factors.csv", 0)
    h = CountFilter("test_files/counted.csv")
    factors = utils.load_csv("test_files/scaling_factors.csv")
    h_norm = h.normalize_with_scaling_factors("test_files/scaling_factors.csv", inplace=False)
    h.normalize_with_scaling_factors(factors)
    assert np.isclose(truth, h.df).all()
    assert h_norm.df.equals(h.df)


def test_filter_low_reads():
    truth = utils.load_csv("test_files/counted_low_rpm_truth.csv", 0)
    h = CountFilter("test_files/counted_low_rpm.csv")
    h.filter_low_reads(threshold=5)
    assert np.isclose(truth, h.df).all()


def test_filter_low_reads_reverse():
    h = CountFilter("test_files/counted.csv")
    low_truth = utils.load_csv(r"test_files/counted_below60_rpm.csv", 0)
    h.filter_low_reads(threshold=60, opposite=True)
    h.df.sort_index(inplace=True)
    low_truth.sort_index(inplace=True)
    print(h.shape)
    print(low_truth.shape)
    print(h.df)
    print(low_truth)

    assert np.all(h.df == low_truth)


def test_deseqfilter_volcano_plot_api():
    d = DESeqFilter("test_files/test_deseq.csv")
    d.volcano_plot()
    d.volcano_plot(alpha=0.000001)
    plt.close('all')


def test_countfilter_pairplot_api():
    c = CountFilter("test_files/counted.csv")
    c.pairplot(log2=False)
    c.pairplot(['cond1', 'cond3'])
    plt.close('all')


def test_countfilter_clustergram_api():
    c = CountFilter("test_files/counted.csv")
    c.clustergram()
    c.clustergram(c.columns[0:2], metric='euclidean', linkage='ward')
    c.clustergram(c.columns[0:2], metric='euclidean', linkage='single')
    with pytest.raises(AssertionError):
        c.clustergram(linkage='invalid')
    with pytest.raises(AssertionError):
        c.clustergram(metric='invalid')
    with pytest.raises(AssertionError):
        c.clustergram(linkage=5)
    plt.close('all')


def test_countfilter_box_plot_api():
    c = CountFilter("test_files/counted.csv")
    c.enhanced_box_plot(ylabel='A different label')
    c.enhanced_box_plot(samples=['cond1', 'cond3'], scatter=True)
    plt.close('all')


def test_countfilter_plot_expression_api():
    c = CountFilter("test_files/counted.csv")
    c.plot_expression('WBGene00007063', {'cond 1 and 2': [0, 1], 'cond3 and 4': ['cond3', 'cond4']})
    c.plot_expression(['WBGene00007064', 'WBGene00044951', 'WBGene00043988', 'WBGene00007066'], {'cond1': ['cond1']})
    c.plot_expression(['WBGene00007064', 'WBGene00044951'], {'cond1': ['cond1'], 'cond2': [1]})
    plt.close('all')


def test_countfilter_scatter_sample_vs_sample_api():
    c = CountFilter("test_files/counted.csv")
    c.scatter_sample_vs_sample('cond1', 'cond2')
    c.scatter_sample_vs_sample('cond3', ['cond2', 'cond1', 'cond4'], highlight={'WBGene00007063', 'WBGene00007064'})
    d = DESeqFilter('test_files/test_deseq.csv').intersection(c, inplace=True)
    c.scatter_sample_vs_sample('cond3', ['cond2', 'cond1', 'cond4'], xlabel='label', title='title', ylabel='ylabel',
                               highlight=d)
    plt.close('all')


def test_countfilter_pca_api():
    c = CountFilter("test_files/counted.csv")
    c.pca()
    c.pca(sample_names=['cond1', 'cond2', 'cond3'], sample_grouping=[1, 1, 2], n_components=2, labels=False)
    with pytest.raises(AssertionError):
        c.pca(n_components=2.0)
    with pytest.raises(AssertionError):
        c.pca(n_components=1)
    plt.close('all')


def test_countfilter_enhanced_box_plot_api():
    c = CountFilter("test_files/counted.csv")
    c.box_plot(notch=True, ylabel='A different label')
    c.box_plot(samples=['cond1', 'cond3'], scatter=True)
    plt.close('all')


def test_countfilter_violin_plot_api():
    c = CountFilter("test_files/counted.csv")
    c.violin_plot(ylabel='A different label')
    c.violin_plot(samples=['cond1', 'cond4'])
    plt.close('all')


def _filter_biotype_tester(filter_obj, truth_protein_coding, truth_pirna):
    protein_coding = filter_obj.filter_biotype(ref=__biotype_ref__, inplace=False)
    pirna = filter_obj.filter_biotype('piRNA', ref=__biotype_ref__, inplace=False)
    pirna.df.sort_index(inplace=True)
    protein_coding.df.sort_index(inplace=True)
    truth_protein_coding.sort_index(inplace=True)
    truth_pirna.sort_index(inplace=True)
    assert np.all(truth_protein_coding == protein_coding.df)
    assert np.all(truth_pirna == pirna.df)


def test_htcount_filter_biotype():
    truth_protein_coding = utils.load_csv('test_files/counted_biotype_protein_coding.csv', 0)
    truth_pirna = utils.load_csv('test_files/counted_biotype_piRNA.csv', 0)
    h = CountFilter("test_files/counted_biotype.csv")
    _filter_biotype_tester(h, truth_protein_coding=truth_protein_coding, truth_pirna=truth_pirna)


def test_htcount_filter_biotype_opposite():
    truth_no_pirna = utils.load_csv(r'test_files/counted_biotype_no_piRNA.csv', 0)
    h = CountFilter("test_files/counted_biotype.csv")
    h.filter_biotype('piRNA', ref=__biotype_ref__, opposite=True, inplace=True)
    h.df.sort_index(inplace=True)
    truth_no_pirna.sort_index(inplace=True)
    assert np.all(h.df == truth_no_pirna)


def test_filter_by_attribute():
    truth = utils.load_csv('test_files/test_deseq_filter_by_attr1.csv', 0)
    d = DESeqFilter('test_files/test_deseq.csv')
    d_notinplace = d.filter_by_attribute('attribute1', ref=__attr_ref__, inplace=False)
    d.filter_by_attribute('attribute1', ref=__attr_ref__)
    truth.sort_index(inplace=True)
    d.df.sort_index(inplace=True)
    d_notinplace.df.sort_index(inplace=True)
    assert np.all(truth == d.df)
    assert np.all(truth == d_notinplace.df)


def test_filter_by_attribute_from_string(monkeypatch):
    monkeypatch.setattr('builtins.input', lambda x: 'attribute1\nattribute2\n')
    union_truth = utils.load_csv('test_files/counted_filter_by_bigtable_union_truth.csv', 0)
    h = CountFilter('test_files/counted_filter_by_bigtable.csv')
    assert np.all(union_truth.sort_index() == h.filter_by_attribute(mode='union',
                                                                    ref=__attr_ref__,
                                                                    inplace=False).df.sort_index())

    monkeypatch.setattr('builtins.input', lambda x: 'attribute1\nattribute2')
    assert np.all(union_truth.sort_index() == h.filter_by_attribute(mode='union',
                                                                    ref=__attr_ref__,
                                                                    inplace=False).df.sort_index())

    monkeypatch.setattr('builtins.input', lambda x: 'attribute1')
    deseq_truth = utils.load_csv('test_files/test_deseq_filter_by_attr1.csv', 0)
    d = DESeqFilter('test_files/test_deseq.csv')
    assert np.all(
        deseq_truth.sort_index() == d.filter_by_attribute(ref=__attr_ref__, inplace=False).df.sort_index())


def test_filter_by_attribute_union():
    union_truth = utils.load_csv('test_files/counted_filter_by_bigtable_union_truth.csv', 0)
    h = CountFilter('test_files/counted_filter_by_bigtable.csv')
    union = h.filter_by_attribute(['attribute1', 'attribute2'], mode='union',
                                  ref=__attr_ref__, inplace=False)
    assert np.all(union.df.sort_index() == union_truth.sort_index())


def test_filter_by_attribute_intersection():
    intersection_truth = utils.load_csv(r'test_files/counted_filter_by_bigtable_intersect_truth.csv', 0)
    h = CountFilter('test_files/counted_filter_by_bigtable.csv')
    intersection = h.filter_by_attribute(['attribute1', 'attribute2'], mode='intersection',
                                         ref=__attr_ref__,
                                         inplace=False)
    intersection.df.sort_index(inplace=True)
    intersection_truth.sort_index(inplace=True)
    assert np.all(intersection.df == intersection_truth)


def test_filter_by_attribute_invalid_mode():
    h = CountFilter('test_files/counted_filter_by_bigtable.csv')
    with pytest.raises(ValueError):
        h.filter_by_attribute(['attribute1', 'attribute2'], mode='difference',
                              ref=__attr_ref__)


def test_split_by_attribute():
    h = CountFilter('test_files/counted_filter_by_bigtable.csv')
    attrs = ['attribute2', 'attribute3', 'attribute4', 'attribute1']
    newobjs = h.split_by_attribute(attrs, ref=__attr_ref__)
    assert len(newobjs) == len(attrs)
    for i, attr in enumerate(attrs):
        assert np.all(
            newobjs[i].df.sort_index() == h.filter_by_attribute(attr, ref=__attr_ref__,
                                                                inplace=False).df.sort_index())


def test_split_by_attribute_multiple():
    f = Filter('test_files/test_deseq.csv')
    attrs = ['attribute2', 'attribute3', 'attribute4', 'attribute1']
    newobjs = f.split_by_attribute(attrs, ref=__attr_ref__)
    assert len(newobjs) == len(attrs)
    for i, attr in enumerate(attrs):
        assert np.all(
            newobjs[i].df.sort_index() == f.filter_by_attribute(attr, ref=__attr_ref__,
                                                                inplace=False).df.sort_index())


def test_split_by_attribute_only_one_attribute():
    f = Filter('test_files/test_deseq.csv')
    newobj = f.split_by_attribute(['attribute1'], ref=__attr_ref__)
    assert len(newobj) == 1
    assert np.all(
        newobj[0].df.sort_index() == f.filter_by_attribute('attribute1', ref=__attr_ref__,
                                                           inplace=False).df.sort_index())
    with pytest.raises(AssertionError):
        f.split_by_attribute('attribute1', ref=__attr_ref__)


def test_split_by_attribute_faulty_attributes():
    f = Filter('test_files/test_deseq.csv')
    with pytest.raises(AssertionError):
        f.split_by_attribute(['attribute1', ['attribute2', 'attribute3']],
                             ref=__attr_ref__)
    with pytest.raises(AssertionError):
        f.split_by_attribute(['attribute1', 2], ref=__attr_ref__)


def test_deseq_filter_significant():
    truth = utils.load_csv("test_files/test_deseq_sig_truth.csv", 0)
    d = DESeqFilter("test_files/test_deseq_sig.csv")
    d.filter_significant(alpha=0.05)
    assert np.all(d.df == truth)


def test_deseq_filter_significant_opposite():
    truth = utils.load_csv(r'test_files/test_deseq_not_sig_truth.csv', 0)
    d = DESeqFilter("test_files/test_deseq_sig.csv")
    d.filter_significant(alpha=0.05, opposite=True)
    d.df.sort_index(inplace=True)
    truth.sort_index(inplace=True)
    truth.fillna(1234567890, inplace=True)
    d.df.fillna(1234567890, inplace=True)
    assert np.all(d.df == truth)


def test_filter_top_n_ascending_number():
    truth = utils.load_csv("test_files/test_deseq_top10.csv", 0)
    d = DESeqFilter("test_files/test_deseq.csv")
    d.filter_top_n('padj', 10)
    d.df.sort_index(inplace=True)
    truth.sort_index(inplace=True)
    assert np.isclose(truth, d.df).all()


def test_filter_top_n_ascending_text():
    truth = utils.load_csv("test_files/test_deseq_top10_text_ascend.csv", 0)
    d = DESeqFilter("test_files/test_deseq_textcol.csv")
    d.filter_top_n('textcol', 10, True)
    d.df.sort_index(inplace=True)
    truth.sort_index(inplace=True)
    assert np.all(truth == d.df)


def test_filter_top_n_descending_number():
    truth = utils.load_csv("test_files/test_deseq_bottom7.csv", 0)
    d = DESeqFilter("test_files/test_deseq.csv")
    d.filter_top_n('log2FoldChange', 7, False)
    d.df.sort_index(inplace=True)
    truth.sort_index(inplace=True)
    assert np.isclose(truth, d.df).all()


def test_filter_top_n_descending_text():
    truth = utils.load_csv("test_files/test_deseq_bottom10_text_descend.csv", 0)
    d = DESeqFilter("test_files/test_deseq_textcol.csv")
    d.filter_top_n('textcol', 10, False)
    d.df.sort_index(inplace=True)
    truth.sort_index(inplace=True)
    assert np.all(truth == d.df)


def test_filter_top_n_nonexisting_column():
    d = DESeqFilter("test_files/test_deseq.csv")
    colname = 'somecol'
    with pytest.raises(AssertionError):
        d.filter_top_n(colname, 5)
        d.filter_top_n([d.df.columns[0], colname])
    assert colname not in d.df.columns


def test_deseq_filter_abs_log2_fold_change():
    truth = utils.load_csv("test_files/test_deseq_fc_4_truth.csv", 0)
    d = DESeqFilter("test_files/test_deseq_fc.csv")
    fc4 = d.filter_abs_log2_fold_change(4, inplace=False)
    fc4.df.sort_index(inplace=True)
    truth.sort_index(inplace=True)
    assert np.all(fc4.df == truth)


def test_deseq_filter_fold_change_direction():
    pos_truth = utils.load_csv("test_files/test_deseq_fc_pos_truth.csv", 0)
    neg_truth = utils.load_csv("test_files/test_deseq_fc_neg_truth.csv", 0)
    d = DESeqFilter("test_files/test_deseq_fc.csv")
    pos = d.filter_fold_change_direction('pos', inplace=False)
    neg = d.filter_fold_change_direction('neg', inplace=False)
    assert np.all(pos.df == pos_truth)
    assert np.all(neg.df == neg_truth)


def test_deseq_split_fold_change():
    d = DESeqFilter("test_files/test_deseq_fc.csv")
    pos_truth = utils.load_csv("test_files/test_deseq_fc_pos_truth.csv", 0)
    neg_truth = utils.load_csv("test_files/test_deseq_fc_neg_truth.csv", 0)
    d = DESeqFilter("test_files/test_deseq_fc.csv")
    pos, neg = d.split_fold_change_direction()
    assert np.all(pos.df == pos_truth)
    assert np.all(neg.df == neg_truth)


def test_intersection():
    intersection_truth = {'WBGene00021375', 'WBGene00044258', 'WBGene00219304', 'WBGene00194708', 'WBGene00018199',
                          'WBGene00019174', 'WBGene00021019', 'WBGene00013816', 'WBGene00045366', 'WBGene00219307',
                          'WBGene00045410', 'WBGene00010100', 'WBGene00077437', 'WBGene00007674', 'WBGene00023036',
                          'WBGene00012648', 'WBGene00022486'}
    set1 = DESeqFilter('test_files/test_deseq_set_ops_1.csv')
    set2 = DESeqFilter('test_files/test_deseq_set_ops_2.csv')

    assert set1.intersection(set2, inplace=False) == intersection_truth


def test_union():
    intersection_truth = {'WBGene00021375', 'WBGene00044258', 'WBGene00219304', 'WBGene00194708', 'WBGene00018199',
                          'WBGene00019174', 'WBGene00021019', 'WBGene00013816', 'WBGene00045366', 'WBGene00219307',
                          'WBGene00045410', 'WBGene00010100', 'WBGene00077437', 'WBGene00007674', 'WBGene00023036',
                          'WBGene00012648', 'WBGene00022486'}
    set2_unique = {'WBGene00018193', 'WBGene00021589', 'WBGene00001118', 'WBGene00010755', 'WBGene00020407',
                   'WBGene00044799', 'WBGene00021654', 'WBGene00012919', 'WBGene00021605'}
    set1_unique = {'WBGene00008447', 'WBGene00021018', 'WBGene00012452', 'WBGene00010507', 'WBGene00022730',
                   'WBGene00012961', 'WBGene00022438', 'WBGene00016635', 'WBGene00044478'}

    set1 = DESeqFilter('test_files/test_deseq_set_ops_1.csv')
    set2 = DESeqFilter('test_files/test_deseq_set_ops_2.csv')
    union_truth = intersection_truth.union(set1_unique.union(set2_unique))
    assert set1.union(set2) == union_truth


def test_difference():
    set2_unique = {'WBGene00018193', 'WBGene00021589', 'WBGene00001118', 'WBGene00010755', 'WBGene00020407',
                   'WBGene00044799', 'WBGene00021654', 'WBGene00012919', 'WBGene00021605'}
    set1_unique = {'WBGene00008447', 'WBGene00021018', 'WBGene00012452', 'WBGene00010507', 'WBGene00022730',
                   'WBGene00012961', 'WBGene00022438', 'WBGene00016635', 'WBGene00044478'}

    set1 = DESeqFilter('test_files/test_deseq_set_ops_1.csv')
    set2 = DESeqFilter('test_files/test_deseq_set_ops_2.csv')

    assert set1.difference(set2, inplace=False) == set1_unique
    assert set2.difference(set1, inplace=False) == set2_unique


def test_symmetric_difference():
    set2_unique = {'WBGene00018193', 'WBGene00021589', 'WBGene00001118', 'WBGene00010755', 'WBGene00020407',
                   'WBGene00044799', 'WBGene00021654', 'WBGene00012919', 'WBGene00021605'}
    set1_unique = {'WBGene00008447', 'WBGene00021018', 'WBGene00012452', 'WBGene00010507', 'WBGene00022730',
                   'WBGene00012961', 'WBGene00022438', 'WBGene00016635', 'WBGene00044478'}

    set1 = DESeqFilter('test_files/test_deseq_set_ops_1.csv')
    set2 = DESeqFilter('test_files/test_deseq_set_ops_2.csv')

    assert set1.symmetric_difference(set2) == set2.symmetric_difference(set1)
    assert set1.symmetric_difference(set2) == set1_unique.union(set2_unique)


def test_deseq_feature_set():
    truth = {'WBGene00008447', 'WBGene00021018', 'WBGene00012452', 'WBGene00010507', 'WBGene00022730',
             'WBGene00012648',
             'WBGene00012961', 'WBGene00022438', 'WBGene00016635', 'WBGene00044478', 'WBGene00021375',
             'WBGene00044258', 'WBGene00219304', 'WBGene00194708', 'WBGene00018199', 'WBGene00022486',
             'WBGene00019174', 'WBGene00021019', 'WBGene00013816', 'WBGene00045366', 'WBGene00219307',
             'WBGene00045410', 'WBGene00010100', 'WBGene00077437', 'WBGene00007674', 'WBGene00023036'}
    d = DESeqFilter('test_files/test_deseq_set_ops_1.csv')
    assert d.index_set == truth


def test_deseq_feature_string():
    truth = {'WBGene00008447', 'WBGene00021018', 'WBGene00012452', 'WBGene00010507', 'WBGene00022730',
             'WBGene00012648',
             'WBGene00012961', 'WBGene00022438', 'WBGene00016635', 'WBGene00044478', 'WBGene00021375',
             'WBGene00044258', 'WBGene00219304', 'WBGene00194708', 'WBGene00018199', 'WBGene00022486',
             'WBGene00019174', 'WBGene00021019', 'WBGene00013816', 'WBGene00045366', 'WBGene00219307',
             'WBGene00045410', 'WBGene00010100', 'WBGene00077437', 'WBGene00007674', 'WBGene00023036'}
    d = DESeqFilter('test_files/test_deseq_set_ops_1.csv')
    assert set(d.index_string.split("\n")) == truth


def test_set_ops_multiple_variable_types():
    set2_unique = {'WBGene00018193', 'WBGene00021589', 'WBGene00001118', 'WBGene00010755', 'WBGene00020407',
                   'WBGene00044799', 'WBGene00021654', 'WBGene00012919', 'WBGene00021605'}
    set1_unique = {'WBGene00008447', 'WBGene00021018', 'WBGene00012452', 'WBGene00010507', 'WBGene00022730',
                   'WBGene00012961', 'WBGene00022438', 'WBGene00016635', 'WBGene00044478'}

    set1 = CountFilter('test_files/test_deseq_set_ops_1.csv')
    set2 = DESeqFilter('test_files/test_deseq_set_ops_2.csv')

    assert set1.symmetric_difference(set2) == set2.symmetric_difference(set1)
    assert set1.symmetric_difference(set2) == set1_unique.union(set2_unique)


def test_htcount_rpm_negative_threshold():
    h = CountFilter("test_files/counted.csv")
    with pytest.raises(AssertionError):
        h.filter_low_reads(threshold=-3)


def test_htcount_threshold_invalid():
    h = CountFilter("test_files/counted.csv")
    with pytest.raises(AssertionError):
        h.filter_low_reads("5")


def test_htcount_split_by_reads():
    h = CountFilter("test_files/counted.csv")
    high_truth = utils.load_csv(r"test_files/counted_above60_rpm.csv", 0)
    low_truth = utils.load_csv(r"test_files/counted_below60_rpm.csv", 0)
    high, low = h.split_by_reads(threshold=60)
    assert np.all(high.df == high_truth)
    assert np.all(low.df == low_truth)


def test_filter_percentile():
    truth = utils.load_csv(r'test_files/test_deseq_percentile_0.25.csv', 0)
    h = DESeqFilter(r'test_files/test_deseq_percentile.csv')
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
    h = DESeqFilter(r'test_files/test_deseq_percentile.csv')
    with pytest.raises(AssertionError):
        h.filter_percentile(-0.2, 'pvalue')
    with pytest.raises(AssertionError):
        h.filter_percentile(1.1, 'baseMean')
    with pytest.raises(AssertionError):
        h.filter_percentile('0.5', 'log2FoldChange')


def test_split_by_percentile():
    truth_below = utils.load_csv(r'test_files/test_deseq_percentile_0.25.csv', 0)
    truth_above = utils.load_csv(r'test_files/test_deseq_percentile_0.75.csv', 0)
    h = DESeqFilter(r'test_files/test_deseq_percentile.csv')
    below, above = h.split_by_percentile(0.25, 'padj')
    for i in [truth_below, truth_above, below.df, above.df]:
        i.sort_index(inplace=True)
    assert np.all(truth_below == below.df)
    assert np.all(truth_above == above.df)


def test_htcount_filter_biotype_multiple():
    truth = utils.load_csv('test_files/counted_biotype_piRNA_protein_coding.csv', 0)
    h = CountFilter("test_files/counted_biotype.csv")
    both = h.filter_biotype(['protein_coding', 'piRNA'], ref=__biotype_ref__,
                            inplace=False)
    both.df.sort_index(inplace=True)
    truth.sort_index(inplace=True)
    assert np.all(truth == both.df)


def test_htcount_filter_biotype_multiple_opposite():
    truth = utils.load_csv('test_files/counted_biotype_piRNA_protein_coding_opposite.csv', 0)
    h = CountFilter("test_files/counted_biotype.csv")
    neither = h.filter_biotype(['protein_coding', 'piRNA'], ref=__biotype_ref__,
                               inplace=False,
                               opposite=True)
    neither.df.sort_index(inplace=True)
    truth.sort_index(inplace=True)
    assert np.all(truth == neither.df)


def test_deseq_filter_biotype():
    truth_protein_coding = utils.load_csv('test_files/test_deseq_biotype_protein_coding.csv', 0)
    truth_pirna = utils.load_csv('test_files/test_deseq_biotype_piRNA.csv', 0)
    d = DESeqFilter("test_files/test_deseq_biotype.csv")
    _filter_biotype_tester(d, truth_protein_coding=truth_protein_coding, truth_pirna=truth_pirna)


def test_deseq_filter_biotype_opposite():
    truth_no_pirna = utils.load_csv(r'test_files/test_deseq_biotype_piRNA_opposite.csv', 0)
    d = DESeqFilter("test_files/test_deseq_biotype.csv")
    d.filter_biotype('piRNA', ref=__biotype_ref__, opposite=True, inplace=True)
    d.df.sort_index(inplace=True)
    truth_no_pirna.sort_index(inplace=True)
    assert np.all(d.df == truth_no_pirna)


def test_deseq_filter_biotype_multiple():
    truth = utils.load_csv('test_files/test_deseq_biotype_piRNA_protein_coding.csv', 0)
    d = DESeqFilter("test_files/test_deseq_biotype.csv")
    both = d.filter_biotype(['protein_coding', 'piRNA'], ref=__biotype_ref__,
                            inplace=False)
    both.df.sort_index(inplace=True)
    truth.sort_index(inplace=True)
    assert np.all(truth == both.df)


def test_deseq_filter_biotype_multiple_opposite():
    truth = utils.load_csv('test_files/test_deseq_biotype_piRNA_protein_coding_opposite.csv', 0)
    d = DESeqFilter("test_files/test_deseq_biotype.csv")
    neither = d.filter_biotype(['protein_coding', 'piRNA'], ref=__biotype_ref__,
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

    set1 = DESeqFilter('test_files/test_deseq_set_ops_1.csv')
    set2 = DESeqFilter('test_files/test_deseq_set_ops_2.csv')
    set3 = {'WBGene00077437', 'WBGene00007674', 'WBGene00023036', 'WBGene00012648', 'WBGene44444444',
            'WBGene99999999',
            'WBGene98765432'}
    union_truth = intersection_truth.union(set1_unique, set2_unique, set3_unique)
    assert set1.union(set2, set3) == union_truth


def test_deseqfilter_intersection_multiple():
    intersection_truth = {'WBGene00077437', 'WBGene00007674', 'WBGene00023036',
                          'WBGene00012648', 'WBGene00022486'}
    set1 = DESeqFilter('test_files/test_deseq_set_ops_1.csv')
    set2 = DESeqFilter('test_files/test_deseq_set_ops_2.csv')
    set3 = {'WBGene00077437', 'WBGene00007674', 'WBGene00023036', 'WBGene00012648', 'WBGene00022486',
            'WBGene99999999',
            'WBGene98765432'}

    assert set1.intersection(set2, set3, inplace=False) == intersection_truth


def test_deseqfilter_difference_multiple():
    set2_unique = {'WBGene00021589', 'WBGene00001118', 'WBGene00010755', 'WBGene00020407',
                   'WBGene00044799', 'WBGene00021654', 'WBGene00012919', 'WBGene00021605'}
    set1_unique = {'WBGene00021018', 'WBGene00012452', 'WBGene00010507', 'WBGene00022730',
                   'WBGene00012961', 'WBGene00022438', 'WBGene00016635', 'WBGene00044478'}

    set1 = DESeqFilter('test_files/test_deseq_set_ops_1.csv')
    set2 = DESeqFilter('test_files/test_deseq_set_ops_2.csv')
    set3 = {'WBGene00018193', 'WBGene00008447', 'WBGene12345678'}

    assert set1.difference(set2, set3, inplace=False) == set1_unique
    assert set2.difference(set3, set1, inplace=False) == set2_unique


def test_intersection_inplace():
    set1_truth = utils.load_csv('test_files/test_deseq_set_ops_1_inplace_intersection.csv', 0)
    set2_truth = utils.load_csv('test_files/test_deseq_set_ops_2_inplace_intersection.csv', 0)
    set1 = DESeqFilter('test_files/test_deseq_set_ops_1.csv')
    set2 = DESeqFilter('test_files/test_deseq_set_ops_2.csv')
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
    set1_truth = utils.load_csv('test_files/test_deseq_set_ops_1_inplace_difference.csv', 0)
    set2_truth = utils.load_csv('test_files/test_deseq_set_ops_2_inplace_difference.csv', 0)
    set1 = DESeqFilter('test_files/test_deseq_set_ops_1.csv')
    set2 = DESeqFilter('test_files/test_deseq_set_ops_2.csv')
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
    truth = utils.load_csv(r'test_files/counted_fold_change_truth.csv', 0)
    truth = truth.squeeze()
    h = CountFilter(r'test_files/counted_fold_change.csv')
    fc = h.fold_change(['cond1_rep1', 'cond1_rep2'], ['cond2_rep1', 'cond2_rep2'])
    assert truth_num_name == fc.numerator
    assert truth_denom_name == fc.denominator
    assert np.all(np.isclose(fc.df, truth))


def test_fcfilter_filter_abs_fc():
    truth = utils.load_csv('test_files/fcfilter_abs_fold_change_truth.csv', 0)
    truth = truth.squeeze()
    truth.sort_index(inplace=True)
    f = FoldChangeFilter('test_files/counted_fold_change_truth.csv', 'numer', 'denom')
    f.filter_abs_log2_fold_change(1)
    f.df.sort_index(inplace=True)
    print(f.df.values)
    print(truth.values)
    assert np.all(np.squeeze(f.df.values) == np.squeeze(truth.values))


def test_fcfilter_fold_change_direction():
    truth_pos = utils.load_csv('test_files/fc_1_pos_fc.csv', 0, squeeze=True)
    truth_neg = utils.load_csv('test_files/fc_1_neg_fc.csv', 0, squeeze=True)
    fc = FoldChangeFilter('test_files/fc_1.csv', 'name', 'name')
    pos = fc.filter_fold_change_direction('pos', inplace=False)
    neg = fc.filter_fold_change_direction('neg', inplace=False)
    assert truth_pos.equals(pos.df)
    assert truth_neg.equals(neg.df)


def test_fcfilter_split_fold_change_direction():
    truth_pos = utils.load_csv('test_files/fc_1_pos_fc.csv', 0, squeeze=True)
    truth_neg = utils.load_csv('test_files/fc_1_neg_fc.csv', 0, squeeze=True)
    fc = FoldChangeFilter('test_files/fc_1.csv', 'name', 'name')
    pos, neg = fc.split_fold_change_direction()
    assert truth_pos.equals(pos.df)
    assert truth_neg.equals(neg.df)


def test_fcfilter_filter_fold_change_direction_bad_input():
    fc = FoldChangeFilter('test_files/fc_1.csv', 'name', 'name')
    with pytest.raises(ValueError):
        fc.filter_fold_change_direction('bad_input')


def test_number_filters_gt():
    truth = utils.load_csv(r'test_files/test_deseq_gt.csv', 0)
    d = DESeqFilter(r'test_files/test_deseq.csv')
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
    truth = utils.load_csv(r'test_files/test_deseq_lt.csv', 0)
    d = DESeqFilter(r'test_files/test_deseq.csv')
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
    truth = utils.load_csv(r'test_files/counted_eq.csv', 0)
    d = CountFilter(r'test_files/counted.csv')
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
    d = CountFilter(r'test_files/counted.csv')
    with pytest.raises(AssertionError):
        d.number_filters('Cond2', 'lt', 5)
    with pytest.raises(AssertionError):
        d.number_filters('cond2', 'contains', 6)
    with pytest.raises(AssertionError):
        d.number_filters('cond2', 'equals', '55')


def test_text_filters_eq():
    truth = utils.load_csv('test_files/text_filters_eq.csv', 0)
    d = CountFilter('test_files/text_filters.csv')
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
    truth = utils.load_csv('test_files/text_filters_ct.csv', 0)
    d = CountFilter('test_files/text_filters.csv')
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
    truth = utils.load_csv('test_files/text_filters_sw.csv', 0)
    d = CountFilter('test_files/text_filters.csv')
    filt_1 = d.text_filters('name', 'sw', '2R', inplace=False)
    filt_2 = d.text_filters('name', 'Starts With', '2R', inplace=False)
    filt_1.df.sort_index(inplace=True)
    filt_2.df.sort_index(inplace=True)
    truth.sort_index(inplace=True)
    print(filt_1.df)
    assert np.all(filt_1.df == filt_2.df)
    assert np.all(np.squeeze(truth) == np.squeeze(filt_1.df))


def test_text_filters_ew():
    truth = utils.load_csv('test_files/text_filters_ew.csv', 0)
    d = CountFilter('test_files/text_filters.csv')
    filt_1 = d.text_filters('name', 'ew', '3', inplace=False)
    filt_2 = d.text_filters('name', 'ends With', '3', inplace=False)
    filt_1.df.sort_index(inplace=True)
    filt_2.df.sort_index(inplace=True)
    truth.sort_index(inplace=True)
    print(filt_1.df)
    assert np.all(filt_1.df == filt_2.df)
    assert np.all(np.squeeze(truth) == np.squeeze(filt_1.df))


def test_text_filters_invalid_input():
    d = CountFilter(r'test_files/counted.csv')
    with pytest.raises(AssertionError):
        d.text_filters('Cond2', 'contains', '5')
    with pytest.raises(AssertionError):
        d.text_filters('cond2', 'lt', '6')
    with pytest.raises(AssertionError):
        d.text_filters('cond2', 'equals', 55)


def test_count_filter_from_folder():
    truth_all_expr = utils.load_csv(r'test_files/test_count_from_folder_all_expr.csv', 0)
    truth_all_feature = utils.load_csv(r'test_files/test_count_from_folder_all_feature.csv', 0)
    truth_norm = utils.load_csv(r'test_files/test_count_from_folder_norm.csv', 0)
    h_notnorm = CountFilter.from_folder('test_files/test_count_from_folder', norm_to_rpm=False, save_csv=True,
                                        counted_fname='__allexpr_temporary_testfile.csv',
                                        uncounted_fname='__allfeature_temporary_testfile.csv')

    os.remove('test_files/test_count_from_folder/__allexpr_temporary_testfile.csv')
    assert np.all(np.isclose(h_notnorm.df, truth_all_expr))

    h_norm = CountFilter.from_folder('test_files/test_count_from_folder', norm_to_rpm=True, save_csv=False)
    assert np.all(np.isclose(h_norm.df, truth_norm))

    all_feature = utils.load_csv('test_files/test_count_from_folder/__allfeature_temporary_testfile.csv', 0)
    all_feature.sort_index(inplace=True)
    truth_all_feature.sort_index(inplace=True)

    os.remove('test_files/test_count_from_folder/__allfeature_temporary_testfile.csv')
    assert np.all(np.isclose(all_feature, truth_all_feature))

    h_nosuffix = CountFilter.from_folder('test_files/test_count_from_folder', norm_to_rpm=False, save_csv=True,
                                         counted_fname='__allexpr_temporary_testfile',
                                         uncounted_fname='__allfeature_temporary_testfile')
    os.remove('test_files/test_count_from_folder/__allexpr_temporary_testfile.csv')
    os.remove('test_files/test_count_from_folder/__allfeature_temporary_testfile.csv')


def test_biotypes():
    truth = utils.load_csv('test_files/biotypes_truth.csv', 0)
    df = CountFilter(r'test_files/counted_biotype.csv').biotypes(ref=__biotype_ref__)
    truth.sort_index(inplace=True)
    df.sort_index(inplace=True)
    assert np.all(df == truth)


def test_biotypes_long_form():
    assert False


def test_filter_by_row_sum():
    truth = utils.load_csv('test_files/test_filter_row_sum.csv', 0)
    h = CountFilter('test_files/counted.csv')
    h.filter_by_row_sum(29)
    h.df.sort_index(inplace=True)
    truth.sort_index(inplace=True)
    assert np.all(h.df == truth)


def test_sort_inplace():
    c = CountFilter('test_files/counted.csv')
    c.sort(by='cond3', ascending=True, inplace=True)
    assert c.df['cond3'].is_monotonic_increasing


def test_sort_not_inplace():
    c = CountFilter('test_files/counted.csv')
    c_copy = utils.load_csv('test_files/counted.csv', 0)
    c_sorted = c.sort(by='cond3', ascending=True, inplace=False)
    assert c_sorted.df['cond3'].is_monotonic_increasing
    assert np.all(c.df == c_copy)


def test_sort_by_multiple_columns():
    truth = utils.load_csv('test_files/counted_sorted_multiple_truth.csv', 0)
    c = CountFilter('test_files/counted.csv')
    c.sort(by=['cond3', 'cond4', 'cond1', 'cond2'], ascending=[True, False, True, False], inplace=True)
    assert np.all(truth == c.df)


def test_sort_with_na_first():
    truth_first = utils.load_csv('test_files/test_deseq_with_nan_sorted_nanfirst_truth.csv', 0)
    truth_last = utils.load_csv('test_files/test_deseq_with_nan_sorted_nanlast_truth.csv', 0)
    c = CountFilter('test_files/test_deseq_with_nan.csv')
    c.sort(by='padj', ascending=True, na_position='first', inplace=True)
    assert truth_first.equals(c.df)
    c.sort(by='padj', ascending=True, na_position='last', inplace=True)
    assert truth_last.equals(c.df)


def test_sort_descending():
    c = CountFilter('test_files/counted.csv')
    c.sort(by='cond3', ascending=False, inplace=True)
    assert c.df['cond3'].is_monotonic_decreasing


def test_filter_missing_values():
    truth = utils.load_csv('test_files/test_deseq_with_nan_all_removed.csv', 0)
    f = Filter('test_files/test_deseq_with_nan.csv')
    f.filter_missing_values()
    assert np.all(f.df.sort_index() == truth.sort_index())


def test_filter_missing_values_foldchangefilter():
    truth = utils.load_csv('test_files/fc_1_nan_removed.csv', 0, squeeze=True)
    f = FoldChangeFilter('test_files/fc_1_nan.csv', 'num', 'denom')
    res_all = f.filter_missing_values(inplace=False)
    assert truth.equals(res_all.df)
    res_foldchange = f.filter_missing_values('Fold Change', inplace=False)
    assert truth.equals(res_foldchange.df)
    with pytest.raises(AssertionError):
        f.filter_missing_values('column that doesnt exist')


def test_filter_missing_values_one_columns():
    truth = utils.load_csv('test_files/test_deseq_with_nan_basemean_removed.csv', 0)
    f = Filter('test_files/test_deseq_with_nan.csv')
    f.filter_missing_values('baseMean')
    print(f.df.sort_index())
    print(truth.sort_index())
    print(f.df.sort_index() == truth.sort_index())
    assert truth.equals(f.df)


def test_filter_missing_values_multiple_columns():
    truth = utils.load_csv('test_files/test_deseq_with_nan_basemean_pvalue_removed.csv', 0)
    f = Filter('test_files/test_deseq_with_nan.csv')
    f.filter_missing_values(['baseMean', 'pvalue'])
    print(f.df.sort_index())
    print(truth.sort_index())
    print(f.df.sort_index() == truth.sort_index())
    assert truth.equals(f.df)


def test_filter_missing_values_invalid_type():
    f = Filter('test_files/test_deseq_with_nan.csv')
    with pytest.raises(TypeError):
        f.filter_missing_values(columns={'baseMean': True, 'log2FolgChange': False})


def test_filter_missing_values_nonexistent_column():
    f = Filter('test_files/test_deseq_with_nan.csv')
    with pytest.raises(AssertionError):
        f.filter_missing_values('pval')
    with pytest.raises(AssertionError):
        f.filter_missing_values(['padj', 'pval'])


def test_pipeline_api():
    pl = Pipeline()
    pl_count = Pipeline('countfilter')
    pl_deseq = Pipeline(DESeqFilter)
    pl = Pipeline(filter_type='FoldChangeFilter')


def test_pipeline_add_function():
    pl = Pipeline()
    pl.add_function(DESeqFilter.filter_biotype, biotype='protein_coding')
    assert len(pl.functions) == 1 and len(pl.params) == 1
    assert pl.functions[0] == DESeqFilter.filter_biotype
    assert pl.params[0] == ((), {'biotype': 'protein_coding'})

    pl = Pipeline()
    pl.add_function('filter_biotype', 'piRNA')
    assert len(pl.functions) == 1 and len(pl.params) == 1
    assert pl.functions[0] == Filter.filter_biotype
    assert pl.params[0] == (('piRNA',), {})

    pl_deseq = Pipeline('DEseqFilter')
    pl_deseq.add_function(Filter.number_filters, 'log2FoldChange', operator='>', value=5)
    assert len(pl_deseq.functions) == 1 and len(pl_deseq.params) == 1
    assert pl_deseq.functions[0] == DESeqFilter.number_filters
    assert pl_deseq.params[0] == (('log2FoldChange',), {'operator': '>', 'value': 5})


def test_pipeline_add_multiple_functions():
    pl_deseq = Pipeline('DEseqFilter')
    pl_deseq.add_function(Filter.number_filters, 'log2FoldChange', operator='>', value=5)
    pl_deseq.add_function(DESeqFilter.filter_significant)
    pl_deseq.add_function('sort', by='log2FoldChange')

    assert len(pl_deseq.functions) == 3 and len(pl_deseq.params) == 3
    assert pl_deseq.functions == [DESeqFilter.number_filters, DESeqFilter.filter_significant, DESeqFilter.sort]
    assert pl_deseq.params == [(('log2FoldChange',), {'operator': '>', 'value': 5}), ((), {}),
                               ((), {'by': 'log2FoldChange'})]


def test_pipeline_remove_last_function():
    pl = Pipeline()
    pl.add_function(DESeqFilter.filter_biotype, biotype='protein_coding',
                    ref=__biotype_ref__)
    pl.remove_last_function()
    assert len(pl.functions) == 0 and len(pl.params) == 0


def test_pipeline_remove_last_from_empty_pipeline():
    pl = Pipeline()
    with pytest.raises(AssertionError):
        pl.remove_last_function()


def test_pipeline_apply_empty_pipeline():
    pl = Pipeline()
    d = DESeqFilter('test_files/test_deseq.csv')
    with pytest.raises(AssertionError):
        pl.apply_to(d)


def test_pipeline_apply_to():
    pl = Pipeline('deseqfilter')
    pl.add_function('filter_significant', 10 ** -70, opposite=True)
    deseq = DESeqFilter('test_files/test_deseq.csv')
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
    pl2.add_function(Filter.filter_biotype, biotype='protein_coding',
                     ref=__biotype_ref__)
    cnt = CountFilter('test_files/counted.csv')
    cnt_truth = cnt.__copy__()
    cnt_truth.filter_biotype('protein_coding', ref=__biotype_ref__)
    cnt_pipelined = pl2.apply_to(cnt, inplace=False)
    pl2.apply_to(cnt, inplace=True)
    cnt.sort(cnt.columns[0])
    cnt_truth.sort(cnt.columns[0])
    cnt_pipelined.sort(cnt.columns[0])
    assert np.all(cnt.df == cnt_truth.df)
    assert np.all(cnt_pipelined.df == cnt_truth.df)


def test_pipeline_apply_to_with_multiple_functions():
    d = DESeqFilter('test_files/test_deseq_with_nan.csv')
    d_copy = d.__copy__()
    p = Pipeline('deseqfilter')
    p.add_function('filter_missing_values')
    p.add_function(Filter.filter_biotype, biotype='protein_coding', ref=__biotype_ref__)
    p.add_function('number_filters', 'log2FoldChange', 'gt', 0.75)
    p.add_function('sort', 'baseMean', ascending=False)
    p.add_function('biotypes')
    p.remove_last_function()

    d_pipelined = p.apply_to(d_copy, inplace=False)
    p.apply_to(d_copy)
    d.filter_missing_values()
    d.filter_biotype('protein_coding', __biotype_ref__)
    d.number_filters('log2FoldChange', 'gt', 0.75)
    d.sort('baseMean', ascending=False)
    assert d.df.equals(d_pipelined.df)
    assert d.df.equals(d_copy.df)


def test_pipeline_apply_to_invalid_object():
    pl = Pipeline('deseqfilter')
    pl.add_function(DESeqFilter.filter_significant, alpha=10 ** -70)
    cnt = utils.load_csv('test_files/counted.csv', 0)
    with pytest.raises(AssertionError):
        pl.apply_to(cnt)


def test_pipeline_init_invalid_filter_type():
    with pytest.raises(AssertionError):
        pl = Pipeline(filter_type='otherFilter')

    class otherFilter:
        def __init__(self):
            self.value = 'value'
            self.othervalue = 'othervalue'

    with pytest.raises(AssertionError):
        pl = Pipeline(filter_type=otherFilter)
    with pytest.raises(AssertionError):
        pl = Pipeline(filter_type=max)
    with pytest.raises(AssertionError):
        pl = Pipeline(filter_type=5)


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
    pl_deseq.add_function(CountFilter.filter_biotype, biotype='protein_coding',
                          ref=__biotype_ref__)
    with pytest.raises(AssertionError):
        pl_deseq.add_function(CountFilter.filter_low_reads, threshold=5)


def _get_pipeline_with_plot(inplace: bool):
    p = Pipeline('DESeqFilter')
    p.add_function('filter_missing_values')
    p.add_function(Filter.filter_top_n, 'baseMean', n=15)
    p.add_function(DESeqFilter.volcano_plot, alpha=0.001)
    p.add_function('biotypes', ref=__biotype_ref__)
    p.add_function('filter_top_n', 'log2FoldChange', 6)
    d = DESeqFilter('test_files/test_deseq_with_nan.csv')
    res_truth = {}
    d_copy = d.__copy__()
    d.filter_missing_values()
    d.filter_top_n('baseMean', n=15)
    res_truth['volcano_plot_1'] = d.volcano_plot(alpha=0.001)
    res_truth['biotypes_1'] = d.biotypes(ref=__biotype_ref__)
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
    assert res['biotypes_1'].equals(res_truth['biotypes_1'])
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
    d = DESeqFilter('test_files/test_deseq_with_nan.csv')
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
    pl_c.add_function(CountFilter.split_hdbscan, min_cluster_size=3, return_prob=True)
    pl_c.add_function(CountFilter.filter_top_n, by='cond1', n=5)
    c = CountFilter('test_files/counted.csv')
    c_pipeline_res, c_pipeline_dict = pl_c.apply_to(c, inplace=False)
    c_res = c.filter_top_n(by='cond2', n=2, opposite=True, inplace=False)
    c_res, prob = c_res.split_hdbscan(min_cluster_size=3, return_prob=True)
    for i in c_res:
        i.filter_top_n(by='cond1', n=5)
    for i, j in zip(c_res, c_pipeline_res):
        assert i == j
    assert np.all(c_pipeline_dict['split_hdbscan_1'] == prob)

    pl_c.remove_last_function()
    pl_c.remove_last_function()
    pl_c.add_function('split_kmeans', k=[2, 3, 4], random_state=42)
    pl_c.add_function(CountFilter.filter_top_n, by='cond1', n=5)
    c = CountFilter('test_files/counted.csv')
    c_pipeline_res = pl_c.apply_to(c, inplace=False)
    c_res = c.filter_top_n(by='cond2', n=2, opposite=True, inplace=False)
    c_res = c_res.split_kmeans(k=[2, 3, 4], random_state=42)
    c_res_cont = []
    for i in c_res:
        c_res_cont.extend(i)
    for i in c_res_cont:
        i.filter_top_n(by='cond1', n=5)
    for i, j in zip(c_res_cont, c_pipeline_res):
        assert i == j


def test_pipeline_apply_to_with_split_function_inplace_raise_error():
    pl = Pipeline('DESeqFilter')
    pl.add_function('filter_missing_values')
    pl.add_function(DESeqFilter.filter_significant, 0.05, opposite=True)
    pl.add_function(DESeqFilter.split_fold_change_direction)
    with pytest.raises(AssertionError):
        pl.apply_to(DESeqFilter('test_files/test_deseq.csv'), inplace=True)


def test_pipeline_apply_to_multiple_splits():
    pl_c = Pipeline('CountFilter')
    pl_c.add_function(CountFilter.filter_top_n, by='cond2', n=2, opposite=True)
    pl_c.add_function(CountFilter.split_hdbscan, min_cluster_size=3, return_prob=True)
    pl_c.add_function(CountFilter.split_kmedoids, k=2, random_state=42)
    pl_c.add_function(CountFilter.split_by_reads, 15)

    c = CountFilter('test_files/counted.csv')
    c_pipeline_res, c_pipeline_dict = pl_c.apply_to(c, inplace=False)
    c_res = c.filter_top_n(by='cond2', n=2, opposite=True, inplace=False)
    c_res, prob = c_res.split_hdbscan(min_cluster_size=3, return_prob=True)
    c_res_cont = []
    for i in c_res:
        c_res_cont.extend(i.split_kmedoids(k=2, random_state=42))
    c_res_cont_fin = []
    for i in c_res_cont:
        c_res_cont_fin.extend(i.split_by_reads(15))
    assert len(c_pipeline_res) == len(c_res_cont_fin)
    for i, j in zip(c_res_cont_fin, c_pipeline_res):
        print(i.df)
        print(j.df)
        assert i == j
    assert np.all(c_pipeline_dict['split_hdbscan_1'] == prob)


def test_pipeline_apply_to_filter_normalize_split_plot():
    scaling_factor_path = 'test_files/big_scaling_factors.csv'
    pl_c = Pipeline('CountFilter')
    pl_c.add_function(CountFilter.normalize_with_scaling_factors, scaling_factor_path)
    pl_c.add_function('biotypes', ref=__biotype_ref__)
    pl_c.add_function(CountFilter.filter_top_n, by='cond3rep1', n=270, opposite=True)
    pl_c.add_function(CountFilter.split_hdbscan, min_cluster_size=100, return_prob=True)
    pl_c.add_function(CountFilter.filter_low_reads, threshold=10)
    pl_c.add_function(CountFilter.clustergram)
    pl_c.add_function(CountFilter.split_kmedoids, k=[2, 3, 7], random_state=42, n_init=1)
    pl_c.add_function(CountFilter.sort, by='cond2rep3')
    pl_c.add_function(CountFilter.biotypes, 'long', __biotype_ref__)

    c = CountFilter('test_files/big_counted.csv')
    c_pipeline_res, c_pipeline_dict = pl_c.apply_to(c, inplace=False)
    c_dict = dict()
    c.normalize_with_scaling_factors(scaling_factor_path)
    c_dict['biotypes_1'] = c.biotypes(ref=__biotype_ref__)
    c_res = c.filter_top_n(by='cond3rep1', n=270, opposite=True, inplace=False)
    c_res, prob = c_res.split_hdbscan(min_cluster_size=100, return_prob=True)
    c_dict['split_hdbscan_1'] = prob
    for obj in c_res:
        obj.filter_low_reads(threshold=10)
    clustergrams = []
    for obj in c_res:
        clustergrams.append(obj.clustergram())
    c_dict['clustergram_1'] = tuple(clustergrams)
    c_res_cont = []
    for obj in c_res:
        k_out = obj.split_kmedoids(k=[2, 3, 7], random_state=42, n_init=1)
        for k in k_out:
            c_res_cont.extend(k)
    for obj in c_res_cont:
        obj.sort(by='cond2rep3')
    biotypes = []
    for obj in c_res_cont:
        biotypes.append(obj.biotypes('long', __biotype_ref__))
    c_dict['biotypes_2'] = tuple(biotypes)

    assert len(c_res_cont) == len(c_pipeline_res)
    for i, j in zip(c_res_cont, c_pipeline_res):
        assert i == j
    assert np.all(c_dict.keys() == c_pipeline_dict.keys())
    assert c_dict['biotypes_1'].equals(c_pipeline_dict['biotypes_1'])
    assert len(c_dict['clustergram_1']) == len(c_pipeline_dict['clustergram_1'])
    for i, j in zip(c_dict['clustergram_1'], c_pipeline_dict['clustergram_1']):
        assert type(i) == type(j)
    assert len(c_dict['biotypes_2']) == len(c_pipeline_dict['biotypes_2'])
    for i, j in zip(c_dict['biotypes_2'], c_pipeline_dict['biotypes_2']):
        assert i.equals(j)


def test_split_kmeans_api():
    c = CountFilter('test_files/big_counted.csv')
    res = c.split_kmeans(k=4)
    assert isinstance(res, tuple)
    assert len(res) == 4
    _test_correct_clustering_split(c, res)
    res2 = c.split_kmeans(k=[2, 3, 7], random_state=42, n_init=1, max_iter=10, plot_style='std_bar')
    assert isinstance(res2, list)
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


def test_split_hdbscan_api():
    c = CountFilter('test_files/big_counted.csv')
    res = c.split_hdbscan(100)
    _test_correct_clustering_split(c, res, True)
    res2 = c.split_hdbscan(4, 5, 'manhattan', 0.2, 'leaf', plot_style='std_area', return_prob=True)
    assert isinstance(res2, list)
    assert isinstance(res2[0], tuple)
    assert isinstance(res2[1], np.ndarray)
    with pytest.raises(AssertionError):
        c.split_hdbscan(min_cluster_size=1)
    with pytest.raises(AssertionError):
        c.split_hdbscan(c.shape[0] + 1)


def test_split_kmedoids():
    c = CountFilter('test_files/big_counted.csv')
    c.filter_top_n(by='cond1rep1', n=2000)
    res = c.split_kmedoids(k=4, n_init=3)
    assert isinstance(res, tuple)
    assert len(res) == 4
    _test_correct_clustering_split(c, res)
    res2 = c.split_kmedoids(k=[2, 3, 7], random_state=42, n_init=1, max_iter=10, plot_style='std_bar')
    assert isinstance(res2, list)
    assert np.all([isinstance(i, tuple) for i in res2])
    [_test_correct_clustering_split(c, i) for i in res2]


def test_gap_statistic():
    assert False


def test_silhouette_method():
    assert False


def test_parse_k(monkeypatch):
    c = CountFilter('test_files/counted.csv')
    monkeypatch.setattr(CountFilter, "_silhouette",
                        lambda self, clusterer_class, random_state, n_init, max_iter, max_clusters: ('success', None))
    assert c._parse_k('silhouette', KMeans, 0, 0, 0, 0) == ['success']

    monkeypatch.setattr(CountFilter, "_gap_statistic",
                        lambda self, clusterer_class, random_state, n_init, max_iter, max_clusters: ('success', None))
    assert c._parse_k('gap', KMeans, 0, 0, 0, 0) == ['success']

    assert list(c._parse_k(10, KMeans, 0, 0, 0, 0)) == [10]
    assert list(c._parse_k([7, 2, 5], KMeans, 0, 0, 0, 0)) == [7, 2, 5]
    assert list(c._parse_k(range(2, 9), KMeans, 0, 0, 0, 0)) == list(range(2, 9))
    with pytest.raises(AssertionError):
        print(list(c._parse_k([5, 2, '3'], KMeans, 0, 0, 0, 0)))
    with pytest.raises(AssertionError):
        c._parse_k('a string', KMeans, 0, 0, 0, 0)
    with pytest.raises(AssertionError):
        c._parse_k(1, KMeans, 0, 0, 0, 0)
    with pytest.raises(AssertionError):
        c._parse_k([3, 5, 1], KMeans, 0, 0, 0, 0)


def test_fc_randomization():
    truth = utils.load_csv('test_files/fc_randomization_truth.csv')
    fc1 = FoldChangeFilter("test_files/fc_1.csv", 'a', 'b')
    fc2 = FoldChangeFilter("test_files/fc_2.csv", "c", "d")
    res = fc1.randomization_test(fc2, random_seed=0)
    assert np.all(truth['significant'] == res['significant'])
    assert np.isclose(truth.iloc[:, :-1], res.iloc[:, :-1]).all()


def test_filter_save_csv():
    assert False


def test_kmedoidsiter_api():
    truth = KMedoids(3, max_iter=300, init='k-medoids++', random_state=42)
    kmeds = _KMedoidsIter(3, max_iter=300, init='k-medoids++', n_init=1, random_state=42)
    c = CountFilter('test_files/counted.csv')
    truth.fit(c.df)
    kmeds.fit(c.df)
    assert np.all(truth.cluster_centers_ == kmeds.cluster_centers_)
    assert np.all(truth.inertia_ == kmeds.inertia_)

    assert np.all(truth.predict(c.df) == kmeds.predict(c.df))
    assert np.all(truth.fit_predict(c.df) == kmeds.fit_predict(c.df))

    kmeds_rand = _KMedoidsIter(3, max_iter=300, init='k-medoids++', n_init=3)
    kmeds_rand.fit(c.df)
    kmeds_rand.predict(c.df)
    kmeds_rand.fit_predict(c.df)


def test_kmedoidsiter_iter():
    kmeds = _KMedoidsIter(3, max_iter=300, init='k-medoids++', n_init=5, random_state=0)
    c = CountFilter('test_files/counted.csv')
    kmeds.fit(c.df)

    inertias = []
    clusterers = []
    for i in range(5):
        clusterers.append(KMedoids(3, max_iter=300, init='k-medoids++', random_state=0).fit(c.df))
        inertias.append(clusterers[i].inertia_)
    truth_inertia = max(inertias)
    truth_kmeds = clusterers[np.argmax(inertias)]
    assert kmeds.inertia_ == truth_inertia
    assert np.all(kmeds.clusterer.predict(c.df) == truth_kmeds.predict(c.df))
