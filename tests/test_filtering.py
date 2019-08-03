import pytest
import numpy as np
import pandas as pd
from pathlib import Path
from rnalysis import general
from rnalysis.filtering import *


def test_deseqfilter_api():
    d = DESeqFilter('test_deseq_biotype.csv')


def test_filter_inplace():
    d = DESeqFilter('test_deseq_no_nans.csv')
    d_copy = DESeqFilter('test_deseq_no_nans.csv')
    truth = general.load_csv('all_expr.csv')
    d_inplace_false = d._inplace(truth, opposite=False, inplace=False, suffix='suffix')
    assert np.all(d_inplace_false.df == truth)
    assert np.all(d.df == d_copy.df)
    d._inplace(truth, opposite=False, inplace=True, suffix='other_suffix')
    assert np.all(d.df == truth)


def test_htcountfilter_api():
    h = HTCountFilter('all_expr.csv')


def test_htcountfilter_norm_reads_to_rpm():
    truth = general.load_csv(r"test_norm_reads_rpm.csv", 0)
    h = HTCountFilter(r"all_expr.csv")
    h.norm_reads_to_rpm(r"all_feature.csv")
    assert np.isclose(truth, h.df).all()


def test_htcountfilter_norm_reads_with_size_factor():
    truth = general.load_csv(r"test_norm_size_factor.csv", 0)
    h = HTCountFilter(r"all_expr.csv")
    h.norm_reads_with_size_factor(r"size_factors.csv")
    assert np.isclose(truth, h.df).all()


def test_filter_low_reads():
    truth = general.load_csv("all_expr_low_rpm_truth.csv", 0)
    h = HTCountFilter("all_expr_low_rpm.csv")
    h.filter_low_reads(threshold=5)
    assert np.isclose(truth, h.df).all()


def test_filter_low_reads_reverse():
    h = HTCountFilter(r"all_expr.csv")
    low_truth = general.load_csv(r"all_expr_below60_rpm.csv", 0)
    h.filter_low_reads(threshold=60, opposite=True)
    h.df.sort_index(inplace=True)
    low_truth.sort_index(inplace=True)
    print(h.shape)
    print(low_truth.shape)
    print(h.df)
    print(low_truth)

    assert np.all(h.df == low_truth)


def test_htcount_filter_biotype():
    truth_protein_coding = general.load_csv('all_expr_biotype_protein_coding.csv', 0)
    truth_pirna = general.load_csv('all_expr_biotype_piRNA.csv', 0)
    h = HTCountFilter("all_expr_biotype.csv")
    protein_coding = h.filter_biotype(ref='big_table_for_tests.csv', inplace=False)
    pirna = h.filter_biotype('piRNA', ref='big_table_for_tests.csv', inplace=False)
    pirna.df.sort_index(inplace=True)
    protein_coding.df.sort_index(inplace=True)
    truth_protein_coding.sort_index(inplace=True)
    truth_pirna.sort_index(inplace=True)
    assert np.all(truth_protein_coding == protein_coding.df)
    assert np.all(truth_pirna == pirna.df)


def test_htcount_filter_biotype_opposite():
    truth_no_pirna = general.load_csv(r'all_expr_biotype_no_piRNA.csv', 0)
    h = HTCountFilter("all_expr_biotype.csv")
    h.filter_biotype('piRNA', opposite=True, inplace=True)
    h.df.sort_index(inplace=True)
    truth_no_pirna.sort_index(inplace=True)
    assert np.all(h.df == truth_no_pirna)


def test_filter_by_bigtable_group_union():
    union_truth = general.load_csv(r'all_expr_filter_by_bigtable_union_truth.csv', 0)
    h = HTCountFilter('all_expr_filter_by_bigtable.csv')
    union = h.filter_by_bigtable_group(['attribute1', 'attribute2'], mode='union',
                                       ref='big_table_for_tests.csv', inplace=False)
    union.df.sort_index(inplace=True)
    union_truth.sort_index(inplace=True)
    assert np.all(union.df == union_truth)


def test_filter_by_bigtable_group_intersection():
    intersection_truth = general.load_csv(r'all_expr_filter_by_bigtable_intersect_truth.csv', 0)
    h = HTCountFilter('all_expr_filter_by_bigtable.csv')
    intersection = h.filter_by_bigtable_group(['attribute1', 'attribute2'], mode='intersection',
                                              ref='big_table_for_tests.csv',
                                              inplace=False)
    intersection.df.sort_index(inplace=True)
    intersection_truth.sort_index(inplace=True)
    assert np.all(intersection.df == intersection_truth)


def test_deseq_filter_significant():
    truth = general.load_csv("test_deseq_sig_truth.csv", 0)
    d = DESeqFilter("test_deseq_sig.csv")
    d.filter_significant(alpha=0.05)
    assert np.all(d.df == truth)


def test_deseq_filter_significant_opposite():
    truth = general.load_csv(r'test_deseq_not_sig_truth.csv', 0)
    d = DESeqFilter("test_deseq_sig.csv")
    d.filter_significant(alpha=0.05, opposite=True)
    d.df.sort_index(inplace=True)
    truth.sort_index(inplace=True)
    truth.fillna(1234567890, inplace=True)
    d.df.fillna(1234567890, inplace=True)
    assert np.all(d.df == truth)


def test_deseq_filter_top_n():
    truth = general.load_csv("test_deseq_top10.csv", 0)
    d = DESeqFilter("test_deseq.csv")
    d.filter_top_n(10)
    assert np.isclose(truth, d.df).all()


def test_deseq_filter_abs_fold_change():
    truth = general.load_csv("test_deseq_fc_4_truth.csv", 0)
    d = DESeqFilter("test_deseq_fc.csv")
    fc4 = d.filter_abs_fold_change(4, inplace=False)
    fc4.df.sort_index(inplace=True)
    truth.sort_index(inplace=True)
    assert np.all(fc4.df == truth)


def test_deseq_filter_fold_change_direction():
    pos_truth = general.load_csv("test_deseq_fc_pos_truth.csv", 0)
    neg_truth = general.load_csv("test_deseq_fc_neg_truth.csv", 0)
    d = DESeqFilter("test_deseq_fc.csv")
    pos = d.filter_fold_change_direction('pos', inplace=False)
    neg = d.filter_fold_change_direction('neg', inplace=False)
    assert np.all(pos.df == pos_truth)
    assert np.all(neg.df == neg_truth)


def test_deseq_split_fold_change():
    d = DESeqFilter("test_deseq_fc.csv")
    pos_truth = general.load_csv("test_deseq_fc_pos_truth.csv", 0)
    neg_truth = general.load_csv("test_deseq_fc_neg_truth.csv", 0)
    d = DESeqFilter("test_deseq_fc.csv")
    pos, neg = d.split_fold_change_direction()
    assert np.all(pos.df == pos_truth)
    assert np.all(neg.df == neg_truth)


def test_deseq_intersection():
    intersection_truth = {'WBGene00021375', 'WBGene00044258', 'WBGene00219304', 'WBGene00194708', 'WBGene00018199',
                          'WBGene00019174', 'WBGene00021019', 'WBGene00013816', 'WBGene00045366', 'WBGene00219307',
                          'WBGene00045410', 'WBGene00010100', 'WBGene00077437', 'WBGene00007674', 'WBGene00023036',
                          'WBGene00012648', 'WBGene00022486'}
    set1 = DESeqFilter('test_deseq_set_ops_1.csv')
    set2 = DESeqFilter('test_deseq_set_ops_2.csv')

    assert set1.intersection(set2) == intersection_truth


def test_deseq_union():
    intersection_truth = {'WBGene00021375', 'WBGene00044258', 'WBGene00219304', 'WBGene00194708', 'WBGene00018199',
                          'WBGene00019174', 'WBGene00021019', 'WBGene00013816', 'WBGene00045366', 'WBGene00219307',
                          'WBGene00045410', 'WBGene00010100', 'WBGene00077437', 'WBGene00007674', 'WBGene00023036',
                          'WBGene00012648', 'WBGene00022486'}
    set2_unique = {'WBGene00018193', 'WBGene00021589', 'WBGene00001118', 'WBGene00010755', 'WBGene00020407',
                   'WBGene00044799', 'WBGene00021654', 'WBGene00012919', 'WBGene00021605'}
    set1_unique = {'WBGene00008447', 'WBGene00021018', 'WBGene00012452', 'WBGene00010507', 'WBGene00022730',
                   'WBGene00012961', 'WBGene00022438', 'WBGene00016635', 'WBGene00044478'}

    set1 = DESeqFilter('test_deseq_set_ops_1.csv')
    set2 = DESeqFilter('test_deseq_set_ops_2.csv')
    union_truth = intersection_truth.union(set1_unique.union(set2_unique))
    assert set1.union(set2) == union_truth


def test_deseq_difference():
    set2_unique = {'WBGene00018193', 'WBGene00021589', 'WBGene00001118', 'WBGene00010755', 'WBGene00020407',
                   'WBGene00044799', 'WBGene00021654', 'WBGene00012919', 'WBGene00021605'}
    set1_unique = {'WBGene00008447', 'WBGene00021018', 'WBGene00012452', 'WBGene00010507', 'WBGene00022730',
                   'WBGene00012961', 'WBGene00022438', 'WBGene00016635', 'WBGene00044478'}

    set1 = DESeqFilter('test_deseq_set_ops_1.csv')
    set2 = DESeqFilter('test_deseq_set_ops_2.csv')

    assert set1.difference(set2) == set1_unique
    assert set2.difference(set1) == set2_unique


def test_deseq_symmetric_difference():
    intersection_truth = {'WBGene00021375', 'WBGene00044258', 'WBGene00219304', 'WBGene00194708', 'WBGene00018199',
                          'WBGene00019174', 'WBGene00021019', 'WBGene00013816', 'WBGene00045366', 'WBGene00219307',
                          'WBGene00045410', 'WBGene00010100', 'WBGene00077437', 'WBGene00007674', 'WBGene00023036',
                          'WBGene00012648', 'WBGene00022486'}
    set2_unique = {'WBGene00018193', 'WBGene00021589', 'WBGene00001118', 'WBGene00010755', 'WBGene00020407',
                   'WBGene00044799', 'WBGene00021654', 'WBGene00012919', 'WBGene00021605'}
    set1_unique = {'WBGene00008447', 'WBGene00021018', 'WBGene00012452', 'WBGene00010507', 'WBGene00022730',
                   'WBGene00012961', 'WBGene00022438', 'WBGene00016635', 'WBGene00044478'}

    set1 = DESeqFilter('test_deseq_set_ops_1.csv')
    set2 = DESeqFilter('test_deseq_set_ops_2.csv')

    assert set1.symmetric_difference(set2) == set2.symmetric_difference(set1)
    assert set1.symmetric_difference(set2) == set1_unique.union(set2_unique)


def test_deseq_feature_set():
    truth = {'WBGene00008447', 'WBGene00021018', 'WBGene00012452', 'WBGene00010507', 'WBGene00022730', 'WBGene00012648',
             'WBGene00012961', 'WBGene00022438', 'WBGene00016635', 'WBGene00044478', 'WBGene00021375',
             'WBGene00044258', 'WBGene00219304', 'WBGene00194708', 'WBGene00018199', 'WBGene00022486',
             'WBGene00019174', 'WBGene00021019', 'WBGene00013816', 'WBGene00045366', 'WBGene00219307',
             'WBGene00045410', 'WBGene00010100', 'WBGene00077437', 'WBGene00007674', 'WBGene00023036'}
    d = DESeqFilter('test_deseq_set_ops_1.csv')
    assert d.features_set() == truth


def test_deseq_feature_string():
    truth = {'WBGene00008447', 'WBGene00021018', 'WBGene00012452', 'WBGene00010507', 'WBGene00022730', 'WBGene00012648',
             'WBGene00012961', 'WBGene00022438', 'WBGene00016635', 'WBGene00044478', 'WBGene00021375',
             'WBGene00044258', 'WBGene00219304', 'WBGene00194708', 'WBGene00018199', 'WBGene00022486',
             'WBGene00019174', 'WBGene00021019', 'WBGene00013816', 'WBGene00045366', 'WBGene00219307',
             'WBGene00045410', 'WBGene00010100', 'WBGene00077437', 'WBGene00007674', 'WBGene00023036'}
    d = DESeqFilter('test_deseq_set_ops_1.csv')
    assert set(d.features_string().split("\n")) == truth


def test_htcount_rpm_negative_threshold():
    h = HTCountFilter("all_expr.csv")
    with pytest.raises(AssertionError):
        h.filter_low_reads(threshold=-3)


def test_htcount_threshold_invalid():
    h = HTCountFilter("all_expr.csv")
    with pytest.raises(AssertionError):
        h.filter_low_reads("5")


def test_htcount_split_by_reads():
    h = HTCountFilter(r"all_expr.csv")
    high_truth = general.load_csv(r"all_expr_above60_rpm.csv", 0)
    low_truth = general.load_csv(r"all_expr_below60_rpm.csv", 0)
    high, low = h.split_by_reads(threshold=60)
    assert np.all(high.df == high_truth)
    assert np.all(low.df == low_truth)


def test_filter_percentile():
    truth = general.load_csv(r'test_deseq_percentile_0.25.csv', 0)
    h = DESeqFilter(r'test_deseq_percentile.csv')
    h.filter_percentile(0.25, 'padj', inplace=True)
    h.df.sort_index(inplace=True)
    truth.sort_index(inplace=True)
    assert np.all(h.df == truth)


def test_split_by_percentile():
    truth_below = general.load_csv(r'test_deseq_percentile_0.25.csv', 0)
    truth_above = general.load_csv(r'test_deseq_percentile_0.75.csv', 0)
    h = DESeqFilter(r'test_deseq_percentile.csv')
    below, above = h.split_by_percentile(0.25, 'padj')
    for i in [truth_below, truth_above, below.df, above.df]:
        i.sort_index(inplace=True)
    assert np.all(truth_below == below.df)
    assert np.all(truth_above == above.df)


def test_htcount_filter_biotype_multiple():
    truth = general.load_csv('all_expr_biotype_piRNA_protein_coding.csv', 0)
    h = HTCountFilter("all_expr_biotype.csv")
    both = h.filter_biotype(['protein_coding', 'piRNA'], ref='big_table_for_tests.csv', inplace=False)
    both.df.sort_index(inplace=True)
    truth.sort_index(inplace=True)
    assert np.all(truth == both.df)


def test_htcount_filter_biotype_multiple_opposite():
    truth = general.load_csv('all_expr_biotype_piRNA_protein_coding_opposite.csv', 0)
    h = HTCountFilter("all_expr_biotype.csv")
    neither = h.filter_biotype(['protein_coding', 'piRNA'], ref='big_table_for_tests.csv', inplace=False, opposite=True)
    neither.df.sort_index(inplace=True)
    truth.sort_index(inplace=True)
    assert np.all(truth == neither.df)


def test_deseq_filter_biotype():
    truth_protein_coding = general.load_csv('test_deseq_biotype_protein_coding.csv', 0)
    truth_pirna = general.load_csv('test_deseq_biotype_piRNA.csv', 0)
    d = DESeqFilter("test_deseq_biotype.csv")
    protein_coding = d.filter_biotype(ref='big_table_for_tests.csv', inplace=False)
    pirna = d.filter_biotype('piRNA', ref='big_table_for_tests.csv', inplace=False)
    pirna.df.sort_index(inplace=True)
    protein_coding.df.sort_index(inplace=True)
    truth_protein_coding.sort_index(inplace=True)
    truth_pirna.sort_index(inplace=True)
    assert np.all(truth_protein_coding == protein_coding.df)
    assert np.all(truth_pirna == pirna.df)


def test_deseq_filter_biotype_opposite():
    truth_no_pirna = general.load_csv(r'test_deseq_biotype_piRNA_opposite.csv', 0)
    d = DESeqFilter("test_deseq_biotype.csv")
    d.filter_biotype('piRNA', opposite=True, inplace=True)
    d.df.sort_index(inplace=True)
    truth_no_pirna.sort_index(inplace=True)
    assert np.all(d.df == truth_no_pirna)


def test_deseq_filter_biotype_multiple():
    truth = general.load_csv('test_deseq_biotype_piRNA_protein_coding.csv', 0)
    d = DESeqFilter("test_deseq_biotype.csv")
    both = d.filter_biotype(['protein_coding', 'piRNA'], ref='big_table_for_tests.csv', inplace=False)
    both.df.sort_index(inplace=True)
    truth.sort_index(inplace=True)
    assert np.all(truth == both.df)


def test_deseq_filter_biotype_multiple_opposite():
    truth = general.load_csv('test_deseq_biotype_piRNA_protein_coding_opposite.csv', 0)
    d = DESeqFilter("test_deseq_biotype.csv")
    neither = d.filter_biotype(['protein_coding', 'piRNA'], ref='big_table_for_tests.csv', inplace=False, opposite=True)
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

    set1 = DESeqFilter('test_deseq_set_ops_1.csv')
    set2 = DESeqFilter('test_deseq_set_ops_2.csv')
    set3 = {'WBGene00077437', 'WBGene00007674', 'WBGene00023036', 'WBGene00012648', 'WBGene44444444', 'WBGene99999999',
            'WBGene98765432'}
    union_truth = intersection_truth.union(set1_unique, set2_unique, set3_unique)
    assert set1.union(set2, set3) == union_truth


def test_deseqfilter_intersection_multiple():
    intersection_truth = {'WBGene00077437', 'WBGene00007674', 'WBGene00023036',
                          'WBGene00012648', 'WBGene00022486'}
    set1 = DESeqFilter('test_deseq_set_ops_1.csv')
    set2 = DESeqFilter('test_deseq_set_ops_2.csv')
    set3 = {'WBGene00077437', 'WBGene00007674', 'WBGene00023036', 'WBGene00012648', 'WBGene00022486', 'WBGene99999999',
            'WBGene98765432'}

    assert set1.intersection(set2, set3) == intersection_truth


def test_deseqfilter_difference_multiple():
    set2_unique = {'WBGene00021589', 'WBGene00001118', 'WBGene00010755', 'WBGene00020407',
                   'WBGene00044799', 'WBGene00021654', 'WBGene00012919', 'WBGene00021605'}
    set1_unique = {'WBGene00021018', 'WBGene00012452', 'WBGene00010507', 'WBGene00022730',
                   'WBGene00012961', 'WBGene00022438', 'WBGene00016635', 'WBGene00044478'}

    set1 = DESeqFilter('test_deseq_set_ops_1.csv')
    set2 = DESeqFilter('test_deseq_set_ops_2.csv')
    set3 = {'WBGene00018193', 'WBGene00008447', 'WBGene12345678'}

    assert set1.difference(set2, set3) == set1_unique
    assert set2.difference(set3, set1) == set2_unique
