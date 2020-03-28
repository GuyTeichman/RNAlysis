import pytest
import numpy as np
import pandas as pd
from pathlib import Path
from rnalysis import general
from rnalysis.filtering import *
import os

def test_deseqfilter_api():
    d = DESeqFilter('test_deseq_biotype.csv')


def test_filter_inplace():
    d = DESeqFilter('test_deseq_no_nans.csv')
    d_copy = DESeqFilter('test_deseq_no_nans.csv')
    truth = general.load_csv('counted.csv')
    d_inplace_false = d._inplace(truth, opposite=False, inplace=False, suffix='suffix')
    assert np.all(d_inplace_false.df == truth)
    assert np.all(d.df == d_copy.df)
    d._inplace(truth, opposite=False, inplace=True, suffix='other_suffix')
    assert np.all(d.df == truth)


def test_countfilter_api():
    h = CountFilter('counted.csv')


def test_countfilter_normalize_to_rpm():
    truth = general.load_csv(r"test_norm_reads_rpm.csv", 0)
    h = CountFilter(r"counted.csv")
    h.normalize_to_rpm(r"uncounted.csv")
    assert np.isclose(truth, h.df).all()


def test_countfilter_norm_reads_with_scaling_factors():
    truth = general.load_csv(r"test_norm_scaling_factors.csv", 0)
    h = CountFilter(r"counted.csv")
    h.normalize_with_scaling_factors(r"scaling_factors.csv")
    assert np.isclose(truth, h.df).all()


def test_filter_low_reads():
    truth = general.load_csv("counted_low_rpm_truth.csv", 0)
    h = CountFilter("counted_low_rpm.csv")
    h.filter_low_reads(threshold=5)
    assert np.isclose(truth, h.df).all()


def test_filter_low_reads_reverse():
    h = CountFilter(r"counted.csv")
    low_truth = general.load_csv(r"counted_below60_rpm.csv", 0)
    h.filter_low_reads(threshold=60, opposite=True)
    h.df.sort_index(inplace=True)
    low_truth.sort_index(inplace=True)
    print(h.shape)
    print(low_truth.shape)
    print(h.df)
    print(low_truth)

    assert np.all(h.df == low_truth)


def test_htcount_filter_biotype():
    truth_protein_coding = general.load_csv('counted_biotype_protein_coding.csv', 0)
    truth_pirna = general.load_csv('counted_biotype_piRNA.csv', 0)
    h = CountFilter("counted_biotype.csv")
    protein_coding = h.filter_biotype(ref='biotype_ref_table_for_tests.csv', inplace=False)
    pirna = h.filter_biotype('piRNA', ref='biotype_ref_table_for_tests.csv', inplace=False)
    pirna.df.sort_index(inplace=True)
    protein_coding.df.sort_index(inplace=True)
    truth_protein_coding.sort_index(inplace=True)
    truth_pirna.sort_index(inplace=True)
    assert np.all(truth_protein_coding == protein_coding.df)
    assert np.all(truth_pirna == pirna.df)


def test_htcount_filter_biotype_opposite():
    truth_no_pirna = general.load_csv(r'counted_biotype_no_piRNA.csv', 0)
    h = CountFilter("counted_biotype.csv")
    h.filter_biotype('piRNA', ref='biotype_ref_table_for_tests.csv', opposite=True, inplace=True)
    h.df.sort_index(inplace=True)
    truth_no_pirna.sort_index(inplace=True)
    assert np.all(h.df == truth_no_pirna)


def test_filter_by_attribute_union():
    union_truth = general.load_csv(r'counted_filter_by_bigtable_union_truth.csv', 0)
    h = CountFilter('counted_filter_by_bigtable.csv')
    union = h.filter_by_attribute(['attribute1', 'attribute2'], mode='union',
                                       ref='attr_ref_table_for_tests.csv', inplace=False)
    union.df.sort_index(inplace=True)
    union_truth.sort_index(inplace=True)
    assert np.all(union.df == union_truth)


def test_filter_by_attribute_intersection():
    intersection_truth = general.load_csv(r'counted_filter_by_bigtable_intersect_truth.csv', 0)
    h = CountFilter('counted_filter_by_bigtable.csv')
    intersection = h.filter_by_attribute(['attribute1', 'attribute2'], mode='intersection',
                                              ref='attr_ref_table_for_tests.csv',
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


def test_filter_top_n_ascending_number():
    truth = general.load_csv("test_deseq_top10.csv", 0)
    d = DESeqFilter("test_deseq.csv")
    d.filter_top_n('padj', 10)
    d.df.sort_index(inplace=True)
    truth.sort_index(inplace=True)
    assert np.isclose(truth, d.df).all()


def test_filter_top_n_ascending_text():
    truth = general.load_csv("test_deseq_top10_text_ascend.csv", 0)
    d = DESeqFilter("test_deseq_textcol.csv")
    d.filter_top_n('textcol', 10, True)
    d.df.sort_index(inplace=True)
    truth.sort_index(inplace=True)
    assert np.all(truth == d.df)


def test_filter_top_n_descending_number():
    truth = general.load_csv("test_deseq_bottom7.csv", 0)
    d = DESeqFilter("test_deseq.csv")
    d.filter_top_n('log2FoldChange', 7, False)
    d.df.sort_index(inplace=True)
    truth.sort_index(inplace=True)
    assert np.isclose(truth, d.df).all()


def test_filter_top_n_descending_text():
    truth = general.load_csv("test_deseq_bottom10_text_descend.csv", 0)
    d = DESeqFilter("test_deseq_textcol.csv")
    d.filter_top_n('textcol', 10, False)
    d.df.sort_index(inplace=True)
    truth.sort_index(inplace=True)
    assert np.all(truth == d.df)


def test_filter_top_n_nonexisting_column():
    d = DESeqFilter("test_deseq.csv")
    colname = 'somecol'
    with pytest.raises(AssertionError):
        d.filter_top_n(colname, 5)
        d.filter_top_n([d.df.columns[0], colname])
    assert colname not in d.df.columns


def test_deseq_filter_abs_log2_fold_change():
    truth = general.load_csv("test_deseq_fc_4_truth.csv", 0)
    d = DESeqFilter("test_deseq_fc.csv")
    fc4 = d.filter_abs_log2_fold_change(4, inplace=False)
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


def test_intersection():
    intersection_truth = {'WBGene00021375', 'WBGene00044258', 'WBGene00219304', 'WBGene00194708', 'WBGene00018199',
                          'WBGene00019174', 'WBGene00021019', 'WBGene00013816', 'WBGene00045366', 'WBGene00219307',
                          'WBGene00045410', 'WBGene00010100', 'WBGene00077437', 'WBGene00007674', 'WBGene00023036',
                          'WBGene00012648', 'WBGene00022486'}
    set1 = DESeqFilter('test_deseq_set_ops_1.csv')
    set2 = DESeqFilter('test_deseq_set_ops_2.csv')

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

    set1 = DESeqFilter('test_deseq_set_ops_1.csv')
    set2 = DESeqFilter('test_deseq_set_ops_2.csv')
    union_truth = intersection_truth.union(set1_unique.union(set2_unique))
    assert set1.union(set2) == union_truth


def test_difference():
    set2_unique = {'WBGene00018193', 'WBGene00021589', 'WBGene00001118', 'WBGene00010755', 'WBGene00020407',
                   'WBGene00044799', 'WBGene00021654', 'WBGene00012919', 'WBGene00021605'}
    set1_unique = {'WBGene00008447', 'WBGene00021018', 'WBGene00012452', 'WBGene00010507', 'WBGene00022730',
                   'WBGene00012961', 'WBGene00022438', 'WBGene00016635', 'WBGene00044478'}

    set1 = DESeqFilter('test_deseq_set_ops_1.csv')
    set2 = DESeqFilter('test_deseq_set_ops_2.csv')

    assert set1.difference(set2, inplace=False) == set1_unique
    assert set2.difference(set1, inplace=False) == set2_unique


def test_symmetric_difference():
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
    assert d.index_set == truth


def test_deseq_feature_string():
    truth = {'WBGene00008447', 'WBGene00021018', 'WBGene00012452', 'WBGene00010507', 'WBGene00022730', 'WBGene00012648',
             'WBGene00012961', 'WBGene00022438', 'WBGene00016635', 'WBGene00044478', 'WBGene00021375',
             'WBGene00044258', 'WBGene00219304', 'WBGene00194708', 'WBGene00018199', 'WBGene00022486',
             'WBGene00019174', 'WBGene00021019', 'WBGene00013816', 'WBGene00045366', 'WBGene00219307',
             'WBGene00045410', 'WBGene00010100', 'WBGene00077437', 'WBGene00007674', 'WBGene00023036'}
    d = DESeqFilter('test_deseq_set_ops_1.csv')
    assert set(d.index_string.split("\n")) == truth


def test_set_ops_multiple_variable_types():
    set2_unique = {'WBGene00018193', 'WBGene00021589', 'WBGene00001118', 'WBGene00010755', 'WBGene00020407',
                   'WBGene00044799', 'WBGene00021654', 'WBGene00012919', 'WBGene00021605'}
    set1_unique = {'WBGene00008447', 'WBGene00021018', 'WBGene00012452', 'WBGene00010507', 'WBGene00022730',
                   'WBGene00012961', 'WBGene00022438', 'WBGene00016635', 'WBGene00044478'}

    set1 = CountFilter('test_deseq_set_ops_1.csv')
    set2 = DESeqFilter('test_deseq_set_ops_2.csv')

    assert set1.symmetric_difference(set2) == set2.symmetric_difference(set1)
    assert set1.symmetric_difference(set2) == set1_unique.union(set2_unique)


def test_htcount_rpm_negative_threshold():
    h = CountFilter("counted.csv")
    with pytest.raises(AssertionError):
        h.filter_low_reads(threshold=-3)


def test_htcount_threshold_invalid():
    h = CountFilter("counted.csv")
    with pytest.raises(AssertionError):
        h.filter_low_reads("5")


def test_htcount_split_by_reads():
    h = CountFilter(r"counted.csv")
    high_truth = general.load_csv(r"counted_above60_rpm.csv", 0)
    low_truth = general.load_csv(r"counted_below60_rpm.csv", 0)
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
    truth = general.load_csv('counted_biotype_piRNA_protein_coding.csv', 0)
    h = CountFilter("counted_biotype.csv")
    both = h.filter_biotype(['protein_coding', 'piRNA'], ref='biotype_ref_table_for_tests.csv', inplace=False)
    both.df.sort_index(inplace=True)
    truth.sort_index(inplace=True)
    assert np.all(truth == both.df)


def test_htcount_filter_biotype_multiple_opposite():
    truth = general.load_csv('counted_biotype_piRNA_protein_coding_opposite.csv', 0)
    h = CountFilter("counted_biotype.csv")
    neither = h.filter_biotype(['protein_coding', 'piRNA'], ref='biotype_ref_table_for_tests.csv', inplace=False,
                               opposite=True)
    neither.df.sort_index(inplace=True)
    truth.sort_index(inplace=True)
    assert np.all(truth == neither.df)


def test_deseq_filter_biotype():
    truth_protein_coding = general.load_csv('test_deseq_biotype_protein_coding.csv', 0)
    truth_pirna = general.load_csv('test_deseq_biotype_piRNA.csv', 0)
    d = DESeqFilter("test_deseq_biotype.csv")
    protein_coding = d.filter_biotype(ref='biotype_ref_table_for_tests.csv', inplace=False)
    pirna = d.filter_biotype('piRNA', ref='biotype_ref_table_for_tests.csv', inplace=False)
    pirna.df.sort_index(inplace=True)
    protein_coding.df.sort_index(inplace=True)
    truth_protein_coding.sort_index(inplace=True)
    truth_pirna.sort_index(inplace=True)
    assert np.all(truth_protein_coding == protein_coding.df)
    assert np.all(truth_pirna == pirna.df)


def test_deseq_filter_biotype_opposite():
    truth_no_pirna = general.load_csv(r'test_deseq_biotype_piRNA_opposite.csv', 0)
    d = DESeqFilter("test_deseq_biotype.csv")
    d.filter_biotype('piRNA', ref='biotype_ref_table_for_tests.csv', opposite=True, inplace=True)
    d.df.sort_index(inplace=True)
    truth_no_pirna.sort_index(inplace=True)
    assert np.all(d.df == truth_no_pirna)


def test_deseq_filter_biotype_multiple():
    truth = general.load_csv('test_deseq_biotype_piRNA_protein_coding.csv', 0)
    d = DESeqFilter("test_deseq_biotype.csv")
    both = d.filter_biotype(['protein_coding', 'piRNA'], ref='biotype_ref_table_for_tests.csv', inplace=False)
    both.df.sort_index(inplace=True)
    truth.sort_index(inplace=True)
    assert np.all(truth == both.df)


def test_deseq_filter_biotype_multiple_opposite():
    truth = general.load_csv('test_deseq_biotype_piRNA_protein_coding_opposite.csv', 0)
    d = DESeqFilter("test_deseq_biotype.csv")
    neither = d.filter_biotype(['protein_coding', 'piRNA'], ref='biotype_ref_table_for_tests.csv', inplace=False,
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

    assert set1.intersection(set2, set3, inplace=False) == intersection_truth


def test_deseqfilter_difference_multiple():
    set2_unique = {'WBGene00021589', 'WBGene00001118', 'WBGene00010755', 'WBGene00020407',
                   'WBGene00044799', 'WBGene00021654', 'WBGene00012919', 'WBGene00021605'}
    set1_unique = {'WBGene00021018', 'WBGene00012452', 'WBGene00010507', 'WBGene00022730',
                   'WBGene00012961', 'WBGene00022438', 'WBGene00016635', 'WBGene00044478'}

    set1 = DESeqFilter('test_deseq_set_ops_1.csv')
    set2 = DESeqFilter('test_deseq_set_ops_2.csv')
    set3 = {'WBGene00018193', 'WBGene00008447', 'WBGene12345678'}

    assert set1.difference(set2, set3, inplace=False) == set1_unique
    assert set2.difference(set3, set1, inplace=False) == set2_unique


def test_intersection_inplace():
    set1_truth = general.load_csv('test_deseq_set_ops_1_inplace_intersection.csv', 0)
    set2_truth = general.load_csv('test_deseq_set_ops_2_inplace_intersection.csv', 0)
    set1 = DESeqFilter('test_deseq_set_ops_1.csv')
    set2 = DESeqFilter('test_deseq_set_ops_2.csv')
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
    set1_truth = general.load_csv('test_deseq_set_ops_1_inplace_difference.csv', 0)
    set2_truth = general.load_csv('test_deseq_set_ops_2_inplace_difference.csv', 0)
    set1 = DESeqFilter('test_deseq_set_ops_1.csv')
    set2 = DESeqFilter('test_deseq_set_ops_2.csv')
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
    truth = general.load_csv(r'counted_fold_change_truth.csv', 0)
    truth = truth.squeeze()
    h = CountFilter(r'counted_fold_change.csv')
    fc = h.fold_change(['cond1_rep1', 'cond1_rep2'], ['cond2_rep1', 'cond2_rep2'])
    assert truth_num_name == fc.numerator
    assert truth_denom_name == fc.denominator
    assert np.all(np.isclose(fc.df, truth))


def test_fc_randomization():
    truth = general.load_csv('fc_randomization_truth.csv')
    fc1 = FoldChangeFilter("fc_1.csv", 'a', 'b')
    fc2 = FoldChangeFilter("fc_2.csv", "c", "d")
    random_state = np.random.get_state()
    res = fc1.randomization_test(fc2)

    try:
        assert np.all(truth['significant'] == res['significant'])
        assert np.isclose(truth.iloc[:, :-1], res.iloc[:, :-1]).all()
    except AssertionError:
        raise AssertionError(f'Enrichment test failed with the numpy.random state: \n{random_state}')


def test_fcfilter_filter_abs_fc():
    truth = general.load_csv('fcfilter_abs_fold_change_truth.csv', 0)
    truth = truth.squeeze()
    truth.sort_index(inplace=True)
    f = FoldChangeFilter('counted_fold_change_truth.csv', 'numer', 'denom')
    f.filter_abs_log2_fold_change(1)
    f.df.sort_index(inplace=True)
    print(f.df.values)
    print(truth.values)
    assert np.all(np.squeeze(f.df.values) == np.squeeze(truth.values))


def test_number_filters_gt():
    truth = general.load_csv(r'test_deseq_gt.csv', 0)
    d = DESeqFilter(r'test_deseq.csv')
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
    truth = general.load_csv(r'test_deseq_lt.csv', 0)
    d = DESeqFilter(r'test_deseq.csv')
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
    truth = general.load_csv(r'counted_eq.csv', 0)
    d = CountFilter(r'counted.csv')
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
    d = CountFilter(r'counted.csv')
    with pytest.raises(AssertionError):
        d.number_filters('Cond2', 'lt', 5)
    with pytest.raises(AssertionError):
        d.number_filters('cond2', 'contains', 6)
    with pytest.raises(AssertionError):
        d.number_filters('cond2', 'equals', '55')


def test_text_filters_eq():
    truth = general.load_csv(r'text_filters_eq.csv', 0)
    d = CountFilter(r'text_filters.csv')
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
    truth = general.load_csv(r'text_filters_ct.csv', 0)
    d = CountFilter(r'text_filters.csv')
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
    truth = general.load_csv(r'text_filters_sw.csv', 0)
    d = CountFilter(r'text_filters.csv')
    filt_1 = d.text_filters('name', 'sw', '2R', inplace=False)
    filt_2 = d.text_filters('name', 'Starts With', '2R', inplace=False)
    filt_1.df.sort_index(inplace=True)
    filt_2.df.sort_index(inplace=True)
    truth.sort_index(inplace=True)
    print(filt_1.df)
    assert np.all(filt_1.df == filt_2.df)
    assert np.all(np.squeeze(truth) == np.squeeze(filt_1.df))


def test_text_filters_ew():
    truth = general.load_csv(r'text_filters_ew.csv', 0)
    d = CountFilter(r'text_filters.csv')
    filt_1 = d.text_filters('name', 'ew', '3', inplace=False)
    filt_2 = d.text_filters('name', 'ends With', '3', inplace=False)
    filt_1.df.sort_index(inplace=True)
    filt_2.df.sort_index(inplace=True)
    truth.sort_index(inplace=True)
    print(filt_1.df)
    assert np.all(filt_1.df == filt_2.df)
    assert np.all(np.squeeze(truth) == np.squeeze(filt_1.df))


def test_text_filters_invalid_input():
    d = CountFilter(r'counted.csv')
    with pytest.raises(AssertionError):
        d.text_filters('Cond2', 'contains', '5')
    with pytest.raises(AssertionError):
        d.text_filters('cond2', 'lt', '6')
    with pytest.raises(AssertionError):
        d.text_filters('cond2', 'equals', 55)


def test_count_filter_from_folder():
    truth_all_expr = general.load_csv(r'test_count_from_folder_all_expr.csv', 0)
    truth_all_feature = general.load_csv(r'test_count_from_folder_all_feature.csv', 0)
    truth_norm = general.load_csv(r'test_count_from_folder_norm.csv', 0)
    h_notnorm = CountFilter.from_folder('test_count_from_folder', norm_to_rpm=False, save_csv=True,
                                        counted_fname='__allexpr_temporary_testfile.csv',
                                        uncounted_fname='__allfeature_temporary_testfile.csv')
    os.remove('test_count_from_folder/__allexpr_temporary_testfile.csv')
    assert np.all(np.isclose(h_notnorm.df, truth_all_expr))

    h_norm = CountFilter.from_folder('test_count_from_folder', norm_to_rpm=True, save_csv=False)
    assert np.all(np.isclose(h_norm.df, truth_norm))

    all_feature = general.load_csv('test_count_from_folder/__allfeature_temporary_testfile.csv', 0)
    all_feature.sort_index(inplace=True)
    truth_all_feature.sort_index(inplace=True)
    try:
        assert np.all(np.isclose(all_feature, truth_all_feature))
    except AssertionError:
        os.remove('test_count_from_folder/__allfeature_temporary_testfile.csv')
        raise AssertionError
    os.remove('test_count_from_folder/__allfeature_temporary_testfile.csv')


def test_biotypes():
    truth = general.load_csv('biotypes_truth.csv', 0)
    df = CountFilter(r'counted_biotype.csv').biotypes(ref='biotype_ref_table_for_tests.csv')
    truth.sort_index(inplace=True)
    df.sort_index(inplace=True)
    assert np.all(df == truth)


def test_filter_by_row_sum():
    truth = general.load_csv('test_filter_row_sum.csv', 0)
    h = CountFilter('counted.csv')
    h.filter_by_row_sum(29)
    h.df.sort_index(inplace=True)
    truth.sort_index(inplace=True)
    assert np.all(h.df == truth)


def test_sort_inplace():
    c = CountFilter('counted.csv')
    c.sort(by='cond3', ascending=True, inplace=True)
    assert c.df['cond3'].is_monotonic_increasing


def test_sort_not_inplace():
    c = CountFilter('counted.csv')
    c_copy = general.load_csv('counted.csv', 0)
    c_sorted = c.sort(by='cond3', ascending=True, inplace=False)
    assert c_sorted.df['cond3'].is_monotonic_increasing
    assert np.all(c.df == c_copy)


def test_sort_by_multiple_columns():
    truth = general.load_csv('counted_sorted_multiple_truth.csv', 0)
    c = CountFilter('counted.csv')
    c.sort(by=['cond3', 'cond4', 'cond1', 'cond2'], ascending=[True, False, True, False], inplace=True)
    assert np.all(truth == c.df)


def test_sort_with_na_first():
    truth_first = general.load_csv('test_deseq_with_nan_sorted_nanfirst_truth.csv', 0)
    truth_last = general.load_csv('test_deseq_with_nan_sorted_nanlast_truth.csv', 0)
    c = CountFilter('test_deseq_with_nan.csv')
    c.sort(by='padj', ascending=True, na_position='first', inplace=True)
    assert truth_first.equals(c.df)
    c.sort(by='padj', ascending=True, na_position='last', inplace=True)
    assert truth_last.equals(c.df)


def test_sort_descending():
    c = CountFilter('counted.csv')
    c.sort(by='cond3', ascending=False, inplace=True)
    assert c.df['cond3'].is_monotonic_decreasing
