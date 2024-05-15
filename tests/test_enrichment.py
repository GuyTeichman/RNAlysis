import os
from unittest.mock import patch, Mock

import pytest
import statsmodels.stats.multitest as multitest

from rnalysis.enrichment import *
from rnalysis.enrichment import _fetch_sets
from tests import __attr_ref__, __biotype_ref__, is_uniprot_available, is_ensembl_available

matplotlib.use('Agg')

ENSEMBL_AVAILABLE = is_ensembl_available()
UNIPROT_AVAILABLE = is_uniprot_available()

up_feature_set = {'WBGene00021187', 'WBGene00195184', 'WBGene00012851', 'WBGene00022486', 'WBGene00011964',
                  'WBGene00012848', 'WBGene00020817', 'WBGene00012452', 'WBGene00016635', 'WBGene00044478',
                  'WBGene00018688', 'WBGene00007489', 'WBGene00019899', 'WBGene00022039', 'WBGene00021188',
                  'WBGene00007523', 'WBGene00195185', 'WBGene00206362', 'WBGene00009453', 'WBGene00017612',
                  'WBGene00018397', 'WBGene00012818', 'WBGene00018204', 'WBGene00018608', 'WBGene00022730',
                  'WBGene00021055', 'WBGene00021189', 'WBGene00007531', 'WBGene00185118', 'WBGene00195186',
                  'WBGene00021019', 'WBGene00001119', 'WBGene00044149', 'WBGene00004120', 'WBGene00013779',
                  'WBGene00044258', 'WBGene00021605', 'WBGene00010067', 'WBGene00017930', 'WBGene00012455',
                  'WBGene00013816', 'WBGene00022728', 'WBGene00206529', 'WBGene00022438', 'WBGene00017631',
                  'WBGene00194708', 'WBGene00018394', 'WBGene00050910', 'WBGene00012909', 'WBGene00018690',
                  'WBGene00007722', 'WBGene00021607', 'WBGene00194982', 'WBGene00206507', 'WBGene00044502',
                  'WBGene00021186', 'WBGene00010769', 'WBGene00008812', 'WBGene00010100', 'WBGene00044439',
                  'WBGene00018252', 'WBGene00022731', 'WBGene00194699', 'WBGene00000443', 'WBGene00010102',
                  'WBGene00012961', 'WBGene00044559', 'WBGene00007674', 'WBGene00011777', 'WBGene00021589',
                  'WBGene00016553', 'WBGene00015321', 'WBGene00019174', 'WBGene00017629', 'WBGene00007091',
                  'WBGene00010507', 'WBGene00008051', 'WBGene00045382', 'WBGene00206492', 'WBGene00006928',
                  'WBGene00009518', 'WBGene00012819', 'WBGene00021375', 'WBGene00015492', 'WBGene00008447',
                  'WBGene00017419', 'WBGene00016520', 'WBGene00017225', 'WBGene00044200', 'WBGene00206390',
                  'WBGene00022523'}


def test_featureset_api():
    FeatureSet(up_feature_set)


def test_featureset_copy():
    s = FeatureSet(up_feature_set, 'name')
    s2 = s.__copy__()
    assert s == s2
    assert s is not s2
    assert s.gene_set is not s2.gene_set


@pytest.mark.parametrize("s1,s2,expected", [
    (FeatureSet(up_feature_set), FeatureSet(up_feature_set), True),
    (FeatureSet(up_feature_set), FeatureSet(up_feature_set.copy()), True),
    (FeatureSet(up_feature_set), FeatureSet(up_feature_set.union({'other'})), False),
    (FeatureSet(up_feature_set, 'name'), FeatureSet(up_feature_set.copy(), 'name2'), False),
    (FeatureSet(set()), FeatureSet(set()), True),
    (FeatureSet(up_feature_set, 'name'), FeatureSet(up_feature_set.union({'other'}), 'name2'), False),
    (FeatureSet(up_feature_set), up_feature_set, False)
])
def test_featureset_eq(s1, s2, expected):
    assert (s1 == s2) == expected


@pytest.mark.parametrize("s", [FeatureSet(set(), 'name'), FeatureSet(up_feature_set)])
def test_featureset_iter(s):
    assert set(s.__iter__()) == s.gene_set


def test_featureset_change_set_name():
    en = FeatureSet(up_feature_set, set_name='up feature set')
    en.change_set_name('different name')
    assert en.set_name == 'different name'
    en.change_set_name('')
    assert en.set_name == ''

    en = FeatureSet(up_feature_set)
    en.change_set_name('different name')
    assert en.set_name == 'different name'

    with pytest.raises(AssertionError):
        en.change_set_name(5)


def test_featureset_contains():
    truth = {'WBGene00017419', 'WBGene00016520', 'WBGene00017225', 'WBGene00044200', 'WBGene00206390',
             'WBGene00022523', 'WBGene00000001', 'WBGene00000002'}
    f = FeatureSet(truth, 'set name')
    for ind in truth:
        assert ind in f
    assert 'set name' not in f
    assert 'WBGene00000003' not in f


def test_featureset_len():
    l = 20
    f = FeatureSet(set([str(i) for i in range(l)]), 'set name')
    assert len(f) == l


def test_featureset_union():
    other = {'WBGene00017419', 'WBGene00016520', 'WBGene00017225', 'WBGene00044200', 'WBGene00206390',
             'WBGene00022523', 'WBGene00000001', 'WBGene00000002'}
    truth = up_feature_set.union(other)
    up = FeatureSet(up_feature_set)
    other_set = FeatureSet(other)
    union_res = up.union(other_set)
    assert np.all(union_res.gene_set == truth)


def test_featureset_intersection():
    other = {'WBGene00017419', 'WBGene00016520', 'WBGene00017225', 'WBGene00044200', 'WBGene00206390',
             'WBGene00022523', 'WBGene00000001', 'WBGene00000002'}
    truth = {'WBGene00017419', 'WBGene00016520', 'WBGene00017225', 'WBGene00044200', 'WBGene00206390',
             'WBGene00022523'}
    up = FeatureSet(up_feature_set)
    other_set = FeatureSet(other)
    intersection_res = up.intersection(other_set)
    assert np.all(intersection_res.gene_set == truth)


def test_featureset_difference():
    other = {'WBGene00017419', 'WBGene00016520', 'WBGene00017225', 'WBGene00044200', 'WBGene00206390',
             'WBGene00022523', 'WBGene00000001', 'WBGene00000002'}
    truth = {'WBGene00000001', 'WBGene00000002'}
    up = FeatureSet(up_feature_set)
    other_set = FeatureSet(other)
    diff_res = other_set.difference(up)
    assert np.all(diff_res.gene_set == truth)


def test_featureset_symmetric_difference():
    first = {'WBGene00016520', 'WBGene00017225', 'WBGene00044200', 'WBGene00206390'}
    second = {'WBGene00044200', 'WBGene00206390',
              'WBGene00022523', 'WBGene00000001', 'WBGene00000002'}
    truth = {'WBGene00016520', 'WBGene00017225', 'WBGene00022523', 'WBGene00000001', 'WBGene00000002'}
    first_set = FeatureSet(first)
    second_set = FeatureSet(second)
    direction1 = second_set.symmetric_difference(first_set)
    direction2 = first_set.symmetric_difference(second_set)
    assert np.all(direction1.gene_set == truth)
    assert np.all(direction2.gene_set == truth)
    with pytest.raises(TypeError):
        first_set.symmetric_difference(second_set, second)


def test_set_operations_invalid_obj():
    first = {'WBGene00016520', 'WBGene00017225', 'WBGene00044200', 'WBGene00206390'}
    first_ep = FeatureSet(first)
    with pytest.raises(TypeError):
        first_ep.intersection(
            ['WBGene00044200', 'WBGene00206390', 'WBGene00022523', 'WBGene00000001', 'WBGene00000002'])


def test_set_operations_with_set():
    first = {'WBGene00016520', 'WBGene00017225', 'WBGene00044200', 'WBGene00206390'}
    second = {'WBGene00044200', 'WBGene00206390', 'WBGene00022523', 'WBGene00000001', 'WBGene00000002'}
    truth = {'WBGene00016520', 'WBGene00017225', 'WBGene00022523', 'WBGene00000001', 'WBGene00000002'}
    first_set = FeatureSet(first)
    symm_diff = first_set.symmetric_difference(second)
    assert np.all(symm_diff.gene_set == truth)


def test_biotypes_from_ref_table():
    truth = io.load_table('tests/test_files/biotypes_truth.csv', 0)
    genes = {'WBGene00048865', 'WBGene00000106', 'WBGene00000137', 'WBGene00199484', 'WBGene00268190', 'WBGene00048864',
             'WBGene00268189', 'WBGene00268195', 'WBGene00255734', 'WBGene00199485', 'WBGene00048863', 'WBGene00000019',
             'WBGene00268191', 'WBGene00000041', 'WBGene00199486', 'WBGene00255735', 'WBGene00000105',
             'index_that_is_not_in_biotype_ref_table'}

    en = FeatureSet(genes)
    df = en.biotypes_from_ref_table(ref=__biotype_ref__)
    df.sort_index(inplace=True)
    truth.sort_index(inplace=True)
    assert np.all(df == truth)


def test_biotypes_from_gtf():
    truth = io.load_table('tests/test_files/biotypes_truth.csv', 0)
    genes = {'WBGene00048865', 'WBGene00000106', 'WBGene00000137', 'WBGene00199484', 'WBGene00268190', 'WBGene00048864',
             'WBGene00268189', 'WBGene00268195', 'WBGene00255734', 'WBGene00199485', 'WBGene00048863', 'WBGene00000019',
             'WBGene00268191', 'WBGene00000041', 'WBGene00199486', 'WBGene00255735', 'WBGene00000105',
             'index_that_is_not_in_biotype_ref_table'}

    en = FeatureSet(genes)
    df = en.biotypes_from_gtf('tests/test_files/test_gtf_for_biotypes.gtf')
    df.sort_index(inplace=True)
    truth.sort_index(inplace=True)
    assert np.all(df == truth)


def test_enrichment_randomization_reliability():
    genes = {'WBGene00000041', 'WBGene00002074', 'WBGene00000019', 'WBGene00000105', 'WBGene00000106', 'WBGene00199484',
             'WBGene00001436', 'WBGene00000137', 'WBGene00001996', 'WBGene00014208'}
    attrs = ['attribute1', 'attribute2', 'attribute4']
    en = FeatureSet(gene_set=genes, set_name='test_set')
    random_seed = 0

    for i in range(2):
        res1 = en.user_defined_enrichment(attrs, statistical_test='randomization', biotype='all',
                                          attr_ref_path=__attr_ref__, biotype_ref_path=__biotype_ref__,
                                          randomization_reps=10000, random_seed=random_seed,
                                          parallel_backend='sequential')
        res2 = en.user_defined_enrichment(attrs, statistical_test='randomization', biotype='all',
                                          attr_ref_path=__attr_ref__, biotype_ref_path=__biotype_ref__,
                                          randomization_reps=10000, random_seed=random_seed + 1,
                                          parallel_backend='sequential')
        res3 = en.user_defined_enrichment(attrs, statistical_test='randomization', biotype='all',
                                          attr_ref_path=__attr_ref__, biotype_ref_path=__biotype_ref__,
                                          randomization_reps=10000, random_seed=random_seed + 2,
                                          parallel_backend='sequential')
        random_seed += 3
        plt.close('all')
        for col in ['samples', 'obs', 'exp', 'log2_fold_enrichment']:
            assert np.all(res1[col] == res2[col])
            assert np.all(res2[col] == res3[col])
        for randcol in ['pval', 'padj']:
            try:
                assert np.isclose(res1[randcol], res2[randcol], atol=4 * 10 ** -4, rtol=0.2).all()
                assert np.isclose(res2[randcol], res3[randcol], atol=4 * 10 ** -4, rtol=0.2).all()
                assert np.isclose(res2[randcol], res1[randcol], atol=4 * 10 ** -4, rtol=0.2).all()
                assert np.isclose(res3[randcol], res2[randcol], atol=4 * 10 ** -4, rtol=0.2).all()
            except AssertionError:
                print(res1)
                print(res2)
                print(res3)
                assert False


def _enrichment_validity(res, truth):
    for col in ['samples', 'obs', 'significant']:
        assert np.all(res[col] == truth[col])
    for closecol in ['exp', 'log2_fold_enrichment']:
        assert np.isclose(res[closecol], truth[closecol], atol=0.0).all()
    for randcol in ['pval']:
        assert np.isclose(res[randcol], truth[randcol], atol=2 * 10 ** -4, rtol=0.25).all()
    pvals = res['pval'].values
    _, padj_truth = multitest.fdrcorrection(pvals, 0.1)
    assert np.isclose(res['padj'].values, padj_truth, atol=0.0).all()


def test_enrichment_randomization_validity():
    truth = io.load_table('tests/test_files/enrichment_randomization_res.csv', 0)
    genes = {'WBGene00000041', 'WBGene00002074', 'WBGene00000105', 'WBGene00000106', 'WBGene00199484',
             'WBGene00001436', 'WBGene00000137', 'WBGene00001996', 'WBGene00014208', 'WBGene00001133'}
    attrs = ['attribute1', 'attribute2']
    en = FeatureSet(gene_set=genes, set_name='test_set')
    res = en.user_defined_enrichment(attrs, statistical_test='randomization', biotype='all', attr_ref_path=__attr_ref__,
                                     biotype_ref_path=__biotype_ref__, randomization_reps=100000, random_seed=0,
                                     parallel_backend='sequential')
    plt.close('all')
    _enrichment_validity(res, truth)


def test_enrichment_randomization_parallel_reliability():
    genes = {'WBGene00000041', 'WBGene00002074', 'WBGene00000019', 'WBGene00000105', 'WBGene00000106', 'WBGene00199484',
             'WBGene00001436', 'WBGene00000137', 'WBGene00001996', 'WBGene00014208'}
    attrs = ['attribute1', 'attribute2', 'attribute4']
    en = FeatureSet(gene_set=genes, set_name='test_set')
    random_seed = 0

    for i in range(2):
        res1 = en.user_defined_enrichment(attrs, statistical_test='randomization', biotype='all',
                                          attr_ref_path=__attr_ref__, biotype_ref_path=__biotype_ref__,
                                          randomization_reps=5000, random_seed=random_seed, parallel_backend='loky')
        res2 = en.user_defined_enrichment(attrs, statistical_test='randomization', biotype='all',
                                          attr_ref_path=__attr_ref__, biotype_ref_path=__biotype_ref__,
                                          randomization_reps=5000, random_seed=random_seed + 1,
                                          parallel_backend='multiprocessing')
        res3 = en.user_defined_enrichment(attrs, statistical_test='randomization', biotype='all',
                                          attr_ref_path=__attr_ref__, biotype_ref_path=__biotype_ref__,
                                          randomization_reps=5000, random_seed=random_seed + 2,
                                          parallel_backend='threading')
        random_seed += 3
        plt.close('all')
        for col in ['samples', 'obs', 'exp', 'log2_fold_enrichment']:
            assert np.all(res1[col] == res2[col])
            assert np.all(res2[col] == res3[col])
        for randcol in ['pval']:
            assert np.isclose(res1[randcol], res2[randcol], atol=4 * 10 ** -4, rtol=0.25).all()
            assert np.isclose(res2[randcol], res3[randcol], atol=4 * 10 ** -4, rtol=0.25).all()
            assert np.isclose(res2[randcol], res1[randcol], atol=4 * 10 ** -4, rtol=0.25).all()
            assert np.isclose(res3[randcol], res2[randcol], atol=4 * 10 ** -4, rtol=0.25).all()


def test_enrichment_parallel_validity():
    truth = io.load_table('tests/test_files/enrichment_randomization_res.csv', 0)
    genes = {'WBGene00000041', 'WBGene00002074', 'WBGene00000105', 'WBGene00000106', 'WBGene00199484',
             'WBGene00001436', 'WBGene00000137', 'WBGene00001996', 'WBGene00014208', 'WBGene00001133'}
    attrs = ['attribute1', 'attribute2']
    en = FeatureSet(gene_set=genes, set_name='test_set')
    res = en.user_defined_enrichment(attrs, statistical_test='randomization', biotype='all', attr_ref_path=__attr_ref__,
                                     biotype_ref_path=__biotype_ref__, randomization_reps=100000, random_seed=0,
                                     parallel_backend="loky")
    plt.close('all')
    _enrichment_validity(res, truth)


def test_enrich_hypergeometric_pvalues():
    truth = io.load_table('tests/test_files/enrichment_hypergeometric_res.csv', 0)
    genes = {'WBGene00000041', 'WBGene00002074', 'WBGene00000105', 'WBGene00000106', 'WBGene00199484',
             'WBGene00001436', 'WBGene00000137', 'WBGene00001996', 'WBGene00014208', 'WBGene00001133'}
    attrs = ['attribute1', 'attribute2']
    en = FeatureSet(gene_set=genes, set_name='test_set')
    res = en.user_defined_enrichment(attrs, statistical_test='hypergeometric', biotype='all',
                                     attr_ref_path=__attr_ref__, biotype_ref_path=__biotype_ref__)
    plt.close('all')
    _enrichment_validity(res, truth)


def test_save_txt():
    try:
        geneset = {'gene1', 'gene2', 'gene3', 'gene5'}
        en = FeatureSet(geneset, 'my gene set')
        en.save_txt('tests/test_files/tmp_enrichment_txt')

        with open('tests/test_files/tmp_enrichment_txt.txt') as f:
            loaded_geneset = {gene.replace('\n', '') for gene in f}
        assert loaded_geneset == geneset
    except Exception as e:
        raise e
    finally:
        try:
            os.remove('tests/test_files/tmp_enrichment_txt.txt')
        except FileNotFoundError:
            pass


def test_featureset_from_string(monkeypatch):
    truth = {'gene1', 'gene2', 'gene 5'}
    monkeypatch.setattr('builtins.input', lambda x: 'gene1\ngene2\ngene 5\n')
    en = FeatureSet(None)
    assert en.gene_set == truth


def test_featureset_repr():
    en = FeatureSet({"1", "2", "4", "5"}, 'my very important set')
    assert repr(en) == "FeatureSet: 'my very important set'"


@pytest.mark.parametrize('kwargs', [{},
                                    dict(plot_horizontal=True, ylim=12),
                                    dict(plot_horizontal=False, ylabel='different ylabel', alpha=0.1, ylim=0.42)])
def test_enrichment_bar_plot(kwargs):
    pth = 'tests/test_files/enrichment_hypergeometric_res.csv'
    enrichment_bar_plot(pth, **kwargs)
    plt.close('all')


def test_rankedset_api():
    en = RankedSet(['1', '9', '4'], 'name')
    assert en.gene_set == {'1', '9', '4'}
    assert (en.ranked_genes == np.array(['1', '9', '4'], dtype='str')).all()
    assert en.set_name == 'name'


def test_rankedset_copy():
    s = RankedSet(['a', 'b', 'd', 'c'], 'name')
    s2 = s.__copy__()
    assert s == s2
    assert s is not s2
    assert s.gene_set is not s2.gene_set
    assert s.ranked_genes is not s2.ranked_genes


@pytest.mark.parametrize("s1,s2,expected", [
    (RankedSet(['a', 'b', 'c']), RankedSet(['a', 'b', 'c']), True),
    (RankedSet(['a', 'b', 'c'], 'name'), RankedSet(['a', 'b', 'c'], 'name'), True),
    (RankedSet(['a', 'b', 'c']), FeatureSet(['a', 'b', 'c']), False),
    (RankedSet(['a', 'b', 'c']), RankedSet(['a', 'b', 'c', 'd']), False),
    (RankedSet(['a', 'b', 'c']), RankedSet(['a', 'c', 'b']), False),
    (RankedSet(['a', 'b', 'c'], 'name'), RankedSet(['a', 'b', 'c'], 'name2'), False),
    (RankedSet([]), RankedSet([]), True),
    (RankedSet(['a', 'b', 'c'], 'name'), RankedSet(['a', 'd', 'b', 'c'], 'name2'), False),
    (RankedSet(['a', 'b', 'c']), ['a', 'b', 'c'], False)
])
def test_rankedset_eq(s1, s2, expected):
    assert (s1 == s2) == expected


@pytest.mark.parametrize("s", [RankedSet([], 'name'), RankedSet(['1', '2', '3', '4', '5']),
                               RankedSet(['a', 'c', 'b', 'e', 'd'], 'othername')])
def test_ranked_iter(s):
    assert list(s.__iter__()) == list(s.ranked_genes)


def test_rankedset_repr():
    en = RankedSet(['1', '9', '4'], 'my very important set')
    assert repr(en) == "RankedSet: 'my very important set'"


def test_rankedset_set_ops_return_type():
    r1 = RankedSet(['1', '2', '3', '4'])
    r2 = RankedSet(['3', '4', '5'])
    f1 = FeatureSet(['3', '5', '6', '7', '9'])
    s1 = {'1', '2', '3', '5'}
    assert isinstance(r1.intersection(r2), FeatureSet)
    assert isinstance(r2.union(f1), FeatureSet)
    assert isinstance(r1.difference(s1), FeatureSet)
    assert isinstance(f1.symmetric_difference(r2), FeatureSet)


def _comp_go_res_df(res, truth):
    res.drop('name', axis=1, inplace=True)
    res.rename_axis('go_id')
    assert res.loc[:, ['n', 'obs']].equals(truth.loc[:, ['n', 'obs']])
    assert np.allclose(res.loc[:, ['exp', 'log2fc']], res.loc[:, ['exp', 'log2fc']])
    assert np.allclose(res['pval'], truth['pval'], atol=0)


def test_rankedset_from_ndarray():
    gene_list = ['gene1', 'gene2', 'gene3', 'gene13']
    arr = np.array(gene_list)
    r = RankedSet(arr)
    for i, j in zip(gene_list, r.ranked_genes):
        assert i == j


def test_rankedset_init_invalid_type():
    with pytest.raises(TypeError):
        _ = RankedSet({1, 2, 3})
    with pytest.raises(TypeError):
        _ = RankedSet({1: 2, 3: 4})
    with pytest.raises(TypeError):
        _ = RankedSet(5)


def test_enrich_non_categorical_api():
    ref_table = 'tests/test_files/attr_ref_table_for_non_categorical.csv'
    s = FeatureSet(
        {'WBGene00000001', 'WBGene00000002', 'WBGene00000003', 'WBGene00000004', 'WBGene00000005', 'WBGene00000006',
         'WBGene00000007', 'WBGene00000008', 'WBGene00000009', 'WBGene00000010', 'WBGene00000011', 'WBGene00000012'})
    res = s.non_categorical_enrichment('attr1', biotype='all', parametric_test=True, attr_ref_path=ref_table,
                                       biotype_ref_path=__biotype_ref__)
    assert isinstance(res, pd.DataFrame)
    res = s.non_categorical_enrichment(['attr3', 'attr2', 'attr4'], biotype='protein_coding', attr_ref_path=ref_table,
                                       biotype_ref_path=__biotype_ref__, plot_log_scale=False, n_bins=30,
                                       plot_style='interleaved', return_fig=True)
    assert isinstance(res, tuple)
    assert isinstance(res[0], pd.DataFrame)
    assert isinstance(res[1], list)
    for fig in res[1]:
        assert isinstance(fig, plt.Figure)
    plt.close('all')


def test_enrich_non_categorial_parametric_test():
    ref_table = 'tests/test_files/attr_ref_table_for_non_categorical.csv'
    res_param_truth = pd.read_csv('tests/test_files/enrich_non_categorical_parametric_truth_fdr15.csv', index_col=0)
    geneset = {'WBGene00000001', 'WBGene00000003', 'WBGene00000004', 'WBGene00000008', 'WBGene00000010',
               'WBGene00000014', 'WBGene00000015', 'WBGene00000020'}

    fs = FeatureSet(geneset)
    res_param = fs.non_categorical_enrichment('all', alpha=0.15, parametric_test=True, attr_ref_path=ref_table,
                                              biotype='all')
    assert res_param.loc[:, ['samples', 'significant']].equals(res_param_truth.loc[:, ['samples', 'significant']])
    assert np.isclose(res_param.loc[:, ['obs', 'exp']].values, res_param_truth.loc[:, ['obs', 'exp']].values,
                      rtol=0.01).all()
    assert np.isclose(res_param.loc[:, ['pval', 'padj']].values, res_param_truth.loc[:, ['pval', 'padj']], atol=0,
                      rtol=0.02).all()


def test_enrich_non_categorial_nonparametric_test():
    ref_table = 'tests/test_files/attr_ref_table_for_non_categorical.csv'
    res_nonparam_truth = pd.read_csv('tests/test_files/enrich_non_categorical_aparametric_truth_fdr15.csv', index_col=0)
    geneset = {'WBGene00000001', 'WBGene00000003', 'WBGene00000004', 'WBGene00000008', 'WBGene00000010',
               'WBGene00000014', 'WBGene00000015', 'WBGene00000020'}

    fs = FeatureSet(geneset)

    res_nonparam = fs.non_categorical_enrichment('all', alpha=0.15, attr_ref_path=ref_table, biotype='all')
    assert res_nonparam.loc[:, ['samples', 'significant']].equals(res_nonparam_truth.loc[:, ['samples', 'significant']])
    assert np.isclose(res_nonparam.loc[:, ['obs', 'exp']].values, res_nonparam_truth.loc[:, ['obs', 'exp']].values,
                      rtol=0.01).all()
    assert np.isclose(res_nonparam.loc[:, ['pval', 'padj']].values, res_nonparam_truth.loc[:, ['pval', 'padj']], atol=0,
                      rtol=0.02).all()
    plt.close('all')


def test_enrich_non_categorical_nan_values():
    ref_table = 'tests/test_files/attr_ref_table_for_non_categorical.csv'
    geneset = {'WBGene00000025', 'WBGene00000023', 'WBGene00000027', 'WBGene00000028', 'WBGene00000029',
               'WBGene00000030', 'WBGene00000015', 'WBGene00000020'}
    truth = pd.read_csv('tests/test_files/enrich_non_categorical_nan_truth.csv', index_col=0)
    truth_param = pd.read_csv('tests/test_files/enrich_non_categorical_nan_parametric_truth.csv', index_col=0)
    fs = FeatureSet(geneset)
    attrs = ['attr1', 'attr5', 'attr2', 'attr3', 'attr4']
    res_param = fs.non_categorical_enrichment(attrs, parametric_test=True, attr_ref_path=ref_table, biotype='all')
    res_nonparam = fs.non_categorical_enrichment(attrs, attr_ref_path=ref_table, biotype='all')
    print(res_param, res_nonparam)

    for df, df_truth in zip((res_nonparam, res_param), (truth, truth_param)):
        assert df.loc[:, ['samples', 'significant']].equals(
            df_truth.loc[:, ['samples', 'significant']])
        assert np.isclose(df.loc[:, ['obs', 'exp']].values, df_truth.loc[:, ['obs', 'exp']].values,
                          rtol=0.01, equal_nan=True).all()
        assert np.isclose(df.loc[:, ['pval', 'padj']].values, df_truth.loc[:, ['pval', 'padj']],
                          atol=0, rtol=0.02, equal_nan=True).all()
    plt.close('all')


def test_enrich_single_set_api():
    genes_ranked = ['WBGene00000019', 'WBGene00000041', 'WBGene00000105', 'WBGene00000106', 'WBGene00000137',
                    'WBGene00001436', 'WBGene00001996', 'WBGene00002074', 'WBGene00003864', 'WBGene00003865',
                    'WBGene00003902', 'WBGene00003915', 'WBGene00000369', 'WBGene00000859', 'WBGene00000860',
                    'WBGene00000861', 'WBGene00000863', 'WBGene00000864', 'WBGene00000865', 'WBGene00001131',
                    'WBGene00001132', 'WBGene00001133', 'WBGene00001134', 'WBGene00048863', 'WBGene00048864',
                    'WBGene00048865', 'WBGene00199484', 'WBGene00199485', 'WBGene00199486', 'WBGene00255734',
                    'WBGene00255735', 'WBGene00268189', 'WBGene00268190', 'WBGene00268191', 'WBGene00268195',
                    'WBGene00004920', 'WBGene00011910', 'WBGene00014208']

    attrs = ['attribute1', 'attribute2']
    en = RankedSet(genes_ranked, set_name='test_set')
    _ = en.single_set_enrichment(attrs, attr_ref_path=__attr_ref__)
    plt.close('all')


@pytest.mark.skipif(not ENSEMBL_AVAILABLE, reason='Ensembl REST API is not available at the moment')
@pytest.mark.skipif(not UNIPROT_AVAILABLE, reason='UniProt REST API is not available at the moment')
@pytest.mark.parametrize('return_nonsignificant', [True, False])
@pytest.mark.parametrize("organism,propagate_annotations", [(6239, 'classic'), ('caenorhabditis elegans', 'elim')])
def test_go_enrichment_single_set_api(organism, propagate_annotations, return_nonsignificant):
    genes_ranked = ['WBGene00000019', 'WBGene00000041', 'WBGene00000105', 'WBGene00000106', 'WBGene00000137',
                    'WBGene00001436', 'WBGene00001996', 'WBGene00002074', 'WBGene00003864', 'WBGene00003865',
                    'WBGene00003902', 'WBGene00003915', 'WBGene00000369', 'WBGene00000859', 'WBGene00000860',
                    'WBGene00000861', 'WBGene00000863', 'WBGene00000864', 'WBGene00000865', 'WBGene00001131',
                    'WBGene00001132', 'WBGene00001133', 'WBGene00001134', 'WBGene00048863', 'WBGene00048864',
                    'WBGene00048865', 'WBGene00199484', 'WBGene00199485', 'WBGene00199486', 'WBGene00255734',
                    'WBGene00255735', 'WBGene00268189', 'WBGene00268190', 'WBGene00268191', 'WBGene00268195',
                    'WBGene00004920', 'WBGene00011910', 'WBGene00014208']

    en = RankedSet(genes_ranked, set_name='test_set')
    _ = en.single_set_go_enrichment(organism, 'WBGene', propagate_annotations=propagate_annotations,
                                    aspects='biological_process', evidence_types='experimental', databases='WB',
                                    return_nonsignificant=return_nonsignificant)
    plt.close('all')


@pytest.mark.skipif(not ENSEMBL_AVAILABLE, reason='Ensembl REST API is not available at the moment')
@pytest.mark.skipif(not UNIPROT_AVAILABLE, reason='UniProt REST API is not available at the moment')
@pytest.mark.parametrize("organism,statistical_test,propagate_annotations,kwargs",
                         [(6239, 'fisher', 'elim', {}),
                          ('caenorhabditis elegans', 'randomization', 'no', dict(randomization_reps=100))])
def test_go_enrichment_api(organism, statistical_test, propagate_annotations, kwargs):
    genes = {'WBGene00048865', 'WBGene00000864', 'WBGene00000105', 'WBGene00001996', 'WBGene00011910', 'WBGene00268195'}
    en = FeatureSet(gene_set=genes, set_name='test_set')
    _ = en.go_enrichment(organism, 'WBGene', statistical_test=statistical_test, biotype='protein_coding',
                         biotype_ref_path=__biotype_ref__, propagate_annotations=propagate_annotations,
                         aspects='biological_process', evidence_types='IMP', databases='WB', **kwargs)
    plt.close('all')


def test_upset_plot_api():
    upset_plot({'obj': {'0', '1', '2'}, 'obj2': {'1', '3'}, 'obj3': {'5', '6', '7', '0'}, 'obj4': {'0', '3', '5'}})
    upset_plot({'obj': {'0', '1', '2'}, 'obj2': FeatureSet({'1', '2', '3'}), 'obj3': RankedSet(['5', '3', '0'])})
    upset_plot(
        {'obj': 'attribute1', 'obj2': 'attribute2', 'obj3': {'WBGene00000001', 'WBGene00001234'}, 'obj4': 'attribute3',
         'obj5': 'attribute4'},
        attr_ref_table_path=__attr_ref__)
    plt.close('all')


def test_venn_diagram_api():
    venn_diagram({'obj': {'0', '1', '2'}, 'obj2': {'1', '3'}, 'obj3': {'5', '6', '7', '0'}})
    venn_diagram({'obj': {'0', '1', '2'}, 'obj2': FeatureSet({'1', '2', '3'}), 'obj3': RankedSet(['5', '3', '0'])},
                 title='title', set_colors=('black', 'purple', 'yellow'), transparency=0.2, linecolor='grey',
                 linestyle='dashed', linewidth=1)
    venn_diagram({'obj': 'attribute1', 'obj2': 'attribute2', 'obj3': {'WBGene00000001', 'WBGene00001234'}},
                 attr_ref_table_path=__attr_ref__, add_outline=False)
    venn_diagram({'obj': {'0', '1', '2'}, 'otherobj': {'2', '3', '5'}})
    venn_diagram({'obj': {'0', '1', '2'}, 'otherobj': {'2', '3', '5'}}, weighted=False)
    plt.close('all')


def test_venn_diagram_invalid_number_of_sets():
    with pytest.raises(ValueError):
        venn_diagram({'obj': {'0', '1', '2'}, 'obj2': {'1', '3'}, 'obj3': {'5', '6', '7', '0'}, 'obj4': {'3', '4'}})
    with pytest.raises(ValueError):
        venn_diagram({'obj1': FeatureSet({'1', '2'})})
    with pytest.raises(ValueError):
        venn_diagram({'obj1': FeatureSet({'1', '2'}), 'obj2': {}, 'obj3': 'attr3', 'obj4': 'attr4', 'obj5': 'attr2'})


@pytest.mark.parametrize('objs,truth', [
    ({'set1': ['one', 'two'], 'set2': {'one', 'two', 'three'}, 'set3': FeatureSet({'two', 'three'})},
     {'set1': {'one', 'two'}, 'set2': {'one', 'two', 'three'}, 'set3': {'two', 'three'}}),
    ({'set1': RankedSet(['one', 'two']), 'set2': Filter('tests/test_files/test_fetch_sets_table.csv'),
      'set3': 'attribute2'},
     {'set1': {'one', 'two'}, 'set2': {'one', 'two', 'three'},
      'set3': {'WBGene00000369', 'WBGene00003864', 'WBGene00003865', 'WBGene00003902', 'WBGene00003915',
               'WBGene00004920', 'WBGene00011910', 'WBGene00014208'}})])
def test_fetch_sets(objs: dict, truth: dict):
    objs_original = objs.copy()
    res = _fetch_sets(objs, __attr_ref__)
    assert res == truth
    assert dict(objs_original) == objs


@pytest.mark.parametrize('objs', [{'set1': {1, 2, 3}, 'set2': True}])
def test_fetch_sets_bad_type(objs: dict):
    with pytest.raises(TypeError):
        _ = _fetch_sets(objs, __attr_ref__)


@pytest.mark.skipif(not ENSEMBL_AVAILABLE, reason='Ensembl REST API is not available at the moment')
@pytest.mark.skipif(not UNIPROT_AVAILABLE, reason='UniProt REST API is not available at the moment')
@pytest.mark.parametrize("organism,statistical_test,kwargs",
                          [(6239, 'fisher', {}),
                          ('caenorhabditis elegans', 'randomization', dict(randomization_reps=100))])
def test_kegg_enrichment_api(organism, statistical_test, kwargs):
    genes = {'WBGene00048865', 'WBGene00000864', 'WBGene00000105', 'WBGene00001996', 'WBGene00011910', 'WBGene00268195',
             'WBGene00000833' 'WBGene00007150', 'WBGene00016995'}
    en = FeatureSet(gene_set=genes, set_name='test_set')
    _ = en.kegg_enrichment(organism, 'WormBase', statistical_test=statistical_test, biotype='protein_coding',
                           biotype_ref_path=__biotype_ref__, **kwargs)
    plt.close('all')


@pytest.mark.skipif(not ENSEMBL_AVAILABLE, reason='Ensembl REST API is not available at the moment')
@pytest.mark.skipif(not UNIPROT_AVAILABLE, reason='UniProt REST API is not available at the moment')
def test_kegg_enrichment_single_list_api():
    genes = ['WBGene00004049', 'WBGene00000864', 'WBGene00004050', 'WBGene00004930', 'WBGene00001686', 'WBGene00268195']
    en = RankedSet(genes, set_name='test_set')
    _ = en.single_set_kegg_enrichment(6239, 'WormBase', pathway_graphs_format='pdf', parallel_backend='sequential')
    plt.close('all')


@pytest.mark.skipif(not UNIPROT_AVAILABLE, reason='UniProt REST API is not available at the moment')
def test_kegg_pathway_graph():
    kegg_pathway_graph('cel04020', None, 'WormBase')
    kegg_pathway_graph('cel04020', {'WBGene00004039', 'WBGene00004039'}, 'WormBase')


def test_gene_ontology_graph():
    gene_ontology_graph('biological_process', 'tests/test_files/go_enrichment_runner_sample_results.csv', 'colName')


@pytest.mark.parametrize("map_to,map_from,remove_unmapped_genes,expected,expected_name", [
    ('UniProtKB AC/ID', 'WormBase', False, 'tests/test_files/counted_translated_with_unmapped.csv',
     'set_name_translateFromWormBasetoUniProtKBACID'),
    ('UniProtKB', 'auto', True, 'tests/test_files/counted_translated_remove_unmapped.csv',
     'set_name_translateFromWormBasetoUniProtKB'),
])
def test_translate_gene_ids(map_to, map_from, remove_unmapped_genes, expected, expected_name, monkeypatch):
    def mock_map_gene_ids(self, ids):
        if self.map_from == 'WormBase':
            return io.GeneIDDict(
                {'WBGene00007063': 'A0A0K3AWR5', 'WBGene00007064': 'A0A2X0T1Z3', 'WBGene00007067': 'D3NQA2',
                 'WBGene00077503': 'H2L2B5', 'WBGene00007071': 'Q17405', 'WBGene00014997': 'Q7JNR0',
                 'WBGene00043988': 'A4F2Z7', 'WBGene00043989': 'G5EFZ2', 'WBGene00007075': 'G5EDW3',
                 'WBGene00007076': 'G5EFZ2'})
        return io.GeneIDDict({})

    monkeypatch.setattr(io.GeneIDTranslator, 'run', mock_map_gene_ids)
    truth = FeatureSet(parsing.data_to_set(io.load_table(expected, index_col=0).index), set_name=expected_name)
    s = FeatureSet(
        {'WBGene00044022', 'WBGene00007075', 'WBGene00007066', 'WBGene00043988', 'WBGene00044951', 'WBGene00077503',
         'WBGene00007076', 'WBGene00007063', 'WBGene00043990', 'WBGene00007074', 'WBGene00007067', 'WBGene00007079',
         'WBGene00077502', 'WBGene00077504', 'WBGene00007064', 'WBGene00043989', 'WBGene00007071', 'WBGene00007078',
         'WBGene00007077', 'WBGene00043987', 'WBGene00007069', 'WBGene00014997'}
        , 'set_name')

    res = s.translate_gene_ids(map_to, map_from, remove_unmapped_genes, inplace=False)
    assert res == truth
    s.translate_gene_ids(map_to, map_from, remove_unmapped_genes, inplace=True)
    assert s == truth


@pytest.mark.parametrize("map_to,map_from,remove_unmapped_genes,expected,expected_name", [
    ('UniProtKB AC/ID', 'WormBase', False,
     ['WBGene00044022', 'G5EDW3', 'WBGene00007066', 'Q7JNR0', 'A4F2Z7', 'WBGene00044951', 'A0A0K3AWR5',
      'WBGene00043990', 'WBGene00007074', 'D3NQA2', 'WBGene00007079', 'WBGene00077502', 'WBGene00077504', 'A0A2X0T1Z3',
      'G5EFZ2', 'Q17405', 'WBGene00007078', 'WBGene00007077', 'WBGene00043987', 'WBGene00007069', 'H2L2B5']
     , 'set_name_translateFromWormBasetoUniProtKBACID'),
    ('UniProtKB', 'auto', True,
     ['G5EDW3', 'Q7JNR0', 'A4F2Z7', 'A0A0K3AWR5', 'D3NQA2', 'A0A2X0T1Z3', 'G5EFZ2', 'Q17405', 'H2L2B5']
     , 'set_name_translateFromWormBasetoUniProtKB'),
])
def test_translate_gene_ids_rankedset(map_to, map_from, remove_unmapped_genes, expected, expected_name, monkeypatch):
    def mock_map_gene_ids(self, ids):
        if self.map_from == 'WormBase':
            return io.GeneIDDict(
                {'WBGene00007063': 'A0A0K3AWR5', 'WBGene00007064': 'A0A2X0T1Z3', 'WBGene00007067': 'D3NQA2',
                 'WBGene00077503': 'H2L2B5', 'WBGene00007071': 'Q17405', 'WBGene00014997': 'Q7JNR0',
                 'WBGene00043988': 'A4F2Z7', 'WBGene00043989': 'G5EFZ2', 'WBGene00007075': 'G5EDW3',
                 'WBGene00007076': 'G5EFZ2'})
        return io.GeneIDDict({})

    monkeypatch.setattr(io.GeneIDTranslator, 'run', mock_map_gene_ids)
    truth = RankedSet(expected, set_name=expected_name)
    s = RankedSet(
        ['WBGene00044022', 'WBGene00007075', 'WBGene00007066', 'WBGene00014997', 'WBGene00043988', 'WBGene00044951',
         'WBGene00007063', 'WBGene00043990', 'WBGene00007074', 'WBGene00007067', 'WBGene00007079',
         'WBGene00077502', 'WBGene00077504', 'WBGene00007064', 'WBGene00043989', 'WBGene00007071', 'WBGene00007078',
         'WBGene00007077', 'WBGene00043987', 'WBGene00007069', 'WBGene00077503']
        , 'set_name')

    res = s.translate_gene_ids(map_to, map_from, remove_unmapped_genes, inplace=False)
    assert res == truth
    s.translate_gene_ids(map_to, map_from, remove_unmapped_genes, inplace=True)
    assert s == truth


def test_filter_by_attribute():
    truth = FeatureSet({'WBGene00011910', 'WBGene00000019', 'WBGene00001436'}, 'set_name_filtbyattrattribute1')
    s = FeatureSet(
        {'WBGene00000027', 'WBGene00000005', 'WBGene00000017', 'WBGene00000025', 'WBGene00000003', 'WBGene00000019',
         'WBGene00000026', 'WBGene00000028', 'WBGene00000004', 'WBGene00000024', 'WBGene00000016', 'WBGene00000010',
         'WBGene00000008', 'WBGene00000009', 'WBGene00000023', 'WBGene00000015', 'WBGene00000006', 'WBGene00000029',
         'WBGene00000012', 'WBGene00001436', 'WBGene00000022', 'WBGene00000013', 'WBGene00000014', 'WBGene00000002',
         'WBGene00000007', 'WBGene00000020', 'WBGene00000011', 'WBGene00011910'}
        , set_name='set_name')
    s_notinplace = s.filter_by_attribute('attribute1', ref=__attr_ref__, inplace=False)
    assert s_notinplace == truth

    s.filter_by_attribute('attribute1', ref=__attr_ref__)
    assert s == truth


def test_filter_by_attribute_rankedset():
    truth = RankedSet(['WBGene00000019', 'WBGene00001436', 'WBGene00011910'], 'set_name_filtbyattrattribute1')
    s = RankedSet(
        ['WBGene00000027', 'WBGene00000005', 'WBGene00000017', 'WBGene00000025', 'WBGene00000003', 'WBGene00000019',
         'WBGene00000026', 'WBGene00000028', 'WBGene00000004', 'WBGene00000024', 'WBGene00000016', 'WBGene00000010',
         'WBGene00000008', 'WBGene00000009', 'WBGene00000023', 'WBGene00000015', 'WBGene00000006', 'WBGene00000029',
         'WBGene00000012', 'WBGene00001436', 'WBGene00000022', 'WBGene00000013', 'WBGene00000014', 'WBGene00000002',
         'WBGene00000007', 'WBGene00000020', 'WBGene00000011', 'WBGene00011910']
        , set_name='set_name')
    s_notinplace = s.filter_by_attribute('attribute1', ref=__attr_ref__, inplace=False)
    assert s_notinplace == truth

    s.filter_by_attribute('attribute1', ref=__attr_ref__)
    assert s == truth


@patch('rnalysis.utils.io.PantherOrthologMapper')
def test_find_paralogs_panther(mock_mapper):
    df_truth = pd.DataFrame({'gene': ['gene1', 'gene1', 'gene2'], 'paralog': ['paralog1', 'paralog2', 'paralog3']})
    mock_mapper_instance = Mock()
    mock_mapper.return_value = mock_mapper_instance
    mock_mapper_instance.get_paralogs.return_value = io.OrthologDict(
        {'gene1': ['paralog1', 'paralog2'], 'gene2': ['paralog3']})

    featureset = FeatureSet(Filter('tests/test_files/test_map_orthologs.csv'))
    result = featureset.find_paralogs_panther(organism='Homo sapiens', gene_id_type='UniProtKB')

    assert result.equals(df_truth)

@pytest.mark.parametrize('filter_percent_identity', [True, False])
@patch('rnalysis.utils.io.EnsemblOrthologMapper')
def test_find_paralogs_ensembl(mock_mapper, filter_percent_identity):
    df_truth = pd.DataFrame({'gene': ['gene1', 'gene1', 'gene2'], 'paralog': ['paralog1', 'paralog2', 'paralog3']})
    mock_mapper_instance = Mock()
    mock_mapper.return_value = mock_mapper_instance
    mock_mapper_instance.get_paralogs.return_value = io.OrthologDict(
        {'gene1': ['paralog1', 'paralog2'], 'gene2': ['paralog3']})

    featureset = FeatureSet(Filter('tests/test_files/test_map_orthologs.csv'))
    result = featureset.find_paralogs_ensembl(organism='Homo sapiens', gene_id_type='UniProtKB',
                                              filter_percent_identity=filter_percent_identity)

    assert result.equals(df_truth)

@pytest.mark.parametrize('ranked_set', [True, False])
@pytest.mark.parametrize('remove_unmapped', [True, False])
@pytest.mark.parametrize('non_unique_mode,filter_least_diverged', [
    ('first', True),
    ('last', False),
])
@patch('rnalysis.utils.io.PantherOrthologMapper')
def test_map_orthologs_panther(mock_mapper, non_unique_mode, filter_least_diverged, remove_unmapped, ranked_set):
    if remove_unmapped:
        filter_obj_path = 'tests/test_files/test_map_orthologs_remove_unmapped_truth.csv'
    else:
        filter_obj_path = 'tests/test_files/test_map_orthologs_truth.csv'
    filter_obj_truth = Filter(filter_obj_path)
    one2many_truth = pd.DataFrame(
        {'gene': ['gene1', 'gene1', 'gene2'], 'ortholog': ['ortholog1', 'ortholog2', 'ortholog3']})
    mock_mapper_instance = Mock()
    mock_mapper.return_value = mock_mapper_instance
    mock_mapper_instance.get_orthologs.return_value = (
        io.OrthologDict({'gene1': 'ortholog2', 'gene2': 'ortholog3'})
        , io.OrthologDict({'gene1': ['ortholog1', 'ortholog2'], 'gene2': ['ortholog3']}))

    if ranked_set:
        featureset = RankedSet(Filter('tests/test_files/test_map_orthologs.csv'))
    else:
        featureset = FeatureSet(Filter('tests/test_files/test_map_orthologs.csv'))
    one2one, one2many = featureset.map_orthologs_panther('Homo sapiens', 'caenorhabditis elegans',
                                                         gene_id_type='UniProtKB', inplace=False,
                                                         non_unique_mode=non_unique_mode,
                                                         filter_least_diverged=filter_least_diverged,
                                                         remove_unmapped_genes=remove_unmapped)

    if ranked_set:
        one2one_truth = RankedSet(filter_obj_truth, one2one.set_name)
    else:
        one2one_truth = FeatureSet(filter_obj_truth, one2one.set_name)

    assert one2one == one2one_truth
    assert one2many.equals(one2many_truth)

@pytest.mark.parametrize('ranked_set', [True, False])
@pytest.mark.parametrize('remove_unmapped', [True, False])
@pytest.mark.parametrize('non_unique_mode,filter_percent_identity', [
    ('first', True),
    ('last', False),
])
@patch('rnalysis.utils.io.EnsemblOrthologMapper')
def test_map_orthologs_ensembl(mock_mapper, non_unique_mode, filter_percent_identity, remove_unmapped, ranked_set):
    if remove_unmapped:
        filter_obj_path = 'tests/test_files/test_map_orthologs_remove_unmapped_truth.csv'
    else:
        filter_obj_path = 'tests/test_files/test_map_orthologs_truth.csv'
    filter_obj_truth = Filter(filter_obj_path)
    one2many_truth = pd.DataFrame(
        {'gene': ['gene1', 'gene1', 'gene2'], 'ortholog': ['ortholog1', 'ortholog2', 'ortholog3']})
    mock_mapper_instance = Mock()
    mock_mapper.return_value = mock_mapper_instance
    mock_mapper_instance.get_orthologs.return_value = (
        io.OrthologDict({'gene1': 'ortholog2', 'gene2': 'ortholog3'})
        , io.OrthologDict({'gene1': ['ortholog1', 'ortholog2'], 'gene2': ['ortholog3']}))

    if ranked_set:
        featureset = RankedSet(Filter('tests/test_files/test_map_orthologs.csv'))
    else:
        featureset = FeatureSet(Filter('tests/test_files/test_map_orthologs.csv'))

    one2one, one2many = featureset.map_orthologs_ensembl('Homo sapiens', 'caenorhabditis elegans',
                                                         gene_id_type='UniProtKB', inplace=False,
                                                         non_unique_mode=non_unique_mode,
                                                         filter_percent_identity=filter_percent_identity,
                                                         remove_unmapped_genes=remove_unmapped)

    if ranked_set:
        one2one_truth = RankedSet(filter_obj_truth, one2one.set_name)
    else:
        one2one_truth = FeatureSet(filter_obj_truth, one2one.set_name)

    assert one2one == one2one_truth
    assert one2many.equals(one2many_truth)

@pytest.mark.parametrize('ranked_set', [True, False])
@pytest.mark.parametrize('remove_unmapped', [True, False])
@pytest.mark.parametrize('non_unique_mode,filter_consistency_score,threshold', [
    ('first', True, 0.5),
    ('last', False, 0.93),
])
@patch('rnalysis.utils.io.PhylomeDBOrthologMapper')
def test_map_orthologs_phylomedb(mock_mapper, non_unique_mode, filter_consistency_score, remove_unmapped, threshold,
                                 ranked_set):
    if remove_unmapped:
        filter_obj_path = 'tests/test_files/test_map_orthologs_remove_unmapped_truth.csv'
    else:
        filter_obj_path = 'tests/test_files/test_map_orthologs_truth.csv'
    filter_obj_truth = Filter(filter_obj_path)
    one2many_truth = pd.DataFrame(
        {'gene': ['gene1', 'gene1', 'gene2'], 'ortholog': ['ortholog1', 'ortholog2', 'ortholog3']})
    mock_mapper_instance = Mock()
    mock_mapper.return_value = mock_mapper_instance
    mock_mapper_instance.get_orthologs.return_value = (
        io.OrthologDict({'gene1': 'ortholog2', 'gene2': 'ortholog3'})
        , io.OrthologDict({'gene1': ['ortholog1', 'ortholog2'], 'gene2': ['ortholog3']}))

    if ranked_set:
        featureset = RankedSet(Filter('tests/test_files/test_map_orthologs.csv'))
    else:
        featureset = FeatureSet(Filter('tests/test_files/test_map_orthologs.csv'))
    one2one, one2many = featureset.map_orthologs_phylomedb('Homo sapiens', 'caenorhabditis elegans',
                                                           gene_id_type='UniProtKB', inplace=False,
                                                           non_unique_mode=non_unique_mode,
                                                           filter_consistency_score=filter_consistency_score,
                                                           consistency_score_threshold=threshold,
                                                           remove_unmapped_genes=remove_unmapped)

    if ranked_set:
        one2one_truth = RankedSet(filter_obj_truth, one2one.set_name)
    else:
        one2one_truth = FeatureSet(filter_obj_truth, one2one.set_name)

    assert one2one == one2one_truth
    assert one2many.equals(one2many_truth)

@pytest.mark.parametrize('ranked_set', [True, False])
@pytest.mark.parametrize('remove_unmapped', [True, False])
@pytest.mark.parametrize('non_unique_mode', ['first', 'last'])
@patch('rnalysis.utils.io.OrthoInspectorOrthologMapper')
def test_map_orthologs_orthoinspector(mock_mapper, non_unique_mode, remove_unmapped, ranked_set):
    if remove_unmapped:
        filter_obj_path = 'tests/test_files/test_map_orthologs_remove_unmapped_truth.csv'
    else:
        filter_obj_path = 'tests/test_files/test_map_orthologs_truth.csv'
    filter_obj_truth = Filter(filter_obj_path)
    one2many_truth = pd.DataFrame(
        {'gene': ['gene1', 'gene1', 'gene2'], 'ortholog': ['ortholog1', 'ortholog2', 'ortholog3']})
    mock_mapper_instance = Mock()
    mock_mapper.return_value = mock_mapper_instance
    mock_mapper_instance.get_orthologs.return_value = (
        io.OrthologDict({'gene1': 'ortholog2', 'gene2': 'ortholog3'})
        , io.OrthologDict({'gene1': ['ortholog1', 'ortholog2'], 'gene2': ['ortholog3']}))

    if ranked_set:
        featureset = RankedSet(Filter('tests/test_files/test_map_orthologs.csv'))
    else:
        featureset = FeatureSet(Filter('tests/test_files/test_map_orthologs.csv'))
    one2one, one2many = featureset.map_orthologs_orthoinspector('Homo sapiens', 'caenorhabditis elegans',
                                                                gene_id_type='UniProtKB', inplace=False,
                                                                non_unique_mode=non_unique_mode,
                                                                remove_unmapped_genes=remove_unmapped)

    if ranked_set:
        one2one_truth = RankedSet(filter_obj_truth, one2one.set_name)
    else:
        one2one_truth = FeatureSet(filter_obj_truth, one2one.set_name)

    assert one2one == one2one_truth
    assert one2many.equals(one2many_truth)
