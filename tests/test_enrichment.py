import pytest
import matplotlib
from rnalysis import filtering
from rnalysis.general import start_parallel_session
from rnalysis.enrichment import *
from tests import __attr_ref__, __biotype_ref__

start_parallel_session()
matplotlib.use('Agg')

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
    up = FeatureSet(up_feature_set)


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


def test_biotypes():
    truth = utils.load_csv('test_files/biotypes_truth.csv', 0)
    genes = {'WBGene00048865', 'WBGene00000106', 'WBGene00000137', 'WBGene00199484', 'WBGene00268190', 'WBGene00048864',
             'WBGene00268189', 'WBGene00268195', 'WBGene00255734', 'WBGene00199485', 'WBGene00048863', 'WBGene00000019',
             'WBGene00268191', 'WBGene00000041', 'WBGene00199486', 'WBGene00255735', 'WBGene00000105',
             'index_that_is_not_in_biotype_ref_table'}

    en = FeatureSet(genes)
    df = en.biotypes(ref=__biotype_ref__)
    df.sort_index(inplace=True)
    truth.sort_index(inplace=True)
    assert np.all(df == truth)


def _enrichment_get_ref_tests_setup(truth, bg_genes):
    genes = {'WBGene00000041', 'WBGene00002074', 'WBGene00000019', 'WBGene00000105', 'WBGene00000106', 'WBGene00199484',
             'WBGene00001436', 'WBGene00000137', 'WBGene00001996', 'WBGene00014208'}
    en = FeatureSet(gene_set=genes, set_name='test_set')

    if not isinstance(bg_genes, str):
        bg_en = FeatureSet(bg_genes, 'background genes')
        res, _ = en._enrichment_get_reference(biotype='all', background_genes=bg_en,
                                              attr_ref_path=__attr_ref__,
                                              biotype_ref_path=__biotype_ref__)
    else:
        res, _ = en._enrichment_get_reference(biotype=bg_genes, background_genes=None,
                                              attr_ref_path=__attr_ref__,
                                              biotype_ref_path=__biotype_ref__)
    truth.sort_index(inplace=True)
    res.sort_index(inplace=True)
    assert np.all(res.index == truth.index)
    assert np.all(res.columns == truth.columns)
    assert np.all(res.attribute1.isna() == truth.attribute1.isna())
    assert np.all(res.attribute2.isna() == truth.attribute2.isna())


def test_enrichment_get_ref_biotype():
    truth = utils.load_csv('test_files/attr_ref_table_for_tests_biotype.csv', 0)
    bg_genes = 'protein_coding'
    _enrichment_get_ref_tests_setup(truth, bg_genes)


def test_enrichment_get_ref_custom_background():
    truth = utils.load_csv('test_files/attr_ref_table_for_tests_specified_bg.csv', 0)
    bg_genes = {'WBGene00003902', 'WBGene00000106', 'WBGene00001436', 'WBGene00000864', 'WBGene00011910',
                'WBGene00000859', 'WBGene00268189', 'WBGene00000865', 'WBGene00003864', 'WBGene00048863',
                'WBGene00000369', 'WBGene00000863', 'WBGene00002074', 'WBGene00000041', 'WBGene00199486',
                'WBGene00000105', 'WBGene00001131'}
    _enrichment_get_ref_tests_setup(truth, bg_genes)


def test_enrichment_get_ref_custom_background_from_featureset_object():
    truth = utils.load_csv('test_files/attr_ref_table_for_tests_specified_bg.csv', 0)
    bg_genes = {'WBGene00003902', 'WBGene00000106', 'WBGene00001436', 'WBGene00000864', 'WBGene00011910',
                'WBGene00000859', 'WBGene00268189', 'WBGene00000865', 'WBGene00003864', 'WBGene00048863',
                'WBGene00000369', 'WBGene00000863', 'WBGene00002074', 'WBGene00000041', 'WBGene00199486',
                'WBGene00000105', 'WBGene00001131'}
    _enrichment_get_ref_tests_setup(truth, bg_genes)


def test_enrichment_get_ref_custom_background_from_filter_object():
    truth = utils.load_csv('test_files/attr_ref_table_for_tests_specified_bg.csv', 0)
    bg_genes = filtering.CountFilter(r'test_files/test_bg_genes_from_filter_object.csv')
    _enrichment_get_ref_tests_setup(truth, bg_genes)


def test_enrichment_get_scales():
    assert False


def tests_enrichment_randomization_api():
    genes = {'WBGene00048865', 'WBGene00000864', 'WBGene00000105', 'WBGene00001996', 'WBGene00011910', 'WBGene00268195',
             'WBGene00255734', 'WBGene00048863', 'WBGene00000369', 'WBGene00000863', 'WBGene00000041', 'WBGene00268190',
             'WBGene00199486', 'WBGene00001131', 'WBGene00003902', 'WBGene00001436', 'WBGene00000865', 'WBGene00001132',
             'WBGene00003864', 'WBGene00000019', 'WBGene00014208', 'WBGene00002074', 'WBGene00000106', 'WBGene00000137',
             'WBGene00000859', 'WBGene00268189'}
    attrs = ['attribute1', 'attribute2']
    en = FeatureSet(gene_set=genes, set_name='test_set')
    _ = en.enrich_randomization(attrs, reps=1, biotype='all', attr_ref_path=__attr_ref__,
                                biotype_ref_path=__biotype_ref__)
    plt.close('all')


def test_enrichment_randomization_reliability():
    genes = {'WBGene00000041', 'WBGene00002074', 'WBGene00000019', 'WBGene00000105', 'WBGene00000106', 'WBGene00199484',
             'WBGene00001436', 'WBGene00000137', 'WBGene00001996', 'WBGene00014208'}
    attrs = ['attribute1', 'attribute2', 'attribute4']
    en = FeatureSet(gene_set=genes, set_name='test_set')
    random_seed = 0

    for i in range(5):
        res1 = en.enrich_randomization(attrs, reps=5000, biotype='all',
                                       attr_ref_path=__attr_ref__,
                                       biotype_ref_path=__biotype_ref__,
                                       random_seed=random_seed)
        res2 = en.enrich_randomization(attrs, reps=5000, biotype='all',
                                       attr_ref_path=__attr_ref__,
                                       biotype_ref_path=__biotype_ref__,
                                       random_seed=random_seed + 1)
        res3 = en.enrich_randomization(attrs, reps=5000, biotype='all',
                                       attr_ref_path=__attr_ref__,
                                       biotype_ref_path=__biotype_ref__,
                                       random_seed=random_seed + 2)
        random_seed += 3
        plt.close('all')
        for col in ['samples', 'obs', 'exp', 'log2_fold_enrichment']:
            assert np.all(res1[col] == res2[col])
            assert np.all(res2[col] == res3[col])
        for randcol in ['pval', 'padj']:
            assert np.isclose(res1[randcol], res2[randcol], atol=4 * 10 ** -4, rtol=0.2).all()
            assert np.isclose(res2[randcol], res3[randcol], atol=4 * 10 ** -4, rtol=0.2).all()
            assert np.isclose(res2[randcol], res1[randcol], atol=4 * 10 ** -4, rtol=0.2).all()
            assert np.isclose(res3[randcol], res2[randcol], atol=4 * 10 ** -4, rtol=0.2).all()


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
    truth = utils.load_csv('test_files/enrichment_randomization_res.csv', 0)
    genes = {'WBGene00000041', 'WBGene00002074', 'WBGene00000105', 'WBGene00000106', 'WBGene00199484',
             'WBGene00001436', 'WBGene00000137', 'WBGene00001996', 'WBGene00014208', 'WBGene00001133'}
    attrs = ['attribute1', 'attribute2']
    en = FeatureSet(gene_set=genes, set_name='test_set')
    res = en.enrich_randomization(attrs, reps=100000, biotype='all',
                                  attr_ref_path=__attr_ref__,
                                  biotype_ref_path=__biotype_ref__, random_seed=0)
    plt.close('all')
    _enrichment_validity(res, truth)


def test_enrichment_randomization_parallel_api():
    genes = {'WBGene00048865', 'WBGene00000864', 'WBGene00000105', 'WBGene00001996', 'WBGene00011910', 'WBGene00268195',
             'WBGene00255734', 'WBGene00048863', 'WBGene00000369', 'WBGene00000863', 'WBGene00000041', 'WBGene00268190',
             'WBGene00199486', 'WBGene00001131', 'WBGene00003902', 'WBGene00001436', 'WBGene00000865', 'WBGene00001132',
             'WBGene00003864', 'WBGene00000019', 'WBGene00014208', 'WBGene00002074', 'WBGene00000106', 'WBGene00000137',
             'WBGene00000859', 'WBGene00268189'}
    attrs = ['attribute1', 'attribute2']
    en = FeatureSet(gene_set=genes, set_name='test_set')
    _ = en.enrich_randomization_parallel(attrs, reps=1, biotype='all',
                                         attr_ref_path=__attr_ref__,
                                         biotype_ref_path=__biotype_ref__)
    plt.close('all')


def test_enrichment_randomization_parallel_reliability():
    genes = {'WBGene00000041', 'WBGene00002074', 'WBGene00000019', 'WBGene00000105', 'WBGene00000106', 'WBGene00199484',
             'WBGene00001436', 'WBGene00000137', 'WBGene00001996', 'WBGene00014208'}
    attrs = ['attribute1', 'attribute2', 'attribute4']
    en = FeatureSet(gene_set=genes, set_name='test_set')
    random_seed = 0

    for i in range(5):
        res1 = en.enrich_randomization_parallel(attrs, reps=5000, biotype='all',
                                                attr_ref_path=__attr_ref__,
                                                biotype_ref_path=__biotype_ref__,
                                                random_seed=random_seed)
        res2 = en.enrich_randomization_parallel(attrs, reps=5000, biotype='all',
                                                attr_ref_path=__attr_ref__,
                                                biotype_ref_path=__biotype_ref__,
                                                random_seed=random_seed + 1)
        res3 = en.enrich_randomization_parallel(attrs, reps=5000, biotype='all',
                                                attr_ref_path=__attr_ref__,
                                                biotype_ref_path=__biotype_ref__,
                                                random_seed=random_seed + 2)
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
    truth = utils.load_csv('test_files/enrichment_randomization_res.csv', 0)
    genes = {'WBGene00000041', 'WBGene00002074', 'WBGene00000105', 'WBGene00000106', 'WBGene00199484',
             'WBGene00001436', 'WBGene00000137', 'WBGene00001996', 'WBGene00014208', 'WBGene00001133'}
    attrs = ['attribute1', 'attribute2']
    en = FeatureSet(gene_set=genes, set_name='test_set')
    res = en.enrich_randomization_parallel(attrs, reps=100000, biotype='all',
                                           attr_ref_path=__attr_ref__,
                                           biotype_ref_path=__biotype_ref__, random_seed=0)
    plt.close('all')
    _enrichment_validity(res, truth)


def test_randomization_int_index_attributes():
    genes = {'WBGene00000041', 'WBGene00002074', 'WBGene00000105', 'WBGene00000106', 'WBGene00199484',
             'WBGene00001436', 'WBGene00000137', 'WBGene00001996', 'WBGene00014208', 'WBGene00001133'}
    en = FeatureSet(gene_set=genes, set_name='test_set')
    attrs_truth = ['attribute1', 'attribute3', 'attribute4']
    attrs = en._enrichment_get_attrs([0, 2, 3], 'test_files/attr_ref_table_for_tests.csv')
    assert attrs == attrs_truth

    attr_truth_single = ['attribute4']
    attr = en._enrichment_get_attrs(3, 'test_files/attr_ref_table_for_tests.csv')
    assert attr == attr_truth_single


def test_randomization_all_attributes():
    genes = {'WBGene00000041', 'WBGene00002074', 'WBGene00000105', 'WBGene00000106', 'WBGene00199484',
             'WBGene00001436', 'WBGene00000137', 'WBGene00001996', 'WBGene00014208', 'WBGene00001133'}
    en = FeatureSet(gene_set=genes, set_name='test_set')
    attrs_truth = ['attribute1', 'attribute2', 'attribute3', 'attribute4']
    attrs = en._enrichment_get_attrs('all', 'test_files/attr_ref_table_for_tests.csv')
    plt.close('all')
    assert attrs == attrs_truth


def test_enrich_statistic_hypergeometric_api():
    genes = {'WBGene00048865', 'WBGene00000864', 'WBGene00000105', 'WBGene00001996', 'WBGene00011910', 'WBGene00268195',
             'WBGene00255734', 'WBGene00048863', 'WBGene00000369', 'WBGene00000863', 'WBGene00000041', 'WBGene00268190',
             'WBGene00199486', 'WBGene00001131', 'WBGene00003902', 'WBGene00001436', 'WBGene00000865', 'WBGene00001132',
             'WBGene00003864', 'WBGene00000019', 'WBGene00014208', 'WBGene00002074', 'WBGene00000106', 'WBGene00000137',
             'WBGene00000859', 'WBGene00268189'}
    attrs = ['attribute1', 'attribute2']
    en = FeatureSet(gene_set=genes, set_name='test_set')
    _ = en.enrich_statistic(attrs, biotype='all', data_scale='boolean', attr_ref_path=__attr_ref__,
                            biotype_ref_path=__biotype_ref__)
    plt.close('all')


def test_enrich_hypergeom_pvalues():
    assert False


def test_calc_hypergeom_pvalues():
    [M, n, N, X] = [13588, 59, 611, 19]
    truth = 4.989682834519698 * 10 ** -12
    pval = FeatureSet._calc_hypergeometric_pval(M, n, N, X)
    assert np.isclose(truth, pval, atol=0, rtol=0.00001)

    [M, n, N, X] = [20000, 430, 700, 6]
    truth = 0.006249179131697138
    pval = FeatureSet._calc_hypergeometric_pval(M, n, N, X)
    assert np.isclose(truth, pval, atol=0, rtol=0.00001)

    [M, n, N, X] = [20000, 43, 300, 3]
    truth = 0.0265186938062861
    pval = FeatureSet._calc_hypergeometric_pval(M, n, N, X)
    assert np.isclose(truth, pval, atol=0, rtol=0.00001)
