import pytest
from rnalysis import general

general.start_parallel_session()
import matplotlib
from rnalysis.enrichment import *

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


def test_enrichment_processing_api():
    up = FeatureSet(up_feature_set)


def test_enrichment_processing_union():
    other = {'WBGene00017419', 'WBGene00016520', 'WBGene00017225', 'WBGene00044200', 'WBGene00206390',
             'WBGene00022523', 'WBGene00000001', 'WBGene00000002'}
    truth = up_feature_set.union(other)
    up = FeatureSet(up_feature_set)
    other_ep = FeatureSet(other)
    up.union(other_ep)
    assert np.all(up.gene_set == truth)


def test_enrichment_processing_intersection():
    other = {'WBGene00017419', 'WBGene00016520', 'WBGene00017225', 'WBGene00044200', 'WBGene00206390',
             'WBGene00022523', 'WBGene00000001', 'WBGene00000002'}
    truth = {'WBGene00017419', 'WBGene00016520', 'WBGene00017225', 'WBGene00044200', 'WBGene00206390',
             'WBGene00022523'}
    up = FeatureSet(up_feature_set)
    other_ep = FeatureSet(other)
    up.intersection(other_ep)
    assert np.all(up.gene_set == truth)


def test_enrichment_processing_difference():
    other = {'WBGene00017419', 'WBGene00016520', 'WBGene00017225', 'WBGene00044200', 'WBGene00206390',
             'WBGene00022523', 'WBGene00000001', 'WBGene00000002'}
    truth = {'WBGene00000001', 'WBGene00000002'}
    up = FeatureSet(up_feature_set)
    other_ep = FeatureSet(other)
    other_ep.difference(up)
    assert np.all(other_ep.gene_set == truth)


def test_enrichment_processing_symmetric_difference():
    first = {'WBGene00016520', 'WBGene00017225', 'WBGene00044200', 'WBGene00206390'}
    second = {'WBGene00044200', 'WBGene00206390',
              'WBGene00022523', 'WBGene00000001', 'WBGene00000002'}
    truth = {'WBGene00016520', 'WBGene00017225', 'WBGene00022523', 'WBGene00000001', 'WBGene00000002'}
    first_ep = FeatureSet(first)
    second_ep = FeatureSet(second)
    direction1 = second_ep.symmetric_difference(first_ep, inplace=False)
    direction2 = first_ep.symmetric_difference(second_ep, inplace=False)
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
    first_ep = FeatureSet(first)
    symm_diff = first_ep.symmetric_difference(second, inplace=False)
    assert np.all(symm_diff.gene_set == truth)


def test_biotypes():
    truth = general.load_csv('biotypes_truth.csv', 0)
    genes = {'WBGene00048865', 'WBGene00000106', 'WBGene00000137', 'WBGene00199484', 'WBGene00268190', 'WBGene00048864',
             'WBGene00268189', 'WBGene00268195', 'WBGene00255734', 'WBGene00199485', 'WBGene00048863', 'WBGene00000019',
             'WBGene00268191', 'WBGene00000041', 'WBGene00199486', 'WBGene00255735', 'WBGene00000105',
             'index_that_is_not_in_biotype_ref_table'}

    en = FeatureSet(genes)
    df = en.biotypes(ref='biotype_ref_table_for_tests.csv')
    df.sort_index(inplace=True)
    truth.sort_index(inplace=True)
    assert np.all(df == truth)


def test_enrichment_get_ref_biotype():
    truth = general.load_csv('attr_ref_table_for_tests_biotype.csv', 0)
    genes = {'WBGene00000041', 'WBGene00002074', 'WBGene00000019', 'WBGene00000105', 'WBGene00000106', 'WBGene00199484',
             'WBGene00001436', 'WBGene00000137', 'WBGene00001996', 'WBGene00014208'}
    en = FeatureSet(gene_set=genes, set_name='test_set')
    res, _ = en._enrichment_get_reference(biotype='protein_coding', background_genes=None,
                                          attr_ref_path='attr_ref_table_for_tests.csv',
                                          biotype_ref_path='biotype_ref_table_for_tests.csv')
    truth.sort_index(inplace=True)
    res.sort_index(inplace=True)

    assert np.all(res.index == truth.index)
    assert np.all(res.columns == truth.columns)
    assert np.all(res.attribute1.isna() == truth.attribute1.isna())
    assert np.all(res.attribute2.isna() == truth.attribute2.isna())
    print(res.int_index)
    print(truth.int_index)
    assert np.all(res.int_index == truth.int_index)


def test_enrichment_get_ref_custom_background():
    truth = general.load_csv('attr_ref_table_for_tests_specified_bg.csv', 0)
    bg_genes = {'WBGene00003902', 'WBGene00000106', 'WBGene00001436', 'WBGene00000864', 'WBGene00011910',
                'WBGene00000859', 'WBGene00268189', 'WBGene00000865', 'WBGene00003864', 'WBGene00048863',
                'WBGene00000369', 'WBGene00000863', 'WBGene00002074', 'WBGene00000041', 'WBGene00199486',
                'WBGene00000105', 'WBGene00001131'}
    genes = {'WBGene00000041', 'WBGene00002074', 'WBGene00000019', 'WBGene00000105', 'WBGene00000106', 'WBGene00199484',
             'WBGene00001436', 'WBGene00000137', 'WBGene00001996', 'WBGene00014208'}
    en = FeatureSet(gene_set=genes, set_name='test_set')

    res, _ = en._enrichment_get_reference(biotype='all', background_genes=bg_genes,
                                          attr_ref_path='attr_ref_table_for_tests.csv',
                                          biotype_ref_path='biotype_ref_table_for_tests.csv')
    truth.sort_index(inplace=True)
    res.sort_index(inplace=True)

    assert np.all(res.index == truth.index)
    assert np.all(res.columns == truth.columns)
    assert np.all(res.attribute1.isna() == truth.attribute1.isna())
    assert np.all(res.attribute2.isna() == truth.attribute2.isna())
    assert np.all(res.int_index == truth.int_index)


def test_enrichment_get_ref_custom_background_from_featureset_object():
    truth = general.load_csv('attr_ref_table_for_tests_specified_bg.csv', 0)
    bg_genes = {'WBGene00003902', 'WBGene00000106', 'WBGene00001436', 'WBGene00000864', 'WBGene00011910',
                'WBGene00000859', 'WBGene00268189', 'WBGene00000865', 'WBGene00003864', 'WBGene00048863',
                'WBGene00000369', 'WBGene00000863', 'WBGene00002074', 'WBGene00000041', 'WBGene00199486',
                'WBGene00000105', 'WBGene00001131'}
    genes = {'WBGene00000041', 'WBGene00002074', 'WBGene00000019', 'WBGene00000105', 'WBGene00000106', 'WBGene00199484',
             'WBGene00001436', 'WBGene00000137', 'WBGene00001996', 'WBGene00014208'}
    en = FeatureSet(gene_set=genes, set_name='test_set')
    bg_en = FeatureSet(bg_genes, 'background genes')
    res, _ = en._enrichment_get_reference(biotype='all', background_genes=bg_en,
                                          attr_ref_path='attr_ref_table_for_tests.csv',
                                          biotype_ref_path='biotype_ref_table_for_tests.csv')
    truth.sort_index(inplace=True)
    res.sort_index(inplace=True)


def test_enrichment_get_ref_custom_background_from_filter_object():
    truth = general.load_csv('attr_ref_table_for_tests_specified_bg.csv', 0)
    bg_genes = filtering.CountFilter(r'test_bg_genes_from_filter_object.csv')
    genes = {'WBGene00000041', 'WBGene00002074', 'WBGene00000019', 'WBGene00000105', 'WBGene00000106', 'WBGene00199484',
             'WBGene00001436', 'WBGene00000137', 'WBGene00001996', 'WBGene00014208'}
    en = FeatureSet(gene_set=genes, set_name='test_set')

    res, _ = en._enrichment_get_reference(biotype='all', background_genes=bg_genes,
                                          attr_ref_path='attr_ref_table_for_tests.csv',
                                          biotype_ref_path='biotype_ref_table_for_tests.csv')
    truth.sort_index(inplace=True)
    res.sort_index(inplace=True)

    assert np.all(res.index == truth.index)
    assert np.all(res.columns == truth.columns)
    assert np.all(res.attribute1.isna() == truth.attribute1.isna())
    assert np.all(res.attribute2.isna() == truth.attribute2.isna())
    assert np.all(res.int_index == truth.int_index)

    assert np.all(res.index == truth.index)
    assert np.all(res.columns == truth.columns)
    assert np.all(res.attribute1.isna() == truth.attribute1.isna())
    assert np.all(res.attribute2.isna() == truth.attribute2.isna())
    assert np.all(res.int_index == truth.int_index)


def tests_enrichment_randomization_api():
    genes = {'WBGene00048865', 'WBGene00000864', 'WBGene00000105', 'WBGene00001996', 'WBGene00011910', 'WBGene00268195',
             'WBGene00255734', 'WBGene00048863', 'WBGene00000369', 'WBGene00000863', 'WBGene00000041', 'WBGene00268190',
             'WBGene00199486', 'WBGene00001131', 'WBGene00003902', 'WBGene00001436', 'WBGene00000865', 'WBGene00001132',
             'WBGene00003864', 'WBGene00000019', 'WBGene00014208', 'WBGene00002074', 'WBGene00000106', 'WBGene00000137',
             'WBGene00000859', 'WBGene00268189'}
    attrs = ['attribute1', 'attribute2']
    en = FeatureSet(gene_set=genes, set_name='test_set')
    _ = en.enrich_randomization(attrs, reps=1, biotype='all', attr_ref_path='attr_ref_table_for_tests.csv',
                                biotype_ref_path='biotype_ref_table_for_tests.csv')


def test_enrichment_randomization_reliability():
    genes = {'WBGene00000041', 'WBGene00002074', 'WBGene00000019', 'WBGene00000105', 'WBGene00000106', 'WBGene00199484',
             'WBGene00001436', 'WBGene00000137', 'WBGene00001996', 'WBGene00014208'}
    attrs = ['attribute1', 'attribute2', 'attribute4']
    en = FeatureSet(gene_set=genes, set_name='test_set')
    random_state = random.getstate()

    for i in range(5):
        res1 = en.enrich_randomization(attrs, reps=5000, biotype='all', attr_ref_path='attr_ref_table_for_tests.csv',
                                       biotype_ref_path='biotype_ref_table_for_tests.csv')
        res2 = en.enrich_randomization(attrs, reps=5000, biotype='all', attr_ref_path='attr_ref_table_for_tests.csv',
                                       biotype_ref_path='biotype_ref_table_for_tests.csv')
        res3 = en.enrich_randomization(attrs, reps=5000, biotype='all', attr_ref_path='attr_ref_table_for_tests.csv',
                                       biotype_ref_path='biotype_ref_table_for_tests.csv')
        try:
            for col in ['samples', 'n obs', 'n exp', 'log2_fold_enrichment']:
                assert np.all(res1[col] == res2[col])
                assert np.all(res2[col] == res3[col])
            for randcol in ['pval', 'padj']:
                assert np.isclose(res1[randcol], res2[randcol], atol=4 * 10 ** -4, rtol=0.2).all()
                assert np.isclose(res2[randcol], res3[randcol], atol=4 * 10 ** -4, rtol=0.2).all()
                assert np.isclose(res2[randcol], res1[randcol], atol=4 * 10 ** -4, rtol=0.2).all()
                assert np.isclose(res3[randcol], res2[randcol], atol=4 * 10 ** -4, rtol=0.2).all()
        except AssertionError:
            raise AssertionError(f'Enrichment test failed with the random state: \n{random_state}')


def test_enrichment_randomization_validity():
    truth = general.load_csv('enrichment_randomization_res.csv', 0)
    genes = {'WBGene00000041', 'WBGene00002074', 'WBGene00000105', 'WBGene00000106', 'WBGene00199484',
             'WBGene00001436', 'WBGene00000137', 'WBGene00001996', 'WBGene00014208', 'WBGene00001133'}
    attrs = ['attribute1', 'attribute2']
    en = FeatureSet(gene_set=genes, set_name='test_set')
    random_state = random.getstate()
    res = en.enrich_randomization(attrs, reps=100000, biotype='all', attr_ref_path='attr_ref_table_for_tests.csv',
                                  biotype_ref_path='biotype_ref_table_for_tests.csv')

    try:
        for col in ['samples', 'n obs', 'significant']:
            assert np.all(res[col] == truth[col])
        for closecol in ['n exp', 'log2_fold_enrichment']:
            assert np.isclose(res[closecol], truth[closecol], atol=0.0).all()
        for randcol in ['pval']:
            assert np.isclose(res[randcol], truth[randcol], atol=2 * 10 ** -4, rtol=0.25).all()
        pvals = res['pval'].values
        _, padj_truth = multitest.fdrcorrection(pvals, 0.1)
        assert np.isclose(res['padj'].values, padj_truth, atol=0.0).all()
    except AssertionError:
        raise AssertionError(f'Enrichment test failed with the random state: \n{random_state}')


def test_enrichment_parallel_api():
    genes = {'WBGene00048865', 'WBGene00000864', 'WBGene00000105', 'WBGene00001996', 'WBGene00011910', 'WBGene00268195',
             'WBGene00255734', 'WBGene00048863', 'WBGene00000369', 'WBGene00000863', 'WBGene00000041', 'WBGene00268190',
             'WBGene00199486', 'WBGene00001131', 'WBGene00003902', 'WBGene00001436', 'WBGene00000865', 'WBGene00001132',
             'WBGene00003864', 'WBGene00000019', 'WBGene00014208', 'WBGene00002074', 'WBGene00000106', 'WBGene00000137',
             'WBGene00000859', 'WBGene00268189'}
    attrs = ['attribute1', 'attribute2']
    en = FeatureSet(gene_set=genes, set_name='test_set')
    _ = en.enrich_randomization_parallel(attrs, reps=1, biotype='all', attr_ref_path='attr_ref_table_for_tests.csv',
                                         biotype_ref_path='biotype_ref_table_for_tests.csv')


def test_enrichment_randomization_parallel_reliability():
    genes = {'WBGene00000041', 'WBGene00002074', 'WBGene00000019', 'WBGene00000105', 'WBGene00000106', 'WBGene00199484',
             'WBGene00001436', 'WBGene00000137', 'WBGene00001996', 'WBGene00014208'}
    attrs = ['attribute1', 'attribute2', 'attribute4']
    en = FeatureSet(gene_set=genes, set_name='test_set')
    random_state = np.random.get_state()

    try:
        for i in range(5):
            res1 = en.enrich_randomization_parallel(attrs, reps=5000, biotype='all',
                                                    attr_ref_path='attr_ref_table_for_tests.csv',
                                                    biotype_ref_path='biotype_ref_table_for_tests.csv')
            res2 = en.enrich_randomization_parallel(attrs, reps=5000, biotype='all',
                                                    attr_ref_path='attr_ref_table_for_tests.csv',
                                                    biotype_ref_path='biotype_ref_table_for_tests.csv')
            res3 = en.enrich_randomization_parallel(attrs, reps=5000, biotype='all',
                                                    attr_ref_path='attr_ref_table_for_tests.csv',
                                                    biotype_ref_path='biotype_ref_table_for_tests.csv')

        for col in ['samples', 'n obs', 'n exp', 'log2_fold_enrichment']:
            assert np.all(res1[col] == res2[col])
            assert np.all(res2[col] == res3[col])
        for randcol in ['pval']:
            assert np.isclose(res1[randcol], res2[randcol], atol=4 * 10 ** -4, rtol=0.25).all()
            assert np.isclose(res2[randcol], res3[randcol], atol=4 * 10 ** -4, rtol=0.25).all()
            assert np.isclose(res2[randcol], res1[randcol], atol=4 * 10 ** -4, rtol=0.25).all()
            assert np.isclose(res3[randcol], res2[randcol], atol=4 * 10 ** -4, rtol=0.25).all()
    except AssertionError:
        raise AssertionError(f'Enrichment test failed with the numpy.random state: \n{random_state}')


def test_enrichment_parallel_validity():
    truth = general.load_csv('enrichment_randomization_res.csv', 0)
    genes = {'WBGene00000041', 'WBGene00002074', 'WBGene00000105', 'WBGene00000106', 'WBGene00199484',
             'WBGene00001436', 'WBGene00000137', 'WBGene00001996', 'WBGene00014208', 'WBGene00001133'}
    attrs = ['attribute1', 'attribute2']
    en = FeatureSet(gene_set=genes, set_name='test_set')
    random_state = np.random.get_state()
    res = en.enrich_randomization_parallel(attrs, reps=100000, biotype='all',
                                           attr_ref_path='attr_ref_table_for_tests.csv',
                                           biotype_ref_path='biotype_ref_table_for_tests.csv')

    try:
        for col in ['samples', 'n obs', 'significant']:
            assert np.all(res[col] == truth[col])
        for closecol in ['n exp', 'log2_fold_enrichment']:
            assert np.isclose(res[closecol], truth[closecol], atol=0.0).all()
        for randcol in ['pval']:
            assert np.isclose(res[randcol], truth[randcol], atol=2 * 10 ** -4, rtol=0.25).all()
        pvals = res['pval'].values
        _, padj_truth = multitest.fdrcorrection(pvals, 0.1)
        assert np.isclose(res['padj'].values, padj_truth, atol=0.0).all()
    except AssertionError:
        raise AssertionError(f'Enrichment test failed with the numpy.random state: \n{random_state}')


def test_randomization_int_index_attributes():
    genes = {'WBGene00000041', 'WBGene00002074', 'WBGene00000105', 'WBGene00000106', 'WBGene00199484',
             'WBGene00001436', 'WBGene00000137', 'WBGene00001996', 'WBGene00014208', 'WBGene00001133'}
    en = FeatureSet(gene_set=genes, set_name='test_set')
    attrs_truth = ['attribute1', 'attribute3', 'attribute4']
    attrs = en._enrichment_get_attrs([0, 2, 3], 'attr_ref_table_for_tests.csv')
    assert attrs == attrs_truth

    attr_truth_single = ['attribute4']
    attr = en._enrichment_get_attrs(3, 'attr_ref_table_for_tests.csv')
    assert attr == attr_truth_single


def test_randomization_all_attributes():
    genes = {'WBGene00000041', 'WBGene00002074', 'WBGene00000105', 'WBGene00000106', 'WBGene00199484',
             'WBGene00001436', 'WBGene00000137', 'WBGene00001996', 'WBGene00014208', 'WBGene00001133'}
    en = FeatureSet(gene_set=genes, set_name='test_set')
    attrs_truth = ['attribute1', 'attribute2', 'attribute3', 'attribute4']
    attrs = en._enrichment_get_attrs('all', 'attr_ref_table_for_tests.csv')
    assert attrs == attrs_truth
