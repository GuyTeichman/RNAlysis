import pytest
from collections import namedtuple

import statsmodels.stats.multitest as multitest

from rnalysis.utils.ontology import *
from rnalysis import filtering
from rnalysis.utils.enrichment_runner import *
from rnalysis.utils.io import *
from tests import __attr_ref__, __biotype_ref__


def test_get_pval_asterisk():
    assert EnrichmentRunner._get_pval_asterisk(0.6) == ('ns', 'normal')
    assert EnrichmentRunner._get_pval_asterisk(0.001, 0.00099) == ('ns', 'normal')
    assert EnrichmentRunner._get_pval_asterisk(0.04) == (u'\u2217', 'bold')
    assert EnrichmentRunner._get_pval_asterisk(0.0099) == (u'\u2217' * 2, 'bold')
    assert EnrichmentRunner._get_pval_asterisk(0) == (u'\u2217' * 4, 'bold')


def test_calc_randomization_pval():
    np.random.seed(42)
    hypergeom_pval = 0.2426153598589023
    avg_pval = 0
    for i in range(5):
        avg_pval += EnrichmentRunner._calc_randomization_pval(500, 1, np.random.random(10000) > 0.1, 100000, 0.11)
    avg_pval /= 5
    assert np.isclose(avg_pval, hypergeom_pval, atol=0.02)


def test_calc_hypergeometric_pvalues():
    [M, n, N, X] = [13588, 59, 611, 19]
    truth = 4.989682834519698 * 10 ** -12
    pval = EnrichmentRunner._calc_hypergeometric_pval(M, n, N, X)
    assert np.isclose(truth, pval, atol=0, rtol=0.00001)

    [M, n, N, X] = [20000, 430, 700, 6]
    truth = 0.006249179131697138
    pval = EnrichmentRunner._calc_hypergeometric_pval(M, n, N, X)
    assert np.isclose(truth, pval, atol=0, rtol=0.00001)

    [M, n, N, X] = [20000, 43, 300, 3]
    truth = 0.0265186938062861
    pval = EnrichmentRunner._calc_hypergeometric_pval(M, n, N, X)
    assert np.isclose(truth, pval, atol=0, rtol=0.00001)


def test_enrichment_get_attrs_int_index_attributes():
    genes = {'WBGene00000041', 'WBGene00002074', 'WBGene00000105', 'WBGene00000106', 'WBGene00199484',
             'WBGene00001436', 'WBGene00000137', 'WBGene00001996', 'WBGene00014208', 'WBGene00001133'}
    e = EnrichmentRunner(genes, [0, 2, 3], 0.05, __attr_ref__, False, '', False, True, 'test_set', False,
                         'hypergeometric', 'all', None, __biotype_ref__)
    e.fetch_annotations()
    e.fetch_attributes()
    e.get_background_set()
    e.update_gene_set()
    e.filter_annotations()
    attrs = e.attributes
    attrs_truth = ['attribute1', 'attribute3', 'attribute4']
    assert attrs == attrs_truth

    e = EnrichmentRunner(genes, 1, 0.05, __attr_ref__, False, '', False, True, 'test_set', False, 'hypergeometric',
                         'all', None, __biotype_ref__)
    e.fetch_annotations()
    e.fetch_attributes()
    e.get_background_set()
    e.update_gene_set()
    e.filter_annotations()
    attr = e.attributes
    attr_truth_single = ['attribute2']
    assert attr == attr_truth_single


def test_enrichment_get_attrs_all_attributes():
    genes = {'WBGene00000041', 'WBGene00002074', 'WBGene00000105', 'WBGene00000106', 'WBGene00199484',
             'WBGene00001436', 'WBGene00000137', 'WBGene00001996', 'WBGene00014208', 'WBGene00001133'}
    e = EnrichmentRunner(genes, 'all', 0.05, __attr_ref__, False, '', False, True, 'test_set', False, 'hypergeometric',
                         'all', None, __biotype_ref__)
    e.fetch_annotations()
    e.fetch_attributes()
    e.get_background_set()
    e.update_gene_set()
    e.filter_annotations()
    attrs = e.attributes
    attrs_truth = ['attribute1', 'attribute2', 'attribute3', 'attribute4']
    plt.close('all')
    assert attrs == attrs_truth


def test_enrichment_get_attrs_from_string(monkeypatch):
    monkeypatch.setattr('builtins.input', lambda x: 'attribute1\nattribute4\n')
    genes = {'WBGene00000041', 'WBGene00002074', 'WBGene00000105', 'WBGene00000106', 'WBGene00199484',
             'WBGene00001436', 'WBGene00000137', 'WBGene00001996', 'WBGene00014208', 'WBGene00001133'}
    e = EnrichmentRunner(genes, None, 0.05, __attr_ref__, False, '', False, True, 'test_set', False, 'hypergeometric',
                         'all', None, __biotype_ref__)
    e.fetch_annotations()
    e.fetch_attributes()
    e.get_background_set()
    e.update_gene_set()
    e.filter_annotations()
    attrs = e.attributes
    attrs_truth = ['attribute1', 'attribute4']
    plt.close('all')
    assert attrs == attrs_truth


def test_enrichment_get_attrs_bad_path():
    e = EnrichmentRunner({'_'}, 'attribute1', 0.05, 'fakepath', False, '', False, True, 'test_set', False,
                         'hypergeometric', biotypes='all', biotype_ref_path=__biotype_ref__)
    with pytest.raises(FileNotFoundError):
        e.fetch_annotations()
        e.fetch_attributes()
        e.get_background_set()
        e.update_gene_set()
        e.filter_annotations()


def _enrichment_get_ref_tests_setup(truth, bg_genes):
    genes = {'WBGene00000041', 'WBGene00002074', 'WBGene00000019', 'WBGene00000105', 'WBGene00000106', 'WBGene00199484',
             'WBGene00001436', 'WBGene00000137', 'WBGene00001996', 'WBGene00014208'}
    biotype = bg_genes if isinstance(bg_genes, str) else 'all'
    background = None if isinstance(bg_genes, str) else bg_genes
    e = EnrichmentRunner(genes, 'all', 0.05, __attr_ref__, False, '', False, True, 'test_set', False, 'hypergeometric',
                         biotype, background, __biotype_ref__)
    e.fetch_annotations()
    e.fetch_attributes()
    e.get_background_set()
    e.update_gene_set()
    e.filter_annotations()
    res = e.annotation_df

    truth.sort_index(inplace=True)
    res.sort_index(inplace=True)
    assert np.all(res.index == truth.index)
    assert np.all(res.columns == truth.columns)
    assert np.all(res.attribute1.isna() == truth.attribute1.isna())
    assert np.all(res.attribute2.isna() == truth.attribute2.isna())


def test_enrichment_get_ref_biotype():
    truth = io.load_csv('tests/test_files/attr_ref_table_for_tests_biotype.csv', 0)
    bg_genes = 'protein_coding'
    _enrichment_get_ref_tests_setup(truth, bg_genes)


def test_enrichment_get_ref_custom_background():
    truth = io.load_csv('tests/test_files/attr_ref_table_for_tests_specified_bg.csv', 0)
    bg_genes = {'WBGene00003902', 'WBGene00000106', 'WBGene00001436', 'WBGene00000864', 'WBGene00011910',
                'WBGene00000859', 'WBGene00268189', 'WBGene00000865', 'WBGene00003864', 'WBGene00048863',
                'WBGene00000369', 'WBGene00000863', 'WBGene00002074', 'WBGene00000041', 'WBGene00199486',
                'WBGene00000105', 'WBGene00001131'}
    _enrichment_get_ref_tests_setup(truth, bg_genes)


def test_enrichment_get_ref_custom_background_from_featureset_object():
    truth = io.load_csv('tests/test_files/attr_ref_table_for_tests_specified_bg.csv', 0)
    bg_genes = {'WBGene00003902', 'WBGene00000106', 'WBGene00001436', 'WBGene00000864', 'WBGene00011910',
                'WBGene00000859', 'WBGene00268189', 'WBGene00000865', 'WBGene00003864', 'WBGene00048863',
                'WBGene00000369', 'WBGene00000863', 'WBGene00002074', 'WBGene00000041', 'WBGene00199486',
                'WBGene00000105', 'WBGene00001131'}
    _enrichment_get_ref_tests_setup(truth, bg_genes)


def test_enrichment_get_ref_custom_background_from_filter_object():
    truth = io.load_csv('tests/test_files/attr_ref_table_for_tests_specified_bg.csv', 0)
    bg_genes = filtering.CountFilter(r'tests/test_files/test_bg_genes_from_filter_object.csv')
    _enrichment_get_ref_tests_setup(truth, bg_genes)


def test_results_to_csv():
    try:
        en = EnrichmentRunner({''}, 'all', 0.05, __attr_ref__, True, 'tests/test_files/tmp_enrichment_csv.csv', False,
                              True, 'test_set', False, 'hypergeometric', 'all', None, __biotype_ref__)
        df = pd.read_csv('tests/test_files/enrichment_hypergeometric_res.csv', index_col=0)
        en.results = df
        en.results_to_csv()
        df_loaded = pd.read_csv('tests/test_files/tmp_enrichment_csv.csv', index_col=0)
        print('\n')
        print(df)
        print(df_loaded)
        assert df.equals(df_loaded)
    except Exception as e:
        raise e
    finally:
        try:
            os.remove('tests/test_files/tmp_enrichment_csv.csv')
        except:
            pass


def _comp_go_res_df(res, truth):
    res.drop('name', axis=1, inplace=True)
    res.rename_axis('go_id')
    assert res.loc[:, ['n', 'obs']].equals(truth.loc[:, ['n', 'obs']])
    assert np.allclose(res.loc[:, ['exp', 'log2fc']], res.loc[:, ['exp', 'log2fc']])
    assert np.allclose(res['pval'], truth['pval'], atol=0)


def test_classic_pvals(monkeypatch):
    goa_df = pd.read_csv('tests/test_files/goa_table.csv', index_col=0).astype('bool')
    gene_set = {'gene1', 'gene2', 'gene5', 'gene12', 'gene13', 'gene17', 'gene19', 'gene25', 'gene27', 'gene28'}
    truth = pd.read_csv('tests/test_files/go_pvalues_classic_truth.csv', index_col=0).sort_index()
    dummy_go_node = namedtuple('DummyGONode', field_names='name')

    class DummyDAGTree:
        def __init__(self):
            self.dummy_node = dummy_go_node('content of name field')

        def __getitem__(self, item):
            return self.dummy_node

    monkeypatch.setattr(io, 'fetch_go_basic', lambda: DummyDAGTree())

    e = GOEnrichmentRunner(gene_set, 'elegans', 'WBGene', 0.05, 'classic', 'any', 'any', None, 'any', None, 'any', None,
                           False, False, '', False, False, False, '', False, 'hypergeometric', 'all')
    e.annotation_df = goa_df
    e.mod_annotation_dfs = goa_df,
    e.attributes = list(goa_df.columns)

    res = pd.DataFrame.from_dict(e._go_classic_over_chunk(tuple(goa_df.columns), 0), orient='index',
                                 columns=['name', 'n', 'obs', 'exp', 'log2fc', 'pval']).sort_index()
    _comp_go_res_df(res, truth)


def test_elim_pvals(monkeypatch):
    goa_df = pd.read_csv('tests/test_files/goa_table.csv', index_col=0).astype('bool')
    threshold = 0.2  # make sure there are both significant and non-significant examples with our small bg size (30)
    gene_set = {'gene1', 'gene2', 'gene5', 'gene12', 'gene13', 'gene17', 'gene19', 'gene25', 'gene27', 'gene28'}
    truth = pd.read_csv('tests/test_files/go_pvalues_elim_truth.csv', index_col=0).sort_index()
    with open('tests/test_files/obo_for_go_tests.obo', 'rb') as f:
        dag_tree = ontology.DAGTree(f, ['is_a'])
    monkeypatch.setattr(io, 'fetch_go_basic', lambda: dag_tree)

    e = GOEnrichmentRunner(gene_set, 'elegans', 'WBGene', threshold, 'classic', 'any', 'any', None, 'any', None, 'any',
                           None, False, False, '', False, False, False, '', False, 'hypergeometric', 'all')
    e.annotation_df = goa_df
    e.mod_annotation_dfs = goa_df.copy(deep=True),
    e.attributes = list(goa_df.columns)
    e.attributes_set = set(e.attributes)

    res = pd.DataFrame.from_dict(e._go_elim_on_aspect('all'), orient='index',
                                 columns=['name', 'n', 'obs', 'exp', 'log2fc', 'pval']).sort_index()

    _comp_go_res_df(res, truth)


def test_weight_pvals(monkeypatch):
    goa_df = pd.read_csv('tests/test_files/goa_table.csv', index_col=0).astype('bool')
    gene_set = {'gene1', 'gene2', 'gene5', 'gene12', 'gene13', 'gene17', 'gene19', 'gene25', 'gene27', 'gene28'}
    truth = pd.read_csv('tests/test_files/go_pvalues_weight_truth.csv', index_col=0).sort_index()
    with open('tests/test_files/obo_for_go_tests.obo', 'rb') as f:
        dag_tree = ontology.DAGTree(f, ['is_a'])

    monkeypatch.setattr(io, 'fetch_go_basic', lambda: dag_tree)

    e = GOEnrichmentRunner(gene_set, 'elegans', 'WBGene', 0.05, 'classic', 'any', 'any', None, 'any', None, 'any',
                           None, False, False, '', False, False, False, '', False, 'hypergeometric', 'all')
    e.annotation_df = goa_df
    e.mod_annotation_dfs = goa_df.copy(deep=True),
    e.attributes = list(goa_df.columns)
    e.attributes_set = set(e.attributes)

    res = pd.DataFrame.from_dict(e._go_weight_on_aspect('all'), orient='index',
                                 columns=['name', 'n', 'obs', 'exp', 'log2fc', 'pval']).sort_index()
    _comp_go_res_df(res, truth)


def test_allm_pvals(monkeypatch):
    goa_df = pd.read_csv('tests/test_files/goa_table.csv', index_col=0).astype('bool')
    threshold = 0.2  # make sure there are both significant and non-significant examples with our small bg size (30)
    gene_set = {'gene1', 'gene2', 'gene5', 'gene12', 'gene13', 'gene17', 'gene19', 'gene25', 'gene27', 'gene28'}
    truth = pd.read_csv('tests/test_files/go_pvalues_allm_truth.csv', index_col=0).sort_index()
    with open('tests/test_files/obo_for_go_tests.obo', 'rb') as f:
        dag_tree = ontology.DAGTree(f, ['is_a'])
    monkeypatch.setattr(io, 'fetch_go_basic', lambda: dag_tree)

    e = GOEnrichmentRunner(gene_set, 'elegans', 'WBGene', threshold, 'classic', 'any', 'any', None, 'any', None, 'any',
                           None, False, False, '', False, False, False, '', False, 'hypergeometric', 'all')
    e.annotation_df = goa_df
    e.mod_annotation_dfs = goa_df.copy(deep=True),
    e.attributes = list(goa_df.columns)
    e.attributes_set = set(e.attributes)

    res = pd.DataFrame.from_dict(e._go_allm_pvalues_parallel(), orient='index',
                                 columns=['name', 'n', 'obs', 'exp', 'log2fc', 'pval']).sort_index()

    _comp_go_res_df(res, truth)
