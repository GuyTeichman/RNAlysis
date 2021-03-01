from collections import namedtuple

import pytest

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


def test_enrichment_runner_update_ranked_genes():
    runner = EnrichmentRunner.__new__(EnrichmentRunner)
    truth = np.array([['111', '000', '222', '777', '333', '888']], dtype='str')
    runner.ranked_genes = np.array(['111', '000', '222', '777', '333', '000', '111', '888', '000'], dtype='str')
    runner.gene_set = {'111', '000', '222', '777', '333', '888'}
    runner._update_ranked_genes()
    assert np.all(truth == runner.ranked_genes)


@pytest.mark.parametrize("test_input,expected", [
    ('fisher', EnrichmentRunner._fisher_enrichment),
    ('HYPERGEOMETRIC', EnrichmentRunner._hypergeometric_enrichment),
    ('Randomization', EnrichmentRunner._randomization_enrichment),
    ('XLmHG', EnrichmentRunner._xlmhg_enrichment)])
def test_enrichment_runner_get_enrichment_func(test_input, expected):
    runner = EnrichmentRunner.__new__(EnrichmentRunner)
    assert runner._get_enrichment_func(test_input).__name__ == expected.__name__
    assert runner._get_enrichment_func(test_input.upper()).__name__ == expected.__name__
    assert runner._get_enrichment_func(test_input.lower()).__name__ == expected.__name__
    assert runner._get_enrichment_func(test_input.capitalize()).__name__ == expected.__name__


@pytest.mark.parametrize("test_input,err", [
    ('fifty', ValueError),
    (50, AssertionError),
    (True, AssertionError),
    (max, AssertionError)])
def test_enrichment_runner_get_enrichment_func_invalid_value(test_input, err):
    runner = EnrichmentRunner.__new__(EnrichmentRunner)
    with pytest.raises(err):
        runner._get_enrichment_func(test_input)


@pytest.mark.parametrize('attr,results',
                         [('attribute1', (38, 6, 11, 3)),
                          ('attribute3', (38, 6, 14, 4)),
                          ('attribute4', (38, 6, 13, 4))])
def test_enrichment_runner_get_hypergeometric_parameters(attr, results):
    runner = EnrichmentRunner.__new__(EnrichmentRunner)
    runner.annotation_df = pd.read_csv('tests/test_files/attr_ref_table_for_tests.csv', index_col=0)
    runner.gene_set = {'WBGene00000019', 'WBGene00000041', 'WBGene00000106', 'WBGene00001133', 'WBGene00003915',
                       'WBGene00268195'}
    assert runner._get_hypergeometric_parameters(attr) == results


@pytest.mark.parametrize('params,truth',
                         [((38, 11, 6, 5), ['attribute', 11, 5, (6 / 38) * 11, np.log2(5 / ((6 / 38) * 11)), 0.05]),
                          ((40, 10, 15, 0), ['attribute', 10, 0, (15 / 40) * 10, -np.inf, 0.05])])
def test_enrichment_runner_hypergeometric_enrichment(monkeypatch, params, truth):
    monkeypatch.setattr(EnrichmentRunner, '_get_hypergeometric_parameters', lambda self, attr: params)

    def alt_calc_pval(self, bg_size, de_size, go_size, go_de_size):
        assert (bg_size, de_size, go_size, go_de_size) == params
        return 0.05

    monkeypatch.setattr(EnrichmentRunner, '_calc_hypergeometric_pval', alt_calc_pval)
    runner = EnrichmentRunner.__new__(EnrichmentRunner)
    assert runner._hypergeometric_enrichment('attribute') == truth


@pytest.mark.parametrize('truth',
                         [(['attribute1', 6, 3, (11 / 38) * 6, np.log2(3 / ((11 / 38) * 6)), 0.05]),
                          (['attribute4', 6, 4, (13 / 38) * 6, np.log2(4 / ((13 / 38) * 6)), 0.05])])
def test_enrichment_runner_randomization_enrichment(monkeypatch, truth):
    reps_truth = 100
    df = pd.read_csv('tests/test_files/attr_ref_table_for_tests.csv', index_col=0)
    bg_array_truth = pd.read_csv('tests/test_files/annotation_df_bg_array_truth.csv', index_col=0)[truth[0]]
    gene_set_truth = {'WBGene00000019', 'WBGene00000041', 'WBGene00000106',
                      'WBGene00001133', 'WBGene00003915', 'WBGene00268195'}

    def alt_calc_pval(self, n: int, log2fc: float, bg_array: np.ndarray, reps: int, obs_frac: float):
        assert reps == reps_truth
        assert n == len(gene_set_truth)
        assert log2fc == truth[4]
        assert isinstance(bg_array, np.ndarray)
        assert obs_frac == truth[2] / truth[1]
        assert np.all(bg_array == bg_array_truth)
        return 0.05

    monkeypatch.setattr(EnrichmentRunner, '_calc_randomization_pval', alt_calc_pval)
    runner = EnrichmentRunner.__new__(EnrichmentRunner)
    runner.gene_set = gene_set_truth
    runner.annotation_df = df
    assert runner._randomization_enrichment(truth[0], reps_truth) == truth


@pytest.mark.parametrize('params,truth',
                         [((38, 11, 6, 5), ['attribute', 11, 5, (6 / 38) * 11, np.log2(5 / ((6 / 38) * 11)), 0.05]),
                          ((40, 10, 15, 0), ['attribute', 10, 0, (15 / 40) * 10, -np.inf, 0.05])])
def test_enrichment_runner_fisher_enrichment(monkeypatch, params, truth):
    monkeypatch.setattr(EnrichmentRunner, '_get_hypergeometric_parameters', lambda self, attr: params)

    def alt_calc_pval(self, bg_size, de_size, go_size, go_de_size):
        assert (bg_size, de_size, go_size, go_de_size) == params
        return 0.05

    monkeypatch.setattr(EnrichmentRunner, '_calc_fisher_pval', alt_calc_pval)
    runner = EnrichmentRunner.__new__(EnrichmentRunner)
    assert runner._fisher_enrichment('attribute') == truth


def test_enrichment_runner_update_gene_set():
    runner = EnrichmentRunner.__new__(EnrichmentRunner)
    runner.single_list = False
    runner.background_set = set(pd.read_csv('tests/test_files/attr_ref_table_for_tests.csv', index_col=0).index)
    runner.gene_set = {'WBGene00000019', 'WBGene00000041', 'WBGene00000106', 'WBGene00001133', 'WBGene00003915',
                       'WBGene99991111'}
    updated_gene_set_truth = {'WBGene00000019', 'WBGene00000041', 'WBGene00000106', 'WBGene00001133', 'WBGene00003915'}
    runner.update_gene_set()
    assert runner.gene_set == updated_gene_set_truth


def test_enrichment_runner_update_gene_set_single_list(monkeypatch):
    monkeypatch.setattr(EnrichmentRunner, '_update_ranked_genes', lambda x: None)
    runner = EnrichmentRunner.__new__(EnrichmentRunner)
    runner.single_list = True
    runner.annotation_df = pd.read_csv('tests/test_files/attr_ref_table_for_tests.csv', index_col=0)
    runner.gene_set = {'WBGene00000019', 'WBGene00000041', 'WBGene00000106', 'WBGene00001133', 'WBGene00003915',
                       'WBGene99991111'}
    updated_gene_set_truth = {'WBGene00000019', 'WBGene00000041', 'WBGene00000106', 'WBGene00001133', 'WBGene00003915'}
    runner.update_gene_set()
    assert runner.gene_set == updated_gene_set_truth


def test_enrichment_runner_api():
    assert False


def test_enrichment_runner_format_results(monkeypatch):
    monkeypatch.setattr(EnrichmentRunner, '_correct_multiple_comparisons', lambda self: None)
    runner = EnrichmentRunner.__new__(EnrichmentRunner)
    results_list = [['name1', 50, 10, 5.5, 2.3, 0.05], ['name2', 17, 0, 3, 0, 1], ['name3', 1, np.nan, -2, -0.7, 0.04]]
    truth = pd.read_csv('tests/test_files/enrichment_runner_format_results_truth.csv', index_col=0)
    runner.en_score_col = 'colName'
    runner.single_list = False

    runner.format_results(results_list)
    assert truth.equals(runner.results)


def test_enrichment_runner_format_results_single_list(monkeypatch):
    monkeypatch.setattr(EnrichmentRunner, '_correct_multiple_comparisons', lambda self: None)
    runner = EnrichmentRunner.__new__(EnrichmentRunner)
    results_list = [['name1', 50, 2.3, 0.05], ['name2', 17, 0, 1], ['name3', 1, -0.7, np.nan]]
    truth = pd.read_csv('tests/test_files/enrichment_runner_single_list_format_results_truth.csv', index_col=0)
    runner.en_score_col = 'colName'
    runner.single_list = True

    runner.format_results(results_list)
    assert truth.equals(runner.results)


def test_enrichment_runner_xlmhg_enrichment():
    runner = EnrichmentRunner.__new__(EnrichmentRunner)
    assert False


@pytest.mark.parametrize('attribute,truth',
                         [('attribute1', np.array([0, 1], dtype='uint16')),
                          ('attribute4', np.array([0, 2], dtype='uint16'))])
def test_enrichment_runner_xlmhg_index_vector(attribute, truth):
    runner = EnrichmentRunner.__new__(EnrichmentRunner)
    runner.ranked_genes = np.array(['WBGene00000106', 'WBGene00000019', 'WBGene00000865', 'WBGene00001131'],
                                   dtype='str')
    runner.annotation_df = pd.read_csv('tests/test_files/attr_ref_table_for_tests.csv', index_col=0)
    assert np.all(runner._xlmhg_index_vector(attribute) == truth)


def test_enrichment_runner_fetch_annotations():
    runner = EnrichmentRunner.__new__(EnrichmentRunner)
    assert False


@pytest.mark.parametrize('attributes,truth',
                         [])
def test_enrichment_runner_fetch_attributes(attributes, truth):
    runner = EnrichmentRunner.__new__(EnrichmentRunner)
    runner.attributes = attributes
    assert False


@pytest.mark.parametrize('attibute_list,all_attrs,is_legal',
                         [(['a', 'b', 'c'], {'c', 'a', 'b', 'd', 'A', 'f'}, True),
                          (['a', 'c', 1], ['c', 'b', 'a'], True),
                          ([3, 2, 1], {'a', 'b', 'c'}, False),
                          (['a', -1, 2], ['a', 'b', 'c'], False),
                          (['A', 'b', 'c'], {'c', 'b', 'a'}, False),
                          (['a', 'b', True], {'a', 'b', True}, False),
                          (['a', 'b', 'c'], 'abc', False),
                          ('abc', {'a', 'b', 'c', 'abc'}, False)])
def test_enrichment_runner_validate_attributes(attibute_list, all_attrs, is_legal):
    runner = EnrichmentRunner.__new__(EnrichmentRunner)
    if is_legal:
        runner._validate_attributes(attibute_list, all_attrs)
    else:
        with pytest.raises(AssertionError):
            runner._validate_attributes(attibute_list, all_attrs)


def test_enrichment_runner_filter_annotations():
    truth = pd.read_csv('tests/test_files/enrichment_runner_filter_annotations_truth.csv', index_col=0)
    runner = EnrichmentRunner.__new__(EnrichmentRunner)
    runner.annotation_df = pd.read_csv('tests/test_files/attr_ref_table_for_tests.csv', index_col=0)
    runner.background_set = {'WBGene00000019', 'WBGene00000106', 'WBGene00000137', 'WBGene00000369', 'WBGene00000860',
                             'WBGene00048865', 'WBGene00268195'}
    runner.attributes = ['attribute1', 'attribute3', 'attribute4']
    runner.single_list = False
    runner.filter_annotations()
    assert truth.equals(runner.annotation_df)


def test_enrichment_runner_filter_annotations_single_list():
    truth = pd.read_csv('tests/test_files/enrichment_runner_filter_annotations_single_list_truth.csv', index_col=0)
    runner = EnrichmentRunner.__new__(EnrichmentRunner)
    runner.annotation_df = pd.read_csv('tests/test_files/attr_ref_table_for_tests.csv', index_col=0)
    runner.attributes = ['attribute1', 'attribute3', 'attribute4']
    runner.single_list = True
    runner.filter_annotations()
    assert truth.equals(runner.annotation_df)


def test_enrichment_runner_calculate_enrichment(monkeypatch):
    random_seed_status = [False]

    def set_seed(self):
        random_seed_status[0] = True

    monkeypatch.setattr(EnrichmentRunner, 'set_random_seed', set_seed)
    monkeypatch.setattr(EnrichmentRunner, '_calculate_enrichment_parallel', lambda self: 'parallel')
    monkeypatch.setattr(EnrichmentRunner, '_calculate_enrichment_serial', lambda self: 'serial')

    runner = EnrichmentRunner.__new__(EnrichmentRunner)
    runner.parallel = True
    assert runner.calculate_enrichment() == 'parallel'
    assert random_seed_status[0]

    random_seed_status[0] = False
    runner.parallel = False
    assert runner.calculate_enrichment() == 'serial'
    assert random_seed_status[0]


@pytest.mark.parametrize('seed,is_legal',
                         [(5, True),
                          (42, True),
                          (-1, False),
                          (0.1, False),
                          ('seed', False)])
def test_enrichment_runner_set_random_seed(seed, is_legal):
    runner = EnrichmentRunner.__new__(EnrichmentRunner)
    runner.random_seed = seed
    if is_legal:
        runner.set_random_seed()
        val = np.random.random()
        np.random.seed(seed)
        truth = np.random.random()
        assert val == truth
    else:
        with pytest.raises(AssertionError):
            runner.set_random_seed()


def _test_enrichment_runner_calculate_enrichment_get_constants():
    runner = EnrichmentRunner.__new__(EnrichmentRunner)
    runner.attributes = ['attribute1', 'attribute2', 'attribute4']
    runner.pvalue_kwargs = {'arg1': 'val1', 'arg2': 'val2'}

    def enrichment_func(attr, **kwargs):
        return [attr, kwargs]

    runner.enrichment_func = enrichment_func
    truth = [[attr, runner.pvalue_kwargs] for attr in runner.attributes]
    return runner, truth


def test_enrichment_runner_calculate_enrichment_parallel():
    runner, truth = _test_enrichment_runner_calculate_enrichment_get_constants()
    assert truth == runner._calculate_enrichment_parallel()


def test_enrichment_runner_calculate_enrichment_serial():
    runner, truth = _test_enrichment_runner_calculate_enrichment_get_constants()
    assert truth == runner._calculate_enrichment_serial()


def test_enrichment_runner_correct_multiple_comparisons():
    runner = EnrichmentRunner.__new__(EnrichmentRunner)

    runner.results = pd.DataFrame([[0, 0, 0, 0, 0, 0.005],
                                   [0, 0, 0, 0, 0, 0.017],
                                   [0, 0, 0, 0, 0, np.nan],
                                   [0, 0, 0, 0, -3, np.nan],
                                   [0, 0, 0, 0, 0, 0.92]],
                                  columns=['name', 'samples', 'obs', 'exp', 'colName', 'pval']).set_index('name')
    runner.alpha = 0.02
    truth = pd.DataFrame([[0, 0, 0, 0, 0, 0.005, 0.015, True],
                          [0, 0, 0, 0, 0, 0.017, 0.0255, False],
                          [0, 0, 0, 0, 0, np.nan, np.nan, np.nan],
                          [0, 0, 0, 0, -3, np.nan, np.nan, np.nan],
                          [0, 0, 0, 0, 0, 0.92, 0.92, False]],
                         columns=['name', 'samples', 'obs', 'exp', 'colName', 'pval', 'padj', 'significant']).set_index(
        'name')

    runner._correct_multiple_comparisons()

    assert np.all(np.isclose(truth['padj'].values, runner.results['padj'].values, atol=0, equal_nan=True))
    for val, val_truth in zip(runner.results['significant'], truth['significant']):
        assert val == val_truth or (np.isnan(val) and np.isnan(val_truth))
    assert truth.loc['name':'pval'].equals(runner.results.loc['name':'pval'])


def test_enrichment_runner_plot_results():
    runner = EnrichmentRunner.__new__(EnrichmentRunner)
    assert False


def test_enrichment_runner_enrichment_bar_plot():
    runner = EnrichmentRunner.__new__(EnrichmentRunner)
    assert False


def test_noncategorical_enrichment_runner_api():
    runner = NonCategoricalEnrichmentRunner.__new__(NonCategoricalEnrichmentRunner)
    assert False


def test_noncategorical_enrichment_runner_get_enrichment_func():
    runner = NonCategoricalEnrichmentRunner.__new__(NonCategoricalEnrichmentRunner)
    assert False


def test_noncategorical_enrichment_runner_sign_test_enrichment():
    runner = NonCategoricalEnrichmentRunner.__new__(NonCategoricalEnrichmentRunner)
    assert False


def test_noncategorical_enrichment_runner_one_sample_t_test_enrichment():
    runner = NonCategoricalEnrichmentRunner.__new__(NonCategoricalEnrichmentRunner)
    assert False


def test_noncategorical_enrichment_runner_format_results():
    runner = NonCategoricalEnrichmentRunner.__new__(NonCategoricalEnrichmentRunner)
    assert False


def test_noncategorical_enrichment_runner_plot_results():
    runner = NonCategoricalEnrichmentRunner.__new__(NonCategoricalEnrichmentRunner)
    assert False


def test_noncategorical_enrichment_runner_enrichment_histogram():
    runner = NonCategoricalEnrichmentRunner.__new__(NonCategoricalEnrichmentRunner)
    assert False


def test_go_enrichment_runner_api():
    runner = GOEnrichmentRunner.__new__(GOEnrichmentRunner)
    assert False


def test_go_enrichment_runner_run():
    runner = GOEnrichmentRunner.__new__(GOEnrichmentRunner)
    assert False


def test_go_enrichment_runner_get_enrichment_func():
    runner = GOEnrichmentRunner.__new__(GOEnrichmentRunner)
    assert False


def test_go_enrichment_runner_get_organism():
    runner = GOEnrichmentRunner.__new__(GOEnrichmentRunner)
    assert False


def test_go_enrichment_runner_fetch_annotations():
    runner = GOEnrichmentRunner.__new__(GOEnrichmentRunner)
    assert False


def test_go_enrichment_runner_get_annotation_iterator():
    runner = GOEnrichmentRunner.__new__(GOEnrichmentRunner)
    assert False


def test_go_enrichment_runner_propagate_annotations():
    runner = GOEnrichmentRunner.__new__(GOEnrichmentRunner)
    assert False


def test_go_enrichment_runner_translate_gene_ids():
    runner = GOEnrichmentRunner.__new__(GOEnrichmentRunner)
    assert False


def test_go_enrichment_runner_get_query_key():
    runner = GOEnrichmentRunner.__new__(GOEnrichmentRunner)
    assert False


def test_go_enrichment_runner_fetch_attributes():
    runner = GOEnrichmentRunner.__new__(GOEnrichmentRunner)
    assert False


def test_go_enrichment_runner_correct_multiple_comparisons():
    runner = GOEnrichmentRunner.__new__(GOEnrichmentRunner)

    runner.results = pd.DataFrame([[0, 0, 0, 0, 0, 0.005],
                                   [0, 0, 0, 0, 0, 0.017],
                                   [0, 0, 0, 0, 0, np.nan],
                                   [0, 0, 0, 0, -3, np.nan],
                                   [0, 0, 0, 0, 0, 0.92]],
                                  columns=['go_id', 'samples', 'obs', 'exp', 'colName', 'pval']).set_index('go_id')
    runner.alpha = 0.04
    truth = pd.DataFrame([[0, 0, 0, 0, 0, 0.005, 0.0275, True],
                          [0, 0, 0, 0, 0, 0.017, 0.04675, False],
                          [0, 0, 0, 0, 0, np.nan, np.nan, np.nan],
                          [0, 0, 0, 0, -3, np.nan, np.nan, np.nan],
                          [0, 0, 0, 0, 0, 0.92, 1.0, False]],
                         columns=['go_id', 'samples', 'obs', 'exp', 'colName', 'pval', 'padj',
                                  'significant']).set_index('go_id')

    runner._correct_multiple_comparisons()

    assert np.all(np.isclose(truth['padj'].values, runner.results['padj'].values, atol=0, equal_nan=True))
    for val, val_truth in zip(runner.results['significant'], truth['significant']):
        assert val == val_truth or (np.isnan(val) and np.isnan(val_truth))
    assert truth.loc['go_id':'pval'].equals(runner.results.loc['go_id':'pval'])


def test_go_enrichment_runner_plot_results():
    runner = GOEnrichmentRunner.__new__(GOEnrichmentRunner)
    assert False


def test_go_enrichment_runner_format_results(monkeypatch):
    monkeypatch.setattr(GOEnrichmentRunner, '_correct_multiple_comparisons', lambda self: None)
    runner = GOEnrichmentRunner.__new__(GOEnrichmentRunner)
    results_dict = {'name1': ['name1', 50, 10, 5.5, 2.3, 0.05], 'name2': ['name2', 17, 0, 3, 0, 1],
                    'name3': ['name3', 1, np.nan, -2, -0.7, 0.04]}
    truth = pd.read_csv('tests/test_files/go_enrichment_runner_format_results_truth.csv', index_col=0)
    runner.en_score_col = 'colName'
    runner.single_list = False
    runner.return_nonsignificant = False

    runner.format_results(results_dict)
    assert truth.equals(runner.results)


@pytest.mark.parametrize('return_nonsignificant,truth_file',
                         [(True, 'tests/test_files/go_enrichment_runner_format_results_with_nonsignificant_truth.csv'),
                          (False, 'tests/test_files/go_enrichment_runner_format_results_truth.csv')])
def test_go_enrichment_runner_format_results(monkeypatch, return_nonsignificant, truth_file):
    truth = pd.read_csv(truth_file, index_col=0)
    results_dict = {'name1': ['desc1', 50, 10, 5, 2.3, 0.05], 'name2': ['desc2', 17, 0, 3, 0, 1],
                    'name3': ['desc3', 1, np.nan, -2, -0.7, 0.04]}

    def add_sig(self):
        self.results['significant'] = [False, True, True]

    monkeypatch.setattr(GOEnrichmentRunner, '_correct_multiple_comparisons', add_sig)

    dag_tree = ontology.DAGTree.__new__(ontology.DAGTree)
    dag_tree.go_terms = {'name1': ontology.GOTerm(), 'name2': ontology.GOTerm(), 'name3': ontology.GOTerm()}
    dag_tree['name1'].level = 2
    dag_tree['name2'].level = 1
    dag_tree['name3'].level = 5

    runner = GOEnrichmentRunner.__new__(GOEnrichmentRunner)
    runner.dag_tree = dag_tree
    runner.en_score_col = 'colName'
    runner.single_list = False
    runner.return_nonsignificant = return_nonsignificant

    runner.format_results(results_dict)

    assert truth.equals(runner.results)


def test_go_enrichment_runner_calculate_enrichment_serial():
    runner = GOEnrichmentRunner.__new__(GOEnrichmentRunner)
    assert False


def test_go_enrichment_runner_calculate_enrichment_parallel():
    runner = GOEnrichmentRunner.__new__(GOEnrichmentRunner)
    assert False


def test_go_enrichment_runner_go_classic_pvalues_serial():
    runner = GOEnrichmentRunner.__new__(GOEnrichmentRunner)
    assert False


def test_go_enrichment_runner_go_classic_pvalues_parallel():
    runner = GOEnrichmentRunner.__new__(GOEnrichmentRunner)
    assert False


def test_go_enrichment_runner_go_elim_pvalues_serial():
    runner = GOEnrichmentRunner.__new__(GOEnrichmentRunner)
    assert False


def test_go_enrichment_runner_go_elim_pvalues_parallel():
    runner = GOEnrichmentRunner.__new__(GOEnrichmentRunner)
    assert False


def test_go_enrichment_runner_go_weight_pvalues_serial():
    runner = GOEnrichmentRunner.__new__(GOEnrichmentRunner)
    assert False


def test_go_enrichment_runner_go_weight_pvalues_parallel():
    runner = GOEnrichmentRunner.__new__(GOEnrichmentRunner)
    assert False


def test_go_enrichment_runner_go_allm_pvalues_serial():
    runner = GOEnrichmentRunner.__new__(GOEnrichmentRunner)
    assert False


def test_go_enrichment_runner_go_allm_pvalues_parallel():
    runner = GOEnrichmentRunner.__new__(GOEnrichmentRunner)
    assert False


def test_go_enrichment_runner_parallel_over_grouping():
    runner = GOEnrichmentRunner.__new__(GOEnrichmentRunner)
    assert False


def test_go_enrichment_runner_calculate_allm():
    runner = GOEnrichmentRunner.__new__(GOEnrichmentRunner)
    assert False


def test_go_enrichment_runner_go_level_iterator():
    runner = GOEnrichmentRunner.__new__(GOEnrichmentRunner)
    assert False


def test_go_enrichment_runner_compute_term_sig():
    runner = GOEnrichmentRunner.__new__(GOEnrichmentRunner)
    assert False


def test_go_enrichment_runner_randomization_enrichment():
    runner = GOEnrichmentRunner.__new__(GOEnrichmentRunner)
    assert False


def test_go_enrichment_runner_xlmhg_enrichment():
    runner = GOEnrichmentRunner.__new__(GOEnrichmentRunner)
    assert False


@pytest.mark.parametrize('attribute,truth',
                         [('attribute1', np.array([0, 1], dtype='uint16')),
                          ('attribute4', np.array([0, 2], dtype='uint16'))])
def test_go_enrichment_runner_xlmhg_index_vector(attribute, truth):
    runner = GOEnrichmentRunner.__new__(GOEnrichmentRunner)
    runner.ranked_genes = np.array(['WBGene00000106', 'WBGene00000019', 'WBGene00000865', 'WBGene00001131'],
                                   dtype='str')
    runner.annotation_df = pd.read_csv('tests/test_files/attr_ref_table_for_tests.csv', index_col=0).notna()
    assert np.all(runner._xlmhg_index_vector(attribute) == truth)


def test_go_enrichment_runner_hypergeometric_enrichment():
    runner = GOEnrichmentRunner.__new__(GOEnrichmentRunner)
    assert False


def test_go_enrichment_runner_fisher_enrichment():
    runner = GOEnrichmentRunner.__new__(GOEnrichmentRunner)
    assert False


def test_go_enrichment_runner_get_hypergeometric_parameters():
    runner = GOEnrichmentRunner.__new__(GOEnrichmentRunner)
    assert False
