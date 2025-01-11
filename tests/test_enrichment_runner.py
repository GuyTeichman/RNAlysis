import copy
from collections import namedtuple

import matplotlib
import pytest

from rnalysis import filtering
from rnalysis.utils import enrichment_runner
from rnalysis.utils.enrichment_runner import *
from rnalysis.utils.io import *
from tests import __attr_ref__

matplotlib.use('Agg')


def _df_to_dict(df, null_mode: bool = True):
    if null_mode:
        return {attr: set(df.filter(pl.col(attr).is_not_null()).select(pl.first()).to_series()) for attr in
                df.columns[1:]}
    return {attr: set(df.filter(pl.col(attr) is True).select(pl.first()).to_series()) for attr in df.columns[1:]}


def test_enrichment_runner_from_results():
    alpha = 0.05
    plot_horizontal = False
    set_name = 'set_name'
    results = io.load_table('tests/test_files/enrichment_hypergeometric_res.csv')
    runner = EnrichmentRunner.from_results_df(results, alpha, plot_horizontal, set_name)

    assert runner.results is results
    assert runner.alpha == 0.05
    assert runner.plot_horizontal == plot_horizontal
    assert runner.set_name == set_name


@pytest.mark.parametrize("pval, alpha, expected", [
    (0.6, None, ('ns', 'normal')),
    (0.001, 0.00099, ('ns', 'normal')),
    (0.04, None, (u'\u2217', 'bold')),
    (0.0099, None, (u'\u2217' * 2, 'bold')),
    (0, None, (u'\u2217' * 4, 'bold'))
])
def test_get_pval_asterisk(pval, alpha, expected):
    if alpha is not None:
        result = EnrichmentPlotter._get_pval_asterisk(pval, alpha)
    else:
        result = EnrichmentPlotter._get_pval_asterisk(pval)
    assert result == expected


def test_calc_randomization_pval():
    np.random.seed(42)
    hypergeom_pval = 0.2426153598589023
    avg_pval = 0
    for i in range(5):
        avg_pval += PermutationTest._calc_permutation_pval(1, 10000, 0.11, 10000, 500, 1000)
    avg_pval /= 5
    assert np.isclose(avg_pval, hypergeom_pval, atol=0.02)


@pytest.mark.parametrize('M,n,N,X,truth', [
    (13588, 59, 611, 19, 4.989682834519698 * 10 ** -12),
    (20000, 430, 700, 6, 0.006249179131697138),
    (20000, 43, 300, 3, 0.0265186938062861)
])
def test_calc_hypergeometric_pvalues(M, n, N, X, truth):
    pval = HypergeometricTest._calc_hg_pval(M, n, N, X)
    assert np.isclose(truth, pval, atol=0, rtol=0.00001)


@pytest.mark.parametrize('attributes,truth', [
    ([0, 2, 3], ['attribute1', 'attribute3', 'attribute4']),
    (1, ['attribute2'])
])
def test_enrichment_get_attrs_int_index_attributes(attributes, truth):
    genes = {'WBGene00000041', 'WBGene00002074', 'WBGene00000105', 'WBGene00000106', 'WBGene00199484',
             'WBGene00001436', 'WBGene00000137', 'WBGene00001996', 'WBGene00014208', 'WBGene00001133'}
    e = EnrichmentRunner(genes, attributes, 0.05, __attr_ref__, True, False, '', False, True, 'test_set', False,
                         HypergeometricTest(), None)
    e.fetch_annotations()
    e.get_attribute_names()
    e.get_background_set()
    e.update_gene_set()
    e.filter_annotations()
    assert e.attributes == truth


def test_enrichment_get_attrs_all_attributes():
    genes = {'WBGene00000041', 'WBGene00002074', 'WBGene00000105', 'WBGene00000106', 'WBGene00199484',
             'WBGene00001436', 'WBGene00000137', 'WBGene00001996', 'WBGene00014208', 'WBGene00001133'}
    e = EnrichmentRunner(genes, 'all', 0.05, __attr_ref__, True, False, '', False, True, 'test_set', False,
                         HypergeometricTest(), None, )
    e.fetch_annotations()
    e.get_attribute_names()
    e.get_background_set()
    e.update_gene_set()
    e.filter_annotations()
    attrs = e.attributes
    attrs_truth = ['attribute1', 'attribute2', 'attribute3', 'attribute4']
    plt.close('all')
    assert attrs == attrs_truth


def test_enrichment_get_attrs_bad_path():
    e = EnrichmentRunner({'_'}, 'attribute1', 0.05, 'fakepath', True, False, '', False, True, 'test_set', False,
                         HypergeometricTest())
    with pytest.raises((FileNotFoundError, AssertionError)):
        e.fetch_annotations()
        e.get_attribute_names()
        e.get_background_set()
        e.update_gene_set()
        e.filter_annotations()


def _enrichment_get_ref_tests_setup(truth, bg_genes, **kwargs):
    genes = {'WBGene00000041', 'WBGene00002074', 'WBGene00000019', 'WBGene00000105', 'WBGene00000106', 'WBGene00199484',
             'WBGene00001436', 'WBGene00000137', 'WBGene00001996', 'WBGene00014208'}
    background = None if isinstance(bg_genes, str) else bg_genes
    e = EnrichmentRunner(genes, 'all', 0.05, __attr_ref__, True, False, '', False, 'test_set', False,
                         HypergeometricTest(), background, **kwargs)
    e.fetch_annotations()
    e.get_attribute_names()
    e.get_background_set()
    e.update_gene_set()
    e.filter_annotations()
    res = e.annotations
    truth = _df_to_dict(truth)
    assert res == truth


def test_enrichment_get_ref():
    truth = io.load_table('tests/test_files/attr_ref_table_for_tests_specified_bg.csv')
    bg_genes = {'WBGene00003902', 'WBGene00000106', 'WBGene00001436', 'WBGene00000864', 'WBGene00011910',
                'WBGene00000859', 'WBGene00268189', 'WBGene00000865', 'WBGene00003864', 'WBGene00048863',
                'WBGene00000369', 'WBGene00000863', 'WBGene00002074', 'WBGene00000041', 'WBGene00199486',
                'WBGene00000105', 'WBGene00001131'}
    _enrichment_get_ref_tests_setup(truth, bg_genes)


def test_enrichment_get_ref_include_unannotated():
    truth = io.load_table('tests/test_files/attr_ref_table_for_tests_specified_bg_unannotated.csv')
    bg_genes = {'WBGene00003902', 'WBGene00000106', 'WBGene00001436', 'WBGene00000864', 'WBGene00011910',
                'WBGene00000859', 'WBGene00268189', 'WBGene00000865', 'WBGene00003864', 'WBGene00048863',
                'WBGene00000369', 'WBGene00000863', 'WBGene00002074', 'WBGene00000041', 'WBGene00199486',
                'WBGene00000105', 'WBGene00001131', 'gene1', 'gene2', 'gene3'}
    _enrichment_get_ref_tests_setup(truth, bg_genes, exclude_unannotated_genes=False)


def test_enrichment_get_ref_from_featureset_object():
    truth = io.load_table('tests/test_files/attr_ref_table_for_tests_specified_bg.csv')
    bg_genes = {'WBGene00003902', 'WBGene00000106', 'WBGene00001436', 'WBGene00000864', 'WBGene00011910',
                'WBGene00000859', 'WBGene00268189', 'WBGene00000865', 'WBGene00003864', 'WBGene00048863',
                'WBGene00000369', 'WBGene00000863', 'WBGene00002074', 'WBGene00000041', 'WBGene00199486',
                'WBGene00000105', 'WBGene00001131'}
    _enrichment_get_ref_tests_setup(truth, bg_genes)


def test_enrichment_get_ref_from_filter_object():
    truth = io.load_table('tests/test_files/attr_ref_table_for_tests_specified_bg.csv')
    bg_genes = filtering.CountFilter(r'tests/test_files/test_bg_genes_from_filter_object.csv').index_set
    _enrichment_get_ref_tests_setup(truth, bg_genes)


def test_results_to_csv():
    try:
        en = EnrichmentRunner({''}, 'all', 0.05, __attr_ref__, True, True, 'tests/test_files/tmp_enrichment_csv.csv',
                              False, True, 'test_set', False, HypergeometricTest(), None)
        df = io.load_table('tests/test_files/enrichment_hypergeometric_res.csv')
        en.results = df
        en.results_to_csv()
        df_loaded = io.load_table('tests/test_files/tmp_enrichment_csv.csv')
        assert df.equals(df_loaded)
    except Exception as e:
        raise e
    finally:
        try:
            os.remove('tests/test_files/tmp_enrichment_csv.csv')
        except:
            pass


def _compare_go_result_dfs(res, truth):
    res = res.sort(pl.first())
    truth = truth.sort(pl.first())
    assert res.select(pl.col('GO ID', 'n', 'obs')).equals(truth.select(pl.col('GO ID', 'n', 'obs')))
    assert np.allclose(res.select(pl.col('exp', 'log2fc', 'pval')), res.select(pl.col('exp', 'log2fc', 'pval')), atol=0)


def test_classic_pvals(monkeypatch):
    annotations = _df_to_dict(io.load_table('tests/test_files/goa_table.csv'), null_mode=False)
    gene_set = {'gene1', 'gene2', 'gene5', 'gene12', 'gene13', 'gene17', 'gene19', 'gene25', 'gene27', 'gene28'}
    background = {f'gene{i + 1}' for i in range(30)}
    truth = io.load_table('tests/test_files/go_pvalues_classic_truth.csv')
    dummy_go_node = namedtuple('DummyGONode', field_names='name')

    class DummyDAGTree:
        def __init__(self):
            self.dummy_node = dummy_go_node('content of name field')

        def __getitem__(self, item):
            return self.dummy_node

    monkeypatch.setattr(ontology, 'fetch_go_basic', lambda: DummyDAGTree())

    e = GOEnrichmentRunner(gene_set, 'elegans', 'WBGene', 0.05, 'classic', 'any', 'any', None, 'any', None, 'any', None,
                           False, False, '', False, False, '', False, HypergeometricTest(), background)
    e.annotations = annotations
    e.mutable_annotations = annotations,
    e.attributes = list(annotations.keys())

    res = pl.DataFrame([i for i in e._go_classic_on_batch(tuple(annotations.keys()), 0).values()],
                       schema=['GO ID', 'n', 'obs', 'exp', 'log2fc', 'pval'])
    _compare_go_result_dfs(res, truth)


def test_elim_pvals(monkeypatch):
    annotations = _df_to_dict(io.load_table('tests/test_files/goa_table.csv'), null_mode=False)
    threshold = 0.2  # make sure there are both significant and non-significant examples with our small bg size (30)
    gene_set = {'gene1', 'gene2', 'gene5', 'gene12', 'gene13', 'gene17', 'gene19', 'gene25', 'gene27', 'gene28'}
    background = {f'gene{i + 1}' for i in range(30)}
    truth = io.load_table('tests/test_files/go_pvalues_elim_truth.csv')
    with open('tests/test_files/obo_for_go_tests.obo', 'r') as f:
        dag_tree = ontology.DAGTree(f, ['is_a'])
    monkeypatch.setattr(ontology, 'fetch_go_basic', lambda: dag_tree)

    e = GOEnrichmentRunner(gene_set, 'elegans', 'WBGene', threshold, 'classic', 'any', 'any', None, 'any', None, 'any',
                           None, False, False, '', False, False, '', False, HypergeometricTest(), background)
    e.annotations = annotations
    e.mutable_annotations = deepcopy(annotations),
    e.attributes = list(annotations.keys())
    e.attributes_set = set(e.attributes)

    res = pl.DataFrame([i for i in e._go_elim_on_aspect('all').values()],
                       schema=['GO ID', 'n', 'obs', 'exp', 'log2fc', 'pval'])

    _compare_go_result_dfs(res, truth)


def test_weight_pvals(monkeypatch):
    goa_df = io.load_table('tests/test_files/goa_table.csv')
    gene_set = {'gene1', 'gene2', 'gene5', 'gene12', 'gene13', 'gene17', 'gene19', 'gene25', 'gene27', 'gene28'}
    truth = io.load_table('tests/test_files/go_pvalues_weight_truth.csv')
    with open('tests/test_files/obo_for_go_tests.obo', 'r') as f:
        dag_tree = ontology.DAGTree(f, ['is_a'])

    monkeypatch.setattr(ontology, 'fetch_go_basic', lambda: dag_tree)

    e = GOEnrichmentRunner(gene_set, 'elegans', 'WBGene', 0.05, 'classic', 'any', 'any', None, 'any', None, 'any', None,
                           False, False, '', False, False, '', False, HypergeometricTest())
    e.annotations = _df_to_dict(goa_df, null_mode=False)
    e.mutable_annotations = {attr: {v: 1 for v in val} for attr, val in e.annotations.items()},
    e.attributes = list(goa_df.columns)
    e.attributes_set = set(e.attributes)
    e.background_set = set.union(*e.annotations.values())

    res = pl.DataFrame([i for i in e._go_weight_on_aspect('all').values()],
                       schema=['GO ID', 'n', 'obs', 'exp', 'log2fc', 'pval'])
    _compare_go_result_dfs(res, truth)


def test_allm_pvals(monkeypatch):
    annotations = _df_to_dict(io.load_table('tests/test_files/goa_table.csv'), null_mode=False)
    threshold = 0.2  # make sure there are both significant and non-significant examples with our small bg size (30)
    gene_set = {'gene1', 'gene2', 'gene5', 'gene12', 'gene13', 'gene17', 'gene19', 'gene25', 'gene27', 'gene28'}
    background = {f'gene{i + 1}' for i in range(30)}
    truth = io.load_table('tests/test_files/go_pvalues_allm_truth.csv')
    with open('tests/test_files/obo_for_go_tests.obo', 'r') as f:
        dag_tree = ontology.DAGTree(f, ['is_a'])
    monkeypatch.setattr(ontology, 'fetch_go_basic', lambda: dag_tree)

    e = GOEnrichmentRunner(gene_set, 'elegans', 'WBGene', threshold, 'classic', 'any', 'any', None, 'any', None, 'any',
                           None, False, False, '', False, False, '', False, HypergeometricTest(), background)
    e.annotations = annotations
    e.mutable_annotations = deepcopy(annotations),
    e.attributes = list(annotations.keys())
    e.attributes_set = set(e.attributes)

    res = pl.DataFrame([val for val in e._go_allm_pvalues_serial().values()],
                       schema=['GO ID', 'n', 'obs', 'exp', 'log2fc', 'pval'])
    _compare_go_result_dfs(res, truth)


def test_enrichment_runner_update_ranked_genes():
    runner = EnrichmentRunner.__new__(EnrichmentRunner)
    truth = np.array([['111', '000', '222', '777', '333', '888']], dtype='str')
    runner.ranked_genes = np.array(['111', '000', '222', '777', '333', '000', '111', '888', '000'], dtype='str')
    runner.gene_set = {'111', '000', '222', '777', '333', '888'}
    runner._update_ranked_genes()
    assert np.all(truth == runner.ranked_genes)


@pytest.mark.parametrize("attr,truth",
                         [('attribute1', (38, 6, 11, 3, 3, 1.736842105263158, 0.7884958948062882)),
                          ('attribute3', (38, 6, 14, 4, 4, 2.2105263157894735, 0.8556100906648252)),
                          ('attribute4', (38, 6, 13, 4, 4, 2.0526315789473686, 0.962525294581337))])
def test_get_hypergeometric_parameters(attr, truth):
    df = io.load_table('tests/test_files/attr_ref_table_for_tests.csv')
    background = parsing.data_to_set(df.select(pl.first()))
    annotation = _df_to_dict(df)
    gene_set = {'WBGene00000019', 'WBGene00000041', 'WBGene00000106', 'WBGene00001133', 'WBGene00003915',
                'WBGene00268195'}

    pval_func = HypergeometricTest()
    res = pval_func.get_hypergeometric_params(annotation[attr], gene_set, background)
    assert res == truth


@pytest.mark.parametrize('params,truth',
                         [((38, 11, 6, 5, 5, (6 / 38) * 11, np.log2(5 / ((6 / 38) * 11))),
                           ['attribute', 11, 5, (6 / 38) * 11, np.log2(5 / ((6 / 38) * 11)), 0.05]),
                          ((40, 10, 15, 0, 0, (15 / 40) * 10, -np.inf),
                           ['attribute', 10, 0, (15 / 40) * 10, -np.inf, 0.05])])
def test_enrichment_runner_hypergeometric_enrichment(monkeypatch, params, truth):
    stats_test = HypergeometricTest()
    monkeypatch.setattr(HypergeometricTest, 'get_hypergeometric_params', lambda self, *args: params)

    def alt_calc_pval(self, bg_size, de_size, go_size, go_de_size):
        assert (bg_size, de_size, go_size, go_de_size) == params[:4]
        return 0.05

    monkeypatch.setattr(HypergeometricTest, '_calc_hg_pval', alt_calc_pval)
    assert stats_test.run('attribute', set(), set(), set()) == truth


@pytest.mark.parametrize('truth', [(['attribute1', 6, 3, (11 / 38) * 6, np.log2(3 / ((11 / 38) * 6)), 0.05]),
                                   (['attribute4', 6, 4, (13 / 38) * 6, np.log2(4 / ((13 / 38) * 6)), 0.05])])
def test_enrichment_runner_randomization_enrichment(monkeypatch, truth):
    reps_truth = 100
    annotations = _df_to_dict(io.load_table('tests/test_files/attr_ref_table_for_tests.csv'))
    bg_array_truth = io.load_table('tests/test_files/annotation_df_bg_array_truth.csv')[truth[0]]
    gene_set_truth = {'WBGene00000019', 'WBGene00000041', 'WBGene00000106',
                      'WBGene00001133', 'WBGene00003915', 'WBGene00268195'}

    def alt_calc_pval(self, log2fc: float, reps: int, obs_frac: float, bg_size: int, en_size: int, attr_size: int):
        assert reps == reps_truth
        assert en_size == len(gene_set_truth)
        assert log2fc == truth[4]
        assert sum(bg_array_truth) == attr_size
        assert len(bg_array_truth) == bg_size
        assert obs_frac == truth[2] / truth[1]
        return 0.05

    monkeypatch.setattr(PermutationTest, '_calc_permutation_pval', alt_calc_pval)
    gene_set = gene_set_truth
    stats_test = PermutationTest(reps_truth)
    background_set = {'WBGene00268190', 'WBGene00048865', 'WBGene00000019', 'WBGene00048864', 'WBGene00003864',
                      'WBGene00000865', 'WBGene00000864', 'WBGene00048863', 'WBGene00003865', 'WBGene00000041',
                      'WBGene00001132', 'WBGene00001133', 'WBGene00004920', 'WBGene00001996', 'WBGene00001134',
                      'WBGene00011910', 'WBGene00199486', 'WBGene00000369', 'WBGene00000861', 'WBGene00014208',
                      'WBGene00003902', 'WBGene00001436', 'WBGene00000105', 'WBGene00268195', 'WBGene00003915',
                      'WBGene00255735', 'WBGene00001131', 'WBGene00000106', 'WBGene00268189', 'WBGene00268191',
                      'WBGene00000860', 'WBGene00002074', 'WBGene00000137', 'WBGene00000863', 'WBGene00199484',
                      'WBGene00000859', 'WBGene00199485', 'WBGene00255734'}

    attr = truth[0]
    res = stats_test.run(attr, annotations[attr], gene_set, background_set)
    print(res)
    assert res == truth


@pytest.mark.parametrize('params,truth',
                         [((38, 11, 6, 5, 5, (6 / 38) * 11, np.log2(5 / ((6 / 38) * 11))),
                           ['attribute', 11, 5, (6 / 38) * 11, np.log2(5 / ((6 / 38) * 11)), 0.05]),
                          ((40, 10, 15, 0, 0, (15 / 40) * 10, -np.inf),
                           ['attribute', 10, 0, (15 / 40) * 10, -np.inf, 0.05])])
def test_enrichment_runner_fisher_enrichment(monkeypatch, params, truth):
    stats_test = FishersExactTest()
    monkeypatch.setattr(FishersExactTest, 'get_hypergeometric_params', lambda self, *args: params)

    def alt_calc_pval(self, bg_size, de_size, go_size, go_de_size):
        assert (bg_size, de_size, go_size, go_de_size) == params[:4]
        return 0.05

    monkeypatch.setattr(FishersExactTest, '_calc_fisher_pval', alt_calc_pval)
    assert stats_test.run('attribute', set(), set(), set()) == truth


@pytest.mark.parametrize('exclude_unannotated', [True, False])
def test_enrichment_runner_update_gene_set(exclude_unannotated):
    runner = EnrichmentRunner.__new__(EnrichmentRunner)
    runner.single_set = False
    runner.exclude_unannotated_genes = exclude_unannotated
    runner.background_set = parsing.data_to_set(
        io.load_table('tests/test_files/attr_ref_table_for_tests.csv').select(pl.first()))
    runner.gene_set = {'WBGene00000019', 'WBGene00000041', 'WBGene00000106', 'WBGene00001133', 'WBGene00003915',
                       'WBGene99991111'}
    updated_gene_set_truth = {'WBGene00000019', 'WBGene00000041', 'WBGene00000106', 'WBGene00001133', 'WBGene00003915'}
    runner.update_gene_set()
    assert runner.gene_set == updated_gene_set_truth


def test_enrichment_runner_update_gene_set_single_list(monkeypatch):
    monkeypatch.setattr(EnrichmentRunner, '_update_ranked_genes', lambda x: None)
    runner = EnrichmentRunner.__new__(EnrichmentRunner)
    runner.single_set = True
    runner.annotations = _df_to_dict(io.load_table('tests/test_files/attr_ref_table_for_tests.csv'))
    runner.gene_set = {'WBGene00000019', 'WBGene00000041', 'WBGene00000106', 'WBGene00001133', 'WBGene00003915',
                       'WBGene99991111'}
    updated_gene_set_truth = {'WBGene00000019', 'WBGene00000041', 'WBGene00000106', 'WBGene00001133', 'WBGene00003915'}
    runner.update_gene_set()
    assert runner.gene_set == updated_gene_set_truth


@pytest.mark.parametrize('exclude_unannotated', [True, False])
@pytest.mark.parametrize('save_csv,', [True, False])
@pytest.mark.parametrize('return_nonsignificant,', [True, False])
@pytest.mark.parametrize('fname', ['fname', None])
@pytest.mark.parametrize(
    'single_list,genes,pval_func,background_set',
    [(True, np.array(['WBGene1', 'WBGene2'], dtype=str), XlmhgTest(), None),
     (False, {'WBGene00000001', 'WBGene00000002'}, PermutationTest(100, 42),
      {'WBGene00000001', 'WBGene00000002', 'EBGene00000003'})])
def test_enrichment_runner_api(exclude_unannotated, return_nonsignificant, save_csv, fname, single_list, genes,
                               pval_func, background_set, monkeypatch):
    monkeypatch.setattr('builtins.input', lambda x: 'fname')

    runner = EnrichmentRunner(genes, ['attr1', 'attr2'], 0.05, 'path/to/attr/ref', return_nonsignificant, save_csv,
                              fname, False, 'set_name', False, pval_func, background_set,
                              exclude_unannotated, single_list)
    if save_csv:
        assert runner.fname == 'fname'
    else:
        with pytest.raises(AttributeError):
            _ = runner.fname


def test_enrichment_runner_format_results(monkeypatch):
    monkeypatch.setattr(EnrichmentRunner, '_correct_multiple_comparisons', lambda self: None)
    runner = EnrichmentRunner.__new__(EnrichmentRunner)
    results_list = [['name1', 50, 10, 5.5, 2.3, 0.05], ['name2', 17, 0, 3, 0, 1], ['name3', 1, np.nan, -2, -0.7, 0.04]]
    truth = io.load_table('tests/test_files/enrichment_runner_format_results_truth.csv')
    runner.en_score_col = 'colName'
    runner.single_set = False
    runner.return_nonsignificant = True

    runner.format_results(results_list)
    assert truth.equals(runner.results.fill_nan(None))


def test_enrichment_runner_format_results_single_list(monkeypatch):
    monkeypatch.setattr(EnrichmentRunner, '_correct_multiple_comparisons', lambda self: None)
    runner = EnrichmentRunner.__new__(EnrichmentRunner)
    results_list = [['name1', 50, 5, 10, 2.3, 0.05], ['name2', 17, 17, 0.56, 0, 1], ['name3', 1, 2, 3.14, -0.7, np.nan]]
    truth = io.load_table('tests/test_files/enrichment_runner_single_list_format_results_truth.csv')
    runner.en_score_col = 'colName'
    runner.single_set = True
    runner.return_nonsignificant = True

    runner.format_results(results_list)
    assert truth.equals(runner.results.fill_nan(None))


@pytest.mark.parametrize('attr,index_vector, pval_kwargs,p_values,e_scores,params_truth',
                         [('attr1', [3, 5, 2], {}, (0.05, 0.9), (2, 1.5),
                           {'L': 1, 'X': 10, 'N': 10, 'table': np.empty((3 + 1, 10 - 3 + 1), dtype=np.longdouble)}),
                          ('attr 2', [1, 2, 3, 4], {'l_fraction': 0.2, 'x': 7}, (0, 0), (0, 0),
                           {'L': 2, 'X': 7, 'N': 10, 'table': np.empty((4 + 1, 10 - 4 + 1), dtype=np.longdouble)}),
                          ('attr,3', [2, 9], {'other_arg': 13}, (0.8, 0.1), (1, 1.2),
                           {'L': 1, 'X': 10, 'N': 10, 'table': np.empty((2 + 1, 10 - 2 + 1), dtype=np.longdouble)})])
def test_xlmhg_enrichment(monkeypatch, attr, index_vector, pval_kwargs, p_values, e_scores, params_truth):
    n_calls_xlmhg_test = [0]
    params_truth['indices'] = index_vector
    obs = 2
    exp = 5 * (0.5 / float(10))

    class ResultObject:
        def __init__(self, pval, escore):
            self.pval = pval
            self.escore = escore
            self.k = 2
            self.K = 5
            self.N = 10
            self.cutoff = 0.5

    monkeypatch.setattr(XlmhgTest, '_generate_xlmhg_index_vectors',
                        lambda self, rankedgenes, ann_set: (index_vector, index_vector))

    def _xlmhg_test_validate_parameters(**kwargs):
        for key in ['X', 'L', 'N', 'indices']:
            assert kwargs[key] == params_truth[key]
        assert kwargs['table'].shape == params_truth['table'].shape
        assert kwargs['table'].dtype == params_truth['table'].dtype

        pval = p_values[n_calls_xlmhg_test[0]]
        escore = e_scores[n_calls_xlmhg_test[0]]
        n_calls_xlmhg_test[0] += 1
        return ResultObject(pval, escore)

    monkeypatch.setattr(xlmhglite, 'get_xlmhg_test_result', _xlmhg_test_validate_parameters)

    ranked_genes = np.array(
        ['gene0', 'gene1', 'gene2', 'gene3', 'gene4', 'gene5', 'gene6', 'gene7', 'gene8', 'gene9'], dtype='str')

    log2fc = np.log2(e_scores[0] if p_values[0] <= p_values[1] else (1 / e_scores[1]))
    output_truth = [attr, len(ranked_genes), obs, exp, log2fc, min(p_values)]

    stats = XlmhgTest(pval_kwargs.get('x', 10), pval_kwargs.get('l_fraction', 0.1))
    result = stats.run(attr, set(), ranked_genes)
    assert result == output_truth


@pytest.mark.parametrize('attribute,truth, truth_rev',
                         [('attribute1', np.array([0, 1], dtype='uint16'), np.array([2, 3], dtype='uint16')),
                          ('attribute4', np.array([0, 2], dtype='uint16'), np.array([1, 3], dtype='uint16'))])
def test_generate_xlmhg_index_vectors(attribute, truth, truth_rev):
    ranked_genes = np.array(['WBGene00000106', 'WBGene00000019', 'WBGene00000865', 'WBGene00001131'], dtype='str')
    annotations = _df_to_dict(io.load_table('tests/test_files/attr_ref_table_for_tests.csv'))
    vec, rev_vec = XlmhgTest()._generate_xlmhg_index_vectors(ranked_genes, annotations[attribute])
    assert np.all(vec == truth)
    assert np.all(rev_vec == truth_rev)


def test_enrichment_runner_fetch_annotations(monkeypatch):
    monkeypatch.setattr(validation, 'validate_attr_table', lambda x: None)
    truth = _df_to_dict(io.load_table('tests/test_files/attr_ref_table_for_tests.csv'))
    runner = EnrichmentRunner.__new__(EnrichmentRunner)
    runner.attr_ref_path = 'tests/test_files/attr_ref_table_for_tests.csv'
    runner.fetch_annotations()
    assert truth == runner.annotations


@pytest.mark.parametrize('attributes,truth,annotation_df_cols',
                         [(['all', ['col1', 'col2', 'col3'], ['col1', 'col2', 'col3']]),
                          (['col2', 'col3'], ['col2', 'col3'], ['col1', 'col2', 'col3']),
                          ([0, 2], ['col1', 'col3'], ['col1', 'col2', 'col3']),
                          ({0, 2}, ['col1', 'col3'], ['col1', 'col2', 'col3']),
                          ({'col2', 'col3'}, ['col2', 'col3'], ['col1', 'col2', 'col3'])])
def test_enrichment_runner_fetch_attributes(attributes, truth, annotation_df_cols, monkeypatch):
    monkeypatch.setattr(EnrichmentRunner, '_validate_attributes', lambda self, attrs, all_attrs: None)
    runner = EnrichmentRunner.__new__(EnrichmentRunner)
    runner.attributes = attributes
    runner.annotations = {attr: set() for attr in annotation_df_cols}
    runner.get_attribute_names()
    assert sorted(runner.attributes) == sorted(truth)


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
    truth = _df_to_dict(io.load_table('tests/test_files/enrichment_runner_filter_annotations_truth.csv'))
    runner = EnrichmentRunner.__new__(EnrichmentRunner)
    runner.annotations = _df_to_dict(io.load_table('tests/test_files/attr_ref_table_for_tests.csv'))
    runner.background_set = {'WBGene00000019', 'WBGene00000106', 'WBGene00000137', 'WBGene00000369', 'WBGene00000860',
                             'WBGene00048865', 'WBGene00268195'}
    runner.attributes = ['attribute1', 'attribute3', 'attribute4']
    runner.single_set = False
    runner.filter_annotations()
    assert truth == runner.annotations


def test_enrichment_runner_filter_annotations_single_list():
    truth = _df_to_dict(io.load_table('tests/test_files/enrichment_runner_filter_annotations_single_list_truth.csv'))
    runner = EnrichmentRunner.__new__(EnrichmentRunner)
    ref_table = io.load_table('tests/test_files/attr_ref_table_for_tests.csv')
    runner.annotations = _df_to_dict(ref_table)
    runner.attributes = ['attribute1', 'attribute3', 'attribute4']
    runner.single_set = True
    runner.filter_annotations()
    assert truth == runner.annotations


def test_enrichment_runner_calculate_enrichment(monkeypatch):
    monkeypatch.setattr(EnrichmentRunner, '_calculate_enrichment_parallel', lambda self: 'parallel')
    monkeypatch.setattr(EnrichmentRunner, '_calculate_enrichment_serial', lambda self: 'serial')

    runner = EnrichmentRunner.__new__(EnrichmentRunner)
    runner.attributes = list(range(10))

    for backend in ['loky', 'threading', 'multiprocessing']:
        runner.parallel_backend = backend
        assert runner.calculate_enrichment() == 'parallel'

    runner.parallel_backend = 'sequential'
    assert runner.calculate_enrichment() == 'serial'


def test_enrichment_runner_correct_multiple_comparisons():
    runner = EnrichmentRunner.__new__(EnrichmentRunner)

    runner.results = pl.DataFrame([[0, 0, 0, 0, 0, 0.005],
                                   [0, 0, 0, 0, 0, 0.017],
                                   [0, 0, 0, 0, 0, np.nan],
                                   [0, 0, 0, 0, -3, np.nan],
                                   [0, 0, 0, 0, 0, 0.92]],
                                  schema=['name', 'samples', 'obs', 'exp', 'colName', 'pval'])
    runner.alpha = 0.02
    truth = pl.DataFrame([[0, 0, 0, 0, 0, 0.005, 0.015, True],
                          [0, 0, 0, 0, 0, 0.017, 0.0255, False],
                          [0, 0, 0, 0, 0, 0.92, 0.92, False],
                          [0, 0, 0, 0, 0, np.nan, None, None],
                          [0, 0, 0, 0, -3, np.nan, None, None],
                          ],
                         schema=['name', 'samples', 'obs', 'exp', 'colName', 'pval', 'padj', 'significant'])
    runner._correct_multiple_comparisons()
    truth = truth.sort(pl.col('padj', 'colName'))
    res = runner.results.sort(pl.col('padj', 'colName'))

    assert np.all(np.isclose(truth['padj'], res['padj'], atol=0, equal_nan=True))
    for val, val_truth in zip(res['significant'], truth['significant']):
        assert val == val_truth
    compared_cols = ('name', 'samples', 'obs', 'exp', 'colName', 'pval')
    assert truth.sort(pl.col('padj', 'colName')).select(compared_cols).equals(
        res.sort(pl.col('padj', 'colName')).select(compared_cols))


@pytest.mark.parametrize('single_list', [True, False])
def test_enrichment_runner_plot_results(monkeypatch, single_list):
    def validate_params(self):
        assert isinstance(self.title, str)
        assert runner.set_name in self.title
        assert isinstance(self.ylabel, str)
        return [plt.Figure()]

    monkeypatch.setattr(BarPlotter, 'run', validate_params)
    runner = EnrichmentRunner.__new__(EnrichmentRunner)
    runner.single_set = single_list
    runner.set_name = 'name of the set'
    runner.results = pl.DataFrame([1])
    runner.en_score_col = 'column_0'
    runner.alpha = 0.05

    res = runner.plot_results().run()
    assert validation.isinstanceiter(res, plt.Figure)


@pytest.mark.parametrize('n_bars', ['all', 2])
@pytest.mark.parametrize('plot_horizontal', [True, False])
@pytest.mark.parametrize('show_expected', [True, False])
@pytest.mark.parametrize('plot_style', ['bar', 'lollipop'])
def test_enrichment_runner_enrichment_bar_plot(plot_horizontal, n_bars, plot_style, show_expected):
    runner = EnrichmentRunner.__new__(EnrichmentRunner)
    runner.results = io.load_table('tests/test_files/enrichment_hypergeometric_res.csv')
    runner.en_score_col = 'log2_fold_enrichment'
    runner.alpha = 0.05
    runner.plot_horizontal = plot_horizontal
    runner.show_expected = show_expected
    runner.single_set = False
    runner.set_name = 'set name'
    runner.plot_style = plot_style
    runner.plot_results().run()


def test_noncategorical_enrichment_runner_api():
    NonCategoricalEnrichmentRunner({'gene1', 'gene2', 'gene4'}, ['attr1', 'attr2'], 0.05,
                                   {'gene1', 'gene2', 'gene3', 'gene4'}, 'path/to/attr/ref',
                                   False, 'fname', False, 'overlap', 5,
                                   'set_name', False, True)


@pytest.mark.parametrize("attr,gene_set,truth",
                         [('attr5', {f'WBGene0000000{i + 1}' for i in range(4)},
                           ['attr5', 4, 16.85912542, 6.986579003, 0.05]),
                          ('attr4', {f'WBGene0000000{i + 1}' for i in range(4)},
                           ['attr4', 4, 15.5, 14, 0.05]),
                          ('attr1', {'WBGene00000001', 'WBGene00000029', 'WBGene00000030'},
                           ['attr1', 3, np.nan, 3, np.nan])])
def test_sign_test_enrichment(monkeypatch, attr, gene_set, truth):
    df = io.load_table('tests/test_files/attr_ref_table_for_non_categorical.csv')
    background_set = {f'WBGene000000{i + 1:02d}' for i in range(30)}

    def validate_params(values, popmean):
        assert np.isclose(popmean, truth[3])
        assert np.all(sorted(values) == sorted(df.filter(pl.first().is_in(gene_set))[attr].to_list()))
        return None, 0.05

    monkeypatch.setattr(enrichment_runner, 'sign_test', validate_params)
    runner = NonCategoricalEnrichmentRunner.__new__(NonCategoricalEnrichmentRunner)
    runner.annotations = df
    runner.gene_set = gene_set
    runner.background_set = background_set
    runner.attr_ref_path = 'tests/test_files/attr_ref_table_for_non_categorical.csv'
    runner.fetch_annotations()
    res = SignTest().run(attr, runner.annotations[attr], gene_set, background_set)
    for res_val, truth_val in zip(res, truth):
        if isinstance(truth_val, float):
            assert np.isclose(truth_val, res_val, equal_nan=True)
        else:
            assert truth_val == res_val


@pytest.mark.parametrize("attr,gene_set,truth",
                         [('attr5', {f'WBGene0000000{i + 1}' for i in range(4)},
                           ['attr5', 4, 14.93335321, 7.709139128, 0.05]),
                          ('attr4', {f'WBGene0000000{i + 1}' for i in range(4)},
                           ['attr4', 4, 14.25, 13.60714286, 0.05]),
                          ('attr1', {'WBGene00000001', 'WBGene00000029', 'WBGene00000030'},
                           ['attr1', 3, np.nan, 4.75, np.nan])])
def test_one_sample_t_test_enrichment(monkeypatch, attr, gene_set, truth):
    df = io.load_table('tests/test_files/attr_ref_table_for_non_categorical.csv')
    background_set = {f'WBGene000000{i + 1:02d}' for i in range(30)}

    def validate_params(values, popmean):
        assert np.isclose(popmean, truth[3])
        assert np.all(sorted(values) == sorted(df.filter(pl.first().is_in(gene_set))[attr].to_list()))
        return None, 0.05

    monkeypatch.setattr(enrichment_runner, 'ttest_1samp', validate_params)
    runner = NonCategoricalEnrichmentRunner.__new__(NonCategoricalEnrichmentRunner)
    runner.annotations = df
    runner.gene_set = gene_set
    runner.background_set = background_set
    runner.attr_ref_path = 'tests/test_files/attr_ref_table_for_non_categorical.csv'
    runner.fetch_annotations()
    res = TTest().run(attr, runner.annotations[attr], gene_set, background_set)
    for res_val, truth_val in zip(res, truth):
        if isinstance(truth_val, float):
            assert np.isclose(truth_val, res_val, equal_nan=True)
        else:
            assert truth_val == res_val


def test_noncategorical_enrichment_runner_format_results(monkeypatch):
    monkeypatch.setattr(NonCategoricalEnrichmentRunner, '_correct_multiple_comparisons', lambda self: None)
    runner = NonCategoricalEnrichmentRunner.__new__(NonCategoricalEnrichmentRunner)
    results_list = [['name1', 50, 10, 5.5, 0.05], ['name2', 17, 0, 3, 1], ['name3', 1, None, -2, None]]
    truth = io.load_table('tests/test_files/non_categorical_enrichment_runner_format_results_truth.csv')
    runner.return_nonsignificant = True
    runner.format_results(results_list)
    print(truth)
    print(runner.results)
    assert truth.equals(runner.results)


def test_noncategorical_enrichment_runner_plot_results(monkeypatch):
    attrs_without_nan_truth = ['attr1', 'attr2', 'attr4']
    attrs_without_nan = []

    def proccess_attr(self, attr):
        attrs_without_nan.append(attr)
        return plt.Figure()

    monkeypatch.setattr(HistogramPlotter, '_plot_histogram', proccess_attr)
    runner = NonCategoricalEnrichmentRunner.__new__(NonCategoricalEnrichmentRunner)
    runner.attributes = ['attr1', 'attr2', 'attr3', 'attr4']
    runner.results = pl.DataFrame({'attributes': runner.attributes, 'padj': [0.05, 0.3, np.nan, 1]})
    runner.annotations = {'attr1': [], 'attr2': [], 'attr3': [], 'attr4': []}
    runner.alpha = 0.05
    runner.gene_set = set()
    runner.background_set = set()
    runner.plot_log_scale = True
    runner.n_bins = 10
    runner.plot_style = 'interleaved'
    runner.set_name = 'set name'
    runner.parametric_test = True

    res = runner.plot_results().run()
    assert attrs_without_nan == attrs_without_nan_truth
    assert len(res) == len(attrs_without_nan)
    for plot in res:
        assert isinstance(plot, plt.Figure)


@pytest.mark.parametrize('plot_style', ['interleaved', 'overlap'])
@pytest.mark.parametrize('plot_log_scale', [True, False])
@pytest.mark.parametrize('parametric_test', [True, False])
@pytest.mark.parametrize('n_bins', [2, 8])
def test_noncategorical_enrichment_runner_enrichment_histogram(plot_style, plot_log_scale, parametric_test, n_bins):
    runner = NonCategoricalEnrichmentRunner.__new__(NonCategoricalEnrichmentRunner)
    runner.plot_style = plot_style
    runner.plot_log_scale = plot_log_scale
    runner.parametric_test = parametric_test
    runner.n_bins = n_bins
    runner.alpha = 0.15
    runner.set_name = 'set_name'
    runner.background_set = {f'WBGene000000{i + 1:02d}' for i in range(30)}
    runner.results = io.load_table('tests/test_files/enrich_non_categorical_nan_parametric_truth.csv')
    runner.attr_ref_path = 'tests/test_files/attr_ref_table_for_non_categorical.csv'
    runner.gene_set = {f'WBGene0000000{i + 1}' for i in range(5)}
    runner.fetch_annotations()
    runner.attributes = ['attr5']

    res = runner.plot_results().run()
    assert validation.isinstanceiter(res, plt.Figure)


@pytest.mark.parametrize('exclude_unannotated', [True, False])
@pytest.mark.parametrize('single_list,genes,pval_func,background_set',
                         [(True, np.array(['WBGene1', 'WBGene2'], dtype=str), XlmhgTest(), None),
                          (False, {'WBGene00000001', 'WBGene00000002'}, PermutationTest(100, 42),
                           {'WBGene00000001', 'WBGene00000002', 'EBGene00000003'},)])
def test_go_enrichment_runner_api(monkeypatch, single_list, genes, pval_func, background_set, exclude_unannotated):
    monkeypatch.setattr(ontology, 'fetch_go_basic', lambda: 'dag_tree')
    runner = GOEnrichmentRunner(genes, 'organism', 'gene_id_type', 0.05, 'elim', 'any', 'any', None, 'any', None, 'any',
                                None, False, False, 'fname', False, False, 'set_name', False, pval_func,
                                background_set, exclude_unannotated, single_list)

    assert runner.dag_tree == 'dag_tree'


def test_go_enrichment_runner_run(monkeypatch):
    organism_truth = 'my_organism'

    def run(self):
        self.results = 'results'

    monkeypatch.setattr(EnrichmentRunner, 'run', run)
    runner = GOEnrichmentRunner.__new__(GOEnrichmentRunner)
    runner.taxon_id, runner.organism = 'taxon_id', organism_truth
    runner.run()
    assert runner.results == 'results'
    assert runner.organism == organism_truth
    assert runner.taxon_id == 'taxon_id'


def test_go_enrichment_runner_fetch_annotations(monkeypatch):
    monkeypatch.setattr(GOEnrichmentRunner, '_get_query_key', lambda self: 'the_query_key')
    monkeypatch.setattr(GOEnrichmentRunner, '_generate_annotation_dict', lambda self: 'goa_df')
    runner = GOEnrichmentRunner.__new__(GOEnrichmentRunner)
    runner.fetch_annotations()
    assert runner.annotations == 'goa_df'
    assert runner.GOA_DF_QUERIES['the_query_key'] == 'goa_df'

    runner.GOA_DF_QUERIES['the_query_key'] = 'another_goa_df'
    runner.fetch_annotations()
    assert runner.annotations == 'another_goa_df'
    assert runner.GOA_DF_QUERIES['the_query_key'] == 'another_goa_df'


def test_go_enrichment_runner_get_annotation_iterator(monkeypatch):
    def alt_init(self, taxon_id, aspects, evidence_types, excluded_evidence_types, databases, excluded_databases,
                 qualifiers, excluded_qualifiers):
        assert taxon_id == 'taxon_id'
        assert aspects == 'aspects'
        assert evidence_types == 'evidence_types'
        assert excluded_evidence_types == 'exc_evidence_types'
        assert databases == 'databases'
        assert excluded_databases == 'exc_databases'
        assert qualifiers == 'qualifiers'
        assert excluded_qualifiers == 'exc_qualifiers'

    monkeypatch.setattr(io.GOlrAnnotationIterator, '__init__', alt_init)
    runner = GOEnrichmentRunner.__new__(GOEnrichmentRunner)
    runner.taxon_id = 'taxon_id'
    runner.aspects = 'aspects'
    runner.evidence_types = 'evidence_types'
    runner.excluded_evidence_types = 'exc_evidence_types'
    runner.databases = 'databases'
    runner.excluded_databases = 'exc_databases'
    runner.qualifiers = 'qualifiers'
    runner.excluded_qualifiers = 'exc_qualifiers'
    res = runner._get_annotation_iterator()
    assert isinstance(res, io.GOlrAnnotationIterator)


@pytest.mark.parametrize('propagate', ['other', 'no'])
@pytest.mark.parametrize("gene_id,go_id,truth",
                         [('gene1', 'GO1',
                           {'GO2': {'gene3', 'gene1'}, 'GO3': {'gene1'}, 'GO5': {'gene1'}, 'GO4': {'gene1'},
                            'GO1': {'gene2', 'gene1'}}),
                          ('gene2', 'GO1',
                           {'GO1': {'gene2', 'gene1'}, 'GO2': {'gene3', 'gene1'}, 'GO3': {'gene2'}, 'GO5': {'gene2'},
                            'GO4': {'gene2'}}
                           ),
                          ('gene1', 'GO2', {'GO1': {'gene2', 'gene1'}, 'GO2': {'gene3', 'gene1'}}),
                          ('gene3', 'GO2', {'GO1': {'gene2', 'gene1'}, 'GO2': {'gene3', 'gene1'}}),
                          ])
def test_go_enrichment_runner_propagate_annotation(monkeypatch, propagate, gene_id, go_id, truth):
    annotation_dict = {'GO1': {'gene1', 'gene2'}, 'GO2': {'gene1', 'gene3'}}
    if propagate == 'no':
        truth = copy.deepcopy(annotation_dict)

    def upper_induced_graph_iter(self, go_id):
        if go_id == 'GO1':
            for upper_go_id in {'GO3', 'GO4', 'GO5'}:
                yield upper_go_id

    monkeypatch.setattr(ontology.DAGTree, 'upper_induced_graph_iter', upper_induced_graph_iter)
    dag = ontology.DAGTree.__new__(ontology.DAGTree)
    runner = GOEnrichmentRunner.__new__(GOEnrichmentRunner)
    runner.dag_tree = dag
    runner.propagate_annotations = propagate
    runner._propagate_annotation(gene_id, go_id, annotation_dict)

    assert annotation_dict == truth


@pytest.mark.parametrize("mapping_dict,truth", [
    ({}, {'GO1': set(), 'GO2': set()}),
    ({'gene1': 'gene1_translated', 'gene3': 'gene3_translated'},
     {'GO1': {'gene1_translated'}, 'GO2': {'gene1_translated', 'gene3_translated'}}),
    ({'gene1': 'gene1_translated', 'gene2': 'gene2_translated', 'gene3': 'gene3_translated'},
     {'GO1': {'gene1_translated', 'gene2_translated'}, 'GO2': {'gene1_translated', 'gene3_translated'}})])
def test_go_enrichment_runner_translate_gene_ids(monkeypatch, mapping_dict, truth):
    monkeypatch.setattr(io.GeneIDTranslator, '__init__', lambda *args: None)
    monkeypatch.setattr(io.GeneIDTranslator, 'run', lambda *args: mapping_dict)
    source_to_gene_id_dict = {'source1': {'gene1', 'gene3'}, 'source2': {'gene2'}}
    sparse_annotation_dict = {'GO1': {'gene1', 'gene2'}, 'GO2': {'gene1', 'gene3'}}

    runner = GOEnrichmentRunner.__new__(GOEnrichmentRunner)
    runner.gene_id_type = 'gene_id_type'

    res = runner._translate_gene_ids(sparse_annotation_dict, source_to_gene_id_dict)
    assert res == truth


@pytest.mark.parametrize('propagate_annotations', ['no', 'elim'])
def test_go_enrichment_runner_get_query_key(propagate_annotations):
    runner = GOEnrichmentRunner.__new__(GOEnrichmentRunner)
    runner.propagate_annotations = propagate_annotations
    runner.taxon_id = 'taxon_id'
    runner.gene_id_type = 'id_type'
    runner.aspects = {'aspect1', 'aspect2'}
    runner.evidence_types = {'ev2', 'ev1', 'ev5'}
    runner.excluded_evidence_types = {'ev3'}
    runner.databases = {'db1'}
    runner.excluded_databases = set()
    runner.qualifiers = {'qual1'}
    runner.excluded_qualifiers = set()

    propagate = True if propagate_annotations != 'no' else False

    truth = (
        'taxon_id', 'id_type', ('aspect1', 'aspect2'), ('ev1', 'ev2', 'ev5'), ('ev3',), ('db1',), tuple(), ('qual1',),
        tuple(), propagate)
    key = runner._get_query_key()
    assert key == truth
    try:
        _ = hash(key)
    except TypeError:
        assert False


def test_go_enrichment_runner_fetch_attributes():
    runner = GOEnrichmentRunner.__new__(GOEnrichmentRunner)
    runner.annotations = _df_to_dict(io.load_table('tests/test_files/attr_ref_table_for_tests.csv'))
    truth_attributes = ['attribute1', 'attribute2', 'attribute3', 'attribute4']
    truth_attributes_set = {'attribute1', 'attribute2', 'attribute3', 'attribute4'}
    runner.get_attribute_names()
    assert runner.attributes == truth_attributes
    assert runner.attributes_set == truth_attributes_set


def test_go_enrichment_runner_correct_multiple_comparisons():
    runner = GOEnrichmentRunner.__new__(GOEnrichmentRunner)

    runner.results = pl.DataFrame([[0, 0, 0, 0, 0, 0.005],
                                   [0, 0, 0, 0, 0, 0.017],
                                   [0, 0, 0, 0, 0, np.nan],
                                   [0, 0, 0, 0, -3, np.nan],
                                   [0, 0, 0, 0, 0, 0.92]],
                                  schema=['GO ID', 'samples', 'obs', 'exp', 'colName', 'pval'])
    runner.alpha = 0.04
    truth = pl.DataFrame([[0, 0, 0, 0, 0, 0.005, 0.0275, True],
                          [0, 0, 0, 0, 0, 0.017, 0.04675, False],
                          [0, 0, 0, 0, 0, np.nan, None, None],
                          [0, 0, 0, 0, -3, np.nan, None, None],
                          [0, 0, 0, 0, 0, 0.92, 1.0, False]],
                         schema=['GO ID', 'samples', 'obs', 'exp', 'colName', 'pval', 'padj', 'significant'])

    runner._correct_multiple_comparisons()
    truth = truth.sort(pl.col('padj', 'colName'))
    res = runner.results.sort(pl.col('padj', 'colName'))

    assert np.all(np.isclose(truth['padj'], res['padj'], atol=0, equal_nan=True))
    for val, val_truth in zip(res['significant'], truth['significant']):
        assert val == val_truth
    compared_cols = ('GO ID', 'samples', 'obs', 'exp', 'colName', 'pval')
    assert truth.sort(pl.col('padj', 'colName')).select(compared_cols).equals(
        res.sort(pl.col('padj', 'colName')).select(compared_cols))


@pytest.mark.parametrize('single_list', [True, False])
@pytest.mark.parametrize('results,n_bars_truth', [([1, 2, 3], 3), (list(range(15)), 10)])
@pytest.mark.parametrize('plot_ontology_graph', [False, True])
def test_go_enrichment_runner_plot_results(monkeypatch, single_list, results, n_bars_truth, plot_ontology_graph):
    dag_plotted = []

    def validate_params_bar(self):
        assert isinstance(self.title, str)
        assert runner.set_name in self.title
        assert isinstance(self.ylabel, str)
        assert isinstance(self.n_bars, int)
        assert self.n_bars == n_bars_truth
        return [plt.Figure()]

    def validate_params_dag(self):
        assert plot_ontology_graph
        dag_plotted.append(True)
        return []

    monkeypatch.setattr(BarPlotter, 'run', validate_params_bar)
    monkeypatch.setattr(GODagPlotter, 'run', validate_params_dag)
    runner = GOEnrichmentRunner.__new__(GOEnrichmentRunner)
    runner.single_set = single_list
    runner.set_name = 'name of the set'
    runner.results = pl.DataFrame(results)
    runner.en_score_col = 'column_0'
    runner.alpha = 0.05
    runner.dag_tree = ontology.DAGTree.__new__(ontology.DAGTree)
    runner.aspects = 'any'
    runner.ontology_graph_format = 'none'
    runner.plot_ontology_graph = plot_ontology_graph
    res = runner.plot_results().run()
    assert validation.isinstanceiter(res, plt.Figure)
    if plot_ontology_graph:
        assert dag_plotted == [True]


@pytest.mark.parametrize('return_nonsignificant,truth_file',
                         [(True, 'tests/test_files/go_enrichment_runner_format_results_with_nonsignificant_truth.csv'),
                          (False, 'tests/test_files/go_enrichment_runner_format_results_truth.csv')])
def test_go_enrichment_runner_format_results(monkeypatch, return_nonsignificant, truth_file):
    truth = io.load_table(truth_file)
    results_dict = {'name1': ['name1', 50, 10, 5, 2.3, None], 'name2': ['name2', 17, 0, 3, 0, 1],
                    'name3': ['name3', 1, 1.3, -2, -0.7, 0.04]}

    def add_sig(self):
        self.results = self.results.with_columns(significant=pl.Series([False, True, True]))

    monkeypatch.setattr(GOEnrichmentRunner, '_correct_multiple_comparisons', add_sig)

    dag_tree = ontology.DAGTree.__new__(ontology.DAGTree)
    dag_tree.go_terms = {'name1': ontology.GOTerm(), 'name2': ontology.GOTerm(), 'name3': ontology.GOTerm()}
    dag_tree['name1'].set_level(2)
    dag_tree['name2'].set_level(1)
    dag_tree['name3'].set_level(5)
    dag_tree['name1'].set_name('desc1')
    dag_tree['name2'].set_name('desc2')
    dag_tree['name3'].set_name('desc3')
    dag_tree.alt_ids = {}

    runner = GOEnrichmentRunner.__new__(GOEnrichmentRunner)
    runner.dag_tree = dag_tree
    runner.en_score_col = 'colName'
    runner.single_set = False
    runner.return_nonsignificant = return_nonsignificant

    runner.format_results(results_dict)
    assert truth.equals(runner.results)


@pytest.mark.parametrize("propagate_annotations,truth",
                         [('no', 'classic'),
                          ('classic', 'classic'),
                          ('elim', 'elim'),
                          ('weight', 'weight'),
                          ('all.m', 'all.m')])
def test_go_enrichment_runner_calculate_enrichment_serial(monkeypatch, propagate_annotations, truth):
    class DAGTreePlaceHolder:
        def __init__(self):
            self.namespaces = ['namespace1', 'namespace2']

    def go_level_iter(self, namespace):
        mapper = {'namespace1': ['attribute1', 'attribute4'], 'namespace2': ['attribute2', 'attribute3']}
        if namespace == 'all':
            return [f'attribute{i + 1}' for i in range(4)]
        return mapper[namespace]

    monkeypatch.setattr(GOEnrichmentRunner, '_go_level_iterator', go_level_iter)
    monkeypatch.setattr(GOEnrichmentRunner, '_go_classic_pvalues_serial', lambda self, desc: 'classic')
    monkeypatch.setattr(GOEnrichmentRunner, '_go_elim_pvalues_serial', lambda self, desc: 'elim')
    monkeypatch.setattr(GOEnrichmentRunner, '_go_weight_pvalues_serial', lambda self, desc: 'weight')
    monkeypatch.setattr(GOEnrichmentRunner, '_go_allm_pvalues_serial', lambda self: 'all.m')
    runner = GOEnrichmentRunner.__new__(GOEnrichmentRunner)
    runner.propagate_annotations = propagate_annotations
    runner.attributes = []
    runner.annotations = _df_to_dict(io.load_table('tests/test_files/attr_ref_table_for_tests.csv'))
    runner.dag_tree = DAGTreePlaceHolder()
    res = runner._calculate_enrichment_serial()

    if propagate_annotations != 'all.m':
        assert len(runner.mutable_annotations) == 1

    if propagate_annotations in {'classic', 'no'}:
        assert runner.annotations is runner.mutable_annotations[0]
    elif propagate_annotations == 'elim':
        assert runner.annotations == runner.mutable_annotations[0]
        assert runner.annotations is not runner.mutable_annotations[0]
    elif propagate_annotations == 'weight':
        assert len(runner.mutable_annotations) == 1
        assert runner.mutable_annotations == ({
                                                  'attribute1': {'WBGene00000105': 1.0, 'WBGene00004920': 1.0,
                                                                 'WBGene00000137': 1.0, 'WBGene00001996': 1.0,
                                                                 'WBGene00014208': 1.0, 'WBGene00000106': 1.0,
                                                                 'WBGene00000019': 1.0, 'WBGene00002074': 1.0,
                                                                 'WBGene00000041': 1.0, 'WBGene00011910': 1.0,
                                                                 'WBGene00001436': 1.0},
                                                  'attribute2': {'WBGene00003902': 1.0, 'WBGene00014208': 1.0,
                                                                 'WBGene00003864': 1.0, 'WBGene00000369': 1.0,
                                                                 'WBGene00003865': 1.0, 'WBGene00011910': 1.0,
                                                                 'WBGene00003915': 1.0, 'WBGene00004920': 1.0},
                                                  'attribute3': {'WBGene00001996': 1.0, 'WBGene00000369': 1.0,
                                                                 'WBGene00000106': 1.0, 'WBGene00000019': 1.0,
                                                                 'WBGene00002074': 1.0, 'WBGene00268189': 1.0,
                                                                 'WBGene00000860': 1.0, 'WBGene00048864': 1.0,
                                                                 'WBGene00000859': 1.0, 'WBGene00001132': 1.0,
                                                                 'WBGene00199484': 1.0, 'WBGene00003915': 1.0,
                                                                 'WBGene00199486': 1.0, 'WBGene00001133': 1.0},
                                                  'attribute4': {'WBGene00000865': 1.0, 'WBGene00003902': 1.0,
                                                                 'WBGene00000137': 1.0, 'WBGene00001134': 1.0,
                                                                 'WBGene00000106': 1.0, 'WBGene00048865': 1.0,
                                                                 'WBGene00003865': 1.0, 'WBGene00268189': 1.0,
                                                                 'WBGene00268195': 1.0, 'WBGene00199485': 1.0,
                                                                 'WBGene00001132': 1.0, 'WBGene00003915': 1.0,
                                                                 'WBGene00001133': 1.0}},)

    assert res == truth


@pytest.mark.parametrize("propagate_annotations,truth",
                         [('no', 'classic'),
                          ('classic', 'classic'),
                          ('elim', 'elim'),
                          ('weight', 'weight'),
                          ('all.m', 'all.m')])
def test_go_enrichment_runner_calculate_enrichment_parallel(monkeypatch, propagate_annotations, truth):
    monkeypatch.setattr(GOEnrichmentRunner, '_go_classic_pvalues_parallel', lambda self, desc: 'classic')
    monkeypatch.setattr(GOEnrichmentRunner, '_go_elim_pvalues_parallel', lambda self, desc: 'elim')
    monkeypatch.setattr(GOEnrichmentRunner, '_go_weight_pvalues_parallel', lambda self, desc: 'weight')
    monkeypatch.setattr(GOEnrichmentRunner, '_go_allm_pvalues_parallel', lambda self: 'all.m')

    def go_level_iter(self, namespace):
        mapper = {'namespace1': ['attribute1', 'attribute4'], 'namespace2': ['attribute2', 'attribute3']}
        if namespace == 'all':
            return [f'attribute{i + 1}' for i in range(4)]
        return mapper[namespace]

    monkeypatch.setattr(GOEnrichmentRunner, '_go_level_iterator', go_level_iter)

    class DAGTreePlaceHolder:
        def __init__(self):
            self.namespaces = ['namespace1', 'namespace2']

    runner = GOEnrichmentRunner.__new__(GOEnrichmentRunner)
    runner.propagate_annotations = propagate_annotations
    runner.attributes = []
    runner.annotations = _df_to_dict(io.load_table('tests/test_files/attr_ref_table_for_tests.csv'))
    runner.dag_tree = DAGTreePlaceHolder()

    res = runner._calculate_enrichment_parallel()

    if propagate_annotations in {'classic', 'no'}:
        assert runner.annotations is runner.mutable_annotations[0]
        assert len(runner.mutable_annotations) == 1
    elif propagate_annotations == 'elim':
        assert len(runner.mutable_annotations) == len(runner.dag_tree.namespaces)
        for namespace, mod_ann in zip(runner.dag_tree.namespaces, runner.mutable_annotations):
            assert {go_id: runner.annotations[go_id] for go_id in go_level_iter(None, namespace)} == mod_ann
    elif propagate_annotations == 'weight':
        assert len(runner.mutable_annotations) == len(runner.dag_tree.namespaces)
        assert runner.mutable_annotations == ({'attribute1': {'WBGene00000105': 1.0, 'WBGene00001436': 1.0,
                                                              'WBGene00004920': 1.0, 'WBGene00000137': 1.0,
                                                              'WBGene00014208': 1.0, 'WBGene00000041': 1.0,
                                                              'WBGene00002074': 1.0, 'WBGene00000019': 1.0,
                                                              'WBGene00001996': 1.0, 'WBGene00000106': 1.0,
                                                              'WBGene00011910': 1.0},
                                               'attribute4': {'WBGene00003865': 1.0, 'WBGene00001134': 1.0,
                                                              'WBGene00268189': 1.0, 'WBGene00048865': 1.0,
                                                              'WBGene00000137': 1.0, 'WBGene00001133': 1.0,
                                                              'WBGene00003902': 1.0, 'WBGene00001132': 1.0,
                                                              'WBGene00000865': 1.0, 'WBGene00199485': 1.0,
                                                              'WBGene00000106': 1.0, 'WBGene00003915': 1.0,
                                                              'WBGene00268195': 1.0}}, {
                                                  'attribute2': {'WBGene00003865': 1.0, 'WBGene00003864': 1.0,
                                                                 'WBGene00004920': 1.0, 'WBGene00014208': 1.0,
                                                                 'WBGene00003902': 1.0, 'WBGene00011910': 1.0,
                                                                 'WBGene00000369': 1.0, 'WBGene00003915': 1.0},
                                                  'attribute3': {'WBGene00268189': 1.0, 'WBGene00199484': 1.0,
                                                                 'WBGene00002074': 1.0, 'WBGene00001133': 1.0,
                                                                 'WBGene00048864': 1.0, 'WBGene00000859': 1.0,
                                                                 'WBGene00000369': 1.0, 'WBGene00000019': 1.0,
                                                                 'WBGene00001132': 1.0, 'WBGene00000860': 1.0,
                                                                 'WBGene00001996': 1.0, 'WBGene00000106': 1.0,
                                                                 'WBGene00003915': 1.0, 'WBGene00199486': 1.0}})

    assert res == truth


def test_go_enrichment_runner_process_annotations_no_annotations(monkeypatch):
    class AnnotationIterator:
        def __init__(self, n_annotations):
            self.n_annotations = n_annotations

        def __iter__(self):
            return [None].__iter__()

    def _get_annotation_iter_zero(self):
        return AnnotationIterator(0)

    monkeypatch.setattr(GOEnrichmentRunner, '_get_annotation_iterator', _get_annotation_iter_zero)
    with pytest.raises(AssertionError) as e:
        runner = GOEnrichmentRunner.__new__(GOEnrichmentRunner)
        runner.propagate_annotations = "no"
        runner.organism = "organism"
        runner.taxon_id = "taxon id"
        runner._process_annotations()
    print(e)


def test_go_enrichment_runner_process_annotations(monkeypatch):
    propagated_annotations = set()
    propagated_annotations_truth = {('gene_id1', 'go_id1'), ('gene_id1', 'go_id2'), ('gene_id2', 'go_id1'),
                                    ('gene_id2', 'go_id3'), ('gene_id3', 'go_id4')}
    annotation_dict_truth = {'go_id1': {'gene_id1', 'gene_id2'}, 'go_id2': {'gene_id1'}, 'go_id3': {'gene_id2'},
                             'go_id4': {'gene_id3'}}
    source_dict_truth = {'source1': {'gene_id1', 'gene_id3'}, 'source2': {'gene_id2'}}

    def get_annotation_iter(self):
        iterator = io.GOlrAnnotationIterator.__new__(io.GOlrAnnotationIterator)
        iterator.n_annotations = 5
        return iterator

    def annotation_iter(self):
        annotations = [{'bioentity_internal_id': 'gene_id1', 'annotation_class': 'go_id1', 'source': 'source1'},
                       {'bioentity_internal_id': 'gene_id1', 'annotation_class': 'go_id2', 'source': 'source1'},
                       {'bioentity_internal_id': 'gene_id2', 'annotation_class': 'go_id1', 'source': 'source2'},
                       {'bioentity_internal_id': 'gene_id2', 'annotation_class': 'go_id3', 'source': 'source2'},
                       {'bioentity_internal_id': 'gene_id3', 'annotation_class': 'go_id4', 'source': 'source1'}]
        for annotation in annotations:
            yield annotation

    def propagate_annotation(self, gene_id, go_id, sparse_annotation_dict):
        propagated_annotations.add((gene_id, go_id))

    monkeypatch.setattr(io.GOlrAnnotationIterator, '_annotation_generator_func', annotation_iter)
    monkeypatch.setattr(GOEnrichmentRunner, '_propagate_annotation', propagate_annotation)
    monkeypatch.setattr(GOEnrichmentRunner, '_get_annotation_iterator', get_annotation_iter)

    dag = ontology.DAGTree.__new__(ontology.DAGTree)
    dag.alt_ids = {}
    dag.go_terms = {'go_id1': ontology.GOTerm(), 'go_id2': ontology.GOTerm(), 'go_id3': ontology.GOTerm(),
                    'go_id4': ontology.GOTerm()}
    dag.go_terms['go_id1'].set_id('go_id1')
    dag.go_terms['go_id2'].set_id('go_id2')
    dag.go_terms['go_id3'].set_id('go_id3')
    dag.go_terms['go_id4'].set_id('go_id4')

    runner = GOEnrichmentRunner.__new__(GOEnrichmentRunner)
    runner.propagate_annotations = 'classic'
    runner.organism = 'organism'
    runner.taxon_id = 'taxon_id'
    runner.dag_tree = dag

    annotation_dict, source_dict = runner._process_annotations()

    assert propagated_annotations == propagated_annotations_truth
    assert annotation_dict == annotation_dict_truth
    assert source_dict == source_dict_truth


def test_go_enrichment_runner_go_classic_pvalues_serial(monkeypatch):
    go_ids = ['attr1', 'attr2', 'attr3']

    def validate_input_params_classic_on_batch(self, go_term_batch, mod_df_index=0):
        try:
            chunk_iterator = iter(go_term_batch)
            assert go_ids == list(chunk_iterator)
        except TypeError:
            assert False

        assert mod_df_index == 0

    monkeypatch.setattr(GOEnrichmentRunner, '_go_classic_on_batch', validate_input_params_classic_on_batch)

    runner = GOEnrichmentRunner.__new__(GOEnrichmentRunner)
    runner.attributes = go_ids
    runner._go_classic_pvalues_serial()


def test_go_enrichment_runner_go_classic_pvalues_parallel(monkeypatch):
    go_ids = [f'id_{i}' for i in range(1500)]
    go_ids_batch_truth = [[f'id_{i}' for i in range(1000)], [f'id_{i}' for i in range(1000, 1500)]]

    def validate_input_params_parallel_over_grouping(self, func, go_term_batches, mod_df_inds, progress_bar_desc):
        assert go_term_batches == go_ids_batch_truth
        assert func.__func__ == GOEnrichmentRunner._go_classic_on_batch
        assert isinstance(progress_bar_desc, str)
        mod_df_inds_lst = [i for i in mod_df_inds]
        for ind in mod_df_inds_lst:
            assert ind == 0
        assert len(mod_df_inds_lst) == len(go_ids_batch_truth)

    monkeypatch.setattr(GOEnrichmentRunner, '_parallel_over_grouping', validate_input_params_parallel_over_grouping)
    runner = GOEnrichmentRunner.__new__(GOEnrichmentRunner)
    runner.attributes = go_ids
    runner._go_classic_pvalues_parallel()


def test_go_enrichment_runner_go_elim_pvalues_serial(monkeypatch):
    def validate_input_params_elim_on_aspect(self, go_aspect, progress_bar_desc):
        assert isinstance(progress_bar_desc, str)
        assert go_aspect == 'all'
        return 'success'

    monkeypatch.setattr(GOEnrichmentRunner, '_go_elim_on_aspect', validate_input_params_elim_on_aspect)

    runner = GOEnrichmentRunner.__new__(GOEnrichmentRunner)
    assert runner._go_elim_pvalues_serial() == 'success'


@pytest.mark.parametrize('aspects', [['aspect0', 'aspect1'], ['aspect0'], ['aspect0', 'aspect1', 'aspect2']])
def test_go_enrichment_runner_go_elim_pvalues_parallel(monkeypatch, aspects):
    class FakeDAG:
        def __init__(self, aspect_list):
            self.namespaces = aspect_list

    @joblib.wrap_non_picklable_objects
    def validate_input_params_elim_on_aspect(aspect, mod_df_ind):
        assert mod_df_ind == int(aspect[-1])
        assert aspect in aspects
        return {aspect: 'success'}

    def is_method_of_class(mthd, cls):
        assert cls == GOEnrichmentRunner
        assert mthd == validate_input_params_elim_on_aspect
        return True

    monkeypatch.setattr(GOEnrichmentRunner, '_go_elim_on_aspect', validate_input_params_elim_on_aspect)
    monkeypatch.setattr(validation, 'is_method_of_class', is_method_of_class)

    runner = GOEnrichmentRunner.__new__(GOEnrichmentRunner)
    runner.dag_tree = FakeDAG(aspects)
    runner.parallel_backend = 'sequential'
    assert runner._go_elim_pvalues_parallel() == {aspect: 'success' for aspect in aspects}


def test_go_enrichment_runner_go_weight_pvalues_serial(monkeypatch):
    def validate_input_params_weight_on_aspect(self, go_aspect, mod_df_ind=0, progress_bar_desc=''):
        assert isinstance(progress_bar_desc, str)
        assert go_aspect == 'all'
        assert mod_df_ind == 0
        return 'success'

    monkeypatch.setattr(GOEnrichmentRunner, '_go_weight_on_aspect', validate_input_params_weight_on_aspect)

    runner = GOEnrichmentRunner.__new__(GOEnrichmentRunner)
    assert runner._go_weight_pvalues_serial() == 'success'


@pytest.mark.parametrize('aspects', [['aspect0', 'aspect1'], ['aspect0'], ['aspect0', 'aspect1', 'aspect2']])
def test_go_enrichment_runner_go_weight_pvalues_parallel(monkeypatch, aspects):
    class FakeDAG:
        def __init__(self, aspect_list):
            self.namespaces = aspect_list

    @joblib.wrap_non_picklable_objects
    def validate_input_params_weight_on_aspect(aspect, mod_df_ind):
        assert mod_df_ind == int(aspect[-1])
        assert aspect in aspects
        return {aspect: 'success'}

    def is_method_of_class(mthd, cls):
        assert cls == GOEnrichmentRunner
        assert mthd == validate_input_params_weight_on_aspect
        return True

    monkeypatch.setattr(GOEnrichmentRunner, '_go_weight_on_aspect', validate_input_params_weight_on_aspect)
    monkeypatch.setattr(validation, 'is_method_of_class', is_method_of_class)

    runner = GOEnrichmentRunner.__new__(GOEnrichmentRunner)
    runner.dag_tree = FakeDAG(aspects)
    runner.parallel_backend = 'sequential'
    assert runner._go_weight_pvalues_parallel() == {aspect: 'success' for aspect in aspects}


def test_go_enrichment_runner_go_allm_pvalues_serial(monkeypatch):
    methods_called = {}
    methods_called_truth = {'classic': True, 'elim': True, 'weight': True}
    outputs_truth = {'classic': 'classic', 'elim': 'elim', 'weight': 'weight'}

    def _calculate_enrichment_serial(self: GOEnrichmentRunner):
        methods_called[self.propagate_annotations] = True
        return self.propagate_annotations

    def _calculate_allm(self, outputs):
        assert outputs == outputs_truth

    monkeypatch.setattr(GOEnrichmentRunner, '_calculate_enrichment_serial', _calculate_enrichment_serial)
    monkeypatch.setattr(GOEnrichmentRunner, '_calculate_allm', _calculate_allm)

    runner = GOEnrichmentRunner.__new__(GOEnrichmentRunner)
    runner._go_allm_pvalues_serial()
    assert methods_called == methods_called_truth
    assert runner.propagate_annotations == 'all.m'


def test_go_enrichment_runner_go_allm_pvalues_serial_error_state(monkeypatch):
    def _calculate_enrichment_serial(self: GOEnrichmentRunner):
        raise AssertionError

    monkeypatch.setattr(GOEnrichmentRunner, '_calculate_enrichment_serial', _calculate_enrichment_serial)

    runner = GOEnrichmentRunner.__new__(GOEnrichmentRunner)
    with pytest.raises(AssertionError):
        runner._go_allm_pvalues_serial()
        assert runner.propagate_annotations == 'all.m'


def test_go_enrichment_runner_go_allm_pvalues_parallel(monkeypatch):
    methods_called = {}
    methods_called_truth = {'classic': True, 'elim': True, 'weight': True}
    outputs_truth = {'classic': 'classic', 'elim': 'elim', 'weight': 'weight'}

    def _calculate_enrichment_parallel(self: GOEnrichmentRunner):
        methods_called[self.propagate_annotations] = True
        return self.propagate_annotations

    def _calculate_allm(self, outputs):
        assert outputs == outputs_truth

    monkeypatch.setattr(GOEnrichmentRunner, '_calculate_enrichment_parallel', _calculate_enrichment_parallel)
    monkeypatch.setattr(GOEnrichmentRunner, '_calculate_allm', _calculate_allm)

    runner = GOEnrichmentRunner.__new__(GOEnrichmentRunner)
    runner._go_allm_pvalues_parallel()
    assert methods_called == methods_called_truth
    assert runner.propagate_annotations == 'all.m'


def test_go_enrichment_runner_go_allm_pvalues_parallel_error_state(monkeypatch):
    def _calculate_enrichment_parallel(self: GOEnrichmentRunner):
        raise AssertionError

    monkeypatch.setattr(GOEnrichmentRunner, '_calculate_enrichment_parallel', _calculate_enrichment_parallel)

    runner = GOEnrichmentRunner.__new__(GOEnrichmentRunner)
    with pytest.raises(AssertionError):
        runner._go_allm_pvalues_parallel()
        assert runner.propagate_annotations == 'all.m'


@pytest.mark.parametrize('attrs,output_dict,truth_dict',
                         [(['attr1', 'attr2'],
                           {'classic': {'attr1': ['classic_val1', 'classic_val2', 0.05],
                                        'attr2': ['classic_val3', 'classic_val4', 1]},
                            'elim': {'attr1': ['elim_val1', 'elim_val2', 0.3],
                                     'attr2': ['elim_val3', 'elim_val4', 0.9999]},
                            'weight': {'attr1': ['weight_val1', 'weight_val2', 0.12],
                                       'attr2': ['weight_val3', 'weight_val4', 0.5]}},
                           {'attr1': ('classic_val1', 'classic_val2', 0.12164404),
                            'attr2': ('classic_val3', 'classic_val4', 0.793674068)})])
def test_go_enrichment_runner_calculate_allm(monkeypatch, attrs, output_dict, truth_dict):
    runner = GOEnrichmentRunner.__new__(GOEnrichmentRunner)
    runner.attributes = attrs
    res = runner._calculate_allm(output_dict)
    assert res.keys() == truth_dict.keys()
    for attr in attrs:
        assert res[attr][:-1] == truth_dict[attr][:-1]
        assert np.isclose(res[attr][-1], truth_dict[attr][-1], atol=0)


def test_go_enrichment_runner_go_level_iterator(monkeypatch):
    this_namespace = 'this_namespace'

    class FakeDAG:
        def __init__(self, go_ids_to_yield: list):
            self.go_ids = go_ids_to_yield

        def level_iter(self, namespace):
            assert namespace == this_namespace
            for go_id in self.go_ids:
                yield go_id

    go_ids_by_level = ['ID1', 'ID5', 'ID4', 'ID10', 'ID2']
    go_ids_in_runner = {'ID1', 'ID2', 'ID3', 'ID4', 'ID10'}

    go_ids_by_level_truth = ['ID1', 'ID4', 'ID10', 'ID2']

    runner = GOEnrichmentRunner.__new__(GOEnrichmentRunner)
    runner.attributes_set = go_ids_in_runner
    runner.dag_tree = FakeDAG(go_ids_by_level)
    assert list(runner._go_level_iterator(this_namespace)) == go_ids_by_level_truth


@pytest.mark.parametrize('grouping,inds,truth',
                         [([['a', 'b', 'c'], ['d', 'e', 'f'], ['g', 'h']], [1, 2, 1],
                           {'a': 1, 'b': 1, 'c': 1, 'd': 2, 'e': 2, 'f': 2, 'g': 1, 'h': 1})])
def test_go_enrichment_runner_parallel_over_grouping(monkeypatch, grouping, inds, truth):
    monkeypatch.setattr(validation, 'is_method_of_class', lambda func, obj_type: True)

    def my_func(group, ind):
        return {obj: ind for obj in group}

    runner = GOEnrichmentRunner.__new__(GOEnrichmentRunner)
    runner.parallel_backend = 'loky'
    res = runner._parallel_over_grouping(my_func, grouping, inds)
    assert res == truth


@pytest.mark.parametrize('obs_truth,exp_truth', [(5, 2.3), (7, 92)])
@pytest.mark.parametrize('pval_fwd,pval_rev,escore_fwd,escore_rev,pval_truth,log2escore_truth',
                         [(0.05, 0.7, 2, 0.5, 0.05, 1), (0.6, 0.3, 1.2, 4, 0.3, -2),
                          (np.nan, np.nan, 2, np.inf, 1, -np.inf)])
def test_extract_xlmhg_results(pval_fwd, pval_rev, escore_fwd, escore_rev, pval_truth,
                               log2escore_truth, obs_truth, exp_truth):
    class XLmHGResultObject:
        def __init__(self, pval, escore):
            self.pval = pval
            self.escore = escore
            self.k = obs_truth
            self.cutoff = 1
            self.K = exp_truth
            self.N = 1

    obs, exp, log2escore, pval = XlmhgTest._extract_xlmhg_results(XLmHGResultObject(pval_fwd, escore_fwd),
                                                                  XLmHGResultObject(pval_rev, escore_rev))

    assert pval == pval_truth
    assert log2escore == log2escore_truth
    assert obs == obs_truth
    assert exp == exp_truth


@pytest.mark.parametrize('exclude_unannotated', [True, False])
@pytest.mark.parametrize(
    'single_list,genes,pval_func,background_set,graph_format',
    [(True, np.array(['WBGene1', 'WBGene2'], dtype=str), XlmhgTest(), None, 'none'),
     (False, {'WBGene00000001', 'WBGene00000002'}, PermutationTest(100, 42),
      {'WBGene00000001', 'WBGene00000002', 'EBGene00000003'}, 'pdf')])
def test_kegg_enrichment_runner_api(monkeypatch, single_list, genes, pval_func, background_set, graph_format,
                                    exclude_unannotated):
    monkeypatch.setattr(io, 'get_taxon_and_id_type', lambda *args: (('a', 'b'), 'gene_id_type'))
    KEGGEnrichmentRunner(genes, 'organism', 'gene_id_type', 0.05, True, False, 'fname', False, False, 'set_name',
                         False, pval_func, background_set, exclude_unannotated, single_list, graph_format)


def test_kegg_enrichment_runner_format_results(monkeypatch):
    def mock_correct_multiple_comparisons(self):
        self.results = self.results.with_columns(significant=pl.Series([False, True, True]))

    monkeypatch.setattr(KEGGEnrichmentRunner, '_correct_multiple_comparisons', mock_correct_multiple_comparisons)
    runner = KEGGEnrichmentRunner.__new__(KEGGEnrichmentRunner)
    results_dict = [['name1', 50, 10, 5.5, 2.3, 0.05], ['name2', 17, 0, 3, 0, 1],
                    ['name3', 1, 2, -2, -0.7, 0.04]]
    truth = io.load_table('tests/test_files/kegg_enrichment_runner_format_results_truth.csv')
    runner.en_score_col = 'colName'
    runner.single_set = False
    runner.return_nonsignificant = False
    runner.pathway_names_dict = {'name1': 'desc1', 'name2': 'desc2', 'name3': 'desc3'}

    runner.format_results(results_dict)
    runner.results.sort(pl.first()).equals(truth.sort(pl.first()))


def test_kegg_enrichment_runner_fetch_annotations(monkeypatch):
    monkeypatch.setattr(KEGGEnrichmentRunner, '_get_query_key', lambda self: 'the_query_key')
    monkeypatch.setattr(KEGGEnrichmentRunner, '_generate_annotation_dict', lambda self: ('kegg_df', 'pathway_names'))
    runner = KEGGEnrichmentRunner.__new__(KEGGEnrichmentRunner)
    runner.fetch_annotations()
    assert runner.annotations == 'kegg_df'
    assert runner.KEGG_DF_QUERIES['the_query_key'] == ('kegg_df', 'pathway_names')
    assert runner.pathway_names_dict == 'pathway_names'

    runner.KEGG_DF_QUERIES['the_query_key'] = ('another_kegg_df', 'other_pathway_names')
    runner.fetch_annotations()
    assert runner.annotations == 'another_kegg_df'
    assert runner.pathway_names_dict == 'other_pathway_names'
    assert runner.KEGG_DF_QUERIES['the_query_key'] == ('another_kegg_df', 'other_pathway_names')


def test_kegg_enrichment_runner_process_annotations(monkeypatch):
    annotation_dict_truth = {'go_id1': {'gene_id1', 'gene_id2'}, 'go_id2': {'gene_id1'}, 'go_id3': {'gene_id2'},
                             'go_id4': {'gene_id3'}}

    name_dict_truth = {'go_id1': 'name1', 'go_id2': 'name2', 'go_id3': 'name3', 'go_id4': 'name4'}

    def get_annotation_iter(self):
        iterator = io.KEGGAnnotationIterator.__new__(io.KEGGAnnotationIterator)
        iterator.n_annotations = 4
        return iterator

    def annotation_iter(self):
        annotations = [['go_id1', 'name1', {'gene_id1', 'gene_id2'}],
                       ['go_id2', 'name2', {'gene_id1'}],
                       ['go_id3', 'name3', {'gene_id2'}],
                       ['go_id4', 'name4', {'gene_id3'}]]
        for annotation in annotations:
            yield annotation

    monkeypatch.setattr(io.KEGGAnnotationIterator, 'get_pathway_annotations', annotation_iter)
    monkeypatch.setattr(KEGGEnrichmentRunner, '_get_annotation_iterator', get_annotation_iter)

    runner = KEGGEnrichmentRunner.__new__(KEGGEnrichmentRunner)
    runner.organism = 'organism'
    runner.taxon_id = 'taxon_id'

    annotation_dict, name_dict = runner._process_annotations()

    assert annotation_dict == annotation_dict_truth
    assert name_dict == name_dict_truth


@pytest.mark.parametrize("mapping_dict,truth", [
    ({}, {'GO1': set(), 'GO2': set()}),
    ({'gene1': 'gene1_translated', 'gene3': 'gene3_translated'},
     {'GO1': {'gene1_translated'}, 'GO2': {'gene1_translated', 'gene3_translated'}}),
    ({'gene1': 'gene1_translated', 'gene2': 'gene2_translated', 'gene3': 'gene3_translated'},
     {'GO1': {'gene1_translated', 'gene2_translated'}, 'GO2': {'gene1_translated', 'gene3_translated'}}
     )])
def test_kegg_enrichment_runner_translate_gene_ids(monkeypatch, mapping_dict, truth):
    monkeypatch.setattr(io.GeneIDTranslator, '__init__', lambda *args: None)
    monkeypatch.setattr(io.GeneIDTranslator, 'run', lambda *args: mapping_dict)

    sparse_annotation_dict = {'GO1': {'gene2', 'gene1'}, 'GO2': {'gene3', 'gene1'}}
    runner = KEGGEnrichmentRunner.__new__(KEGGEnrichmentRunner)
    runner.gene_id_type = 'gene_id_type'

    res = runner._translate_gene_ids(sparse_annotation_dict)
    assert res == truth


def test_kegg_enrichment_runner_get_query_key():
    runner = KEGGEnrichmentRunner.__new__(KEGGEnrichmentRunner)
    runner.taxon_id = 'taxon_id'
    runner.gene_id_type = 'id_type'

    truth = ('taxon_id', 'id_type')
    key = runner._get_query_key()
    assert key == truth
    try:
        _ = hash(key)
    except TypeError:
        assert False


def test_kegg_enrichment_runner_fetch_attributes():
    runner = KEGGEnrichmentRunner.__new__(KEGGEnrichmentRunner)
    runner.annotations = _df_to_dict(io.load_table('tests/test_files/attr_ref_table_for_tests.csv'))
    truth_attributes = ['attribute1', 'attribute2', 'attribute3', 'attribute4']
    truth_attributes_set = {'attribute1', 'attribute2', 'attribute3', 'attribute4'}
    runner.get_attribute_names()
    assert runner.attributes == truth_attributes
    assert runner.attributes_set == truth_attributes_set


@pytest.mark.parametrize('single_set,graph_format', [(True, 'png'), (False, 'none')])
def test_kegg_enrichment_runner_pathway_plot(single_set, graph_format, monkeypatch):
    pathway_id_truth = 'pth_id'
    translator_truth = 'translator obj'
    plotted = []

    class MockPathway:
        def __init__(self, pathway_id, translator):
            assert pathway_id == pathway_id_truth
            assert translator == translator_truth

        def plot_pathway(self, significant, ylabel, graph_format):
            assert significant == runner.gene_set
            if single_set:
                assert ylabel == KEGGEnrichmentRunner.SINGLE_SET_ENRICHMENT_SCORE_YLABEL
            else:
                assert ylabel == KEGGEnrichmentRunner.ENRICHMENT_SCORE_YLABEL
            assert graph_format == graph_format
            plotted.append(True)

    monkeypatch.setattr(ontology, 'fetch_kegg_pathway', lambda *args, **kwargs: MockPathway(*args, **kwargs))
    monkeypatch.setattr(BarPlotter, 'run', lambda *args: list())

    runner = KEGGEnrichmentRunner.__new__(KEGGEnrichmentRunner)
    runner.gene_set = {'a', 'b', 'c'}
    runner.single_set = single_set
    runner.pathway_graphs_format = graph_format
    runner.gene_id_map = translator_truth
    runner.alpha = 0.05
    runner.set_name = 'set name'
    runner.pathway_graphs_format = 'none'
    runner.results = pl.DataFrame({'column_0': [pathway_id_truth], 'padj': [True]})
    runner.en_score_col = 'column_0'
    runner.plot_pathway_graphs = True
    runner.pathway_names_dict = {pathway_id_truth: 'name'}
    runner.plot_results().run()
