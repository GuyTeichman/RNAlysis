import pytest
import numpy as np

from rnalysis.utils.clustering import *
from rnalysis.utils.parallel import *
from rnalysis.utils.validation import *
from rnalysis.utils.preprocessing import *
from rnalysis.utils.ref_tables import *
from rnalysis.utils.parsing import *
from rnalysis.utils.io import *
from rnalysis.utils.io import _format_ids_iter
from rnalysis import __attr_file_key__, __biotype_file_key__


class MockResponse(object):
    def __init__(self, status_code: int = 200, url: str = 'http://httpbin.org/get', headers: dict = {'blaa': '1234'},
                 text: str = ''):
        self.status_code = status_code
        self.url = url
        self.headers = headers
        self.text = text
        self.ok = self.status_code == 200

    def raise_for_status(self):
        if not self.ok:
            raise ConnectionError('request not ok')


def test_is_df_dataframe():
    my_df = pd.DataFrame()
    assert (check_is_df_like(my_df))


def test_is_df_str():
    correct_path = "myfile.csv"
    correct_path_2 = r"myfolder\anotherfile.csv"
    assert (not check_is_df_like(correct_path))
    assert (not check_is_df_like(correct_path_2))


def test_is_df_pathobj():
    correct_pathobj = Path("all_feature_96_new.csv")
    assert (not check_is_df_like(correct_pathobj))


def test_is_df_str_notcsv():
    incorrect_pth = "myfile.xlsx"
    incorrect_pth2 = "all_feature_96_new"
    with pytest.raises(ValueError):
        check_is_df_like(incorrect_pth)
    with pytest.raises(ValueError):
        check_is_df_like(incorrect_pth2)


def test_is_df_pathobj_notcsv():
    incorrect_pathobj = Path("test_general.py")
    with pytest.raises(ValueError):
        check_is_df_like(incorrect_pathobj)


def test_is_df_invalid_type():
    invalid_type = 67
    with pytest.raises(ValueError):
        check_is_df_like(invalid_type)


def test_load_csv_bad_input():
    invalid_input = 2349
    with pytest.raises(AssertionError):
        a = load_csv(invalid_input)


def test_load_csv():
    truth = pd.DataFrame({'idxcol': ['one', 'two', 'three'], 'othercol': [4, 5, 6]})
    truth.set_index('idxcol', inplace=True)
    pth = "test_files/test_load_csv.csv"
    assert (truth == load_csv(pth, 0)).all().all()


def test_load_csv_drop_columns():
    loaded = load_csv('test_files/counted.csv', 0, drop_columns='cond1')
    print(loaded)
    assert list(loaded.columns) == ['cond2', 'cond3', 'cond4']

    loaded = load_csv('test_files/counted.csv', 0, drop_columns=['cond2', 'cond4'])
    assert list(loaded.columns) == ['cond1', 'cond3']

    with pytest.raises(IndexError):
        load_csv('test_files/counted.csv', 0, drop_columns=['cond1', 'cond6'])


def test_biotype_table_assertions():
    assert False


def test_attr_table_assertions():
    assert False


def test_get_settings_file_path():
    make_temp_copy_of_settings_file()

    set_temp_copy_of_settings_file_as_default()
    remove_temp_copy_of_settings_file()
    assert False


def test_get_attr_ref_path():
    make_temp_copy_of_settings_file()
    update_settings_file('path/to/attr/ref/file', __attr_file_key__)
    with get_settings_file_path().open() as f:
        success = False
        for line in f.readlines():
            if line.startswith(__attr_file_key__):
                success = line.startswith('attribute_reference_table: path/to/attr/ref/file')
                if not success:
                    print(f'failiure at: {line}')

    set_temp_copy_of_settings_file_as_default()
    remove_temp_copy_of_settings_file()
    assert success


def test_get_biotype_ref_path():
    make_temp_copy_of_settings_file()
    update_settings_file('path/to/biotype/ref/file', __biotype_file_key__)
    with get_settings_file_path().open() as f:
        success = False
        for line in f.readlines():
            if line.startswith(__biotype_file_key__):
                success = line.startswith('biotype_reference_table: path/to/biotype/ref/file')
                if not success:
                    print(f'failiure at: {line}')

    set_temp_copy_of_settings_file_as_default()
    remove_temp_copy_of_settings_file()
    assert success


def test_update_ref_table_attributes():
    make_temp_copy_of_settings_file()
    try:
        get_settings_file_path().unlink()
    except FileNotFoundError:
        pass
    update_settings_file('path/to/biotype/ref/file', __biotype_file_key__)
    with get_settings_file_path().open() as f:
        counter = 0
        success = False
        for line in f.readlines():
            if line.startswith(__biotype_file_key__):
                counter += 1
                if counter > 1:
                    print(f'Failure with counter={counter}')
                success = line.startswith('biotype_reference_table: path/to/biotype/ref/file')
                if not success:
                    print(f'Failure at: {line}')

    update_settings_file('a/new/path/to/biotype', __biotype_file_key__)
    with get_settings_file_path().open() as f:
        counter_2 = 0
        success_2 = False
        for line in f.readlines():
            if line.startswith(__biotype_file_key__):
                counter_2 += 1
                if counter_2 > 1:
                    print(f'Failure with counter={counter_2}')
                success_2 = line.startswith('biotype_reference_table: a/new/path/to/biotype')
                if not success_2:
                    print(f'Failure at: {line}')

    set_temp_copy_of_settings_file_as_default()
    remove_temp_copy_of_settings_file()
    assert success
    assert counter == 1
    assert counter_2 == 1


def test_save_csv():
    try:
        df = pd.read_csv('test_files/enrichment_hypergeometric_res.csv', index_col=0)
        save_csv(df, 'test_files/tmp_test_save_csv.csv')
        df_loaded = pd.read_csv('test_files/tmp_test_save_csv.csv', index_col=0)
        assert df.equals(df_loaded)
        df = pd.read_csv('test_files/enrichment_hypergeometric_res.csv')
        save_csv(df, 'test_files/tmp_test_save_csv.csv', '_2', index=False)
        df_loaded = pd.read_csv('test_files/tmp_test_save_csv_2.csv', index_col=0)
        df = pd.read_csv('test_files/enrichment_hypergeometric_res.csv', index_col=0)
        assert df.equals(df_loaded)

    except Exception as e:
        raise e
    finally:
        try:
            os.remove('test_files/tmp_test_save_csv.csv')
            os.remove('test_files/tmp_test_save_csv_2.csv')
        except:
            pass


def test_standardize():
    np.random.seed(42)
    data = np.random.randint(-200, 100000, (100, 5))
    res = standardize(data)
    assert res.shape == data.shape
    assert np.isclose(res.mean(axis=0), 0).all()
    assert np.isclose(res.std(axis=0), 1).all()


def test_standard_box_cos():
    np.random.seed(42)
    data = np.random.randint(-200, 100000, (100, 5))
    res = standard_box_cox(data)
    assert res.shape == data.shape
    assert np.isclose(res.mean(axis=0), 0).all()
    assert np.isclose(res.std(axis=0), 1).all()
    assert not np.isclose(res, standardize(data)).all()


def test_kmedoidsiter_api():
    truth = KMedoids(3, max_iter=300, init='k-medoids++', random_state=42)
    kmeds = KMedoidsIter(3, max_iter=300, init='k-medoids++', n_init=1, random_state=42)
    df = load_csv('test_files/counted.csv', 0)
    truth.fit(df)
    kmeds.fit(df)
    assert np.all(truth.cluster_centers_ == kmeds.cluster_centers_)
    assert np.all(truth.inertia_ == kmeds.inertia_)

    assert np.all(truth.predict(df) == kmeds.predict(df))
    assert np.all(truth.fit_predict(df) == kmeds.fit_predict(df))

    kmeds_rand = KMedoidsIter(3, max_iter=300, init='k-medoids++', n_init=3)
    kmeds_rand.fit(df)
    kmeds_rand.predict(df)
    kmeds_rand.fit_predict(df)

    assert repr(kmeds) == repr(kmeds.clusterer)
    assert str(kmeds) == str(kmeds.clusterer)


def test_kmedoidsiter_iter():
    kmeds = KMedoidsIter(3, max_iter=300, init='k-medoids++', n_init=5, random_state=0)
    df = load_csv('test_files/counted.csv', 0)
    kmeds.fit(df)

    inertias = []
    clusterers = []
    for i in range(5):
        clusterers.append(KMedoids(3, max_iter=300, init='k-medoids++', random_state=0).fit(df))
        inertias.append(clusterers[i].inertia_)
    truth_inertia = max(inertias)
    truth_kmeds = clusterers[np.argmax(inertias)]
    assert kmeds.inertia_ == truth_inertia
    assert np.all(kmeds.clusterer.predict(df) == truth_kmeds.predict(df))


def test_parse_go_id():
    line = b"is_a: 123456 GO1234567 GO:123456 GO:7654321! primary alcohol metabolic process"
    truth = "GO:7654321"
    assert parsing.parse_go_id(line) == truth
    line_2 = b"is_a: GO:1234567 GO:0034308 ! primary alcohol metabolic process"
    truth_2 = "GO:1234567"
    assert parsing.parse_go_id(line_2) == truth_2


def test_dag_tree_parser_api():
    file = 'test_files/go_mini.obo'
    with open(file, 'rb') as f:
        dag_tree = parsing.DAGTreeParser(f, ['is_a', 'part_of', 'regulates'])
    with open(file, 'rb') as f:
        dag_tree._parse_file(f)


def test_dag_tree_parser_construction():
    file = 'test_files/go_mini.obo'
    levels_truth = [{'GO:0034308': 0, 'GO:0051125': 0},
                    {'GO:0034315': 0, 'GO:0006040': 0, 'GO:0006793': 0, 'GO:0009225': 0, 'GO:0009226': 0},
                    {'GO:2001315': 0, 'GO:2001313': 0}]
    data_version_truth = 'releases/2020-07-16'
    with open(file, 'rb') as f:
        dag_tree = parsing.DAGTreeParser(f, ['is_a', 'part_of', 'regulates'])

    assert dag_tree.data_version == data_version_truth
    assert dag_tree.parent_relationship_types == ['is_a', 'part_of', 'regulates']

    assert dag_tree.alt_ids == {'GO:0034619': 'GO:0006040'}
    assert dag_tree['GO:0034619'] == dag_tree['GO:0006040']

    assert len(dag_tree.levels) == len(levels_truth)
    for level, truth_level in zip(dag_tree.levels, levels_truth):
        assert level.keys() == truth_level.keys()
        for item in truth_level:
            assert item in dag_tree

    assert 'GO:0003840' not in dag_tree


def test_go_term_get_parents():
    file = 'test_files/go_mini.obo'
    with open(file, 'rb') as f:
        dag_tree = parsing.DAGTreeParser(f, ['is_a', 'part_of', 'regulates'])
    parents_truth = {'GO:0034308': [], 'GO:0051125': [], 'GO:0034315': ['GO:0051125'], 'GO:0006040': ['GO:0034308'],
                     'GO:0006793': ['GO:0034308'], 'GO:0009225': ['GO:0034308'], 'GO:0009226': ['GO:0034308'],
                     'GO:2001315': ['GO:0009226', 'GO:0006793'],
                     'GO:2001313': ['GO:0009225', 'GO:0006793', 'GO:0006040']}
    for node in dag_tree.go_terms:
        assert node in parents_truth
        assert dag_tree[node].get_parents(dag_tree.parent_relationship_types).sort() == parents_truth[node].sort()


def test_go_term_get_children():
    file = 'test_files/go_mini.obo'
    with open(file, 'rb') as f:
        dag_tree = parsing.DAGTreeParser(f, ['is_a', 'part_of', 'regulates'])
    children_truth = {'GO:0034308': ['GO:0006040', 'GO:0006793', 'GO:0009225', 'GO:0009226'],
                      'GO:0051125': ['GO:0034315'], 'GO:0034315': [], 'GO:0006040': ['GO:2001313'],
                      'GO:0006793': ['GO:2001313', 'GO:2001315'], 'GO:0009225': ['GO:2001313'],
                      'GO:0009226': ['GO:2001315'], 'GO:2001315': [], 'GO:2001313': []}
    for node in dag_tree.go_terms:
        assert node in children_truth
        assert dag_tree[node].get_children(dag_tree.parent_relationship_types).sort() == children_truth[node].sort()


def test_dag_tree_parser_level_iterator():
    levels_truth = [{'GO:0034308', 'GO:0051125'},
                    {'GO:0034315', 'GO:0006040', 'GO:0006793', 'GO:0009225', 'GO:0009226'},
                    {'GO:2001315', 'GO:2001313'}]
    levels_truth.reverse()
    file = 'test_files/go_mini.obo'
    with open(file, 'rb') as f:
        dag_tree = parsing.DAGTreeParser(f, ['is_a', 'part_of', 'regulates'])

    level_iter = dag_tree.level_iter()
    for truth_level in levels_truth:
        level = set()
        for i in range(len(truth_level)):
            level.add(next(level_iter))
        assert truth_level == level


def test_dag_tree_parser_upper_induced_tree_iterator():
    parents_truth = {'GO:0034308': [], 'GO:0051125': [], 'GO:0034315': ['GO:0051125'], 'GO:0006040': ['GO:0034308'],
                     'GO:0006793': ['GO:0034308'], 'GO:0009225': ['GO:0034308'], 'GO:0009226': ['GO:0034308'],
                     'GO:2001315': ['GO:0009226', 'GO:0006793', 'GO:0034308'],
                     'GO:2001313': ['GO:0009225', 'GO:0006793', 'GO:0006040', 'GO:0034308']}
    file = 'test_files/go_mini.obo'
    with open(file, 'rb') as f:
        dag_tree = parsing.DAGTreeParser(f, ['is_a', 'part_of', 'regulates'])

    for node in parents_truth:
        parents_truth[node].sort()
        ui_tree = list(dag_tree.upper_induced_graph_iter(node))
        ui_tree.sort()
        assert ui_tree == parents_truth[node]


def test_data_to_list():
    assert parsing.data_to_list([1, 2, 'hi']) == [1, 2, 'hi']
    assert parsing.data_to_list((1, 2, 'hi')) == [1, 2, 'hi']
    assert parsing.data_to_list('fifty seven brave men') == ['fifty seven brave men']
    assert sorted(parsing.data_to_list({'three', 'different', 'elements'})) == sorted(
        ['three', 'different', 'elements'])
    assert parsing.data_to_list(np.array([6, 9, 2])) == [6, 9, 2]
    assert parsing.data_to_list(67.2) == [67.2]


def test_data_to_set():
    assert parsing.data_to_set([1, 2, 'hi']) == {1, 2, 'hi'}
    assert parsing.data_to_set((1, 2, 'hi')) == {1, 2, 'hi'}
    assert parsing.data_to_set('fifty seven brave men') == {'fifty seven brave men'}
    assert parsing.data_to_set({'three', 'different', 'elements'}) == {'three', 'different', 'elements'}
    assert parsing.data_to_set(np.array([6, 9, 2])) == {6, 9, 2}
    assert parsing.data_to_set(67.2) == {67.2}


def test_data_to_tuple():
    assert parsing.data_to_tuple([1, 2, 'hi']) == (1, 2, 'hi')
    assert parsing.data_to_tuple('fifty seven brave men') == ('fifty seven brave men',)
    assert sorted(parsing.data_to_tuple({'three', 'different', 'elements'})) == sorted(
        ('three', 'different', 'elements'))
    assert parsing.data_to_tuple(np.array([6, 9, 2])) == (6, 9, 2)
    assert parsing.data_to_tuple(67.2) == (67.2,)


def test_from_string(monkeypatch):
    monkeypatch.setattr('builtins.input', lambda x: 'one\ntwo \nthree;\n')
    assert from_string() == ['one', 'two ', 'three;']


def test_intersection_nonempty():
    assert intersection_nonempty({1, 2, 3}, {2, 5, 7, 1}) == {1, 2}
    assert intersection_nonempty({1, 3, 7}, set()) == {1, 3, 7}
    assert intersection_nonempty({1, 2, 4, 5}, set(), {1, 3, 4, 7}) == {1, 4}


def test_fetch_go_basic_connectivity():
    _ = fetch_go_basic()


def test_map_gene_ids_connectivity():
    ids_uniprot = ['P34544', 'Q27395', 'P12844']
    ids_wormbase = ['WBGene00019883', 'WBGene00023497', 'WBGene00003515']
    entrez_to_wb_truth = {'176183': 'WBGene00019883', '173203': 'WBGene00012343'}
    wb_to_entrez_truth = {val: key for key, val in zip(entrez_to_wb_truth.keys(), entrez_to_wb_truth.values())}
    mapped_ids_truth = {uniprot: wb for uniprot, wb in zip(ids_uniprot, ids_wormbase)}
    mapped_ids_truth_rev = {b: a for a, b in zip(mapped_ids_truth.keys(), mapped_ids_truth.values())}

    mapped_ids = map_gene_ids(ids_uniprot, map_from='UniProtKB', map_to='WormBase')
    for geneid in ids_uniprot:
        assert geneid in mapped_ids
        assert mapped_ids[geneid] == mapped_ids_truth[geneid]
    assert mapped_ids.mapping_dict == mapped_ids_truth

    mapped_ids = map_gene_ids(ids_wormbase, map_from='WormBase', map_to='UniProtKB')
    for geneid in ids_wormbase:
        assert geneid in mapped_ids
        assert mapped_ids[geneid] == mapped_ids_truth_rev[geneid]
    assert mapped_ids.mapping_dict == mapped_ids_truth_rev

    mapped_ids = map_gene_ids(entrez_to_wb_truth.keys(), 'Entrez Gene ID', 'WormBase')
    for geneid in entrez_to_wb_truth:
        assert mapped_ids[geneid] == entrez_to_wb_truth[geneid]

    mapped_ids = map_gene_ids(wb_to_entrez_truth.keys(), 'WormBase', 'Entrez Gene ID')
    for geneid in wb_to_entrez_truth:
        assert mapped_ids[geneid] == wb_to_entrez_truth[geneid]


def test_map_gene_ids_request(monkeypatch):
    def mock_get(url):
        return MockResponse()

    monkeypatch.setattr(requests, 'get', mock_get)
    assert False


def test_map_gene_ids_to_same_set():
    mapper = map_gene_ids(['it', 'doesnt', 'matter', 'what', 'is', 'in', 'here'], 'UniProtKB', 'UniProtKB')
    assert mapper.mapping_dict is None
    for i in ['it', 'not', False, 42, 3.14]:
        assert i in mapper
        assert mapper[i] == i


def test_map_gene_ids_parsing(monkeypatch):
    def mock_get(url, params):
        return MockResponse(text='From\tTo\na\tA\nQ27395\tb\nB\tWBGene00003515\n')

    monkeypatch.setattr(requests, 'get', mock_get)
    assert False


def test_format_ids_iter():
    assert list(_format_ids_iter('one two three')) == ['one two three']
    assert list(_format_ids_iter(123)) == ['123']
    assert list(_format_ids_iter(['one', ' two', 'three; ', 'four'])) == ['one  two three;  four']
    assert list(_format_ids_iter(['1', 'two', '3', '4', 'five', '6', '7'], 3)) == ['1 two 3', '4 five 6', '7']


def _load_id_abbreviation_dict():
    assert False


def test_gene_id_translator_api():
    _ = GeneIDTranslator({1: 2, 3: 4})
    _ = GeneIDTranslator()


def test_gene_id_translator_getitem():
    translator = GeneIDTranslator({1: 2, 3: 4})
    assert translator[1] == 2
    assert translator[3] == 4
    translator = GeneIDTranslator(None)
    for something in [2, 3, '1', False, True, {}, 3.141592]:
        assert translator[something] == something


def test_gene_id_translator_contains():
    translator = GeneIDTranslator({1: 2, 3: 4})
    assert 1 in translator and 3 in translator
    for invalid in [2, 4, '1', False]:
        assert invalid not in translator
    translator = GeneIDTranslator(None)
    for something in [2, 3, '1', False, True, {}, 3.141592]:
        assert something in translator


def test_map_taxon_id_connectivity():
    assert map_taxon_id(6239) == (6239, 'Caenorhabditis elegans')
    assert map_taxon_id('canis lupus familiaris') == (9615, 'Canis lupus familiaris')
    with pytest.raises(ValueError):
        map_taxon_id('Lorem ipsum dolor sit amet')


def test_map_taxon_id_parsing(monkeypatch):
    assert False


def test_map_taxon_id_no_results(monkeypatch):
    assert False


def test_map_taxon_id_multiple_results(monkeypatch):
    assert False


def test_map_taxon_id_no_connection(monkeypatch):
    assert False


def test_golr_annotation_iterator_connectivity():
    assert False


def test_golr_annotation_iterator_parsing(monkeypatch):
    assert False


def test_stop_ipcluster():
    stop_ipcluster()


def test_start_ipcluster():
    start_ipcluster(1)


def test_is_iterable():
    assert isiterable(range(10))
    assert isiterable([])
    assert isiterable('fortytwo')
    assert not isiterable(42)
    assert not isiterable(3.14)
    assert not isiterable(True)
    assert isiterable((i for i in range(10)))
    assert isiterable({})
    assert isiterable({1, 2, 4})


def test_is_iterable_not_emptying_generator():
    gen = (i for i in range(10) if i != 7)
    assert isiterable(gen)
    assert list(gen) == [i for i in range(10) if i != 7]


def test_validate_uniprot_dataset_name():
    validate_uniprot_dataset_name({'one': 1, 'two': 2, 'three': 3}, 'one', 'two', 'three')
    with pytest.raises(AssertionError):
        validate_uniprot_dataset_name({'one': 1, 'two': 2, 'three': 3}, 'one', 2, 'three')
    with pytest.raises(AssertionError):
        validate_uniprot_dataset_name({'one': 1, 'two': 2, 'three': 3}, 'one', 'Two', 'three')


def test_validate_hdbscan_parameters():
    with pytest.raises(AssertionError):
        validate_hdbscan_parameters(1, 'euclidean', 'eom', 1)

    with pytest.raises(AssertionError):
        validate_hdbscan_parameters(2.0, 'euclidean', 'eom', 13)

    with pytest.raises(AssertionError):
        validate_hdbscan_parameters(14, 'euclidean', 'eom', 13)

    validate_hdbscan_parameters(13, 'euclidean', 'eom', 13)

    with pytest.raises(AssertionError):
        validate_hdbscan_parameters(2, 'euclidean', 5, 13)

    with pytest.raises(AssertionError):
        validate_hdbscan_parameters(2, 5, 'EOM', 13)
