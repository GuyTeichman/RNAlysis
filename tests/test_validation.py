import pytest

from rnalysis.utils.validation import *


class DummyClass:
    def __init__(self):
        pass


class DummyClassChild(DummyClass):
    def mthd(self):
        pass


class DummyClassNotChild(dict):
    pass


@pytest.mark.parametrize('obj,cls,expected', [
    ([2, 3], list, True),
    ([2, 3], dict, False),
    ([2, 3], (dict, list), True),
    ([2, 3], (dict, tuple), False),
    (DummyClass(), DummyClass, True),
    (DummyClass(), DummyClassChild, False),
    (DummyClassChild(), DummyClass, True),
    (DummyClassChild(), (DummyClass, int), True),
    (DummyClassChild(), DummyClassChild, True),
    (DummyClassNotChild(), DummyClass, False),
    (DummyClassNotChild(), (DummyClass, DummyClassChild), False)])
def test_isinstanceinh(obj, cls, expected):
    assert isinstanceinh(obj, cls) == expected


def test_is_df_dataframe():
    my_df = pl.DataFrame()
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


@pytest.mark.parametrize("from_dict,to_dict,from_names,to_names,is_legal", [
    ({'one': 1, 'two': 2, 'three': 3}, {'one': 1, 'two': 2, 'three': 3}, ['one', 'two'], ['three'], True),
    ({'one': 1, 'two': 2, 'three': 3}, {'one': 1, 'two': 2, 'three': 3}, ['one'], [2, 'three'], False),
    ({'one': 1, 'two': 2, 'three': 3}, {'one': 1, 'two': 2, 'three': 3}, ['one', 'two', 'three'],
     ['one', 'two', 'three'], True),
    ({'one': 1, 'two': 2, 'three': 3}, {'one': 1, 'two': 2, 'three': 3}, ['one', 'Two', 'three'],
     ['one', 'two', 'three'], False),
    ({'one': 1, 'two': 2, 'three': 3}, {'two': 2, 'three': 3}, ['one', 'two', 'three'],
     ['one', 'two', 'three'], False),
    ({'one': 1, 'three': 3}, {'two': 2, 'three': 3}, ['one', 'two', 'three'],
     ['two', 'three'], False),
])
def test_validate_uniprot_dataset_name(from_dict, to_dict, from_names, to_names, is_legal):
    if is_legal:
        validate_uniprot_dataset_name((to_dict, from_dict), to_names, from_names)
    else:
        with pytest.raises(AssertionError):
            validate_uniprot_dataset_name((to_dict, from_dict), to_names, from_names)


@pytest.mark.parametrize('min_cluster_size,metric,cluster_selection_method,n_features,expected_to_pass',
                         [(1, 'euclidean', 'eom', 1, False),
                          (2.0, 'euclidean', 'eom', 13, False),
                          (14, 'euclidean', 'eom', 13, False),
                          (13, 'euclidean', 'eom', 13, True),
                          (2, 'euclidean', 5, 13, False),
                          (2, 5, 'EOM', 13, False)])
def test_validate_hdbscan_parameters(min_cluster_size, metric, cluster_selection_method, n_features, expected_to_pass):
    if expected_to_pass:
        validate_hdbscan_parameters(min_cluster_size, metric, cluster_selection_method, n_features)
    else:
        with pytest.raises(AssertionError):
            validate_hdbscan_parameters(min_cluster_size, metric, cluster_selection_method, n_features)


@pytest.mark.parametrize("test_input,expected_type,expected", [
    (["3", "hello world", "", ], str, True),
    ({5, 3, 7, 8}, int, True),
    ({-0.1, 3, 7.2, 8, 0}, (int, float), True),
    ((5, 3, '7', 8), int, False),
    ((5, 3, [7, '8'], 8), (int, list), True),
    ([[], [1], [1, 2]], int, False),
    ([], DummyClass, True),
    ({'one': DummyClass(), 'two': DummyClass()}.values(), DummyClass, True)
])
def test_isinstanceiter(test_input, expected_type, expected):
    assert isinstanceiter(test_input, expected_type) == expected


@pytest.mark.parametrize("test_input,expected_type,expected", [
    (["3", "hello world", "", ], str, True),
    ({5, 3, '7', 8}, int, True),
    ({-0.1, 3, 7.2, 8, 0}, (int, float), True),
    ((5, 3, '7', 8), float, False),
    ((5, 3, [7, '8'], 8), list, True),
    ([[], [1], [1, 2]], int, False),
    ([1, [], [1], [1, 2]], (str, int), True),
    ([], DummyClass, False),
    ({'one': DummyClass(), 'two': DummyClass()}.values(), DummyClass, True)
])
def test_isinstanceiter_any(test_input, expected_type, expected):
    assert isinstanceiter_any(test_input, expected_type) == expected


class TestClass:
    def cls_method(self):
        pass

    def other_cls_method(self):
        pass


@pytest.mark.parametrize('mthd,cls,truth',
                         [(str.lower, str, True),
                          (str.lower, int, False),
                          (lambda x: x + 1, str, False),
                          (TestClass.cls_method, TestClass, True),
                          (TestClass().cls_method, TestClass, True),
                          (str.lower, TestClass, False)])
def test_is_method_of_class(mthd, cls, truth):
    assert is_method_of_class(mthd, cls) == truth


@pytest.mark.parametrize('metric,linkage,expected_to_pass', [
    ('euclidean', 'single', True),
    (1, 'single', False),
    ('euclidean', 1, False),
    ('Spearman', 'COMPLETE', True),
    ('sportman', 'complete', False),
    ('spearman', 'compete', False),
    ('jackknife', None, True),
    (None, None, False),
    (None, 'single', False),
    ('sportman', None, False)])
def test_validate_clustering_parameters(metric, linkage, expected_to_pass):
    legal_metrics = {'euclidean', 'spearman', 'jackknife', 'other'}
    if expected_to_pass:
        truth = (metric.lower(), linkage.lower()) if linkage is not None else metric.lower()
        assert validate_clustering_parameters(legal_metrics, metric, linkage) == truth
    else:
        with pytest.raises(AssertionError):
            validate_clustering_parameters(legal_metrics, metric, linkage)


@pytest.mark.parametrize('table_path,expected_to_pass',
                         [('attr_ref_table_for_tests.csv', True),
                          ('attr_ref_table_one_col.csv', False),
                          ('attr_ref_table_one_row.csv', False),
                          ('attr_ref_table_no_index_title.csv', True),
                          ('attr_ref_table_wrong_index_title.csv', True)])
def test_validate_attr_table(table_path, expected_to_pass):
    df = pl.read_csv('tests/test_files/' + table_path)
    if expected_to_pass:
        validate_attr_table(df)
    else:
        with pytest.raises(AssertionError):
            validate_attr_table(df)


@pytest.mark.parametrize('table_path,expected_to_pass',
                         [('biotype_ref_table_for_tests.csv', True),
                          ('biotype_ref_table_one_col.csv', False),
                          ('biotype_ref_table_too_many_cols.csv', False),
                          ('biotype_ref_table_one_row.csv', False),
                          ('biotype_ref_table_no_titles.csv', True),
                          ('biotype_ref_table_wrong_titles.csv', True)])
def test_validate_biotype_table(table_path, expected_to_pass):
    df = pl.read_csv('tests/test_files/' + table_path)
    if expected_to_pass:
        validate_biotype_table(df)
        assert df.columns[0] == 'gene'
        assert df.columns[1] == 'biotype'
        assert len(df.columns) == 2
    else:
        with pytest.raises(AssertionError):
            validate_biotype_table(df)


@pytest.mark.parametrize('threshold,expected_to_pass', [
    (5, True),
    (0.1, True),
    ('5', False),
    ([5, 6], False),
    (0, True),
    (0.0, True),
    (-2, False),
    (-3.14, False)])
def test_validate_threshold(threshold, expected_to_pass):
    if expected_to_pass:
        validate_threshold(threshold)
    else:
        with pytest.raises(AssertionError):
            validate_threshold(threshold)


@pytest.mark.parametrize('path,is_legal_truth', [
    ('tests/test_files/test_deseq.csv', True),
    ('test_deseq.csv', False),
    ('tests/test_files', False),
    ('tests/test_files/test_gene_set.txt', True)
])
def test_is_legal_file_path(path, is_legal_truth):
    assert is_legal_file_path(path) == is_legal_truth


@pytest.mark.parametrize('path,is_legal_truth', [
    ('tests', True),
    ('test_deseq.csv', False),
    ('tests/test_files/test_count_from_folder', True),
    ('tests/test_files/fake_dir', False)
])
def test_is_legal_dir_path(path, is_legal_truth):
    assert is_legal_dir_path(path) == is_legal_truth


@pytest.mark.parametrize('path,accept_gtf,accept_gff3,err_exp,truth', [
    ('tests/test_files/counted.csv', True, True, True, None),
    ('tests/test_files/kallisto_tests/transcripts.gtf', True, True, False, 'gtf'),
    ('tests/test_files/doesnt_exist.gtf', True, True, True, None),
    ('tests/test_files/test_gff_wormbase.gff3', True, True, False, 'gff3'),
    ('tests/test_files/kallisto_tests/transcripts.gtf', True, False, False, 'gtf'),
    ('tests/test_files/kallisto_tests/transcripts.gtf', False, True, True, None),
    ('tests/test_files/test_gff_wormbase.gff3', False, True, False, 'gff3'),
    ('tests/test_files/test_gff_wormbase.gff3', True, False, True, None),

])
def test_validate_genome_annotation_file(path, accept_gtf, accept_gff3, err_exp, truth):
    if err_exp:
        with pytest.raises((AssertionError, ValueError)):
            validate_genome_annotation_file(path, accept_gtf, accept_gff3)
    else:
        res = validate_genome_annotation_file(path, accept_gtf, accept_gff3)
        assert res == truth
