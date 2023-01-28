import pytest

from rnalysis.utils.parsing import *
from rnalysis.utils.parsing import _parse_r_arg


class DummyClass:
    def __init__(self):
        pass


@pytest.mark.parametrize('is_input_path', [True, False])
@pytest.mark.parametrize('path,expected_path', [
    ('file.csv', 'file'),
    ('path/to/file.exe', 'path/to/file'),
    ('name.fastq.gz', 'name'),
    ('path/to/a/filename.sam.tar.gz', 'path/to/a/filename'),
])
def test_remove_suffixes(path, expected_path, is_input_path):
    if is_input_path:
        assert remove_suffixes(Path(path)) == Path(expected_path)
    else:
        assert remove_suffixes(path) == expected_path


@pytest.mark.parametrize('input_val,sort,truth',
                         [([1, 2, 'hi'], False, [1, 2, 'hi']),
                          ((1, 2, 'hi'), False, [1, 2, 'hi']),
                          ('fifty seven brave men', False, ['fifty seven brave men']),
                          ('fifty seven brave men', True, ['fifty seven brave men']),
                          ({'three', 'different', 'elements'}, True, ['different', 'elements', 'three']),
                          (np.array([6, 9, 2]), False, [6, 9, 2]),
                          (np.array([6, 9, 2]), True, [2, 6, 9]),
                          (None, False, [None])])
def test_data_to_list(input_val, sort, truth):
    assert data_to_list(input_val, sort=sort) == truth


@pytest.mark.parametrize('input_val,sort,truth',
                         [([1, 2, 'hi'], False, (1, 2, 'hi')),
                          ((1, 2, 'hi'), False, (1, 2, 'hi')),
                          ('fifty seven brave men', False, ('fifty seven brave men',)),
                          ('fifty seven brave men', True, ('fifty seven brave men',)),
                          ({'three', 'different', 'elements'}, True, ('different', 'elements', 'three')),
                          (np.array([6, 9, 2]), False, (6, 9, 2)),
                          (np.array([6, 9, 2]), True, (2, 6, 9)),
                          (67.2, False, (67.2,)),
                          ((67.2,), False, (67.2,)),
                          ((67.5,), True, (67.5,)),
                          (None, False, (None,))])
def test_data_to_tuple(input_val, sort, truth):
    assert data_to_tuple(input_val, sort=sort) == truth


def test_data_to_list_invalid_type():
    with pytest.raises(TypeError):
        data_to_list(DummyClass())


def test_data_to_tuple_invalid_type():
    with pytest.raises(TypeError):
        data_to_tuple(DummyClass())


def test_data_to_set_invalid_type():
    with pytest.raises(TypeError):
        data_to_set(DummyClass())


def test_data_to_set():
    assert data_to_set([1, 2, 'hi']) == {1, 2, 'hi'}
    assert data_to_set((1, 2, 'hi')) == {1, 2, 'hi'}
    assert data_to_set('fifty seven brave men') == {'fifty seven brave men'}
    assert data_to_set({'three', 'different', 'elements'}) == {'three', 'different', 'elements'}
    assert data_to_set(np.array([6, 9, 2])) == {6, 9, 2}
    assert data_to_set(67.2) == {67.2}
    assert data_to_set(None) == {None}


def test_from_string(monkeypatch):
    monkeypatch.setattr('builtins.input', lambda x: 'one\t\ntwo \nthree; and four\n')
    assert from_string() == ['one\t', 'two ', 'three; and four']
    assert from_string(del_spaces=True) == ['one\t', 'two', 'three;andfour']


def test_uniprot_tab_to_dict():
    tab = 'From\tTo\nWBGene00019883\tP34544\nWBGene00023497\tQ27395\nWBGene00003515\tP12844\nWBGene00000004\t' \
          'A0A0K3AVL7\nWBGene00000004\tO17395\n'
    tab_rev = 'From\tTo\nP34544\tWBGene00019883\nQ27395\tWBGene00023497\nP12844\tWBGene00003515\nA0A0K3AVL7\t' \
              'WBGene00000004\nO17395\tWBGene00000004\n'

    truth = ({'WBGene00019883': 'P34544', 'WBGene00023497': 'Q27395', 'WBGene00003515': 'P12844'},
             ['A0A0K3AVL7', 'O17395'])
    truth_rev = ({'P34544': 'WBGene00019883', 'Q27395': 'WBGene00023497', 'P12844': 'WBGene00003515',
                  'A0A0K3AVL7': 'WBGene00000004', 'O17395': 'WBGene00000004'}, [])
    assert truth == uniprot_tab_to_dict(tab)

    assert truth_rev == uniprot_tab_to_dict(tab_rev)


def test_uniprot_tab_to_dict_empty():
    tab = 'From\tTo\n'
    assert uniprot_tab_to_dict(tab) == ({}, [])


def test_uniprot_tab_with_score_to_dict_empty():
    tab = 'Entry\tAnnotation\tyourlist:M20200816216DA2B77BFBD2E6699CA9B6D1C41EB2A5FE6AF\n'
    assert uniprot_tab_with_score_to_dict(tab) == {}
    assert uniprot_tab_with_score_to_dict(tab, True) == {}


def test_uniprot_tab_with_score_to_dict():
    tab = 'Entry\tAnnotation\tyourlist:M20200816216DA2B77BFBD2E6699CA9B6D1C41EB2A5FE6AF\nP34544\t5 out of 5\t' \
          'WBGene00019883\nQ27395\t4 out of 5\tWBGene00023497\nP12844\t5 out of 5\tWBGene00003515\nA0A0K3AVL7\t' \
          '1 out of 5\tWBGene00000004\nO17395\t2 out of 5\tWBGene00000004\n'

    truth = {'WBGene00019883': 'P34544', 'WBGene00023497': 'Q27395', 'WBGene00003515': 'P12844',
             'WBGene00000004': 'O17395'}
    truth_rev = {'P34544': 'WBGene00019883', 'Q27395': 'WBGene00023497', 'P12844': 'WBGene00003515',
                 'A0A0K3AVL7': 'WBGene00000004', 'O17395': 'WBGene00000004'}
    assert truth == uniprot_tab_with_score_to_dict(tab)

    assert truth_rev == uniprot_tab_with_score_to_dict(tab, True)


def test_sparse_dict_to_bool_df():
    truth = pd.read_csv('tests/test_files/sparse_dict_to_df_truth.csv', index_col=0).sort_index(axis=1)
    sparse_dict = {'gene1': {'a', 'c'}, 'gene2': {'b'}, 'gene3': {'b', 'g', 'h'}, 'gene4': {'e', 'd', 'a', 'c', 'f'},
                   'gene5': set()}
    res = sparse_dict_to_bool_df(sparse_dict).sort_index(axis=1)
    assert res.equals(truth)


@pytest.mark.parametrize("lst,chunk_size,truth", [([1, 2, 3, 4, 5], 1, [[1], [2], [3], [4], [5]]),
                                                  ([1, 2, 3, 4], 2, [[1, 2], [3, 4]]),
                                                  ([1, 2, 3, 4, 5], 2, [[1, 2], [3, 4], [5]]),
                                                  ([], 3, [[]]),
                                                  ((1, 2, 3, 4), 3, [(1, 2, 3), (4,)]),
                                                  ((1, 2, 3, 4, 5), 5, [(1, 2, 3, 4, 5)]),
                                                  ((1, 2, 3, 4), 5, [(1, 2, 3, 4)]),
                                                  (tuple(), 1, [tuple()])])
def test_partition_list(lst, chunk_size, truth):
    res = partition_list(lst, chunk_size)
    assert res == truth


@pytest.mark.parametrize("lst,chunk_size", [([1, 2, 3, 4], 0),
                                            ([1, 2], -1),
                                            ({1, 2}, 2),
                                            ([1, 2, 3], 1.0),
                                            ((i for i in range(3)), 2)])
def test_partition_list_invalid_input(lst, chunk_size):
    with pytest.raises(AssertionError):
        _ = partition_list(lst, chunk_size)


@pytest.mark.parametrize("lst,truth", [([1, 2, 3, 4, 5], [1, 2, 3, 4, 5]),
                                       ([], []),
                                       ([1, 2, [3, 4], [5]], [1, 2, 3, 4, 5]),
                                       ([1, 2, [3, [], 4], 5, [], 6], [1, 2, 3, 4, 5, 6]),
                                       ([[1, 2], [3, 4, [5, 6], 7, [], [8, [9]]]], [1, 2, 3, 4, 5, 6, 7, 8, 9])])
def test_flatten(lst, truth):
    assert flatten(lst) == truth


def test_generate_upset_series():
    names = ('a', 'b', 'c')
    tuples = [(True, True, True), (True, False, True), (True, True, False), (True, False, False),
              (False, True, True), (False, True, False), (False, False, True)]
    multi_index_truth = pd.MultiIndex.from_tuples(tuples, names=names).sort_values()
    srs_truth = pd.Series(index=multi_index_truth, dtype='uint32').sort_index()
    srs_truth.loc[True, True, True] = 1
    srs_truth.loc[True, True, False] = 2
    srs_truth.loc[True, False, True] = 1
    srs_truth.loc[True, False, False] = 0
    srs_truth.loc[False, True, True] = 1
    srs_truth.loc[False, True, False] = 1
    srs_truth.loc[False, False, True] = 0
    srs = generate_upset_series(
        {'a': {'1', '2', '3', '6'}, 'b': {'2', '3', '4', '5', '6'}, 'c': {'1', '5', '6'}}).sort_index()
    assert srs.index.sort_values().equals(multi_index_truth)
    assert srs.sort_index().equals(srs_truth.sort_index())


@pytest.mark.parametrize('version,expected', [('3.0.0', [3, 0, 0]), ('0.1.3', [0, 1, 3]), ('2.0.5', [2, 0, 5])])
def test_parse_version(version, expected):
    assert parse_version(version) == expected


@pytest.mark.parametrize('attribute_str,expected', [
    ('attr1 "Val1"; attr_2 "val_2"; ', {'attr1': 'Val1', 'attr_2': 'val_2'}),
    ('gene_id "ENSG00000186092"; transcript_id "ENST00000335137"; '
     'gene_name "OR4F5"; gene_source "ensembl_havana"; gene_biotype "protein_coding"; transcript_name "OR4F5-001"; '
     'transcript_source "ensembl_havana"; tag "CCDS"; ccds_id "CCDS30547";',
     {'gene_id': 'ENSG00000186092', 'transcript_id': 'ENST00000335137', 'gene_name': 'OR4F5',
      'gene_source': 'ensembl_havana', 'gene_biotype': 'protein_coding', 'transcript_name': 'OR4F5-001',
      'transcript_source': 'ensembl_havana', 'tag': 'CCDS', 'ccds_id': 'CCDS30547'})
])
def test_parse_gtf_attributes(attribute_str, expected):
    assert parse_gtf_attributes(attribute_str) == expected


@pytest.mark.parametrize('arg,expected', [
    (None, 'NULL'),
    (35, '35'),
    (7.22, '7.22'),
    (-3.52, '-3.52'),
    ('my str', '"my str"'),
    (True, 'TRUE'),
    (False, 'FALSE'),
    (['a', 'ab', 'cde'], 'c("a", "ab", "cde")'),
    ({}, 'c()')
])
def test_parse_r_arg(arg, expected):
    res = _parse_r_arg(arg)
    assert res == expected


@pytest.mark.parametrize('kwargs,expected', [
    ({'arg1': False},'arg1 = FALSE'),
    ({'arg1': 'val1', 'arg2': -0.25, 'arg3': True},'arg1 = "val1", \narg2 = -0.25, \narg3 = TRUE'),
    ({'arg1': None, 'arg2': ['b', 'dc', 'ad f']},'arg1 = NULL, \narg2 = c("b", "dc", "ad f")')
])
def test_python_to_r_kwargs(kwargs, expected):
    res = python_to_r_kwargs(kwargs)
    assert res == expected
