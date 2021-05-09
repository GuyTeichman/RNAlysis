import pytest

from rnalysis.utils.parsing import *


class DummyClass:
    def __init__(self):
        pass


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
