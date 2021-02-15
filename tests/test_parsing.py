import pytest
from rnalysis.utils.parsing import *


class DummyClass:
    def __init__(self):
        pass


def test_data_to_list():
    assert data_to_list([1, 2, 'hi']) == [1, 2, 'hi']
    assert data_to_list((1, 2, 'hi')) == [1, 2, 'hi']
    assert data_to_list('fifty seven brave men') == ['fifty seven brave men']
    assert sorted(data_to_list({'three', 'different', 'elements'})) == sorted(
        ['three', 'different', 'elements'])
    assert data_to_list(np.array([6, 9, 2])) == [6, 9, 2]
    assert data_to_list(67.2) == [67.2]


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


def test_data_to_tuple():
    assert data_to_tuple([1, 2, 'hi']) == (1, 2, 'hi')
    assert data_to_tuple('fifty seven brave men') == ('fifty seven brave men',)
    assert sorted(data_to_tuple({'three', 'different', 'elements'})) == sorted(
        ('three', 'different', 'elements'))
    assert data_to_tuple(np.array([6, 9, 2])) == (6, 9, 2)
    assert data_to_tuple((67.2,)) == (67.2,)
    assert data_to_tuple(67.2) == (67.2,)


def test_from_string(monkeypatch):
    monkeypatch.setattr('builtins.input', lambda x: 'one\t\ntwo \nthree; and four\n')
    assert from_string() == ['one\t', 'two ', 'three; and four']
    assert from_string(del_spaces=True) == ['one\t', 'two', 'three;andfour']


def test_parse_go_id():
    line = b"is_a: 123456 GO1234567 GO:123456 GO:7654321! primary alcohol metabolic process"
    truth = "GO:7654321"
    assert parse_go_id(line) == truth
    line_2 = b"is_a: GO:1234567 GO:0034308 ! primary alcohol metabolic process"
    truth_2 = "GO:1234567"
    assert parse_go_id(line_2) == truth_2


def test_dag_tree_parser_api():
    file = 'tests/test_files/go_mini.obo'
    with open(file, 'rb') as f:
        dag_tree = DAGTree(f, ('is_a', 'part_of', 'regulates'))
    with open(file, 'rb') as f:
        dag_tree._parse_file(f)


def test_dag_tree_parser_construction():
    file = 'tests/test_files/go_mini.obo'
    levels_truth = [{'GO:0034308': 0, 'GO:0051125': 0},
                    {'GO:0034315': 0, 'GO:0006040': 0, 'GO:0006793': 0, 'GO:0009225': 0, 'GO:0009226': 0},
                    {'GO:2001315': 0, 'GO:2001313': 0}]
    data_version_truth = 'releases/2020-07-16'
    with open(file, 'rb') as f:
        dag_tree = DAGTree(f, ['is_a', 'part_of', 'regulates'])

    assert dag_tree.data_version == data_version_truth
    assert dag_tree.parent_relationship_types == ('is_a', 'part_of', 'regulates')

    assert dag_tree.alt_ids == {'GO:0034619': 'GO:0006040'}
    assert dag_tree['GO:0034619'] == dag_tree['GO:0006040']

    assert len(dag_tree.levels) == len(levels_truth)
    for level, truth_level in zip(dag_tree.levels, levels_truth):
        assert level.keys() == truth_level.keys()
        for item in truth_level:
            assert item in dag_tree

    assert 'GO:0003840' not in dag_tree


def test_go_term_get_parents():
    file = 'tests/test_files/go_mini.obo'
    with open(file, 'rb') as f:
        dag_tree = DAGTree(f, ['is_a', 'part_of', 'regulates'])
    parents_truth = {'GO:0034308': [], 'GO:0051125': [], 'GO:0034315': ['GO:0051125'], 'GO:0006040': ['GO:0034308'],
                     'GO:0006793': ['GO:0034308'], 'GO:0009225': ['GO:0034308'], 'GO:0009226': ['GO:0034308'],
                     'GO:2001315': ['GO:0009226', 'GO:0006793'],
                     'GO:2001313': ['GO:0009225', 'GO:0006793', 'GO:0006040']}
    for node in dag_tree.go_terms:
        assert node in parents_truth
        assert dag_tree[node].get_parents(dag_tree.parent_relationship_types).sort() == parents_truth[node].sort()


def test_go_term_get_children():
    file = 'tests/test_files/go_mini.obo'
    with open(file, 'rb') as f:
        dag_tree = DAGTree(f, ['is_a', 'part_of', 'regulates'])
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
    levels_truth_bioprocess = [{'GO:0034308'},
                               {'GO:0034315', 'GO:0006040', 'GO:0006793', 'GO:0009225', 'GO:0009226'},
                               {'GO:2001315', 'GO:2001313'}]
    levels_truth.reverse()
    levels_truth_bioprocess.reverse()
    file = 'tests/test_files/go_mini.obo'
    with open(file, 'rb') as f:
        dag_tree = DAGTree(f, ['is_a', 'part_of', 'regulates'])

    level_iter = dag_tree.level_iter()
    for truth_level in levels_truth:
        level = set()
        for i in range(len(truth_level)):
            level.add(next(level_iter))
        assert truth_level == level

    level_iter_bioprocess = dag_tree.level_iter('biological_process')
    for truth_level in levels_truth_bioprocess:
        level = set()
        for i in range(len(truth_level)):
            level.add(next(level_iter_bioprocess))
        assert truth_level == level


def test_dag_tree_parser_upper_induced_tree_iterator():
    parents_truth = {'GO:0034308': [], 'GO:0051125': [], 'GO:0034315': ['GO:0051125'], 'GO:0006040': ['GO:0034308'],
                     'GO:0006793': ['GO:0034308'], 'GO:0009225': ['GO:0034308'], 'GO:0009226': ['GO:0034308'],
                     'GO:2001315': ['GO:0009226', 'GO:0006793', 'GO:0034308'],
                     'GO:2001313': ['GO:0009225', 'GO:0006793', 'GO:0006040', 'GO:0034308']}
    file = 'tests/test_files/go_mini.obo'
    with open(file, 'rb') as f:
        dag_tree = DAGTree(f, ['is_a', 'part_of', 'regulates'])

    for node in parents_truth:
        parents_truth[node].sort()
        ui_tree = list(dag_tree.upper_induced_graph_iter(node))
        ui_tree.sort()
        assert ui_tree == parents_truth[node]

    file = 'tests/test_files/obo_for_go_tests.obo'
    with open(file, 'rb') as f:
        dag_tree = DAGTree(f, ['is_a'])
    parents_truth_file_2 = {'GO:0008150': [], 'GO:0050896': ['GO:0008150'], 'GO:0007610': ['GO:0008150'],
                            'GO:0042221': ['GO:0008150', 'GO:0050896'], 'GO:0009605': ['GO:0050896', 'GO:0008150'],
                            'GO:0009991': ['GO:0008150', 'GO:0050896', 'GO:0009605'],
                            'GO:0031667': ['GO:0008150', 'GO:0050896', 'GO:0009605', 'GO:0009991'],
                            'GO:0007584': ['GO:0042221', 'GO:0031667', 'GO:0008150', 'GO:0050896', 'GO:0009605',
                                           'GO:0009991'],
                            'GO:0051780': ['GO:0007584', 'GO:0007610', 'GO:0042221', 'GO:0031667', 'GO:0008150',
                                           'GO:0050896', 'GO:0009605', 'GO:0009991']}
    for node in parents_truth_file_2:
        parents_truth_file_2[node].sort()
        ui_tree = list(dag_tree.upper_induced_graph_iter(node))
        ui_tree.sort()
        print(node)
        assert ui_tree == parents_truth_file_2[node]


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


def test_parse_go_aspects():
    go_dict = {'a': 'aspect a', 'b': 'aspect b', 'c': 'aspect c', '_a_': 'aspect a'}
    assert parse_go_aspects('any', go_dict) == {'aspect a', 'aspect b', 'aspect c'}
    assert parse_go_aspects('B', go_dict) == {'aspect b'}
    assert parse_go_aspects(['b'], go_dict) == {'aspect b'}
    assert parse_go_aspects(['a', 'z', 'c'], go_dict) == {'aspect a', 'aspect c', 'z'}
    assert parse_go_aspects(['b', 'c', 'A'], go_dict) == {'aspect a', 'aspect b', 'aspect c'}
    assert parse_go_aspects(['aspect z'], go_dict) == {'aspect z'}


def test_parse_evidence_types():
    ev_dict = {'a': 'eva', 'b': 'evb', 'c': 'evc', 'ab': {'eva', 'evb'}, 'bc': {'evb', 'evc'}, 'de': {'evd', 'eve'}}
    assert parse_evidence_types('any', ev_dict) == {'eva', 'evb', 'evc', 'evd', 'eve'}
    assert parse_evidence_types('bc', ev_dict) == {'evb', 'evc'}
    assert parse_evidence_types('c', ev_dict) == {'evc'}
    assert parse_evidence_types({'a', 'bc', 'f'}, ev_dict) == {'eva', 'evb', 'evc', 'f'}
    assert parse_evidence_types({'a', 'ab'}, ev_dict) == {'eva', 'evb'}
    assert parse_evidence_types({'z', 'v'}, ev_dict) == {'z', 'v'}


def test_sparse_dict_to_bool_df():
    truth = pd.read_csv('tests/test_files/sparse_dict_to_df_truth.csv', index_col=0).sort_index(axis=1)
    sparse_dict = {'gene1': {'a', 'c'}, 'gene2': {'b'}, 'gene3': {'b', 'g', 'h'}, 'gene4': {'e', 'd', 'a', 'c', 'f'},
                   'gene5': set()}
    res = sparse_dict_to_bool_df(sparse_dict).sort_index(axis=1)
    assert res.equals(truth)
