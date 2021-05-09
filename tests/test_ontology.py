import pytest
from rnalysis.utils.ontology import *


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


def test_go_term_with_properties():
    term = GOTerm.with_properties('go_id', 'name', 'namespace', 'level')
    assert isinstance(term, GOTerm)
    assert term.id == 'go_id'
    assert term.name == 'name'
    assert term.namespace == 'namespace'
    assert term.level == 'level'


def test_go_term_set_get_id():
    term = GOTerm()
    assert term.id is None
    term.set_id('id')
    assert term.id == 'id'


def test_go_term_set_get_name():
    term = GOTerm()
    assert term.name is None
    term.set_name('name')
    assert term.name == 'name'


def test_go_term_set_get_namespace():
    term = GOTerm()
    assert term.namespace is None
    term.set_namespace('namespace')
    assert term.namespace == 'namespace'


def test_go_term_set_get_level():
    term = GOTerm()
    assert term.level is None
    term.set_level('level')
    assert term.level == 'level'


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

        last_node = node
        ui_tree = list(dag_tree.upper_induced_graph_iter(last_node))
        ui_tree.sort()
        assert ui_tree == parents_truth[last_node]

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
