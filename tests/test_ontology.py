import os
import matplotlib
import matplotlib.pyplot as plt
import pytest
from rnalysis.utils.ontology import *

matplotlib.use('Agg')


def test_kegg_pathway_parser_construction():
    empty_rel = {'activation': set(), 'binding/association': set(), 'compound': set(), 'dephosphorylation': set(),
                 'dissociation': set(), 'expression': set(), 'glycosylation': set(), 'indirect effect': set(),
                 'inhibition': set(), 'irreversible reaction': set(), 'methylation': set(),
                 'missing interaction': set(),
                 'phosphorylation': set(), 'repression': set(), 'reversible reaction': set(), 'state change': set(),
                 'ubiquitination': set(), 'unknown': set()}
    truth_attrs = {1: dict(id=1, name='hsa:1029', display_name='CDKN2A', type='gene',
                           children_relationships=empty_rel,
                           relationships={'activation': set(), 'binding/association': set(), 'compound': set(),
                                          'dephosphorylation': set(), 'dissociation': set(), 'expression': {(2, 'A')},
                                          'glycosylation': set(), 'indirect effect': set(), 'inhibition': {(3, 'B')},
                                          'irreversible reaction': set(), 'methylation': set(),
                                          'missing interaction': set(), 'phosphorylation': {(3, 'C')},
                                          'repression': set(), 'reversible reaction': set(), 'state change': set(),
                                          'ubiquitination': set(), 'unknown': set()}),
                   2: dict(id=2, name='hsa:1030', display_name='ARF', type='gene', relationships=empty_rel,
                           children_relationships={'activation': set(), 'binding/association': set(),
                                                   'compound': {(4, '->')}, 'dephosphorylation': set(),
                                                   'dissociation': set(), 'expression': {(1, 'A')},
                                                   'glycosylation': set(), 'indirect effect': set(),
                                                   'inhibition': set(), 'irreversible reaction': set(),
                                                   'methylation': set(), 'missing interaction': set(),
                                                   'phosphorylation': set(), 'repression': set(),
                                                   'reversible reaction': set(), 'state change': set(),
                                                   'ubiquitination': set(), 'unknown': set()}),
                   3: dict(id=3, name='hsa:51343', display_name='FZR1', type='gene', relationships=empty_rel,
                           children_relationships={'activation': set(), 'binding/association': set(),
                                                   'compound': set(), 'dephosphorylation': set(),
                                                   'dissociation': set(), 'expression': set(),
                                                   'glycosylation': set(), 'indirect effect': set(),
                                                   'inhibition': {(1, 'B')}, 'irreversible reaction': set(),
                                                   'methylation': set(), 'missing interaction': set(),
                                                   'phosphorylation': {(1, 'C')}, 'repression': set(),
                                                   'reversible reaction': set(), 'state change': set(),
                                                   'ubiquitination': set(), 'unknown': set()}),
                   4: dict(id=4, name='hsa:4171', display_name='MCM2',
                           type='gene', children_relationships=empty_rel,
                           relationships={'activation': set(), 'binding/association': set(), 'compound': {(2, '->')},
                                          'dephosphorylation': set(), 'dissociation': set(), 'expression': set(),
                                          'glycosylation': set(), 'indirect effect': set(), 'inhibition': set(),
                                          'irreversible reaction': set(), 'methylation': set(),
                                          'missing interaction': set(), 'phosphorylation': set(), 'repression': set(),
                                          'reversible reaction': set(), 'state change': set(), 'ubiquitination': set(),
                                          'unknown': set()}),
                   5: dict(id=5, name='undefined', display_name=['3', '4'], type='group', relationships=empty_rel,
                           children_relationships=empty_rel, )}
    with open('tests/test_files/test_kgml.xml') as f:
        tree = ElementTree.parse(f)

    pathway = KEGGPathway(tree, {})
    assert len(pathway.entries) == 5
    assert 'GO:0003840' not in pathway
    for this_id, attrs in truth_attrs.items():
        for attr, val in attrs.items():
            assert getattr(pathway[this_id], attr) == val


def test_plot_kegg_pathway_api():
    with open('tests/test_files/test_kgml.xml') as f:
        tree = ElementTree.parse(f)

    pathway = KEGGPathway(tree, {})
    pathway.plot_pathway()
    plt.close('all')


def test_kegg_entry_with_properties():
    entry = KEGGEntry.with_properties(32, 'name', 'type', 'display_name')
    assert isinstance(entry, KEGGEntry)
    assert entry.id == 32
    assert entry.name == 'name'
    assert entry.type == 'type'
    assert entry.display_name == 'display_name'


def test_kegg_entry_set_get_id():
    entry = KEGGEntry()
    assert entry.id is None
    entry.set_id(12)
    assert entry.id == 12


def test_kegg_entry_set_get_name():
    entry = KEGGEntry()
    assert entry.name is None
    entry.set_name('name')
    assert entry.name == 'name'


def test_kegg_entry_set_get_type():
    entry = KEGGEntry()
    assert entry.type is None
    entry.set_type('type')
    assert entry.type == 'type'


def test_kegg_entry_set_get_display_name():
    entry = KEGGEntry()
    assert entry.display_name is None
    entry.set_display_name('disp')
    assert entry.display_name == 'disp'


def test_parse_go_id():
    line = "is_a: 123456 GO1234567 GO:123456 GO:7654321! primary alcohol metabolic process"
    truth = "GO:7654321"
    assert parse_go_id(line) == truth
    line_2 = "is_a: GO:1234567 GO:0034308 ! primary alcohol metabolic process"
    truth_2 = "GO:1234567"
    assert parse_go_id(line_2) == truth_2


def test_dag_tree_parser_api():
    file = 'tests/test_files/go_mini.obo'
    with open(file, 'r') as f:
        dag_tree = DAGTree(f, ('is_a', 'part_of', 'regulates'))
    with open(file, 'r') as f:
        dag_tree._parse_file(f)


def test_dag_tree_parser_construction():
    file = 'tests/test_files/go_mini.obo'
    levels_truth = [{'GO:0034308': 0, 'GO:0051125': 0},
                    {'GO:0034315': 0, 'GO:0006040': 0, 'GO:0006793': 0, 'GO:0009225': 0, 'GO:0009226': 0},
                    {'GO:2001315': 0, 'GO:2001313': 0}]
    data_version_truth = 'releases/2020-07-16'
    with open(file, 'r') as f:
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
    term = GOTerm.with_properties('go_id', 'name', 'namespace', 42)
    assert isinstance(term, GOTerm)
    assert term.id == 'go_id'
    assert term.name == 'name'
    assert term.namespace == 'namespace'
    assert term.level == 42


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
    term.set_level(42)
    assert term.level == 42


def test_go_term_get_parents():
    file = 'tests/test_files/go_mini.obo'
    with open(file, 'r') as f:
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
    with open(file, 'r') as f:
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
    with open(file, 'r') as f:
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
    with open(file, 'r') as f:
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
    with open(file, 'r') as f:
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
        assert ui_tree == parents_truth_file_2[node]


def test_dag_plot_ontology():
    results = io.load_csv('tests/test_files/go_enrichment_runner_sample_results.csv', index_col=0)
    dag_tree = fetch_go_basic()
    en_score_col = 'colName'
    ontology_graph_format = 'png'
    dag_tree.plot_ontology('biological_process', results, en_score_col, 'title', 'ylabel', ontology_graph_format)
