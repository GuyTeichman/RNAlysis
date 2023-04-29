import pytest

from rnalysis.utils.parsing import *
from rnalysis.utils.parsing import _parse_r_arg


class DummyClass:
    def __init__(self):
        pass

    def __eq__(self, other):
        return type(self) == type(other)


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
                          (DummyClass(), False, [DummyClass()]),
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
                          (DummyClass(), False, (DummyClass(),)),
                          (None, False, (None,))])
def test_data_to_tuple(input_val, sort, truth):
    assert data_to_tuple(input_val, sort=sort) == truth


@pytest.mark.parametrize('val,expected', [
    ([1, 2, 'hi'], {1, 2, 'hi'}),
    ((1, 2, 'hi'), {1, 2, 'hi'}),
    ('fifty seven brave men', {'fifty seven brave men'}),
    ({'three', 'different', 'elements'}, {'three', 'different', 'elements'}),
    (np.array([6, 9, 2]), {6, 9, 2}),
    (67.2, {67.2}),
    (None, {None}),
])
def test_data_to_set(val, expected):
    assert data_to_set(val) == expected


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
    ({'arg1': False}, 'arg1 = FALSE'),
    ({'arg1': 'val1', 'arg2': -0.25, 'arg3': True}, 'arg1 = "val1",\narg2 = -0.25,\narg3 = TRUE'),
    ({'arg1': None, 'arg2': ['b', 'dc', 'ad f']}, 'arg1 = NULL,\narg2 = c("b", "dc", "ad f")')
])
def test_python_to_r_kwargs(kwargs, expected):
    res = python_to_r_kwargs(kwargs)
    assert res == expected


@pytest.mark.parametrize('attrs,expected', [
    ('ID=mrna0001;Name=sonichedgehog', {'ID': 'mrna0001', 'Name': 'sonichedgehog'}),
    ('ID=exon00003;Parent=mrna0001', {'ID': 'exon00003', 'Parent': 'mrna0001'}),
    ('ID=B0348.6a.1;Parent=WBGene00002061', {'ID': 'B0348.6a.1', 'Parent': 'WBGene00002061'}),
    ('ID=exon00004;Parent=mRNA00001,mRNA00002,mRNA00003', {'ID': 'exon00004', 'Parent':
        ['mRNA00001', 'mRNA00002', 'mRNA00003']}),
    ('ID=exon00005;Parent=mRNA00001,mRNA00003', {'ID': 'exon00005', 'Parent': ['mRNA00001', 'mRNA00003']}),
    ('ID=mRNA00003;Parent=gene00001;Name=EDEN.3', {'ID': 'mRNA00003', 'Parent': 'gene00001', 'Name': 'EDEN.3'}),
])
def test_parse_gff3_attributes(attrs, expected):
    res = parse_gff3_attributes(attrs)
    assert res == expected


@pytest.mark.parametrize('allow_unicode,string,expected', [
    (False, 'text', 'text'),
    (False, 'text with spaces', 'text-with-spaces'),
    (True, 'text_without_spaces', 'text_without_spaces'),
    (False, 'text---without--spaces', 'text-without-spaces'),
    (False, 'dash_', 'dash'),
    (True, 'dash-', 'dash'),
    (False, 'more!characters?:&%a', 'morecharactersa')
])
def test_slugify(allow_unicode, string, expected):
    res = slugify(string, allow_unicode)
    assert res == expected


@pytest.mark.parametrize('regex,repl,item,expected', [
    ('cat', 'dog', 'The cat sat on the cat', 'The cat sat on the dog'),
    ('1', 'one', '1234 1 5678 1', '1234 1 5678 one'),
    ('([A-Za-z]+)', r'\1\1', 'hello world', 'hello worldworld'),
    ('abc', 'xyz', 'abcdefabc', 'abcdefxyz'),
    ('[0-9]+', '', 'abcd1234', 'abcd'),
    ('(abc)+', 'def', 'abcabcabc', 'def'),
    ('', 'xyz', 'Hello, world!', 'Hello, world!xyz'),
    ('[A-Z][a-z]+', 'Mr.', 'Hello John Smith', 'Hello John Mr.'),
    ('[0-9]{2,3}', 'one', '12345', '123one'),
    ('e', 'replaced', 'escaped', 'escapreplacedd')
])
def test_replace_last_occurrence(regex, repl, item, expected):
    res = replace_last_occurrence(regex, repl, item)
    assert res == expected


@pytest.fixture
def sample_df():
    return pd.DataFrame({
        "A": [1.0, 2, 3, 4, 5],
        "B": [6.00111, 7, 8, 9, 10],
        "C": [11, 12, 13, 14, 15],
        "D": [16, 17, 18, 19, 20],
        "E": [21, 22, 23, 24, 25]
    })


def test_df_to_html_returns_string(sample_df):
    html = df_to_html(sample_df)
    assert isinstance(html, str)


@pytest.mark.parametrize('n_rows', [2, 3, 4])
def test_df_to_html_limits_number_of_rows(sample_df, n_rows):
    html = df_to_html(sample_df, max_rows=n_rows)
    assert html.count("<tr>") == n_rows + 2  # includes header row and '...'


def test_df_to_html_formats_floats_with_two_decimals(sample_df):
    html = df_to_html(sample_df)
    assert "1.00" in html
    assert "6.00" in html


def test_df_to_html_does_not_truncate_values(sample_df):
    long_string = "a" * 100
    df = pd.DataFrame({
        "A": [1, 2, 3],
        "B": [long_string, long_string, long_string],
    })
    html = df_to_html(df)
    assert long_string in html


def test_df_to_html_removes_redundant_dots(sample_df):
    df = pd.DataFrame({
        "A": ["A", "B", "C"],
        "B": ["...", "...", "..."],
    })
    html = df_to_html(df, max_rows=2, max_cols=2)
    assert "...<br>" not in html


@pytest.mark.parametrize("items, expected_output", [
    ([], '<table border="1" class="dataframe">\n</table>'),
    (['item'],
     '<table border="1" class="dataframe">\n<tr><td style="border: 1px solid black; border-collapse: collapse;"><b>item</b></td></tr>\n</table>'),
    (['item1', 'item2', 'item3'],
     '<table border="1" class="dataframe">\n<tr><td style="border: 1px solid black; border-collapse: collapse;"><b>item1</b></td></tr>\n<tr><td style="border: 1px solid black; border-collapse: collapse;"><b>item2</b></td></tr>\n<tr><td style="border: 1px solid black; border-collapse: collapse;"><b>item3</b></td></tr>\n</table>'),
])
def test_items_to_html_table(items, expected_output):
    assert items_to_html_table(items) == expected_output


@pytest.mark.parametrize("d, expected", [
    ({}, ''),
    ({'a': 1}, "a = 1"),
    ({'a': 1, 'b': 'hello'}, "a = 1, \nb = 'hello'"),
    ({'a': [1, 2, 3], 'b': {'x': 1}}, "a = [1, 2, 3], \nb = {'x': 1}")
])
def test_format_dict_for_display(d, expected):
    assert format_dict_for_display(d) == expected
