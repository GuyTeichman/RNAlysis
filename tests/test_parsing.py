import inspect
import re

import pytest

from rnalysis.exceptions import InternalError
from rnalysis.utils.parsing import *
from rnalysis.utils.parsing import _parse_r_arg


class DummyClass:
    def __init__(self):
        pass

    def __eq__(self, other):
        return type(self) == type(other)


@pytest.mark.parametrize('is_input_path', [True, False])
@pytest.mark.parametrize(
    'path,expected_path',
    [
        ('file.csv', 'file'),
        ('path/to/file.exe', 'path/to/file'),
        ('name.fastq.gz', 'name'),
        ('path/to/a/filename.sam.tar.gz', 'path/to/a/filename'),
    ],
)
def test_remove_suffixes(path, expected_path, is_input_path):
    if is_input_path:
        assert remove_suffixes(Path(path)) == Path(expected_path)
    else:
        assert remove_suffixes(path) == expected_path


@pytest.mark.parametrize(
    'input_val,sort,truth',
    [
        ([1, 2, 'hi'], False, [1, 2, 'hi']),
        ((1, 2, 'hi'), False, [1, 2, 'hi']),
        ('fifty seven brave men', False, ['fifty seven brave men']),
        ('fifty seven brave men', True, ['fifty seven brave men']),
        ({'three', 'different', 'elements'}, True, ['different', 'elements', 'three']),
        (np.array([6, 9, 2]), False, [6, 9, 2]),
        (np.array([6, 9, 2]), True, [2, 6, 9]),
        (DummyClass(), False, [DummyClass()]),
        (None, False, [None]),
    ],
)
def test_data_to_list(input_val, sort, truth):
    assert data_to_list(input_val, sort=sort) == truth


@pytest.mark.parametrize(
    'input_val,sort,truth',
    [
        ([1, 2, 'hi'], False, (1, 2, 'hi')),
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
        (None, False, (None,)),
    ],
)
def test_data_to_tuple(input_val, sort, truth):
    assert data_to_tuple(input_val, sort=sort) == truth


@pytest.mark.parametrize(
    'val,expected',
    [
        ([1, 2, 'hi'], {1, 2, 'hi'}),
        ((1, 2, 'hi'), {1, 2, 'hi'}),
        ('fifty seven brave men', {'fifty seven brave men'}),
        ({'three', 'different', 'elements'}, {'three', 'different', 'elements'}),
        (np.array([6, 9, 2]), {6, 9, 2}),
        (67.2, {67.2}),
        (None, {None}),
    ],
)
def test_data_to_set(val, expected):
    assert data_to_set(val) == expected


def test_from_string(monkeypatch):
    monkeypatch.setattr('builtins.input', lambda x: 'one\t\ntwo \nthree; and four\n')
    assert from_string() == ['one\t', 'two ', 'three; and four']
    assert from_string(del_spaces=True) == ['one\t', 'two', 'three;andfour']


def test_uniprot_tab_to_dict():
    tab = (
        'From\tTo\nWBGene00019883\tP34544\nWBGene00023497\tQ27395\nWBGene00003515\tP12844\nWBGene00000004\t'
        'A0A0K3AVL7\nWBGene00000004\tO17395\n'
    )
    tab_rev = (
        'From\tTo\nP34544\tWBGene00019883\nQ27395\tWBGene00023497\nP12844\tWBGene00003515\nA0A0K3AVL7\t'
        'WBGene00000004\nO17395\tWBGene00000004\n'
    )

    truth = (
        {'WBGene00019883': 'P34544', 'WBGene00023497': 'Q27395', 'WBGene00003515': 'P12844'},
        ['A0A0K3AVL7', 'O17395'],
    )
    truth_rev = (
        {
            'P34544': 'WBGene00019883',
            'Q27395': 'WBGene00023497',
            'P12844': 'WBGene00003515',
            'A0A0K3AVL7': 'WBGene00000004',
            'O17395': 'WBGene00000004',
        },
        [],
    )
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
    tab = (
        'Entry\tAnnotation\tyourlist:M20200816216DA2B77BFBD2E6699CA9B6D1C41EB2A5FE6AF\nP34544\t5 out of 5\t'
        'WBGene00019883\nQ27395\t4 out of 5\tWBGene00023497\nP12844\t5 out of 5\tWBGene00003515\nA0A0K3AVL7\t'
        '1 out of 5\tWBGene00000004\nO17395\t2 out of 5\tWBGene00000004\n'
    )

    truth = {
        'WBGene00019883': 'P34544',
        'WBGene00023497': 'Q27395',
        'WBGene00003515': 'P12844',
        'WBGene00000004': 'O17395',
    }
    truth_rev = {
        'P34544': 'WBGene00019883',
        'Q27395': 'WBGene00023497',
        'P12844': 'WBGene00003515',
        'A0A0K3AVL7': 'WBGene00000004',
        'O17395': 'WBGene00000004',
    }
    assert truth == uniprot_tab_with_score_to_dict(tab)

    assert truth_rev == uniprot_tab_with_score_to_dict(tab, True)


@pytest.mark.parametrize(
    'lst,chunk_size,truth',
    [
        ([1, 2, 3, 4, 5], 1, [[1], [2], [3], [4], [5]]),
        ([1, 2, 3, 4], 2, [[1, 2], [3, 4]]),
        ([1, 2, 3, 4, 5], 2, [[1, 2], [3, 4], [5]]),
        ([], 3, [[]]),
        ((1, 2, 3, 4), 3, [(1, 2, 3), (4,)]),
        ((1, 2, 3, 4, 5), 5, [(1, 2, 3, 4, 5)]),
        ((1, 2, 3, 4), 5, [(1, 2, 3, 4)]),
        (tuple(), 1, [tuple()]),
    ],
)
def test_partition_list(lst, chunk_size, truth):
    res = partition_list(lst, chunk_size)
    assert res == truth


@pytest.mark.parametrize(
    'lst,chunk_size', [([1, 2, 3, 4], 0), ([1, 2], -1), ({1, 2}, 2), ([1, 2, 3], 1.0), ((i for i in range(3)), 2)]
)
def test_partition_list_invalid_input(lst, chunk_size):
    with pytest.raises(InternalError):
        _ = partition_list(lst, chunk_size)


@pytest.mark.parametrize(
    'lst,truth',
    [
        ([1, 2, 3, 4, 5], [1, 2, 3, 4, 5]),
        ([], []),
        ([1, 2, [3, 4], [5]], [1, 2, 3, 4, 5]),
        ([1, 2, [3, [], 4], 5, [], 6], [1, 2, 3, 4, 5, 6]),
        ([[1, 2], [3, 4, [5, 6], 7, [], [8, [9]]]], [1, 2, 3, 4, 5, 6, 7, 8, 9]),
    ],
)
def test_flatten(lst, truth):
    assert flatten(lst) == truth


def test_generate_upset_series():
    names = ('a', 'b', 'c')
    tuples = [
        (True, True, True),
        (True, False, True),
        (True, True, False),
        (True, False, False),
        (False, True, True),
        (False, True, False),
        (False, False, True),
    ]
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
        {'a': {'1', '2', '3', '6'}, 'b': {'2', '3', '4', '5', '6'}, 'c': {'1', '5', '6'}}
    ).sort_index()
    assert srs.index.sort_values().equals(multi_index_truth)
    assert srs.sort_index().equals(srs_truth.sort_index())


@pytest.mark.parametrize('version,expected', [('3.0.0', [3, 0, 0]), ('0.1.3', [0, 1, 3]), ('2.0.5', [2, 0, 5])])
def test_parse_version(version, expected):
    assert parse_version(version) == expected


@pytest.mark.parametrize(
    'attribute_str,expected',
    [
        ('attr1 "Val1"; attr_2 "val_2"; ', {'attr1': 'Val1', 'attr_2': 'val_2'}),
        (
            'gene_id "ENSG00000186092"; transcript_id "ENST00000335137"; '
            'gene_name "OR4F5"; gene_source "ensembl_havana"; gene_biotype "protein_coding"; transcript_name "OR4F5-001"; '
            'transcript_source "ensembl_havana"; tag "CCDS"; ccds_id "CCDS30547";',
            {
                'gene_id': 'ENSG00000186092',
                'transcript_id': 'ENST00000335137',
                'gene_name': 'OR4F5',
                'gene_source': 'ensembl_havana',
                'gene_biotype': 'protein_coding',
                'transcript_name': 'OR4F5-001',
                'transcript_source': 'ensembl_havana',
                'tag': 'CCDS',
                'ccds_id': 'CCDS30547',
            },
        ),
    ],
)
def test_parse_gtf_attributes(attribute_str, expected):
    assert parse_gtf_attributes(attribute_str) == expected


@pytest.mark.parametrize(
    'arg,expected',
    [
        (None, 'NULL'),
        (35, '35'),
        (7.22, '7.22'),
        (-3.52, '-3.52'),
        ('my str', '"my str"'),
        (True, 'TRUE'),
        (False, 'FALSE'),
        (['a', 'ab', 'cde'], 'c("a", "ab", "cde")'),
        ({}, 'c()'),
    ],
)
def test_parse_r_arg(arg, expected):
    res = _parse_r_arg(arg)
    assert res == expected


@pytest.mark.parametrize(
    'kwargs,expected',
    [
        ({'arg1': False}, 'arg1 = FALSE'),
        ({'arg1': 'val1', 'arg2': -0.25, 'arg3': True}, 'arg1 = "val1",\narg2 = -0.25,\narg3 = TRUE'),
        ({'arg1': None, 'arg2': ['b', 'dc', 'ad f']}, 'arg1 = NULL,\narg2 = c("b", "dc", "ad f")'),
    ],
)
def test_python_to_r_kwargs(kwargs, expected):
    res = python_to_r_kwargs(kwargs)
    assert res == expected


@pytest.mark.parametrize(
    'attrs,expected',
    [
        ('ID=mrna0001;Name=sonichedgehog', {'ID': 'mrna0001', 'Name': 'sonichedgehog'}),
        ('ID=exon00003;Parent=mrna0001', {'ID': 'exon00003', 'Parent': 'mrna0001'}),
        ('ID=B0348.6a.1;Parent=WBGene00002061', {'ID': 'B0348.6a.1', 'Parent': 'WBGene00002061'}),
        (
            'ID=exon00004;Parent=mRNA00001,mRNA00002,mRNA00003',
            {'ID': 'exon00004', 'Parent': ['mRNA00001', 'mRNA00002', 'mRNA00003']},
        ),
        ('ID=exon00005;Parent=mRNA00001,mRNA00003', {'ID': 'exon00005', 'Parent': ['mRNA00001', 'mRNA00003']}),
        ('ID=mRNA00003;Parent=gene00001;Name=EDEN.3', {'ID': 'mRNA00003', 'Parent': 'gene00001', 'Name': 'EDEN.3'}),
    ],
)
def test_parse_gff3_attributes(attrs, expected):
    res = parse_gff3_attributes(attrs)
    assert res == expected


@pytest.mark.parametrize(
    'input_string,expected',
    [
        ('P34544', 'P34544'),  # standard 6-character accession, letter outside the O/P/Q range
        ('Q27395', 'Q27395'),  # standard 6-character accession with an O/P/Q prefix
        ('O17395', 'O17395'),  # standard 6-character accession with an O/P/Q prefix
        ('A0A0K3AVL7', 'A0A0K3AVL7'),  # extended 10-character accession format
        ('sp|P34544|GENE_HUMAN', 'P34544'),  # accession embedded in a UniProt FASTA header
        ('UniProtKB:Q27395', 'Q27395'),  # accession embedded after a database prefix
        ('P34544-2', 'P34544'),  # isoform suffix is not part of the match and is dropped
        ('P34544.1', 'P34544'),  # version suffix is not part of the match and is dropped
        ('P34544-2.3', 'P34544'),  # combined isoform + version suffix is dropped
        ('not a uniprot id', None),  # no match anywhere in the string
        ('', None),  # empty string
        ('p34544', None),  # malformed: lowercase is not a valid accession (case-sensitive)
        ('P3454', None),  # malformed: one character too short to be a valid accession
        ('unicode text café P34544 café', 'P34544'),  # valid id surrounded by unicode text
    ],
)
def test_parse_uniprot_id(input_string, expected):
    assert parse_uniprot_id(input_string) == expected


@pytest.mark.parametrize(
    'allow_unicode,string,expected',
    [
        (False, 'text', 'text'),
        (False, 'text with spaces', 'text-with-spaces'),
        (True, 'text_without_spaces', 'text_without_spaces'),
        (False, 'text---without--spaces', 'text-without-spaces'),
        (False, 'dash_', 'dash'),
        (True, 'dash-', 'dash'),
        (False, 'more!characters?:&%a', 'morecharactersa'),
    ],
)
def test_slugify(allow_unicode, string, expected):
    res = slugify(string, allow_unicode)
    assert res == expected


@pytest.mark.parametrize(
    's1,s2,expected',
    [
        ('hello_world', 'goodbye_world', '_world'),  # shared suffix is the longest common substring
        ('abcdef', 'xabcdefy', 'abcdef'),  # one string fully contained within the other
        ('abc', 'xyz', ''),  # no common substring at all
        ('same', 'same', 'same'),  # identical strings: LCS is the whole string
        ('', '', ''),  # both strings empty
        ('', 'abc', ''),  # one string empty
        ('cafe123', 'cafe456', 'cafe'),  # shared prefix
        ('café123', 'café456', 'café'),  # unicode shared prefix
    ],
)
def test_longest_common_substring(s1, s2, expected):
    assert longest_common_substring(s1, s2) == expected


@pytest.mark.parametrize(
    'strings,expected',
    [
        ([], ''),  # empty input
        (['hello'], 'hello'),  # single item: the whole string is returned as-is
        (['file_A.txt', 'file_B.txt'], '.txt'),  # a real shared suffix
        (['abc', 'abc', 'abc'], 'abc'),  # all-identical inputs: the whole string is the suffix
        (['apple', 'banana'], ''),  # no common suffix
        (['', ''], ''),  # all strings empty
        (['a', 'ba', 'cba'], 'a'),  # strings of different lengths, suffix bounded by the shortest one
        (['日本語abc', '中国語abc'], '語abc'),  # unicode suffix
        (['ab', 'aab'], 'ab'),  # lexicographically-smaller string is longer than the larger one (issue #180)
        (['aac', 'ab'], ''),  # same as above, but with no shared suffix at all
        # issue #195: a string that is neither the lexicographic min nor max must still be able to break
        # the suffix. 'ax' and 'zx' (min/max) both end in 'x', but 'cy' does not -> no common suffix.
        (['zx', 'ax', 'cy'], ''),
        (['abc', 'xbc', 'ybc'], 'bc'),  # issue #195: a genuine shared suffix across 3+ strings
    ],
)
def test_common_suffix(strings, expected):
    assert common_suffix(strings) == expected


@pytest.mark.parametrize(
    's,suffix,expected',
    [
        ('sample_trimmed', '_trimmed', 'sample'),  # single suffix given as a plain string
        ('sample_trimmed', ['_sorted', '_trimmed'], 'sample'),  # suffix found within a list of candidates
        ('sample', '_trimmed', 'sample'),  # no matching suffix: string is returned unchanged
        ('sample_trimmed', ['med', '_trimmed'], 'sample'),  # longest matching suffix wins over a shorter one
        ('sample_med', ['_trimmed', 'med'], 'sample_'),  # falls back to a shorter suffix when the longest doesn't match
        ('sample', [], 'sample'),  # empty suffix list: string is returned unchanged
        ('café_trimmed', '_trimmed', 'café'),  # unicode string
        ('', '_trimmed', ''),  # empty string
        ('_trimmed', '_trimmed', ''),  # the whole string is the suffix
    ],
)
def test_remove_suffix(s, suffix, expected):
    assert remove_suffix(s, suffix) == expected


@pytest.mark.parametrize(
    'file_pairs,expected',
    [
        ([], []),  # empty input
        ([('sample1_R1', 'sample1_R2')], ['sample1_R']),  # single pair
        ([('sample1_R1', 'sample1_R2'), ('sample2_R1', 'sample2_R2')], ['sample1', 'sample2']),  # two pairs
        (
            [('ctrl1_R1', 'ctrl1_R2'), ('ctrl2_R1', 'ctrl2_R2'), ('treat1_R1', 'treat1_R2')],
            ['ctrl1', 'ctrl2', 'treat1'],
        ),  # issue #195: the genuinely shared "_R" suffix is now trimmed evenly from every pair. Before
        # the common_suffix fix this ignored the non-min/max 'ctrl2_R' and mis-trimmed "1_R", yielding
        # the inconsistent ['ctrl', 'ctrl2_R', 'treat'].
        ([('abc', 'xyz')], ['abcxyz']),  # no common substring: names are concatenated
        ([('same', 'same')], ['same']),  # all-identical inputs
        ([('sample1_R1_trimmed', 'sample1_R2_trimmed')], ['sample1_R']),  # known suffixes are stripped first
        (
            [('x_R1', 'x_R2'), ('x_R1', 'x_R2')],
            ['x_R1x_R2', 'x_R1x_R2'],
        ),  # duplicate pairs: trimming to empty is avoided
        ([('café_R1', 'café_R2')], ['café_R']),  # unicode
        # regression test for issue #181: a scope leak used the LAST pair's s1/s2 (from the loop
        # variable) to filter every pair's result instead of each pair's own s1/s2. Since the last
        # pair here is much shorter, the earlier (longer) pair's valid LCS result was wrongly
        # excluded from the common-suffix computation, so the shared trailing "_R" was never trimmed.
        ([('sampleA_long_name_R1', 'sampleA_long_name_R2'), ('sB_R1', 'sB_R2')], ['sampleA_long_name', 'sB']),
        # regression test for issue #195: the third pair shares nothing beyond the global "common" LCS, so
        # it falls back to concatenation ('commonPPcommonQQ'). The old tautological `lcs_only` filter kept
        # that fallback in the common-suffix computation, so its trailing 'QQ' masked the real "_S" suffix
        # shared by the genuine LCS results and it was never trimmed. Now the suffix is computed over the
        # genuine LCS results only, so "_S" is correctly stripped from the first two names.
        (
            [('commonAB_S1', 'commonAB_S2'), ('commonCD_S1', 'commonCD_S2'), ('commonPP', 'commonQQ')],
            ['commonAB', 'commonCD', 'commonPPcommonQQ'],
        ),
    ],
)
def test_generate_common_name(file_pairs, expected):
    assert generate_common_name(file_pairs) == expected


@pytest.mark.parametrize(
    'regex,repl,item,expected',
    [
        ('cat', 'dog', 'The cat sat on the cat', 'The cat sat on the dog'),
        ('1', 'one', '1234 1 5678 1', '1234 1 5678 one'),
        ('([A-Za-z]+)', r'\1\1', 'hello world', 'hello worldworld'),
        ('abc', 'xyz', 'abcdefabc', 'abcdefxyz'),
        ('[0-9]+', '', 'abcd1234', 'abcd'),
        ('(abc)+', 'def', 'abcabcabc', 'def'),
        ('', 'xyz', 'Hello, world!', 'Hello, world!xyz'),
        ('[A-Z][a-z]+', 'Mr.', 'Hello John Smith', 'Hello John Mr.'),
        ('[0-9]{2,3}', 'one', '12345', '123one'),
        ('e', 'replaced', 'escaped', 'escapreplacedd'),
    ],
)
def test_replace_last_occurrence(regex, repl, item, expected):
    res = replace_last_occurrence(regex, repl, item)
    assert res == expected


@pytest.fixture
def sample_df():
    return pl.DataFrame(
        {
            'A': [1.0, 2, 3, 4, 5],
            'B': [6.00111, 7, 8, 9, 10],
            'C': [11, 12, 13, 14, 15],
            'D': [16, 17, 18, 19, 20],
            'E': [21, 22, 23, 24, 25],
        }
    )


def test_df_to_html_returns_string(sample_df):
    html = df_to_html(sample_df)
    assert isinstance(html, str)


@pytest.mark.parametrize('n_rows', [2, 3, 4])
def test_df_to_html_limits_number_of_rows(sample_df, n_rows):
    html = df_to_html(sample_df, max_rows=n_rows)
    assert html.count('<tr>') == n_rows + 1  # includes header row and '...'


def test_df_to_html_formats_floats_with_two_decimals(sample_df):
    html = df_to_html(sample_df)
    assert '1.00' in html
    assert '6.00' in html


def test_df_to_html_does_not_truncate_values(sample_df):
    long_string = 'a' * 100
    df = pl.DataFrame(
        {
            'A': [1, 2, 3],
            'B': [long_string, long_string, long_string],
        }
    )
    html = df_to_html(df)
    assert long_string in html


def test_df_to_html_removes_redundant_dots(sample_df):
    df = pl.DataFrame(
        {
            'A': ['A', 'B', 'C'],
            'B': ['...', '...', '...'],
        }
    )
    html = df_to_html(df, max_rows=2, max_cols=2)
    assert '...<br>' not in html


@pytest.mark.parametrize(
    'items, expected_output',
    [
        ([], '<table border="1" class="dataframe">\n</table>'),
        (
            ['item'],
            '<table border="1" class="dataframe">\n<tr><td style="border: 1px solid black; border-collapse: collapse;"><b>item</b></td></tr>\n</table>',
        ),
        (
            ['item1', 'item2', 'item3'],
            '<table border="1" class="dataframe">\n<tr><td style="border: 1px solid black; border-collapse: collapse;"><b>item1</b></td></tr>\n<tr><td style="border: 1px solid black; border-collapse: collapse;"><b>item2</b></td></tr>\n<tr><td style="border: 1px solid black; border-collapse: collapse;"><b>item3</b></td></tr>\n</table>',
        ),
    ],
)
def test_items_to_html_table(items, expected_output):
    assert items_to_html_table(items) == expected_output


@pytest.mark.parametrize(
    'd, expected',
    [
        ({}, ''),
        ({'a': 1}, 'a = 1'),
        ({'a': 1, 'b': 'hello'}, "a = 1, \nb = 'hello'"),
        ({'a': [1, 2, 3], 'b': {'x': 1}}, "a = [1, 2, 3], \nb = {'x': 1}"),
    ],
)
def test_format_dict_for_display(d, expected):
    assert format_dict_for_display(d) == expected


class TestRMakeNames:
    def test_valid_names(self):
        names = ['validname', 'another.validname', 'yet.another.valid.name']
        expected_output = ['validname', 'another.validname', 'yet.another.valid.name']
        assert r_make_names(names) == expected_output

    def test_names_with_invalid_characters(self):
        names = ['invalid!name', 'another@invalid#name', 'yet$another%invalid^name']
        expected_output = ['invalid.name', 'another.invalid.name', 'yet.another.invalid.name']
        assert r_make_names(names) == expected_output

    def test_names_with_underscore(self):
        names = ['valid_name', 'another_valid_name_', 'yet_another__valid_name', '_invalid']
        expected_output = ['valid_name', 'another_valid_name_', 'yet_another__valid_name', '._invalid']
        assert r_make_names(names) == expected_output

    def test_names_starting_with_invalid_characters(self):
        names = ['!invalidname', '@anotherinvalidname', '$yetanotherinvalidname']
        expected_output = ['.invalidname', '.anotherinvalidname', '.yetanotherinvalidname']
        assert r_make_names(names) == expected_output

    def test_names_starting_with_dot_followed_by_number(self):
        names = ['.1invalidname', '.2anotherinvalidname', '.3yetanotherinvalidname']
        expected_output = ['..1invalidname', '..2anotherinvalidname', '..3yetanotherinvalidname']
        assert r_make_names(names) == expected_output

    def test_duplicate_names(self):
        names = ['duplicate', 'duplicate', 'anotherduplicate', 'anotherduplicate', 'unique']
        expected_output = ['duplicate', 'duplicate.1', 'anotherduplicate', 'anotherduplicate.1', 'unique']
        assert r_make_names(names) == expected_output

    def test_duplicate_names_with_invalid_characters(self):
        names = ['duplicate!', 'duplicate@', 'another#duplicate', 'another$duplicate', 'unique%']
        expected_output = ['duplicate.', 'duplicate..1', 'another.duplicate', 'another.duplicate.1', 'unique.']
        assert r_make_names(names) == expected_output

    def test_empty_list(self):
        names = []
        expected_output = []
        assert r_make_names(names) == expected_output


class TestSanitizeFilenameComponent:
    """issue #281: user values are spliced into output file names, so they must be made filename-safe."""

    @pytest.mark.parametrize(
        'text',
        [
            'cond1',
            'rep_2',
            'KO-1.5',
            'wild type',
            "KO'1",
            'poly(age, degree = 2)',
            'a+b=c,d[e]{f}@g#h%i&j!k~l$m;n^o`p',
            'Wildtyp-\u00fc',
            '\u05d0\u05d1',
            # non-ASCII symbols and punctuation are legal file name characters on every OS, so they are kept
            '37\u00b0C',
            'WT\u2013KO',
            'KO\u2019s',
            '\u00b1 2\u00d7',
            'non\u00a0breaking',
            '',
        ],
    )
    def test_filename_safe_values_unchanged(self, text):
        # values that are already safe must keep producing byte-identical file names (no arbitrary changes)
        assert sanitize_filename_component(text) == text

    @pytest.mark.parametrize(
        'text,expected',
        [
            ('KO/1', 'KO_1'),
            ('KO:1', 'KO_1'),
            ('back\\slash', 'back_slash'),
            ('q"uote', 'q_uote'),
            ('a<b>c|d?e*f', 'a_b_c_d_e_f'),
            ('tab\there\nnewline', 'tab_here_newline'),
            ('line\u2028separator\x7fdel', 'line_separator_del'),
        ],
    )
    def test_unsafe_characters_replaced(self, text, expected):
        assert sanitize_filename_component(text) == expected

    @pytest.mark.parametrize(
        'text,replacements,expected',
        [
            # the established LRT interaction-term mappings must not change
            ('factor1:factor2', {':': '_x_'}, 'factor1_x_factor2'),
            ('factor1:factor2', {':': 'X'}, 'factor1Xfactor2'),
            ('a:b/c', {':': '_x_'}, 'a_x_b_c'),
        ],
    )
    def test_replacements(self, text, replacements, expected):
        assert sanitize_filename_component(text, replacements) == expected

    def test_unsafe_replacement_is_rejected(self):
        with pytest.raises(ValueError):
            sanitize_filename_component('a:b', {':': '/'})


class TestRQuote:
    """issue #281: Python strings entering R source must be escaped deliberately, not by coincidence."""

    @pytest.mark.parametrize(
        'text,expected',
        [
            ('cond1', '"cond1"'),
            ('', '""'),
            ("KO'1", '"KO\'1"'),
            ('q"uote', '"q\\"uote"'),
            ('back\\slash', '"back\\\\slash"'),
            ('C:/path/to/file.csv', '"C:/path/to/file.csv"'),
            ('\\"', '"\\\\\\""'),
            ('Wildtyp-\u00fc', '"Wildtyp-\u00fc"'),
        ],
    )
    def test_r_quote(self, text, expected):
        assert r_quote(text) == expected

    def test_r_quote_non_string(self):
        assert r_quote(5) == '"5"'


@pytest.mark.parametrize(
    'sample_grouping, mode, expected_output',
    [
        (['group1', 'group2'], 'display', ['group1', 'group2']),
        (['group1', 'group2'], 'auto', ['group1', 'group2']),
        ([['group1', 'subgroup1'], 'group2'], 'display', ['group1\nsubgroup1', 'group2']),
        ([['group1', 'subgroup1'], 'group2'], 'auto', ['group1;subgroup1', 'group2']),
        ([('group1', 'subgroup1'), 'group2'], 'display', ['group1\nsubgroup1', 'group2']),
        ([('group1', 'subgroup1'), 'group2'], 'auto', ['group1;subgroup1', 'group2']),
        (
            [['sample1', 'sample2'], ['sample3'], ['sample4', 'sample5', 'sample6']],
            'display',
            ['sample1\nsample2', 'sample3', 'sample4\nsample5\nsample6'],
        ),
        (
            [['sample1', 'sample2'], ['sample3'], ['sample4', 'sample5', 'sample6']],
            'auto',
            ['sample1;sample2', 'sample3', 'sample4;sample5;sample6'],
        ),
        (
            [('sample1', 'sample2'), ('sample3',), ('sample4', 'sample5', 'sample6')],
            'display',
            ['sample1\nsample2', 'sample3', 'sample4\nsample5\nsample6'],
        ),
        (
            [('sample1', 'sample2'), ('sample3',), ('sample4', 'sample5', 'sample6')],
            'auto',
            ['sample1;sample2', 'sample3', 'sample4;sample5;sample6'],
        ),
    ],
)
def test_make_group_names(sample_grouping, mode, expected_output):
    assert make_group_names(sample_grouping, mode) == expected_output


def _legacy_parse_docstring(docstring: str):
    # verbatim copy of parse_docstring() before the fix for issue #339, kept as a back-compat reference
    docstring = re.sub(' +', ' ', docstring)
    split = docstring.split('\n\n')
    desc = split[0]
    params_str = split[1]
    free_text_match = r'["#\w\s\.\(\)\-\:%\*=@!\?\+\_,/`><\[\]' + "'" + ']'
    params_matches = list(
        re.finditer(r'^:param ([a-zA-Z_0-9]+):(' + free_text_match + r'+?)(?=^:.*:)', params_str, re.MULTILINE)
    )
    params = {match.group(1): match.group(2).replace('. ', '. \n') for match in params_matches}
    return desc, params


def _gui_exposed_functions():
    from rnalysis import enrichment, fastq, filtering

    funcs = {}
    for cls in (
        filtering.Filter,
        filtering.CountFilter,
        filtering.DESeqFilter,
        filtering.FoldChangeFilter,
        enrichment.FeatureSet,
        enrichment.RankedSet,
    ):
        # mirror gui.TabPage.get_all_actions(): every public callable attribute of the class.
        # EXCLUDED_FUNCS are deliberately kept (a stricter superset, and avoids importing Qt here).
        # Methods inherited from builtins (e.g. set.update) are not RNAlysis functions and have no reST params.
        for name in dir(cls):
            if name.startswith('_') or not callable(getattr(cls, name)):
                continue
            func = getattr(cls, name)
            if not getattr(func, '__module__', None) or not func.__module__.startswith('rnalysis'):
                continue
            funcs[f'{cls.__name__}.{name}'] = func
    for module in (fastq, enrichment):
        for name, func in inspect.getmembers(module, inspect.isfunction):
            if not name.startswith('_') and func.__module__ == module.__name__:
                funcs[f'{module.__name__}.{name}'] = func
    return funcs


def _raw_docstring(func):
    # the same preprocessing io.get_method_docstring() applies before calling parse_docstring()
    return inspect.cleandoc(inspect.getdoc(func) or '')


def test_parse_docstring_basic():
    docstring = (
        'Short description of the function.\n\n'
        ':param a: first param. It has two sentences.\n'
        ':type a: int\n'
        ':param b: second param\n'
        ':type b: str\n'
        ':return: something\n'
        ':rtype: int'
    )
    desc, params = parse_docstring(docstring)
    assert desc == 'Short description of the function.'
    assert params == {'a': ' first param. \nIt has two sentences.\n', 'b': ' second param\n'}


@pytest.mark.parametrize(
    'text',
    [
        "If 'first', puts NaNs at the beginning; if 'last', puts NaNs at the end",
        "for example {'method':'k-medoids', 'n_clusters':[3,5]}",
        'either a | b',
        'approximately ~5 reads',
        'x^2 + y^2',
        'costs $5',
        'concentration in µM — roughly',
        'a & b',
        'path\\to\\file',
        'semicolon; at the end;',
        'p < 0.05 and q > 0.1',
    ],
)
def test_parse_docstring_special_characters(text):
    docstring = f'Description.\n\n:param a: {text}\n:type a: str\n:param b: plain text\n:type b: int'
    _, params = parse_docstring(docstring)
    assert params == {'a': f' {text}\n', 'b': ' plain text\n'}


def test_parse_docstring_multiline_param():
    docstring = (
        'Description.\n\n'
        ':param a: first line; \n'
        'second line, continued. Third sentence.\n'
        'a real new line\n'
        ':type a: str\n'
        ':return: nothing'
    )
    _, params = parse_docstring(docstring)
    assert params == {'a': ' first line; \nsecond line, continued. \nThird sentence.\na real new line\n'}


class _BackslashContinuationDocstring:
    def method(self, a):
        """
        Description.

        :param a: first line; \
        second line
        :type a: str
        :return: nothing
        """


def test_parse_docstring_backslash_continuation():
    # the way real RNAlysis docstrings continue a line: Python consumes the trailing backslash and newline
    _, params = parse_docstring(_raw_docstring(_BackslashContinuationDocstring.method))
    assert params == {'a': ' first line; second line\n'}


@pytest.mark.parametrize(
    'tail,expected',
    [
        ('\n:return: a value', ' last; param\n'),
        ('\n:rtype: int', ' last; param\n'),
        ('\n:Examples:\n>>> f()', ' last; param\n'),
        ('', ' last; param'),
        ('\n\n:Examples:\n>>> f(a=1)\n1', ' last; param'),
        ('\n\n.. figure:: /figures/fig.png\n :align: center', ' last; param'),
    ],
)
def test_parse_docstring_last_param(tail, expected):
    docstring = 'Description.\n\n:type a: str\n:param a: last; param' + tail
    _, params = parse_docstring(docstring)
    assert params == {'a': expected}


def test_parse_docstring_params_in_later_paragraph():
    docstring = (
        'Description.\n\n'
        ':param a: first\n'
        ':type a: int\n\n'
        ':param b: second; after a blank line\n'
        ':type b: int\n'
        ':return: nothing'
    )
    desc, params = parse_docstring(docstring)
    assert desc == 'Description.'
    assert params == {'a': ' first\n', 'b': ' second; after a blank line\n'}


@pytest.mark.parametrize('docstring', ['Only a description.', 'Return self|value.', ''])
def test_parse_docstring_no_params_section(docstring):
    assert parse_docstring(docstring) == (docstring, {})


@pytest.mark.parametrize(
    'func_name,param_name,expected',
    [
        (
            'CountFilter.filter_low_reads',
            'threshold',
            ' The minimal number of reads (counts, rpm, rpkm, tpm, etc) a feature should have in at least '
            '`n_samples` samples in order not to be filtered out.\n',
        ),
        ('CountFilter.split_kmeans', 'n_clusters', ' The number of clusters the algorithm will seek.\n'),
        (
            'FeatureSet.go_enrichment',
            'organism',
            ' organism name or NCBI taxon ID for which the function will fetch GO annotations.\n',
        ),
    ],
)
def test_parse_docstring_real_snapshot(func_name, param_name, expected):
    _, params = parse_docstring(_raw_docstring(_gui_exposed_functions()[func_name]))
    assert params[param_name] == expected


def test_parse_docstring_backwards_compatible():
    # every parameter description that parsed before the fix for issue #339 must stay byte-identical
    n_compared = 0
    for qualname, func in _gui_exposed_functions().items():
        raw = _raw_docstring(func)
        if '\n\n' not in raw:
            continue
        legacy_desc, legacy_params = _legacy_parse_docstring(raw)
        desc, params = parse_docstring(raw)
        assert desc == legacy_desc, qualname
        for name, text in legacy_params.items():
            assert params.get(name) == text, f'{qualname}: {name}'
            n_compared += 1
    assert n_compared > 1000


def test_parse_docstring_covers_every_gui_param():
    # contract test (issue #339): every ':param name:' in a GUI-exposed function's docstring
    # must reach the GUI with non-empty help text
    missing = []
    for qualname, func in _gui_exposed_functions().items():
        raw = _raw_docstring(func)
        _, params = parse_docstring(raw)
        for name in re.findall(r'^[ \t]*:param ([a-zA-Z_0-9]+):', raw, re.MULTILINE):
            if not params.get(name, '').strip():
                missing.append(f'{qualname}: {name}')
    assert not missing, f'{len(missing)} parameters have no help text:\n' + '\n'.join(missing)
