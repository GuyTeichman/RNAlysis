import pytest
from rnalysis.utils.io import *
from rnalysis.utils.io import _format_ids_iter


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


def test_load_csv_bad_input():
    invalid_input = 2349
    with pytest.raises(AssertionError):
        a = load_csv(invalid_input)


def test_load_csv():
    truth = pd.DataFrame({'idxcol': ['one', 'two', 'three'], 'othercol': [4, 5, 6]})
    truth.set_index('idxcol', inplace=True)
    pth = "tests/test_files/test_load_csv.csv"
    assert (truth == load_csv(pth, 0)).all().all()


def test_load_csv_drop_columns():
    loaded = load_csv('tests/test_files/counted.csv', 0, drop_columns='cond1')
    print(loaded)
    assert list(loaded.columns) == ['cond2', 'cond3', 'cond4']

    loaded = load_csv('tests/test_files/counted.csv', 0, drop_columns=['cond2', 'cond4'])
    assert list(loaded.columns) == ['cond1', 'cond3']

    with pytest.raises(IndexError):
        load_csv('tests/test_files/counted.csv', 0, drop_columns=['cond1', 'cond6'])


def test_save_csv():
    try:
        df = pd.read_csv('tests/test_files/enrichment_hypergeometric_res.csv', index_col=0)
        save_csv(df, 'tests/test_files/tmp_test_save_csv.csv')
        df_loaded = pd.read_csv('tests/test_files/tmp_test_save_csv.csv', index_col=0)
        assert df.equals(df_loaded)
        df = pd.read_csv('tests/test_files/enrichment_hypergeometric_res.csv')
        save_csv(df, 'tests/test_files/tmp_test_save_csv.csv', '_2', index=False)
        df_loaded = pd.read_csv('tests/test_files/tmp_test_save_csv_2.csv', index_col=0)
        df = pd.read_csv('tests/test_files/enrichment_hypergeometric_res.csv', index_col=0)
        assert df.equals(df_loaded)

    except Exception as e:
        raise e
    finally:
        try:
            os.remove('tests/test_files/tmp_test_save_csv.csv')
            os.remove('tests/test_files/tmp_test_save_csv_2.csv')
        except:
            pass


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


def test_format_ids_iter():
    assert list(_format_ids_iter('one two three')) == ['one two three']
    assert list(_format_ids_iter(123)) == ['123']
    assert list(_format_ids_iter(['one', ' two', 'three; ', 'four'])) == ['one  two three;  four']
    assert list(_format_ids_iter(['1', 'two', '3', '4', 'five', '6', '7'], 3)) == ['1 two 3', '4 five 6', '7']


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


def _load_id_abbreviation_dict():
    assert False
