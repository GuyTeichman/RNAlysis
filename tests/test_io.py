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


@pytest.mark.parametrize("test_input,expected", [
    ('any', {'aspect a', 'aspect b', 'aspect c'}),
    ('B', {'aspect b'}),
    (['b'], {'aspect b'}),
    (['a', 'z', 'c'], {'aspect a', 'aspect c', 'z'}),
    (['b', 'c', 'A'], {'aspect a', 'aspect b', 'aspect c'}),
    (['aspect z'], {'aspect z'})
])
def test_golr_annotation_iterator_parse_go_aspects(monkeypatch, test_input, expected):
    go_dict = {'a': 'aspect a', 'b': 'aspect b', 'c': 'aspect c', '_a_': 'aspect a'}
    monkeypatch.setattr(GOlrAnnotationIterator, '_ASPECTS_DICT', go_dict)
    assert GOlrAnnotationIterator._parse_go_aspects(test_input) == expected


def test_golr_annotation_iterator_api(monkeypatch):
    def null_method(self):
        pass

    def parse_method(self, param):
        return set()

    monkeypatch.setattr(GOlrAnnotationIterator, '_generate_query', null_method)
    monkeypatch.setattr(GOlrAnnotationIterator, '_get_n_annotations', null_method)
    monkeypatch.setattr(GOlrAnnotationIterator, '_parse_go_aspects', parse_method)
    monkeypatch.setattr(GOlrAnnotationIterator, '_parse_evidence_types', parse_method)

    golr = GOlrAnnotationIterator(1234)


@pytest.mark.parametrize("test_input,expected", [
    ('any', {'eva', 'evb', 'evc', 'evd', 'eve'}),
    ('bc', {'evb', 'evc'}),
    ('c', {'evc'}),
    ({'a', 'bc', 'f'}, {'eva', 'evb', 'evc', 'f'}),
    ({'a', 'ab'}, {'eva', 'evb'}),
    ({'z', 'v'}, {'z', 'v'})
])
def test_golr_annotation_iterator_parse_evidence_types(monkeypatch, test_input, expected):
    ev_dict = {'a': 'eva', 'b': 'evb', 'c': 'evc', 'ab': {'eva', 'evb'}, 'bc': {'evb', 'evc'}, 'de': {'evd', 'eve'}}
    monkeypatch.setattr(GOlrAnnotationIterator, '_EVIDENCE_TYPE_DICT', ev_dict)
    assert GOlrAnnotationIterator._parse_evidence_types(test_input) == expected


def test_golr_annotation_iterator_get_n_annotations(monkeypatch):
    num_found_truth = 126311

    def fake_request(self, params):
        assert isinstance(self, GOlrAnnotationIterator)
        assert isinstance(params, dict)
        with open(f'tests/test_files/golr_header.txt') as f:
            return f.readline()

    monkeypatch.setattr(GOlrAnnotationIterator, '_golr_request', fake_request)
    golr = GOlrAnnotationIterator.__new__(GOlrAnnotationIterator)
    golr.default_params = {}
    assert golr._get_n_annotations() == num_found_truth


def test_golr_annotation_iterator_generate_query():
    golr = GOlrAnnotationIterator.__new__(GOlrAnnotationIterator)
    golr.aspects = {'P', 'C'}
    golr.databases = {'DB1', 'DB2'}
    golr.evidence_types = {'IEA', 'IMP'}
    golr.excluded_databases = set()
    golr.excluded_evidence_types = {'EXP', 'IDA'}
    golr.excluded_qualifiers = {'not_a'}
    golr.qualifiers = set()
    golr.taxon_id = 6239

    aspects_iter = iter(golr.aspects)
    db_iter = iter(golr.databases)
    evidence_iter = iter(golr.evidence_types)

    query_truth = ['document_category:"annotation"', 'taxon:"NCBITaxon:6239"',
                   f'source:"{next(db_iter)}" OR source:"{next(db_iter)}"',
                   f'evidence_type:"{next(evidence_iter)}" OR evidence_type:"{next(evidence_iter)}"',
                   '-qualifier:"not_a"', '-evidence_type:"EXP"',
                   '-evidence_type:"IDA"', f'aspect:"{next(aspects_iter)}" OR aspect:"{next(aspects_iter)}"']

    assert sorted(golr._generate_query()) == sorted(query_truth)


def test_golr_annotation_iterator_golr_request_connectivity(monkeypatch):
    fake_params = {'param': 'value', 'other_param': 'other_value'}
    assert isinstance(GOlrAnnotationIterator._golr_request(fake_params), str)


def test_golr_annotation_iterator_golr_request(monkeypatch):
    correct_url = GOlrAnnotationIterator.URL
    correct_params = {'param': 'value', 'other_param': 'other_value'}

    def mock_get(url, params: dict):
        assert url == correct_url
        assert params == correct_params
        return MockResponse(text='the correct text')

    monkeypatch.setattr(requests, 'get', mock_get)
    assert GOlrAnnotationIterator._golr_request(correct_params) == 'the correct text'

    def mock_get_failed(url, params: dict):
        assert url == correct_url
        assert params == correct_params
        return MockResponse(text='the correct text', status_code=404)

    monkeypatch.setattr(requests, 'get', mock_get_failed)
    with pytest.raises(ConnectionError):
        _ = GOlrAnnotationIterator._golr_request(correct_params)


def test_golr_annotation_iterator_parsing(monkeypatch):
    truth_params = {
        "q": "*:*",
        "wt": "json",  # return format
        "rows": 5,  # how many rows to return
        # how many annotations to fetch (fetch 0 to find n_annotations, then fetch in iter_size increments
        "start": 0,  # from which annotation number to start fetching
        "fq": [f'document_category:"annotation"', 'taxon:"NCBITaxon:6239"'],  # search query
        "fl": "source,bioentity_internal_id,annotation_class",  # fields
        "omitHeader": 'true'}

    records_truth = [{"source": "WB", "bioentity_internal_id": "WBGene00011482", "annotation_class": "GO:0003923"},
                     {"source": "WB", "bioentity_internal_id": "WBGene00011482", "annotation_class": "GO:0016255"},
                     {"source": "WB", "bioentity_internal_id": "WBGene00011481", "annotation_class": "GO:0004190"},
                     {"source": "WB", "bioentity_internal_id": "WBGene00011481", "annotation_class": "GO:0005783"},
                     {"source": "WB", "bioentity_internal_id": "WBGene00011481", "annotation_class": "GO:0005789"}]

    def fake_request(self, params):
        assert isinstance(self, GOlrAnnotationIterator)
        assert params == truth_params
        with open(f'tests/test_files/golr_response.txt') as f:
            return f.readline()

    monkeypatch.setattr(GOlrAnnotationIterator, '_golr_request', fake_request)

    request_params = {
        "q": "*:*",
        "wt": "json",  # return format
        "rows": 5,  # how many rows to return
        # how many annotations to fetch (fetch 0 to find n_annotations, then fetch in iter_size increments
        "fq": [f'document_category:"annotation"', 'taxon:"NCBITaxon:6239"'],  # search query
        "fl": "source,bioentity_internal_id,annotation_class"}  # fields

    golr = GOlrAnnotationIterator.__new__(GOlrAnnotationIterator)
    golr.default_params = request_params
    golr.iter_size = 5
    golr.n_annotations = 5
    records = [i for i in golr]
    assert len(records) == len(records_truth)
    for record, true_record in zip(records, records_truth):
        assert record == true_record


def test_map_taxon_id_parsing(monkeypatch):
    assert False


def test_map_taxon_id_no_results(monkeypatch):
    assert False


def test_map_taxon_id_multiple_results(monkeypatch):
    assert False


def test_map_taxon_id_no_connection(monkeypatch):
    assert False


def _load_id_abbreviation_dict():
    assert False
