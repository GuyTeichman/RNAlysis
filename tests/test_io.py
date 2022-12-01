import pytest

import rnalysis.utils.io
from rnalysis.utils import io
from rnalysis.utils.io import *
from rnalysis.utils.io import _format_ids_iter, _ensmbl_lookup_post_request


class MockResponse(object):
    def __init__(self, status_code: int = 200, url: str = 'http://httpbin.org/get', headers: dict = 'default',
                 text: str = '', json_output: dict = dict(), content: str = ''):
        self.status_code = status_code
        self.url = url
        self.headers = {'default': 'default'} if headers == 'default' else headers
        self.text = text
        self.ok = self.status_code == 200
        self._json = json_output
        self.content = bytes(content, 'utf8')

    def raise_for_status(self):
        if not self.ok:
            raise ConnectionError('request not ok')

    def json(self):
        return self._json

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        return


def test_load_csv_bad_input():
    invalid_input = 2349
    with pytest.raises(AssertionError):
        a = load_csv(invalid_input)


@pytest.mark.parametrize('pth', ("tests/test_files/test_load_csv.csv", "tests/test_files/test_load_csv.tsv",
                                 "tests/test_files/test_load_csv_tabs.txt",
                                 "tests/test_files/test_load_csv_other_sep.txt"))
def test_load_csv(pth):
    truth = pd.DataFrame({'idxcol': ['one', 'two', 'three'], 'othercol': [4, 5, 6]})
    truth.set_index('idxcol', inplace=True)
    loaded = load_csv(pth, 0)
    assert loaded.equals(truth)


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


@pytest.mark.parametrize("test_input,expected", [
    ('any', {'aspect a', 'aspect b', 'aspect c'}),
    ('Asp_B', {'aspect b'}),
    (['asp_b'], {'aspect b'}),
    (['a', 'z', 'c'], {'aspect a', 'aspect c', 'z'}),
    (['asp_b', 'c', 'A'], {'aspect a', 'aspect b', 'aspect c'}),
    (['aspect z'], {'aspect z'})
])
def test_golr_annotation_iterator_parse_go_aspects(monkeypatch, test_input, expected):
    go_dict = {'a': 'aspect a', 'asp_b': 'aspect b', 'c': 'aspect c', '_a_': 'aspect a'}
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
    ({'z', 'v'}, {'z', 'v'}),
    (None, set())
])
def test_golr_annotation_iterator_parse_evidence_types(monkeypatch, test_input, expected):
    ev_dict = {'a': 'eva', 'b': 'evb', 'c': 'evc', 'ab': {'eva', 'evb'}, 'bc': {'evb', 'evc'}, 'de': {'evd', 'eve'}}
    monkeypatch.setattr(GOlrAnnotationIterator, '_EVIDENCE_TYPE_DICT', ev_dict)
    assert GOlrAnnotationIterator._parse_evidence_types(test_input) == expected


def test_golr_annotation_iterator_get_n_annotations(monkeypatch):
    num_found_truth = 126311

    def fake_request(self, params, cached_filename):
        assert isinstance(self, GOlrAnnotationIterator)
        assert isinstance(params, dict)
        assert cached_filename == 'test.json'
        with open(f'tests/test_files/golr_header.txt') as f:
            return f.readline()

    monkeypatch.setattr(GOlrAnnotationIterator, '_golr_request', fake_request)
    monkeypatch.setattr(GOlrAnnotationIterator, '_generate_cached_filename', lambda self, start: 'test.json')
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


def remove_cached_test_file(cached_filename: str):
    try:
        os.remove(io.get_todays_cache_dir().joinpath(cached_filename))
    except FileNotFoundError:
        pass


def test_golr_annotation_iterator_golr_request(monkeypatch):
    cached_filename = 'test.json'
    remove_cached_test_file(cached_filename)

    correct_url = GOlrAnnotationIterator.URL
    correct_params = {'param': 'value', 'other_param': 'other_value'}

    def mock_get(url, params: dict):
        assert url == correct_url
        assert params == correct_params
        return MockResponse(text='the correct text')

    monkeypatch.setattr(requests, 'get', mock_get)
    monkeypatch.setattr(GOlrAnnotationIterator, '_generate_cached_filename', lambda self, start: 'test.json')
    assert GOlrAnnotationIterator._golr_request(correct_params, cached_filename) == 'the correct text'

    def mock_get_uncached(url, params: dict):
        raise AssertionError("This function should not be called if a cached file was found!")

    monkeypatch.setattr(requests, 'get', mock_get_uncached)
    try:
        assert GOlrAnnotationIterator._golr_request(correct_params, cached_filename) == 'the correct text'
    finally:
        remove_cached_test_file(cached_filename)

    def mock_get_failed(url, params: dict):
        assert url == correct_url
        assert params == correct_params
        return MockResponse(text='the correct text', status_code=404)

    monkeypatch.setattr(requests, 'get', mock_get_failed)
    try:
        with pytest.raises(ConnectionError):
            _ = GOlrAnnotationIterator._golr_request(correct_params)
    finally:
        remove_cached_test_file(cached_filename)


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

    def fake_request(self, params, cached_filename):
        assert isinstance(self, GOlrAnnotationIterator)
        assert params == truth_params
        assert cached_filename == 'test.json'
        with open(f'tests/test_files/golr_response.txt') as f:
            return f.readline()

    monkeypatch.setattr(GOlrAnnotationIterator, '_golr_request', fake_request)
    monkeypatch.setattr(GOlrAnnotationIterator, '_generate_cached_filename', lambda self, start: 'test.json')

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


def test_map_taxon_id_connectivity():
    assert map_taxon_id(6239) == (6239, 'Caenorhabditis elegans')
    assert map_taxon_id('canis lupus familiaris') == (9615, 'Canis lupus familiaris')
    with pytest.raises(ValueError):
        map_taxon_id('Lorem ipsum dolor sit amet')


def test_map_taxon_id(monkeypatch):
    taxon_name = 'c elegans'

    def mock_requests_get(url, params):
        assert url == 'https://rest.uniprot.org/taxonomy/search?'
        assert params == {'format': 'tsv', 'query': taxon_name}
        return MockResponse(text='Taxon Id\tMnemonic\tScientific name\tCommon name\tSynonym\tOther Names\tReviewed\t'
                                 'Rank\tLineage\tParent\tVirus hosts\n6239\tCAEEL\tCaenorhabditis elegans\t\t\t'
                                 'Caenorhabditis elegans (Maupas, 1900); Rhabditis elegans; Rhabditis elegans Maupas, '
                                 '1900; roundworm\treviewed\tSpecies\tEukaryota; Metazoa; Ecdysozoa; Nematoda; '
                                 'Chromadorea; Rhabditida; Rhabditina; Rhabditomorpha; Rhabditoidea; Rhabditidae; '
                                 'Peloderinae; Caenorhabditis\t6237\t\n')

    monkeypatch.setattr(requests, 'get', mock_requests_get)
    assert map_taxon_id(taxon_name) == (6239, 'Caenorhabditis elegans')


def test_map_taxon_id_no_results(monkeypatch):
    def mock_requests_get(url, params):
        return MockResponse(text='')

    monkeypatch.setattr(requests, 'get', mock_requests_get)
    with pytest.raises(ValueError):
        map_taxon_id('')


def test_map_taxon_id_multiple_results(monkeypatch):
    def mock_requests_get(url, params):
        return MockResponse(
            text='Taxon Id\tScientific name\n9615\tCanis lupus familiaris\n2509620\t'
                 'Wlobachia endosymbiont of Canis lupus familiaris\n990119\tCanis lupus x Canis lupus familiaris')

    monkeypatch.setattr(requests, 'get', mock_requests_get)
    assert map_taxon_id('') == (9615, 'Canis lupus familiaris')


def test_map_taxon_id_no_connection(monkeypatch):
    def mock_requests_get(url, params):
        return MockResponse(status_code=100)

    monkeypatch.setattr(requests, 'get', mock_requests_get)
    with pytest.raises(ConnectionError):
        map_taxon_id('name')


def test_ensmbl_lookup_post_request(monkeypatch):
    ids = ('id1', 'id2', 'id3')

    def mock_post_request(url, headers, data):
        assert url == 'https://rest.ensembl.org/lookup/id'
        assert headers == {"Content-Type": "application/json", "Accept": "application/json"}
        assert isinstance(data, str)
        assert json.loads(data) == {'ids': list(ids)}

        return MockResponse(json_output={this_id: {} for this_id in ids})

    monkeypatch.setattr(requests, 'post', mock_post_request)
    assert _ensmbl_lookup_post_request(ids) == {'id1': {}, 'id2': {}, 'id3': {}}


@pytest.mark.parametrize("gene_id_info,truth", [
    ({'id1': {'source': 'src1'}, 'id2': {'source': 'src1'}}, {'src1': {'id1', 'id2'}}),
    ({'id1': {'source': 'src1'}, 'id2': {'source': 'src1'}, 'id3': None}, {'src1': {'id1', 'id2'}}),
    ({'id1': {'source': 'src1'}, 'id2': {'source': 'src1'}, 'id3': {'source': 'src2'}},
     {'src1': {'id1', 'id2'}, 'src2': {'id3'}}),
    ({'id1': None, 'id2': None}, {}),
    ({}, {})
])
def test_infer_sources_from_gene_ids(monkeypatch, gene_id_info, truth):
    monkeypatch.setattr(io, '_ensmbl_lookup_post_request', lambda x: gene_id_info)
    assert infer_sources_from_gene_ids([]) == truth


@pytest.mark.parametrize("gene_id_info,truth", [
    ({'id1': {'species': 'c_elegans'}, 'id2': {'species': 'c_elegans'}, 'id3': None}, 'c elegans'),
    ({'id1': {'species': 'c_elegans'}, 'id2': {'species': 'm_musculus'}, 'id3': {'species': 'm_musculus'}},
     'm musculus')])
def test_infer_taxon_from_gene_ids(monkeypatch, gene_id_info, truth):
    monkeypatch.setattr(io, 'map_taxon_id', lambda x: x)
    monkeypatch.setattr(io, '_ensmbl_lookup_post_request', lambda x: gene_id_info)
    assert infer_taxon_from_gene_ids([])[0] == truth


def test_infer_taxon_from_gene_ids_no_species(monkeypatch):
    gene_id_info = {'id1': None, 'id2': None}
    monkeypatch.setattr(io, '_ensmbl_lookup_post_request', lambda x: gene_id_info)
    with pytest.raises(ValueError):
        infer_taxon_from_gene_ids([])


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


@pytest.mark.parametrize('id_type', ['UniProtKB', 'Entrez', 'WormBase'])
def test_map_gene_ids_to_same_set(id_type):
    mapper = map_gene_ids(['it', 'doesnt', 'matter', 'what', 'is', 'in', 'here'], id_type, id_type)
    assert mapper.mapping_dict is None
    for i in ['it', 'not', False, 42, 3.14]:
        assert i in mapper
        assert mapper[i] == i


@pytest.mark.parametrize('ids,map_from,map_to,req_from,req_to,req_query,txt,truth',
                         [(['P34544', 'Q27395', 'P12844'], 'UniProtKB', 'WormBase', 'ACC', 'WORMBASE_ID',
                           'P34544 Q27395 P12844',
                           'From\tTo\nP34544\tWBGene00019883\nQ27395\tWBGene00023497\nP12844\tWBGene00003515\n',
                           {'P34544': 'WBGene00019883', 'Q27395': 'WBGene00023497', 'P12844': 'WBGene00003515'}
                           )])
def test_map_gene_ids_request(monkeypatch, ids, map_from, map_to, req_from, req_to, req_query, txt, truth):
    def mock_get(url, params):
        assert url == 'https://www.uniprot.org/uploadlists/'
        assert params == {'from': req_from,
                          'to': req_to,
                          'format': 'tab',
                          'query': req_query,
                          'columns': 'id'}
        return MockResponse(text=txt)

    monkeypatch.setattr(requests, 'get', mock_get)
    res = map_gene_ids(ids, map_from, map_to)
    for gene_id in truth:
        assert res[gene_id] == truth[gene_id]


@pytest.mark.parametrize('ids,map_from,map_to,txt,rev_txt,truth',
                         [(['WBGene00000003', 'WBGene00000004'], 'WormBase', 'UniProtKB',
                           'From\tEntry\tAnnotation\nWBGene00000003\tQ19151\t110\nWBGene00000004\tA0A0K3AVL7\t57\nWBGene00000004\tO17395\t137.2\n'
                           , '', {'WBGene00000003': 'Q19151', 'WBGene00000004': 'O17395'}),
                          (
                              ['id1', 'id2'], 'UniProtKB', 'WormBase',
                              'From\tTo\nid1\tWBID1\nid2\tWBID2.2\nid2\tWBID2.1\n',
                              'From\tEntry\tAnnotation\nWBID1\tid1\t112.5\nWBID2.1\tid2\t112.5\nWBID2.2\tid2\t235\n'
                              , {'id1': 'WBID1', 'id2': 'WBID2.2'})
                          ])
def test_map_gene_ids_with_duplicates(monkeypatch, ids, map_from, map_to, txt, rev_txt, truth):
    def mock_abbrev_dict():
        d = {'WormBase': 'WormBase',
             'UniProtKB_to': 'UniProtKB',
             'UniProtKB_from': 'UniProtKB_AC-ID',
             'UniProtKB': 'UniProtKB'}
        return d, d

    def mock_get_mapping_results(api_url: str, from_db: str, to_db: str, ids: List[str], polling_interval: float,
                                 session, verbose):
        mock_abbrev_dict_to, mock_abbrev_dict_from = mock_abbrev_dict()
        if to_db == mock_abbrev_dict_to['UniProtKB_to']:
            return_txt = txt if map_to == 'UniProtKB' else rev_txt
        elif from_db == mock_abbrev_dict_from['UniProtKB_from']:
            return_txt = txt if map_from == 'UniProtKB' else rev_txt
        else:
            raise ValueError(to_db, from_db)
        return return_txt.split('\n')

    monkeypatch.setattr(io, '_get_id_abbreviation_dicts', mock_abbrev_dict)
    monkeypatch.setattr(io, 'get_mapping_results', mock_get_mapping_results)
    res = map_gene_ids(ids, map_from, map_to)
    for gene_id in truth:
        assert res[gene_id] == truth[gene_id]


def test_get_todays_cache_dir():
    today = date.today()
    today_str = str(today.year) + '_' + str(today.month).zfill(2) + '_' + str(today.day).zfill(2)
    cache_dir_truth = os.path.join(appdirs.user_cache_dir('RNAlysis'), today_str)
    assert cache_dir_truth == str(io.get_todays_cache_dir())


def test_load_cached_file():
    cached_filename = 'test.txt'
    remove_cached_test_file(cached_filename)

    cache_content_truth = "testing\n123"
    cache_dir = io.get_todays_cache_dir()
    path = os.path.join(cache_dir, cached_filename)

    assert load_cached_file(cached_filename) is None

    with open(path, 'x') as f:
        f.write(cache_content_truth)

    try:
        assert load_cached_file(cached_filename) == cache_content_truth
    finally:
        remove_cached_test_file(cached_filename)


def test_cache_file():
    cached_filename = 'test.txt'
    remove_cached_test_file(cached_filename)

    cache_content_truth = "testing\n123"
    cache_dir = io.get_todays_cache_dir()
    path = os.path.join(cache_dir, cached_filename)

    cache_file(cache_content_truth, cached_filename)
    try:
        with open(path, 'r') as f:
            assert f.read() == cache_content_truth
    finally:
        remove_cached_test_file(cached_filename)


@pytest.mark.parametrize("gene_set,expected_split", [
    ({1, 2, 3}, ['1', '2', '3']),
    ({'geneA', 'geneB', 'geneC', 'geneD'}, ["geneA", "geneB", "geneC", "geneD"])
])
def test_save_gene_set(gene_set, expected_split):
    pth = 'tests/test_files/tmp_saved_gene_set.txt'
    try:
        save_gene_set(gene_set, pth)
        with open(pth) as f:
            split = f.read().split('\n')
        assert sorted(split) == sorted(expected_split)
    finally:
        try:
            os.unlink(pth)
        except FileNotFoundError:
            pass


@pytest.mark.parametrize("args", [(6239,), (6239, 'all'), (6239, ['id1', 'id2'])])
def test_kegg_annotation_iterator_api(args):
    _ = KEGGAnnotationIterator(*args)


@pytest.mark.parametrize('arguments,url_truth', [('argument', 'https://rest.kegg.jp/operation/argument'),
                                                 (['arg1', 'arg2', 'arg3'],
                                                  'https://rest.kegg.jp/operation/arg1/arg2/arg3'),
                                                 (['argument'], 'https://rest.kegg.jp/operation/argument')], )
def test_kegg_annotation_iterator_kegg_request(monkeypatch, arguments, url_truth):
    truth = '{"sample": "json", "lorem": "ipsum"}'
    cached = []

    def mock_get_cached_file(filename):
        return None

    def mock_cache_file(content, filename):
        assert filename == 'cached_filename.csv'
        assert content == truth
        cached.append(True)

    def mock_get(url):
        assert url == url_truth
        return MockResponse(text=truth)

    monkeypatch.setattr(io, 'load_cached_file', mock_get_cached_file)
    monkeypatch.setattr(requests, 'get', mock_get)
    monkeypatch.setattr(io, 'cache_file', mock_cache_file)
    assert KEGGAnnotationIterator._kegg_request('operation', arguments, 'cached_filename.csv') == (
        truth, False)
    assert cached == [True]


def test_kegg_annotation_iterator_kegg_request_cached(monkeypatch):
    truth = {"sample": "json", "lorem": "ipsum"}

    def mock_get_cached_file(filename):
        return truth

    monkeypatch.setattr(io, 'load_cached_file', mock_get_cached_file)
    assert KEGGAnnotationIterator._kegg_request('operation', 'argument', 'cached_filename.csv') == (truth, True)


PATHWAY_NAMES_TRUTH = {'path:cel00010': 'Glycolysis / Gluconeogenesis - Caenorhabditis elegans (nematode)',
                       'path:cel00020': 'Citrate cycle (TCA cycle) - Caenorhabditis elegans (nematode)',
                       'path:cel00030': 'Pentose phosphate pathway - Caenorhabditis elegans (nematode)',
                       'path:cel00040': 'Pentose and glucuronate interconversions - Caenorhabditis elegans (nematode)',
                       'path:cel00051': 'Fructose and mannose metabolism - Caenorhabditis elegans (nematode)',
                       'path:cel00052': 'Galactose metabolism - Caenorhabditis elegans (nematode)', }


def test_kegg_annotation_iterator_get_pathways(monkeypatch):
    truth = (PATHWAY_NAMES_TRUTH, 6)
    organism_code = 'cel'

    def mock_kegg_request(self, operation, arguments):
        assert operation == 'list'
        assert arguments == ['pathway', organism_code]
        with open('tests/test_files/kegg_pathways.txt') as f:
            return f.read(), False

    monkeypatch.setattr(KEGGAnnotationIterator, '_kegg_request', mock_kegg_request)

    kegg = KEGGAnnotationIterator.__new__(KEGGAnnotationIterator)
    kegg.organism_code = organism_code
    assert kegg.get_pathways() == truth


def test_kegg_annotation_iterator_get_custom_pathways(monkeypatch):
    truth = {'path:cel00010': None, 'path:cel00030': None, 'path:cel00051': None}

    def mock_kegg_request(self, operation, arguments):
        return 'cel', False

    monkeypatch.setattr(KEGGAnnotationIterator, '_kegg_request', mock_kegg_request)

    kegg = KEGGAnnotationIterator(6239, [i for i in truth.keys()])
    assert kegg.pathway_names == truth


def test_kegg_annotation_iterator_get_pathway_annotations(monkeypatch):
    truth = {'path:cel00010': ['Glycolysis / Gluconeogenesis - Caenorhabditis elegans (nematode)',
                               {'cel:CELE_F14B4.2', 'cel:CELE_Y87G2A.8', 'cel:CELE_C50F4.2', 'cel:CELE_Y71H10A.1'}],
             'path:cel00020': ['Citrate cycle (TCA cycle) - Caenorhabditis elegans (nematode)',
                               {'cel:CELE_T20G5.2', 'cel:CELE_B0365.1', 'cel:CELE_D1005.1'}],
             'path:cel00030': ['Pentose phosphate pathway - Caenorhabditis elegans (nematode)',
                               {'cel:CELE_Y87G2A.8', 'cel:CELE_B0035.5'}],
             'path:cel00040': ['Pentose and glucuronate interconversions - Caenorhabditis elegans (nematode)',
                               {'cel:CELE_Y105E8B.9', 'cel:CELE_B0310.5', 'cel:CELE_T04H1.7', 'cel:CELE_T04H1.8'}],
             'path:cel00051': ['Fructose and mannose metabolism - Caenorhabditis elegans (nematode)',
                               {'cel:CELE_C05C8.7', 'cel:CELE_ZK632.4'}],
             'path:cel00052': ['Galactose metabolism - Caenorhabditis elegans (nematode)', {'cel:CELE_C01B4.6'}]}
    args_truth = ['path:cel00010+path:cel00020+path:cel00030', 'path:cel00040+path:cel00051+path:cel00052']

    def mock_kegg_request(self, operation, arguments, fname):
        assert operation == 'get'
        assert arguments == args_truth[0] or arguments == args_truth[1]
        if arguments == args_truth[0]:
            pth = 'tests/test_files/kegg_annotation_1of2.txt'
        else:
            pth = 'tests/test_files/kegg_annotation_2of2.txt'
        with open(pth) as f:
            return f.read(), False

    monkeypatch.setattr(KEGGAnnotationIterator, '_kegg_request', mock_kegg_request)
    monkeypatch.setattr(KEGGAnnotationIterator, 'REQ_MAX_ENTRIES', 3)
    kegg = KEGGAnnotationIterator.__new__(KEGGAnnotationIterator)
    kegg.pathway_annotations = None
    kegg.taxon_id = 6239
    kegg.organism_code = 'cel'
    kegg.pathway_names = PATHWAY_NAMES_TRUTH
    assert {key: [name, ann] for key, name, ann in kegg.get_pathway_annotations()} == truth


def test_kegg_annotation_iterator_get_pathway_annotations_cached():
    annotation = {'a': 1, 'b': 2, 'c': 3}
    names = {'a': 'namea', 'b': 'nameb', 'c': 'namec'}
    truth = {'a': ['namea', 1], 'b': ['nameb', 2], 'c': ['namec', 3]}
    kegg = KEGGAnnotationIterator.__new__(KEGGAnnotationIterator)
    kegg.pathway_annotations = annotation
    kegg.pathway_names = names
    assert {key: [name, ann] for key, name, ann in kegg.get_pathway_annotations()} == truth


@pytest.mark.parametrize('taxon_id,truth', [
    (6239, 'cel'),
    (270351, 'maqu'),
    (9606, 'hsa'),
    (6238, 'cbr')
])
def test_kegg_annotation_iterator_get_kegg_organism_code(monkeypatch, taxon_id, truth):
    def mock_get_tree():
        with open('tests/test_files/kegg_taxon_tree_truth.json') as f:
            return json.loads(f.read())

    monkeypatch.setattr(KEGGAnnotationIterator, '_get_taxon_tree', mock_get_tree)

    assert KEGGAnnotationIterator.get_kegg_organism_code(taxon_id) == truth


def test_kegg_annotation_iterator_get_taxon_tree(monkeypatch):
    truth = {"sample": "json", "lorem": "ipsum"}
    truth_text = '{"sample": "json", "lorem": "ipsum"}'
    cached = []

    def mock_get_cached_file(filename):
        return None

    def mock_cache_file(content, filename):
        assert filename == 'kegg_taxon_tree.json'
        assert content == truth_text
        cached.append(True)

    def mock_get(url, params):
        assert url == 'https://www.genome.jp/kegg-bin/download_htext?htext=br08610'
        assert params == {'format': 'json'}
        return MockResponse(content=truth_text)

    monkeypatch.setattr(io, 'load_cached_file', mock_get_cached_file)
    monkeypatch.setattr(requests, 'get', mock_get)
    monkeypatch.setattr(io, 'cache_file', mock_cache_file)

    assert KEGGAnnotationIterator._get_taxon_tree() == truth
    assert cached == [True]


def test_kegg_annotation_iterator_get_taxon_tree_cached(monkeypatch):
    truth = {"sample": "json", "lorem": "ipsum"}

    def mock_get_cached_file(filename):
        return '{"sample":"json","lorem":"ipsum"}'

    monkeypatch.setattr(io, 'load_cached_file', mock_get_cached_file)
    assert KEGGAnnotationIterator._get_taxon_tree() == truth


def test_get_gui_videos(monkeypatch):
    mock_dir = Path('tests/.pytest_cache/video_test_cache')
    video_files = ('first_vid.webp', 'second_vid.webp', 'third_vid.webp')
    visited = []

    def mock_get(url, params):
        base_url = 'https://github.com/GuyTeichman/RNAlysis/blob/master/rnalysis/gui/videos/'
        assert url.startswith(base_url)
        visited.append(url.split(base_url)[1])
        return MockResponse(content='content')

    monkeypatch.setattr(requests, 'get', mock_get)
    monkeypatch.setattr(rnalysis.utils.io, 'get_tutorial_videos_dir',
                        lambda: mock_dir)

    try:
        list(get_gui_videos(video_files))
        assert sorted(visited) == ['first_vid.webp', 'second_vid.webp', 'third_vid.webp']
        assert mock_dir.exists()
        file_count = 0
        for file in mock_dir.iterdir():
            assert file.name in visited
            with open(file) as f:
                assert f.read() == 'content'
            file_count += 1
        assert file_count == len(video_files)
    finally:
        shutil.rmtree(mock_dir, ignore_errors=True)


class MockProcess:
    def __init__(self, returncode: int):
        self.stdout = [b'things', b'to', b'print']
        self.returncode = returncode

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        return


@pytest.mark.parametrize('r_path,expected', [
    ('auto', ['Rscript', "tests/test_files/test_r_script.R"]),
    ('D:/Program Files/R', ["D:/Program Files/R/bin/Rscript", "tests/test_files/test_r_script.R"])
])
def test_run_r_script(monkeypatch, r_path, expected):
    script_path = 'tests/test_files/test_r_script.R'

    def mock_popen(process, stdout, stderr, shell=False):
        assert process == expected
        return MockProcess(0)

    monkeypatch.setattr(subprocess, 'Popen', mock_popen)
    monkeypatch.setattr(subprocess, 'run', lambda *args, **kwargs: MockProcess(0))

    run_r_script(script_path, r_path)


def test_run_r_script_not_installed(monkeypatch):
    def mock_run(process, stdout, stderr, shell=False):
        return MockProcess(1)

    monkeypatch.setattr(subprocess, 'run', mock_run)
    with pytest.raises(FileNotFoundError):
        run_r_script('tests/test_files/test_r_script.R')


@pytest.mark.parametrize("stdout", [True, False])
@pytest.mark.parametrize("stderr", [True, False])
def test_run_subprocess(stdout, stderr):
    run_subprocess(['python3', '--version'], stdout, stderr)


@pytest.mark.parametrize('this_version,response,expected', [
    ('3.0.0', MockResponse(status_code=503), False),
    ('3.0.0', MockResponse(text='{"info":{"author":"Guy Teichman","version":"3.1.0"}}'), True),
    ('3.1.1', MockResponse(text='{"info":{"author":"Guy Teichman","version":"3.1.0"}}'), False),
    ('3.1.1', MockResponse(text='{"info":{"author":"Guy Teichman","version":"3.1.1"}}'), False),
])
def test_is_rnalysis_outdated(monkeypatch, this_version, response, expected):
    monkeypatch.setattr(requests, 'get', lambda *args, **kwargs: response)
    monkeypatch.setattr(rnalysis.utils.io, '__version__', this_version)
    assert is_rnalysis_outdated() == expected


@pytest.mark.parametrize('feature_type,truth', [
    ('gene', {
        'ENSG00000168671': 'protein_coding',
        'ENSG00000249641': 'antisense_RNA',
        'ENSG00000123364': 'protein_coding',
        'ENSG00000123407': 'protein_coding',
        'ENSG00000123388': 'protein_coding',
        'ENSG00000180818': 'protein_coding',
        'ENSG00000197757': 'protein_coding',
        'ENSG00000180806': 'protein_coding',
        'ENSG00000037965': 'protein_coding',
        'ENSG00000198353': 'protein_coding',
        'ENSG00000172789': 'other'}),
    ('transcript', {
        "ENST00000282507": "protein_coding",
        "ENST00000513300": "protein_coding",
        "ENST00000504685": "nonsense_mediated_decay",
        "ENST00000504954": "processed_transcript",
        "ENST00000515131": "protein_coding",
        "ENST00000512916": "antisense_RNA",
        "ENST00000243056": "protein_coding",
        "ENST00000243103": "protein_coding",
        "ENST00000546378": "protein_coding",
        "ENST00000243082": "protein_coding",
        "ENST00000515593": "protein_coding",
        "ENST00000303460": "protein_coding",
        "ENST00000511575": "processed_transcript",
        "ENST00000514415": "processed_transcript",
        "ENST00000504315": "protein_coding",
        "ENST00000509328": "protein_coding",
        "ENST00000394331": "protein_coding",
        "ENST00000243108": "protein_coding",
        "ENST00000513413": "processed_transcript",
        "ENST00000504557": "processed_transcript",
        "ENST00000508190": "protein_coding",
        "ENST00000303450": "protein_coding",
        "ENST00000040584": "protein_coding",
        "ENST00000303406": "protein_coding",
        "ENST00000507650": "processed_transcript",
        "ENST00000430889": "protein_coding",
        "ENST00000312492": "other", })
])
def test_map_gene_to_attr(feature_type, truth):
    gtf_path = 'tests/test_files/kallisto_tests/transcripts.gtf'
    res = map_gene_to_attr(gtf_path, feature_type + '_biotype', feature_type, False, False, False)
    assert res == truth
