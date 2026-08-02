import gzip

import pytest

from rnalysis.utils.genome_annotation import *


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


@pytest.mark.parametrize('gtf_path,feature_type,len_method,truth', [
    ('tests/test_files/test_gtf_wormbase.gtf', 'gene', 'mean',
     {"WBGene00197333": 150,
      "WBGene00198386": 150,
      "WBGene00015153": 1177.5,
      "WBGene00002061": 6137 / 6}),
    ('tests/test_files/test_gtf_wormbase.gtf', 'gene', 'median', {
        "WBGene00197333": 150,
        "WBGene00198386": 150,
        "WBGene00015153": 1177.5,
        "WBGene00002061": 1021.5}),
    ('tests/test_files/test_gtf_wormbase.gtf', 'gene', 'geometric_mean', {
        "WBGene00197333": 150,
        "WBGene00198386": 150,
        "WBGene00015153": 1167.9708900482055,
        "WBGene00002061": 1019.8567195006932}),
    ('tests/test_files/test_gtf_wormbase.gtf', 'gene', 'merged_exons',
     {"WBGene00197333": 150,
      "WBGene00198386": 150,
      "WBGene00015153": 1327,
      "WBGene00002061": 1108}
     ),
    ('tests/test_files/test_gtf_wormbase.gtf', 'gene', 'min',
     {"WBGene00197333": 150,
      "WBGene00198386": 150,
      "WBGene00015153": 1028,
      "WBGene00002061": 940}),
    ('tests/test_files/test_gtf_wormbase.gtf', 'transcript', 'mean',
     {'cTel3X.2': 150,
      'cTel3X.3': 150,
      'B0348.5a': 1327,
      'B0348.5b': 1028,
      'B0348.6b.2': 1108,
      'B0348.6a.2': 940,
      'B0348.6b.1': 949,
      'B0348.6c.2': 946,
      'B0348.6a.1': 1094,
      'B0348.6c.1': 1100,
      }),
    ('tests/test_files/test_gff_wormbase.gff3', 'gene', 'geometric_mean',
     {"WBGene00015153": 1167.9708900482055,
      "WBGene00002061": 1019.8567195006932}
     ),
    ('tests/test_files/test_gff_wormbase.gff3', 'gene', 'max',
     {"WBGene00015153": 1327,
      "WBGene00002061": 1108}
     ),
    ('tests/test_files/test_gff_wormbase.gff3', 'gene', 'merged_exons',
     {"WBGene00015153": 1327,
      "WBGene00002061": 1108}),
])
def test_get_genomic_feature_lengths(gtf_path, feature_type, len_method, truth):
    res = get_genomic_feature_lengths(gtf_path, feature_type, len_method)
    assert res.keys() == truth.keys()
    for key in res:
        assert np.isclose(res[key], truth[key])


# ======================================================================================================================
# Characterization tests (epic #67 / PR-A, issue #150).
#
# These pin the EXACT current behaviour of the three annotation functions so the planned Polars rewrite (PR-B, #151)
# can be held to byte-identical results. Several assertions deliberately lock in known quirks/bugs that are only
# scheduled to be FIXED in PR-C (#152) -- they are marked below. Do not "fix" them here: if the behaviour changes,
# that belongs in #152 with a HISTORY.rst note, and these tests are what will flag it.
#
# Small, hand-verifiable fixtures (two genes; gene A has two transcripts of length 200 & 100 so averaging methods
# diverge; gene B has one transcript of length 100). GENCODE inline versions mirror the Ensembl separate-version
# attributes, giving the two formats deliberately parallel outputs.
# ======================================================================================================================
ENSEMBL_GTF = 'tests/test_files/test_gtf_ensembl.gtf'  # bare IDs + gene_version/transcript_version, gene_biotype
GENCODE_GTF = 'tests/test_files/test_gtf_gencode.gtf'  # inline-versioned IDs, NO version attrs, gene_type
MULTIPARENT_GFF3 = 'tests/test_files/test_gff_multiparent.gff3'  # mRNA feature, biotype=, comma-separated Parent=a,b
WORMBASE_GFF3 = 'tests/test_files/test_gff_wormbase.gff3'  # bare WBGene/dotted mRNA IDs, no type: prefix


def _malformed_gtf(tmp_path):
    """Valid header + gene/transcript/exon lines, plus one line with the wrong number of columns."""
    pth = tmp_path / 'malformed.gtf'
    pth.write_text(
        '#comment\n'
        '1\tsrc\tgene\t100\t200\t.\t+\t.\tgene_id "G1"; gene_name "GENEA";\n'
        '1\tsrc\ttranscript\t100\t200\t.\t+\t.\tgene_id "G1"; transcript_id "T1";\n'
        '1\tsrc\texon\t100\t200\t.\t+\t.\tgene_id "G1"; transcript_id "T1";\n'
        '1\tsrc\texon\t100\t200\t.\t+\n',  # only 7 columns -> malformed (a valid line has 9)
        encoding='utf-8', newline='\n')
    return pth


# ------------------------------------------------------------------ map_transcripts_to_genes (currently uncovered)
@pytest.mark.parametrize('use_name,use_version,split_ids,truth', [
    # Ensembl: bare IDs -> split_ids is a no-op; use_version drives the version suffix.
    (False, False, False, {'ENST00000000001': 'ENSG00000000001', 'ENST00000000002': 'ENSG00000000001',
                           'ENST00000000003': 'ENSG00000000002'}),
    (False, False, True, {'ENST00000000001': 'ENSG00000000001', 'ENST00000000002': 'ENSG00000000001',
                          'ENST00000000003': 'ENSG00000000002'}),
    (False, True, False, {'ENST00000000001.2': 'ENSG00000000001.3', 'ENST00000000002.1': 'ENSG00000000001.3',
                          'ENST00000000003.4': 'ENSG00000000002.1'}),
    (False, True, True, {'ENST00000000001.2': 'ENSG00000000001.3', 'ENST00000000002.1': 'ENSG00000000001.3',
                         'ENST00000000003.4': 'ENSG00000000002.1'}),
    (True, False, False, {'ENST00000000001': 'GENEA', 'ENST00000000002': 'GENEA', 'ENST00000000003': 'GENEB'}),
    (True, False, True, {'ENST00000000001': 'GENEA', 'ENST00000000002': 'GENEA', 'ENST00000000003': 'GENEB'}),
    (True, True, False, {'ENST00000000001.2': 'GENEA', 'ENST00000000002.1': 'GENEA', 'ENST00000000003.4': 'GENEB'}),
    (True, True, True, {'ENST00000000001.2': 'GENEA', 'ENST00000000002.1': 'GENEA', 'ENST00000000003.4': 'GENEB'}),
])
def test_map_transcripts_to_genes_ensembl(use_name, use_version, split_ids, truth):
    assert map_transcripts_to_genes(ENSEMBL_GTF, use_name, use_version, split_ids) == truth


@pytest.mark.parametrize('use_name,use_version,split_ids,truth', [
    # GENCODE: inline-versioned IDs -> use_version is inert (no *_version attrs); split_ids strips the inline version.
    (False, False, False, {'ENST00000000001.2': 'ENSG00000000001.3', 'ENST00000000002.1': 'ENSG00000000001.3',
                           'ENST00000000003.4': 'ENSG00000000002.1'}),
    (False, False, True, {'ENST00000000001': 'ENSG00000000001', 'ENST00000000002': 'ENSG00000000001',
                          'ENST00000000003': 'ENSG00000000002'}),
    (False, True, False, {'ENST00000000001.2': 'ENSG00000000001.3', 'ENST00000000002.1': 'ENSG00000000001.3',
                          'ENST00000000003.4': 'ENSG00000000002.1'}),
    (False, True, True, {'ENST00000000001': 'ENSG00000000001', 'ENST00000000002': 'ENSG00000000001',
                         'ENST00000000003': 'ENSG00000000002'}),
    (True, False, False, {'ENST00000000001.2': 'GENEA', 'ENST00000000002.1': 'GENEA', 'ENST00000000003.4': 'GENEB'}),
    (True, False, True, {'ENST00000000001': 'GENEA', 'ENST00000000002': 'GENEA', 'ENST00000000003': 'GENEB'}),
    (True, True, False, {'ENST00000000001.2': 'GENEA', 'ENST00000000002.1': 'GENEA', 'ENST00000000003.4': 'GENEB'}),
    (True, True, True, {'ENST00000000001': 'GENEA', 'ENST00000000002': 'GENEA', 'ENST00000000003': 'GENEB'}),
])
def test_map_transcripts_to_genes_gencode(use_name, use_version, split_ids, truth):
    assert map_transcripts_to_genes(GENCODE_GTF, use_name, use_version, split_ids) == truth


def test_map_transcripts_to_genes_wormbase_split_ids_collapse():
    # QUIRK (fix in #152): the DEFAULT split_ids=True splits WormBase IDs on '.', but those dots are part of the ID
    # (cTel3X.2, B0348.5a), so transcripts collapse onto their prefix and first-wins keeps only 2 of 10 entries.
    assert map_transcripts_to_genes('tests/test_files/test_gtf_wormbase.gtf') == {
        'cTel3X': 'WBGene00197333', 'B0348': 'WBGene00015153'}
    # With split_ids=False every transcript is preserved (the intended mapping for this file).
    assert map_transcripts_to_genes('tests/test_files/test_gtf_wormbase.gtf',
                                    use_name=False, use_version=False, split_ids=False) == {
               'cTel3X.2': 'WBGene00197333', 'cTel3X.3': 'WBGene00198386',
               'B0348.5a': 'WBGene00015153', 'B0348.5b': 'WBGene00015153',
               'B0348.6b.2': 'WBGene00002061', 'B0348.6a.2': 'WBGene00002061', 'B0348.6b.1': 'WBGene00002061',
               'B0348.6c.2': 'WBGene00002061', 'B0348.6a.1': 'WBGene00002061', 'B0348.6c.1': 'WBGene00002061'}
    # use_name=True maps the (collapsed) key to the gene NAME.
    assert map_transcripts_to_genes('tests/test_files/test_gtf_wormbase.gtf', use_name=True) == {
        'cTel3X': 'cTel3X.2', 'B0348': 'B0348.5'}


def test_map_transcripts_to_genes_kallisto_realfile():
    # Ensembl-style human file: split_ids does not collapse (IDs have no internal dots); use_version toggles the
    # inline .N suffix. Spot-check size + representative entries rather than the full 27-entry dict.
    versioned = map_transcripts_to_genes('tests/test_files/kallisto_tests/transcripts.gtf')
    assert len(versioned) == 27
    assert versioned['ENST00000282507.7'] == 'ENSG00000168671.9'
    unversioned = map_transcripts_to_genes('tests/test_files/kallisto_tests/transcripts.gtf', use_version=False)
    assert len(unversioned) == 27
    assert unversioned['ENST00000282507'] == 'ENSG00000168671'


def test_map_transcripts_to_genes_malformed_skips_and_warns(tmp_path):
    # #152: previously RAISED ValueError on a bad-column line; now unified on skip-and-warn like the other
    # two functions. The one 7-column line is dropped (with a warning) and the valid transcript still maps.
    with pytest.warns(UserWarning, match='malformed line'):
        assert map_transcripts_to_genes(_malformed_gtf(tmp_path), use_version=False, split_ids=False) == {'T1': 'G1'}


# ------------------------------------------------------------------ get_genomic_feature_lengths (extra coverage)
@pytest.mark.parametrize('gtf_path,gene_keys', [(ENSEMBL_GTF, ('ENSG00000000001', 'ENSG00000000002')),
                                                (GENCODE_GTF, ('ENSG00000000001.3', 'ENSG00000000002.1'))])
@pytest.mark.parametrize('method,gene_a_len', [('mean', 150), ('median', 150), ('max', 200), ('min', 100),
                                               ('geometric_mean', 141.42135623730945), ('merged_exons', 200)])
def test_get_genomic_feature_lengths_gene_all_methods(gtf_path, gene_keys, method, gene_a_len):
    res = get_genomic_feature_lengths(gtf_path, 'gene', method)
    truth = {gene_keys[0]: gene_a_len, gene_keys[1]: 100}
    assert res.keys() == truth.keys()
    for key in res:
        assert np.isclose(res[key], truth[key])


@pytest.mark.parametrize('gtf_path,truth', [
    (ENSEMBL_GTF, {'ENST00000000001': 200, 'ENST00000000002': 100, 'ENST00000000003': 100}),
    (GENCODE_GTF, {'ENST00000000001.2': 200, 'ENST00000000002.1': 100, 'ENST00000000003.4': 100}),
])
def test_get_genomic_feature_lengths_transcript_method_ignored(gtf_path, truth):
    with pytest.warns(UserWarning, match='method parameter is ignored'):
        res = get_genomic_feature_lengths(gtf_path, 'transcript', 'mean')
    assert res == truth


@pytest.mark.parametrize('method,gene_len', [('mean', 175), ('median', 175), ('max', 200), ('min', 150),
                                             ('geometric_mean', 173.2050807568876), ('merged_exons', 250)])
def test_get_genomic_feature_lengths_gff3_multiparent_gene(method, gene_len):
    # #152: GFF3 SO-term prefixes (gene:/transcript:) are stripped so IDs match the user's count table.
    # merged_exons attributes the shared exon to the FIRST parent's gene only, so the union is 100+50+100=250.
    res = get_genomic_feature_lengths(MULTIPARENT_GFF3, 'gene', method)
    assert res.keys() == {'ENSG1'}
    assert np.isclose(res['ENSG1'], gene_len)


def test_get_genomic_feature_lengths_gff3_multiparent_transcript():
    # The comma-separated Parent=a,b shared exon (len 100) is added to BOTH transcripts (prefixes stripped, #152).
    with pytest.warns(UserWarning, match='method parameter is ignored'):
        res = get_genomic_feature_lengths(MULTIPARENT_GFF3, 'transcript', 'mean')
    assert res == {'ENST1': 150, 'ENST2': 200}


def test_get_genomic_feature_lengths_malformed_skips_and_warns(tmp_path):
    # #152: the bad-column line is skipped WITH a warning now (previously silent), matching the unified behaviour.
    pth = _malformed_gtf(tmp_path)
    with pytest.warns(UserWarning, match='malformed line'):
        assert get_genomic_feature_lengths(pth, 'transcript', 'mean') == {'T1': 101}
    with pytest.warns(UserWarning, match='malformed line'):
        assert get_genomic_feature_lengths(pth, 'gene', 'mean') == {'G1': 101}


# ------------------------------------------------------------------ map_gene_to_attr (extra coverage + known bugs)
@pytest.mark.parametrize('use_name,use_version,split_ids,truth', [
    (False, False, False, {'ENSG00000000001': 'protein_coding', 'ENSG00000000002': 'lincRNA'}),
    (False, False, True, {'ENSG00000000001': 'protein_coding', 'ENSG00000000002': 'lincRNA'}),
    (False, True, False, {'ENSG00000000001.3': 'protein_coding', 'ENSG00000000002.1': 'lincRNA'}),
    (False, True, True, {'ENSG00000000001.3': 'protein_coding', 'ENSG00000000002.1': 'lincRNA'}),
    (True, False, False, {'GENEA': 'protein_coding', 'GENEB': 'lincRNA'}),
    (True, True, True, {'GENEA': 'protein_coding', 'GENEB': 'lincRNA'}),
])
def test_map_gene_to_attr_ensembl_gene(use_name, use_version, split_ids, truth):
    assert map_gene_to_attr(ENSEMBL_GTF, 'gene_biotype', 'gene', use_name, use_version, split_ids) == truth


@pytest.mark.parametrize('use_name,use_version,split_ids,truth', [
    (False, False, False, {'ENST00000000001': 'protein_coding', 'ENST00000000002': 'protein_coding',
                           'ENST00000000003': 'lincRNA'}),
    (False, True, False, {'ENST00000000001.2': 'protein_coding', 'ENST00000000002.1': 'protein_coding',
                          'ENST00000000003.4': 'lincRNA'}),
    # BUG #1 (fix in #152): with use_name on a transcript, the value is gated on gene_name but READ from
    # transcript_name, so results are keyed by TRANSCRIPT name (not gene name, not transcript id).
    (True, False, False, {'GENEA-201': 'protein_coding', 'GENEA-202': 'protein_coding', 'GENEB-201': 'lincRNA'}),
    (True, True, True, {'GENEA-201': 'protein_coding', 'GENEA-202': 'protein_coding', 'GENEB-201': 'lincRNA'}),
])
def test_map_gene_to_attr_ensembl_transcript(use_name, use_version, split_ids, truth):
    assert map_gene_to_attr(ENSEMBL_GTF, 'transcript_biotype', 'transcript', use_name, use_version, split_ids) == truth


@pytest.mark.parametrize('feature_type,attribute,split_ids,truth', [
    ('gene', 'gene_type', False, {'ENSG00000000001.3': 'protein_coding', 'ENSG00000000002.1': 'lincRNA'}),
    ('gene', 'gene_type', True, {'ENSG00000000001': 'protein_coding', 'ENSG00000000002': 'lincRNA'}),
    ('transcript', 'transcript_type', False, {'ENST00000000001.2': 'protein_coding',
                                              'ENST00000000002.1': 'protein_coding', 'ENST00000000003.4': 'lincRNA'}),
    ('transcript', 'transcript_type', True, {'ENST00000000001': 'protein_coding', 'ENST00000000002': 'protein_coding',
                                             'ENST00000000003': 'lincRNA'}),
])
def test_map_gene_to_attr_gencode_no_version(feature_type, attribute, split_ids, truth):
    assert map_gene_to_attr(GENCODE_GTF, attribute, feature_type, False, False, split_ids) == truth


@pytest.mark.parametrize('feature_type,attribute,truth', [
    ('gene', 'gene_type', {'ENSG00000000001.3': 'protein_coding', 'ENSG00000000002.1': 'lincRNA'}),
    ('transcript', 'transcript_type', {'ENST00000000001.2': 'protein_coding', 'ENST00000000002.1': 'protein_coding',
                                       'ENST00000000003.4': 'lincRNA'}),
])
def test_map_gene_to_attr_use_version_tolerates_missing(feature_type, attribute, truth):
    # BUG #2 FIXED (#152): use_version no longer raises KeyError on a file without *_version attrs; the version
    # suffix is simply omitted, so the result matches the use_version=False case (inline-versioned GENCODE IDs).
    assert map_gene_to_attr(GENCODE_GTF, attribute, feature_type, False, True, False) == truth


def _transcript_missing_name_gtf(tmp_path):
    """One transcript has gene_name but NO transcript_name; the other has both."""
    pth = tmp_path / 'missing_tx_name.gtf'
    pth.write_text(
        '1\tsrc\ttranscript\t100\t200\t.\t+\t.\t'
        'gene_id "G1"; transcript_id "T1"; gene_name "GENEA"; transcript_name "GENEA-201"; gene_biotype "pc";\n'
        '1\tsrc\ttranscript\t300\t400\t.\t+\t.\t'
        'gene_id "G1"; transcript_id "T2"; gene_name "GENEA"; gene_biotype "pc";\n',  # no transcript_name
        encoding='utf-8', newline='\n')
    return pth


def test_map_gene_to_attr_transcript_use_name_gates_on_transcript_name(tmp_path):
    # BUG #1 FIXED (#152): with a transcript + use_name, the gate and the read are now the SAME attribute
    # (transcript_name). The transcript lacking transcript_name is dropped rather than producing a None key.
    res = map_gene_to_attr(_transcript_missing_name_gtf(tmp_path), 'gene_biotype', 'transcript', True, False, False)
    assert res == {'GENEA-201': 'pc'}
    assert None not in res


# ------------------------------------------------------------------ anchored-key regex (PR-B, #151, item #3)
# The Polars rewrite extracts attributes with an anchored regex so a key is never matched inside a longer
# compound key (e.g. gene_id must not be read from havana_gene_id). These are the tripwire for that: a naive
# unanchored pattern would capture the WRONG value below.
def _compound_key_gtf(tmp_path):
    pth = tmp_path / 'compound.gtf'
    pth.write_text(
        '1\thavana\tgene\t100\t200\t.\t+\t.\t'
        'havana_gene_id "WRONGGENE"; gene_id "RIGHTGENE"; gene_name "GENEA";\n'
        '1\thavana\ttranscript\t100\t200\t.\t+\t.\t'
        'havana_gene_id "WRONGGENE"; gene_id "RIGHTGENE"; '
        'havana_transcript_id "WRONGTX"; transcript_id "RIGHTTX"; gene_name "GENEA";\n',
        encoding='utf-8', newline='\n')
    return pth


def test_map_transcripts_to_genes_anchored_key_ignores_compound_key(tmp_path):
    # gene_id / transcript_id must come from the real keys, not havana_gene_id / havana_transcript_id.
    assert map_transcripts_to_genes(_compound_key_gtf(tmp_path), use_version=False, split_ids=False) == {
        'RIGHTTX': 'RIGHTGENE'}


def test_map_gene_to_attr_anchored_key_ignores_compound_key(tmp_path):
    # Same guarantee via the map_gene_to_attr extraction path: the gene id key is RIGHTGENE, not WRONGGENE.
    assert map_gene_to_attr(_compound_key_gtf(tmp_path), 'gene_name', 'gene', False, False, False) == {
        'RIGHTGENE': 'GENEA'}


# ====================================================================================================================
# #152 (PR-C): GFF3 parity + robustness (formats, gzip, content-sniff, aliases). New behaviour, not characterization.
# ====================================================================================================================

# ------------------------------------------------------------------ GFF3 parity for map_transcripts_to_genes
@pytest.mark.parametrize('use_version,split_ids', [(False, False), (True, False), (False, True), (True, True)])
def test_map_transcripts_to_genes_gff3_multiparent(use_version, split_ids):
    # GFF3 now accepted: ID->transcript, Parent->gene, SO-term prefixes stripped. GFF3 has no *_version, so
    # use_version is inert; ENST1/ENST2 have no internal dot, so split_ids is a no-op here.
    assert map_transcripts_to_genes(MULTIPARENT_GFF3, use_name=False, use_version=use_version,
                                    split_ids=split_ids) == {'ENST1': 'ENSG1', 'ENST2': 'ENSG1'}


def test_map_transcripts_to_genes_gff3_wormbase():
    res = map_transcripts_to_genes(WORMBASE_GFF3, use_version=False, split_ids=False)
    assert res['B0348.5a'] == 'WBGene00015153'
    assert res['B0348.6c.1'] == 'WBGene00002061'
    assert len(res) == 8  # 8 mRNA rows, all with distinct IDs


def test_map_transcripts_to_genes_gff3_split_ids_collapses_dotted():
    # split_ids strips at the first dot, collapsing the WormBase isoform ids onto their base (first-wins).
    res = map_transcripts_to_genes(WORMBASE_GFF3, use_version=False, split_ids=True)
    assert res == {'B0348': 'WBGene00015153'}


# ------------------------------------------------------------------ GFF3 parity for map_gene_to_attr
def test_map_gene_to_attr_gff3_gene():
    # biotype= lives on the gene row; feature_type='gene' keys by the (prefix-stripped) gene ID.
    assert map_gene_to_attr(MULTIPARENT_GFF3, 'biotype', 'gene', False, False, False) == {'ENSG1': 'protein_coding'}


def test_map_gene_to_attr_gff3_transcript():
    # biotype= is also on the mRNA rows; feature_type='transcript' keys by the (prefix-stripped) transcript ID.
    assert map_gene_to_attr(MULTIPARENT_GFF3, 'biotype', 'transcript', False, False, False) == {
        'ENST1': 'protein_coding', 'ENST2': 'protein_coding'}


def test_map_gene_to_attr_gff3_use_name():
    # Name= on the gene row is used as the key when use_name=True.
    assert map_gene_to_attr(MULTIPARENT_GFF3, 'biotype', 'gene', True, False, False) == {'GENEA': 'protein_coding'}


# ------------------------------------------------------------------ formats: .gff extension, gzip, content-sniff
def test_map_transcripts_to_genes_accepts_gff_extension(tmp_path):
    pth = tmp_path / 'ann.gff'
    pth.write_text(Path(MULTIPARENT_GFF3).read_text(encoding='utf-8'), encoding='utf-8', newline='\n')
    assert map_transcripts_to_genes(pth, use_version=False, split_ids=False) == {'ENST1': 'ENSG1', 'ENST2': 'ENSG1'}


def test_map_transcripts_to_genes_reads_gzip(tmp_path):
    pth = tmp_path / 'ensembl.gtf.gz'
    with gzip.open(pth, 'wt', encoding='utf-8', newline='\n') as fh:
        fh.write(Path(ENSEMBL_GTF).read_text(encoding='utf-8'))
    assert map_transcripts_to_genes(pth, use_version=False, split_ids=False) == {
        'ENST00000000001': 'ENSG00000000001', 'ENST00000000002': 'ENSG00000000001',
        'ENST00000000003': 'ENSG00000000002'}


def test_get_genomic_feature_lengths_reads_gzip(tmp_path):
    pth = tmp_path / 'ensembl.gtf.gz'
    with gzip.open(pth, 'wt', encoding='utf-8', newline='\n') as fh:
        fh.write(Path(ENSEMBL_GTF).read_text(encoding='utf-8'))
    assert get_genomic_feature_lengths(pth, 'gene', 'mean') == {'ENSG00000000001': 150, 'ENSG00000000002': 100}


def test_content_sniff_overrides_extension(tmp_path):
    # GFF3 content in a file misnamed '.gtf' must still parse as GFF3 (content-sniff wins over the extension hint).
    pth = tmp_path / 'mislabeled.gtf'
    pth.write_text(Path(MULTIPARENT_GFF3).read_text(encoding='utf-8'), encoding='utf-8', newline='\n')
    assert map_transcripts_to_genes(pth, use_version=False, split_ids=False) == {'ENST1': 'ENSG1', 'ENST2': 'ENSG1'}


def test_no_malformed_warning_for_tabless_lines(tmp_path, recwarn):
    # A trailing GFF3 ##FASTA block (and blank lines) have no tabs; they must be skipped SILENTLY, not counted
    # as malformed (otherwise a real genome would emit a warning about thousands of "malformed" lines).
    pth = tmp_path / 'with_fasta.gff3'
    pth.write_text(
        '##gff-version 3\n'
        '1\ttest\tgene\t100\t499\t.\t+\t.\tID=gene:ENSG1;biotype=protein_coding;Name=GENEA\n'
        '1\ttest\tmRNA\t100\t499\t.\t+\t.\tID=transcript:ENST1;Parent=gene:ENSG1\n'
        '\n'
        '##FASTA\n'
        '>1\n'
        'ACGTACGTACGTACGT\n'
        'ACGTACGTACGTACGT\n',
        encoding='utf-8', newline='\n')
    assert map_transcripts_to_genes(pth, use_version=False, split_ids=False) == {'ENST1': 'ENSG1'}
    assert not any('malformed' in str(w.message) for w in recwarn)


# ------------------------------------------------------------------ biotype attribute auto-resolution
def test_map_gene_to_attr_biotype_autoresolves_family():
    # Requesting the Ensembl name 'gene_biotype' on a GENCODE file (which uses 'gene_type') auto-resolves across
    # the biotype family, so non-programmer defaults keep working regardless of the annotation's convention.
    assert map_gene_to_attr(GENCODE_GTF, 'gene_biotype', 'gene', False, False, True) == {
        'ENSG00000000001': 'protein_coding', 'ENSG00000000002': 'lincRNA'}


def test_map_gene_to_attr_non_biotype_attribute_not_resolved(tmp_path):
    # Auto-resolution is limited to the biotype family; a non-family attribute stays an exact lookup.
    assert map_gene_to_attr(_compound_key_gtf(tmp_path), 'gene_name', 'gene', False, False, False) == {
        'RIGHTGENE': 'GENEA'}


def _gene_type_with_transcript_biotype_gtf(tmp_path):
    """Gene-level biotype lives in gene_type ('lincRNA'); a transcript carries a DIFFERENT transcript_biotype."""
    pth = tmp_path / 'mixed_biotype.gtf'
    pth.write_text(
        '1\tsrc\tgene\t100\t400\t.\t+\t.\tgene_id "G1"; gene_type "lincRNA";\n'
        '1\tsrc\ttranscript\t100\t400\t.\t+\t.\tgene_id "G1"; transcript_id "T1"; gene_type "lincRNA"; '
        'transcript_biotype "protein_coding";\n',
        encoding='utf-8', newline='\n')
    return pth


def test_map_gene_to_attr_biotype_resolution_stays_within_feature_level(tmp_path):
    # Review regression guard: a gene-level 'gene_biotype' request must NOT fall back to a transcript-level
    # attribute (transcript_biotype). It resolves to the same-level gene_type -> the gene's true biotype 'lincRNA',
    # never the transcript's 'protein_coding'.
    assert map_gene_to_attr(_gene_type_with_transcript_biotype_gtf(tmp_path), 'gene_biotype', 'gene',
                            False, False, False) == {'G1': 'lincRNA'}
