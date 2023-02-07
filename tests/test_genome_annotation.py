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
    assert res == truth
