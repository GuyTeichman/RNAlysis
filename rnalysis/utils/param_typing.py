import collections.abc
import ftplib
import functools
import typing
import warnings

import requests
import tenacity
import typing_extensions

from rnalysis import FROZEN_ENV
from rnalysis.utils import io, parsing

PARALLEL_BACKENDS = ('multiprocessing', 'sequential') if FROZEN_ENV else (
    'multiprocessing', 'loky', 'threading', 'sequential')
QUANTILE_INTERPOLATION_METHODS = ("nearest", "higher", "lower", "midpoint", "linear")
SUMMATION_METHODS = ('scaled_tpm', 'raw')
LIMMA_NORM = ('scale', 'quantile', 'cyclicloess')
K_CRITERIA = ('gap', 'silhouette', 'calinski_harabasz', 'davies_bouldin', 'bic')

LEGAL_GENE_LENGTH_METHODS = ('mean', 'median', 'max', 'min', 'geometric_mean', 'merged_exons')

LEGAL_FASTQ_SUFFIXES = ('.fastq', '.fastq.gz', '.fq', '.fq.gz')
LEGAL_ALIGNMENT_SUFFIXES = ('.sam', '.bam', '.cram',)
LEGAL_BOWTIE2_PRESETS = ('very-fast', 'fast', 'sensitive', 'very-sensitive')
LEGAL_BOWTIE2_MODES = ('end-to-end', 'local')
LEGAL_QUAL_SCORE_TYPES = ('phred33', 'phred64', 'solexa-quals', 'int-quals')

GRAPHVIZ_FORMATS = ('pdf', 'png', 'svg', 'none')

BIOTYPES = ('protein_coding', 'pseudogene', 'lincRNA', 'miRNA', 'ncRNA', 'piRNA', 'rRNA', 'snoRNA', 'snRNA', 'tRNA')
BIOTYPE_ATTRIBUTE_NAMES = ('biotype', 'gene_biotype', 'transcript_biotype', 'gene_type', 'transcript_type')
# Feature-type keywords (column 3) that denote a transcript, across GTF and GFF3 conventions. Used to identify
# transcript rows when parsing annotation files (mRNA/primary_transcript are the GFF3/GTF equivalents of 'transcript').
TRANSCRIPT_FEATURE_NAMES = ('transcript', 'mRNA', 'primary_transcript')

GO_ASPECTS = ('biological_process', 'cellular_component', 'molecular_function')
GO_EVIDENCE_TYPES = ('experimental', 'phylogenetic', 'computational', 'author', 'curator', 'electronic')
GO_QUALIFIERS = ('not', 'contributes_to', 'colocalizes_with')

DEFAULT_ORGANISMS = tuple(sorted(['Caenorhabditis elegans',
                                  'Mus musculus',
                                  'Drosophila melanogaster',
                                  'Homo sapiens',
                                  'Arabidopsis thaliana',
                                  'Danio rerio',
                                  'Escherichia coli',
                                  'Saccharomyces cerevisiae',
                                  'Schizosaccharomyces pombe']))

ORTHOLOG_NON_UNIQUE_MODES = ('first', 'last', 'random', 'none')

Fraction = typing.NewType('Fraction', float)
PositiveInt = typing.NewType('PositiveInt', int)
NegativeInt = typing.NewType('NegativeInt', int)
NonNegativeInt = typing.NewType('NonNegativeInt', int)

ColumnName = typing.NewType('ColumnName', str)
ColumnNames = typing.NewType('ColumnNames', typing.Union[ColumnName, typing.Iterable[ColumnName]])
GroupedColumns = typing.NewType('GroupedColumns', typing.List[typing.Iterable[ColumnName]])

Color = typing.NewType('Color', typing.Union[str, typing.Tuple[float, float, float]])
ColorList = typing.NewType('ColorList', typing.Union[typing.List[Color], typing.Tuple[Color, ...]])
ColorMap = typing.NewType('ColorMap', str)


def type_to_supertype(this_type):
    if hasattr(this_type, '__supertype__'):
        return type_to_supertype(this_type.__supertype__)

    args = typing_extensions.get_args(this_type)
    if isinstance(args, tuple) and len(args) > 0:
        origin_map = {list: typing.List, set: typing.Set, dict: typing.Dict, collections.abc.Iterable: typing.Iterable,
                      frozenset: typing.FrozenSet, typing.Union: typing.Union,
                      typing_extensions.Literal: typing_extensions.Literal}
        components = []
        for component in args:
            components.append(type_to_supertype(component))
        origin = origin_map[typing_extensions.get_origin(this_type)]
        return origin[parsing.data_to_tuple(components)]

    return this_type


# Errors that mean "couldn't fetch or parse the remote list of legal values". These getters run at
# import time to populate Literal[...] type annotations, so an uncaught error here crashes
# `import rnalysis.filtering` for every user and every test-collection run — degrade to an empty tuple
# plus a warning instead. RequestException covers connection/HTTP/timeout failures (the underlying
# io.get_legal_* calls now pass a request timeout); ValueError covers JSON-decode errors from a
# 200-with-garbage body; Key/Type/Index errors cover otherwise-parseable-but-unexpected shapes.
# These getters are lru_cached, so a degraded empty result is intentionally pinned for the process:
# the Literal[...] annotations are evaluated exactly once at import, and the warning tells the user to
# restart RNAlysis to retry.
_LEGAL_VALUE_FETCH_ERRORS = (requests.exceptions.RequestException, tenacity.RetryError,
                             ValueError, KeyError, TypeError, IndexError)


@functools.lru_cache(maxsize=2)
def get_gene_id_types() -> typing.Tuple[str, ...]:
    try:
        gene_id_types = parsing.data_to_tuple(io.get_legal_gene_id_types()[0].keys())
    except _LEGAL_VALUE_FETCH_ERRORS:
        gene_id_types = tuple()
        warnings.warn('Failed to retrieve gene ID mapping data from UniProtKB. '
                      'Some features may not work as intended. '
                      'To fix this issue, make sure your computer has internet connection, '
                      'and restart RNAlysis. ')
    return gene_id_types


@functools.lru_cache(maxsize=2)
def get_panther_taxons() -> typing.Tuple[str, ...]:
    try:
        taxons = io.get_legal_panther_taxons()
    except _LEGAL_VALUE_FETCH_ERRORS:
        taxons = tuple()
        warnings.warn('Failed to retrieve legal taxons from PantherDB. '
                      'Some features may not work as intended. '
                      'To fix this issue, make sure your computer has internet connection, '
                      'and restart RNAlysis. ')
    return taxons


@functools.lru_cache(maxsize=2)
def get_phylomedb_taxons() -> typing.Tuple[str, ...]:
    try:
        taxons = io.get_legal_phylomedb_taxons()
    except (*_LEGAL_VALUE_FETCH_ERRORS, *ftplib.all_errors):
        taxons = tuple()
        warnings.warn('Failed to retrieve legal taxons from PhylomeDB. '
                      'Some features may not work as intended. '
                      'To fix this issue, make sure your computer has internet connection, '
                      'and restart RNAlysis. ')
    return taxons


@functools.lru_cache(maxsize=2)
def get_ensembl_taxons() -> typing.Tuple[str, ...]:
    try:
        taxons = parsing.data_to_tuple(io.get_legal_ensembl_taxons())
    except _LEGAL_VALUE_FETCH_ERRORS:
        taxons = tuple()
        warnings.warn('Failed to retrieve legal taxons from Ensembl. '
                      'Some features may not work as intended. '
                      'To fix this issue, make sure your computer has internet connection, '
                      'and restart RNAlysis. ')
    return taxons
