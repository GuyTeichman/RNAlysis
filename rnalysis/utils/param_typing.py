import collections.abc
import ftplib
import functools
import typing
import warnings

import requests
import tenacity
import typing_extensions
from polars._typing import RollingInterpolationMethod

from rnalysis import FROZEN_ENV
from rnalysis.utils import io, parsing

PARALLEL_BACKENDS = ('multiprocessing', 'sequential') if FROZEN_ENV else (
'multiprocessing', 'loky', 'threading', 'sequential')
QUANTILE_INTERPOLATION_METHODS = RollingInterpolationMethod
SUMMATION_METHODS = ('scaled_tpm', 'raw')

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

GO_ASPECTS = ('biological_process', 'cellular_component', 'molecular_function')
GO_EVIDENCE_TYPES = ('experimental', 'phylogenetic', 'computational', 'author', 'curator', 'electronic')
GO_QUALIFIERS = ('not', 'contributes_to', 'colocalizes_with')

DEFAULT_ORGANISMS = tuple(sorted(['Caenorhabditis elegans',
                                  'Mus musculus',
                                  'Drosophila melanogaster',
                                  'Homo sapiens',
                                  'Arabodopsis thaliana',
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


@functools.lru_cache(maxsize=2)
def get_gene_id_types() -> typing.Tuple[str, ...]:
    try:
        gene_id_types = parsing.data_to_tuple(io.get_legal_gene_id_types()[0].keys())
    except requests.exceptions.ConnectionError:
        gene_id_types = tuple()
        warnings.warn('Failed to retreive gene ID mapping data from UniProtKB. '
                      'Some features may not work as intended. '
                      'To fix this issue, make sure your computer has internet connection, '
                      'and restart RNAlysis. ')
    return gene_id_types


@functools.lru_cache(maxsize=2)
def get_panther_taxons() -> typing.Tuple[str, ...]:
    try:
        taxons = io.get_legal_panther_taxons()
    except requests.exceptions.ConnectionError:
        taxons = tuple()
        warnings.warn('Failed to retreive legal taxons from PantherDB. '
                      'Some features may not work as intended. '
                      'To fix this issue, make sure your computer has internet connection, '
                      'and restart RNAlysis. ')
    return taxons


@functools.lru_cache(maxsize=2)
def get_phylomedb_taxons() -> typing.Tuple[str, ...]:
    try:
        taxons = io.get_legal_phylomedb_taxons()
    except ftplib.all_errors:
        taxons = tuple()
        warnings.warn('Failed to retreive legal taxons from PhylomeDB. '
                      'Some features may not work as intended. '
                      'To fix this issue, make sure your computer has internet connection, '
                      'and restart RNAlysis. ')
    return taxons


@functools.lru_cache(maxsize=2)
def get_ensembl_taxons() -> typing.Tuple[str, ...]:
    try:
        taxons = parsing.data_to_tuple(io.get_legal_ensembl_taxons())
    except (requests.exceptions.ConnectionError, tenacity.RetryError):
        taxons = tuple()
        warnings.warn('Failed to retreive legal taxons from Ensembl. '
                      'Some features may not work as intended. '
                      'To fix this issue, make sure your computer has internet connection, '
                      'and restart RNAlysis. ')
    return taxons
