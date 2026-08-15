import collections.abc
import functools
import json
import typing
import warnings
from pathlib import Path

import typing_extensions

from rnalysis import FROZEN_ENV
from rnalysis.utils import parsing

PARALLEL_BACKENDS = ('multiprocessing', 'sequential') if FROZEN_ENV else (
    'multiprocessing', 'loky', 'threading', 'sequential')
QUANTILE_INTERPOLATION_METHODS = ("nearest", "higher", "lower", "midpoint", "linear")
SUMMATION_METHODS = ('scaled_tpm', 'raw')
LIMMA_NORM = ('scale', 'quantile', 'cyclicloess')
K_CRITERIA = ('gap', 'silhouette', 'calinski_harabasz', 'davies_bouldin', 'bic')
# The transforms applied to a count matrix before PCA/clustering ('power_transform'). This parameter used to be
# a bool, and the legacy True ('box-cox') / False ('none') are still accepted by the API -- see
# generic.parse_power_transform -- but the GUI only ever offers the named transforms.
POWER_TRANSFORMS = ('box-cox', 'log', 'none')

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
# Attribute names suggested in the GUI for filtering/annotating a table from a GTF/GFF file. Free text is also
# allowed (used as Union[Literal[GTF_ATTRIBUTE_NAMES], str]). 'chromosome'/'source'/'strand' are reserved names that
# read the fixed GTF/GFF columns; everything else is looked up as a column-9 key=value attribute.
GTF_ATTRIBUTE_NAMES = ('gene_biotype', 'transcript_biotype', 'biotype', 'gene_name', 'gene_id', 'transcript_id',
                       'chromosome', 'source', 'strand')

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


# The legal values below (UniProtKB gene ID types, PantherDB/PhylomeDB/Ensembl taxons) come from
# remote services, but they are needed to build Literal[...] type annotations on public
# Filter/FeatureSet methods — which are evaluated when those class bodies run, i.e. during
# `import rnalysis.filtering`. Fetching them live would therefore put four web requests on the import
# path: slow for everyone, and on an offline machine (a cluster, a locked-down lab PC) it would leave
# every one of those dropdowns empty. Instead they are read from a snapshot packaged with RNAlysis and
# regenerated at release time by packaging/generate_api_vocabularies.py, which makes the legal values
# a property of the RNAlysis *version* — identical on every machine, reproducible, and available
# offline. Live operations (actual ID mapping, ortholog lookups, enrichment) still query the services
# at call time through rnalysis.utils.io; only this vocabulary is snapshotted.
API_VOCABULARIES_PATH = Path(__file__).parent.parent.joinpath('data_files/api_vocabularies.json')

# Errors that mean "the packaged snapshot is missing or unreadable": OSError for a missing/unreadable
# file, ValueError for a JSON-decode failure, Key/Type/Attribute/Index errors for a file that parses
# but has an unexpected shape. Because these getters run inside class bodies, an uncaught error here
# would crash the import for every user and every test-collection run — degrade to an empty tuple plus
# a warning instead. Unlike the old live-fetch version, a failure now means a broken installation
# rather than a network problem, so the warning says so; it is still pinned by lru_cache for the
# process, since the annotations are evaluated exactly once at import.
_SNAPSHOT_READ_ERRORS = (OSError, ValueError, KeyError, TypeError, AttributeError, IndexError)


@functools.lru_cache(maxsize=1)
def _load_api_vocabularies() -> typing.Dict[str, dict]:
    """Read the packaged vocabulary snapshot. Raises one of _SNAPSHOT_READ_ERRORS if it is unusable."""
    with open(API_VOCABULARIES_PATH, encoding='utf-8') as snapshot_file:
        vocabularies = json.load(snapshot_file)['vocabularies']
    if not isinstance(vocabularies, dict):
        raise TypeError(f"expected a dict of vocabularies, got {type(vocabularies)}")
    return vocabularies


def _get_snapshot_vocabulary(key: str, description: str) -> typing.Tuple[str, ...]:
    try:
        values = _load_api_vocabularies()[key]['values']
        # Literal[...] can only be built from hashable values, so anything but strings here would
        # crash the import instead of emptying a single dropdown.
        if not isinstance(values, (list, tuple)) or not all(isinstance(value, str) for value in values):
            raise TypeError(f"expected a list of strings under vocabulary '{key}'")
        return parsing.data_to_tuple(values)
    except _SNAPSHOT_READ_ERRORS:
        warnings.warn(f'Failed to load the list of {description} from the packaged data file '
                      f'"{API_VOCABULARIES_PATH}". '
                      f'Some features may not work as intended, and some drop-down menus will be empty. '
                      f'This usually means your RNAlysis installation is incomplete or corrupted - '
                      f'to fix this issue, re-install RNAlysis. ')
        return tuple()


@functools.lru_cache(maxsize=2)
def get_gene_id_types() -> typing.Tuple[str, ...]:
    return _get_snapshot_vocabulary('gene_id_types', 'gene ID types from UniProtKB')


@functools.lru_cache(maxsize=2)
def get_panther_taxons() -> typing.Tuple[str, ...]:
    return _get_snapshot_vocabulary('panther_taxons', 'legal taxons from PantherDB')


@functools.lru_cache(maxsize=2)
def get_phylomedb_taxons() -> typing.Tuple[str, ...]:
    return _get_snapshot_vocabulary('phylomedb_taxons', 'legal taxons from PhylomeDB')


@functools.lru_cache(maxsize=2)
def get_ensembl_taxons() -> typing.Tuple[str, ...]:
    return _get_snapshot_vocabulary('ensembl_taxons', 'legal taxons from Ensembl')
