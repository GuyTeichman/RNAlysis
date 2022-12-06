import collections.abc
import typing
import warnings

import requests
import typing_extensions

from rnalysis.utils import io, parsing

BIOTYPES = ('protein_coding', 'pseudogene', 'lincRNA', 'miRNA', 'ncRNA', 'piRNA', 'rRNA', 'snoRNA', 'snRNA', 'tRNA')
BIOTYPE_ATTRIBUTE_NAMES = ('biotype', 'gene_biotype', 'transcript_biotype', 'gene_type', 'transcript_type')
GO_ASPECTS = ('biological_process', 'molecular function', 'cellular component')
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

ColumnName = typing.NewType('ColumnName', str)
ColumnNames = typing.NewType('ColumnNames', typing.Union[ColumnName, typing.Iterable[ColumnName]])
GroupedColumns = typing.NewType('GroupedColumns', typing.List[typing.Iterable[ColumnName]])

Color = typing.NewType('Color', typing.Union[str, typing.Tuple[float, float, float]])
ColorList = typing.NewType('ColorList', typing.Union[typing.List[Color], typing.Tuple[Color, ...]])


def type_to_supertype(this_type):
    if hasattr(this_type, '__supertype__'):
        return type_to_supertype(this_type.__supertype__)

    args = typing_extensions.get_args(this_type)
    if isinstance(args, tuple) and len(args) > 0:
        origin_map = {list: typing.List, set: typing.Set, dict: typing.Dict, collections.abc.Iterable: typing.Iterable,
                      frozenset: typing.FrozenSet, typing.Union: typing.Union, typing.Literal: typing.Literal}
        components = []
        for component in args:
            components.append(type_to_supertype(component))
        origin = origin_map[typing_extensions.get_origin(this_type)]
        return origin[parsing.data_to_tuple(components)]

    return this_type


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
