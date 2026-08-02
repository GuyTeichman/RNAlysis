"""
Parsing and length/attribute extraction for genome-annotation files (GTF and GFF3).

This module is the single home for all GTF/GFF parsing in RNAlysis. The three public helpers
(:func:`map_transcripts_to_genes`, :func:`get_genomic_feature_lengths`, :func:`map_gene_to_attr`) are
consumed from ``filtering`` and ``fastq``; their signatures and dict return-types are a stable contract.

Internally the pipeline is Polars-native: the file is scanned with :func:`pl.scan_csv`, one row per line,
and attribute values are pulled out column-wise with anchored regular expressions (see
:func:`_extract_gtf_attribute` / :func:`_extract_gff3_attribute`). The per-string dict parsers
:func:`parse_gtf_attributes` / :func:`parse_gff3_attributes` are kept here (and re-exported from
``utils.parsing`` for backwards compatibility) but are no longer on the hot path.
"""
import re
import warnings
from pathlib import Path
from typing import Literal, Union

import numpy as np
import polars as pl
from scipy.stats.mstats import gmean

from rnalysis.utils import generic, validation
from rnalysis.utils.param_typing import BIOTYPE_ATTRIBUTE_NAMES, LEGAL_GENE_LENGTH_METHODS, TRANSCRIPT_FEATURE_NAMES

# Read whole annotation lines as a single Polars column: pick a field separator that never occurs in a
# GTF/GFF file, so scan_csv yields one row per line and we split on the real '\t' ourselves. This lets us
# enforce the 9-column rule ourselves: malformed (wrong-column-count) lines are dropped, with a single warning,
# uniformly across all three public functions (see :func:`_load_annotation_fields`).
_LINE_SEP = '\x01'
_N_COLUMNS = 9
_FEATURE_IDX, _START_IDX, _END_IDX, _ATTR_IDX = 2, 3, 4, 8
# Ensembl-style GFF3 gives IDs an SO-term prefix (``ID=gene:ENSG1``, ``Parent=transcript:ENST1``). Strip it so the
# emitted gene/transcript IDs match a user's count table (WormBase-style bare IDs have no colon and are untouched).
_SO_PREFIX_RE = r'^[A-Za-z_]+:'


def parse_gtf_attributes(attr_str: str):
    attributes_dict = {}
    for this_attr in attr_str.split('; '):
        this_attr = this_attr.strip()
        key_end = this_attr.find(' ')
        if key_end == -1:
            continue
        key = this_attr[:key_end]
        val_start = this_attr.find('"', key_end)
        val_end = this_attr.find('"', val_start + 1)
        val = this_attr[val_start + 1:val_end]
        attributes_dict[key] = val
    return attributes_dict


def parse_gff3_attributes(attr_str: str):
    attributes_dict = {}
    for this_attr in attr_str.rstrip().split(';'):
        this_attr = this_attr.strip()
        if len(this_attr) == 0:
            continue
        key, val = this_attr.split('=')
        val = val.split(',')
        if len(val) == 1:
            val = val[0]
        attributes_dict[key] = val
    return attributes_dict


def _extract_gtf_attribute(attr_col: pl.Expr, key: str) -> pl.Expr:
    """Column-wise extraction of a GTF attribute value (``key "value";``).

    The key is anchored on a start-of-string / ``;`` boundary so that, e.g., extracting ``gene_id`` does
    **not** match the tail of a compound key like ``havana_gene_id`` (correctness item #3). Returns null
    where the key is absent, mirroring ``key not in parse_gtf_attributes(...)``. If a key occurs more than
    once in a line (non-standard), the first occurrence wins; standard GTF attribute keys are unique.
    """
    pattern = rf'(?:^|;)\s*{re.escape(key)}\s+"([^"]*)"'
    return attr_col.str.extract(pattern, 1)


def _extract_gff3_attribute(attr_col: pl.Expr, key: str) -> pl.Expr:
    """Column-wise extraction of a GFF3 attribute value (``key=value``), anchored like the GTF variant."""
    pattern = rf'(?:^|;)\s*{re.escape(key)}=([^;]*)'
    return attr_col.str.extract(pattern, 1)


def _strip_so_prefix(expr: pl.Expr) -> pl.Expr:
    """Strip a leading GFF3 SO-term ``type:`` prefix from an ID expression (``transcript:ENST1`` -> ``ENST1``)."""
    return expr.str.replace(_SO_PREFIX_RE, '')


def _scan_annotation_lines(path: Union[str, Path]) -> pl.LazyFrame:
    """Scan an annotation file into one row per (non-comment) line, split into its tab-separated fields.

    Yields columns ``_fields`` (List[str]) and ``_n`` (field count) so callers can enforce the 9-column
    rule themselves. ``quote_char=None`` keeps the literal ``"`` characters of GTF values intact. ``scan_csv``
    transparently reads gzip-compressed (``.gz``) files.
    """
    return (pl.scan_csv(path, separator=_LINE_SEP, has_header=False, comment_prefix='#',
                        quote_char=None, schema={'raw': pl.String})
            .with_columns(pl.col('raw').str.split('\t').alias('_fields'))
            .with_columns(pl.col('_fields').list.len().alias('_n')))


def _load_annotation_fields(path: Union[str, Path]) -> pl.DataFrame:
    """Scan the file to per-line tab-fields and drop malformed (``_n != 9``) lines, warning once with the count.

    This is the single choke-point that unifies malformed-line handling across all three public functions
    (previously they diverged: some raised, some silently skipped). Returns a DataFrame with the ``raw``,
    ``_fields`` (List[str]) and ``_n`` columns, restricted to well-formed 9-column rows.
    """
    df = _scan_annotation_lines(path).collect()
    # A record line has 8 tabs (9 fields); count as malformed only lines that DO contain tabs but the wrong number
    # of them. Tab-less lines (blank lines, or a GFF3 ``##FASTA`` block's sequence lines) are silently skipped so
    # they don't inflate the warning to thousands of "malformed" lines on a real genome annotation.
    n_malformed = df.filter((pl.col('_n') != _N_COLUMNS) & (pl.col('_n') > 1)).height
    if n_malformed:
        warnings.warn(f'Skipped {n_malformed} malformed line(s) in the annotation file '
                      f'(expected {_N_COLUMNS} tab-separated columns).')
    return df.filter(pl.col('_n') == _N_COLUMNS)


def _read_annotation_records(path: Union[str, Path], file_type: Literal['gtf', 'gff3']) -> pl.DataFrame:
    """Materialize the feature records needed for length aggregation.

    Malformed lines are dropped (with a warning) by :func:`_load_annotation_fields`. Columns: ``feature``,
    ``start``, ``end``, ``_transcript`` (this feature's own id, for transcript/mRNA rows), ``_gene`` (its parent
    gene, for transcript/mRNA rows) and ``_parents`` (List[str] of parent transcript ids, for exon rows). For GFF3,
    SO-term ``type:`` prefixes are stripped from IDs so they match a user's count table.
    """
    base = _load_annotation_fields(path).lazy().select(
        pl.col('_fields').list.get(_FEATURE_IDX).alias('feature'),
        pl.col('_fields').list.get(_START_IDX).cast(pl.Int64).alias('start'),
        pl.col('_fields').list.get(_END_IDX).cast(pl.Int64).alias('end'),
        pl.col('_fields').list.get(_ATTR_IDX).alias('_attr'),
    )
    if file_type == 'gtf':
        records = base.with_columns(
            _extract_gtf_attribute(pl.col('_attr'), 'transcript_id').alias('_transcript'),
            _extract_gtf_attribute(pl.col('_attr'), 'gene_id').alias('_gene'),
        ).with_columns(pl.concat_list(pl.col('_transcript')).alias('_parents'))
    else:  # gff3: a transcript's id is its own ID, its gene is Parent; an exon's Parent may be comma-separated
        records = base.with_columns(
            _strip_so_prefix(_extract_gff3_attribute(pl.col('_attr'), 'ID')).alias('_transcript'),
            _strip_so_prefix(_extract_gff3_attribute(pl.col('_attr'), 'Parent')).alias('_gene'),
        ).with_columns(_extract_gff3_attribute(pl.col('_attr'), 'Parent').str.split(',')
                       .list.eval(_strip_so_prefix(pl.element())).alias('_parents'))
    return records.collect()


# Biotype attribute names by feature level. Fallback stays WITHIN a level so a gene-level request can never
# silently resolve to a transcript-level biotype (and vice versa); ``biotype`` is level-agnostic (GFF3 convention).
_GENE_BIOTYPE_FAMILY = ('gene_biotype', 'gene_type', 'biotype')
_TRANSCRIPT_BIOTYPE_FAMILY = ('transcript_biotype', 'transcript_type', 'biotype')


def _resolve_attribute_name(fields: pl.DataFrame, file_type: Literal['gtf', 'gff3'], attribute: str,
                            feature_type: str) -> str:
    """Resolve a biotype attribute across its naming-convention family (#152 hybrid alias resolution).

    If ``attribute`` is a member of :data:`BIOTYPE_ATTRIBUTE_NAMES` (e.g. the default ``gene_biotype``) but is
    absent from the file, fall back to the first *same-level* family member that is present (e.g. GENCODE's
    ``gene_type`` for a gene request), so non-programmer defaults keep working regardless of the annotation's
    convention. The fallback never crosses feature levels, so a gene request cannot pick up a transcript-level
    biotype. An explicitly requested attribute that is present, and any non-family attribute, is returned unchanged.
    """
    if attribute not in BIOTYPE_ATTRIBUTE_NAMES:
        return attribute
    family = _GENE_BIOTYPE_FAMILY if feature_type == 'gene' else _TRANSCRIPT_BIOTYPE_FAMILY
    candidates = [attribute] + [name for name in family if name != attribute]
    extract = _extract_gtf_attribute if file_type == 'gtf' else _extract_gff3_attribute
    attr_col = pl.col('_fields').list.get(_ATTR_IDX)
    present = fields.lazy().select(
        [extract(attr_col, name).is_not_null().any().alias(name) for name in candidates]).collect()
    for name in candidates:
        if present[name][0]:
            return name
    return attribute


def _sum_exon_lengths_per_parent(exons: pl.DataFrame) -> dict:
    """Sum exon lengths onto each parent transcript. A shared exon (multiple parents) is added to each,
    preserving the GFF3 comma-separated ``Parent=a,b`` semantics (correctness item #4)."""
    # exon _parents lists always hold >=1 parent, so empty_as_null is immaterial here; set it explicitly to
    # silence the Polars 1.x deprecation warning and pin behaviour across the 2.0 default change.
    grouped = exons.explode('_parents', empty_as_null=False).group_by('_parents').agg(pl.col('_length').sum())
    return dict(zip(grouped['_parents'].to_list(), grouped['_length'].to_list()))


def map_transcripts_to_genes(gtf_path: Union[str, Path], use_name: bool = False, use_version: bool = True,
                             split_ids: bool = True):
    file_type = validation.validate_genome_annotation_file(gtf_path)

    fields = _load_annotation_fields(gtf_path)  # drops malformed lines, warning once
    base = fields.lazy().select(
        pl.col('_fields').list.get(_FEATURE_IDX).alias('feature'),
        pl.col('_fields').list.get(_ATTR_IDX).alias('_attr'),
    )
    if file_type == 'gtf':
        transcripts = base.filter(pl.col('feature') == 'transcript').with_columns(
            _extract_gtf_attribute(pl.col('_attr'), 'transcript_id').alias('transcript_id'),
            _extract_gtf_attribute(pl.col('_attr'), 'gene_id').alias('gene_id'),
            _extract_gtf_attribute(pl.col('_attr'), 'transcript_version').alias('transcript_version'),
            _extract_gtf_attribute(pl.col('_attr'), 'gene_version').alias('gene_version'),
            _extract_gtf_attribute(pl.col('_attr'), 'gene_name').alias('gene_name'),
        )
    else:  # gff3: ID -> transcript, Parent -> gene; no *_version attrs (versions are inline in the ID)
        # NOTE: use_name on GFF3 reads the transcript row's own ``Name`` (GFF3 puts the gene's Name on the gene
        # row, not the transcript row), so it maps transcript -> transcript-name rather than -> gene-name. Callers
        # (fastq count-summation, biotype ref-table) try use_name=False first, so this fallback is rarely reached.
        transcripts = base.filter(pl.col('feature').is_in(TRANSCRIPT_FEATURE_NAMES)).with_columns(
            _strip_so_prefix(_extract_gff3_attribute(pl.col('_attr'), 'ID')).alias('transcript_id'),
            _strip_so_prefix(_extract_gff3_attribute(pl.col('_attr'), 'Parent')).alias('gene_id'),
            pl.lit(None, dtype=pl.String).alias('transcript_version'),
            pl.lit(None, dtype=pl.String).alias('gene_version'),
            _extract_gff3_attribute(pl.col('_attr'), 'Name').alias('gene_name'),
        )
    transcripts = transcripts.filter(pl.col('transcript_id').is_not_null() & pl.col('gene_id').is_not_null())

    if split_ids:
        transcripts = transcripts.with_columns(
            pl.col('transcript_id').str.split('.').list.first().alias('transcript_id'),
            pl.col('gene_id').str.split('.').list.first().alias('gene_id'),
        )
    if use_version:
        # a version suffix is appended only when BOTH transcript_version and gene_version are present
        both_versions = pl.col('transcript_version').is_not_null() & pl.col('gene_version').is_not_null()
        transcripts = transcripts.with_columns(
            pl.when(both_versions).then(pl.col('transcript_id') + '.' + pl.col('transcript_version'))
            .otherwise(pl.col('transcript_id')).alias('transcript_id'),
            pl.when(both_versions).then(pl.col('gene_id') + '.' + pl.col('gene_version'))
            .otherwise(pl.col('gene_id')).alias('gene_id'),
        )
    if use_name:
        transcripts = transcripts.filter(pl.col('gene_name').is_not_null())
        value = pl.col('gene_name')
    else:
        value = pl.col('gene_id')

    # first-wins on duplicate transcript ids (e.g. WormBase ids collapsed by split_ids)
    result = transcripts.select(pl.col('transcript_id').alias('key'), value.alias('val')).unique(
        subset=['key'], keep='first', maintain_order=True).collect()
    return dict(zip(result['key'].to_list(), result['val'].to_list()))


def get_genomic_feature_lengths(gtf_path: Union[str, Path], feature_type: Literal['gene', 'transcript'],
                                method: Literal[LEGAL_GENE_LENGTH_METHODS]):
    method_funcs = {'mean': np.mean, 'median': np.median, 'max': np.max, 'min': np.min, 'geometric_mean': gmean,
                    'merged_exons': None}
    assert feature_type in ('gene', 'transcript'), f"Illegal feature_type: '{feature_type}'"
    assert method in method_funcs, f"Illegal method: '{method}'."
    file_type = validation.validate_genome_annotation_file(gtf_path)

    records = _read_annotation_records(gtf_path, file_type)
    # map each transcript/mRNA to its parent gene (dict assignment -> last-wins on duplicate ids)
    transcript_rows = records.filter(pl.col('feature').is_in(('transcript', 'mRNA')))
    transcripts_to_genes = dict(zip(transcript_rows['_transcript'].to_list(), transcript_rows['_gene'].to_list()))
    exons = records.filter(pl.col('feature') == 'exon').with_columns(
        (pl.col('end') - pl.col('start') + 1).alias('_length'))

    if feature_type == 'transcript':
        warnings.warn("Since feature_type='transcript', the method parameter is ignored. ")
        return _sum_exon_lengths_per_parent(exons)

    if method == 'merged_exons':
        # attribute each exon to the FIRST parent transcript's gene, then union the exon intervals per gene
        gene_intervals = {}
        for parents, start, end in exons.select('_parents', 'start', 'end').iter_rows():
            gene_id = transcripts_to_genes[parents[0]]
            gene_intervals.setdefault(gene_id, []).append((start, end))
        return {gene_id: generic.sum_intervals_inclusive(intervals)
                for gene_id, intervals in gene_intervals.items()}

    # average the per-transcript lengths within each gene, per the chosen method
    transcript_lengths = _sum_exon_lengths_per_parent(exons)
    gene_lengths = {}
    for transcript, gene_id in transcripts_to_genes.items():
        gene_lengths.setdefault(gene_id, []).append(transcript_lengths[transcript])

    avg_func = method_funcs[method]
    avg_gene_lengths = {}
    for gene in gene_lengths:
        if len(gene_lengths[gene]) == 1:
            avg_gene_lengths[gene] = gene_lengths[gene][0]
        else:
            avg_gene_lengths[gene] = avg_func(gene_lengths[gene])
    return avg_gene_lengths


def map_gene_to_attr(gtf_path: Union[str, Path], attribute: str, feature_type: str, use_name: bool, use_version: bool,
                     split_ids: bool):
    file_type = validation.validate_genome_annotation_file(gtf_path)
    assert feature_type in {'gene', 'transcript'}, f"Invalid feature_type: '{feature_type}'"

    fields = _load_annotation_fields(gtf_path)  # drops malformed lines, warning once
    attribute = _resolve_attribute_name(fields, file_type, attribute, feature_type)  # #152 biotype-family resolution
    if file_type == 'gtf':
        return _map_gene_to_attr_gtf(fields, attribute, feature_type, use_name, use_version, split_ids)
    return _map_gene_to_attr_gff3(fields, attribute, feature_type, use_name, split_ids)


def _map_gene_to_attr_gtf(fields: pl.DataFrame, attribute: str, feature_type: str, use_name: bool, use_version: bool,
                          split_ids: bool) -> dict:
    # extract everything the reduction below reads, one row per line, in file order
    rows = fields.lazy().with_columns(
        pl.col('_fields').list.get(_ATTR_IDX).alias('_attr')).select(
        _extract_gtf_attribute(pl.col('_attr'), attribute).alias('attr_value'),
        _extract_gtf_attribute(pl.col('_attr'), 'gene_id').alias('gene_id'),
        _extract_gtf_attribute(pl.col('_attr'), 'transcript_id').alias('transcript_id'),
        _extract_gtf_attribute(pl.col('_attr'), 'gene_version').alias('gene_version'),
        _extract_gtf_attribute(pl.col('_attr'), 'transcript_version').alias('transcript_version'),
        _extract_gtf_attribute(pl.col('_attr'), 'gene_name').alias('gene_name'),
        _extract_gtf_attribute(pl.col('_attr'), 'transcript_name').alias('transcript_name'),
        _extract_gtf_attribute(pl.col('_attr'), 'name').alias('name'),
    ).collect()

    mapping = {}
    for row in rows.iter_rows(named=True):
        if row['attr_value'] is None:  # this line does not carry the requested attribute
            continue
        feature_name = None
        if feature_type == 'gene':
            if row['gene_id'] is None:
                raise KeyError('gene_id')
            feature_id = row['gene_id'].split('.')[0] if split_ids else row['gene_id']
            if use_version and row['gene_version'] is not None:  # #152 bug #2: tolerate a missing version suffix
                feature_id += '.' + row['gene_version']
            if use_name:
                if row['gene_name'] is not None:
                    feature_name = row['gene_name']
                elif row['name'] is not None:
                    feature_name = row['name']
                else:
                    continue
        else:
            if row['transcript_id'] is None:
                raise KeyError('transcript_id')
            feature_id = row['transcript_id'].split('.')[0] if split_ids else row['transcript_id']
            if use_version and row['transcript_version'] is not None:  # #152 bug #2
                feature_id += '.' + row['transcript_version']
            if use_name:
                # #152 bug #1 fixed: gate on and read the SAME attribute (transcript_name)
                if row['transcript_name'] is not None:
                    feature_name = row['transcript_name']
                elif row['name'] is not None:
                    feature_name = row['name']
                else:
                    continue

        if feature_id in mapping:
            continue
        mapping[feature_name if use_name else feature_id] = row['attr_value']

    return mapping


def _map_gene_to_attr_gff3(fields: pl.DataFrame, attribute: str, feature_type: str, use_name: bool,
                           split_ids: bool) -> dict:
    # In GFF3 the attribute lives on the feature it describes: gene-level attrs on 'gene' rows, transcript-level
    # attrs on transcript rows. The feature id is that row's own (prefix-stripped) ID. GFF3 has no *_version attrs.
    role_features = ('gene',) if feature_type == 'gene' else TRANSCRIPT_FEATURE_NAMES
    rows = fields.lazy().select(
        pl.col('_fields').list.get(_FEATURE_IDX).alias('feature'),
        pl.col('_fields').list.get(_ATTR_IDX).alias('_attr'),
    ).filter(pl.col('feature').is_in(role_features)).select(
        _extract_gff3_attribute(pl.col('_attr'), attribute).alias('attr_value'),
        _strip_so_prefix(_extract_gff3_attribute(pl.col('_attr'), 'ID')).alias('feature_id'),
        _extract_gff3_attribute(pl.col('_attr'), 'Name').alias('name'),
    ).collect()

    mapping = {}
    for row in rows.iter_rows(named=True):
        if row['attr_value'] is None:
            continue
        if row['feature_id'] is None:
            raise KeyError('ID')
        feature_id = row['feature_id'].split('.')[0] if split_ids else row['feature_id']
        if use_name:
            if row['name'] is None:
                continue
            feature_name = row['name']
        else:
            feature_name = None

        if feature_id in mapping:
            continue
        mapping[feature_name if use_name else feature_id] = row['attr_value']

    return mapping
