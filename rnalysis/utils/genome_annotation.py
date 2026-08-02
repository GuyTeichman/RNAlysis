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
from rnalysis.utils.param_typing import LEGAL_GENE_LENGTH_METHODS

# Read whole annotation lines as a single Polars column: pick a field separator that never occurs in a
# GTF/GFF file, so scan_csv yields one row per line and we split on the real '\t' ourselves. This lets us
# reproduce the two functions' divergent handling of malformed (wrong-column-count) lines exactly:
# map_transcripts_to_genes raises, while get_genomic_feature_lengths silently skips.
_LINE_SEP = '\x01'
_N_COLUMNS = 9
_FEATURE_IDX, _START_IDX, _END_IDX, _ATTR_IDX = 2, 3, 4, 8


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


def _scan_annotation_lines(path: Union[str, Path]) -> pl.LazyFrame:
    """Scan an annotation file into one row per (non-comment) line, split into its tab-separated fields.

    Yields columns ``_fields`` (List[str]) and ``_n`` (field count) so callers can enforce the 9-column
    rule themselves. ``quote_char=None`` keeps the literal ``"`` characters of GTF values intact.
    """
    return (pl.scan_csv(path, separator=_LINE_SEP, has_header=False, comment_prefix='#',
                        quote_char=None, schema={'raw': pl.String})
            .with_columns(pl.col('raw').str.split('\t').alias('_fields'))
            .with_columns(pl.col('_fields').list.len().alias('_n')))


def _read_annotation_records(path: Union[str, Path], file_type: Literal['gtf', 'gff3']) -> pl.DataFrame:
    """Materialize the feature records needed for length aggregation.

    Malformed lines (``_n != 9``) are silently dropped, reproducing ``get_genomic_feature_lengths``'s
    original skip-and-continue behaviour. Columns: ``feature``, ``start``, ``end``, ``_transcript`` (this
    feature's own id, for transcript/mRNA rows), ``_gene`` (its parent gene, for transcript/mRNA rows) and
    ``_parents`` (List[str] of parent transcript ids, for exon rows).
    """
    base = _scan_annotation_lines(path).filter(pl.col('_n') == _N_COLUMNS).select(
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
            _extract_gff3_attribute(pl.col('_attr'), 'ID').alias('_transcript'),
            _extract_gff3_attribute(pl.col('_attr'), 'Parent').alias('_gene'),
        ).with_columns(_extract_gff3_attribute(pl.col('_attr'), 'Parent').str.split(',').alias('_parents'))
    return records.collect()


def _sum_exon_lengths_per_parent(exons: pl.DataFrame) -> dict:
    """Sum exon lengths onto each parent transcript. A shared exon (multiple parents) is added to each,
    preserving the GFF3 comma-separated ``Parent=a,b`` semantics (correctness item #4)."""
    # exon _parents lists always hold >=1 parent, so empty_as_null is immaterial here; set it explicitly to
    # silence the Polars 1.x deprecation warning and pin behaviour across the 2.0 default change.
    grouped = exons.explode('_parents', empty_as_null=False).group_by('_parents').agg(pl.col('_length').sum())
    return dict(zip(grouped['_parents'].to_list(), grouped['_length'].to_list()))


def map_transcripts_to_genes(gtf_path: Union[str, Path], use_name: bool = False, use_version: bool = True,
                             split_ids: bool = True):
    validation.validate_genome_annotation_file(gtf_path, accept_gff3=False)

    df = _scan_annotation_lines(gtf_path).collect()
    malformed = df.filter(pl.col('_n') != _N_COLUMNS)
    if malformed.height:
        raise ValueError(f'Invalid GTF format: expected 9 tab-separated columns but got '
                         f'{malformed["_n"][0]} in line: "{malformed["raw"][0]}"')

    transcripts = df.lazy().select(
        pl.col('_fields').list.get(_FEATURE_IDX).alias('feature'),
        pl.col('_fields').list.get(_ATTR_IDX).alias('_attr'),
    ).filter(pl.col('feature') == 'transcript').with_columns(
        _extract_gtf_attribute(pl.col('_attr'), 'transcript_id').alias('transcript_id'),
        _extract_gtf_attribute(pl.col('_attr'), 'gene_id').alias('gene_id'),
        _extract_gtf_attribute(pl.col('_attr'), 'transcript_version').alias('transcript_version'),
        _extract_gtf_attribute(pl.col('_attr'), 'gene_version').alias('gene_version'),
        _extract_gtf_attribute(pl.col('_attr'), 'gene_name').alias('gene_name'),
    ).filter(pl.col('transcript_id').is_not_null() & pl.col('gene_id').is_not_null())

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
    validation.validate_genome_annotation_file(gtf_path, accept_gff3=False)
    assert feature_type in {'gene', 'transcript'}, f"Invalid feature_type: '{feature_type}'"

    # extract everything the (bug-for-bug) reduction below reads, one row per line, in file order
    rows = _scan_annotation_lines(gtf_path).filter(pl.col('_n') == _N_COLUMNS).with_columns(
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
            if use_version:
                if row['gene_version'] is None:
                    raise KeyError('gene_version')
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
            if use_version:
                if row['transcript_version'] is None:
                    raise KeyError('transcript_version')
                feature_id += '.' + row['transcript_version']
            if use_name:
                # BUG #1 (to be fixed in #152): gated on gene_name but the value read is transcript_name
                if row['gene_name'] is not None:
                    feature_name = row['transcript_name']
                elif row['name'] is not None:
                    feature_name = row['name']
                else:
                    continue

        if feature_id in mapping:
            continue
        mapping[feature_name if use_name else feature_id] = row['attr_value']

    return mapping
