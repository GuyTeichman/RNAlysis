import warnings
from pathlib import Path
from typing import Literal, Union

import numpy as np
from scipy.stats.mstats import gmean

from rnalysis.utils import generic, parsing, validation
from rnalysis.utils.param_typing import LEGAL_GENE_LENGTH_METHODS


def map_transcripts_to_genes(gtf_path: Union[str, Path], use_name: bool = False, use_version: bool = True,
                             split_ids: bool = True):
    validation.validate_genome_annotation_file(gtf_path, accept_gff3=False)

    mapping = {}
    with open(gtf_path, errors="ignore") as f:
        for i, line in enumerate(f.readlines()):
            if len(line) == 0 or line[0] == '#':
                continue
            line_split = line.strip().split('\t')
            if len(line_split) != 9:
                raise ValueError(f'Invalid GTF format at line #{i + 1}: "{line}"')
            if line_split[2] == 'transcript':
                attributes = line_split[8]
                attributes_dict = parsing.parse_gtf_attributes(attributes)
                if 'transcript_id' not in attributes_dict or 'gene_id' not in attributes_dict:
                    continue

                transcript_id = attributes_dict['transcript_id'].split(".")[0] if split_ids else attributes_dict[
                    'transcript_id']
                gene_id = attributes_dict['gene_id'].split(".")[0] if split_ids else attributes_dict['gene_id']
                if use_version:
                    if 'transcript_version' in attributes_dict and 'gene_version' in attributes_dict:
                        transcript_id += '.' + attributes_dict['transcript_version']
                        gene_id += '.' + attributes_dict['gene_version']
                gene_name = None
                if use_name:
                    if 'gene_name' not in attributes_dict:
                        continue
                    gene_name = attributes_dict['gene_name']

                if transcript_id in mapping:
                    continue

                mapping[transcript_id] = gene_name if use_name else gene_id
    return mapping


def get_genomic_feature_lengths(gtf_path: Union[str, Path], feature_type: Literal['gene', 'transcript'],
                                method: Literal[LEGAL_GENE_LENGTH_METHODS]):
    method_funcs = {'mean': np.mean, 'median': np.median, 'max': np.max, 'min': np.min, 'geometric_mean': gmean,
                    'merged_exons': None}
    assert feature_type in ('gene', 'transcript'), f"Illegal feature_type: '{feature_type}'"
    assert method in method_funcs, f"Illegal method: '{method}'."
    file_type = validation.validate_genome_annotation_file(gtf_path)

    transcript_lengths = {}
    transcripts_to_genes = {}
    with open(gtf_path, errors="ignore") as f:
        for line in f.readlines():
            if line.startswith('#'):  # skip comments
                continue
            split = line.rstrip().split('\t')
            if len(split) != 9:  # skip invalid lines
                continue
            this_feature_type = split[2]
            attributes_dict = parsing.parse_gtf_attributes(split[8]) if file_type == 'gtf' else \
                parsing.parse_gff3_attributes(split[8])
            # map transcripts/mRNAs to their parent gene
            if this_feature_type == 'transcript' or this_feature_type == 'mRNA':
                transcript = attributes_dict['ID'] if file_type == 'gff3' else attributes_dict['transcript_id']
                gene_id = attributes_dict['Parent'] if file_type == 'gff3' else attributes_dict['gene_id']
                transcripts_to_genes[transcript] = gene_id
            # get exon length info
            elif this_feature_type == 'exon':
                start, end = int(split[3]), int(split[4])
                exon_length = end - start + 1
                transcript_ids = attributes_dict['Parent'] if file_type == 'gff3' else \
                    attributes_dict['transcript_id']
                # for 'merged exons' method, map the exons directly to the gene IDs, and
                # keep track of the interval of this exon. it will be used later.
                # No need to add each exon multiple times if it belongs to multiple transcripts.
                if method == 'merged_exons' and feature_type == 'gene':
                    gene_id = transcripts_to_genes[parsing.data_to_list(transcript_ids)[0]]
                    if gene_id not in transcript_lengths:
                        transcript_lengths[gene_id] = []
                    transcript_lengths[gene_id].append((start, end))
                else:
                    # map the exon lengths to their parent transcript/transcripts
                    for this_id in parsing.data_to_list(transcript_ids):
                        transcript_lengths[this_id] = transcript_lengths.get(this_id, 0) + exon_length

    if feature_type == 'transcript':
        warnings.warn("Since feature_type='transcript', the method parameter is ignored. ")
        return transcript_lengths

    elif method == 'merged_exons':
        # if the method is 'merged_exons', sum the intervals of the exons for each gene.
        # No need to map to transcripts since it was done earlier.
        summed_gene_lengths = dict()
        for gene_id, intervals in transcript_lengths.items():
            summed_gene_lengths[gene_id] = generic.sum_intervals_inclusive(intervals)
        return summed_gene_lengths

    # map transcripts to genes
    gene_lengths = {}
    for transcript_ids, gene_id in transcripts_to_genes.items():
        if gene_id not in gene_lengths:
            gene_lengths[gene_id] = []
        gene_lengths[gene_id].append(transcript_lengths[transcript_ids])

    # average transcript lengths per gene according to the specified function (method)
    avg_gene_lengths = {}
    avg_func = method_funcs[method]
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

    mapping = {}
    with open(gtf_path, errors="ignore") as f:
        for line in f.readlines():
            if len(line) == 0 or line[0] == '#':
                continue
            line_split = line.strip().split('\t')
            attributes = line_split[8]
            attributes_dict = parsing.parse_gtf_attributes(attributes)
            if attribute not in attributes_dict:
                continue

            feature_name = None
            if feature_type == 'gene':
                feature_id = attributes_dict['gene_id'].split(".")[0] if split_ids else attributes_dict['gene_id']
                if use_version:
                    feature_id += '.' + attributes_dict['gene_version']

                if use_name:
                    if 'gene_name' in attributes_dict:
                        feature_name = attributes_dict['gene_name']
                    elif 'name' in attributes_dict:
                        feature_name = attributes_dict['name']
                    else:
                        continue

            else:
                feature_id = attributes_dict['transcript_id'].split(".")[0] if split_ids else attributes_dict[
                    'transcript_id']
                if use_version:
                    feature_id += '.' + attributes_dict['transcript_version']

                if use_name:
                    if 'gene_name' in attributes_dict:
                        feature_name = attributes_dict['transcript_name']
                    elif 'name' in attributes_dict:
                        feature_name = attributes_dict['name']
                    else:
                        continue

            if feature_id in mapping:
                continue
            mapping[feature_name if use_name else feature_id] = attributes_dict[attribute]

    return mapping
