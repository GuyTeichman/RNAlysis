import warnings
from pathlib import Path
from typing import Union, List

import pandas as pd
from tqdm.auto import tqdm

from rnalysis import filtering
from rnalysis.utils import parsing, io

try:
    from typing import Literal
except ImportError:
    from typing_extensions import Literal
try:
    import cutadapt

    HAS_CUTADAPT = True

except ImportError:
    HAS_CUTADAPT = False

LEGAL_FASTQ_SUFFIXES = ['.fastq', '.fastq.gz', '.fq', '.fq.gz']


def kallisto_create_index(transcriptome_fasta: Union[str, Path],
                          kallisto_installation_folder: Union[str, Path, Literal['auto']] = 'auto',
                          kmer_length: int = 31, make_unique: bool = False):
    """
    Create an

    :param transcriptome_fasta:
    :type transcriptome_fasta:
    :param kallisto_installation_folder:
    :type kallisto_installation_folder:
    :param kmer_length:
    :type kmer_length:
    :param make_unique:
    :type make_unique:
    """
    assert isinstance(kmer_length, int), f"parameter 'kmer_length' must be an integer. Instead, got {type(kmer_length)}"
    assert 0 < kmer_length <= 31 and kmer_length % 2 == 1, f"'kmer_length' must be an odd integer between 1 and 31"

    if kallisto_installation_folder == 'auto':
        call = ['kallisto']
    else:
        call = [Path(kallisto_installation_folder).joinpath("kallisto").as_posix()]
    call.append('index')

    transcriptome_fasta = Path(transcriptome_fasta)
    assert Path(transcriptome_fasta).exists(), 'the transcriptome FASTA file does not exist!'

    index_filename = parsing.remove_suffixes(transcriptome_fasta).with_suffix('.idx').as_posix()
    call.extend(['-i', index_filename])

    call.extend(['-k', str(kmer_length)])

    if make_unique:
        call.append('--unique')

    call.append(transcriptome_fasta.as_posix())

    print(f"Running command: \n{' '.join(call)}")
    io.run_subprocess(call)
    print("Done")


def kallisto_quantify_single_end(fastq_folder: Union[str, Path], output_folder: Union[str, Path],
                                 index_file: Union[str, Path], gtf_file: Union[str, Path],
                                 average_fragment_length: float, stdev_fragment_length: float,
                                 kallisto_installation_folder: Union[str, Path, Literal['auto']] = 'auto',
                                 new_sample_names: Union[List[str], Literal['auto']] = 'auto',
                                 stranded: Literal['no', 'forward', 'reverse'] = 'no',
                                 learn_bias: bool = False, seek_fusion_genes: bool = False,
                                 bootstrap_samples: Union[int, None] = None):
    """
    Quantify transcript abundance in single-end mRNA sequencing data using kallisto.

    :param fastq_folder: Path to the folder containing the FASTQ files you want to quantify
    :type fastq_folder: str or Path
    :param output_folder: Path to a folder in which the quantified results, as well as the log files, will be saved. \
    The individual output of each pair of FASTQ files will reside in a different sub-folder within the output folder, \
    and a summarized results table will be saved in the output folder itself.
    :type output_folder: str/Path to an existing folder
    :param index_file: Path to a pre-built kallisto index of the target transcriptome. \
    Can either be downloaded from  the \
    `kallisto transcriptome indices site <https://github.com/pachterlab/kallisto-transcriptome-indices/releases>`_, \
    or generated manually from a FASTA file using the function `kallisto_create_index`.
    :type index_file: str or Path
    :param gtf_file: Path to a GTF annotation file. This file will be used to map per-transcript abundances to \
    per-gene estimated counts. The transcript names in the GTF files should match the ones in the index file - \
    we recommend downloading cDNA FASTA/index files and GTF files from the same data source.
    :type gtf_file: str or Path
    :param average_fragment_length: Estimated average fragment length. Typical Illumina libraries produce fragment \
    lengths ranging from 180–200bp, but it’s best to determine this from a library quantification with an instrument \
    such as an Agilent Bioanalyzer.
    :type average_fragment_length: float > 0
    :param stdev_fragment_length: Estimated standard deviation of fragment length. Typical Illumina libraries \
    produce fragment lengths ranging from 180–200bp, but it’s best to determine this from a library quantification \
    with an instrument such as an Agilent Bioanalyzer.
    :type stdev_fragment_length: float > 0
    :param kallisto_installation_folder: Path to the installation folder of kallisto. For example: \
    'C:/Program Files/kallisto'
    :type kallisto_installation_folder: str, Path, or 'auto' (default='auto')
    :param new_sample_names: Give a new name to each quantified sample (optional). \
    If sample_names='auto', sample names \
    will be given automatically. Otherwise, sample_names should be a list of new names, with the order of the names \
    matching the order of the file pairs.
    :type new_sample_names: list of str or 'auto' (default='auto')
    :param stranded: Indicates the strandedness of the data. 'no' indicates the data is not stranded. \
    'forward' indicates the data is stranded, where the first read in the pair pseudoaligns \
    to the forward strand of a transcript. 'reverse' indicates the data is stranded, where the first read in the pair \
    pseudoaligns to the reverse strand of a transcript.
    :type stranded: 'no', 'forward', 'reverse' (default='no')
    :param learn_bias: if True, kallisto learns parameters for a model of sequences specific bias \
    and corrects the abundances accordlingly.
    :type learn_bias: bool (default=False)
    :param seek_fusion_genes: if True, does normal quantification, but additionally looks for reads that do not \
    pseudoalign because they are potentially from fusion genes. \
    All output is written to the file fusion.txt in the output folder.
    :type seek_fusion_genes: bool (default=False)
    :param bootstrap_samples: Number of bootstrap samples to be generated. Bootstrap samples do not affect the \
    estimated count values, but generates an additional .hdf5 output file which contains \
    uncertainty estimates for the expression levels. This is required if you use tools such as Sleuth for downstream \
    differential expression analysis, but not for more traditional tools such as DESeq2 and edgeR.
    :type bootstrap_samples: int >=0 or None (default=None)
    """
    if new_sample_names != 'auto':
        new_sample_names = parsing.data_to_list(new_sample_names)
    output_folder = Path(output_folder)

    call = _parse_kallisto_misc_args(output_folder, index_file, kallisto_installation_folder, stranded, learn_bias,
                                     seek_fusion_genes, bootstrap_samples)

    call.append("--single")
    call.extend(["-s", str(stdev_fragment_length)])
    call.extend(["-l", str(average_fragment_length)])

    legal_samples = []
    for item in Path(fastq_folder).iterdir():
        if item.is_file():
            name = item.name
            if any([name.endswith(suffix) for suffix in LEGAL_FASTQ_SUFFIXES]):
                legal_samples.append(item)

    assert (new_sample_names == 'auto') or (len(new_sample_names) == len(legal_samples)), \
        f'Number of samples ({len(legal_samples)}) does not match number of sample names ({len(new_sample_names)})!'

    calls = []
    for i, item in enumerate(sorted(legal_samples)):
        this_call = call.copy()
        if new_sample_names == 'auto':
            this_name = parsing.remove_suffixes(item).stem
        else:
            this_name = new_sample_names[i]

        this_call[-6] = Path(this_call[-6]).joinpath(this_name).as_posix()
        this_call.append(item.as_posix())
        calls.append(this_call)

    for kallisto_call in tqdm(calls, 'quantifying transcript abundance', unit='files'):
        io.run_subprocess(kallisto_call)
        print(f"File saved successfully at {kallisto_call[-7]}")

    return _process_kallisto_outputs(output_folder, gtf_file)


def kallisto_quantify_paired_end(r1_files: List[str], r2_files: List[str], output_folder: Union[str, Path],
                                 index_file: Union[str, Path], gtf_file: Union[str, Path],
                                 kallisto_installation_folder: Union[str, Path, Literal['auto']] = 'auto',
                                 new_sample_names: Union[List[str], Literal['auto']] = 'auto',
                                 stranded: Literal['no', 'forward', 'reverse'] = 'no', learn_bias: bool = False,
                                 seek_fusion_genes: bool = False, bootstrap_samples: Union[int, None] = None):
    """
    Quantify transcript abundance in paired-end mRNA sequencing data using kallisto.

    :param r1_files: a list of paths to your Read#1 files. The files should be sorted in tandem with r2_files, \
    so that they line up to form pairs of R1 and R2 files.
    :type r1_files: list of str/Path to existing FASTQ files
    :param r2_files: a list of paths to your Read#2 files. The files should be sorted in tandem with r1_files, \
    so that they line up to form pairs of R1 and R2 files.
    :type r2_files: list of str/Path to existing FASTQ files
    :param output_folder: Path to a folder in which the quantified results, as well as the log files, will be saved. \
    The individual output of each pair of FASTQ files will reside in a different sub-folder within the output folder, \
    and a summarized results table will be saved in the output folder itself.
    :type output_folder: str/Path to an existing folder
    :param index_file: Path to a pre-built kallisto index of the target transcriptome. \
    Can either be downloaded from  the \
    `kallisto transcriptome indices site <https://github.com/pachterlab/kallisto-transcriptome-indices/releases>`_, \
    or generated manually from a FASTA file using the function `kallisto_create_index`.
    :type index_file: str or Path
    :param gtf_file: Path to a GTF annotation file. This file will be used to map per-transcript abundances to \
    per-gene estimated counts. The transcript names in the GTF files should match the ones in the index file - \
    we recommend downloading cDNA FASTA/index files and GTF files from the same data source.
    :type gtf_file: str or Path
    :param kallisto_installation_folder: Path to the installation folder of kallisto. For example: \
        'C:/Program Files/kallisto'
    :type kallisto_installation_folder: str, Path, or 'auto' (default='auto')
    :param new_sample_names: Give a new name to each quantified sample (optional). If sample_names='auto', sample names \
    will be given automatically. Otherwise, sample_names should be a list of new names, with the order of the names \
    matching the order of the file pairs.
    :type new_sample_names: list of str or 'auto' (default='auto')
    :param stranded: Indicates the strandedness of the data. 'no' indicates the data is not stranded. \
    'forward' indicates the data is stranded, where the first read in the pair pseudoaligns \
    to the forward strand of a transcript. 'reverse' indicates the data is stranded, where the first read in the pair \
    pseudoaligns to the reverse strand of a transcript.
    :type stranded: 'no', 'forward', 'reverse' (default='no')
    :param learn_bias: if True, kallisto learns parameters for a model of sequences specific bias \
    and corrects the abundances accordlingly.
    :type learn_bias: bool (default=False)
    :param seek_fusion_genes: if True, does normal quantification, but additionally looks for reads that do not \
    pseudoalign because they are potentially from fusion genes. \
    All output is written to the file fusion.txt in the output folder.
    :type seek_fusion_genes: bool (default=False)
    :param bootstrap_samples: Number of bootstrap samples to be generated. Bootstrap samples do not affect the \
    estimated count values, but generates an additional .hdf5 output file which contains \
    uncertainty estimates for the expression levels. This is required if you use tools such as Sleuth for downstream \
    differential expression analysis, but not for more traditional tools such as DESeq2 and edgeR.
    :type bootstrap_samples: int >=0 or None (default=None)
    """
    if new_sample_names != 'auto':
        new_sample_names = parsing.data_to_list(new_sample_names)
    assert len(r1_files) == len(r2_files), f"Got an uneven number of R1 and R2 files: " \
                                           f"{len(r1_files)} and {len(r2_files)} respectively"
    assert (new_sample_names == 'auto') or (len(new_sample_names) == len(
        r1_files)), f'Number of samples ({len(r1_files)}) does not match number of sample names ({len(new_sample_names)})!'

    output_folder = Path(output_folder)
    call = _parse_kallisto_misc_args(output_folder, index_file, kallisto_installation_folder, stranded, learn_bias,
                                     seek_fusion_genes, bootstrap_samples)

    calls = []
    for i, (file1, file2) in enumerate(zip(r1_files, r2_files)):
        file1 = Path(file1)
        file2 = Path(file2)
        this_call = call.copy()
        if new_sample_names == 'auto':
            this_name = f"{parsing.remove_suffixes(file1).stem}_{parsing.remove_suffixes(file2).stem}"
        else:
            this_name = new_sample_names[i]

        this_call[-1] = Path(this_call[-1]).joinpath(this_name).as_posix()
        this_call.extend([file1.as_posix(), file2.as_posix()])
        calls.append(this_call)

    for kallisto_call in tqdm(calls, 'Quantifying transcript abundance', unit='file pairs'):
        io.run_subprocess(kallisto_call)
        print(f"Files saved successfully at {kallisto_call[-3]}")

    return _process_kallisto_outputs(output_folder, gtf_file)


def _process_kallisto_outputs(output_folder, gtf_file):
    counts, tpm = _merge_kallisto_outputs(output_folder)
    genes_scaled_tpm = _sum_transcripts_to_genes(tpm, counts, gtf_file)

    io.save_csv(counts, output_folder.joinpath('transcript_counts.csv'))
    io.save_csv(tpm, output_folder.joinpath('transcript_tpm.csv'))
    io.save_csv(genes_scaled_tpm, output_folder.joinpath('kallisto_output_scaled_per_gene.csv'))

    print("Done")
    return filtering.CountFilter.from_dataframe(genes_scaled_tpm, 'kallisto_output_scaled_per_gene',
                                                is_normalized=False)


def _parse_kallisto_misc_args(output_folder, index_file: str, kallisto_installation_folder: Union[str, Path],
                              stranded: Literal['no', 'forward', 'reverse'] = 'no', learn_bias: bool = False,
                              seek_fusion_genes: bool = False, bootstrap_samples: Union[int, None] = None):
    output_folder = Path(output_folder)
    index_file = Path(index_file)
    assert output_folder.exists(), "supplied 'output_folder' does not exist!"
    assert index_file.exists(), f"supplied 'index_file' does not exist!"
    assert isinstance(stranded, str) and stranded.lower() in ["no", "forward", "reverse"], \
        f"invalid value for parameter 'stranded': {stranded}"

    if kallisto_installation_folder == 'auto':
        call = ['kallisto']
    else:
        call = [Path(kallisto_installation_folder).joinpath("kallisto").as_posix()]
    call.append('quant')
    call.extend(["-i", index_file.as_posix()])

    if learn_bias:
        call.append("--bias")

    if seek_fusion_genes:
        call.append("--fusion")

    stranded = stranded.lower()
    if stranded == "forward":
        call.append("--fr-stranded")
    elif stranded == "reverse":
        call.append("--rf-stranded")

    if bootstrap_samples is not None:
        assert isinstance(bootstrap_samples,
                          int) and bootstrap_samples >= 0, "'bootstrap_samples' must be a non-negative integer!"
        call.extend(['-b', str(bootstrap_samples)])

    call.extend(["-o", output_folder.as_posix()])
    return call


def _merge_kallisto_outputs(output_folder: Union[str, Path]):
    """
    output a merged csv file of transcript estimated counts, and a merged csv file of transcript estimated TPMs.

    :param output_folder:
    :type output_folder:
    :return:
    :rtype:
    """
    counts = pd.DataFrame()
    tpm = pd.DataFrame()
    for item in Path(output_folder).iterdir():
        if item.is_dir():
            abundance_path = item.joinpath('abundance.tsv')
            if abundance_path.exists():
                this_df = io.load_csv(abundance_path, index_col=0)
                sample_name = item.name
                counts[sample_name] = this_df['est_counts']
                tpm[sample_name] = this_df['tpm']
    return counts, tpm


def _sum_transcripts_to_genes(tpm: pd.DataFrame, counts: pd.DataFrame, gtf_path: Union[str, Path]):
    transcript_to_gene_map = _map_transcripts_to_genes(gtf_path)
    library_sizes = counts.sum(axis=0) / (10 ** 6)
    tpm_cpy = tpm.copy()
    tpm_cpy['Gene ID'] = transcript_to_gene_map
    tpm_by_gene = tpm_cpy.groupby('Gene ID').sum()
    scaled_tpm = tpm_by_gene.multiply(library_sizes, axis=1)

    return scaled_tpm


def _map_transcripts_to_genes(gtf_file: Union[str, Path], use_name: bool = False, use_version: bool = True):
    mapping = {}
    with open(gtf_file, errors="ignore") as f:
        for line in f.readlines():
            if len(line) == 0 or line[0] == '#':
                continue
            line_split = line.strip().split('\t')
            if line_split[2] == 'transcript':
                info = line_split[8]
                d = {}
                for x in info.split('; '):
                    x = x.strip()
                    p = x.find(' ')
                    if p == -1:
                        continue
                    k = x[:p]
                    p = x.find('"', p)
                    p2 = x.find('"', p + 1)
                    v = x[p + 1:p2]
                    d[k] = v

                if 'transcript_id' not in d or 'gene_id' not in d:
                    continue

                transcript_id = d['transcript_id'].split(".")[0]
                gene_id = d['gene_id'].split(".")[0]
                if use_version:
                    if 'transcript_version' in d and 'gene_version' in d:
                        transcript_id += '.' + d['transcript_version']
                        gene_id += '.' + d['gene_version']
                gene_name = None
                if use_name:
                    if 'gene_name' not in d:
                        continue
                    gene_name = d['gene_name']

                if transcript_id in mapping:
                    continue

                mapping[transcript_id] = gene_name if use_name else gene_id
    return mapping


def trim_adapters_single_end(fastq_folder: Union[str, Path], output_folder: Union[str, Path],
                             three_prime_adapters: Union[None, str, List[str]],
                             five_prime_adapters: Union[None, str, List[str]] = None,
                             any_position_adapters: Union[None, str, List[str]] = None,
                             quality_trimming: Union[int, None] = 20, trim_n: bool = True,
                             minimum_read_length: Union[int, None] = 10, maximum_read_length: Union[int, None] = None,
                             discard_untrimmed_reads: bool = True, error_tolerance: float = 0.1,
                             minimum_overlap: int = 3, allow_indels: bool = True, parallel: bool = True):
    """
    Trim adapters from single-end reads using `CutAdapt <https://cutadapt.readthedocs.io/>`_.

    :param fastq_folder: Path to the folder containing your untrimmed FASTQ files
    :type fastq_folder: str/Path to an existing folder
    :param output_folder: Path to a folder in which the trimmed FASTQ files, as well as the log files, will be saved.
    :type output_folder: str/Path to an existing folder
    :param three_prime_adapters: the sequence of the adapter/adapters to trim from the 3' end of the reads.
    :type three_prime_adapters: str, list of str, or None
    :param five_prime_adapters: the sequence of the adapter/adapters to trim from the 5' end of the reads.
    :type five_prime_adapters: str, list of str, or None (default=None)
    :param any_position_adapters: the sequence of the adapter/adapters to trim from either end \
    (or from the middle) of the reads.
    :type any_position_adapters: str, list of str, or None (default=None)
    :param quality_trimming: if specified, trim low-quality 3' end from the reads. Any bases with quality score below \
    the specified value will be trimmed from the 3' end of the read.
    :type quality_trimming: int or None (default=20)
    :param trim_n: if True, removem flanking N bases from each read. For example, a read with the sequence \
    'NNACGTACGTNNNN' will be trimmed down to 'ACGTACGT'. This occurs after adapter trimming.
    :type trim_n: bool (default=True)
    :param minimum_read_length: if specified (default), \
    discard processed reads that are shorter than minimum_read_length.
    :type minimum_read_length: int or None (default=10)
    :param maximum_read_length: if specified, discard processed reads that are shorter than minimum_read_length.
    :type maximum_read_length: int or None (default=None)
    :param discard_untrimmed_reads: if True, discards reads in which no adapter was found.
    :type discard_untrimmed_reads: bool (default=True)
    :param error_tolerance: The level of error tolerance permitted when searching for adapters, with the lowest \
    value being 0 (no error tolerance) and the maximum being 1 (100% error tolerance). \
    Allowed errors are mismatches, insertions and deletions.
    :type error_tolerance: float between 0 and 1 (default=0.1)
    :param minimum_overlap: the minimum number of nucleotides that must match exactly to the adapter sequence \
    in order to trim it.
    :type minimum_overlap: int >= 0 (default=3)
    :param allow_indels: if False, insertions and deletions in the adapter sequence are not allowed - only mismatches.
    :type allow_indels: bool (default=True)
    :param parallel: if True, runs CutAdapt on all available cores in parallel. \
    Otherwise, run CutAdapt on a single processor only.
    :type parallel: bool (default=True)
    """
    if not HAS_CUTADAPT:
        warnings.warn(f"Python package 'cutadapt' is not installed. \n"
                      f"If you want to use the adapter trimming feature, "
                      f"please install python package 'cutadapt' and try again. ")
        return

    call = ['cutadapt']

    for adapter_group, prefix in zip([three_prime_adapters, five_prime_adapters, any_position_adapters],
                                     ['--adapter', '--front', '--anywhere']):
        if adapter_group is not None:
            for adapter in parsing.data_to_list(adapter_group):
                assert isinstance(adapter, str), f"The following adapter is invalid: {adapter}"
                call.extend([prefix, adapter])

    call.extend(_parse_cutadapt_misc_args(quality_trimming, trim_n, minimum_read_length, maximum_read_length,
                                          discard_untrimmed_reads, error_tolerance, minimum_overlap, allow_indels,
                                          parallel))

    calls = []
    for item in Path(fastq_folder).iterdir():
        if item.is_file():
            name = item.name
            if any([name.endswith(suffix) for suffix in LEGAL_FASTQ_SUFFIXES]):
                stem = parsing.remove_suffixes(item).stem
                output_path = Path(output_folder).joinpath(stem + '_trimmed.fastq.gz')
                this_call = call.copy()
                this_call.extend(['--output', output_path.as_posix(), item.as_posix()])
                calls.append(this_call)

    for cutadapt_call in tqdm(calls, 'Trimming adapters', unit='files'):
        infile_stem = parsing.remove_suffixes(Path(cutadapt_call[-1])).stem
        log_filename = Path(output_folder).joinpath(f'cutadapt_log_{infile_stem}.log').absolute().as_posix()
        io.run_subprocess(cutadapt_call, log_filename=log_filename)
        print(f"File saved successfully at {cutadapt_call[-2]}")
    print("Done")


def trim_adapters_paired_end(r1_files: List[Union[str, Path]], r2_files: List[Union[str, Path]],
                             output_folder: Union[str, Path],
                             three_prime_adapters_r1: Union[None, str, List[str]],
                             three_prime_adapters_r2: Union[None, str, List[str]],
                             five_prime_adapters_r1: Union[None, str, List[str]] = None,
                             five_prime_adapters_r2: Union[None, str, List[str]] = None,
                             any_position_adapters_r1: Union[None, str, List[str]] = None,
                             any_position_adapters_r2: Union[None, str, List[str]] = None,
                             quality_trimming: Union[int, None] = 20, trim_n: bool = True,
                             minimum_read_length: Union[int, None] = 10, maximum_read_length: Union[int, None] = None,
                             discard_untrimmed_reads: bool = True,
                             pair_filter_if: Literal['both', 'any', 'first'] = 'both',
                             error_tolerance: float = 0.1, minimum_overlap: int = 3, allow_indels: bool = True,
                             parallel: bool = True):
    """
    Trim adapters from paired-end reads using `CutAdapt <https://cutadapt.readthedocs.io/>`_.

    :param r1_files: a list of paths to your Read#1 files. The files should be sorted in tandem with r2_files, \
    so that they line up to form pairs of R1 and R2 files.
    :type r1_files: list of str/Path to existing FASTQ files
    :param r2_files: a list of paths to your Read#2 files. The files should be sorted in tandem with r1_files, \
    so that they line up to form pairs of R1 and R2 files.
    :type r2_files: list of str/Path to existing FASTQ files
    :param output_folder: Path to a folder in which the trimmed FASTQ files, as well as the log files, will be saved.
    :type output_folder: str/Path to an existing folder
    :param three_prime_adapters_r1: the sequence of the adapter/adapters \
    to trim from the 3' end of the reads in Read#1 files.
    :type three_prime_adapters_r1: str, list of str, or None
    :param three_prime_adapters_r2: the sequence of the adapter/adapters \
    to trim from the 3' end of the reads in Read#2 files.
    :type three_prime_adapters_r2: str, list of str, or None
    :param five_prime_adapters_r1: the sequence of the adapter/adapters \
    to trim from the 5' end of the reads in Read#1 files.
    :type five_prime_adapters_r1: str, list of str, or None (default=None)
    :param five_prime_adapters_r2: the sequence of the adapter/adapters \
    to trim from the 5' end of the reads in Read#2 files.
    :type five_prime_adapters_r2: str, list of str, or None (default=None)
    :param any_position_adapters_r1: the sequence of the adapter/adapters \
    to trim from either end (or the middle) of the reads in Read#1 files.
    :type any_position_adapters_r1: str, list of str, or None (default=None)
    :param any_position_adapters_r2: the sequence of the adapter/adapters \
    to trim from either end (or the middle) of the reads in Read#2 files.
    :type any_position_adapters_r2: str, list of str, or None (default=None)
    :param quality_trimming: if specified, trim low-quality 3' end from the reads. Any bases with quality score below \
    the specified value will be trimmed from the 3' end of the read.
    :type quality_trimming: int or None (default=20)
    :param trim_n: if True, removem flanking N bases from each read. For example, a read with the sequence \
    'NNACGTACGTNNNN' will be trimmed down to 'ACGTACGT'. This occurs after adapter trimming.
    :type trim_n: bool (default=True)
    :param minimum_read_length: if specified (default), \
    discard processed reads that are shorter than minimum_read_length.
    :type minimum_read_length: int or None (default=10)
    :param maximum_read_length: if specified, discard processed reads that are shorter than minimum_read_length.
    :type maximum_read_length: int or None (default=None)
    :param discard_untrimmed_reads: if True, discards reads in which no adapter was found.
    :type discard_untrimmed_reads: bool (default=True)
    :param pair_filter_if: Cutadapt always discards both reads of a pair if it determines that the pair \
    should be discarded. This parameter determines how to combine the filters for Read#1 and Read#2 into a \
    single decision about the read pair. When the value is 'both', you require that filtering criteria must apply \
    to both reads in order for a read pair to be discarded. When the value is 'any', you require that \
    at least one of the reads (R1 or R2) fulfills the filtering criterion in order to discard them. \
    When the value is 'first', only the first read in each pair determines whether to discard the pair or not.
    :type pair_filter_if: 'both', 'any', or 'first' (default='both')
    :param error_tolerance: The level of error tolerance permitted when searching for adapters, with the lowest \
    value being 0 (no error tolerance) and the maximum being 1 (100% error tolerance). \
    Allowed errors are mismatches, insertions and deletions.
    :type error_tolerance: float between 0 and 1 (default=0.1)
    :param minimum_overlap: the minimum number of nucleotides that must match exactly to the adapter sequence \
    in order to trim it.
    :type minimum_overlap: int >= 0 (default=3)
    :param allow_indels: if False, insertions and deletions in the adapter sequence are not allowed - only mismatches.
    :type allow_indels: bool (default=True)
    :param parallel: if True, runs CutAdapt on all available cores in parallel. \
    Otherwise, run CutAdapt on a single processor only.
    :type parallel: bool (default=True)
    """
    if not HAS_CUTADAPT:
        warnings.warn(f"Python package 'cutadapt' is not installed. \n"
                      f"If you want to use the adapter trimming feature, "
                      f"please install python package 'cutadapt' and try again. ")
        return
    assert len(r1_files) == len(r2_files), f"Got an uneven number of R1 and R2 files: " \
                                           f"{len(r1_files)} and {len(r2_files)} respectively"

    call = ['cutadapt']

    for r1_group, r2_group, prefix in zip(
        [three_prime_adapters_r1, five_prime_adapters_r1, any_position_adapters_r1],
        [three_prime_adapters_r2, five_prime_adapters_r2, any_position_adapters_r2],
        ['-a', '-g', '-b']):

        if r1_group is not None:
            for adapter in parsing.data_to_list(r1_group):
                assert isinstance(adapter, str), f"The following R1 adapter is invalid: {adapter}"
                call.extend([prefix, adapter])

        if r2_group is not None:
            for adapter in parsing.data_to_list(r2_group):
                assert isinstance(adapter, str), f"The following R2 adapter is invalid: {adapter}"
                call.extend([prefix.upper(), adapter])

    call.extend(_parse_cutadapt_misc_args(quality_trimming, trim_n, minimum_read_length, maximum_read_length,
                                          discard_untrimmed_reads, error_tolerance, minimum_overlap, allow_indels,
                                          parallel))
    call.append(f'--pair-filter={pair_filter_if}')
    calls = []
    for file1, file2 in zip(r1_files, r2_files):
        file1 = Path(file1)
        file2 = Path(file2)
        output_path_r1 = Path(output_folder).joinpath(parsing.remove_suffixes(file1).stem + '_trimmed.fastq.gz')
        output_path_r2 = Path(output_folder).joinpath(parsing.remove_suffixes(file2).stem + '_trimmed.fastq.gz')

        this_call = call.copy()
        this_call.extend(['--output', output_path_r1.as_posix()])
        this_call.extend(['--paired-output', output_path_r2.as_posix()])
        this_call.extend([file1.as_posix(), file2.as_posix()])
        calls.append(this_call)

    for cutadapt_call in tqdm(calls, 'Trimming adapters', unit='file pairs'):
        infile1_stem = parsing.remove_suffixes(Path(cutadapt_call[-2])).stem
        infile2_stem = parsing.remove_suffixes(Path(cutadapt_call[-1])).stem
        log_filename = Path(output_folder).joinpath(
            f'cutadapt_log_{infile1_stem}_{infile2_stem}.log').absolute().as_posix()
        io.run_subprocess(cutadapt_call, log_filename=log_filename)
        print(f"Files saved successfully at {cutadapt_call[-2]} and  {cutadapt_call[-1]}")
    print("Done")


def _parse_cutadapt_misc_args(quality_trimming: Union[int, None], trim_n: bool,
                              minimum_read_length: Union[int, None],
                              maximum_read_length: Union[int, None],
                              discard_untrimmed_reads: bool, error_tolerance: float, minimum_overlap: int,
                              allow_indels: bool, parallel: bool) -> List[str]:
    call = []
    if quality_trimming is not None:
        assert isinstance(quality_trimming, int), f"'quality_trimming' must be an integer. " \
                                                  f"Instead, got type {type(quality_trimming)}"
        call.extend(['--quality-cutoff', str(quality_trimming)])

    if minimum_read_length is not None:
        assert isinstance(minimum_read_length, int), f"'minimum_read_length' must be an integer. " \
                                                     f"Instead, got type {type(minimum_read_length)}"
        call.extend(['--minimum-length', str(minimum_read_length)])

    if maximum_read_length is not None:
        assert isinstance(maximum_read_length, int), f"'maximum_read_length' must be an integer. " \
                                                     f"Instead, got type {type(maximum_read_length)}"
        call.extend(['--maximum-length', str(maximum_read_length)])
    if trim_n:
        call.append('--trim-n')

    if discard_untrimmed_reads:
        call.append('--discard-untrimmed')

    if parallel:
        call.extend(['--cores', '0'])

    if not allow_indels:
        call.append('--no-indels')

    call.extend(['--overlap', str(minimum_overlap)])
    call.extend(['--error-rate', str(error_tolerance)])

    return call
