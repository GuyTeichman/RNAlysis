import warnings
from pathlib import Path
from typing import Union, List
from tqdm.auto import tqdm
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

    legal_suffixes = ['.fastq', '.fastq.gz', '.fq', '.fq.gz']
    calls = []
    for item in Path(fastq_folder).iterdir():
        if item.is_file():
            name = item.name
            if any([name.endswith(suffix) for suffix in legal_suffixes]):
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
