import filecmp
import gzip
import os.path
import shutil

import pytest

from rnalysis import fastq
from rnalysis.fastq import *
from rnalysis.utils import io


def are_dir_trees_equal(dir1, dir2):
    """
    Compare two directories recursively. Files in each directory are \
    assumed to be equal if their names and contents are equal.\
    credit: bhttps://stackoverflow.com/a/6681395

    :param dir1: First directory path
    :param dir2: Second directory path

    :return: True if the dir trees are the same and there were no errors while accessing the directories or files, \
    False otherwise.
   """

    dirs_cmp = filecmp.dircmp(dir1, dir2)
    if len(dirs_cmp.left_only) > 0 or len(dirs_cmp.right_only) > 0 or \
        len(dirs_cmp.funny_files) > 0:
        return False
    (_, mismatch, errors) = filecmp.cmpfiles(
        dir1, dir2, dirs_cmp.common_files, shallow=False)
    if len(mismatch) > 0 or len(errors) > 0:
        return False
    for common_dir in dirs_cmp.common_dirs:
        new_dir1 = Path(dir1).joinpath(common_dir).as_posix()
        new_dir2 = Path(dir2).joinpath(common_dir).as_posix()
        if not are_dir_trees_equal(new_dir1, new_dir2):
            return False
    return True


@pytest.mark.parametrize("transcriptome_fasta,kallisto_installation_folder,kmer_length,make_unique,expected_command", [
    ('tests/test_files/kallisto_tests/transcripts.fasta', 'path/to/kallisto', 5, True,
     ['path/to/kallisto/kallisto', 'index', '-i', 'tests/test_files/kallisto_tests/transcripts.idx', '-k', '5',
      '--unique', 'tests/test_files/kallisto_tests/transcripts.fasta']),
    ('tests/test_files/kallisto_tests/transcripts.fasta', 'auto', 3, False,
     ['kallisto', 'index', '-i', 'tests/test_files/kallisto_tests/transcripts.idx', '-k', '3',
      'tests/test_files/kallisto_tests/transcripts.fasta']),
])
def test_kallisto_create_index_command(monkeypatch, transcriptome_fasta, kallisto_installation_folder, kmer_length,
                                       make_unique, expected_command):
    index_created = []

    def mock_run_subprocess(args, print_stdout=True, print_stderr=True, log_filename: str = None):
        assert args == expected_command
        assert print_stdout
        assert print_stderr
        index_created.append(True)

    monkeypatch.setattr(io, 'run_subprocess', mock_run_subprocess)

    kallisto_create_index(transcriptome_fasta, kallisto_installation_folder, kmer_length, make_unique)
    assert index_created == [True]


def test_kallisto_create_index():
    out_path = 'tests/test_files/kallisto_tests/transcripts.idx'
    truth_path = 'tests/test_files/kallisto_tests/transcripts_truth.idx'
    try:
        kallisto_create_index('tests/test_files/kallisto_tests/transcripts.fasta')
        with open(truth_path, 'rb') as truth, open(out_path, 'rb') as out:
            assert truth.read() == out.read()
    finally:
        if Path(out_path).exists():
            Path(out_path).unlink()


def test_kallisto_quantify_single_end():
    in_dir = 'tests/test_files/kallisto_tests'
    gtf_file = 'tests/test_files/kallisto_tests/transcripts.gtf'
    index_file = 'tests/test_files/kallisto_tests/transcripts_truth.idx'
    out_dir = 'tests/test_files/kallisto_tests/outdir'
    truth_dir = 'tests/test_files/kallisto_tests/truth/single'
    try:
        kallisto_quantify_single_end(in_dir, out_dir, index_file, gtf_file, 175.0, 25.0)
        for item in Path(out_dir).rglob('run_info.json'):
            item.unlink()
        assert are_dir_trees_equal(out_dir, truth_dir)
    finally:
        for item in Path(out_dir).iterdir():
            if 'gitignore' in item.name:
                continue
            if item.is_file():
                item.unlink()
            else:
                shutil.rmtree(item)


def test_kallisto_quantify_paired_end():
    in1_path = 'tests/test_files/kallisto_tests/reads_1.fastq.gz'
    in2_path = 'tests/test_files/kallisto_tests/reads_2.fastq.gz'
    gtf_file = 'tests/test_files/kallisto_tests/transcripts.gtf'
    index_file = 'tests/test_files/kallisto_tests/transcripts_truth.idx'
    out_dir = 'tests/test_files/kallisto_tests/outdir'
    truth_dir = 'tests/test_files/kallisto_tests/truth/paired'
    try:
        kallisto_quantify_paired_end([in1_path], [in2_path], out_dir, index_file, gtf_file)
        for item in Path(out_dir).rglob('run_info.json'):
            item.unlink()
        assert are_dir_trees_equal(out_dir, truth_dir)

    finally:
        for item in Path(out_dir).iterdir():
            if 'gitignore' in item.name:
                continue
            if item.is_file():
                item.unlink()
            else:
                shutil.rmtree(item)


@pytest.mark.parametrize(
    "fastq_folder,output_folder,index_file,gtf_file,average_fragment_length,stdev_fragment_length,"
    "kallisto_installation_folder,"
    "new_sample_names,stranded,learn_bias,seek_fusion_genes,bootstrap_samples,expected_command", [
        ('tests/test_files/kallisto_tests', 'tests/test_files/kallisto_tests/outdir',
         'tests/test_files/kallisto_tests/transcripts_truth.idx', 'tests/test_files/kallisto_tests/transcripts.gtf',
         125, 14, 'path/to/kallisto', 'auto', 'no', False, False, None,
         ['path/to/kallisto/kallisto', 'quant', '-i', 'tests/test_files/kallisto_tests/transcripts_truth.idx',
          '-o', 'outfolder', '--single', '-s', '14', '-l', '125']),
        ('tests/test_files/kallisto_tests', 'tests/test_files/kallisto_tests/outdir',
         'tests/test_files/kallisto_tests/transcripts_truth.idx', 'tests/test_files/kallisto_tests/transcripts.gtf',
         8.5, 0.2, 'path/to/kallisto', ['new_name_1', 'new_name_2'], 'reverse', True, True, 3,
         ['path/to/kallisto/kallisto', 'quant', '-i', 'tests/test_files/kallisto_tests/transcripts_truth.idx',
          '--bias', '--fusion', '--rf-stranded', '-b', '3', '-o', 'outfolder', '--single', '-s', '0.2', '-l', '8.5']),
    ])
def test_kallisto_quantify_single_end_command(monkeypatch, fastq_folder, output_folder, index_file, gtf_file,
                                              average_fragment_length, stdev_fragment_length,
                                              kallisto_installation_folder, new_sample_names, stranded, learn_bias,
                                              seek_fusion_genes, bootstrap_samples, expected_command):
    files_to_cover = ['reads_1.fastq.gz', 'reads_2.fastq.gz']
    file_stems = ['reads_1', 'reads_2']
    files_covered = []
    output_processed = []

    def mock_run_subprocess(args, print_stdout=True, print_stderr=True, log_filename: str = None):
        for i in range(len(files_to_cover)):
            if files_to_cover[i] in args[-1]:
                exp = expected_command.copy()
                if new_sample_names == 'auto':
                    exp[-6] = f'{output_folder}/{file_stems[i]}'
                else:
                    exp[-6] = f'{output_folder}/{new_sample_names[i]}'

                assert args == exp + [f'tests/test_files/kallisto_tests/{files_to_cover[i]}']
                files_covered.append(files_to_cover[i])
                break
        assert print_stdout
        assert print_stderr

    def mock_process_outputs(out_folder, gtf):
        assert Path(gtf_file) == Path(gtf)
        assert Path(out_folder) == Path(output_folder)
        output_processed.append(True)

    monkeypatch.setattr(io, 'run_subprocess', mock_run_subprocess)
    monkeypatch.setattr(fastq, '_process_kallisto_outputs', mock_process_outputs)

    kallisto_quantify_single_end(fastq_folder, output_folder, index_file, gtf_file, average_fragment_length,
                                 stdev_fragment_length,
                                 kallisto_installation_folder, new_sample_names, stranded, learn_bias,
                                 seek_fusion_genes, bootstrap_samples)
    assert sorted(files_covered) == sorted(files_to_cover)
    assert output_processed == [True]


@pytest.mark.parametrize(
    "r1_files,r2_files,output_folder,index_file,gtf_file,kallisto_installation_folder,"
    "new_sample_names,stranded,learn_bias,seek_fusion_genes,bootstrap_samples,expected_command", [
        (['tests/test_files/kallisto_tests/reads_1.fastq.gz'], ['tests/test_files/kallisto_tests/reads_2.fastq.gz'],
         'tests/test_files/kallisto_tests/outdir',
         'tests/test_files/kallisto_tests/transcripts_truth.idx', 'tests/test_files/kallisto_tests/transcripts.gtf',
         'path/to/kallisto', 'auto', 'no', False, False, None,
         ['path/to/kallisto/kallisto', 'quant', '-i', 'tests/test_files/kallisto_tests/transcripts_truth.idx',
          '-o', ]),
        (['tests/test_files/kallisto_tests/reads_1.fastq.gz'], ['tests/test_files/kallisto_tests/reads_2.fastq.gz'],
         'tests/test_files/kallisto_tests/outdir',
         'tests/test_files/kallisto_tests/transcripts_truth.idx', 'tests/test_files/kallisto_tests/transcripts.gtf',
         'path/to/kallisto', ['new_pair_name'], 'reverse', True, True, 3,
         ['path/to/kallisto/kallisto', 'quant', '-i', 'tests/test_files/kallisto_tests/transcripts_truth.idx',
          '--bias', '--fusion', '--rf-stranded', '-b', '3', '-o', ]),
    ])
def test_kallisto_quantify_paired_end_command(monkeypatch, r1_files, r2_files, output_folder, index_file, gtf_file,
                                              kallisto_installation_folder, new_sample_names, stranded, learn_bias,
                                              seek_fusion_genes, bootstrap_samples, expected_command):
    pairs_to_cover = [('reads_1.fastq.gz', 'reads_2.fastq.gz'), ]
    pair_stems = [('reads_1', 'reads_2'), ]
    pairs_covered = []
    output_processed = []

    def mock_run_subprocess(args, print_stdout=True, print_stderr=True, log_filename: str = None):
        for i in range(len(pairs_to_cover)):
            if pairs_to_cover[i][-1] in args[-1]:
                if new_sample_names == 'auto':
                    assert args == expected_command + [f'{output_folder}/{pair_stems[i][0]}_{pair_stems[i][1]}',
                                                       f'tests/test_files/kallisto_tests/{pairs_to_cover[i][0]}',
                                                       f'tests/test_files/kallisto_tests/{pairs_to_cover[i][1]}', ]
                else:
                    assert args == expected_command + [f'{output_folder}/{new_sample_names[i]}',
                                                       f'tests/test_files/kallisto_tests/{pairs_to_cover[i][0]}',
                                                       f'tests/test_files/kallisto_tests/{pairs_to_cover[i][1]}', ]
                pairs_covered.append(pairs_to_cover[i])
                break
        assert print_stdout
        assert print_stderr

    def mock_process_outputs(out_folder, gtf):
        assert Path(gtf_file) == Path(gtf)
        assert Path(out_folder) == Path(output_folder)
        output_processed.append(True)

    monkeypatch.setattr(io, 'run_subprocess', mock_run_subprocess)
    monkeypatch.setattr(fastq, '_process_kallisto_outputs', mock_process_outputs)

    kallisto_quantify_paired_end(r1_files, r2_files, output_folder, index_file, gtf_file,
                                 kallisto_installation_folder, new_sample_names, stranded, learn_bias,
                                 seek_fusion_genes, bootstrap_samples)
    assert sorted(pairs_covered) == sorted(pairs_to_cover)
    assert output_processed == [True]


def test_trim_adapters_single_end():
    in_dir = 'tests/test_files/test_fastqs/dir3'
    out_dir = 'tests/test_files/test_fastqs/outdir'
    out_path = 'tests/test_files/test_fastqs/outdir/test_fastq_trimmed.fastq.gz'
    truth_path = 'tests/test_files/test_fastqs/outdir/test_fastq_trimmed_truth.fastq.gz'
    adapter_seq = 'AACTTCTTA'
    try:
        trim_adapters_single_end(in_dir, out_dir, adapter_seq)
        with gzip.open(truth_path) as truth, gzip.open(out_path) as out:
            assert truth.read() == out.read()
    finally:
        Path(out_path).unlink()
        for file in Path(out_dir).iterdir():
            if file.is_file() and file.suffix == '.log':
                file.unlink()


def test_trim_adapters_paired_end():
    in1_path = 'tests/test_files/test_fastqs/dir4/paired_1.fastq'
    in2_path = 'tests/test_files/test_fastqs/dir4/paired_2.fastq'
    out_dir = 'tests/test_files/test_fastqs/outdir'
    out1_path = 'tests/test_files/test_fastqs/outdir/paired_1_trimmed.fastq.gz'
    out2_path = 'tests/test_files/test_fastqs/outdir/paired_2_trimmed.fastq.gz'
    truth1_path = 'tests/test_files/test_fastqs/outdir/paired_1_trimmed_truth.fastq.gz'
    truth2_path = 'tests/test_files/test_fastqs/outdir/paired_2_trimmed_truth.fastq.gz'
    adapter1_seq = 'TTA'
    adapter2_seq = 'TGT'
    try:
        trim_adapters_paired_end([in1_path], [in2_path], out_dir, adapter1_seq, adapter2_seq, minimum_read_length=5)

        with gzip.open(truth1_path) as truth, gzip.open(out1_path) as out:
            assert truth.read() == out.read()
        with gzip.open(truth2_path) as truth, gzip.open(out2_path) as out:
            assert truth.read() == out.read()
    finally:

        Path(out1_path).unlink()
        Path(out2_path).unlink()
        for file in Path(out_dir).iterdir():
            if file.is_file() and file.suffix == '.log':
                file.unlink()


@pytest.mark.parametrize('fastq_folder,output_folder,three_prime_adapters,five_prime_adapters,any_position_adapters,'
                         'quality_trimming,trim_n,minimum_read_length,maximum_read_length,discard_untrimmed_reads,'
                         'error_tolerance,minimum_overlap,allow_indels,parallel,expected_command', [
                             (
                                 'tests/test_files/test_fastqs/dir1', 'out/folder', 'ATGGAC', None, 'CCTGA', None,
                                 False, None, None, False, 0.1, 5, True, False,
                                 ['cutadapt', '--adapter', 'ATGGAC', '--anywhere', 'CCTGA', '--overlap', '5',
                                  '--error-rate', '0.1', '--output']),
                             (
                                 'tests/test_files/test_fastqs/dir1', 'out/folder', None, 'ATGGGG', ['CCTGA', 'ATGC'],
                                 15, True, 12, 100, True, 0.999, 2, False, True,
                                 ['cutadapt', '--front', 'ATGGGG', '--anywhere', 'CCTGA', '--anywhere', 'ATGC',
                                  '--quality-cutoff', '15', '--minimum-length', '12', '--maximum-length', '100',
                                  '--trim-n', '--discard-untrimmed', '--cores', '0', '--no-indels', '--overlap', '2',
                                  '--error-rate', '0.999', '--output'])
                         ])
def test_trim_adapters_single_end_command(monkeypatch, fastq_folder, output_folder, three_prime_adapters,
                                          five_prime_adapters, any_position_adapters, quality_trimming, trim_n,
                                          minimum_read_length, maximum_read_length, discard_untrimmed_reads,
                                          error_tolerance, minimum_overlap, allow_indels, parallel, expected_command):
    files_to_cover = ['fq1.fastq', 'fq2.fastq.gz']
    file_stems = ['fq1', 'fq2']
    files_covered = []

    def mock_run_subprocess(args, print_stdout=True, print_stderr=True, log_filename: str = None):
        for i in range(len(files_to_cover)):
            if files_to_cover[i] in args[-1]:
                assert args == expected_command + [f'{output_folder}/{file_stems[i]}_trimmed.fastq.gz',
                                                   f'tests/test_files/test_fastqs/dir1/{files_to_cover[i]}']
                files_covered.append(files_to_cover[i])
                break
        assert print_stdout
        assert print_stderr

    monkeypatch.setattr(io, 'run_subprocess', mock_run_subprocess)

    trim_adapters_single_end(fastq_folder, output_folder, three_prime_adapters,
                             five_prime_adapters, any_position_adapters, quality_trimming, trim_n,
                             minimum_read_length, maximum_read_length, discard_untrimmed_reads,
                             error_tolerance, minimum_overlap, allow_indels, parallel)
    assert sorted(files_covered) == sorted(files_to_cover)


@pytest.mark.parametrize(
    'fastq_1, fastq_2, output_folder, three_prime_r1, three_prime_r2,five_prime_r1, five_prime_r2, any_position_r1, '
    'any_position_r2,quality_trimming, trim_n, minimum_read_length, maximum_read_length,discard_untrimmed_reads, '
    'pair_filter_if,error_tolerance, minimum_overlap, allow_indels, parallel, expected_command',
    [
        (
            ['tests/test_files/test_fastqs/dir1/fq1.fastq', 'tests/test_files/test_fastqs/dir1/fq2.fastq.gz'],
            ['tests/test_files/test_fastqs/dir2/fq4.fq.gz', 'tests/test_files/test_fastqs/dir2/fq3.fq'],
            'out/folder', 'ATGGAC', 'AAATTTT', None, None, 'CCTGA', 'GTGGAA', None, False,
            None, None, False, 'any', 0.1, 5, True, False,
            ['cutadapt', '-a', 'ATGGAC', '-A', 'AAATTTT', '-b', 'CCTGA', '-B', 'GTGGAA', '--overlap', '5',
             '--error-rate', '0.1', '--pair-filter=any', '--output']),
        (
            ['tests/test_files/test_fastqs/dir1/fq1.fastq', 'tests/test_files/test_fastqs/dir1/fq2.fastq.gz'],
            ['tests/test_files/test_fastqs/dir2/fq4.fq.gz', 'tests/test_files/test_fastqs/dir2/fq3.fq'],
            'out/folder', None, None, None, 'ATGGGG', ['CCTGA', 'ATGC'], ['AAAA', 'GGGG'], 15,
            True, 12, 100, True, 'both', 0.999, 2, False, True,
            ['cutadapt', '-G', 'ATGGGG', '-b', 'CCTGA', '-b', 'ATGC', '-B', 'AAAA', '-B', 'GGGG',
             '--quality-cutoff', '15', '--minimum-length', '12', '--maximum-length', '100',
             '--trim-n', '--discard-untrimmed', '--cores', '0', '--no-indels', '--overlap', '2',
             '--error-rate', '0.999', '--pair-filter=both', '--output'])
    ])
def test_trim_adapters_paired_end_command(monkeypatch, fastq_1, fastq_2, output_folder, three_prime_r1, three_prime_r2,
                                          five_prime_r1, five_prime_r2, any_position_r1, any_position_r2,
                                          quality_trimming, trim_n, minimum_read_length, maximum_read_length,
                                          discard_untrimmed_reads, pair_filter_if: Literal['any', 'both', 'first'],
                                          error_tolerance, minimum_overlap, allow_indels, parallel, expected_command):
    pairs_to_cover = [('dir1/fq1.fastq', 'dir2/fq4.fq.gz'), ('dir1/fq2.fastq.gz', 'dir2/fq3.fq')]
    pair_stems = [('fq1', 'fq4'), ('fq2', 'fq3')]
    pairs_covered = []

    def mock_run_subprocess(args, print_stdout=True, print_stderr=True, log_filename: str = None):
        for i in range(len(pairs_to_cover)):
            if pairs_to_cover[i][-1] in args[-1]:
                assert args == expected_command + [f'{output_folder}/{pair_stems[i][0]}_trimmed.fastq.gz',
                                                   '--paired-output',
                                                   f'{output_folder}/{pair_stems[i][1]}_trimmed.fastq.gz',
                                                   f'tests/test_files/test_fastqs/{pairs_to_cover[i][0]}',
                                                   f'tests/test_files/test_fastqs/{pairs_to_cover[i][1]}', ]
                pairs_covered.append(pairs_to_cover[i])
                break
        assert print_stdout
        assert print_stderr

    monkeypatch.setattr(io, 'run_subprocess', mock_run_subprocess)

    trim_adapters_paired_end(fastq_1, fastq_2, output_folder, three_prime_r1, three_prime_r2,
                             five_prime_r1, five_prime_r2, any_position_r1, any_position_r2,
                             quality_trimming, trim_n, minimum_read_length, maximum_read_length,
                             discard_untrimmed_reads, pair_filter_if,
                             error_tolerance, minimum_overlap, allow_indels, parallel)
    assert sorted(pairs_covered) == sorted(pairs_to_cover)
