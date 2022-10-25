import pytest
from rnalysis.fastq import *
from rnalysis.utils import io


@pytest.mark.parametrize('fastq_folder,outout_folder,three_prime_adapters,five_prime_adapters,any_position_adapters,'
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
def test_trim_adapters_single_end_command(monkeypatch, fastq_folder, outout_folder, three_prime_adapters,
                                          five_prime_adapters, any_position_adapters, quality_trimming, trim_n,
                                          minimum_read_length, maximum_read_length, discard_untrimmed_reads,
                                          error_tolerance, minimum_overlap, allow_indels, parallel, expected_command):
    files_to_cover = ['fq1.fastq', 'fq2.fastq.gz']
    file_stems = ['fq1', 'fq2']
    files_covered = []

    def mock_run_subprocess(args, print_stdout=True, print_stderr=True):
        for i in range(len(files_to_cover)):
            if files_to_cover[i] in args[-1]:
                assert args == expected_command + [f'{outout_folder}/{file_stems[i]}_trimmed.fastq.gz',
                                                   f'tests/test_files/test_fastqs/dir1/{files_to_cover[i]}']
                files_covered.append(files_to_cover[i])
                break
        assert print_stdout
        assert print_stderr

    monkeypatch.setattr(io, 'run_subprocess', mock_run_subprocess)

    trim_adapters_single_end(fastq_folder, outout_folder, three_prime_adapters,
                             five_prime_adapters, any_position_adapters, quality_trimming, trim_n,
                             minimum_read_length, maximum_read_length, discard_untrimmed_reads,
                             error_tolerance, minimum_overlap, allow_indels, parallel)
    assert files_covered == files_to_cover


@pytest.mark.parametrize(
    'fastq_1, fastq_2, outout_folder, three_prime_r1, three_prime_r2,five_prime_r1, five_prime_r2, any_position_r1, '
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
def test_trim_adapters_paired_end_command(monkeypatch, fastq_1, fastq_2, outout_folder, three_prime_r1, three_prime_r2,
                                          five_prime_r1, five_prime_r2, any_position_r1, any_position_r2,
                                          quality_trimming, trim_n, minimum_read_length, maximum_read_length,
                                          discard_untrimmed_reads, pair_filter_if: Literal['any', 'both', 'first'],
                                          error_tolerance, minimum_overlap, allow_indels, parallel, expected_command):
    pairs_to_cover = [('dir1/fq1.fastq', 'dir2/fq4.fq.gz'), ('dir1/fq2.fastq.gz', 'dir2/fq3.fq')]
    pair_stems = [('fq1', 'fq4'), ('fq2', 'fq3')]
    pairs_covered = []

    def mock_run_subprocess(args, print_stdout=True, print_stderr=True):
        for i in range(len(pairs_to_cover)):
            if pairs_to_cover[i][-1] in args[-1]:
                assert args == expected_command + [f'{outout_folder}/{pair_stems[i][0]}_trimmed.fastq.gz',
                                                   '--paired-output',
                                                   f'{outout_folder}/{pair_stems[i][1]}_trimmed.fastq.gz',
                                                   f'tests/test_files/test_fastqs/{pairs_to_cover[i][0]}',
                                                   f'tests/test_files/test_fastqs/{pairs_to_cover[i][1]}', ]
                pairs_covered.append(pairs_to_cover[i])
                break
        assert print_stdout
        assert print_stderr

    monkeypatch.setattr(io, 'run_subprocess', mock_run_subprocess)

    trim_adapters_paired_end(fastq_1, fastq_2, outout_folder, three_prime_r1, three_prime_r2,
                             five_prime_r1, five_prime_r2, any_position_r1, any_position_r2,
                             quality_trimming, trim_n, minimum_read_length, maximum_read_length,
                             discard_untrimmed_reads, pair_filter_if,
                             error_tolerance, minimum_overlap, allow_indels, parallel)
    assert pairs_covered == pairs_to_cover
