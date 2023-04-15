import filecmp
import gzip
import platform

import pytest

from rnalysis import fastq
from rnalysis.fastq import *
from rnalysis.utils import io


def unlink_tree(dir):
    for item in Path(dir).iterdir():
        if 'gitignore' in item.name:
            continue
        if item.is_file():
            item.unlink()
        else:
            shutil.rmtree(item)


def are_dir_trees_equal(dir1, dir2, compare_contents: bool = True):
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
        print(f"mismatch between {dir1} and {dir2} with left_only={dirs_cmp.left_only}, "
              f"right_only={dirs_cmp.right_only}, funny={dirs_cmp.funny_files}")
        return False
    (_, mismatch, errors) = filecmp.cmpfiles(
        dir1, dir2, dirs_cmp.common_files, shallow=False)
    if (len(mismatch) > 0 or len(errors) > 0) and compare_contents:
        print(f"mismatch between {dir1} and {dir2} in the files {mismatch} with errors {errors}")
        for item in mismatch:
            items = []
            for this_dir in [dir1, dir2]:
                pth = Path(this_dir).joinpath(item)
                with open(pth) as f:
                    items.append(f.read())
            if items[0] != items[1]:
                for i in items:
                    print(i)
                    print('---------------------------')
                return False
    for common_dir in dirs_cmp.common_dirs:
        new_dir1 = Path(dir1).joinpath(common_dir).as_posix()
        new_dir2 = Path(dir2).joinpath(common_dir).as_posix()
        if not are_dir_trees_equal(new_dir1, new_dir2, compare_contents):
            return False
    return True


@pytest.mark.parametrize(
    'genome_fastas,output_folder,index_name,bowtie2_installation_folder,random_seed,threads,expected_command', [
        ('tests/test_files/bowtie2_tests/transcripts.fasta', 'tests/test_files/bowtie2_tests/index', 'transcripts',
         'auto', 0, 1,
         ['bowtie2-build-s', '--wrapper', 'basic-0', '--seed', '0', '--threads', '1',
          'tests/test_files/bowtie2_tests/transcripts.fasta',
          'tests/test_files/bowtie2_tests/index/transcripts']),

        (['tests/test_files/bowtie2_tests/transcripts.fasta', 'tests/test_files/bowtie2_tests/transcripts.fasta'],
         'tests/test_files/bowtie2_tests/index', 'transcripts', 'path/to/bt2', 42, 10,
         ['path/to/bt2/bowtie2-build-s', '--wrapper', 'basic-0', '--seed', '42', '--threads', '10',
          'tests/test_files/bowtie2_tests/transcripts.fasta,tests/test_files/bowtie2_tests/transcripts.fasta',
          'tests/test_files/bowtie2_tests/index/transcripts']),
    ])
def test_bowtie2_create_index_command(monkeypatch, genome_fastas, output_folder, index_name,
                                      bowtie2_installation_folder, random_seed, threads, expected_command):
    index_created = []

    def mock_run_subprocess(args, print_stdout=True, print_stderr=True, log_filename: str = None, shell: bool = False):
        assert shell
        if args[-1] == '--version':
            return 0
        assert args == expected_command
        assert print_stdout
        assert print_stderr
        index_created.append(True)

    monkeypatch.setattr(io, 'run_subprocess', mock_run_subprocess)

    bowtie2_create_index(genome_fastas, output_folder, index_name,
                         bowtie2_installation_folder, random_seed, threads)
    assert index_created == [True]


@pytest.mark.parametrize(
    "r1_files,r2_files,output_folder,index_file,bowtie2_installation_folder,"
    "new_sample_names,mode,settings_preset,ignore_qualities,quality_score_type,mate_orientations,"
    "min_fragment_length,max_fragment_length, allow_individual_alignment,allow_disconcordant_alignment,"
    "random_seed,threads,expected_command", [
        (['tests/test_files/kallisto_tests/reads_1.fastq'], ['tests/test_files/kallisto_tests/reads_2.fastq'],
         'tests/test_files/bowtie2_tests/outdir', 'tests/test_files/bowtie2_tests/index/transcripts', 'auto',
         'auto', 'end-to-end', 'very-sensitive', False, 'phred33', 'fwd-rev', 50, 500, True, False, 0, 1,
         ['bowtie2', '--end-to-end', '--very-sensitive', '--phred33', '--seed', '0', '--threads', '1', '-x',
          'tests/test_files/bowtie2_tests/index/transcripts', '-I', '50', '-X', '500', '--no-discorcordant', '--fr',
          '-1', 'tests/test_files/kallisto_tests/reads_1.fastq', '-2',
          'tests/test_files/kallisto_tests/reads_2.fastq', '-S',
          'tests/test_files/bowtie2_tests/outdir/reads_1_reads_2.sam']
         ),
        (['tests/test_files/kallisto_tests/reads_1.fastq'], ['tests/test_files/kallisto_tests/reads_2.fastq'],
         'tests/test_files/bowtie2_tests/outdir', 'tests/test_files/bowtie2_tests/index/transcripts.1.bt2',
         'path/to/bowtie2inst', ['newName'], 'local', 'fast', True, 'phred64', 'fwd-fwd', 0, 250, False, True, 42, 12,
         ['path/to/bowtie2inst/bowtie2', '--local', '--fast', '--phred64', '--ignore-quals', '--seed', '42',
          '--threads', '12', '-x', 'tests/test_files/bowtie2_tests/index/transcripts', '-I', '0', '-X', '250',
          '--no-mixed', '--ff', '-1', 'tests/test_files/kallisto_tests/reads_1.fastq', '-2',
          'tests/test_files/kallisto_tests/reads_2.fastq', '-S', 'tests/test_files/bowtie2_tests/outdir/newName.sam']
         ),
    ])
def test_bowtie2_align_paired_end_command(monkeypatch, r1_files, r2_files, output_folder, index_file,
                                          bowtie2_installation_folder, new_sample_names, mode, settings_preset,
                                          ignore_qualities, quality_score_type,
                                          mate_orientations, min_fragment_length, max_fragment_length,
                                          allow_individual_alignment,
                                          allow_disconcordant_alignment, random_seed, threads, expected_command):
    pairs_to_cover = [('reads_1.fastq', 'reads_2.fastq')]
    pairs_covered = []

    def mock_run_subprocess(args, print_stdout=True, print_stderr=True, log_filename: str = None, shell: bool = False):
        assert shell
        if args[1] == '--version':
            return 0

        assert args == expected_command
        pairs_covered.append(pairs_to_cover[0])
        assert print_stdout
        assert print_stderr

    monkeypatch.setattr(io, 'run_subprocess', mock_run_subprocess)

    bowtie2_align_paired_end(r1_files, r2_files, output_folder, index_file, bowtie2_installation_folder,
                             new_sample_names, mode, settings_preset, ignore_qualities, quality_score_type,
                             mate_orientations, min_fragment_length, max_fragment_length, allow_individual_alignment,
                             allow_disconcordant_alignment, random_seed, threads)
    assert sorted(pairs_covered) == sorted(pairs_to_cover)


@pytest.mark.parametrize(
    "fastq_folder,output_folder,index_file,bowtie2_installation_folder,new_sample_names,mode,settings_preset,"
    "ignore_qualities,quality_score_type,random_seed,threads,expected_command", [
        ('tests/test_files/kallisto_tests',
         'tests/test_files/bowtie2_tests/outdir', 'tests/test_files/bowtie2_tests/index/transcripts', 'auto',
         'auto', 'end-to-end', 'very-sensitive', False, 'phred33', 0, 1,
         ['bowtie2', '--end-to-end', '--very-sensitive', '--phred33', '--seed', '0', '--threads', '1', '-x',
          'tests/test_files/bowtie2_tests/index/transcripts', '-U', '', '-S', '']
         ),
        ('tests/test_files/kallisto_tests',
         'tests/test_files/bowtie2_tests/outdir', 'tests/test_files/bowtie2_tests/index/transcripts.1.bt2',
         'path/to/bowtie2inst', ['newName1', 'newName2', ], 'local', 'fast', True, 'phred64', 42, 12,
         ['path/to/bowtie2inst/bowtie2', '--local', '--fast', '--phred64', '--ignore-quals', '--seed', '42',
          '--threads', '12', '-x', 'tests/test_files/bowtie2_tests/index/transcripts', '-U', '', '-S', '']
         ),
    ])
def test_bowtie2_align_single_end_command(monkeypatch, fastq_folder, output_folder, index_file,
                                          bowtie2_installation_folder, new_sample_names, mode, settings_preset,
                                          ignore_qualities, quality_score_type, random_seed, threads, expected_command):
    files_to_cover = ['reads_1.fastq', 'reads_2.fastq']
    file_stems = ['reads_1', 'reads_2']
    files_covered = []

    def mock_run_subprocess(args, print_stdout=True, print_stderr=True, log_filename: str = None, shell: bool = False):
        assert shell
        if args[1] == '--version':
            return 0
        for i in range(len(files_to_cover)):
            if files_to_cover[i] in args[-3]:
                exp = expected_command.copy()
                exp[-3] = fastq_folder + '/' + files_to_cover[i]
                if new_sample_names == 'auto':
                    exp[-1] = f'{output_folder}/{file_stems[i]}.sam'
                else:
                    exp[-1] = f'{output_folder}/{new_sample_names[i]}.sam'

                assert args == exp
                files_covered.append(files_to_cover[i])
                break
        assert print_stdout
        assert print_stderr

    monkeypatch.setattr(io, 'run_subprocess', mock_run_subprocess)

    bowtie2_align_single_end(fastq_folder, output_folder, index_file,
                             bowtie2_installation_folder, new_sample_names, mode, settings_preset,
                             ignore_qualities, quality_score_type, random_seed, threads)
    assert sorted(files_covered) == sorted(files_to_cover)


@pytest.mark.parametrize(
    "fastq_folder,output_folder,genome_fasta,shortstack_installation_folder,new_sample_names,known_rnas,trim_adapter,"
    "autotrim_key,multimap_mode,align_only,show_secondary,dicer_min_length,dicer_max_length,loci_file,locus,"
    "search_microrna,strand_cutoff,min_coverage,pad,threads,expected_command", [
        ('tests/test_files/kallisto_tests',
         'tests/test_files/shortstack_tests/outdir', 'tests/test_files/shortstack_tests/transcripts.fasta',
         'auto',
         'auto', None, None, 'autokey', 'fractional', False, True, 10, 30, None, 'locus string', None, 0.8, 2, 75, 1,
         ['ShortStack', '--genomefile', 'tests/test_files/shortstack_tests/transcripts.fasta', '--mmap', 'f',
          '--show_secondaries', '--locus', 'locus string', '--nohp', '--dicermin', '10', '--dicermax', '30',
          '--strand_cutoff', '0.8', '--mincov', '2', '--pad', '75', '--threads', '1', '--readfile',
          'tests/test_files/kallisto_tests/reads_1.fastq', '--outdir',
          'tests/test_files/shortstack_tests/outdir/reads_1']

         ),
        ('tests/test_files/kallisto_tests',
         'tests/test_files/shortstack_tests/outdir', 'tests/test_files/shortstack_tests/transcripts.fasta',
         'path/to/shortstackinst', ['newName1', 'newName2', ], 'tests/test_files/test_deseq.csv',
         'autotrim', 'autokey', 'random', True, False, 21, 24, 'tests/test_files/counted.csv',
         None, 'known-rnas', 0.9, 1.2, 3, 12,
         ['path/to/shortstackinst/ShortStack', '--genomefile', 'tests/test_files/shortstack_tests/transcripts.fasta',
          '--mmap', 'r', '--adapter', 'autotrim', '--align_only', '--locifile', 'tests/test_files/counted.csv',
          '--knownRNAs', 'tests/test_files/test_deseq.csv',
          '--dicermin', '21', '--dicermax', '24', '--strand_cutoff',
          '0.9', '--mincov', '1.2', '--pad', '3', '--threads', '12', '--readfile',
          'tests/test_files/kallisto_tests/reads_1.fastq', '--outdir',
          'tests/test_files/shortstack_tests/outdir/newName1']

         ),
    ])
def test_shortstack_command(monkeypatch, fastq_folder, output_folder, genome_fasta, shortstack_installation_folder,
                            new_sample_names, known_rnas, trim_adapter, autotrim_key, multimap_mode, align_only,
                            show_secondary, dicer_min_length, dicer_max_length, loci_file, locus, search_microrna,
                            strand_cutoff, min_coverage, pad, threads, expected_command):
    files_to_cover = ['reads_1.fastq', 'reads_2.fastq']
    file_stems = ['reads_1', 'reads_2']
    files_covered = []

    def mock_run_subprocess(args, print_stdout=True, print_stderr=True, log_filename: str = None, shell: bool = False):
        assert shell
        if args[1] == '--version':
            return 0
        for i in range(len(files_to_cover)):
            if files_to_cover[i] in args[-3]:
                exp = expected_command.copy()
                exp[-3] = fastq_folder + '/' + files_to_cover[i]
                if new_sample_names == 'auto':
                    exp[-1] = f'{output_folder}/{file_stems[i]}'
                else:
                    exp[-1] = f'{output_folder}/{new_sample_names[i]}'

                assert args == exp
                files_covered.append(files_to_cover[i])
                break
        assert print_stdout
        assert print_stderr

    monkeypatch.setattr(io, 'run_subprocess', mock_run_subprocess)

    shortstack_align_smallrna(fastq_folder, output_folder, genome_fasta, shortstack_installation_folder,
                              new_sample_names, known_rnas, trim_adapter, autotrim_key, multimap_mode, align_only,
                              show_secondary, dicer_min_length, dicer_max_length, loci_file, locus, search_microrna,
                              strand_cutoff, min_coverage, pad, threads)
    assert sorted(files_covered) == sorted(files_to_cover)


def test_bowtie2_create_index():
    out_path = 'tests/test_files/bowtie2_tests/outdir'
    truth_path = 'tests/test_files/bowtie2_tests/index'
    try:
        bowtie2_create_index(['tests/test_files/bowtie2_tests/transcripts.fasta'], out_path, random_seed=0)
        assert are_dir_trees_equal(out_path, truth_path, compare_contents=False)
    except Exception as e:
        if platform.system() == 'Windows':
            pytest.xfail('bowtie2 PATH is not defined properly on GitHub Actions Windows')
        raise e
    finally:
        unlink_tree(out_path)


def test_bowtie2_align_single_end():
    in_dir = 'tests/test_files/kallisto_tests'
    index_file = 'tests/test_files/bowtie2_tests/index/transcripts.1.bt2'
    out_dir = 'tests/test_files/bowtie2_tests/outdir'
    truth_dir = 'tests/test_files/bowtie2_tests/truth/single'
    try:
        bowtie2_align_single_end(in_dir, out_dir, index_file)
        assert are_dir_trees_equal(out_dir, truth_dir, compare_contents=False)
    finally:
        unlink_tree(out_dir)


def test_bowtie2_align_paired_end():
    in_1 = 'tests/test_files/kallisto_tests/reads_1.fastq'
    in_2 = 'tests/test_files/kallisto_tests/reads_2.fastq'
    index_file = 'tests/test_files/bowtie2_tests/index/transcripts.1.bt2'
    out_dir = 'tests/test_files/bowtie2_tests/outdir'
    truth_dir = 'tests/test_files/bowtie2_tests/truth/paired'
    try:
        bowtie2_align_paired_end([in_1], [in_2], out_dir, index_file, random_seed=42)
        assert are_dir_trees_equal(out_dir, truth_dir, compare_contents=False)
    finally:
        unlink_tree(out_dir)


@pytest.mark.parametrize("transcriptome_fasta,kallisto_installation_folder,kmer_length,make_unique,expected_command", [
    ('tests/test_files/kallisto_tests/transcripts.fasta', 'auto', 5, True,
     ['kallisto', 'index', '-i', 'tests/test_files/kallisto_tests/transcripts.idx', '-k', '5',
      '--make-unique', 'tests/test_files/kallisto_tests/transcripts.fasta']),
    ('tests/test_files/kallisto_tests/transcripts.fasta', 'pth/to/kallisto', 3, False,
     ['pth/to/kallisto/kallisto', 'index', '-i', 'tests/test_files/kallisto_tests/transcripts.idx', '-k', '3',
      'tests/test_files/kallisto_tests/transcripts.fasta']),
])
def test_kallisto_create_index_command(monkeypatch, transcriptome_fasta, kallisto_installation_folder, kmer_length,
                                       make_unique, expected_command):
    index_created = []

    def mock_run_subprocess(args, print_stdout=True, print_stderr=True, log_filename: str = None, shell: bool = False):
        assert not shell
        if args[-1] == 'version':
            return 0
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

    log_path = 'tests/test_files/kallisto_tests/kallisto-index_transcripts.log'
    log_truth_path = 'tests/test_files/kallisto_tests/kallisto-index_transcripts_truth.log'
    try:
        kallisto_create_index('tests/test_files/kallisto_tests/transcripts.fasta')
        with open(truth_path, 'rb') as truth, open(out_path, 'rb') as out:
            assert truth.read() == out.read()
        with open(log_truth_path) as truth, open(log_path) as out:
            assert truth.read().replace('\n', '') == out.read().replace('\n', '')
    finally:
        if Path(out_path).exists():
            Path(out_path).unlink()
        if Path(log_path).exists():
            Path(log_path).unlink()


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
        unlink_tree(out_dir)


def test_kallisto_quantify_paired_end():
    in1_path = 'tests/test_files/kallisto_tests/reads_1.fastq'
    in2_path = 'tests/test_files/kallisto_tests/reads_2.fastq'
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
        unlink_tree(out_dir)


@pytest.mark.parametrize(
    "fastq_folder,output_folder,index_file,gtf_file,average_fragment_length,stdev_fragment_length,"
    "kallisto_installation_folder,"
    "new_sample_names,stranded,learn_bias,seek_fusion_genes,bootstrap_samples,expected_command", [
        ('tests/test_files/kallisto_tests', 'tests/test_files/kallisto_tests/outdir',
         'tests/test_files/kallisto_tests/transcripts_truth.idx', 'tests/test_files/kallisto_tests/transcripts.gtf',
         125, 14, 'auto', 'auto', 'no', False, False, None,
         ['kallisto', 'quant', '-i', 'tests/test_files/kallisto_tests/transcripts_truth.idx',
          '-o', 'outfolder', '--single', '-s', '14', '-l', '125']),
        ('tests/test_files/kallisto_tests', 'tests/test_files/kallisto_tests/outdir',
         'tests/test_files/kallisto_tests/transcripts_truth.idx', 'tests/test_files/kallisto_tests/transcripts.gtf',
         8.5, 0.2, 'inst/folder', ['new_name_1', 'new_name_2'], 'reverse', True, True, 3,
         ['inst/folder/kallisto', 'quant', '-i', 'tests/test_files/kallisto_tests/transcripts_truth.idx',
          '--bias', '--fusion', '--rf-stranded', '-b', '3', '-o', 'outfolder', '--single', '-s', '0.2', '-l', '8.5']),
    ])
def test_kallisto_quantify_single_end_command(monkeypatch, fastq_folder, output_folder, index_file, gtf_file,
                                              average_fragment_length, stdev_fragment_length,
                                              kallisto_installation_folder, new_sample_names, stranded, learn_bias,
                                              seek_fusion_genes, bootstrap_samples, expected_command):
    files_to_cover = ['reads_1.fastq', 'reads_2.fastq']
    file_stems = ['reads_1', 'reads_2']
    files_covered = []
    output_processed = []

    def mock_run_subprocess(args, print_stdout=True, print_stderr=True, log_filename: str = None, shell: bool = False):
        assert not shell
        if args[1] == 'version':
            return 0
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
        (['tests/test_files/kallisto_tests/reads_1.fastq'], ['tests/test_files/kallisto_tests/reads_2.fastq'],
         'tests/test_files/kallisto_tests/outdir',
         'tests/test_files/kallisto_tests/transcripts_truth.idx', 'tests/test_files/kallisto_tests/transcripts.gtf',
         'auto', 'auto', 'no', False, False, None,
         ['kallisto', 'quant', '-i', 'tests/test_files/kallisto_tests/transcripts_truth.idx',
          '-o', ]),
        (['tests/test_files/kallisto_tests/reads_1.fastq'], ['tests/test_files/kallisto_tests/reads_2.fastq'],
         'tests/test_files/kallisto_tests/outdir',
         'tests/test_files/kallisto_tests/transcripts_truth.idx', 'tests/test_files/kallisto_tests/transcripts.gtf',
         'kallisto/inst/folder', ['new_pair_name'], 'reverse', True, True, 3,
         ['kallisto/inst/folder/kallisto', 'quant', '-i', 'tests/test_files/kallisto_tests/transcripts_truth.idx',
          '--bias', '--fusion', '--rf-stranded', '-b', '3', '-o', ]),
    ])
def test_kallisto_quantify_paired_end_command(monkeypatch, r1_files, r2_files, output_folder, index_file, gtf_file,
                                              kallisto_installation_folder, new_sample_names, stranded, learn_bias,
                                              seek_fusion_genes, bootstrap_samples, expected_command):
    pairs_to_cover = [('reads_1.fastq', 'reads_2.fastq'), ]
    pair_stems = [('reads_1', 'reads_2'), ]
    pairs_covered = []
    output_processed = []

    def mock_run_subprocess(args, print_stdout=True, print_stderr=True, log_filename: str = None, shell: bool = False):
        assert not shell
        if args[1] == 'version':
            return 0

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
        trim_adapters_single_end(in_dir, out_dir, adapter_seq, gzip_output=True)
        with gzip.open(truth_path) as truth, gzip.open(out_path) as out:
            assert truth.read() == out.read()
    finally:
        if Path(out_path).exists():
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
        trim_adapters_paired_end([in1_path], [in2_path], out_dir, adapter1_seq, adapter2_seq, minimum_read_length=5,
                                 gzip_output=True)

        with gzip.open(truth1_path) as truth, gzip.open(out1_path) as out:
            assert truth.read() == out.read()
        with gzip.open(truth2_path) as truth, gzip.open(out2_path) as out:
            assert truth.read() == out.read()
    finally:
        for pth in [out1_path, out2_path]:
            if Path(pth).exists():
                Path(pth).unlink()
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
    monkeypatch.setattr(io, 'generate_base_call', lambda *args, **kwargs: ['cutadapt'])

    trim_adapters_single_end(fastq_folder, output_folder, three_prime_adapters,
                             five_prime_adapters, any_position_adapters, quality_trimming, trim_n,
                             minimum_read_length, maximum_read_length, discard_untrimmed_reads,
                             error_tolerance, minimum_overlap, allow_indels, parallel, gzip_output=True)
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
    monkeypatch.setattr(io, 'generate_base_call', lambda *args, **kwargs: ['cutadapt'])

    trim_adapters_paired_end(fastq_1, fastq_2, output_folder, three_prime_r1, three_prime_r2,
                             five_prime_r1, five_prime_r2, any_position_r1, any_position_r2,
                             quality_trimming, trim_n, minimum_read_length, maximum_read_length,
                             discard_untrimmed_reads, pair_filter_if,
                             error_tolerance, minimum_overlap, allow_indels, parallel, gzip_output=True)
    assert sorted(pairs_covered) == sorted(pairs_to_cover)


def test_featurecounts_single_end():
    counts_truth = filtering.CountFilter('tests/test_files/featurecounts_tests/truth/single/featureCounts_counts.csv')
    annotations_truth = io.load_csv('tests/test_files/featurecounts_tests/truth/single/featureCounts_annotation.csv', 0)
    stats_truth = io.load_csv('tests/test_files/featurecounts_tests/truth/single/featureCounts_stats.csv', 0)
    truth_outdir = 'tests/test_files/featurecounts_tests/truth/single'
    outdir = 'tests/test_files/featurecounts_tests/outdir'
    gtf_file = 'tests/test_files/featurecounts_tests/single/bamfile_no_qualities.gtf'
    new_sample_names = ['sample1_new']
    try:
        counts, annotations, stats = featurecounts_single_end('tests/test_files/featurecounts_tests/single', outdir,
                                                              gtf_file, new_sample_names=new_sample_names)
        assert counts == counts_truth
        assert annotations.equals(annotations_truth)
        assert stats.equals(stats_truth)
        assert are_dir_trees_equal(outdir, truth_outdir)
    finally:
        unlink_tree(outdir)


def test_featurecounts_paired_end():
    counts_truth = filtering.CountFilter('tests/test_files/featurecounts_tests/truth/paired/featureCounts_counts.csv')
    annotations_truth = io.load_csv('tests/test_files/featurecounts_tests/truth/paired/featureCounts_annotation.csv', 0)
    stats_truth = io.load_csv('tests/test_files/featurecounts_tests/truth/paired/featureCounts_stats.csv', 0)
    truth_outdir = 'tests/test_files/featurecounts_tests/truth/paired'
    outdir = 'tests/test_files/featurecounts_tests/outdir'
    try:
        counts, annotations, stats = featurecounts_paired_end('tests/test_files/featurecounts_tests/paired', outdir,
                                                              'tests/test_files/featurecounts_tests/test-minimum.gtf')
        assert counts == counts_truth
        assert annotations.equals(annotations_truth)
        assert stats.equals(stats_truth)
        assert are_dir_trees_equal(outdir, truth_outdir)
    finally:
        unlink_tree(outdir)


def test_SingleEndPipeline_add_function():
    param_dict1 = {'allow_indels': True, 'any_position_adapters': None, 'discard_untrimmed_reads': False,
                   'error_tolerance': 0.1}
    p = SingleEndPipeline()
    assert p.functions == [] and p.params == []
    p.add_function(trim_adapters_single_end, 'ATGGG',
                   **param_dict1)

    assert p.functions == [trim_adapters_single_end] and p.params == [(('ATGGG',), param_dict1)]

    p.add_function('bowtie2_align_single_end', 'arg', 'arg', kw1='val', kw2='val2')
    assert p.functions == [trim_adapters_single_end, bowtie2_align_single_end] and \
           p.params == [(('ATGGG',), param_dict1), (('arg', 'arg'), {'kw1': 'val', 'kw2': 'val2'})]


def test_SingleEndPipeline_import():
    truth = SingleEndPipeline()
    truth.add_function(trim_adapters_single_end,
                       **{'allow_indels': True, 'any_position_adapters': None, 'discard_untrimmed_reads': False,
                          'error_tolerance': 0.1, 'five_prime_adapters': None, 'maximum_read_length': None,
                          'minimum_overlap': 3, 'minimum_read_length': 10, 'parallel': True, 'quality_trimming': 20,
                          'three_prime_adapters': 'ATGGGTATATGGGT', 'trim_n': True, 'gzip_output': False})
    truth.add_function(bowtie2_align_single_end, **{'bowtie2_installation_folder': 'auto', 'ignore_qualities': False,
                                                    'index_file': 'tests/test_files/bowtie2_tests/index/transcripts.1.bt2',
                                                    'mode': 'end-to-end', 'new_sample_names': 'auto',
                                                    'quality_score_type': 'phred33', 'random_seed': 0,
                                                    'settings_preset': 'very-sensitive', 'threads': 1})
    truth.add_function(featurecounts_single_end, **{'count_fractionally': False, 'count_multi_mapping_reads': False,
                                                    'count_multi_overlapping_reads': False, 'gtf_attr_name': 'gene_id',
                                                    'gtf_feature_type': 'exon',
                                                    'gtf_file': 'tests/test_files/kallisto_tests/transcripts.gtf',
                                                    'ignore_secondary': True, 'is_long_read': False,
                                                    'min_mapping_quality': 0, 'new_sample_names': 'auto',
                                                    'r_installation_folder': 'auto', 'report_read_assignment': None,
                                                    'stranded': 'no', 'threads': 1})

    pth = 'tests/test_files/test_singleend_pipeline.yaml'
    p = SingleEndPipeline.import_pipeline(pth)
    assert p == truth


def test_SingleEndPipeline_apply_to():
    in_dir = 'tests/test_files/fastq_pipeline_tests/in'
    out_dir = Path('tests/test_files/fastq_pipeline_tests/single/outdir')
    truth_dir = Path('tests/test_files/fastq_pipeline_tests/single/truth')
    pth = 'tests/test_files/test_singleend_pipeline.yaml'
    p = SingleEndPipeline.import_pipeline(pth)

    try:
        p.apply_to(in_dir, out_dir)
        for file in Path(out_dir).joinpath('01_trim_adapters_single_end').iterdir():
            if file.is_file() and file.suffix == '.log':
                file.unlink()

        assert are_dir_trees_equal(out_dir.joinpath('01_trim_adapters_single_end'),
                                   truth_dir.joinpath('01_trim_adapters_single_end'))
        assert are_dir_trees_equal(out_dir.joinpath('02_bowtie2_align_single_end'),
                                   truth_dir.joinpath('02_bowtie2_align_single_end'),compare_contents=False)
        assert are_dir_trees_equal(out_dir.joinpath('03_featurecounts_single_end'),
                                   truth_dir.joinpath('03_featurecounts_single_end'))
    finally:
        unlink_tree(out_dir)
