"""
The *fastq* module provides a unified programmatic interface to external tools that process FASTQ files.
Those currently include the *CutAdapt* adapter-trimming tool, the *kallisto* RNA-sequencing quantification tool,
the *bowtie2* alignment tool, and the *featureCounts* feature counting tool.
"""

import abc
import itertools
import os
import shutil
import sys
import types
import typing
import warnings
from pathlib import Path
from typing import Union, List, Tuple, Literal

import pandas as pd
from tqdm.auto import tqdm

from rnalysis import filtering
from rnalysis.utils import validation, parsing, io, feature_counting, genome_annotation, generic, installs
from rnalysis.utils.generic import readable_name
from rnalysis.utils.param_typing import PositiveInt, NonNegativeInt, Fraction, LEGAL_FASTQ_SUFFIXES, \
    LEGAL_BOWTIE2_PRESETS, LEGAL_BOWTIE2_MODES, LEGAL_QUAL_SCORE_TYPES, LEGAL_ALIGNMENT_SUFFIXES

try:
    from cutadapt.cli import main as cutadapt_main

    HAS_CUTADAPT = True

except ImportError:  # pragma: no cover
    HAS_CUTADAPT = False


    def cutadapt_main(x):
        return None


def _func_type(func_type: Literal['single', 'paired', 'both']):
    def decorator(item):
        item.func_type = func_type
        return item

    return decorator


class _FASTQPipeline(generic.GenericPipeline, abc.ABC):
    def __str__(self):
        string = ''
        if len(self) > 0:
            string += ":\n\t" + '\n\t'.join(
                self._readable_func_signature(func, params[0], params[1]) for func, params in
                zip(self.functions, self.params))
        return string

    def __repr__(self):
        string = ''
        if len(self) > 0:
            string += ": " + "-->".join(
                self._func_signature(func, params[0], params[1]) for func, params in zip(self.functions, self.params))
        return string

    @property
    def is_paired_end(self):
        raise NotImplementedError

    @staticmethod
    def _is_paired_end_func(func: typing.Callable):
        if func.__name__.endswith('paired_end'):
            return True
        elif func.__name__.endswith('single_end'):
            return False
        raise ValueError(f"Cannot determine whether function '{func}' is single-end or paired-end.")

    def add_function(self, func: Union[str, types.FunctionType], *args, **kwargs):
        if isinstance(func, str):
            thismodule = sys.modules[__name__]
            func = getattr(thismodule, func)

        assert self._is_paired_end_func(func) == self.is_paired_end, \
            f"{'paired' * (not self.is_paired_end) + 'single' * self.is_paired_end}-end function " \
            f"cannot be added to {'single' * (not self.is_paired_end) + 'paired' * self.is_paired_end}-end Pipeline!"

        super().add_function(func, *args, **kwargs)

    def _init_from_dict(self, pipeline_dict: dict):
        self.params = [(parsing.data_to_tuple(p[0]), p[1]) for p in pipeline_dict['params']]
        self.functions = [getattr(sys.modules[__name__], func) for func in pipeline_dict['functions']]


class SingleEndPipeline(_FASTQPipeline):
    def __str__(self):
        return "Pipeline for sequence files (single-end)" + super().__str__()

    def __repr__(self):
        return "SingleEndPipeline()" + super().__repr__()

    def _get_pipeline_dict(self):
        d = super()._get_pipeline_dict()
        d['metadata']['pipeline_type'] = 'single'
        return d

    @readable_name('Apply Pipeline to single-end sequencing data')
    def apply_to(self, input_folder: Union[str, Path], output_folder: Union[str, Path]):
        input_folder = Path(input_folder)
        output_folder = Path(output_folder)
        return_values = []

        assert input_folder.exists() and input_folder.is_dir(), "input_folder does not exist!"
        assert output_folder.exists() and output_folder.is_dir(), "output_folder does not exist!"

        current_in_dir = input_folder
        for i, (func, (args, kwargs)) in enumerate(zip(self.functions, self.params)):
            params = generic.get_method_signature(func)
            current_out_dir = output_folder.joinpath(f'{i + 1:02d}_{func.__name__}')
            try:
                current_out_dir.mkdir(parents=True)
            except OSError:
                pass

            res = func(current_in_dir, current_out_dir, *args, **kwargs)
            if res is not None:  # TODO: make uniform with Pipeline output dict?
                if isinstance(res, (list, tuple)):
                    return_values.extend(res)
                else:
                    return_values.append(res)

            if 'new_sample_names' in params:
                current_in_dir = current_out_dir

        return parsing.data_to_tuple(return_values)

    @property
    def is_paired_end(self):
        return False


class PairedEndPipeline(_FASTQPipeline):
    def __str__(self):
        return "Pipeline for sequence files (paired-end)" + super().__str__()

    def __repr__(self):
        return "PairedEndPipeline()" + super().__repr__()

    def _get_pipeline_dict(self):
        d = super()._get_pipeline_dict()
        d['metadata']['pipeline_type'] = 'paired'
        return d

    @property
    def is_paired_end(self):
        return True

    @readable_name('Apply Pipeline to paired-end sequencing data')
    def apply_to(self, r1_files: List[str], r2_files: List[str], output_folder: Union[str, Path]):
        output_folder = Path(output_folder)
        return_values = []

        assert output_folder.exists() and output_folder.is_dir(), "output_folder does not exist!"

        current_r1, current_r2 = r1_files, r2_files
        current_in_dir = output_folder.joinpath('00_input')
        # if the first function in the Pipeline takes a single input folder as input,
        # create an 'input' folder and copy the starting files into it.
        params = generic.get_method_signature(self.functions[0])
        if 'input_folder' in params:
            try:
                current_in_dir.mkdir(parents=True)
            except OSError:
                pass
            with tqdm(total=len(r1_files) + len(r2_files), desc='Copying files', unit='files') as pbar:
                for file in (itertools.chain(r1_files, r2_files)):
                    shutil.copy(file, current_in_dir)
                    pbar.update(1)

        for i, (func, (args, kwargs)) in enumerate(zip(self.functions, self.params)):
            current_out_dir = output_folder.joinpath(f'{i + 1:02d}_{func.__name__}')
            try:
                current_out_dir.mkdir(parents=True)
            except OSError:
                pass

            params = generic.get_method_signature(func)

            if 'return_new_filenames' in params:
                kwargs['return_new_filenames'] = True

            if 'r1_files' in params and 'r2_files' in params:
                res = func(current_r1, current_r2, current_out_dir, *args, **kwargs)
            elif 'input_folder' in params:
                res = func(current_in_dir, current_out_dir, *args, **kwargs)
            else:
                raise NotImplementedError(f"Cannot run function '{func.__name__}' in this Pipeline.")

            if res is not None and not (validation.isiterable(res) and validation.isinstanceiter(res,
                                                                                                 str)):  # TODO: make uniform with Pipeline output dict?
                if isinstance(res, (list, tuple)):
                    return_values.extend(res)
                else:
                    return_values.append(res)

            if 'new_sample_names' in params:
                current_in_dir = current_out_dir
            if 'return_new_filenames' in params:
                current_r1, current_r2 = res

        return parsing.data_to_tuple(return_values)


def _generate_picard_basecall(command: str, picard_installation_folder: Union[str, Path, Literal['auto']]) -> List[str]:
    if picard_installation_folder == 'auto':
        installs.install_picard()
        picard_installation_folder = installs.PICARD_JAR
    args = ['-jar', picard_installation_folder, command]
    base_call = io.generate_base_call('java', installs.JDK_PATH, args=args)
    return base_call


def _run_picard_calls(calls: List[List[str]], script_name: str, output_folder: Path):
    with tqdm(total=len(calls), desc=f'Applying "{script_name}"', unit='files') as pbar:
        for picard_call in calls:
            print(f"Running command: \n{' '.join(picard_call)}")
            log_filename = Path(output_folder).joinpath(
                f'picard_{script_name}_{Path(picard_call[-1]).stem}.log').absolute().as_posix()
            io.run_subprocess(picard_call, shell=True, log_filename=log_filename)
            print(f"File saved successfully at {picard_call[-1]}")
            pbar.update(1)


def _parse_sam2fastq_misc_args(picard_installation_folder: Union[str, Path, Literal['auto']] = 'auto',
                               re_reverse_reads: bool = True, include_non_primary_alignments: bool = False,
                               quality_trim: Union[PositiveInt, None] = None):
    script_name = 'SamToFastq'
    base_call = _generate_picard_basecall(script_name, picard_installation_folder)
    base_call.append(f"RE_REVERSE={'true' if re_reverse_reads else 'false'}")
    base_call.append(f"INCLUDE_NON_PRIMARY_ALIGNMENTS ={'true' if include_non_primary_alignments else 'false'}")
    if quality_trim is not None:
        base_call.append(f"QUALITY={quality_trim}")
    return base_call, script_name


@_func_type('single')
@readable_name('Convert SAM/BAM files to FASTQ files (single-end)')
def sam_to_fastq_single(input_folder: Union[str, Path], output_folder: Union[str, Path],
                        picard_installation_folder: Union[str, Path, Literal['auto']] = 'auto',
                        new_sample_names: Union[List[str], Literal['auto']] = 'auto',
                        re_reverse_reads: bool = True, include_non_primary_alignments: bool = False,
                        quality_trim: Union[PositiveInt, None] = None):
    """
    Convert SAM/BAM files to FASTQ files using \
    `Picard SamToFastq <https://broadinstitute.github.io/picard/command-line-overview.html#SamToFastq>`_.

    :param input_folder: Path to the folder containing the SAM/BAM files you want to convert.
    :type input_folder: str or Path
    :param output_folder: Path to a folder in which the converted FASTQ files will be saved.
    :type output_folder: str or Path
    :param picard_installation_folder: Path to the installation folder of Picard. For example: \
    'C:/Program Files/Picard'
    :type picard_installation_folder: str, Path, or 'auto' (default='auto')
    :param new_sample_names: Give a new name to each converted sample (optional). \
    If sample_names='auto', sample names \
    will be given automatically. Otherwise, sample_names should be a list of new names, with the order of the names \
    matching the order of the files in the directory.
    :type new_sample_names: list of str or 'auto' (default='auto')
    :param re_reverse_reads: Re-reverse bases and qualities of reads with the negative-strand flag \
    before writing them to FASTQ.
    :type re_reverse_reads: bool (default=True)
    :param include_non_primary_alignments: If true, include non-primary alignments in the output. \
    Support of non-primary alignments in SamToFastq is not comprehensive, \
    so there may be exceptions if this is set to true and there are paired reads with non-primary alignments.
    :type include_non_primary_alignments: bool (default=False)
    :param quality_trim: If enabled, End-trim reads using the phred/bwa quality trimming algorithm and this quality.
    :type quality_trim: positive int or None (default=None)
    :return: a list of the paths to the generated FASTQ files.
    :rtype: list of str
    """
    input_folder = Path(input_folder)
    output_folder = Path(output_folder)
    assert input_folder.exists() and input_folder.is_dir(), "input_folder does not exist!"
    assert output_folder.exists() and output_folder.is_dir(), "output_folder does not exist!"

    base_call, script_name = _parse_sam2fastq_misc_args(picard_installation_folder, re_reverse_reads,
                                                        include_non_primary_alignments, quality_trim)
    calls = []

    legal_samples = _get_legal_samples(input_folder, 'alignment')
    assert (new_sample_names == 'auto') or (len(new_sample_names) == len(legal_samples)), \
        f'Number of samples ({len(legal_samples)}) does not match number of sample names ({len(new_sample_names)})!'

    for i, sam_file in enumerate(sorted(legal_samples)):
        this_call = base_call.copy()
        if new_sample_names == 'auto':
            this_name = parsing.remove_suffixes(sam_file).stem + '_sam2fastq'
        else:
            this_name = new_sample_names[i]
        this_call.append(f"INPUT={sam_file.as_posix()}")
        this_call.append(f"FASTQ={output_folder.joinpath(f'{this_name}.fastq').as_posix()}")
        calls.append(this_call)

    _run_picard_calls(calls, script_name, output_folder)


@_func_type('paired')
@readable_name('Convert SAM/BAM files to FASTQ files (paired-end)')
def sam_to_fastq_paired(input_folder: Union[str, Path], output_folder: Union[str, Path],
                        picard_installation_folder: Union[str, Path, Literal['auto']] = 'auto',
                        new_sample_names: Union[List[str], Literal['auto']] = 'auto',
                        re_reverse_reads: bool = True, include_non_primary_alignments: bool = False,
                        quality_trim: Union[PositiveInt, None] = None, return_new_filenames: bool = False):
    """
    Convert SAM/BAM files to FASTQ files using \
    `Picard SamToFastq <https://broadinstitute.github.io/picard/command-line-overview.html#SamToFastq>`_.

    :param input_folder: Path to the folder containing the SAM/BAM files you want to convert.
    :type input_folder: str or Path
    :param output_folder: Path to a folder in which the converted FASTQ files will be saved.
    :type output_folder: str or Path
    :param picard_installation_folder: Path to the installation folder of Picard. For example: \
    'C:/Program Files/Picard'
    :type picard_installation_folder: str, Path, or 'auto' (default='auto')
    :param new_sample_names: Give a new name to each converted sample (optional). \
    If sample_names='auto', sample names \
    will be given automatically. Otherwise, sample_names should be a list of new names, with the order of the names \
    matching the order of the files in the directory.
    :type new_sample_names: list of str or 'auto' (default='auto')
    before writing them to FASTQ.
    :type re_reverse_reads: bool (default=True)
    :param include_non_primary_alignments: If true, include non-primary alignments in the output. \
    Support of non-primary alignments in SamToFastq is not comprehensive, \
    so there may be exceptions if this is set to true and there are paired reads with non-primary alignments.
    :type include_non_primary_alignments: bool (default=False)
    :param quality_trim: If enabled, End-trim reads using the phred/bwa quality trimming algorithm and this quality.
    :type quality_trim: positive int or None (default=None)
    :return: a list of the paths to the generated FASTQ files.
    :return: a list of the paths to the generated FASTQ files.
    :rtype: list of str
    """
    input_folder = Path(input_folder)
    output_folder = Path(output_folder)
    assert input_folder.exists() and input_folder.is_dir(), "input_folder does not exist!"
    assert output_folder.exists() and output_folder.is_dir(), "output_folder does not exist!"

    base_call, script_name = _parse_sam2fastq_misc_args(picard_installation_folder, re_reverse_reads,
                                                        include_non_primary_alignments, quality_trim)

    calls = []

    legal_samples = _get_legal_samples(input_folder, 'alignment')
    assert (new_sample_names == 'auto') or (len(new_sample_names) == len(legal_samples)), \
        f'Number of samples ({len(legal_samples)}) does not match number of sample names ({len(new_sample_names)})!'

    r1_out = []
    r2_out = []
    for i, sam_file in enumerate(sorted(legal_samples)):
        this_call = base_call.copy()
        if new_sample_names == 'auto':
            this_name = parsing.remove_suffixes(sam_file).stem + '_sam2fastq'
        else:
            this_name = new_sample_names[i]
        this_call.append(f"INPUT={sam_file.as_posix()}")
        this_r1 = output_folder.joinpath(f'{this_name}_R1.fastq').as_posix()
        r1_out.append(this_r1)
        this_call.append(f"FASTQ={this_r1}")

        this_r2 = output_folder.joinpath(f'{this_name}_R2.fastq').as_posix()
        r2_out.append(this_r2)
        this_call.append(f"SECOND_END_FASTQ={this_r2}")

        calls.append(this_call)

    _run_picard_calls(calls, script_name, output_folder)

    if return_new_filenames:
        return r1_out, r2_out


def _parse_fastq2sam_misc_args(picard_installation_folder: Union[str, Path, Literal['auto']] = 'auto',
                               quality_score_type: Union[Literal['auto'], Literal[LEGAL_QUAL_SCORE_TYPES]] = 'auto'):
    script_name = 'FastqToSam'
    base_call = _generate_picard_basecall(script_name, picard_installation_folder)
    if quality_score_type == 'auto':
        pass
    elif quality_score_type == 'solexa-quals':
        base_call.append("QUALITY_FORMAT=Solexa")
    elif quality_score_type == 'phred64':
        base_call.append("QUALITY_FORMAT=Illumina")
    elif quality_score_type == 'phred33':
        base_call.append("QUALITY_FORMAT=Standard")
    else:
        warnings.warn(f"Quality format {quality_score_type} is not supported by Picard. "
                      f"Using default quality format.")
    return base_call, script_name


@_func_type('single')
@readable_name('Convert FASTQ files to SAM/BAM files (single-end reads)')
def fastq_to_sam_single(input_folder: Union[str, Path], output_folder: Union[str, Path],
                        picard_installation_folder: Union[str, Path, Literal['auto']] = 'auto',
                        new_sample_names: Union[List[str], Literal['auto']] = 'auto',
                        output_format: Literal['sam', 'bam'] = 'bam',
                        quality_score_type: Union[Literal['auto'], Literal[LEGAL_QUAL_SCORE_TYPES]] = 'auto'):
    """
    Convert SAM/BAM files to FASTQ files using \
    `Picard SamToFastq <https://broadinstitute.github.io/picard/command-line-overview.html#SamToFastq>`_.

    :param input_folder: Path to the folder containing the SAM/BAM files you want to convert.
    :type input_folder: str or Path
    :param output_folder: Path to a folder in which the converted FASTQ files will be saved.
    :type output_folder: str or Path
    :param picard_installation_folder: Path to the installation folder of Picard. For example: \
    'C:/Program Files/Picard'
    :type picard_installation_folder: str, Path, or 'auto' (default='auto')
    :param new_sample_names: Give a new name to each converted sample (optional). \
    If sample_names='auto', sample names \
    will be given automatically. Otherwise, sample_names should be a list of new names, with the order of the names \
    matching the order of the files in the directory.
    :type new_sample_names: list of str or 'auto' (default='auto')
    :return: a list of the paths to the generated FASTQ files.
    :rtype: list of str
    """
    input_folder = Path(input_folder)
    output_folder = Path(output_folder)
    assert input_folder.exists() and input_folder.is_dir(), "input_folder does not exist!"
    assert output_folder.exists() and output_folder.is_dir(), "output_folder does not exist!"

    base_call, script_name = _parse_fastq2sam_misc_args(picard_installation_folder, quality_score_type)
    calls = []

    legal_samples = _get_legal_samples(input_folder, 'alignment')
    assert (new_sample_names == 'auto') or (len(new_sample_names) == len(legal_samples)), \
        f'Number of samples ({len(legal_samples)}) does not match number of sample names ({len(new_sample_names)})!'

    for i, fastq_file in enumerate(sorted(legal_samples)):
        this_call = base_call.copy()
        if new_sample_names == 'auto':
            this_name = parsing.remove_suffixes(fastq_file).stem + '_sam2fastq'
        else:
            this_name = new_sample_names[i]
        this_call.append(f"FASTQ={fastq_file.as_posix()}")
        this_call.append(f"OUTPUT={output_folder.joinpath(f'{this_name}.{output_format}').as_posix()}")
        calls.append(this_call)

    _run_picard_calls(calls, script_name, output_folder)


@_func_type('paired')
@readable_name('Convert FASTQ files to SAM/BAM files (paired-end reads)')
def fastq_to_sam_paired(r1_files: List[str], r2_files: List[str], output_folder: Union[str, Path],
                        picard_installation_folder: Union[str, Path, Literal['auto']] = 'auto',
                        new_sample_names: Union[List[str], Literal['auto']] = 'auto',
                        output_format: Literal['sam', 'bam'] = 'bam',
                        quality_score_type: Union[Literal['auto'], Literal[LEGAL_QUAL_SCORE_TYPES]] = 'auto'):
    """
    Convert SAM/BAM files to FASTQ files using \
    `Picard SamToFastq <https://broadinstitute.github.io/picard/command-line-overview.html#SamToFastq>`_.

    :param input_folder: Path to the folder containing the SAM/BAM files you want to convert.
    :type input_folder: str or Path
    :param output_folder: Path to a folder in which the converted FASTQ files will be saved.
    :type output_folder: str or Path
    :param picard_installation_folder: Path to the installation folder of Picard. For example: \
    'C:/Program Files/Picard'
    :type picard_installation_folder: str, Path, or 'auto' (default='auto')
    :param new_sample_names: Give a new name to each converted sample (optional). \
    If sample_names='auto', sample names \
    will be given automatically. Otherwise, sample_names should be a list of new names, with the order of the names \
    matching the order of the files in the directory.
    :type new_sample_names: list of str or 'auto' (default='auto')
    :return: a list of the paths to the generated FASTQ files.
    :rtype: list of str
    """
    output_folder = Path(output_folder)
    assert output_folder.exists() and output_folder.is_dir(), "output_folder does not exist!"

    base_call, script_name = _parse_fastq2sam_misc_args(picard_installation_folder, quality_score_type)

    calls = []

    for i, (file1, file2) in enumerate(zip(r1_files, r2_files)):
        this_call = base_call.copy()
        file1 = Path(file1)
        file2 = Path(file2)
        if new_sample_names == 'auto':
            this_name = f"{parsing.remove_suffixes(file1).stem}_{parsing.remove_suffixes(file2).stem}_sam2fastq"
        else:
            this_name = new_sample_names[i]
        this_call.append(f"FASTQ={file1.as_posix()}")
        this_call.append(f"FASTQ2={file2.as_posix()}")
        this_call.append(f"OUTPUT={output_folder.joinpath(f'{this_name}.{output_format}').as_posix()}")
        calls.append(this_call)

    _run_picard_calls(calls, script_name, output_folder)


@_func_type('both')
@readable_name('Validate SAM/BAM files')
def validate_sam(input_folder: Union[str, Path], output_folder: Union[str, Path],
                 picard_installation_folder: Union[str, Path, Literal['auto']] = 'auto',
                 verbose: bool = True, is_bisulfite_sequenced: bool = False):
    """
    Validate SAM/BAM files using \
    `Picard ValidateSamFile <https://broadinstitute.github.io/picard/command-line-overview.html#ValidateSamFile>`_.

    :param input_folder: Path to the folder containing the SAM/BAM files you want to validate.
    :type input_folder: str or Path
    :param output_folder: Path to a folder in which the validation reports will be saved.
    :type output_folder: str or Path
    :param picard_installation_folder: Path to the installation folder of Picard. For example: \
    'C:/Program Files/Picard'
    :type picard_installation_folder: str, Path, or 'auto' (default='auto')
    :param new_sample_names: Give a new name to each converted sample (optional). \
    If sample_names='auto', sample names \
    will be given automatically. Otherwise, sample_names should be a list of new names, with the order of the names \
    matching the order of the files in the directory.
    :type new_sample_names: list of str or 'auto' (default='auto')
    :param verbose: If True, the validation report will be verbose. If False, the validation report will be a summary.
    :type verbose: bool (default=True)
    :param is_bisulfite_sequenced: Indicates whether the SAM/BAM file consists of bisulfite sequenced reads. \
    If so, C->T is not counted as en error in computer the value of the NM tag.
    :type is_bisulfite_sequenced: bool (default=False)
    :return: a list of the paths to the generated FASTQ files.
    :rtype: list of str
    """
    script_name = 'ValidateSamFile'
    input_folder = Path(input_folder)
    output_folder = Path(output_folder)
    assert input_folder.exists() and input_folder.is_dir(), "input_folder does not exist!"
    assert output_folder.exists() and output_folder.is_dir(), "output_folder does not exist!"

    base_call = _generate_picard_basecall(script_name, picard_installation_folder)
    base_call.append(f"MODE={'VERBOSE' if verbose else 'SUMMARY'}")
    base_call.append(f"IS_BISULFITE_SEQUENCED={'true' if is_bisulfite_sequenced else 'false'}")

    legal_samples = _get_legal_samples(input_folder, 'alignment')

    calls = []
    for i, sam_file in enumerate(sorted(legal_samples)):
        this_call = base_call.copy()

        this_call.append(f"INPUT={sam_file.as_posix()}")
        this_call.append(f"FASTQ={output_folder.joinpath(f'{sam_file.stem}_report.txt').as_posix()}")
        calls.append(this_call)

    _run_picard_calls(calls, script_name, output_folder)


@_func_type('both')
@readable_name('Sort SAM/BAM files')
def sort_sam(input_folder: Union[str, Path], output_folder: Union[str, Path],
             picard_installation_folder: Union[str, Path, Literal['auto']] = 'auto',
             new_sample_names: Union[List[str], Literal['auto']] = 'auto',
             sort_order: Literal['coordinate', 'queryname', 'duplicate'] = 'coordinate'):
    """
    Sort SAM/BAM files using \
    `Picard SortSam <https://broadinstitute.github.io/picard/command-line-overview.html#SortSam>`_.

    :param input_folder: Path to the folder containing the SAM/BAM files you want to sort.
    :type input_folder: str or Path
    :param output_folder: Path to a folder in which the sorted SAM/BAM files will be saved.
    :type output_folder: str or Path
    :param picard_installation_folder: Path to the installation folder of Picard. For example: \
    'C:/Program Files/Picard'
    :type picard_installation_folder: str, Path, or 'auto' (default='auto')
    :param new_sample_names: Give a new name to each converted sample (optional). \
    If sample_names='auto', sample names \
    will be given automatically. Otherwise, sample_names should be a list of new names, with the order of the names \
    matching the order of the files in the directory.
    :type new_sample_names: list of str or 'auto' (default='auto')
    :param sort_order: The order in which the alignments should be sorted.
    :type sort_order: 'coordinate', 'queryname', or 'duplicate' (default='coordinate')
    """
    script_name = 'SortSam'
    input_folder = Path(input_folder)
    output_folder = Path(output_folder)
    assert input_folder.exists() and input_folder.is_dir(), "input_folder does not exist!"
    assert output_folder.exists() and output_folder.is_dir(), "output_folder does not exist!"

    base_call = _generate_picard_basecall(script_name, picard_installation_folder)

    if sort_order == 'coordinate':
        base_call.append("SORT_ORDER=coordinate")
    elif sort_order == 'queryname':
        base_call.append("SORT_ORDER=queryname")
    elif sort_order == 'duplicate':
        base_call.append("SORT_ORDER=duplicate")
    else:
        raise ValueError(f"Sort order '{sort_order}' is not supported by Picard.")

    legal_samples = _get_legal_samples(input_folder, 'alignment')
    assert (new_sample_names == 'auto') or (len(new_sample_names) == len(legal_samples)), \
        f'Number of samples ({len(legal_samples)}) does not match number of sample names ({len(new_sample_names)})!'

    calls = []
    for i, sam_file in enumerate(sorted(legal_samples)):
        this_call = base_call.copy()
        if new_sample_names == 'auto':
            this_name = parsing.remove_suffixes(sam_file).stem + '_sorted'
        else:
            this_name = new_sample_names[i]
        this_call.append(f"INPUT={sam_file.as_posix()}")
        this_call.append(f"OUTPUT={output_folder.joinpath(f'{this_name}.{sam_file.suffix}').as_posix()}")
        calls.append(this_call)

    _run_picard_calls(calls, script_name, output_folder)


@readable_name('Find PCR/optical duplicates in SAM/BAM files')
def find_duplicates(input_folder: Union[str, Path], output_folder: Union[str, Path],
                    picard_installation_folder: Union[str, Path, Literal['auto']] = 'auto',
                    new_sample_names: Union[List[str], Literal['auto']] = 'auto',
                    output_format: Literal['sam', 'bam'] = 'bam',
                    duplicate_handling: Literal['mark', 'remove_optical', 'remove_all'] = 'remove_all',
                    duplicate_scoring_strategy: Literal[
                        'reference_length', 'sum_of_base_qualities', 'random'] = 'sum_of_base_qualities',
                    optical_duplicate_pixel_distance: int = 100):
    """
    Find duplicate reads in SAM/BAM files using \
    `Picard MarkDuplicates <https://broadinstitute.github.io/picard/command-line-overview.html#MarkDuplicates>`_.

    :param input_folder: Path to the folder containing the SAM/BAM files you want to sort.
    :type input_folder: str or Path
    :param output_folder: Path to a folder in which the sorted SAM/BAM files will be saved.
    :type output_folder: str or Path
    :param picard_installation_folder: Path to the installation folder of Picard. For example: \
    'C:/Program Files/Picard'
    :type picard_installation_folder: str, Path, or 'auto' (default='auto')
    :param new_sample_names: Give a new name to each converted sample (optional). \
    If sample_names='auto', sample names \
    will be given automatically. Otherwise, sample_names should be a list of new names, with the order of the names \
    matching the order of the files in the directory.
    :type new_sample_names: list of str or 'auto' (default='auto')
    :param output_format: Format of the output file.
    :type output_format: 'sam' or 'bam' (default='bam')
    :param duplicate_handling: How to handle detected duplicate reads. If 'mark', duplicate reads will be marked \
    with a 1024 flag. If 'remove_optical', 'optical' duplicates and other duplicates that appear to have arisen \
    from the sequencing process instead of the library preparation process will be removed. \
    If 'remove_all', all duplicate reads will be removed.
    :type duplicate_handling: 'mark', 'remove_optical', or 'remove_all' (default='remove_all')
    :param duplicate_scoring_strategy: How to score duplicate reads. If 'reference_length', the length of the \
    reference sequence will be used. If 'sum_of_base_qualities', the sum of the base qualities will be used.
    :type duplicate_scoring_strategy: 'reference_length', 'sum_of_base_qualities', or 'random' \
    (default='sum_of_base_qualities')
    :param optical_duplicate_pixel_distance: The maximum offset between two duplicate clusters in order to consider \
    them optical duplicates. The default (100) is appropriate for unpatterned versions of the Illumina platform. \
    For the patterned flowcell models, 2500 is moreappropriate. \
    For other platforms and models, users should experiment to find what works best.
    :type optical_duplicate_pixel_distance: int (default=100)
    """
    script_name = 'MarkDuplicates'

    input_folder = Path(input_folder)
    output_folder = Path(output_folder)
    assert input_folder.exists() and input_folder.is_dir(), "input_folder does not exist!"
    assert output_folder.exists() and output_folder.is_dir(), "output_folder does not exist!"

    base_call = _generate_picard_basecall(script_name, picard_installation_folder)

    base_call.append(f"OPTICAL_DUPLICATE_PIXEL_DISTANCE={optical_duplicate_pixel_distance}")

    if duplicate_handling == 'mark':
        base_call.append("REMOVE_DUPLICATES=false")
        base_call.append("REMOVE_SEQUENCING_DUPLICATES=false")
    elif duplicate_handling == 'remove_optical':
        base_call.append("REMOVE_SEQUENCING_DUPLICATES=true")
        base_call.append("REMOVE_DUPLICATES=false")
    elif duplicate_handling == 'remove_all':
        base_call.append("REMOVE_DUPLICATES=true")
    else:
        raise ValueError(f"Duplicate handling method '{duplicate_handling}' is not supported by Picard.")

    if duplicate_scoring_strategy == 'sum_of_base_qualities':
        base_call.append("DUPLICATE_SCORING_STRATEGY=SUM_OF_BASE_QUALITIES")
    elif duplicate_scoring_strategy == 'reference_length':
        base_call.append("DUPLICATE_SCORING_STRATEGY=TOTAL_MAPPED_REFERENCE_LENGTH")
    elif duplicate_scoring_strategy == 'random':
        base_call.append("DUPLICATE_SCORING_STRATEGY=RANDOM")
    else:
        raise ValueError(f"Duplicate scoring strategy '{duplicate_scoring_strategy}' is not supported by Picard.")

    calls = []
    legal_samples = _get_legal_samples(input_folder, 'alignment')
    assert (new_sample_names == 'auto') or (len(new_sample_names) == len(legal_samples)), \
        f'Number of samples ({len(legal_samples)}) does not match number of sample names ({len(new_sample_names)})!'

    for i, sam_file in enumerate(sorted(legal_samples)):
        this_call = base_call.copy()
        if new_sample_names == 'auto':
            this_name = parsing.remove_suffixes(sam_file).stem + '_sam2fastq'
        else:
            this_name = new_sample_names[i]
        this_call.append(f"INPUT={sam_file.as_posix()}")
        this_call.append(f"OUTPUT={output_folder.joinpath(f'{this_name}.{output_format}').as_posix()}")
        this_call.append(f"METRICS_FILE ={output_folder.joinpath(f'{this_name}_metrics.txt').as_posix()}")
        calls.append(this_call)

    _run_picard_calls(calls, script_name, output_folder)


@_func_type('both')
@readable_name('Convert SAM/BAM files to BAM/SAM files')
def convert_sam_format(input_folder: Union[str, Path], output_folder: Union[str, Path],
                       picard_installation_folder: Union[str, Path, Literal['auto']] = 'auto',
                       new_sample_names: Union[List[str], Literal['auto']] = 'auto',
                       output_format: Literal['sam', 'bam'] = 'bam'):
    """
    Convert SAM files to BAM files or vice versa using \
    `Picard SamFormatConverter <https://broadinstitute.github.io/picard/command-line-overview.html#SamFormatConverter>`_.

    :param input_folder: Path to the folder containing the SAM/BAM files you want to convert.
    :type input_folder: str or Path
    :param output_folder: Path to a folder in which the converted FASTQ files will be saved.
    :type output_folder: str or Path
    :param picard_installation_folder: Path to the installation folder of Picard. For example: \
    'C:/Program Files/Picard'
    :type picard_installation_folder: str, Path, or 'auto' (default='auto')
    :param new_sample_names: Give a new name to each converted sample (optional). \
    If sample_names='auto', sample names \
    will be given automatically. Otherwise, sample_names should be a list of new names, with the order of the names \
    matching the order of the files in the directory.
    :type new_sample_names: list of str or 'auto' (default='auto')
    :param output_format: format to convert the files into.
    :type output_format: 'sam' or 'bam' (default='bam')
    """
    script_name = 'SamFormatConverter'

    input_folder = Path(input_folder)
    output_folder = Path(output_folder)
    assert input_folder.exists() and input_folder.is_dir(), "input_folder does not exist!"
    assert output_folder.exists() and output_folder.is_dir(), "output_folder does not exist!"

    base_call = _generate_picard_basecall(script_name, picard_installation_folder)

    legal_samples = _get_legal_samples(input_folder, 'alignment')
    assert (new_sample_names == 'auto') or (len(new_sample_names) == len(legal_samples)), \
        f'Number of samples ({len(legal_samples)}) does not match number of sample names ({len(new_sample_names)})!'

    calls = []
    for i, sam_file in enumerate(sorted(legal_samples)):
        this_call = base_call.copy()
        if new_sample_names == 'auto':
            this_name = parsing.remove_suffixes(sam_file).stem
        else:
            this_name = new_sample_names[i]
        this_call.append(f"INPUT={sam_file.as_posix()}")
        this_call.append(f"OUTPUT={output_folder.joinpath(f'{this_name}.{output_format}').as_posix()}")
        calls.append(this_call)

    _run_picard_calls(calls, script_name, output_folder)


@_func_type('single')
@readable_name('featureCounts count (single-end reads)')
def featurecounts_single_end(input_folder: Union[str, Path], output_folder: Union[str, Path],
                             gtf_file: Union[str, Path],
                             gtf_feature_type: str = 'exon', gtf_attr_name: str = 'gene_id',
                             r_installation_folder: Union[str, Path, Literal['auto']] = 'auto',
                             new_sample_names: Union[List[str], Literal['auto']] = 'auto',
                             stranded: Literal['no', 'forward', 'reverse'] = 'no', min_mapping_quality: int = 0,
                             count_multi_mapping_reads: bool = False, count_multi_overlapping_reads: bool = False,
                             ignore_secondary: bool = True,
                             count_fractionally: bool = False, is_long_read: bool = False,
                             report_read_assignment: Union[Literal['bam', 'sam', 'core'], None] = None,
                             threads: PositiveInt = 1) -> Tuple[filtering.CountFilter, pd.DataFrame, pd.DataFrame]:
    """
    Assign mapped single-end sequencing reads to specified genomic features using \
    `RSubread featureCounts <https://doi.org/10.1093/bioinformatics/btt656>`_.
    Returns a count matrix (CountFilter) containing feature counts for all input files, \
    a DataFrame summarizing the features reads were aligned to, and a DataFrame summarizing the alignment statistics.

    :param input_folder: Path to the folder containing the SAM/BAM files you want to quantfy.
    :type input_folder: str or Path
    :param output_folder: Path to a folder in which the quantified results, as well as the log files, will be saved.
    :type output_folder: str or Path
    :param gtf_file: Path to a GTF annotation file. This file will be used to map reads to features. \
    The chromosome names in the GTF files should match the ones in the index file with which you aligned the reads.
    :type gtf_file: str or Path
    :param gtf_feature_type: the feature type or types used to select rows in the GTF annotation \
    which will be used for read summarization.
    :type gtf_feature_type: str (default='exon')
    :param gtf_attr_name: the attribute type in the GTF annotation which will be used to group features (eg. exons) \
    into meta-features (eg. genes).
    :type gtf_attr_name: str (default='gene_id')
    :param r_installation_folder: Path to the installation folder of R. For example: \
    'C:/Program Files/R/R-4.2.1'
    :type r_installation_folder: str, Path, or 'auto' (default='auto')
    :param new_sample_names: Give a new name to each quantified sample (optional). \
    If sample_names='auto', sample names \
    will be given automatically. Otherwise, sample_names should be a list of new names, with the order of the names \
    matching the order of the alphabetical order of the files in the directory.
    :type new_sample_names: list of str or 'auto' (default='auto')
    :param stranded: Indicates the strandedness of the data. 'no' indicates the data is not stranded. \
    'forward' indicates the data is stranded, where the reads align \
    to the forward strand of a transcript. 'reverse' indicates the data is stranded, where the reads align \
    to the reverse strand of a transcript.
    :type stranded: 'no', 'forward', 'reverse' (default='no')
    :param min_mapping_quality: the minimum mapping quality score a read must satisfy in order to be counted.
    :type min_mapping_quality: int >= 0 (default=0)
    :param count_multi_mapping_reads:  indicating if multi-mapping reads/fragments should be counted \
    ('NH' tag in BAM/SAM files).
    :type count_multi_mapping_reads: bool (default=True)
    :param count_multi_overlapping_reads:  indicating if a read is allowed to be assigned to more than one feature \
    (or meta-feature) if it is found to overlap with more than one feature (or meta-feature).
    :type count_multi_overlapping_reads: bool (default=False)
    :param ignore_secondary: indicating if only primary alignments should be counted. Primary and secondary alignments \
    are identified using bit 0x100 in the Flag field of SAM/BAM files. \
    If True, all primary alignments in a dataset will be counted no matter they are from multi-mapping reads or not.
    :type ignore_secondary: bool (default=True)
    :param count_fractionally: indicating if fractional counts are produced for \
    multi-mapping reads and/or multi-overlapping reads.
    :type count_fractionally: bool (default=False)
    :param is_long_read: indicating if input data contain long reads. \
    This option should be set to True if counting Nanopore or PacBio long reads.
    :type is_long_read: bool (default=False)
    :param report_read_assignment: if not None, featureCounts will generated detailed read assignment results \
    for each read. These results can be saved in one of three formats: BAM, SAM, or CORE.
    :type report_read_assignment: 'bam', 'sam', 'core', or None (default=None)
    :param threads: number of threads to run bowtie2-build on. More threads will generally make index building faster.
    :type threads: int > 0 (default=1)
    :return: a count matrix (CountFilter) containing feature counts for all input files, \
    a DataFrame summarizing the features reads were aligned to, and a DataFrame summarizing the alignment statistics.
    :rtype: (filtering.CountFilter, pd.DataFrame, pd.DataFrame)
    """
    output_folder = Path(output_folder)
    assert output_folder.exists(), 'Output folder does not exist!'
    kwargs = _parse_featurecounts_misc_args(input_folder, output_folder, gtf_file, gtf_feature_type, gtf_attr_name,
                                            stranded, min_mapping_quality, count_multi_mapping_reads,
                                            count_multi_overlapping_reads, ignore_secondary, count_fractionally,
                                            is_long_read, report_read_assignment, threads)

    new_sample_names = _featurecounts_get_sample_names(kwargs['files'], new_sample_names)

    feature_counting.run_featurecounts_analysis(kwargs, output_folder, r_installation_folder)
    counts, annotation, stats = _process_featurecounts_output(output_folder, new_sample_names)
    return counts, annotation, stats


@_func_type('paired')
@readable_name('featureCounts count (paired-end reads)')
def featurecounts_paired_end(input_folder: Union[str, Path], output_folder: Union[str, Path],
                             gtf_file: Union[str, Path], gtf_feature_type: str = 'exon', gtf_attr_name: str = 'gene_id',
                             r_installation_folder: Union[str, Path, Literal['auto']] = 'auto',
                             new_sample_names: Union[List[str], Literal['auto']] = 'auto',
                             stranded: Literal['no', 'forward', 'reverse'] = 'no', min_mapping_quality: int = 0,
                             count_multi_mapping_reads: bool = False, count_multi_overlapping_reads: bool = False,
                             ignore_secondary: bool = True, count_fractionally: bool = False,
                             is_long_read: bool = False, require_both_mapped: bool = True,
                             count_chimeric_fragments: bool = False, min_fragment_length: NonNegativeInt = 50,
                             max_fragment_length: Union[PositiveInt, None] = 600,
                             report_read_assignment: Union[Literal['bam', 'sam', 'core'], None] = None,
                             threads: PositiveInt = 1) -> Tuple[filtering.CountFilter, pd.DataFrame, pd.DataFrame]:
    """
    Assign mapped paired-end sequencing reads to specified genomic features using \
    `RSubread featureCounts <https://doi.org/10.1093/bioinformatics/btt656>`_. \
    Returns a count matrix (CountFilter) containing feature counts for all input files, \
    a DataFrame summarizing the features reads were aligned to, and a DataFrame summarizing the alignment statistics.

    :param input_folder: Path to the folder containing the SAM/BAM files you want to quantfy.
    :type input_folder: str or Path
    :param output_folder: Path to a folder in which the quantified results, \
    as well as the log files and R script used to generate them, will be saved.
    :type output_folder: str or Path
    :param gtf_file: Path to a GTF annotation file. This file will be used to map reads to features. \
    The chromosome names in the GTF files should match the ones in the index file with which you aligned the reads.
    :type gtf_file: str or Path
    :param gtf_feature_type: the feature type or types used to select rows in the GTF annotation \
    which will be used for read summarization.
    :type gtf_feature_type: str (default='exon')
    :param gtf_attr_name: the attribute type in the GTF annotation which will be used to group features (eg. exons) \
    into meta-features (eg. genes).
    :type gtf_attr_name: str (default='gene_id')
    :param r_installation_folder: Path to the installation folder of R. For example: \
    'C:/Program Files/R/R-4.2.1'
    :type r_installation_folder: str, Path, or 'auto' (default='auto')
    :param new_sample_names: Give a new name to each quantified sample (optional). \
    If sample_names='auto', sample names \
    will be given automatically. Otherwise, sample_names should be a list of new names, with the order of the names \
    matching the order of the file pairs.
    :type new_sample_names: list of str or 'auto' (default='auto')
    :param stranded: Indicates the strandedness of the data. 'no' indicates the data is not stranded. \
    'forward' indicates the data is stranded, where the first read in the pair aligns \
    to the forward strand of a transcript. 'reverse' indicates the data is stranded, where the first read in the pair \
    aligns to the reverse strand of a transcript.
    :type stranded: 'no', 'forward', 'reverse' (default='no')
    :param min_mapping_quality: the minimum mapping quality score a read must satisfy in order to be counted. \
    For paired-end reads, at least one end should satisfy this criteria.
    :type min_mapping_quality: int >= 0 (default=0)
    :param count_multi_mapping_reads:  indicating if multi-mapping reads/fragments should be counted \
    ('NH' tag in BAM/SAM files).
    :type count_multi_mapping_reads: bool (default=True)
    :param count_multi_overlapping_reads:  indicating if a read is allowed to be assigned to more than one feature \
    (or meta-feature) if it is found to overlap with more than one feature (or meta-feature).
    :type count_multi_overlapping_reads: bool (default=False)
    :param ignore_secondary: indicating if only primary alignments should be counted. Primary and secondary alignments \
    are identified using bit 0x100 in the Flag field of SAM/BAM files. \
    If True, all primary alignments in a dataset will be counted no matter they are from multi-mapping reads or not.
    :type ignore_secondary: bool (default=True)
    :param count_fractionally: indicating if fractional counts are produced for \
    multi-mapping reads and/or multi-overlapping reads.
    :type count_fractionally: bool (default=False)
    :param is_long_read: indicating if input data contain long reads. \
    This option should be set to True if counting Nanopore or PacBio long reads.
    :type is_long_read: bool (default=False)
    :param report_read_assignment: if not None, featureCounts will generated detailed read assignment results \
    for each read pair. These results can be saved in one of three formats: BAM, SAM, or CORE.
    :type report_read_assignment: 'bam', 'sam', 'core', or None (default=None)
    :param require_both_mapped: indicating if both ends from the same fragment are required to be successfully aligned \
    before the fragment can be assigned to a feature or meta-feature.
    :type require_both_mapped: bool (default=True)
    :param count_chimeric_fragments: indicating whether a chimeric fragment, \
    which has its two reads mapped to different chromosomes, should be counted or not.
    :type count_chimeric_fragments: bool(default=False)
    :param min_fragment_length: The minimum fragment length for valid paired-end alignments. \
    Read pairs with shorter fragments will not be counted.
    :type min_fragment_length: int >= 0 (default=50)
    :param max_fragment_length: The maximum fragment length for valid paired-end alignments. \
    Read pairs with longer fragments will not be counted.
    :type max_fragment_length: int > 0 or None (default=600)
    :param threads: number of threads to run bowtie2-build on. More threads will generally make index building faster.
    :type threads: int > 0 (default=1)
    :return: a count matrix (CountFilter) containing feature counts for all input files, \
    a DataFrame summarizing the features reads were aligned to, and a DataFrame summarizing the alignment statistics.
    :rtype: (filtering.CountFilter, pd.DataFrame, pd.DataFrame)
    """
    output_folder = Path(output_folder)
    assert output_folder.exists(), 'Output folder does not exist!'

    kwargs = _parse_featurecounts_misc_args(input_folder, output_folder, gtf_file, gtf_feature_type, gtf_attr_name,
                                            stranded, min_mapping_quality, count_multi_mapping_reads,
                                            count_multi_overlapping_reads, ignore_secondary, count_fractionally,
                                            is_long_read, report_read_assignment, threads)
    paired_kwargs = {'isPairedEnd': True, 'requireBothEndsMapped': require_both_mapped,
                     'countChimericFragments': count_chimeric_fragments, 'minFragLength': min_fragment_length,
                     'maxFragLength': max_fragment_length, 'countReadPairs': True}
    kwargs.update(paired_kwargs)

    new_sample_names = _featurecounts_get_sample_names(kwargs['files'], new_sample_names)

    feature_counting.run_featurecounts_analysis(kwargs, output_folder, r_installation_folder)
    counts, annotation, stats = _process_featurecounts_output(output_folder, new_sample_names)
    return counts, annotation, stats


def _parse_featurecounts_misc_args(input_folder: Union[str, Path], output_folder: Path, gtf_file: Union[str, Path],
                                   gtf_feature_type: str, gtf_attr_name: str,
                                   stranded: Literal['no', 'forward', 'reverse'], min_mapping_quality: int,
                                   count_multi_mapping_reads: bool, count_multi_overlapping_reads: bool,
                                   ignore_secondary: bool, count_fractionally: bool, is_long_read: bool,
                                   report_read_assignment: Union[Literal['bam', 'sam', 'core'], None],
                                   threads: PositiveInt):
    strand_dict = {'no': 0, 'forward': 1, 'reverse': 2}
    assert stranded in strand_dict, f"Invalid value for 'stranded': '{stranded}'."
    read_assignment_formats = {'bam': 'BAM', 'sam': 'SAM', 'core': 'CORE', None: None}
    assert report_read_assignment in read_assignment_formats, \
        f"Invalid value for 'report_read_asignment': {report_read_assignment}"

    gtf_file = Path(gtf_file)
    assert gtf_file.exists() and gtf_file.is_file(), "'gtf_file' does not exist!"

    input_folder = Path(input_folder)
    assert input_folder.exists() and input_folder.is_dir(), "input_folder does not exist!"
    files = []
    for item in _get_legal_samples(input_folder, 'alignment'):
        files.append(item.as_posix())
    assert len(files) > 0, f"No legal input files were find in input_folder '{input_folder.as_posix()}'!"

    kwargs = {'files': files, 'annot.ext': gtf_file.as_posix(),
              'isGTFAnnotationFile': gtf_file.suffix.lower() != '.saf',
              'GTF.featureType': gtf_feature_type, 'GTF.attrType': gtf_attr_name,
              'allowMultiOverlap': count_multi_overlapping_reads, 'countMultiMappingReads': count_multi_mapping_reads,
              'fraction': count_fractionally, 'isLongRead': is_long_read, 'minMQS': min_mapping_quality,
              'primaryOnly': ignore_secondary, 'strandSpecific': strand_dict[stranded], 'nthreads': threads}

    if report_read_assignment is not None:
        kwargs.update({'reportReads': read_assignment_formats[report_read_assignment],
                       'reportReadsPath': output_folder.as_posix()})
    return kwargs


def _featurecounts_get_sample_names(files: list, new_sample_names):
    if new_sample_names == 'auto':
        new_sample_names = [Path(fname).stem for fname in files]
    else:
        new_sample_names = parsing.data_to_list(new_sample_names)
        assert len(new_sample_names) == len(files), f"The number of samples {len(files)} does not match the number of " \
                                                    f"new sample names ({len(new_sample_names)})!"
    return new_sample_names


def _process_featurecounts_output(output_folder, new_sample_names):
    counts_path = Path(output_folder).joinpath('featureCounts_counts.csv')
    annotation_path = Path(output_folder).joinpath('featureCounts_annotation.csv')
    stats_path = Path(output_folder).joinpath('featureCounts_stats.csv')

    counts = filtering.CountFilter(counts_path)
    counts.df.columns = new_sample_names
    io.save_table(counts.df, counts_path)  # re-save to reflect changes in column names

    annotation = io.load_table(annotation_path, 0)
    io.save_table(annotation, annotation_path)  # re-save to reflect changes in column names

    stats = io.load_table(stats_path, 0)
    stats.columns = ['Status'] + new_sample_names
    io.save_table(stats, stats_path)  # re-save to reflect changes in column names

    return counts, annotation, stats


@readable_name('Bowtie2 build index')
def bowtie2_create_index(genome_fastas: List[Union[str, Path]], output_folder: Union[str, Path],
                         index_name: Union[str, Literal['auto']] = 'auto',
                         bowtie2_installation_folder: Union[str, Path, Literal['auto']] = 'auto',
                         random_seed: Union[NonNegativeInt, None] = None, threads: PositiveInt = 1):
    """
    builds a bowtie index from FASTA formatted files of target sequences (genome). \
    The index files will be saved in the same folder as your first FASTA file, with the `.bt2` suffix. \
    Be aware that there are pre-built bowtie2 indices for popular model organisms. \
    These can be downloaded from  the \
    `bowtie2 website <https://bowtie-bio.sourceforge.net/bowtie2/manual.shtml>`_ (from menu on the right).

    :param genome_fastas: Path to the FASTA file/files which contain reference sequences to be aligned to.
    :type genome_fastas: list of str or Path
    :param output_folder: Path to the folder in which the bowtie2 index files will be saved.
    :type output_folder: str or Path
    :param index_name: The basename of the index files. \
    bowtie2 will create files named index_name.1.bt2, index_name.2.bt2, index_name.3.bt2, index_name.4.bt2, \
    index_name.rev.1.bt2, and index_name.rev.2.bt2. \
    if index_name='auto', the index name used will be the stem of the first supplied genome FASTA file \
    (for example: if the first genome FASTA file is 'path/to/genome.fa.gz', the index name will be 'genome').
    :type index_name: str or 'auto' (default='auto')
    :param bowtie2_installation_folder: Path to the installation folder of bowtie2. For example:
    'C:/Program Files/bowtie2-2.5.1'. if installation folder is set to 'auto', \
    RNAlysis will attempt to find it automatically.
    :type bowtie2_installation_folder: str, Path, or 'auto' (default='auto')
    :param random_seed: if specified, determines the seed for pseudo-random number generator.
    :type random_seed: int >=0 or None (default=None)
    :param threads: number of threads to run bowtie2-build on. More threads will generally make index building faster.
    :type threads: int > 0 (default=1)
    """
    # determine where to create a small or large bowtie2 index
    small_index_max_size = 4 * 1024 ** 3 - 200
    total_size = 0
    for fasta in parsing.data_to_list(genome_fastas):
        fasta = Path(fasta)
        if fasta.exists():
            if fasta.suffix in ['.gz', '.Z']:
                total_size += io.get_gunzip_size(fasta)
            else:
                statinfo = os.stat(fasta)
                total_size += statinfo.st_size

    command = 'bowtie2-build-l' if total_size > small_index_max_size else 'bowtie2-build-s'
    call = io.generate_base_call(command, bowtie2_installation_folder, shell=True, args=['--wrapper', 'basic-0'])

    if random_seed is not None:
        assert isinstance(random_seed, int) and random_seed >= 0, "'random_seed' must be an integer >=0 !"
        call.extend(['--seed', str(random_seed)])

    assert isinstance(threads, int) and threads > 0, "'threads' must be an integer >0 !"
    call.extend(['--threads', str(threads)])

    genome_fastas = [Path(pth) for pth in parsing.data_to_list(genome_fastas)]
    for fasta in genome_fastas:
        assert fasta.exists(), f"the FASTA file '{fasta.as_posix()}' does not exist!"
    call.append(','.join([fasta.as_posix() for fasta in genome_fastas]))

    output_folder = Path(output_folder)
    assert output_folder.exists(), "output_folder does not exist!"

    if index_name == 'auto':
        index_name = parsing.remove_suffixes(genome_fastas[0]).stem
    else:
        assert isinstance(index_name, str), f"'index_name' must be a string, instead got {type(index_name)}."
    call.append(output_folder.joinpath(index_name).as_posix())

    print(f"Running command: \n{' '.join(call)}")
    with tqdm(total=1, desc='Building bowtie2 index', unit='index') as pbar:
        log_filename = Path(output_folder).joinpath(f'bowtie2-build_{index_name}.log').absolute().as_posix()
        io.run_subprocess(call, shell=True, log_filename=log_filename)
        pbar.update()


def _get_legal_samples(in_dir: Union[str, Path], file_type: Literal['sequence', 'alignment'] = 'sequence') -> List[
    Path]:
    in_dir = Path(in_dir)
    assert in_dir.exists(), "'fastq_folder' does not exist!"

    if file_type == 'sequence':
        legal_suffixes = LEGAL_FASTQ_SUFFIXES
    elif file_type == 'alignment':
        legal_suffixes = LEGAL_ALIGNMENT_SUFFIXES
    else:
        raise ValueError(f"Invalid file type '{file_type}'.")

    legal_samples = []
    for item in sorted(in_dir.iterdir()):
        if item.is_file():
            name = item.name
            if any([name.lower().endswith(suffix.lower()) for suffix in legal_suffixes]):
                legal_samples.append(item)
        else:
            legal_samples.extend(_get_legal_samples(item, file_type))
    return legal_samples


@_func_type('single')
@readable_name('ShortStack align (small RNAs)')
def shortstack_align_smallrna(fastq_folder: Union[str, Path], output_folder: Union[str, Path],
                              genome_fasta: Union[str, Path],
                              shortstack_installation_folder: Union[str, Path, Literal['auto']] = 'auto',
                              new_sample_names: Union[List[str], Literal['auto']] = 'auto',
                              known_rnas: Union[str, Path, None] = None,
                              trim_adapter: Union[str, Literal['autotrim'], None] = None,
                              autotrim_key: str = 'TCGGACCAGGCTTCATTCCCC',
                              multimap_mode: Literal['fractional', 'unique', 'random'] = 'fractional',
                              align_only: bool = False, show_secondary_alignments: bool = False,
                              dicer_min_length: PositiveInt = 21, dicer_max_length: PositiveInt = 24,
                              loci_file: Union[str, Path, None] = None, locus: Union[str, None] = None,
                              search_microrna: Union[None, Literal['de-novo', 'known-rnas']] = 'known-rnas',
                              strand_cutoff: Fraction = 0.8, min_coverage: float = 2, pad: PositiveInt = 75,
                              threads: PositiveInt = 1):
    """
    Align small RNA single-end reads from FASTQ files to a reference sequence using the \
    `ShortStack <https://github.com/MikeAxtell/ShortStack>`_ aligner (version 4). \
    ShortStack is currently not supported on computers running Windows.

    :param fastq_folder: Path to the folder containing the FASTQ files you want to quantify
    :type fastq_folder: str or Path
    :param output_folder: Path to a folder in which the aligned reads, as well as the log files, will be saved.
    :type output_folder: str/Path to an existing folder
    :param genome_fasta: Path to the FASTA file which contain the reference sequences to be aligned to.
    :type genome_fasta: str or Path
    :param shortstack_installation_folder: Path to the installation folder of ShortStack. For example: \
    '/home/myuser/anaconda3/envs/myenv/bin'. if installation folder is set to 'auto', \
    RNAlysis will attempt to find it automatically.
    :type shortstack_installation_folder: str, Path, or 'auto' (default='auto')
    :param new_sample_names: Give a new name to each quantified sample (optional). \
    If sample_names='auto', sample names \
    will be given automatically. Otherwise, sample_names should be a list of new names, with the order of the names \
    matching the order of the file pairs.
    :type new_sample_names: list of str or 'auto' (default='auto')
    :param known_rnas: Path to FASTA-formatted file of known small RNAs. \
    FASTA must be formatted such that a single RNA sequence is on one line only. \
    ATCGUatcgu characters are acceptable. \
    These RNAs are typically the sequences of known microRNAs. \
    For instance, a FASTA file of mature miRNAs pulled from https://www.mirbase.org. \
    Providing these data increases the accuracy of MIRNA locus identification.
    :type known_rnas: str, Path, or None (default=None)
    :param trim_adapter: Determines whether ShortStack will attempt to trim the supplied reads. \
    If `trim_adapter` is not provided (default), no trimming will be run. \
    If `trim_adapter` is set to 'autotrim', ShortStack will automatically infer the 3' adapter sequence \
    of the untrimmed reads, and the uses that to coordinate read trimming. \
    If `trim_adapter` is a DNA sequence, ShortStack will trim the reads using the given DNA sequence as the 3' adapter.
    :type trim_adapter: str, 'autotrim', or None (default=None)
    :param autotrim_key: A DNA sequence to use as a known suffix during the autotrim procedure. \
    This parameter is used only if `trim_adapter` is set to 'autotrim'. \
    ShortStack's autotrim discovers the 3' adapter by scanning for reads that begin with the sequence given by \
    `autotrim_key`. This should be the sequence of a small RNA that is known to be highly abundant \
    in all the libraries. The default sequence is for miR166, \
    a microRNA that is present in nearly all plants at high levels. \
    For non-plant experiments, or if the default is not working well, \
    consider providing an alternative to the default.
    :type autotrim_key: str (default="TCGGACCAGGCTTCATTCCCC" (miR166))
    :param multimap_mode: Sets the mode by which multi-mapped reads are handled. \
    These modes are described in Johnson et al. (2016). The default mode ('fractional') has the best performance. \
    In 'fractional' mode, ShortStack will use a fractional weighting scheme for placement of multi-mapped reads. \
    In 'unique' mode, only uniquely-aligned reads are used as weights for placement of multi-mapped reads. \
    In 'random' mode, multi-mapped read placement is random.
    :type multimap_mode: 'fractional', 'unique', or 'random' (default='fractional')
    :param align_only: if True, ShortStack will terminate after the alignment phase; \
    no additional analysis will occur.
    :type align_only: bool (default=False)
    :param show_secondary_alignments: if True, ShortStack will retain secondary alignments for multimapped reads. \
    This will increase bam file size, possibly by a lot.
    :type show_secondary_alignments: bool (default=False)
    :param dicer_min_length: the minimum size (in nucleotides) of a valid small RNA. \
    Together with `dicer_max_length`, \
    this option sets the bounds to discriminate Dicer-derived small RNA loci from other loci. \
    At least 80% of the reads in a given cluster \
    must be in the range indicated by `dicer_min_length` and `dicer_max_length`.
    :type dicer_min_length: positive int (default=21)
    :param dicer_max_length: the maximun size (in nucleotides) of a valid small RNA. \
    Together with `dicer_min_length`, \
    this option sets the bounds to discriminate Dicer-derived small RNA loci from other loci. \
    At least 80% of the reads in a given cluster \
    must be in the range indicated by `dicer_min_length` and `dicer_max_length`.
    :type dicer_max_length: positive int (default=24)
    :param loci_file: Path to a file of pre-determined loci to analyze. \
    This will prevent de-novo discovery of small RNA loci. \
    The file may be in gff3, bed, or simple tab-delimited format (Chr:Start-Stop[tab]Name). \
    Mutually exclusive with `locus`.
    :type loci_file: str, Path, or None (default=None)
    :param locus: A single locus to analyze, \
    given as a string in the format Chr:Start-Stop (using one-based, inclusive numbering). \
    This will prevent de novo discovery of small RNA loci. Mutually exclusive with `loci_file`.
    :type locus: str or None (default=None)
    :param search_microrna: determines whether and how search for microRNAs will be performed. \
    if `search_microrna` is None, ShortStack will not search for microRNAs. \
    This saves computational time, but MIRNA loci will not be differentiated from other types of small RNA clusters. \
    if `search_microrna` is 'known-rnas', t ShortStack will confine MIRNA analysis to loci where one or more \
    queries from the `known_rnas` file are aligned to the genome. \
    if `search_microrna` is 'de-novo', ShortStack will run  a more comprehensive genome-wide scan for MIRNA loci. \
    Discovered loci that do not overlap already known microRNAs should be treated with caution.
    :type search_microrna: 'de-novo', 'known-rnas', or None (default='known-rnas')
    :param strand_cutoff: Floating point number that sets the cutoff for standedness. \
    By default (strand_cutoff=0.8), loci with >80% reads on the top genomic strand are '+' stranded, \
    loci with <20% reads on the top genomic strand are '-' stranded, and all others are unstranded '.'.
    :type strand_cutoff: float between 0.5 and 1 (default=0.8)
    :param min_coverage: Minimum alignment depth, in units of 'reads per million', \
    required to nucleate a small RNA cluster during de-novo cluster search.
    :type min_coverage: float > 0 (default=2)
    :param pad: initial peaks (continuous regions with depth exceeding the argument `min_coverage`) are merged \
    if they are this distance or less from each other.
    :type pad: positive int (default=75)
    :param threads: number of threads to run ShortStack on. More threads will generally make index building faster.
    :type threads: int > 0 (default=1)
    """
    available_multimap_modes = {'fractional': 'f', 'unique': 'u', 'random': 'r'}

    output_folder = Path(output_folder)
    assert output_folder.exists(), "supplied 'output_folder' does not exist!"

    call = io.generate_base_call('ShortStack', shortstack_installation_folder, shell=True)

    genome_fasta = Path(genome_fasta)
    assert genome_fasta.exists(), f"file 'genome_fasta' at {genome_fasta.as_posix()} does not exist!"
    call.extend(['--genomefile', genome_fasta.as_posix()])

    assert multimap_mode in available_multimap_modes, f"Illegal value for 'multimap_mode': '{multimap_mode}'"
    call.extend(['--mmap', available_multimap_modes[multimap_mode]])

    if trim_adapter is not None:
        if trim_adapter == 'auto':
            assert isinstance(autotrim_key, str), f"'autotrim_key' must be a string, instead got {type(autotrim_key)}!"
            call.extend(['--autotrim', '--autotrim_key', autotrim_key])
        elif isinstance(trim_adapter, str):
            call.extend(['--adapter', trim_adapter])
        else:
            raise TypeError(f"Illegal value for 'trim_adapter': '{trim_adapter}'.")

    if align_only:
        call.append('--align_only')

    if show_secondary_alignments:
        call.append('--show_secondaries')

    if loci_file is not None:
        assert locus is None, "Cannot specify both 'loci_file' and 'locus'!"
        loci_file = Path(loci_file)
        assert loci_file.exists() and loci_file.is_file(), \
            f"File 'loci_file' at {loci_file.as_posix()} does not exist!"
        call.extend(['--locifile', loci_file.as_posix()])

    if locus is not None:
        assert isinstance(locus, str), f"'locus' must be a string, instead got {type(locus)}!"
        call.extend(['--locus', locus])

    if search_microrna == 'known-rnas':
        if known_rnas is not None:
            known_rnas = Path(known_rnas)
            assert known_rnas.exists() and known_rnas.is_file(), \
                f"File 'known_rnas' at {known_rnas.as_posix()} does not exist!"
            call.extend(['--knownRNAs', known_rnas.as_posix()])
    else:
        if known_rnas is not None:
            warnings.warn(f"'search_microrna' was set to '{search_microrna}', "
                          f"and therefore parameter 'known_rnas' is ignored")
        if search_microrna == 'de-novo':
            call.append('--dn_mirna')
        elif search_microrna is None:
            call.append('--nohp')
        else:
            raise TypeError(f"Invalid type for 'search_microrna': {type(search_microrna)}. ")

    assert isinstance(dicer_min_length, int) and dicer_min_length > 0, "'dicer_min_length' must be a positive integer!"
    assert isinstance(dicer_max_length, int) and dicer_max_length > 0, "'dicer_max_length' must be a positive integer!"
    assert dicer_min_length <= dicer_max_length, "'dicer_min_length' must be <= 'dicer_max_length'!"
    call.extend(['--dicermin', str(dicer_min_length)])
    call.extend(['--dicermax', str(dicer_max_length)])

    assert isinstance(strand_cutoff,
                      float) and 0.5 < strand_cutoff < 1, \
        "'strand_cutoff' must be a fraction beetween 0.5 and 1 (non-inclusive)!"
    call.extend(['--strand_cutoff', str(strand_cutoff)])

    assert isinstance(min_coverage, (int, float)) and min_coverage > 0, "'min_coverage' must be a positive number!"
    call.extend(['--mincov', str(min_coverage)])

    assert isinstance(pad, int) and pad > 0, "'pad' must be a positive integer!"
    call.extend(['--pad', str(pad)])

    assert isinstance(threads, int) and threads >= 0, "'threads' must be a non-negative int!"
    call.extend(['--threads', str(threads)])

    legal_samples = _get_legal_samples(fastq_folder)

    assert (new_sample_names == 'auto') or (len(new_sample_names) == len(legal_samples)), \
        f'Number of samples ({len(legal_samples)}) does not match number of sample names ({len(new_sample_names)})!'

    calls = []
    for i, item in enumerate(sorted(legal_samples)):
        this_call = call.copy()
        if new_sample_names == 'auto':
            this_name = parsing.remove_suffixes(item).stem
        else:
            this_name = new_sample_names[i]

        this_call.extend(['--readfile', item.as_posix()])
        this_call.extend(['--outdir', output_folder.joinpath(f'{this_name}').as_posix()])
        calls.append(this_call)

    with tqdm(total=len(calls), desc='Aligning reads', unit='files') as pbar:
        for shortstack_call in calls:
            print(f"Running command: \n{' '.join(shortstack_call)}")
            log_filename = Path(output_folder).joinpath(
                f'ShortStack_{Path(shortstack_call[-1]).stem}.log').absolute().as_posix()
            io.run_subprocess(shortstack_call, shell=True, log_filename=log_filename)
            print(f"Files saved successfully at {shortstack_call[-1]}")
            pbar.update(1)


@_func_type('single')
@readable_name('Bowtie2 align (single-end reads)')
def bowtie2_align_single_end(fastq_folder: Union[str, Path], output_folder: Union[str, Path],
                             index_file: Union[str, Path],
                             bowtie2_installation_folder: Union[str, Path, Literal['auto']] = 'auto',
                             new_sample_names: Union[List[str], Literal['auto']] = 'auto',
                             mode: Literal[LEGAL_BOWTIE2_MODES] = 'end-to-end',
                             settings_preset: Literal[LEGAL_BOWTIE2_PRESETS] = 'very-sensitive',
                             ignore_qualities: bool = False,
                             quality_score_type: Literal[LEGAL_QUAL_SCORE_TYPES] = 'phred33',
                             random_seed: NonNegativeInt = 0, threads: PositiveInt = 1):
    """
    Align single-end reads from FASTQ files to a reference sequence using the \
    `bowtie2 <https://bowtie-bio.sourceforge.net/bowtie2>`_ aligner. \
    The FASTQ files will be individually aligned, and the aligned SAM files will be saved in the output folder. \
    You can read more about how bowtie2 works in the \
    `bowtie2 manual <https://bowtie-bio.sourceforge.net/bowtie2/manual.shtml>`_.

    :param fastq_folder: Path to the folder containing the FASTQ files you want to quantify
    :type fastq_folder: str or Path
    :param output_folder: Path to a folder in which the aligned reads, as well as the log files, will be saved.
    :type output_folder: str/Path to an existing folder
    :param index_file: Path to a pre-built bowtie2 index of the target genome. \
    Can either be downloaded from the `bowtie2 website <https://bowtie-bio.sourceforge.net/bowtie2>`_ \
    (menu on the right), or generated manually from FASTA files using the function 'bowtie2_create_index'. \
    Note that bowtie2 indices are composed of multiple files ending with the '.bt2' suffix. \
    All of those files should be in the same location. It is enough to specify the path to one of those files \
    (for example, 'path/to/index.1.bt2'), or to the main name of the index (for example, 'path/to/index').
    :type index_file: str or Path
    :param bowtie2_installation_folder: Path to the installation folder of bowtie2. For example: \
    'C:/Program Files/bowtie2-2.5.1'. if installation folder is set to 'auto', \
    RNAlysis will attempt to find it automatically.
    :type bowtie2_installation_folder: str, Path, or 'auto' (default='auto')
    :param new_sample_names: Give a new name to each quantified sample (optional). \
    If sample_names='auto', sample names \
    will be given automatically. Otherwise, sample_names should be a list of new names, with the order of the names \
    matching the order of the file pairs.
    :type new_sample_names: list of str or 'auto' (default='auto')
    :param mode: determines the alignment mode of bowtie2. \
    end-to-end mode will look for alignments involving all the read characters. \
    local mode will allow 'clipping' of nucleotides from both sides of the read, if that maximizes the alignment score.
    :type mode: 'end-to-end' or 'local' (default='end-to-end')
    :param settings_preset: determines the alignment sensitivity preset. Higher sensitivity will result in more \
    accurate alignments, but will take longer to calculate. You can read more about the settings presets in the \
    `bowtie2 manual`_.
    :type settings_preset: 'very-sensitive', 'sensitive', 'fast', or 'very-fast' (default='very-sensitive')
    :param ignore_qualities: if True, bowtie2 will ignore the qualities of the reads and treat them all as \
    maximum quality.
    :type ignore_qualities: bool (default=False)
    :param quality_score_type: determines the encoding type of the read quality scores. \
    Most modern sequencing setups use phred+33.
    :type quality_score_type: 'phred33', 'phred64', 'solexa-quals', or 'int-quals' (default='phred33')
    :param random_seed:  determines the seed for pseudo-random number generator.
    :type random_seed: int >=0 (default=0)
    :param threads: number of threads to run bowtie2-build on. More threads will generally make index building faster.
    :type threads: int > 0 (default=1)
    """
    output_folder = Path(output_folder)
    call = _parse_bowtie2_misc_args(output_folder, index_file, bowtie2_installation_folder, mode, settings_preset,
                                    ignore_qualities, quality_score_type, random_seed, threads)

    legal_samples = _get_legal_samples(fastq_folder)

    assert (new_sample_names == 'auto') or (len(new_sample_names) == len(legal_samples)), \
        f'Number of samples ({len(legal_samples)}) does not match number of sample names ({len(new_sample_names)})!'

    calls = []
    for i, item in enumerate(sorted(legal_samples)):
        this_call = call.copy()
        if new_sample_names == 'auto':
            this_name = parsing.remove_suffixes(item).stem
        else:
            this_name = new_sample_names[i]

        this_call.extend(['-U', item.as_posix()])
        this_call.extend(['-S', output_folder.joinpath(f'{this_name}.sam').as_posix()])
        calls.append(this_call)

    with tqdm(total=len(calls), desc='Aligning reads', unit='files') as pbar:
        for bt2_call in calls:
            print(f"Running command: \n{' '.join(bt2_call)}")
            log_filename = Path(output_folder).joinpath(
                f'bowtie2-align_{Path(bt2_call[-1]).stem}.log').absolute().as_posix()
            io.run_subprocess(bt2_call, shell=True, log_filename=log_filename)
            print(f"File saved successfully at {bt2_call[-1]}")
            pbar.update(1)


@_func_type('paired')
@readable_name('Bowtie2 align (paired-end reads)')
def bowtie2_align_paired_end(r1_files: List[str], r2_files: List[str], output_folder: Union[str, Path],
                             index_file: Union[str, Path],
                             bowtie2_installation_folder: Union[str, Path, Literal['auto']] = 'auto',
                             new_sample_names: Union[List[str], Literal['auto']] = 'auto',
                             mode: Literal[LEGAL_BOWTIE2_MODES] = 'end-to-end',
                             settings_preset: Literal[LEGAL_BOWTIE2_PRESETS] = 'very-sensitive',
                             ignore_qualities: bool = False,
                             quality_score_type: Literal[LEGAL_QUAL_SCORE_TYPES] = 'phred33',
                             mate_orientations: Literal['fwd-rev', 'rev-fwd', 'fwd-fwd'] = 'fwd-rev',
                             min_fragment_length: NonNegativeInt = 0,
                             max_fragment_length: PositiveInt = 500,
                             allow_individual_alignment: bool = True,
                             allow_disconcordant_alignment: bool = True,
                             random_seed: NonNegativeInt = 0, threads: PositiveInt = 1):
    """
    Align paired-end reads from FASTQ files to a reference sequence using the \
    `bowtie2 <https://bowtie-bio.sourceforge.net/bowtie2>`_ aligner. \
    The FASTQ file pairs will be individually aligned, and the aligned SAM files will be saved in the output folder. \
    You can read more about how bowtie2 works in the \
    `bowtie2 manual <https://bowtie-bio.sourceforge.net/bowtie2/manual.shtml>`_.

    :param r1_files: a list of paths to your Read#1 files. The files should be sorted in tandem with r2_files, \
    so that they line up to form pairs of R1 and R2 files.
    :type r1_files: list of str/Path to existing FASTQ files
    :param r2_files: a list of paths to your Read#2 files. The files should be sorted in tandem with r1_files, \
    so that they line up to form pairs of R1 and R2 files.
    :type r2_files: list of str/Path to existing FASTQ files
    :param output_folder: Path to a folder in which the aligned reads, as well as the log files, will be saved.
    :type output_folder: str/Path to an existing folder
    :param index_file: Path to a pre-built bowtie2 index of the target genome. \
    Can either be downloaded from the `bowtie2 website <https://bowtie-bio.sourceforge.net/bowtie2>`_ \
    (menu on the right), or generated manually from FASTA files using the function 'bowtie2_create_index'. \
    Note that bowtie2 indices are composed of multiple files ending with the '.bt2' suffix. \
    All of those files should be in the same location. It is enough to specify the path to one of those files \
    (for example, 'path/to/index.1.bt2'), or to the main name of the index (for example, 'path/to/index').
    :type index_file: str or Path
    :param bowtie2_installation_folder: Path to the installation folder of bowtie2. For example: \
    'C:/Program Files/bowtie2-2.5.1'. if installation folder is set to 'auto', \
    RNAlysis will attempt to find it automatically.
    :type bowtie2_installation_folder: str, Path, or 'auto' (default='auto')
    :param new_sample_names: Give a new name to each quantified sample (optional). \
    If sample_names='auto', sample names \
    will be given automatically. Otherwise, sample_names should be a list of new names, with the order of the names \
    matching the order of the file pairs.
    :type new_sample_names: list of str or 'auto' (default='auto')
    :param mode: determines the alignment mode of bowtie2. \
    end-to-end mode will look for alignments involving all the read characters. \
    local mode will allow 'clipping' of nucleotides from both sides of the read, if that maximizes the alignment score.
    :type mode: 'end-to-end' or 'local' (default='end-to-end')
    :param settings_preset: determines the alignment sensitivity preset. Higher sensitivity will result in more \
    accurate alignments, but will take longer to calculate. You can read more about the settings presets in the \
    `bowtie2 manual`_.
    :type settings_preset: 'very-sensitive', 'sensitive', 'fast', or 'very-fast' (default='very-sensitive')
    :param ignore_qualities: if True, bowtie2 will ignore the qualities of the reads and treat them all as \
    maximum quality.
    :type ignore_qualities: bool (default=False)
    :param quality_score_type: determines the encoding type of the read quality scores. \
    Most modern sequencing setups use phred+33.
    :type quality_score_type: 'phred33', 'phred64', 'solexa-quals', or 'int-quals' (default='phred33')
    :param mate_orientations:
    :type mate_orientations: 'fwd-rev', 'rev-fwd', or 'fwd-fwd' (default='fwd-rev')
    :param min_fragment_length: The minimum fragment length for valid paired-end alignments.
    :type min_fragment_length: int >= 0 (default=0)
    :param max_fragment_length: The maximum fragment length for valid paired-end alignments.
    :type max_fragment_length: int > 0 (default=500)
    :param allow_individual_alignment:
    :type allow_individual_alignment: bool (default=
    :param allow_disconcordant_alignment:
    :type allow_disconcordant_alignment: bool (default=
    :param random_seed:  determines the seed for pseudo-random number generator.
    :type random_seed: int >=0 (default=0)
    :param threads: number of threads to run bowtie2-build on. More threads will generally make index building faster.
    :type threads: int > 0 (default=1)
    """
    output_folder = Path(output_folder)
    call = _parse_bowtie2_misc_args(output_folder, index_file, bowtie2_installation_folder, mode, settings_preset,
                                    ignore_qualities, quality_score_type, random_seed, threads)

    assert len(r1_files) == len(r2_files), f"Got an uneven number of R1 and R2 files: " \
                                           f"{len(r1_files)} and {len(r2_files)} respectively"
    assert (new_sample_names == 'auto') or (len(new_sample_names) == len(r1_files)), \
        f'Number of samples ({len(r1_files)}) does not match number of sample names ({len(new_sample_names)})!'

    assert isinstance(min_fragment_length, int) and min_fragment_length >= 0, \
        "'min_fragment_len' must be a non-negative int!"
    assert isinstance(max_fragment_length, int) and max_fragment_length >= 0, \
        "'max_fragment_len' must be a non-negative int!"
    call.extend(['-I', str(min_fragment_length)])
    call.extend(['-X', str(max_fragment_length)])

    if not allow_disconcordant_alignment:
        call.append('--no-discorcordant')
    if not allow_individual_alignment:
        call.append('--no-mixed')

    if mate_orientations == 'fwd-rev':
        call.append('--fr')
    elif mate_orientations == 'rev-fwd':
        call.append('--rf')
    elif mate_orientations == 'fwd-fwd':
        call.append('--ff')
    else:
        raise ValueError(f"Invalid value for 'mate_orientations': '{mate_orientations}'.")

    calls = []
    for i, (file1, file2) in enumerate(zip(r1_files, r2_files)):
        file1 = Path(file1)
        file2 = Path(file2)
        this_call = call.copy()
        if new_sample_names == 'auto':
            this_name = f"{parsing.remove_suffixes(file1).stem}_{parsing.remove_suffixes(file2).stem}.sam"
        else:
            this_name = Path(new_sample_names[i]).with_suffix('.sam').name

        this_call.extend(['-1', file1.as_posix(), '-2', file2.as_posix()])
        this_call.extend(['-S', output_folder.joinpath(this_name).as_posix()])
        calls.append(this_call)

    with tqdm(total=len(calls), desc='Aligning reads', unit='file pairs') as pbar:
        for bt2_call in calls:
            print(f"Running command: \n{' '.join(bt2_call)}")
            log_filename = Path(output_folder).joinpath(
                f'bowtie2-align_{Path(bt2_call[-1]).stem}.log').absolute().as_posix()
            io.run_subprocess(bt2_call, shell=True, log_filename=log_filename)
            print(f"File saved successfully at {bt2_call[-1]}")
            pbar.update(1)


def _parse_bowtie2_misc_args(output_folder, index_file: str, bowtie2_installation_folder: Union[str, Path],
                             mode: Literal[LEGAL_BOWTIE2_MODES], settings_preset: Literal[LEGAL_BOWTIE2_PRESETS],
                             ignore_qualities: bool, quality_score_type: Literal[LEGAL_QUAL_SCORE_TYPES],
                             random_seed: NonNegativeInt, threads: PositiveInt):
    output_folder = Path(output_folder)
    index_file = parsing.remove_suffixes(Path(index_file))
    assert output_folder.exists(), "supplied 'output_folder' does not exist!"

    call = io.generate_base_call('bowtie2', bowtie2_installation_folder, shell=True)

    assert mode in LEGAL_BOWTIE2_MODES, f"Invalid value for 'mode': '{mode}'."
    call.append(f'--{mode}')
    assert settings_preset in LEGAL_BOWTIE2_PRESETS, f"Invalid value for 'settings_preset': '{settings_preset}'."
    call.append(f'--{settings_preset}')
    assert quality_score_type in LEGAL_QUAL_SCORE_TYPES, \
        f"Invalid value for 'quality_score_type': '{quality_score_type}'."
    call.append(f'--{quality_score_type}')

    if ignore_qualities:
        call.append('--ignore-quals')

    assert isinstance(random_seed, int) and random_seed >= 0, "'random_seed' must be a non-negative int!"
    call.extend(['--seed', str(random_seed)])
    assert isinstance(threads, int) and threads >= 0, "'threads' must be a non-negative int!"
    call.extend(['--threads', str(threads)])

    call.extend(['-x', index_file.as_posix()])

    return call


@readable_name('Kallisto build index')
def kallisto_create_index(transcriptome_fasta: Union[str, Path],
                          kallisto_installation_folder: Union[str, Path, Literal['auto']] = 'auto',
                          kmer_length: PositiveInt = 31, make_unique: bool = False):
    """
    builds a kallisto index from a FASTA formatted file of target sequences (transcriptome). \
    The index file will be saved in the same folder as your FASTA file, with the `.idx` suffix. \
    Be aware that there are pre-built kallisto indices for popular model organisms. \
    These can be downloaded from the \
    `kallisto transcriptome indices site <https://github.com/pachterlab/kallisto-transcriptome-indices/releases>`_.

    :param transcriptome_fasta: Path to the FASTA file of your desired transcriptome.
    :type transcriptome_fasta: str or Path
    :param kallisto_installation_folder: Path to the installation folder of kallisto. For example: \
    'C:/Program Files/kallisto'. if installation folder is set to 'auto', \
    RNAlysis will attempt to find it automatically.
    :type kallisto_installation_folder: str, Path, or 'auto' (default='auto')
    :param kmer_length: k-mer length of the index.
    :type kmer_length: an odd int between 1 and 31 (default=31)
    :param make_unique: if True, replace repeated target names with unique names.
    :type make_unique: bool (default=False)
    """
    assert isinstance(kmer_length, int), f"parameter 'kmer_length' must be an integer. Instead, got {type(kmer_length)}"
    assert 0 < kmer_length <= 31 and kmer_length % 2 == 1, "'kmer_length' must be an odd integer between 1 and 31"

    call = io.generate_base_call('kallisto', kallisto_installation_folder, 'version')
    call.append('index')

    transcriptome_fasta = Path(transcriptome_fasta)
    assert transcriptome_fasta.exists(), 'the transcriptome FASTA file does not exist!'

    index_filename = parsing.remove_suffixes(transcriptome_fasta).with_suffix('.idx').as_posix()
    call.extend(['-i', index_filename])

    call.extend(['-k', str(kmer_length)])

    if make_unique:
        call.append('--make-unique')

    call.append(transcriptome_fasta.as_posix())

    print(f"Running command: \n{' '.join(call)}")
    with tqdm(total=1, desc='Building kallisto index', unit='index') as pbar:
        log_filename = transcriptome_fasta.parent.joinpath(
            f'kallisto-index_{transcriptome_fasta.stem}.log').absolute().as_posix()
        io.run_subprocess(call, log_filename=log_filename)
        pbar.update()


@_func_type('single')
@readable_name('Kallisto quantify (single-end reads)')
def kallisto_quantify_single_end(fastq_folder: Union[str, Path], output_folder: Union[str, Path],
                                 index_file: Union[str, Path], gtf_file: Union[str, Path],
                                 average_fragment_length: float, stdev_fragment_length: float,
                                 kallisto_installation_folder: Union[str, Path, Literal['auto']] = 'auto',
                                 new_sample_names: Union[List[str], Literal['auto']] = 'auto',
                                 stranded: Literal['no', 'forward', 'reverse'] = 'no',
                                 bootstrap_samples: Union[PositiveInt, None] = None,
                                 **legacy_args) -> filtering.CountFilter:
    """
    Quantify transcript abundance in single-end mRNA sequencing data using \
    `kallisto <https://pachterlab.github.io/kallisto/>`_. \
    The FASTQ files will be individually quantified and saved in the output folder, each in its own sub-folder. \
    Alongside these files, three .csv files will be saved: a per-transcript count estimate table, \
    a per-transcript TPM estimate table, and a per-gene scaled output table. \
    The per-gene scaled output table is generated using the *scaledTPM* method \
    (scaling the TPM estimates up to the library size) as described by \
    `Soneson et al 2015 <https://doi.org/10.12688/f1000research.7563.2>`_ and used in the \
    `tximport <https://ycl6.gitbook.io/guide-to-rna-seq-analysis/differential-expression-analysis/tximport#scaling>`_ \
    R package. This table format is considered un-normalized for library size, \
    and can therefore be used directly by count-based statistical inference tools such as DESeq2.
    *RNAlysis* will return this table once the analysis is finished.

    :param fastq_folder: Path to the folder containing the FASTQ files you want to quantify
    :type fastq_folder: str or Path
    :param output_folder: Path to a folder in which the quantified results, as well as the log files, will be saved. \
    The individual output of each pair of FASTQ files will reside in a different sub-folder within the output folder, \
    and a summarized results table will be saved in the output folder itself.
    :type output_folder: str/Path to an existing folder
    :param index_file: Path to a pre-built kallisto index of the target transcriptome. \
    Can either be downloaded from the \
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
    'C:/Program Files/kallisto'. if installation folder is set to 'auto', \
    RNAlysis will attempt to find it automatically.
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
    and corrects the abundances accordlingly. \
    Note that this feature is not supported by kallisto versions beyond 0.48.0.
    :type learn_bias: bool (default=False)
    :param seek_fusion_genes: if True, does normal quantification, but additionally looks for reads that do not \
    pseudoalign because they are potentially from fusion genes. \
    All output is written to the file fusion.txt in the output folder. \
    Note that this feature is not supported by kallisto versions beyond 0.48.0.
    :type seek_fusion_genes: bool (default=False)
    :param bootstrap_samples: Number of bootstrap samples to be generated. Bootstrap samples do not affect the \
    estimated count values, but generates an additional .hdf5 output file which contains \
    uncertainty estimates for the expression levels. This is required if you use tools such as Sleuth for downstream \
    differential expression analysis, but not for more traditional tools such as DESeq2 and edgeR.
    :type bootstrap_samples: int >0 or None (default=None)
    """
    # handle legacy arguments
    learn_bias = legacy_args.get('learn_bias', False)
    seek_fusion_genes = legacy_args.get('seek_fusion_genes', False)
    if new_sample_names != 'auto':
        new_sample_names = parsing.data_to_list(new_sample_names)
    output_folder = Path(output_folder)

    call = _parse_kallisto_misc_args(output_folder, index_file, kallisto_installation_folder, stranded, learn_bias,
                                     seek_fusion_genes, bootstrap_samples)

    call.append("--single")
    call.extend(["-s", str(stdev_fragment_length)])
    call.extend(["-l", str(average_fragment_length)])

    legal_samples = _get_legal_samples(fastq_folder)

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

    with tqdm(total=len(calls), desc='Quantifying transcript abundance', unit='files') as pbar:
        for kallisto_call in calls:
            io.run_subprocess(kallisto_call)
            print(f"File saved successfully at {kallisto_call[-7]}")
            pbar.update(1)

    return _process_kallisto_outputs(output_folder, gtf_file)


@_func_type('paired')
@readable_name('Kallisto quantify (paired-end reads)')
def kallisto_quantify_paired_end(r1_files: List[str], r2_files: List[str], output_folder: Union[str, Path],
                                 index_file: Union[str, Path], gtf_file: Union[str, Path],
                                 kallisto_installation_folder: Union[str, Path, Literal['auto']] = 'auto',
                                 new_sample_names: Union[List[str], Literal['auto']] = 'auto',
                                 stranded: Literal['no', 'forward', 'reverse'] = 'no',
                                 bootstrap_samples: Union[PositiveInt, None] = None,
                                 **legacy_args) -> filtering.CountFilter:
    """
    Quantify transcript abundance in paired-end mRNA sequencing data using \
    `kallisto <https://pachterlab.github.io/kallisto/>`_. \
    The FASTQ file pairs will be individually quantified and saved in the output folder, each in its own sub-folder. \
    Alongside these files, three .csv files will be saved: a per-transcript count estimate table, \
    a per-transcript TPM estimate table, and a per-gene scaled output table. \
    The per-gene scaled output table is generated using the *scaledTPM* method \
    (scaling the TPM estimates up to the library size) as described by \
    `Soneson et al 2015 <https://doi.org/10.12688/f1000research.7563.2>`_ and used in the \
    `tximport <https://ycl6.gitbook.io/guide-to-rna-seq-analysis/differential-expression-analysis/tximport#scaling>`_ \
    R package. This table format is considered un-normalized for library size, \
    and can therefore be used directly by count-based statistical inference tools such as DESeq2.
    *RNAlysis* will return this table once the analysis is finished.

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
        'C:/Program Files/kallisto'. if installation folder is set to 'auto', \
    RNAlysis will attempt to find it automatically.
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
    and corrects the abundances accordlingly. \
    Note that this feature is not supported by kallisto versions beyond 0.48.0.
    :type learn_bias: bool (default=False)
    :param seek_fusion_genes: if True, does normal quantification, but additionally looks for reads that do not \
    pseudoalign because they are potentially from fusion genes. \
    All output is written to the file fusion.txt in the output folder. \
    Note that this feature is not supported by kallisto versions beyond 0.48.0.
    :type seek_fusion_genes: bool (default=False)
    :param bootstrap_samples: Number of bootstrap samples to be generated. Bootstrap samples do not affect the \
    estimated count values, but generates an additional .hdf5 output file which contains \
    uncertainty estimates for the expression levels. This is required if you use tools such as Sleuth for downstream \
    differential expression analysis, but not for more traditional tools such as DESeq2 and edgeR.
    :type bootstrap_samples: int >0 or None (default=None)
    """
    if new_sample_names != 'auto':
        new_sample_names = parsing.data_to_list(new_sample_names)
    assert len(r1_files) == len(r2_files), f"Got an uneven number of R1 and R2 files: " \
                                           f"{len(r1_files)} and {len(r2_files)} respectively"
    assert (new_sample_names == 'auto') or (len(new_sample_names) == len(r1_files)), \
        f'Number of samples ({len(r1_files)}) does not match number of sample names ({len(new_sample_names)})!'

    # handle legacy arguments
    learn_bias = legacy_args.get('learn_bias', False)
    seek_fusion_genes = legacy_args.get('seek_fusion_genes', False)

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

    with tqdm(total=len(calls), desc='Quantifying transcript abundance', unit='file pairs') as pbar:
        for kallisto_call in calls:
            io.run_subprocess(kallisto_call)
            print(f"Files saved successfully at {kallisto_call[-3]}")
            pbar.update(1)

    return _process_kallisto_outputs(output_folder, gtf_file)


def _process_kallisto_outputs(output_folder, gtf_file):
    counts, tpm = _merge_kallisto_outputs(output_folder)
    genes_scaled_tpm = _sum_transcripts_to_genes(tpm, counts, gtf_file)

    io.save_table(counts, output_folder.joinpath('transcript_counts.csv'))
    io.save_table(tpm, output_folder.joinpath('transcript_tpm.csv'))
    io.save_table(genes_scaled_tpm, output_folder.joinpath('kallisto_output_scaled_per_gene.csv'))

    return filtering.CountFilter.from_dataframe(genes_scaled_tpm, 'kallisto_output_scaled_per_gene',
                                                is_normalized=False)


def _parse_kallisto_misc_args(output_folder, index_file: str, kallisto_installation_folder: Union[str, Path],
                              stranded: Literal['no', 'forward', 'reverse'] = 'no', learn_bias: bool = False,
                              seek_fusion_genes: bool = False, bootstrap_samples: Union[int, None] = None):
    # handle legacy arguments
    if learn_bias:
        warnings.warn("The 'learn_bias' argument is no longer supported by kallisto versions beyond 0.48.0. ")
    if seek_fusion_genes:
        warnings.warn("The 'seek_fusion_genes' argument is no longer supported by kallisto versions beyond 0.48.0. ")

    output_folder = Path(output_folder)
    index_file = Path(index_file)
    assert output_folder.exists(), "supplied 'output_folder' does not exist!"
    assert index_file.exists(), "supplied 'index_file' does not exist!"
    assert isinstance(stranded, str) and stranded.lower() in ["no", "forward", "reverse"], \
        f"invalid value for parameter 'stranded': {stranded}"

    call = io.generate_base_call('kallisto', kallisto_installation_folder, 'version')
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
    for item in sorted(Path(output_folder).iterdir()):
        if item.is_dir():
            abundance_path = item.joinpath('abundance.tsv')
            if abundance_path.exists():
                this_df = io.load_table(abundance_path, index_col=0)
                sample_name = item.name
                counts[sample_name] = this_df['est_counts']
                tpm[sample_name] = this_df['tpm']
    return counts, tpm


def _sum_transcripts_to_genes(tpm: pd.DataFrame, counts: pd.DataFrame, gtf_path: Union[str, Path]):
    with tqdm(desc='Parsing GTF file', total=8) as pbar:
        for use_name in [False, True]:
            for use_version in [True, False]:
                for split_ids in [True, False]:
                    transcript_to_gene_map = genome_annotation.map_transcripts_to_genes(gtf_path, use_name, use_version,
                                                                                        split_ids)
                    pbar.update(1)
                    if len(transcript_to_gene_map) == 0:
                        continue

                    library_sizes = counts.sum(axis=0) / (10 ** 6)
                    tpm_cpy = tpm.copy()
                    tpm_cpy['Gene ID'] = pd.Series(transcript_to_gene_map)
                    tpm_by_gene = tpm_cpy.groupby('Gene ID').sum()

                    if tpm_by_gene.shape[0] == 0:
                        continue
                    scaled_tpm = tpm_by_gene.multiply(library_sizes, axis=1)
                    pbar.update(8)
                    if isinstance(scaled_tpm, pd.Series):
                        scaled_tpm = scaled_tpm.to_frame()
                    return scaled_tpm

    raise ValueError("Failed to map transcripts to genes with the given GTF file!")


@_func_type('single')
@readable_name('CutAdapt (single-end reads)')
def trim_adapters_single_end(fastq_folder: Union[str, Path], output_folder: Union[str, Path],
                             three_prime_adapters: Union[None, str, List[str]],
                             five_prime_adapters: Union[None, str, List[str]] = None,
                             any_position_adapters: Union[None, str, List[str]] = None,
                             new_sample_names: Union[List[str], Literal['auto']] = 'auto',
                             quality_trimming: Union[NonNegativeInt, None] = 20, trim_n: bool = True,
                             minimum_read_length: NonNegativeInt = 10,
                             maximum_read_length: Union[PositiveInt, None] = None, discard_untrimmed_reads: bool = True,
                             error_tolerance: Fraction = 0.1, minimum_overlap: NonNegativeInt = 3,
                             allow_indels: bool = True, parallel: bool = True, gzip_output: bool = False):
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
    :param gzip_output: if True, gzips the output FASTQ files.
    :type gzip_output: bool (default=False)
    :param new_sample_names: Give a new name to each trimmed sample (optional). \
    If sample_names='auto', sample names \
    will be given automatically. Otherwise, sample_names should be a list of new names, with the order of the names \
    matching the alphabetical order of the input files.
    :type new_sample_names: list of str or 'auto' (default='auto')
    """
    if not HAS_CUTADAPT:
        warnings.warn("Python package 'cutadapt' is not installed. \n"
                      "If you want to use the adapter trimming feature, "
                      "please install python package 'cutadapt' and try again. ")
        return

    try:
        call = io.generate_base_call('cutadapt', 'auto')
        found_cli = True
    except FileNotFoundError:
        call = []
        found_cli = False

    for adapter_group, prefix in zip([three_prime_adapters, five_prime_adapters, any_position_adapters],
                                     ['--adapter', '--front', '--anywhere']):
        if adapter_group is not None:
            for adapter in parsing.data_to_list(adapter_group):
                assert isinstance(adapter, str), f"The following adapter is invalid: {adapter}"
                call.extend([prefix, adapter])

    call.extend(_parse_cutadapt_misc_args(quality_trimming, trim_n, minimum_read_length, maximum_read_length,
                                          discard_untrimmed_reads, error_tolerance, minimum_overlap, allow_indels,
                                          parallel))

    legal_samples = _get_legal_samples(fastq_folder)
    assert (new_sample_names == 'auto') or (len(new_sample_names) == len(legal_samples)), \
        f'Number of samples ({len(legal_samples)}) does not match number of sample names ({len(new_sample_names)})!'

    calls = []
    for i, item in enumerate(legal_samples):
        if item.is_file():
            name = item.name
            if any([name.endswith(suffix) for suffix in LEGAL_FASTQ_SUFFIXES]):
                if new_sample_names == 'auto':
                    suffix = '_trimmed.fastq.gz' if gzip_output else '_trimmed.fastq'
                    this_name = parsing.remove_suffixes(item).stem + suffix
                else:
                    suffix = '.fastq.gz' if gzip_output else '.fastq'
                    this_name = new_sample_names[i] + suffix

                output_path = Path(output_folder).joinpath(this_name)
                this_call = call.copy()
                this_call.extend(['--output', output_path.as_posix(), item.as_posix()])
                calls.append(this_call)

    with tqdm(total=len(calls), desc='Trimming adapters', unit='files') as pbar:
        for cutadapt_call in calls:
            infile_stem = parsing.remove_suffixes(Path(cutadapt_call[-1])).stem
            log_filename = Path(output_folder).joinpath(f'cutadapt_log_{infile_stem}.log').absolute().as_posix()

            if found_cli:
                io.run_subprocess(cutadapt_call, log_filename=log_filename)
            else:
                cutadapt_main(cutadapt_call)
            print(f"File saved successfully at {cutadapt_call[-2]}")
            pbar.update(1)


@_func_type('paired')
@readable_name('CutAdapt (paired-end reads)')
def trim_adapters_paired_end(r1_files: List[Union[str, Path]], r2_files: List[Union[str, Path]],
                             output_folder: Union[str, Path],
                             three_prime_adapters_r1: Union[None, str, List[str]],
                             three_prime_adapters_r2: Union[None, str, List[str]],
                             five_prime_adapters_r1: Union[None, str, List[str]] = None,
                             five_prime_adapters_r2: Union[None, str, List[str]] = None,
                             any_position_adapters_r1: Union[None, str, List[str]] = None,
                             any_position_adapters_r2: Union[None, str, List[str]] = None,
                             new_sample_names: Union[List[str], Literal['auto']] = 'auto',
                             quality_trimming: Union[NonNegativeInt, None] = 20, trim_n: bool = True,
                             minimum_read_length: NonNegativeInt = 10,
                             maximum_read_length: Union[PositiveInt, None] = None,
                             discard_untrimmed_reads: bool = True,
                             pair_filter_if: Literal['both', 'any', 'first'] = 'both',
                             error_tolerance: Fraction = 0.1, minimum_overlap: NonNegativeInt = 3,
                             allow_indels: bool = True, parallel: bool = True, gzip_output: bool = False,
                             return_new_filenames: bool = False):
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
    :param gzip_output: if True, gzips the output FASTQ files.
    :type gzip_output: bool (default=False)
    :param new_sample_names: Give a new name to each trimmed sample (optional). \
    If sample_names='auto', sample names \
    will be given automatically. Otherwise, sample_names should be a list of new names, with the order of the names \
    matching the order of the file pairs.
    :type new_sample_names: list of str or 'auto' (default='auto')
    """
    if not HAS_CUTADAPT:
        warnings.warn("Python package 'cutadapt' is not installed. \n"
                      "If you want to use the adapter trimming feature, "
                      "please install python package 'cutadapt' and try again. ")
        return
    assert len(r1_files) == len(r2_files), f"Got an uneven number of R1 and R2 files: " \
                                           f"{len(r1_files)} and {len(r2_files)} respectively"
    assert (new_sample_names == 'auto') or (len(new_sample_names) == len(r1_files)), \
        f'Number of samples ({len(r1_files)}) does not match number of sample names ({len(new_sample_names)})!'

    try:
        call = io.generate_base_call('cutadapt', 'auto')
        found_cli = True
    except FileNotFoundError:
        call = []
        found_cli = False

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
    r1_out = []
    r2_out = []
    for i, (file1, file2) in enumerate(zip(r1_files, r2_files)):
        file1 = Path(file1)
        file2 = Path(file2)
        if new_sample_names == 'auto':
            suffix = '_trimmed.fastq.gz' if gzip_output else '_trimmed.fastq'
            base_name_r1 = f"{parsing.remove_suffixes(file1).stem}{suffix}"
            base_name_r2 = f"{parsing.remove_suffixes(file2).stem}{suffix}"
        else:
            suffix = '.fastq.gz' if gzip_output else '.fastq'
            base_name_r1 = f"{new_sample_names[i]}_R1{suffix}"
            base_name_r2 = f"{new_sample_names[i]}_R2{suffix}"

        output_path_r1 = Path(output_folder).joinpath(base_name_r1)
        output_path_r2 = Path(output_folder).joinpath(base_name_r2)
        r1_out.append(output_path_r1)
        r2_out.append(output_path_r2)

        this_call = call.copy()
        this_call.extend(['--output', output_path_r1.as_posix()])
        this_call.extend(['--paired-output', output_path_r2.as_posix()])
        this_call.extend([file1.as_posix(), file2.as_posix()])
        calls.append(this_call)

    with tqdm(total=len(calls), desc='Trimming adapters', unit='file pairs') as pbar:
        for cutadapt_call in calls:
            infile1_stem = parsing.remove_suffixes(Path(cutadapt_call[-2])).stem
            infile2_stem = parsing.remove_suffixes(Path(cutadapt_call[-1])).stem
            log_filename = Path(output_folder).joinpath(
                f'cutadapt_log_{infile1_stem}_{infile2_stem}.log').absolute().as_posix()
            if found_cli:
                io.run_subprocess(cutadapt_call, log_filename=log_filename)
            else:
                cutadapt_main(cutadapt_call)
            print(f"Files saved successfully at {cutadapt_call[-2]} and  {cutadapt_call[-1]}")
            pbar.update(1)

    if return_new_filenames:
        return r1_out, r2_out


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
