import shutil
import subprocess
import sys

import numpy as np
import polars.selectors as cs
import pytest

from rnalysis.utils.differential_expression import *


class TestLimmaVoomRunner:
    @pytest.mark.parametrize(
        'data,design_matrix,comparisons,random_effect,expected_path',
        [
            (
                'tests/test_files/big_counted.csv',
                'tests/test_files/test_design_matrix.csv',
                [('condition', 'cond2', 'cond1')],
                None,
                'tests/test_files/limma_tests/case1/expected_limma_script_1.R',
            ),
            (
                'counted.csv',
                'tests/test_files/test_design_matrix.csv',
                [('condition', 'cond3', 'cond2'), ('replicate', 'rep2', 'rep1'), ('condition', 'cond1', 'cond2')],
                None,
                'tests/test_files/limma_tests/case2/expected_limma_script_2.R',
            ),
            (
                'tests/test_files/big_counted.csv',
                'tests/test_files/test_design_matrix.csv',
                [('condition', 'cond2', 'cond1')],
                'replicate',
                'tests/test_files/limma_tests/case3/expected_limma_script_3.R',
            ),
        ],
    )
    def test_create_limma_script(self, data, design_matrix, comparisons, random_effect, expected_path):
        with open(expected_path) as f:
            expected = f.read()
        runner = LimmaVoomRunner(data, design_matrix, comparisons, random_effect=random_effect)
        out_path = runner.create_script()
        assert Path(out_path).exists()
        with open(out_path) as f:
            out = f.read()

        expected = expected.replace('*PLACEHOLDERPATH*', out_path.parent.as_posix())
        assert out == expected

    @pytest.mark.parametrize(
        'comparisons,expected_paths',
        [
            (
                [('replicate', 'rep2', 'rep3')],
                ['tests/test_files/limma_tests/case1/LimmaVoom_replicate_rep2_vs_rep3_truth.csv'],
            ),
            (
                [('condition', 'cond2', 'cond1'), ('condition', 'cond3', 'cond2')],
                [
                    'tests/test_files/limma_tests/case2/LimmaVoom_condition_cond2_vs_cond1_truth.csv',
                    'tests/test_files/limma_tests/case2/LimmaVoom_condition_cond3_vs_cond2_truth.csv',
                ],
            ),
        ],
    )
    def test_run_limma_analysis(self, comparisons, expected_paths):
        data_path = 'tests/test_files/big_counted.csv'
        design_mat_path = 'tests/test_files/test_design_matrix.csv'
        runner = LimmaVoomRunner(data_path, design_mat_path, comparisons)
        output_dir = runner.run()

        dfs = []
        for item in output_dir.iterdir():
            if item.is_file() and item.suffix == '.csv':
                dfs.append(io.load_table(item))

        expected_dfs = []
        for item in expected_paths:
            expected_dfs.append(io.load_table(item))

        for out, truth in zip(dfs, expected_dfs):
            out = out.sort(pl.first())
            truth = truth.sort(pl.first())
            assert out.shape == truth.shape
            assert np.all(out.columns == truth.columns)
            assert out.select(pl.first()).equals(truth.select(pl.first()))
            if sys.platform == 'win32':  # running DESeq in linux gives slightly different results
                assert np.allclose(out.drop(cs.first()), truth.drop(cs.first()), equal_nan=True, atol=1 * 10 ** (-4))

    @pytest.mark.parametrize(
        'data,design_matrix,comparisons,covariates,expected_path',
        [
            (
                'tests/test_files/big_counted.csv',
                'tests/test_files/test_design_matrix_advanced.csv',
                [('condition', 'cond2', 'cond1')],
                ['covariate1', 'covariate2'],
                'tests/test_files/limma_tests/case4/expected_limma_script_covariates.R',
            ),
        ],
    )
    def test_create_limma_script_with_covariates(self, data, design_matrix, comparisons, covariates, expected_path):
        with open(expected_path) as f:
            expected = f.read()
        runner = LimmaVoomRunner(data, design_matrix, comparisons, covariates=covariates)
        out_path = runner.create_script()
        assert Path(out_path).exists()
        with open(out_path) as f:
            out = f.read()

        expected = expected.replace('*PLACEHOLDERPATH*', out_path.parent.as_posix())
        assert out == expected

    @pytest.mark.parametrize(
        'data,design_matrix,comparisons,lrt_factors,expected_path',
        [
            (
                'tests/test_files/big_counted.csv',
                'tests/test_files/test_design_matrix_advanced.csv',
                [('condition', 'cond2', 'cond1')],
                ['factor1', 'factor2'],
                'tests/test_files/limma_tests/case5/expected_limma_script_lrt_factors.R',
            ),
        ],
    )
    def test_create_limma_script_with_lrt_factors(self, data, design_matrix, comparisons, lrt_factors, expected_path):
        with open(expected_path) as f:
            expected = f.read()
        runner = LimmaVoomRunner(data, design_matrix, comparisons, lrt_factors=lrt_factors)
        out_path = runner.create_script()
        assert Path(out_path).exists()
        with open(out_path) as f:
            out = f.read()

        expected = expected.replace('*PLACEHOLDERPATH*', out_path.parent.as_posix())
        assert out == expected

    @pytest.mark.parametrize(
        'data,design_matrix,comparisons,covariates,lrt_factors,expected_path',
        [
            (
                'tests/test_files/big_counted.csv',
                'tests/test_files/test_design_matrix_advanced.csv',
                [('condition', 'cond2', 'cond1')],
                ['covariate1', 'covariate2'],
                ['factor1', 'factor2'],
                'tests/test_files/limma_tests/case6/expected_limma_script_combined.R',
            ),
        ],
    )
    def test_create_limma_script_combined(
        self, data, design_matrix, comparisons, covariates, lrt_factors, expected_path
    ):
        with open(expected_path) as f:
            expected = f.read()
        runner = LimmaVoomRunner(data, design_matrix, comparisons, covariates=covariates, lrt_factors=lrt_factors)
        out_path = runner.create_script()
        assert Path(out_path).exists()
        with open(out_path) as f:
            out = f.read()

        expected = expected.replace('*PLACEHOLDERPATH*', out_path.parent.as_posix())
        assert out == expected

    @pytest.mark.parametrize(
        'model_factors,expected_coef_names',
        [
            (
                ['factor1', 'factor2', 'factor3'],
                {'factor1': {'factor1B', 'factor1C'}, 'factor2': {'factor2'}, 'factor3': {'factor3Y', 'factor3Z'}},
            ),
            (
                ['factor1', 'poly(factor2, degree = 2)', 'factor3'],
                {
                    'factor1': {'factor1B', 'factor1C'},
                    'poly(factor2, degree = 2)': {'poly.factor2..degree...2.1', 'poly.factor2..degree...2.2'},
                    'factor3': {'factor3Y', 'factor3Z'},
                },
            ),
            (
                ['factor1', 'factor2', 'factor1:factor2'],
                {
                    'factor1': {'factor1B', 'factor1C'},
                    'factor2': {'factor2'},
                    'factor1:factor2': {'factor1B.factor2', 'factor1C.factor2'},
                },
            ),
            (
                ['factor1', 'factor2', 'factor3', 'factor1:factor2:factor3'],
                {
                    'factor1': {'factor1B', 'factor1C'},
                    'factor2': {'factor2'},
                    'factor3': {'factor3Y', 'factor3Z'},
                    'factor1:factor2:factor3': {
                        'factor1B.factor2.factor3Y',
                        'factor1B.factor2.factor3Z',
                        'factor1C.factor2.factor3Y',
                        'factor1C.factor2.factor3Z',
                    },
                },
            ),
        ],
    )
    def test_get_coef_names(self, model_factors, expected_coef_names):
        design_matrix = 'tests/test_files/limma_tests/test_design_mat.csv'
        runner = LimmaVoomRunner('data.csv', design_matrix, [('factor1', 'B', 'A')], model_factors=model_factors)
        coef_names = runner._get_coef_names()
        assert coef_names == expected_coef_names

    @pytest.mark.parametrize(
        'interaction,expected_terms',
        [
            ('factor1:factor2', {'factor1B:factor2', 'factor1C:factor2'}),
            ('factor1:factor3', {'factor1B:factor3Y', 'factor1B:factor3Z', 'factor1C:factor3Y', 'factor1C:factor3Z'}),
            ('factor2:factor3', {'factor2:factor3Y', 'factor2:factor3Z'}),
            (
                'factor1:factor2:factor3',
                {
                    'factor1B:factor2:factor3Y',
                    'factor1B:factor2:factor3Z',
                    'factor1C:factor2:factor3Y',
                    'factor1C:factor2:factor3Z',
                },
            ),
        ],
    )
    def test_generate_interaction_terms(self, interaction, expected_terms):
        design_matrix = 'tests/test_files/limma_tests/test_design_mat.csv'
        runner = LimmaVoomRunner('data.csv', design_matrix, [('factor1', 'B', 'A')])
        interaction_terms = runner._generate_interaction_terms(interaction)
        assert interaction_terms == expected_terms


class TestDESeqRunner:
    @pytest.mark.parametrize(
        'data,design_matrix,comparisons,cooks_cutoff,expected_path',
        [
            (
                'tests/test_files/big_counted.csv',
                'tests/test_files/test_design_matrix.csv',
                [('condition', 'cond2', 'cond1')],
                False,
                'tests/test_files/deseq2_tests/case1/expected_deseq_script_1.R',
            ),
            (
                'counted.csv',
                'tests/test_files/test_design_matrix.csv',
                [('condition', 'cond3', 'cond2'), ('replicate', 'rep2', 'rep1'), ('condition', 'cond1', 'cond2')],
                True,
                'tests/test_files/deseq2_tests/case2/expected_deseq_script_2.R',
            ),
        ],
    )
    def test_create_deseq2_script(self, data, design_matrix, comparisons, cooks_cutoff, expected_path):
        runner = DESeqRunner(data, design_matrix, comparisons, cooks_cutoff=cooks_cutoff)
        with open(expected_path) as f:
            expected = f.read()

        out_path = runner.create_script()
        assert Path(out_path).exists()
        with open(out_path) as f:
            out = f.read()

        expected = expected.replace('*PLACEHOLDERPATH*', out_path.parent.as_posix())
        assert out == expected

    @pytest.mark.parametrize(
        'comparisons,expected_paths',
        [
            (
                [('replicate', 'rep2', 'rep3')],
                ['tests/test_files/deseq2_tests/case1/DESeq2_replicate_rep2_vs_rep3_truth.csv'],
            ),
            (
                [('condition', 'cond2', 'cond1'), ('condition', 'cond3', 'cond2')],
                [
                    'tests/test_files/deseq2_tests/case2/DESeq2_condition_cond2_vs_cond1_truth.csv',
                    'tests/test_files/deseq2_tests/case2/DESeq2_condition_cond3_vs_cond2_truth.csv',
                ],
            ),
        ],
    )
    def test_run_deseq2_analysis(self, comparisons, expected_paths):
        data_path = 'tests/test_files/big_counted.csv'
        design_mat_path = 'tests/test_files/test_design_matrix.csv'
        runner = DESeqRunner(data_path, design_mat_path, comparisons)
        output_dir = runner.run()

        dfs = []
        for item in output_dir.iterdir():
            if item.is_file() and item.suffix == '.csv':
                dfs.append(io.load_table(item))

        expected_dfs = []
        for item in expected_paths:
            expected_dfs.append(io.load_table(item))

        for out, truth in zip(dfs, expected_dfs):
            out = out.sort(pl.first())
            truth = truth.sort(pl.first())
            assert out.shape == truth.shape
            assert np.all(out.columns == truth.columns)
            assert out.select(pl.first()).equals(truth.select(pl.first()))
            if sys.platform == 'win32':  # running DESeq in linux gives slightly different results
                assert np.allclose(out.drop(cs.first()), truth.drop(cs.first()), equal_nan=True, atol=1 * 10 ** (-4))

    @pytest.mark.parametrize(
        'data,design_matrix,comparisons,covariates,expected_path',
        [
            (
                'tests/test_files/big_counted.csv',
                'tests/test_files/test_design_matrix.csv',
                [('condition', 'cond2', 'cond1')],
                ['covariate1', 'covariate2'],
                'tests/test_files/deseq2_tests/case3/expected_deseq_script_covariates.R',
            ),
        ],
    )
    def test_create_deseq2_script_with_covariates(self, data, design_matrix, comparisons, covariates, expected_path):
        runner = DESeqRunner(data, design_matrix, comparisons, covariates=covariates)
        with open(expected_path) as f:
            expected = f.read()

        out_path = runner.create_script()
        assert Path(out_path).exists()
        with open(out_path) as f:
            out = f.read()

        expected = expected.replace('*PLACEHOLDERPATH*', out_path.parent.as_posix())
        assert out == expected

    @pytest.mark.parametrize(
        'data,design_matrix,comparisons,lrt_factors,expected_path',
        [
            (
                'tests/test_files/big_counted.csv',
                'tests/test_files/test_design_matrix.csv',
                [('condition', 'cond2', 'cond1')],
                ['factor1', 'factor2'],
                'tests/test_files/deseq2_tests/case4/expected_deseq_script_lrt_factors.R',
            ),
        ],
    )
    def test_create_deseq2_script_with_lrt_factors(self, data, design_matrix, comparisons, lrt_factors, expected_path):
        runner = DESeqRunner(data, design_matrix, comparisons, lrt_factors=lrt_factors)
        with open(expected_path) as f:
            expected = f.read()

        out_path = runner.create_script()
        assert Path(out_path).exists()
        with open(out_path) as f:
            out = f.read()

        expected = expected.replace('*PLACEHOLDERPATH*', out_path.parent.as_posix())
        assert out == expected

    @pytest.mark.parametrize(
        'data,design_matrix,comparisons,covariates,lrt_factors,expected_path',
        [
            (
                'tests/test_files/big_counted.csv',
                'tests/test_files/test_design_matrix.csv',
                [('condition', 'cond2', 'cond1')],
                ['covariate1', 'covariate2'],
                ['factor1', 'factor2'],
                'tests/test_files/deseq2_tests/case5/expected_deseq_script_combined.R',
            ),
        ],
    )
    def test_create_deseq2_script_combined(
        self, data, design_matrix, comparisons, covariates, lrt_factors, expected_path
    ):
        runner = DESeqRunner(data, design_matrix, comparisons, covariates=covariates, lrt_factors=lrt_factors)
        with open(expected_path) as f:
            expected = f.read()

        out_path = runner.create_script()
        assert Path(out_path).exists()
        with open(out_path) as f:
            out = f.read()

        expected = expected.replace('*PLACEHOLDERPATH*', out_path.parent.as_posix())
        assert out == expected

    @pytest.mark.parametrize(
        'data,design_matrix,comparisons,scale_factor_path,scale_factors_ndims,expected_path',
        [
            (
                'tests/test_files/big_counted.csv',
                'tests/test_files/test_design_matrix.csv',
                [('condition', 'cond2', 'cond1')],
                'tests/test_files/deseq2_tests/scaling_factors_1d.csv',
                1,
                'tests/test_files/deseq2_tests/case6/expected_deseq_script_normfactors_1d.R',
            ),
            (
                'tests/test_files/big_counted.csv',
                'tests/test_files/test_design_matrix.csv',
                [('condition', 'cond2', 'cond1')],
                'tests/test_files/deseq2_tests/scaling_factors_2d.csv',
                2,
                'tests/test_files/deseq2_tests/case7/expected_deseq_script_normfactors_2d.R',
            ),
        ],
    )
    def test_create_deseq2_script_with_scale_factors(
        self, data, design_matrix, comparisons, scale_factor_path, scale_factors_ndims, expected_path
    ):
        runner = DESeqRunner(
            data,
            design_matrix,
            comparisons,
            scale_factor_path=scale_factor_path,
            scale_factors_ndims=scale_factors_ndims,
        )
        with open(expected_path) as f:
            expected = f.read()

        out_path = runner.create_script()
        assert Path(out_path).exists()
        with open(out_path) as f:
            out = f.read()

        expected = expected.replace('*PLACEHOLDERPATH*', out_path.parent.as_posix())
        assert out == expected

    @pytest.mark.parametrize(
        'data,design_matrix,comparisons,expected_path',
        [
            (
                'tests/test_files/big_counted.csv',
                'tests/test_files/test_design_matrix.csv',
                [('condition', 'cond2', 'cond1')],
                'tests/test_files/deseq2_tests/case1/expected_deseq_script_1.R',
            ),
        ],
    )
    def test_create_deseq2_script_without_scale_factors_strips_normfactors_placeholder(
        self, data, design_matrix, comparisons, expected_path
    ):
        # cooks_cutoff=False matches the golden file used by test_create_deseq2_script's case1
        runner = DESeqRunner(data, design_matrix, comparisons, scale_factor_path=None, cooks_cutoff=False)
        with open(expected_path) as f:
            expected = f.read()

        out_path = runner.create_script()
        assert Path(out_path).exists()
        with open(out_path) as f:
            out = f.read()

        assert '$NORMFACTORS' not in out
        expected = expected.replace('*PLACEHOLDERPATH*', out_path.parent.as_posix())
        assert out == expected


class _RCalled(Exception):
    pass


class TestDiffExpPreflightValidation:
    """
    comparisons/covariates/lrt_factors must be validated against the design matrix before any R process
    starts (issue #280), so a typo fails instantly with a plain-language error instead of after a full R fit.
    """

    DESIGN_MAT = (
        'sample,condition,replicate,age,batch,genotype\n'
        's1,KO,rep1,10,1,mutA\n'
        's2,KO,rep2,15,2,mutB\n'
        's3,KO,rep3,12.5,1,mutA\n'
        's4,WT,rep1,8,2,mutB\n'
        's5,WT,rep2,11,1,mutA\n'
        's6,WT,rep3,9,2,mutB\n'
    )

    @pytest.fixture
    def design_mat_path(self, tmp_path):
        pth = tmp_path.joinpath('design_mat.csv')
        pth.write_text(self.DESIGN_MAT)
        return pth

    @pytest.fixture
    def no_r(self, monkeypatch):
        # any attempt to touch R (installing packages or running a script) means validation came too late
        def fail(*args, **kwargs):
            raise _RCalled('R was invoked')

        monkeypatch.setattr(installs, 'install_deseq2', fail)
        monkeypatch.setattr(installs, 'install_limma', fail)
        monkeypatch.setattr(io, 'run_r_script', fail)

    @pytest.mark.parametrize('runner_cls', [DESeqRunner, LimmaVoomRunner])
    @pytest.mark.parametrize(
        'comparisons,covariates,lrt_factors,expected_in_msg',
        [
            # comparison factor is not a design-matrix column
            ([('conditon', 'KO', 'WT')], [], [], ["'conditon'", "'condition'"]),
            # numerator is not an observed value of the factor (case typo)
            ([('condition', 'Ko', 'WT')], [], [], ["'Ko'", "'KO'", "'WT'"]),
            # denominator has a trailing space
            ([('condition', 'KO', 'WT ')], [], [], ["'WT '", "'KO'", "'WT'"]),
            # denominator belongs to a different factor
            ([('condition', 'KO', 'rep1')], [], [], ["'rep1'"]),
            # numerator == denominator
            ([('condition', 'KO', 'KO')], [], [], ["'KO'"]),
            # comparisons on a numeric column (R treats it as continuous, so it has no levels to compare)
            ([('batch', '2', '1')], [], [], ["'batch'", 'numeric']),
            ([('batch', 2, 1)], [], [], ["'batch'", 'numeric']),
            # malformed comparison
            ([('condition', 'KO')], [], [], ['three']),
            # covariate is not a column
            ([('condition', 'KO', 'WT')], ['agee'], [], ["'agee'", "'age'"]),
            # covariate is not numeric
            ([('condition', 'KO', 'WT')], ['genotype'], [], ["'genotype'", 'numeric']),
            # LRT factor is not a column
            ([('condition', 'KO', 'WT')], [], ['genotyp'], ["'genotyp'", "'genotype'"]),
            # one part of an interaction term is not a column
            ([('condition', 'KO', 'WT')], [], ['condition:genotyp'], ["'genotyp'", "'condition:genotyp'"]),
            # polynomial LRT term on a factor that is not a column
            ([('condition', 'KO', 'WT')], [], ['poly(agee, degree = 2)'], ["'agee'"]),
        ],
    )
    def test_invalid_inputs_fail_before_r(
        self, design_mat_path, no_r, runner_cls, comparisons, covariates, lrt_factors, expected_in_msg
    ):
        runner = runner_cls('counts.csv', design_mat_path, comparisons, covariates=covariates, lrt_factors=lrt_factors)
        with pytest.raises(InvalidValueError) as e:
            runner.run()
        for expected in expected_in_msg:
            assert expected in str(e.value)

    @pytest.mark.parametrize('runner_cls', [DESeqRunner, LimmaVoomRunner])
    @pytest.mark.parametrize(
        'comparisons,covariates,lrt_factors',
        [
            ([('condition', 'KO', 'WT')], [], []),
            ([('condition', 'WT', 'KO'), ('replicate', 'rep2', 'rep1')], [], []),
            ([('genotype', 'mutB', 'mutA')], ['age'], []),
            ([('condition', 'KO', 'WT')], ['age', 'batch'], ['genotype']),
            ([('condition', 'KO', 'WT')], [], ['condition:genotype']),
            ([('condition', 'KO', 'WT')], [], ['poly(age, degree = 2)']),
            ([], [], []),
        ],
    )
    def test_valid_inputs_pass_validation(
        self, design_mat_path, no_r, runner_cls, comparisons, covariates, lrt_factors
    ):
        runner = runner_cls('counts.csv', design_mat_path, comparisons, covariates=covariates, lrt_factors=lrt_factors)
        # validation passes, so the run proceeds to the (mocked) R step
        with pytest.raises(_RCalled):
            runner.run()


def _r_string_literals(script: str) -> List[List[str]]:
    """
    Minimal lexer for the R scripts generated by RNAlysis: return the (unescaped) double/single-quoted string \
    literals on each line, ignoring comments. Raise AssertionError if a string literal is left unterminated.
    """
    literals_per_line = []
    for line_num, line in enumerate(script.splitlines()):
        literals = []
        i = 0
        while i < len(line):
            char = line[i]
            if char == '#':
                break
            if char in '"\'':
                quote = char
                i += 1
                current = ''
                while True:
                    assert i < len(line), f'unterminated R string literal on line {line_num + 1}: {line}'
                    if line[i] == '\\':
                        current += line[i + 1]
                        i += 2
                        continue
                    if line[i] == quote:
                        break
                    current += line[i]
                    i += 1
                literals.append(current)
            i += 1
        literals_per_line.append(literals)
    return literals_per_line


def _assert_r_parses(script_path: Path):
    rscript = shutil.which('Rscript')
    if rscript is None:
        return
    proc = subprocess.run(
        [rscript, '-e', f'invisible(parse(file="{Path(script_path).as_posix()}", encoding="UTF-8"))'],
        capture_output=True,
    )
    assert proc.returncode == 0, proc.stderr.decode('utf8', errors='replace')


def _read_script(script_path: Path) -> str:
    with open(script_path, encoding='utf-8') as f:
        return f.read()


def _written_files(script: str) -> List[Path]:
    paths = []
    for line, literals in zip(script.splitlines(), _r_string_literals(script)):
        if line.startswith('write.csv('):
            assert len(literals) == 1
            paths.append(Path(literals[0]))
    return paths


class TestHostileDesignValues:
    """
    issue #281: design-matrix values are spliced into R source and into output file names. \
    Values containing quotes, backslashes, path separators, or other characters that are illegal in file names \
    must produce valid R source and filename-safe output names, while safe values keep their exact names.
    """

    HOSTILE = {
        # value: expected file-name component
        'KO/1': 'KO_1',
        'a:b': 'a_b',
        'back\\slash': 'back_slash',
        'q"uote': 'q_uote',
        "ap'os": "ap'os",
        'wild type': 'wild type',
        'Wildtyp-ü': 'Wildtyp-ü',
    }
    # R's make.names() of the limma coefficient name for each value
    LIMMA_COEFS = {
        'a:b': 'conditiona.b',
        'back\\slash': 'conditionback.slash',
        'q"uote': 'conditionq.uote',
        "ap'os": 'conditionap.os',
        'wild type': 'conditionwild.type',
        'Wildtyp-ü': 'conditionWildtyp.ü',
    }

    @pytest.fixture
    def hostile_design_path(self, tmp_path):
        values = list(self.HOSTILE) + ['WT']
        df = pl.DataFrame(
            {
                'sample': [f's{i}' for i in range(len(values))],
                'condition': values,
                'replicate': [f'rep{i % 2}' for i in range(len(values))],
            }
        )
        pth = tmp_path.joinpath('design.csv')
        df.write_csv(pth)
        return pth

    @pytest.fixture
    def no_r(self, monkeypatch):
        def fail(*args, **kwargs):
            raise _RCalled('R was invoked')

        monkeypatch.setattr(installs, 'install_deseq2', fail)
        monkeypatch.setattr(installs, 'install_limma', fail)
        monkeypatch.setattr(io, 'run_r_script', fail)

    def test_deseq2_script_hostile_values(self, hostile_design_path):
        comparisons = [('condition', value, 'WT') for value in self.HOSTILE]
        runner = DESeqRunner('counts.csv', hostile_design_path, comparisons)
        runner.validate_design_inputs()
        script_path = runner.create_script()
        script = _read_script(script_path)
        literals = _r_string_literals(script)  # raises if any string literal is malformed

        contrasts = [lits for line, lits in zip(script.splitlines(), literals) if 'contrast=c(' in line]
        assert contrasts == [list(comp) for comp in comparisons]

        written = _written_files(script)
        assert [pth.parent for pth in written] == [script_path.parent] * len(comparisons)
        assert [pth.name for pth in written] == [f'DESeq2_condition_{name}_vs_WT.csv' for name in self.HOSTILE.values()]
        _assert_r_parses(script_path)

    def test_limma_script_hostile_values(self, hostile_design_path):
        comparisons = [('condition', value, 'WT') for value in self.HOSTILE]
        runner = LimmaVoomRunner('counts.csv', hostile_design_path, comparisons)
        runner.validate_design_inputs()
        script_path = runner.create_script()
        script = _read_script(script_path)
        literals = _r_string_literals(script)

        levels = [lits for line, lits in zip(script.splitlines(), literals) if 'design_matrix$condition <-' in line]
        assert levels == [sorted(list(self.HOSTILE) + ['WT'])]

        contrasts = [lits for line, lits in zip(script.splitlines(), literals) if 'makeContrasts(' in line]
        # 'KO/1' is the baseline level (alphabetically first), so it has no coefficient of its own
        expected_contrasts = [[' - conditionWT']] + [
            [f'{self.LIMMA_COEFS[value]} - conditionWT'] for value in list(self.HOSTILE)[1:]
        ]
        assert contrasts == expected_contrasts

        written = _written_files(script)
        assert [pth.parent for pth in written] == [script_path.parent] * len(comparisons)
        assert [pth.name for pth in written] == [
            f'LimmaVoom_condition_{name}_vs_WT.csv' for name in self.HOSTILE.values()
        ]
        _assert_r_parses(script_path)

    @pytest.mark.parametrize(
        'runner_cls,expected_names',
        [
            (
                DESeqRunner,
                [
                    'DESeq2_condition_cond2_vs_cond1.csv',
                    'DESeq2_factor1_B_vs_A.csv',
                    'DESeq2_covariate1_covariate.csv',
                    'DESeq2_factor1_x_factor2_LRT.csv',
                    'DESeq2_poly(covariate2, degree = 2)_LRT.csv',
                ],
            ),
            (
                LimmaVoomRunner,
                [
                    'LimmaVoom_condition_cond2_vs_cond1.csv',
                    'LimmaVoom_factor1_B_vs_A.csv',
                    'LimmaVoom_covariate1_covariate.csv',
                    'LimmaVoom_factor1Xfactor2_LRT.csv',
                    'LimmaVoom_poly(covariate2, degree = 2)_LRT.csv',
                ],
            ),
        ],
    )
    def test_safe_values_keep_their_file_names(self, runner_cls, expected_names):
        # rule 5: names that were already filename-safe must not change (byte-identical to RNAlysis 4.3.0)
        runner = runner_cls(
            'counts.csv',
            'tests/test_files/test_design_matrix_advanced.csv',
            [('condition', 'cond2', 'cond1'), ('factor1', 'B', 'A')],
            covariates=['covariate1'],
            lrt_factors=['factor1:factor2', 'poly(covariate2, degree = 2)'],
        )
        script_path = runner.create_script()
        names = [pth.name for pth in _written_files(_read_script(script_path))]
        assert names == expected_names

    @pytest.mark.parametrize('runner_cls', [DESeqRunner, LimmaVoomRunner])
    @pytest.mark.parametrize(
        'values,comparisons,expected_in_msg',
        [
            # two distinct values that sanitize to the same file name
            (['KO/1', 'KO:1', 'WT'], [('condition', 'KO/1', 'WT'), ('condition', 'KO:1', 'WT')], ["'KO/1'", "'KO:1'"]),
            # names that differ only by case collide on case-insensitive file systems (Windows, macOS)
            (['KO', 'ko', 'WT'], [('condition', 'KO', 'WT'), ('condition', 'ko', 'WT')], ["'KO'", "'ko'"]),
        ],
    )
    def test_file_name_collisions_fail_before_r(self, tmp_path, no_r, runner_cls, values, comparisons, expected_in_msg):
        pth = tmp_path.joinpath('design.csv')
        pl.DataFrame({'sample': [f's{i}' for i in range(len(values))], 'condition': values}).write_csv(pth)
        runner = runner_cls('counts.csv', pth, comparisons)
        with pytest.raises(InvalidValueError) as e:
            runner.run()
        for expected in expected_in_msg:
            assert expected in str(e.value)
        # the same file names are also refused when the script is created directly
        with pytest.raises(InvalidValueError):
            runner.create_script()

    def test_duplicate_comparison_is_not_a_collision(self, tmp_path):
        # listing the very same comparison twice already worked before #281 (it writes the same file twice)
        pth = tmp_path.joinpath('design.csv')
        pl.DataFrame({'sample': ['s0', 's1'], 'condition': ['KO', 'WT']}).write_csv(pth)
        runner = DESeqRunner('counts.csv', pth, [('condition', 'KO', 'WT'), ('condition', 'KO', 'WT')])
        runner.validate_design_inputs()
        runner.create_script()

    def test_distinct_names_that_casefold_alike_are_not_a_collision(self, tmp_path):
        # 'ß' and 'ss' are different file names even on case-insensitive file systems
        pth = tmp_path.joinpath('design.csv')
        pl.DataFrame({'sample': ['s0', 's1', 's2'], 'condition': ['Straße', 'Strasse', 'WT']}).write_csv(pth)
        runner = DESeqRunner('counts.csv', pth, [('condition', 'Straße', 'WT'), ('condition', 'Strasse', 'WT')])
        runner.validate_design_inputs()

    @pytest.mark.parametrize(
        'name,expected',
        [
            ('conditionKO-1', 'conditionKO.1'),
            ('conditionwild type', 'conditionwild.type'),
            ('cell.typecond2', 'cell.typecond2'),
            ('conditionü', 'conditionü'),
            ('condition²', 'condition.'),
            ('_x', 'X_x'),
            ('.1x', 'X.1x'),
        ],
    )
    def test_r_syntactic_name(self, name, expected):
        # mirrors R's make.names() (checked against R 4.4.2)
        assert LimmaVoomRunner._r_syntactic_name(name) == expected

    def test_limma_contrast_with_dotted_factor_name(self, tmp_path):
        # model.matrix() names coefficients after the factor as written in the formula ('cell.type'),
        # so the contrast must too (the output file name keeps its established sanitized spelling)
        pth = tmp_path.joinpath('design.csv')
        pl.DataFrame({'sample': ['s0', 's1', 's2'], 'cell.type': ['cond1', 'cond2', 'cond3']}).write_csv(pth)
        comparisons = [('cell.type', 'cond2', 'cond1'), ('cell.type', 'cond1', 'cond3')]
        runner = LimmaVoomRunner('counts.csv', pth, comparisons)
        runner.validate_design_inputs()
        script = _read_script(runner.create_script())
        literals = _r_string_literals(script)
        contrasts = [lits for line, lits in zip(script.splitlines(), literals) if 'makeContrasts(' in line]
        assert contrasts == [['cell.typecond2'], [' - cell.typecond3']]
        assert [pth.name for pth in _written_files(script)] == [
            'LimmaVoom_cell_type_cond2_vs_cond1.csv',
            'LimmaVoom_cell_type_cond1_vs_cond3.csv',
        ]

    def test_limma_coefficient_collision_fails_before_r(self, tmp_path, no_r):
        # 'KO.1' and 'KO-1' are distinct file names, but R's make.names() maps both limma coefficients
        # to 'conditionKO.1', so the contrast would be ambiguous
        pth = tmp_path.joinpath('design.csv')
        pl.DataFrame({'sample': ['s0', 's1', 's2'], 'condition': ['KO-1', 'KO.1', 'WT']}).write_csv(pth)
        runner = LimmaVoomRunner('counts.csv', pth, [('condition', 'KO-1', 'WT')])
        with pytest.raises(InvalidValueError) as e:
            runner.run()
        assert "'KO-1'" in str(e.value) and "'KO.1'" in str(e.value)

    @staticmethod
    def _hostile_big_design(tmp_path) -> Path:
        # the standard test design matrix, with 'cond1' and 'cond2' renamed to values that contain path separators
        # ('wt-2/x' is also renamed by R's make.names(), so it exercises limma's coefficient naming)
        design = io.load_table('tests/test_files/test_design_matrix.csv')
        design = design.with_columns(pl.col('condition').replace({'cond1': 'KO/1', 'cond2': 'wt-2/x'}))
        pth = tmp_path.joinpath('design.csv')
        design.write_csv(pth)
        return pth

    @staticmethod
    def _assert_matches_truth(out_path: Path, truth_path: str):
        out = io.load_table(out_path).sort(pl.first())
        truth = io.load_table(truth_path).sort(pl.first())
        assert out.shape == truth.shape
        assert out.columns == truth.columns
        assert out.select(pl.first()).equals(truth.select(pl.first()))
        if sys.platform == 'win32':  # running DESeq in linux gives slightly different results
            assert np.allclose(out.drop(cs.first()), truth.drop(cs.first()), equal_nan=True, atol=1 * 10 ** (-4))

    def test_run_deseq2_hostile_value(self, tmp_path):
        runner = DESeqRunner(
            'tests/test_files/big_counted.csv', self._hostile_big_design(tmp_path), [('condition', 'wt-2/x', 'KO/1')]
        )
        output_dir = runner.run()
        out_path = output_dir.joinpath('DESeq2_condition_wt-2_x_vs_KO_1.csv')
        assert out_path.exists()
        # renaming the values must not change the results
        self._assert_matches_truth(
            out_path, 'tests/test_files/deseq2_tests/case2/DESeq2_condition_cond2_vs_cond1_truth.csv'
        )

    def test_run_limma_hostile_value(self, tmp_path):
        runner = LimmaVoomRunner(
            'tests/test_files/big_counted.csv', self._hostile_big_design(tmp_path), [('condition', 'wt-2/x', 'KO/1')]
        )
        output_dir = runner.run()
        out_path = output_dir.joinpath('LimmaVoom_condition_wt-2_x_vs_KO_1.csv')
        assert out_path.exists()
        self._assert_matches_truth(
            out_path, 'tests/test_files/limma_tests/case2/LimmaVoom_condition_cond2_vs_cond1_truth.csv'
        )
