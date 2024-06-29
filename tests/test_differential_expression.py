import sys

import numpy as np
import polars.selectors as cs
import pytest

from rnalysis.utils.differential_expression import *


class TestLimmaVoomRunner:
    @pytest.mark.parametrize("data,design_matrix,comparisons,random_effect,expected_path", [
        ('tests/test_files/big_counted.csv', 'tests/test_files/test_design_matrix.csv',
         [('condition', 'cond2', 'cond1')],
         None, 'tests/test_files/limma_tests/case1/expected_limma_script_1.R'),
        ('counted.csv', 'tests/test_files/test_design_matrix.csv',
         [('condition', 'cond3', 'cond2'), ('replicate', 'rep2', 'rep1'), ('condition', 'cond1', 'cond2')],
         None, 'tests/test_files/limma_tests/case2/expected_limma_script_2.R'),
        ('tests/test_files/big_counted.csv', 'tests/test_files/test_design_matrix.csv',
         [('condition', 'cond2', 'cond1')],
         'replicate', 'tests/test_files/limma_tests/case3/expected_limma_script_3.R'),
    ])
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

    @pytest.mark.parametrize('comparisons,expected_paths', [
        (
            [('replicate', 'rep2', 'rep3')],
            ['tests/test_files/limma_tests/case1/LimmaVoom_replicate_rep2_vs_rep3_truth.csv']),
        ([('condition', 'cond2', 'cond1'), ('condition', 'cond3', 'cond2')],
         ['tests/test_files/limma_tests/case2/LimmaVoom_condition_cond2_vs_cond1_truth.csv',
          'tests/test_files/limma_tests/case2/LimmaVoom_condition_cond3_vs_cond2_truth.csv'])
    ])
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
                assert np.allclose(out.drop(cs.first()), truth.drop(cs.first()), equal_nan=True, atol=1 * 10 ** (- 4))

    @pytest.mark.parametrize("data,design_matrix,comparisons,covariates,expected_path", [
        ('tests/test_files/big_counted.csv', 'tests/test_files/test_design_matrix_advanced.csv',
         [('condition', 'cond2', 'cond1')], ['covariate1', 'covariate2'],
         'tests/test_files/limma_tests/case4/expected_limma_script_covariates.R'),
    ])
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

    @pytest.mark.parametrize("data,design_matrix,comparisons,lrt_factors,expected_path", [
        ('tests/test_files/big_counted.csv', 'tests/test_files/test_design_matrix_advanced.csv',
         [('condition', 'cond2', 'cond1')], ['factor1', 'factor2'],
         'tests/test_files/limma_tests/case5/expected_limma_script_lrt_factors.R'),
    ])
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

    @pytest.mark.parametrize("data,design_matrix,comparisons,covariates,lrt_factors,expected_path", [
        ('tests/test_files/big_counted.csv', 'tests/test_files/test_design_matrix_advanced.csv',
         [('condition', 'cond2', 'cond1')], ['covariate1', 'covariate2'], ['factor1', 'factor2'],
         'tests/test_files/limma_tests/case6/expected_limma_script_combined.R'),
    ])
    def test_create_limma_script_combined(self, data, design_matrix, comparisons, covariates, lrt_factors,
                                          expected_path):
        with open(expected_path) as f:
            expected = f.read()
        runner = LimmaVoomRunner(data, design_matrix, comparisons, covariates=covariates, lrt_factors=lrt_factors)
        out_path = runner.create_script()
        assert Path(out_path).exists()
        with open(out_path) as f:
            out = f.read()

        expected = expected.replace('*PLACEHOLDERPATH*', out_path.parent.as_posix())
        assert out == expected

    @pytest.mark.parametrize("model_factors,expected_coef_names", [
        (['factor1', 'factor2', 'factor3'],
         {'factor1': {'factor1B', 'factor1C'},
          'factor2': {'factor2'},
          'factor3': {'factor3Y', 'factor3Z'}}),
        (['factor1', 'poly(factor2, degree = 2)', 'factor3'],
         {'factor1': {'factor1B', 'factor1C'},
          'poly(factor2, degree = 2)': {'poly.factor2..degree...2.1', 'poly.factor2..degree...2.2'},
          'factor3': {'factor3Y', 'factor3Z'}}),
        (['factor1', 'factor2', 'factor1:factor2'],
         {'factor1': {'factor1B', 'factor1C'},
          'factor2': {'factor2'},
          'factor1:factor2': {'factor1B.factor2', 'factor1C.factor2'}}),
        (['factor1', 'factor2', 'factor3', 'factor1:factor2:factor3'],
         {'factor1': {'factor1B', 'factor1C'},
          'factor2': {'factor2'},
          'factor3': {'factor3Y', 'factor3Z'},
          'factor1:factor2:factor3': {'factor1B.factor2.factor3Y', 'factor1B.factor2.factor3Z',
                                      'factor1C.factor2.factor3Y', 'factor1C.factor2.factor3Z'}}),
    ])
    def test_get_coef_names(self, model_factors, expected_coef_names):
        design_matrix = 'tests/test_files/limma_tests/test_design_mat.csv'
        runner = LimmaVoomRunner('data.csv', design_matrix, [('factor1', 'B', 'A')],
                                 model_factors=model_factors)
        coef_names = runner._get_coef_names()
        assert coef_names == expected_coef_names

    @pytest.mark.parametrize("interaction,expected_terms", [
        ('factor1:factor2', {'factor1B:factor2', 'factor1C:factor2'}),
        ('factor1:factor3', {'factor1B:factor3Y', 'factor1B:factor3Z', 'factor1C:factor3Y', 'factor1C:factor3Z'}),
        ('factor2:factor3', {'factor2:factor3Y', 'factor2:factor3Z'}),
        ('factor1:factor2:factor3', {'factor1B:factor2:factor3Y', 'factor1B:factor2:factor3Z',
                                     'factor1C:factor2:factor3Y', 'factor1C:factor2:factor3Z'}),
    ])
    def test_generate_interaction_terms(self, interaction, expected_terms):
        design_matrix = 'tests/test_files/limma_tests/test_design_mat.csv'
        runner = LimmaVoomRunner('data.csv', design_matrix, [('factor1', 'B', 'A')])
        interaction_terms = runner._generate_interaction_terms(interaction)
        assert interaction_terms == expected_terms


class TestDESeqRunner:
    @pytest.mark.parametrize("data,design_matrix,comparisons,cooks_cutoff,expected_path", [
        ('tests/test_files/big_counted.csv', 'tests/test_files/test_design_matrix.csv',
         [('condition', 'cond2', 'cond1')], False,
         'tests/test_files/deseq2_tests/case1/expected_deseq_script_1.R'),
        ('counted.csv', 'tests/test_files/test_design_matrix.csv',
         [('condition', 'cond3', 'cond2'), ('replicate', 'rep2', 'rep1'), ('condition', 'cond1', 'cond2')], True,
         'tests/test_files/deseq2_tests/case2/expected_deseq_script_2.R')
    ])
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

    @pytest.mark.parametrize('comparisons,expected_paths', [
        ([('replicate', 'rep2', 'rep3')],
         ['tests/test_files/deseq2_tests/case1/DESeq2_replicate_rep2_vs_rep3_truth.csv']),
        ([('condition', 'cond2', 'cond1'), ('condition', 'cond3', 'cond2')],
         ['tests/test_files/deseq2_tests/case2/DESeq2_condition_cond2_vs_cond1_truth.csv',
          'tests/test_files/deseq2_tests/case2/DESeq2_condition_cond3_vs_cond2_truth.csv'])
    ])
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
                assert np.allclose(out.drop(cs.first()), truth.drop(cs.first()), equal_nan=True, atol=1 * 10 ** (- 4))

    @pytest.mark.parametrize("data,design_matrix,comparisons,covariates,expected_path", [
        ('tests/test_files/big_counted.csv', 'tests/test_files/test_design_matrix.csv',
         [('condition', 'cond2', 'cond1')], ['covariate1', 'covariate2'],
         'tests/test_files/deseq2_tests/case3/expected_deseq_script_covariates.R'),
    ])
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

    @pytest.mark.parametrize("data,design_matrix,comparisons,lrt_factors,expected_path", [
        ('tests/test_files/big_counted.csv', 'tests/test_files/test_design_matrix.csv',
         [('condition', 'cond2', 'cond1')], ['factor1', 'factor2'],
         'tests/test_files/deseq2_tests/case4/expected_deseq_script_lrt_factors.R'),
    ])
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

    @pytest.mark.parametrize("data,design_matrix,comparisons,covariates,lrt_factors,expected_path", [
        ('tests/test_files/big_counted.csv', 'tests/test_files/test_design_matrix.csv',
         [('condition', 'cond2', 'cond1')], ['covariate1', 'covariate2'], ['factor1', 'factor2'],
         'tests/test_files/deseq2_tests/case5/expected_deseq_script_combined.R'),
    ])
    def test_create_deseq2_script_combined(self, data, design_matrix, comparisons, covariates, lrt_factors,
                                           expected_path):
        runner = DESeqRunner(data, design_matrix, comparisons, covariates=covariates, lrt_factors=lrt_factors)
        with open(expected_path) as f:
            expected = f.read()

        out_path = runner.create_script()
        assert Path(out_path).exists()
        with open(out_path) as f:
            out = f.read()

        expected = expected.replace('*PLACEHOLDERPATH*', out_path.parent.as_posix())
        assert out == expected
