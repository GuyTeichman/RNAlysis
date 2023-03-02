import platform
import sys

import numpy as np
import pytest

from rnalysis.utils.differential_expression import *


def test_install_limma():
    install_limma()


@pytest.mark.parametrize("data,design_matrix,comparisons,expected_path", [
    ('tests/test_files/big_counted.csv', 'tests/test_files/test_design_matrix.csv', [('condition', 'cond2', 'cond1')],
     'tests/test_files/limma_tests/case1/expected_limma_script_1.R'),
    ('counted.csv', 'tests/test_files/test_design_matrix.csv',
     [('condition', 'cond3', 'cond2'), ('replicate', 'rep2', 'rep1'), ('condition', 'cond1', 'cond2')],
     'tests/test_files/limma_tests/case2/expected_limma_script_2.R')
])
def test_create_limma_script(data, design_matrix, comparisons, expected_path):
    with open(expected_path) as f:
        expected = f.read()

    out_path = create_limma_script(data, design_matrix, comparisons)
    assert Path(out_path).exists()
    with open(out_path) as f:
        out = f.read()

    expected = expected.replace('*PLACEHOLDERPATH*', out_path.parent.as_posix())
    assert out == expected


@pytest.mark.parametrize('comparisons,expected_paths', [
    (
    [('replicate', 'rep2', 'rep3')], ['tests/test_files/limma_tests/case1/LimmaVoom_replicate_rep2_vs_rep3_truth.csv']),
    ([('condition', 'cond2', 'cond1'), ('condition', 'cond3', 'cond2')],
     ['tests/test_files/limma_tests/case2/LimmaVoom_condition_cond2_vs_cond1_truth.csv',
      'tests/test_files/limma_tests/case2/LimmaVoom_condition_cond3_vs_cond2_truth.csv'])
])
def test_run_limma_analysis(comparisons, expected_paths):
    data_path = 'tests/test_files/big_counted.csv'
    design_mat_path = 'tests/test_files/test_design_matrix.csv'
    output_dir = run_limma_analysis(data_path, design_mat_path, comparisons)

    dfs = []
    for item in output_dir.iterdir():
        if item.is_file() and item.suffix == '.csv':
            dfs.append(io.load_csv(item, index_col=0))

    expected_dfs = []
    for item in expected_paths:
        expected_dfs.append(io.load_csv(item, index_col=0))

    for out, truth in zip(dfs, expected_dfs):
        assert out.shape == truth.shape
        assert np.all(out.columns == truth.columns)
        assert np.all(sorted(out.index) == sorted(truth.index))
        if sys.platform == 'win32':  # running DESeq in linux gives slightly different results
            assert np.allclose(out, truth, equal_nan=True, atol=1 * 10 ** (- 5))


def test_install_deseq2():
    install_deseq2()


@pytest.mark.parametrize("data,design_matrix,comparisons,expected_path", [
    ('tests/test_files/big_counted.csv', 'tests/test_files/test_design_matrix.csv', [('condition', 'cond2', 'cond1')],
     'tests/test_files/deseq2_tests/case1/expected_deseq_script_1.R'),
    ('counted.csv', 'tests/test_files/test_design_matrix.csv',
     [('condition', 'cond3', 'cond2'), ('replicate', 'rep2', 'rep1'), ('condition', 'cond1', 'cond2')],
     'tests/test_files/deseq2_tests/case2/expected_deseq_script_2.R')
])
def test_create_deseq2_script(data, design_matrix, comparisons, expected_path):
    with open(expected_path) as f:
        expected = f.read()

    out_path = create_deseq2_script(data, design_matrix, comparisons)
    assert Path(out_path).exists()
    with open(out_path) as f:
        out = f.read()

    expected = expected.replace('*PLACEHOLDERPATH*', out_path.parent.as_posix())
    assert out == expected


@pytest.mark.parametrize('comparisons,expected_paths', [
    ([('replicate', 'rep2', 'rep3')], ['tests/test_files/deseq2_tests/case1/DESeq2_replicate_rep2_vs_rep3_truth.csv']),
    ([('condition', 'cond2', 'cond1'), ('condition', 'cond3', 'cond2')],
     ['tests/test_files/deseq2_tests/case2/DESeq2_condition_cond2_vs_cond1_truth.csv',
      'tests/test_files/deseq2_tests/case2/DESeq2_condition_cond3_vs_cond2_truth.csv'])
])
def test_run_deseq2_analysis(comparisons, expected_paths):
    data_path = 'tests/test_files/big_counted.csv'
    design_mat_path = 'tests/test_files/test_design_matrix.csv'
    output_dir = run_deseq2_analysis(data_path, design_mat_path, comparisons)

    dfs = []
    for item in output_dir.iterdir():
        if item.is_file() and item.suffix == '.csv':
            dfs.append(io.load_csv(item, index_col=0))

    expected_dfs = []
    for item in expected_paths:
        expected_dfs.append(io.load_csv(item, index_col=0))

    for out, truth in zip(dfs, expected_dfs):
        assert out.shape == truth.shape
        assert np.all(out.columns == truth.columns)
        assert np.all(sorted(out.index) == sorted(truth.index))
        if sys.platform == 'win32':  # running DESeq in linux gives slightly different results
            assert np.allclose(out, truth, equal_nan=True, atol=1 * 10 ** (- 5))
