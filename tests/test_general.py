import pytest
import numpy as np
import pandas as pd
from pathlib import Path
from rnalysis.general import *


def test_is_df_dataframe():
    my_df = pd.DataFrame()
    assert (check_is_df(my_df))


def test_is_df_str():
    correct_path = "myfile.csv"
    correct_path_2 = r"myfolder\anotherfile.csv"
    assert (not check_is_df(correct_path))
    assert (not check_is_df(correct_path_2))


def test_is_df_pathobj():
    correct_pathobj = Path("all_feature_96_new.csv")
    assert (not check_is_df(correct_pathobj))


def test_is_df_str_notcsv():
    incorrect_pth = "myfile.xlsx"
    incorrect_pth2 = "all_feature_96_new"
    with pytest.raises(ValueError):
        check_is_df(incorrect_pth)
    with pytest.raises(ValueError):
        check_is_df(incorrect_pth2)


def test_is_df_pathobj_notcsv():
    incorrect_pathobj = Path("test_general.py")
    with pytest.raises(ValueError):
        check_is_df(incorrect_pathobj)


def test_is_df_invalid_type():
    invalid_type = 67
    with pytest.raises(ValueError):
        check_is_df(invalid_type)


def test_load_csv_bad_input():
    invalid_input = 2349
    with pytest.raises(AssertionError):
        a = load_csv(invalid_input)


def test_load_csv():
    truth = pd.DataFrame({'idxcol': ['one', 'two', 'three'], 'othercol': [4, 5, 6]})
    truth.set_index('idxcol', inplace=True)
    pth = "test_load_csv.csv"
    assert (truth == load_csv(pth, 0)).all().all()


def test_norm_reads_to_rpm():
    truth = load_csv(r"test_norm_reads_rpm.csv", 0)
    expr = load_csv(r"all_expr.csv", 0)
    feature = load_csv("all_feature.csv", 0)
    norm = norm_reads_to_rpm(expr, feature)
    assert np.isclose(truth, norm).all()


def test_remove_unindexed_rows():
    truth = load_csv("all_expr_missing_rows_deleted.csv", 0)
    missing = load_csv("all_expr_missing_rows.csv", 0)
    assert (truth == remove_unindexed_rows(missing)).all().all()


def test_filter_low_rpm():
    truth = load_csv("all_expr_low_rpm_truth.csv", 0)
    lowrpm = load_csv("all_expr_low_rpm.csv", 0)
    assert np.isclose(truth, filter_low_rpm(lowrpm, threshold=5)).all()
