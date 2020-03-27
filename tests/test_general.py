import pytest
import numpy as np
import pandas as pd
from pathlib import Path
from rnalysis.general import *
from rnalysis.general import _check_is_df,_remove_unindexed_rows


def test_is_df_dataframe():
    my_df = pd.DataFrame()
    assert (_check_is_df(my_df))


def test_is_df_str():
    correct_path = "myfile.csv"
    correct_path_2 = r"myfolder\anotherfile.csv"
    assert (not _check_is_df(correct_path))
    assert (not _check_is_df(correct_path_2))


def test_is_df_pathobj():
    correct_pathobj = Path("all_feature_96_new.csv")
    assert (not _check_is_df(correct_pathobj))


def test_is_df_str_notcsv():
    incorrect_pth = "myfile.xlsx"
    incorrect_pth2 = "all_feature_96_new"
    with pytest.raises(ValueError):
        _check_is_df(incorrect_pth)
    with pytest.raises(ValueError):
        _check_is_df(incorrect_pth2)


def test_is_df_pathobj_notcsv():
    incorrect_pathobj = Path("test_general.py")
    with pytest.raises(ValueError):
        _check_is_df(incorrect_pathobj)


def test_is_df_invalid_type():
    invalid_type = 67
    with pytest.raises(ValueError):
        _check_is_df(invalid_type)


def test_load_csv_bad_input():
    invalid_input = 2349
    with pytest.raises(AssertionError):
        a = load_csv(invalid_input)


def test_load_csv():
    truth = pd.DataFrame({'idxcol': ['one', 'two', 'three'], 'othercol': [4, 5, 6]})
    truth.set_index('idxcol', inplace=True)
    pth = "test_load_csv.csv"
    assert (truth == load_csv(pth, 0)).all().all()


def test_remove_unindexed_rows():
    truth = load_csv("counted_missing_rows_deleted.csv", 0)
    missing = load_csv("counted_missing_rows.csv", 0)
    assert (truth == _remove_unindexed_rows(missing)).all().all()


def test_parse_wbgene_string():
    string = ' WBGene WBGenes WBGene12345678, WBGene98765432WBGene00000000\n the geneWBGene44444444 ' \
             'daf-16 A5gHB.5 /// WBGene55555555'
    truth = {'WBGene12345678', 'WBGene98765432', 'WBGene00000000', 'WBGene44444444', 'WBGene55555555'}
    assert truth == parse_wbgene_string(string)


def test_parse_sequence_name_string():
    string = 'CELE_Y55D5A.5T23G5.6 /// WBGene00000000 daf-16\nZK662.4 '
    truth = {'Y55D5A.5', 'T23G5.6', 'ZK662.4'}
    assert truth == parse_sequence_name_string(string)


def test_parse_gene_name_string():
    string = 'saeg-2 \\\ lin-15B cyp-23A1lin-15A WBGene12345678\n GHF5H.3'
    truth = {'saeg-2', 'lin-15B', 'cyp-23A1', 'lin-15A'}
    assert truth == parse_gene_name_string(string)
