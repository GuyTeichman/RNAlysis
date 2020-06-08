import pytest
import os
from rnalysis.general import *
from rnalysis import __biotype_file_key__, __attr_file_key__


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
    string = 'saeg-2 /// \tlin-15B cyp-23A1lin-15A WBGene12345678\n GHF5H.3'
    truth = {'saeg-2', 'lin-15B', 'cyp-23A1', 'lin-15A'}
    assert truth == parse_gene_name_string(string)


def test_save_to_csv():
    assert False


def test_print_settings_file(capfd):
    utils.make_temp_copy_of_settings_file()
    set_biotype_ref_table_path('temp_path')
    set_attr_ref_table_path('temp_path')

    print_settings_file()

    utils.set_temp_copy_of_settings_file_as_default()
    utils.remove_temp_copy_of_settings_file()


def test_set_attr_ref_path():
    utils.make_temp_copy_of_settings_file()
    if utils.get_settings_file_path().exists():
        utils.get_settings_file_path().unlink()

    set_attr_ref_table_path('old/path')
    success_1 = utils.get_attr_ref_path('predefined') == 'old/path'

    set_attr_ref_table_path('new/path')
    success_2 = utils.get_attr_ref_path('predefined') == 'new/path'

    utils.set_temp_copy_of_settings_file_as_default()
    utils.remove_temp_copy_of_settings_file()
    assert success_1
    assert success_2


def test_set_biotype_ref_path():
    utils.make_temp_copy_of_settings_file()
    if utils.get_settings_file_path().exists():
        utils.get_settings_file_path().unlink()

    set_biotype_ref_table_path('old/path')
    success_1 = utils.get_biotype_ref_path('predefined') == 'old/path'

    set_biotype_ref_table_path('new/path')
    success_2 = utils.get_biotype_ref_path('predefined') == 'new/path'

    utils.set_temp_copy_of_settings_file_as_default()
    utils.remove_temp_copy_of_settings_file()
    assert success_1
    assert success_2


def test_reset_settings_file():
    utils.make_temp_copy_of_settings_file()
    try:
        utils.get_settings_file_path().unlink()
    except FileNotFoundError:
        pass
    utils.update_settings_file('path/to/biotype/ref/file', __biotype_file_key__)
    utils.update_settings_file('path/to/attr/ref/file', __attr_file_key__)
    reset_settings_file()

    success = not utils.get_settings_file_path().exists()

    utils.set_temp_copy_of_settings_file_as_default()
    utils.remove_temp_copy_of_settings_file()
    assert success


def test_start_parallel_session():
    assert False
