from rnalysis.general import *
from rnalysis import __biotype_file_key__, __attr_file_key__
import pytest

from rnalysis.filtering import Filter


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


def test_print_settings_file(capfd):
    ref_tables.make_temp_copy_of_settings_file()
    set_biotype_ref_table_path('temp_path')
    set_attr_ref_table_path('temp_path')

    print_settings_file()

    ref_tables.set_temp_copy_of_settings_file_as_default()
    ref_tables.remove_temp_copy_of_settings_file()


def test_set_attr_ref_path():
    ref_tables.make_temp_copy_of_settings_file()
    if ref_tables.get_settings_file_path().exists():
        ref_tables.get_settings_file_path().unlink()

    set_attr_ref_table_path('old/path')
    success_1 = ref_tables.get_attr_ref_path('predefined') == 'old/path'

    set_attr_ref_table_path('new/path')
    success_2 = ref_tables.get_attr_ref_path('predefined') == 'new/path'

    ref_tables.set_temp_copy_of_settings_file_as_default()
    ref_tables.remove_temp_copy_of_settings_file()
    assert success_1
    assert success_2


def test_set_biotype_ref_path():
    ref_tables.make_temp_copy_of_settings_file()
    if ref_tables.get_settings_file_path().exists():
        ref_tables.get_settings_file_path().unlink()

    set_biotype_ref_table_path('old/path')
    success_1 = ref_tables.get_biotype_ref_path('predefined') == 'old/path'

    set_biotype_ref_table_path('new/path')
    success_2 = ref_tables.get_biotype_ref_path('predefined') == 'new/path'

    ref_tables.set_temp_copy_of_settings_file_as_default()
    ref_tables.remove_temp_copy_of_settings_file()
    assert success_1
    assert success_2


def test_reset_settings_file():
    ref_tables.make_temp_copy_of_settings_file()
    try:
        ref_tables.get_settings_file_path().unlink()
    except FileNotFoundError:
        pass
    ref_tables.update_settings_file('path/to/biotype/ref/file', __biotype_file_key__)
    ref_tables.update_settings_file('path/to/attr/ref/file', __attr_file_key__)
    reset_settings_file()

    success = not ref_tables.get_settings_file_path().exists()

    ref_tables.set_temp_copy_of_settings_file_as_default()
    ref_tables.remove_temp_copy_of_settings_file()
    assert success


def test_save_to_csv(monkeypatch):
    monkeypatch.setattr('rnalysis.utils.io.save_csv',
                        lambda x, y: None if isinstance(x, pd.DataFrame) else (_ for _ in ()).throw(AssertionError))
    df = pd.DataFrame([[1, 2, 3], [4, 5, 6]])
    save_to_csv(df, 'filename')
    f = Filter(('filename', df))
    save_to_csv(f, 'filename')
    with pytest.raises(TypeError):
        save_to_csv(5, 'filename')
