import pytest
from rnalysis.utils.ref_tables import *
from rnalysis import __attr_file_key__, __biotype_file_key__


def test_get_settings_file_path():
    make_temp_copy_of_settings_file()

    set_temp_copy_of_settings_file_as_default()
    remove_temp_copy_of_settings_file()
    assert False


def test_get_attr_ref_path():
    make_temp_copy_of_settings_file()
    update_settings_file('path/to/attr/ref/file', __attr_file_key__)
    with get_settings_file_path().open() as f:
        success = False
        for line in f.readlines():
            if line.startswith(__attr_file_key__):
                success = line.startswith('attribute_reference_table: path/to/attr/ref/file')
                if not success:
                    print(f'failiure at: {line}')

    set_temp_copy_of_settings_file_as_default()
    remove_temp_copy_of_settings_file()
    assert success


def test_get_biotype_ref_path():
    make_temp_copy_of_settings_file()
    update_settings_file('path/to/biotype/ref/file', __biotype_file_key__)
    with get_settings_file_path().open() as f:
        success = False
        for line in f.readlines():
            if line.startswith(__biotype_file_key__):
                success = line.startswith('biotype_reference_table: path/to/biotype/ref/file')
                if not success:
                    print(f'failiure at: {line}')

    set_temp_copy_of_settings_file_as_default()
    remove_temp_copy_of_settings_file()
    assert success


def test_update_ref_table_attributes():
    make_temp_copy_of_settings_file()
    try:
        get_settings_file_path().unlink()
    except FileNotFoundError:
        pass
    update_settings_file('path/to/biotype/ref/file', __biotype_file_key__)
    with get_settings_file_path().open() as f:
        counter = 0
        success = False
        for line in f.readlines():
            if line.startswith(__biotype_file_key__):
                counter += 1
                if counter > 1:
                    print(f'Failure with counter={counter}')
                success = line.startswith('biotype_reference_table: path/to/biotype/ref/file')
                if not success:
                    print(f'Failure at: {line}')

    update_settings_file('a/new/path/to/biotype', __biotype_file_key__)
    with get_settings_file_path().open() as f:
        counter_2 = 0
        success_2 = False
        for line in f.readlines():
            if line.startswith(__biotype_file_key__):
                counter_2 += 1
                if counter_2 > 1:
                    print(f'Failure with counter={counter_2}')
                success_2 = line.startswith('biotype_reference_table: a/new/path/to/biotype')
                if not success_2:
                    print(f'Failure at: {line}')

    set_temp_copy_of_settings_file_as_default()
    remove_temp_copy_of_settings_file()
    assert success
    assert counter == 1
    assert counter_2 == 1


def test_biotype_table_assertions():
    assert False


def test_attr_table_assertions():
    assert False
