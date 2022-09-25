import pytest
import pathlib
from rnalysis.utils.settings import *
from rnalysis import __attr_file_key__, __biotype_file_key__


@pytest.fixture
def use_temp_settings_file(request):
    make_temp_copy_of_settings_file()
    request.addfinalizer(remove_temp_copy_of_settings_file)
    request.addfinalizer(set_temp_copy_of_settings_file_as_default)


def test_get_settings_file_path(use_temp_settings_file):
    true_path = pathlib.Path(os.path.join(os.getcwd(), 'rnalysis', 'settings.yaml'))
    assert true_path == get_settings_file_path()


def test_update_attr_ref_path(use_temp_settings_file):
    update_settings_file('path/to/attr/ref/file', __attr_file_key__)
    with get_settings_file_path().open() as f:
        success = False
        for line in f.readlines():
            if line.startswith(__attr_file_key__):
                success = line.startswith('attribute_reference_table: path/to/attr/ref/file')
                if not success:
                    print(f'failiure at: {line}')
    assert success


def test_update_biotype_ref_path(use_temp_settings_file):
    update_settings_file('path/to/biotype/ref/file', __biotype_file_key__)
    with get_settings_file_path().open() as f:
        success = False
        for line in f.readlines():
            if line.startswith(__biotype_file_key__):
                success = line.startswith('biotype_reference_table: path/to/biotype/ref/file')
                if not success:
                    print(f'failiure at: {line}')
    assert success


def test_update_ref_table_attributes(use_temp_settings_file):
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

    assert success
    assert counter == 1
    assert counter_2 == 1
    assert success_2


def test_read_get_attr_ref_path_no_settings(use_temp_settings_file, monkeypatch):
    true_path = 'new/path/file.csv'
    monkeypatch.setattr('builtins.input', lambda __prompt: true_path)
    try:
        get_settings_file_path().unlink()
    except FileNotFoundError:
        pass
    assert get_attr_ref_path('predefined') == true_path


def test_read_get_biotype_ref_path_no_settings(use_temp_settings_file, monkeypatch):
    true_path = 'new/path/file.csv'
    monkeypatch.setattr('builtins.input', lambda __prompt: true_path)
    try:
        get_settings_file_path().unlink()
    except FileNotFoundError:
        pass
    assert get_biotype_ref_path('predefined') == true_path


def test_get_biotype_ref_path(use_temp_settings_file):
    pth = 'path/to/biotype/ref/file'
    update_settings_file(pth, __biotype_file_key__)
    assert get_biotype_ref_path('predefined') == pth


def test_get_attr_ref_path(use_temp_settings_file):
    pth = 'path/to/attr/ref/file'
    update_settings_file(pth, __attr_file_key__)
    assert get_attr_ref_path('predefined') == pth


def test_get_biotype_ref_path_not_predefined():
    pth = 'some/biotype/path'
    assert get_biotype_ref_path(pth) == pth


def test_get_attr_ref_path_not_predefined():
    pth = 'some/attr/ref/path'
    assert get_biotype_ref_path(pth) == pth


def test_reset_settings(use_temp_settings_file):
    with open(get_settings_file_path(), 'w') as f:
        f.write('abcde\nfgh')
    reset_settings()
    assert not get_settings_file_path().exists()


def test_reset_settings_no_previous_file(use_temp_settings_file):
    reset_settings()
    assert not get_settings_file_path().exists()
