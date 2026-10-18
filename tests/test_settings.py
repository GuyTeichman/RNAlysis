import warnings

import pytest
import yaml

from rnalysis import (
    __attr_file_key__,
    __biotype_file_key__,
    __databases_key__,
    __font_key__,
    __font_size_key__,
    __report_gen_key__,
    __show_tutorial_key__,
    __stylesheet_key__,
)
from rnalysis.utils import io
from rnalysis.utils.settings import *


@pytest.fixture
def use_temp_settings_file(request):
    make_temp_copy_of_settings_file()
    request.addfinalizer(remove_temp_copy_of_settings_file)
    request.addfinalizer(set_temp_copy_of_settings_file_as_default)


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


@pytest.fixture
def isolated_settings_dir(tmp_path, monkeypatch):
    """Point the settings file at a temporary directory, so the user's real settings file is never touched."""
    monkeypatch.setattr(io, 'get_data_dir', lambda: tmp_path)
    return tmp_path


@pytest.mark.parametrize(
    'content',
    [
        'font: Arial\nfont_size: [10\n',  # malformed YAML (unterminated flow sequence)
        'font: "Arial\n',  # truncated mid-string
        '\tfont: Arial\n',  # tab indentation is invalid YAML
    ],
)
def test_load_settings_file_malformed_yaml(isolated_settings_dir, content):
    pth = get_settings_file_path()
    pth.write_text(content)
    with pytest.warns(UserWarning, match=r'settings\.yaml\.bad'):
        settings = load_settings_file()
    assert settings == {}
    backup = isolated_settings_dir.joinpath('settings.yaml.bad')
    assert backup.exists()
    assert backup.read_text() == content
    assert not pth.exists()


def test_load_settings_file_undecodable_bytes(isolated_settings_dir):
    pth = get_settings_file_path()
    content = b'font: \xff\xfe\x00\x81garbage\n'
    pth.write_bytes(content)
    with pytest.warns(UserWarning, match=r'settings\.yaml\.bad'):
        settings = load_settings_file()
    assert settings == {}
    assert isolated_settings_dir.joinpath('settings.yaml.bad').read_bytes() == content


@pytest.mark.parametrize('content', ['just a truncated strin', '- a\n- list\n', '42\n'])
def test_load_settings_file_non_dict_content(isolated_settings_dir, content):
    pth = get_settings_file_path()
    pth.write_text(content)
    with pytest.warns(UserWarning, match=r'settings\.yaml\.bad'):
        settings = load_settings_file()
    assert settings == {}
    assert isolated_settings_dir.joinpath('settings.yaml.bad').read_text() == content
    assert not pth.exists()


def test_load_settings_file_corrupt_overwrites_previous_backup(isolated_settings_dir):
    backup = isolated_settings_dir.joinpath('settings.yaml.bad')
    backup.write_text('old backup')
    get_settings_file_path().write_text('font: [\n')
    with pytest.warns(UserWarning):
        load_settings_file()
    assert backup.read_text() == 'font: [\n'


@pytest.mark.parametrize('content', ['', '\n', '# just a comment\n'])
def test_load_settings_file_empty_file(isolated_settings_dir, content):
    get_settings_file_path().write_text(content)
    with warnings.catch_warnings():
        warnings.simplefilter('error')
        assert load_settings_file() == {}
    assert not isolated_settings_dir.joinpath('settings.yaml.bad').exists()
    assert get_settings_file_path().exists()


def test_load_settings_file_no_file(isolated_settings_dir):
    with warnings.catch_warnings():
        warnings.simplefilter('error')
        assert load_settings_file() == {}
    assert not isolated_settings_dir.joinpath('settings.yaml.bad').exists()


def test_get_gui_settings_corrupt_file_falls_back_to_defaults(isolated_settings_dir):
    get_settings_file_path().write_text('font: "Arial\nfont_size: 1')
    with pytest.warns(UserWarning, match=r'settings\.yaml\.bad'):
        font, font_size, stylesheet, databases, show_tutorial, prompt_report_gen = get_gui_settings()
    assert font == DEFAULT_VALUES[__font_key__]
    assert font_size == DEFAULT_VALUES[__font_size_key__]
    assert stylesheet == DEFAULT_VALUES[__stylesheet_key__]
    assert databases == DEFAULT_VALUES[__databases_key__]
    assert show_tutorial == DEFAULT_VALUES[__show_tutorial_key__]
    assert prompt_report_gen == DEFAULT_VALUES[__report_gen_key__]
    # the defaults are written back to a fresh, valid settings file
    assert load_settings_file()[__font_key__] == DEFAULT_VALUES[__font_key__]


def test_update_settings_file_failed_dump_keeps_original(isolated_settings_dir, monkeypatch):
    pth = get_settings_file_path()
    update_settings_file('Arial', __font_key__)
    original_bytes = pth.read_bytes()

    def failing_dump(data, stream=None, **kwargs):
        stream.write('font: Tim')  # simulate a partial write before the failure
        raise RuntimeError('simulated crash mid-write')

    monkeypatch.setattr(yaml, 'safe_dump', failing_dump)
    with pytest.raises(RuntimeError):
        update_settings_file('Comic Sans', __font_key__)
    assert pth.read_bytes() == original_bytes
    assert sorted(p.name for p in isolated_settings_dir.iterdir()) == ['settings.yaml']


def test_corrupt_settings_warning_can_be_caught_by_the_gui(isolated_settings_dir):
    # the GUI records this warning category while building its window, and shows it in a message box
    get_settings_file_path().write_text('font: [\n')
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter('always', CorruptSettingsWarning)
        load_settings_file()
    corrupt = [w for w in caught if issubclass(w.category, CorruptSettingsWarning)]
    assert len(corrupt) == 1
    assert issubclass(CorruptSettingsWarning, UserWarning)
    assert 'settings.yaml.bad' in str(corrupt[0].message)


def test_corrupt_settings_backup_falls_back_to_copy(isolated_settings_dir, monkeypatch):
    pth = get_settings_file_path()
    pth.write_text('font: [\n')
    real_replace = os.replace

    def replace_fails_for_backup(src, dst):
        if str(dst).endswith('.bad'):
            raise PermissionError('file is locked')
        return real_replace(src, dst)

    monkeypatch.setattr(os, 'replace', replace_fails_for_backup)
    with pytest.warns(UserWarning, match=r'settings\.yaml\.bad'):
        assert load_settings_file() == {}
    assert isolated_settings_dir.joinpath('settings.yaml.bad').read_text() == 'font: [\n'


def test_update_settings_file_retries_when_target_locked(isolated_settings_dir, monkeypatch):
    update_settings_file('Arial', __font_key__)
    real_replace = os.replace
    calls = []

    def flaky_replace(src, dst):
        calls.append(src)
        if len(calls) < 3:
            raise PermissionError('[WinError 5] Access is denied')
        return real_replace(src, dst)

    monkeypatch.setattr(os, 'replace', flaky_replace)
    monkeypatch.setattr('rnalysis.utils.settings.time.sleep', lambda s: None)
    update_settings_file('Comic Sans', __font_key__)
    assert len(calls) == 3
    assert load_settings_file()[__font_key__] == 'Comic Sans'
    assert sorted(p.name for p in isolated_settings_dir.iterdir()) == ['settings.yaml']


def test_update_settings_file_falls_back_to_direct_write_when_target_stays_locked(isolated_settings_dir, monkeypatch):
    update_settings_file('Arial', __font_key__)

    def always_locked(src, dst):
        raise PermissionError('[WinError 5] Access is denied')

    monkeypatch.setattr(os, 'replace', always_locked)
    monkeypatch.setattr('rnalysis.utils.settings.time.sleep', lambda s: None)
    update_settings_file('Comic Sans', __font_key__)
    assert load_settings_file()[__font_key__] == 'Comic Sans'
    assert sorted(p.name for p in isolated_settings_dir.iterdir()) == ['settings.yaml']


def test_update_settings_file_uses_unique_temp_files(isolated_settings_dir, monkeypatch):
    real_replace = os.replace
    temp_names = []

    def recording_replace(src, dst):
        temp_names.append(Path(src).name)
        return real_replace(src, dst)

    monkeypatch.setattr(os, 'replace', recording_replace)
    update_settings_file('Arial', __font_key__)
    update_settings_file('Comic Sans', __font_key__)
    assert len(temp_names) == 2
    assert temp_names[0] != temp_names[1]
    assert all(name.startswith('settings.yaml.') and name.endswith('.tmp') for name in temp_names)


def test_update_settings_file_leaves_no_temp_file(isolated_settings_dir):
    update_settings_file('Arial', __font_key__)
    update_settings_file(12, __font_size_key__)
    assert sorted(p.name for p in isolated_settings_dir.iterdir()) == ['settings.yaml']


def test_update_settings_file_format_unchanged(isolated_settings_dir):
    values = {
        __font_key__: 'Times New Roman',
        __font_size_key__: '10',
        __stylesheet_key__: 'base',
        __databases_key__: ['Google', 'NCBI Genes', 'UniProtKB'],
        __show_tutorial_key__: True,
        __report_gen_key__: None,
        __attr_file_key__: 'path/to/attr ref.csv',
    }
    for key, val in values.items():
        update_settings_file(val, key)
    # the on-disk format must be exactly what a plain open('w') + yaml.safe_dump produces
    expected = isolated_settings_dir.joinpath('expected.yaml')
    with expected.open('w') as f:
        yaml.safe_dump(values, f)
    assert get_settings_file_path().read_bytes() == expected.read_bytes()
    assert load_settings_file() == values
