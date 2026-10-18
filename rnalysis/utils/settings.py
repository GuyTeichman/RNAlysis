import os
import shutil
import time
import uuid
import warnings
from pathlib import Path
from typing import Any, Union

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

DEFAULT_VALUES = {
    __font_key__: 'Times New Roman',
    __font_size_key__: 10,
    __stylesheet_key__: 'base',
    __show_tutorial_key__: True,
    __databases_key__: ['Google', 'NCBI Genes', 'UniProtKB'],
    __report_gen_key__: None,
}


def get_settings_file_path():
    """
    Generates the full path of the 'settings.yaml' file.
    :returns: the path of the settings.yaml file.
    :rtype: pathlib.Path
    """
    directory = io.get_data_dir()
    if not directory.exists():
        directory.mkdir(parents=True)
    return directory.joinpath('settings.yaml')


def reset_settings():
    pth = get_settings_file_path()
    if pth.exists():
        pth.unlink()


def load_settings_file():
    """
    loads and parses the settings.yaml file into a dictionary.
    :rtype: dict
    """
    settings_pth = get_settings_file_path()
    if not settings_pth.exists():
        return dict()
    try:
        with settings_pth.open() as f:
            settings = yaml.safe_load(f)
    except (yaml.YAMLError, UnicodeDecodeError):
        return _handle_corrupt_settings_file(settings_pth)
    if settings is None:  # an empty file (or one with only comments) is a valid, empty settings file
        return dict()
    if not isinstance(settings, dict):
        return _handle_corrupt_settings_file(settings_pth)
    return settings


_REPLACE_RETRIES = 5
_REPLACE_RETRY_DELAY = 0.05  # seconds; doubled after every failed attempt


class CorruptSettingsWarning(UserWarning):
    """
    Warns that the settings file could not be read and was reset to the default settings. \
    The GUI catches this warning while building its window (before its log console exists) and shows it in a message box.
    """


def _handle_corrupt_settings_file(settings_pth: Path) -> dict:
    """
    Move an unreadable settings file aside to 'settings.yaml.bad' (overwriting any previous backup), \
    warn the user, and return empty settings so that the default values are used instead.
    A corrupt preferences file should never prevent RNAlysis from starting.
    """
    backup_pth = settings_pth.with_name(settings_pth.name + '.bad')
    try:
        os.replace(settings_pth, backup_pth)
        backup_msg = f'A copy of the unreadable file was saved to: {backup_pth}'
    except OSError:
        try:  # e.g. the file is locked by another program on Windows: keep a copy instead of moving it
            shutil.copyfile(settings_pth, backup_pth)
            backup_msg = f'A copy of the unreadable file was saved to: {backup_pth}'
        except OSError:
            backup_msg = f'RNAlysis could not save a copy of the unreadable file ({settings_pth}).'
    message = (
        'Your RNAlysis settings file could not be read (it may be damaged or incomplete), '
        'so RNAlysis is using its default settings instead. '
        'Any settings you changed earlier (such as font or Reference Table paths) need to be set again. '
        f'{backup_msg}'
    )
    warnings.warn(message, CorruptSettingsWarning)
    return dict()


def _replace_with_retries(src: Path, dst: Path):
    """
    os.replace() src over dst. On Windows, replacing a file fails with PermissionError while another program \
    (e.g. an antivirus scanner, or another RNAlysis process reading the settings) briefly holds it open, \
    so retry a few times with a short backoff before giving up.
    """
    delay = _REPLACE_RETRY_DELAY
    for attempt in range(_REPLACE_RETRIES):
        try:
            os.replace(src, dst)
            return
        except PermissionError:
            if attempt == _REPLACE_RETRIES - 1:
                raise
            time.sleep(delay)
            delay *= 2


def update_settings_file(value: Any, key: str):
    """
    Receives a key and a value, and updates/adds the key and value to the settings.yaml file.

    :param value: the value to be added/updated (such as Reference Table path)
    :type value: str
    :param key: the key to be added/updated (such as __attr_file_key__)
    :type key: str
    """
    settings_pth = get_settings_file_path()
    out = load_settings_file()
    out[key] = value
    # write to a temporary file in the same directory and then atomically replace the settings file,
    # so that a crash or kill mid-write can never leave a truncated (corrupt) settings file behind.
    # the temporary file's name is unique, so that concurrent writers never write into each other's temporary file.
    # (the file is opened with plain text-mode 'x', so its encoding, newlines and permissions match a plain open('w'))
    tmp_pth = settings_pth.with_name(f'{settings_pth.name}.{os.getpid()}.{uuid.uuid4().hex[:8]}.tmp')
    try:
        with tmp_pth.open('x') as f:
            yaml.safe_dump(out, f)
            f.flush()
            os.fsync(f.fileno())
        try:
            _replace_with_retries(tmp_pth, settings_pth)
        except PermissionError:
            # the settings file stays locked by another program; fall back to writing it in place
            # (the pre-4.3.1 behavior), which Windows allows even while another program has the file open
            with settings_pth.open('w') as f:
                yaml.safe_dump(out, f)
    finally:
        # after a successful os.replace() the temporary file no longer exists; otherwise (a failed dump,
        # or the in-place fallback above) remove it so that no stray temporary files are left behind
        try:
            tmp_pth.unlink(missing_ok=True)
        except OSError:
            pass


def is_setting_in_file(key) -> bool:
    """
    Returns True if a settings value is defined for given key, and False otherwise.
    :type key: str
    :param key: the key in the settings file whose value to read.
    :rtype: bool
    """
    settings = load_settings_file()
    return key in settings


def read_value_from_settings(key):
    """
    Attempt to read the value corresponding to a given key from the settings.yaml file. \
    If the key was not previously defined, the user will be prompted to define it.

    :type key: str
    :param key: the key in the settings file whose value to read.

    :return:
    The path of the reference table.
    """
    settings = load_settings_file()
    if key not in settings:
        if key in DEFAULT_VALUES:
            val = DEFAULT_VALUES[key]
        else:
            val = input(f'Please insert the value of {key}:\n')
        update_settings_file(val, key)
        settings = load_settings_file()
    return settings[key]


def get_biotype_ref_path(ref: Union[str, Path]):
    """
    Returns the predefined Biotype Reference Table path from the settings file if ref='predefined', \
    otherwise returns 'ref' unchanged.

    :param ref: the 'ref' argument from a filtering module/enrichment module function
    :type ref: str, pathlib.Path or 'predefined'
    :returns: if ref is 'predefined', returns the predefined Biotype Reference Table path from the same file. \
    Otherwise, returns 'ref'.
    :rtype: str
    """
    if ref == 'predefined':
        pth = read_value_from_settings(__biotype_file_key__)
        print(f'Biotype Reference Table used: {pth}')

        return pth
    else:
        if ref is not None and ref != '':
            print(f'Biotype Reference Table used: {ref}')
        return ref


def get_attr_ref_path(ref):
    """
    Returns the predefined Attribute Reference Table path from the settings file if ref='predefined', \
    otherwise returns 'ref' unchanged.

    :param ref: the 'ref' argument from a filtering module/enrichment module function
    :type ref: str, pathlib.Path or 'predefined'
    :returns: if ref is 'predefined', returns the predefined Attribute Reference Table path from the same file. \
    Otherwise, returns 'ref'.
    :rtype: str
    """
    if ref == 'predefined':
        pth = read_value_from_settings(__attr_file_key__)
        print(f'Attribute Reference Table used: {pth}')
        return pth
    else:
        if ref is not None and ref != '':
            print(f'Attribute Reference Table used: {ref}')
        return ref


def get_gui_settings():
    font = read_value_from_settings(__font_key__)
    font_size = int(read_value_from_settings(__font_size_key__))
    stylesheet = read_value_from_settings(__stylesheet_key__)
    databases = read_value_from_settings(__databases_key__)
    show_tutorial = read_value_from_settings(__show_tutorial_key__)
    prompt_report_gen = read_value_from_settings(__report_gen_key__)
    return font, font_size, stylesheet, databases, show_tutorial, prompt_report_gen


def set_gui_settings(
    font: str, font_size: int, stylesheet: str, databases, show_tutorial: bool, prompt_report_gen: Union[bool, None]
):
    update_settings_file(font, __font_key__)
    update_settings_file(str(font_size), __font_size_key__)
    update_settings_file(stylesheet, __stylesheet_key__)
    update_settings_file(databases, __databases_key__)
    update_settings_file(show_tutorial, __show_tutorial_key__)
    update_settings_file(prompt_report_gen, __report_gen_key__)


def get_databases_settings():
    return read_value_from_settings(__databases_key__)


def get_show_tutorial_settings():
    return read_value_from_settings(__show_tutorial_key__)


def set_show_tutorial_settings(show_tutorial: bool):
    update_settings_file(show_tutorial, __show_tutorial_key__)


def get_report_gen_settings():
    return read_value_from_settings(__report_gen_key__)


def set_report_gen_settings(prompt_report_gen: Union[bool, None]):
    update_settings_file(prompt_report_gen, __report_gen_key__)


def set_table_settings(attr_ref_path: str, biotype_ref_path: str):
    if attr_ref_path:
        update_settings_file(attr_ref_path, __attr_file_key__)
    if biotype_ref_path:
        update_settings_file(biotype_ref_path, __biotype_file_key__)


def make_temp_copy_of_settings_file():
    """
    Make a temporary copy ('temp_settings.yaml') of the default settings file ('settings.yaml').
    """
    pth = get_settings_file_path()
    try:
        remove_temp_copy_of_settings_file()
    except FileNotFoundError:
        print('no previous temporary settings file existed')
    if not pth.exists():
        print('no previous settings file exists')
        return
    with open(os.path.join(str(pth.parent), 'temp_settings.yaml'), 'w') as tempfile, pth.open() as originfile:
        tempfile.writelines(originfile.readlines())


def set_temp_copy_of_settings_file_as_default():
    """
    Copy the contents of the temporary settings file ('temp_settings.yaml') into the default settings file \
    ('settings.yaml'), if a temporary settings file exists.
    """
    pth = get_settings_file_path()
    if pth.exists():
        pth.unlink()
    if not Path(os.path.join(str(pth.parent), 'temp_settings.yaml')).exists():
        print('no temporary settings file exists')
        return
    with open(os.path.join(str(pth.parent), 'temp_settings.yaml')) as temp_file, pth.open('w') as original_file:
        original_file.writelines(temp_file.readlines())


def remove_temp_copy_of_settings_file():
    """
    Remove the temporary copy of the settings file ('temp_settings.yaml'), if such file exists.
    """
    pth = get_settings_file_path()
    tmp_pth = Path(os.path.join(str(pth.parent), 'temp_settings.yaml'))
    if tmp_pth.exists():
        tmp_pth.unlink()
