import os
from pathlib import Path
from typing import Union

import yaml

from rnalysis import __attr_file_key__, __biotype_file_key__, __path__


def get_settings_file_path():
    """
    Generates the full path of the 'settings.yaml' file.
    :returns: the path of the settings.yaml file.
    :rtype: pathlib.Path
    """
    # return Path(os.path.join(os.path.dirname(__file__), 'settings.yaml'))
    return Path(os.path.join(__path__[0], 'settings.yaml'))


def load_settings_file():
    """
    loads and parses the settings.yaml file into a dictionary.
    :rtype: dict
    """
    settings_pth = get_settings_file_path()
    if not settings_pth.exists():
        return dict()
    with settings_pth.open() as f:
        settings = yaml.safe_load(f)
        if settings is None:
            settings = dict()
        return settings


def update_settings_file(value: str, key: str):
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
    with settings_pth.open('w') as f:
        yaml.safe_dump(out, f)


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
        update_settings_file(input(f'Please insert the full path of {key}:\n'), key)
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
        print(f'Attribute Reference Table used: {ref}')
        return ref


def make_temp_copy_of_settings_file():
    """
    Make a temporary copy ('temp_settings.yaml') of the default settings file ('settings.yaml').
    """
    pth = get_settings_file_path()
    try:
        remove_temp_copy_of_settings_file()
    except FileNotFoundError:
        print("no previous temporary settings file existed")
    if not pth.exists():
        print("no previous settings file exists")
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
        print("no temporary settings file exists")
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
