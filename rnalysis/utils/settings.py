import os
from pathlib import Path
from typing import Union, Any

import yaml

from rnalysis import __attr_file_key__, __biotype_file_key__, __font_key__, __font_size_key__, __stylesheet_key__, \
    __show_tutorial_key__, __report_gen_key__, __databases_key__
from rnalysis.utils import io

DEFAULT_VALUES = {__font_key__: 'Times New Roman', __font_size_key__: 10, __stylesheet_key__: 'base',
                  __show_tutorial_key__: True, __databases_key__: ['Google', 'NCBI Genes', 'UniProtKB'],
                  __report_gen_key__: None}


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
    with settings_pth.open() as f:
        settings = yaml.safe_load(f)
        if settings is None:
            settings = dict()
        return settings


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
    with settings_pth.open('w') as f:
        yaml.safe_dump(out, f)


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


def set_gui_settings(font: str, font_size: int, stylesheet: str, databases, show_tutorial: bool,
                     prompt_report_gen: Union[bool, None]):
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
