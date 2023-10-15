import filecmp
import os
import shutil
import socket
from pathlib import Path

import requests

__attr_ref__ = 'tests/test_files/attr_ref_table_for_tests.csv'
__biotype_ref__ = 'tests/test_files/biotype_ref_table_for_tests.csv'


def unlink_tree(dir):
    for item in Path(dir).iterdir():
        if 'gitignore' in item.name:
            continue
        if item.is_file():
            item.unlink()
        else:
            shutil.rmtree(item)


def are_dir_trees_equal(dir1, dir2, compare_contents: bool = True):
    """
    Compare two directories recursively. Files in each directory are \
    assumed to be equal if their names and contents are equal.\
    credit: bhttps://stackoverflow.com/a/6681395

    :param dir1: First directory path
    :param dir2: Second directory path

    :return: True if the dir trees are the same and there were no errors while accessing the directories or files, \
    False otherwise.
   """

    dirs_cmp = filecmp.dircmp(dir1, dir2)
    if len(dirs_cmp.left_only) > 0 or len(dirs_cmp.right_only) > 0 or \
        len(dirs_cmp.funny_files) > 0:
        print(f"mismatch between {dir1} and {dir2} with left_only={dirs_cmp.left_only}, "
              f"right_only={dirs_cmp.right_only}, funny={dirs_cmp.funny_files}")
        return False
    (_, mismatch, errors) = filecmp.cmpfiles(
        dir1, dir2, dirs_cmp.common_files, shallow=False)
    if (len(mismatch) > 0 or len(errors) > 0) and compare_contents:
        print(f"mismatch between {dir1} and {dir2} in the files {mismatch} with errors {errors}")
        for item in mismatch:
            items = []
            for this_dir in [dir1, dir2]:
                pth = Path(this_dir).joinpath(item)
                try:
                    with open(pth) as f:
                        txt = f.read()
                except UnicodeDecodeError:
                    with open(pth, 'rb') as f:
                        txt = f.read()
                items.append(txt)
            if items[0] != items[1]:
                for i in items:
                    print(i)
                    print('---------------------------')
                return False
    for common_dir in dirs_cmp.common_dirs:
        new_dir1 = Path(dir1).joinpath(common_dir).as_posix()
        new_dir2 = Path(dir2).joinpath(common_dir).as_posix()
        if not are_dir_trees_equal(new_dir1, new_dir2, compare_contents):
            return False
    return True


def is_uniprot_available():
    req = requests.get('https://rest.uniprot.org/uniprotkb/search?size=1&query=P53&fields=accession%2Cgene_names')
    if str(req.status_code)[0] == '5':
        return False
    return True


def is_ensembl_available():
    req = requests.get('https://rest.ensembl.org/lookup/id')
    if str(req.status_code)[0] == '5':
        return False
    return True


def is_phylomedb_available():
    try:
        host = socket.gethostbyname('ftp.phylomedb.org')
        s = socket.create_connection((host, 21), timeout=5)
        s.close()
        return True
    except (socket.gaierror, ConnectionError):
        return False


if os.getcwd().endswith('tests'):
    try:
        os.chdir('../../RNAlysis')
    except FileNotFoundError:
        pass
