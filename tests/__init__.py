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


def are_dir_trees_equal(dir1, dir2, compare_contents: bool = True, ignore: list = None):
    """
    Compare two directories recursively. Files in each directory are \
    assumed to be equal if their names and contents are equal.\
    credit: bhttps://stackoverflow.com/a/6681395

    :param dir1: First directory path
    :param dir2: Second directory path

    :return: True if the dir trees are the same and there were no errors while accessing the directories or files, \
    False otherwise.
   """
    ignore = [] if ignore is None else ignore
    ignore = filecmp.DEFAULT_IGNORES + ignore
    dirs_cmp = filecmp.dircmp(dir1, dir2, ignore=ignore)
    if len(dirs_cmp.left_only) > 0 or len(dirs_cmp.right_only) > 0 or len(dirs_cmp.funny_files) > 0:
        print(f"mismatch between {dir1} and {dir2} with left_only={dirs_cmp.left_only}, "
              f"right_only={dirs_cmp.right_only}, funny={dirs_cmp.funny_files}")
        return False
    files_to_cmp = [item for item in dirs_cmp.common_files if Path(item).suffix not in ignore]
    (_, mismatch, errors) = filecmp.cmpfiles(dir1, dir2, files_to_cmp, shallow=False)
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
        if not are_dir_trees_equal(new_dir1, new_dir2, compare_contents, ignore):
            return False
    return True


def is_pantherdb_available():
    # Probe the exact capability the live TestPantherOrthologMapper suite needs: a real
    # `ortholog/matchortho` query for a well-conserved gene (C. elegans G5EDF7 -> human) that a
    # healthy PantherDB answers with a non-empty ortholog mapping. An earlier probe hit the param-
    # free `supportedgenomes` endpoint, but PantherDB degrades unevenly -- it can answer that with
    # HTTP 200 while the ortholog endpoint returns empty/invalid 200 bodies, so the suite ran and
    # its `assert len(...) > 0` checks *failed* instead of *skipping*. Mirror the OrthoInspector
    # probe: exercise the same endpoint the tests do and treat an empty/malformed mapping as
    # "unavailable". The probe hits PantherDB raw (not via the mapper), so a genuine mapper
    # regression still surfaces as a real test failure rather than a skip. Catch RequestException
    # (not just TimeoutError/HTTPError) so an unreachable service can't crash collection at import.
    url = 'https://www.pantherdb.org/services/oai/pantherdb/ortholog/matchortho'
    params = dict(geneInputList='G5EDF7', organism='6239', targetOrganism='9606', orthologType='all')
    try:
        req = requests.post(url, params=params, headers={'accept': 'application/json'}, timeout=(10, 30))
    except requests.exceptions.RequestException:
        return False
    if str(req.status_code)[0] in ['4', '5']:
        return False
    try:
        mapped = req.json()['search']['mapping']['mapped']
    except (requests.exceptions.JSONDecodeError, KeyError, TypeError):
        return False
    return bool(mapped)


def is_uniprot_available():
    try:
        req = requests.get('https://rest.uniprot.org/uniprotkb/search?size=1&query=P53&fields=accession%2Cgene_names')
    except TimeoutError:
        return False
    if str(req.status_code)[0] == '5':
        return False
    return True


def is_ensembl_available():
    try:
        req = requests.get('https://rest.ensembl.org/lookup/id')
    except TimeoutError:
        return False
    if str(req.status_code)[0] == '5':
        return False
    return True


def is_orthoinspector_available():
    # OrthoInspector's per-database `orthologs` endpoints are the flaky part -- the `databases`/`species`
    # endpoints can respond while these time out (the real-world outage this guard exists for). Probe the
    # exact capability the live tests need, mirroring OrthoInspectorOrthologMapper's URL and timeout.
    url = 'https://api.bigest-icube.fr/orthoinspector/Eukaryota2016/species/6239/orthologs/6238'
    try:
        req = requests.get(url, timeout=(10, 30))
    except requests.exceptions.RequestException:
        return False
    if str(req.status_code)[0] in ['4', '5']:
        return False
    return True


def is_phylomedb_available():
    try:
        host = socket.gethostbyname('ftp.phylomedb.org')
        s = socket.create_connection((host, 21), timeout=5)
        s.close()
        return True
    except (socket.gaierror, socket.timeout, ConnectionError, TimeoutError):
        return False


if os.getcwd().endswith('tests'):
    try:
        os.chdir('../../RNAlysis')
    except FileNotFoundError:
        pass
