import os
import requests

__attr_ref__ = 'tests/test_files/attr_ref_table_for_tests.csv'
__biotype_ref__ = 'tests/test_files/biotype_ref_table_for_tests.csv'

def is_uniprot_available():
    req = requests.get('https://rest.uniprot.org/idmapping/run')
    if str(req.status_code)[0] == '5':
        return False
    return True

def is_ensembl_available():
    req = requests.get('https://rest.ensembl.org/lookup/id')
    if str(req.status_code)[0] == '5':
        return False
    return True

if os.getcwd().endswith('tests'):
    try:
        os.chdir('../../RNAlysis')
    except FileNotFoundError:
        pass
