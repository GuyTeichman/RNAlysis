import os

__attr_ref__ = 'tests/test_files/attr_ref_table_for_tests.csv'
__biotype_ref__ = 'tests/test_files/biotype_ref_table_for_tests.csv'

if os.getcwd().endswith('tests'):
    try:
        os.chdir('../../RNAlysis')
    except FileNotFoundError:
        pass
