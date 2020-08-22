import os

try:
    os.chdir('tests/')
except FileNotFoundError:
    pass

__attr_ref__ = 'test_files/attr_ref_table_for_tests.csv'
__biotype_ref__ = 'test_files/biotype_ref_table_for_tests.csv'
