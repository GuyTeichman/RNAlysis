from rnalysis.utils.feature_counting import *


def test_install_rsubread():
    install_rsubread()


def test_create_featurecounts_script():
    expected_path = 'tests/test_files/expected_featurecounts_script.R'
    output_dir = Path('C:/specified/output/dir')
    kwargs = {'arg': True, 'arg2': None, 'arg3': 'string', 'arg4': ['str1', 'str2', 'str 3']}
    with open(expected_path) as f:
        expected = f.read()

    out_path = create_featurecounts_script(kwargs, output_dir)
    assert Path(out_path).exists()
    with open(out_path) as f:
        out = f.read()

    assert out == expected


# def test_run_featurecounts_analysis():
#     assert False
