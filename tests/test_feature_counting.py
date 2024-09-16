from rnalysis.utils.feature_counting import *



def test_create_featurecounts_script():
    expected_path = 'tests/test_files/featurecounts_tests/expected_featurecounts_script.R'
    output_dir = Path('tests/test_files/featurecounts_tests/outdir')
    kwargs = {'arg': True, 'arg2': None, 'arg3': 'string', 'arg4': ['str1', 'str2', 'str 3']}
    with open(expected_path) as f:
        expected = f.read()

    runner = FeatureCountsRunner(kwargs, output_dir)
    out_path = runner._create_script()
    assert Path(out_path).exists()
    with open(out_path) as f:
        out = f.read()

    try:
        assert out == expected
    finally:
        out_path.unlink()
