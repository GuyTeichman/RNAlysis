from pathlib import Path

DOCS_DIR = Path('docs/source')
USER_GUIDE = DOCS_DIR / 'user_guide.rst'
GUI_GUIDE = DOCS_DIR / 'user_guide_gui.rst'


def test_user_guide_examples_use_current_api_names_and_paths():
    # The guides contain non-ASCII characters (e.g. Polars table-repr box-drawing glyphs), so the
    # encoding must be pinned to UTF-8 -- otherwise read_text() uses the platform default (cp1252 on
    # Windows), which fails to decode them and breaks this test on Windows only.
    user_guide = USER_GUIDE.read_text(encoding='utf-8')
    gui_guide = GUI_GUIDE.read_text(encoding='utf-8')

    assert 'normalize_quantile' not in user_guide
    assert 'normalize_quantile' not in gui_guide
    assert 'normalize_to_quantile' in user_guide
    assert 'normalize_to_quantile' in gui_guide

    assert 'normalize_mrn' not in user_guide
    assert 'normalize_mrn' not in gui_guide
    assert 'normalize_median_of_ratios' in user_guide
    assert 'normalize_median_of_ratios' in gui_guide

    assert 'test_files/test_files/' not in user_guide
    assert '>>> counts = CountFilter(' not in user_guide
    assert '>>> pipe = Pipeline()' not in user_guide
    assert '>>> counts = filtering.CountFilter(' in user_guide
    assert '>>> pipe = filtering.Pipeline()' in user_guide

    assert 'FeatureSet.enrich_randomization' not in user_guide
    assert 'enrich_randomization' not in user_guide
    assert 'FeatureSet.user_defined_enrichment' in user_guide

    assert 'read_attr_ref_table_path' not in user_guide
    assert 'read_biotype_ref_table_path' not in user_guide
