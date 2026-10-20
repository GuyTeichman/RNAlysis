"""Guards for ``.gitattributes``: LF text everywhere, and no line-ending conversion of byte-exact content.

Without these attributes a clone with ``core.autocrlf=true`` (Git for Windows' system default) gets CRLF working
copies, which the ``mixed-line-ending --fix=lf`` pre-commit hook then keeps "fixing", and tests may read fixtures
with different bytes than CI does.
"""

import shutil
import subprocess
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parent.parent

pytestmark = pytest.mark.skipif(
    shutil.which('git') is None or not (REPO_ROOT / '.git').exists(), reason='requires a git checkout'
)


def _attrs(path: str) -> dict:
    out = subprocess.run(
        ['git', 'check-attr', 'text', 'eol', 'binary', '--', path],
        cwd=REPO_ROOT,
        capture_output=True,
        text=True,
        check=True,
    ).stdout
    return {line.split(': ')[1]: line.split(': ')[2] for line in out.splitlines()}


@pytest.mark.parametrize(
    'path', ['rnalysis/filtering.py', 'README.rst', '.pre-commit-config.yaml', 'docs/source/a.rst', 'make.bat']
)
def test_text_files_are_lf(path):
    # make.bat included: a CRLF rule would fight the `mixed-line-ending --fix=lf` pre-commit hook
    attrs = _attrs(path)
    assert attrs['text'] == 'auto'
    assert attrs['eol'] == 'lf'


@pytest.mark.parametrize(
    'path',
    [
        'tests/test_files/counted.csv',
        'tests/test_files/deseq2_tests/case1/expected_deseq_script_1.R',
        'rnalysis/data_files/r_templates/logging.R',
        'rnalysis/data_files/report_templates/vis-network.min.js',
    ],
)
def test_byte_exact_content_is_not_converted(path):
    assert _attrs(path)['text'] == 'unset'


@pytest.mark.parametrize('path', ['docs/source/logo.png', 'tests/test_files/test_session.rnal', 'x/y.fastq.gz'])
def test_binaries_are_binary(path):
    assert _attrs(path)['binary'] == 'set'
