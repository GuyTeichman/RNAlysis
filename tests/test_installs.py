from unittest.mock import patch, MagicMock

import pytest

from rnalysis.utils.installs import *

# Constants for mocking
MOCK_JDK_DIRS = ['jdk-17.0.1', 'jdk-16', 'jdk-15']
MOCK_JDK_VERSION_OUTPUT = "java version 17"


# Test get_jdk_path
def test_get_jdk_path_success():
    with patch('os.listdir', return_value=MOCK_JDK_DIRS), \
        patch('os.path.isdir', return_value=True):
        expected_path = JDK_ROOT.joinpath(f'{MOCK_JDK_DIRS[0]}/bin')
        assert get_jdk_path() == expected_path


def test_get_jdk_path_failure():
    with patch('os.listdir', return_value=[]):
        with pytest.raises(FileNotFoundError):
            get_jdk_path()


# Test is_jdk_installed
def test_is_jdk_installed_success():
    with patch('pathlib.Path.exists', return_value=True), \
        patch('subprocess.check_output', return_value=MOCK_JDK_VERSION_OUTPUT):
        assert is_jdk_installed()


def test_is_jdk_installed_failure():
    with patch('pathlib.Path.exists', return_value=False):
        assert not is_jdk_installed()


def test_is_jdk_installed_subprocess_error():
    with patch('pathlib.Path.exists', return_value=True), \
        patch('subprocess.check_output', side_effect=subprocess.CalledProcessError(1, 'cmd')):
        assert not is_jdk_installed()


# Test install_jdk
def test_install_jdk_not_installed():
    with patch('rnalysis.utils.installs.is_jdk_installed', return_value=False), \
        patch('pathlib.Path.exists', return_value=False), \
        patch('rnalysis.utils.installs.jdk.install') as mock_install:
        install_jdk()
        mock_install.assert_called_once()


def test_install_jdk_already_installed():
    with patch('rnalysis.utils.installs.is_jdk_installed', return_value=True), \
        patch('rnalysis.utils.installs.jdk.install') as mock_install:
        install_jdk()
        mock_install.assert_not_called()


# Test is_picard_installed
def test_is_picard_installed_success():
    with patch('rnalysis.utils.installs.io.run_subprocess', return_value=(None, ['Version'])):
        assert is_picard_installed()


def test_is_picard_installed_failure():
    with patch('rnalysis.utils.installs.io.run_subprocess', return_value=(None, [])):
        assert not is_picard_installed()


def test_is_picard_installed_file_not_found_error():
    with patch('rnalysis.utils.installs.io.run_subprocess', side_effect=FileNotFoundError):
        assert not is_picard_installed()

    # Test install_picard


def test_install_picard_not_installed():
    with patch('rnalysis.utils.installs.is_picard_installed', return_value=False), \
        patch('rnalysis.utils.installs.install_jdk'), \
        patch('urllib.request.urlopen'), \
        patch('builtins.open', new_callable=MagicMock):
        install_picard()

#
# def test_install_limma():
#     install_limma()
#
#
# def test_install_deseq2():
#     install_deseq2()
