from unittest.mock import MagicMock

import pytest

from rnalysis.utils import installs
from rnalysis.utils.installs import *

# Constants for mocking
MOCK_JDK_DIRS = [f'jdk-{installs.JDK_VERSION}.0.1', 'jdk-16', 'jdk-15']
MOCK_JDK_VERSION_OUTPUT = f"java version {installs.JDK_VERSION}"


@pytest.fixture
def mock_jdk_install(monkeypatch):
    monkeypatch.setattr('pathlib.Path.exists', lambda self: True)
    monkeypatch.setattr('subprocess.check_output', lambda *args, **kwargs: MOCK_JDK_VERSION_OUTPUT)
    monkeypatch.setattr(installs, 'get_jdk_path', lambda: JDK_ROOT.joinpath(f'{MOCK_JDK_DIRS[0]}/bin'))


# Test get_jdk_path
def test_get_jdk_path_success(monkeypatch):
    monkeypatch.setattr('os.listdir', lambda x: MOCK_JDK_DIRS)
    monkeypatch.setattr('os.path.isdir', lambda x: True)
    monkeypatch.setattr('pathlib.Path.rglob', lambda self, pattern: [JDK_ROOT.joinpath(f'{MOCK_JDK_DIRS[0]}/bin/java')])
    expected_path = JDK_ROOT.joinpath(f'{MOCK_JDK_DIRS[0]}/bin')
    assert get_jdk_path() == expected_path


def test_get_jdk_path_failure(monkeypatch):
    monkeypatch.setattr('os.listdir', lambda x: [])
    assert get_jdk_path() == ""


# Test is_jdk_installed
def test_is_jdk_installed_success(monkeypatch, mock_jdk_install):
    assert is_jdk_installed()


def test_is_jdk_installed_failure_no_dir(monkeypatch):
    monkeypatch.setattr('pathlib.Path.exists', lambda self: False)
    assert not is_jdk_installed()


def test_is_jdk_installed_failure_wrong_version(monkeypatch):
    monkeypatch.setattr('pathlib.Path.exists', lambda self: True)
    monkeypatch.setattr('subprocess.check_output',
                        lambda *args, **kwargs: MOCK_JDK_VERSION_OUTPUT.replace(str(installs.JDK_VERSION), '15'))
    assert not is_jdk_installed()


def test_is_jdk_installed_subprocess_error(monkeypatch):
    def mock_check_output(*args, **kwargs):
        raise subprocess.CalledProcessError(1, 'cmd')

    monkeypatch.setattr('pathlib.Path.exists', lambda self: True)
    monkeypatch.setattr('subprocess.check_output', mock_check_output)
    assert not is_jdk_installed()

    # Test install_jdk


def test_install_jdk_not_installed(monkeypatch):
    monkeypatch.setattr(installs, 'is_jdk_installed', lambda: False)
    monkeypatch.setattr('pathlib.Path.exists', lambda self: False)
    mock_install = MagicMock()
    monkeypatch.setattr('jdk.install', mock_install)
    install_jdk()
    mock_install.assert_called_once()


def test_install_jdk_already_installed(monkeypatch):
    monkeypatch.setattr(installs, 'is_jdk_installed', lambda: True)
    mock_install = MagicMock()
    monkeypatch.setattr('jdk.install', mock_install)
    install_jdk()
    mock_install.assert_not_called()


# Test is_picard_installed
def test_is_picard_installed_success(monkeypatch, mock_jdk_install):
    monkeypatch.setattr(io, 'run_subprocess', lambda *args, **kwargs: (None, ['Version']))
    assert is_picard_installed()


def test_is_picard_installed_failure(monkeypatch):
    monkeypatch.setattr(io, 'run_subprocess', lambda *args, **kwargs: (None, []))
    assert not is_picard_installed()


def test_is_picard_installed_file_not_found_error(monkeypatch):
    def mock_run_subprocess(*args, **kwargs):
        raise FileNotFoundError()

    monkeypatch.setattr(io, 'run_subprocess', mock_run_subprocess)
    assert not is_picard_installed()


# Test install_picard
def test_install_picard_not_installed(monkeypatch):
    monkeypatch.setattr(installs, 'is_picard_installed', lambda: False)
    monkeypatch.setattr(installs, 'install_jdk', MagicMock())
    monkeypatch.setattr('urllib.request.urlopen', MagicMock())
    monkeypatch.setattr('builtins.open', MagicMock())
    install_picard()


def test_install_picard_found_installation(monkeypatch, mock_jdk_install):
    monkeypatch.setattr(installs, 'is_picard_installed', lambda: True)
    mock_open = MagicMock()
    monkeypatch.setattr(installs, 'urlopen', mock_open)
    monkeypatch.setattr(installs, 'copyfileobj', mock_open)
    install_jdk()
    mock_open.assert_not_called()


def test_install_limma():
    install_limma()


def test_install_deseq2():
    install_deseq2()


def test_install_rsubread():
    install_rsubread()
