import os
import platform
import shutil
import subprocess
import warnings
from pathlib import Path
from shutil import copyfileobj
from typing import Literal, Union
from urllib.request import urlopen

from rnalysis.utils import io

try:
    import jdk
except ImportError:  # pragma: no cover
    warnings.warn("'install-jdk' not found. To use PicardTools, you would need to install JDK manually.")


    class jdk:  # pragma: no cover
        @staticmethod
        def install(version):
            warnings.warn("Cannot install JDK.")

PICARD_JAR = Path(os.environ.get('PICARDTOOLS_JAR', io.get_data_dir().joinpath('picard.jar')))
JDK_VERSION = "21"
JDK_ROOT = io.get_data_dir().joinpath(f'jdk{JDK_VERSION}')


def get_jdk_path()->Path:
    try:
        if platform.system() == 'Windows':
            pattern = "java.exe"
        else:
            pattern = "java"
        # List all files and directories in the given directory
        items = os.listdir(JDK_ROOT)
        # Filter for directories starting with 'jdk-'
        jdk_directories = sorted([item for item in items if
                                  item.startswith('jdk-') and os.path.isdir(os.path.join(JDK_ROOT, item))],
                                 reverse=True)
        if len(jdk_directories) == 0:
            raise FileNotFoundError('No JDK directory found')
        base_dir = JDK_ROOT.joinpath(f'{jdk_directories[0]}')

        # Filter for files starting with 'java' or 'java.exe'
        matches = list(base_dir.rglob(pattern))
        if len(matches) == 0:
            raise FileNotFoundError(f'No java executable found in {base_dir}')
        return matches[0].parent
    except FileNotFoundError:
        # find jdk installation in PATH
        return Path()

def is_jdk_installed():
    try:
        if not JDK_ROOT.exists():
            return False
        # Run the "java -version" command and capture the output
        jdk_path = get_jdk_path()
        output = subprocess.check_output(
            [Path(jdk_path).joinpath("java").as_posix(), "-version"],
            stderr=subprocess.STDOUT, text=True)
        # Check if the output contains "version X" or "X." (for Java X)
        if f"version {JDK_VERSION}" in output or f" {JDK_VERSION}." in output:
            return True
    except (subprocess.CalledProcessError, FileNotFoundError):
        # If the "java -version" command returns an error, Java X is not installed.
        return False
    return False


def install_jdk():
    if not is_jdk_installed():
        print("Installing Java...")
        if JDK_ROOT.exists():
            shutil.rmtree(JDK_ROOT)
        JDK_ROOT.mkdir(parents=True, exist_ok=True)
        jdk.install(str(JDK_VERSION), path=JDK_ROOT.as_posix())
        print('Done')


def is_picard_installed():
    try:
        _, stderr = io.run_subprocess([f'{get_jdk_path()}/java', '-jar', PICARD_JAR, 'SortVcf', '--version'])
        return len(stderr) >= 1 and 'Version' in stderr[0]
    except FileNotFoundError:
        return False


def install_picard():
    install_jdk()
    if is_picard_installed():
        print("Found Picard installation")
        return
    picard_url = 'https://github.com/broadinstitute/picard/releases/latest/download/picard.jar'
    print(f'downloading picard.jar from {picard_url} to {PICARD_JAR}...')
    with urlopen(picard_url) as (response), open(PICARD_JAR, 'wb') as (f):
        copyfileobj(response, f)
    print('Done')


def _run_r_installer(script_name: str, package_name: str,
                     r_installation_folder: Union[str, Path, Literal['auto']] = 'auto'):
    """
    Run one of the bundled R package-installation scripts, and turn a failed installation into
    actionable guidance.

    ``io.run_r_script`` raises :class:`ChildProcessError` when the R script itself fails, which on
    fresh machines and CI is most often a missing write permission to R's library folder. The
    original error is chained (``from``) so the R output stays visible. Failures that indicate an
    RNAlysis bug rather than an environment problem (for example a missing bundled script, which
    raises :class:`~rnalysis.exceptions.InternalError`) are deliberately left to propagate as-is.
    """
    script_path = Path(__file__).parent.parent.joinpath(f'data_files/r_templates/{script_name}')
    try:
        io.run_r_script(script_path, r_installation_folder)
    except ChildProcessError as e:
        raise ChildProcessError(f"Failed to install {package_name}. "
                                f"Please make sure you have write permission to R's library folder, "
                                f"or try to install {package_name} manually.") from e


def install_limma(r_installation_folder: Union[str, Path, Literal['auto']] = 'auto'):
    _run_r_installer('limma_install.R', 'limma', r_installation_folder)


def install_deseq2(r_installation_folder: Union[str, Path, Literal['auto']] = 'auto'):
    _run_r_installer('deseq2_install.R', 'DESeq2', r_installation_folder)


def install_rsubread(r_installation_folder: Union[str, Path, Literal['auto']] = 'auto'):
    _run_r_installer('rsubread_install.R', 'RSubread', r_installation_folder)
