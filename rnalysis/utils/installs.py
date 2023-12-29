import os
import os.path
import shutil
import subprocess
from pathlib import Path
from shutil import copyfileobj
from typing import Union, Literal
from urllib.request import urlopen

import jdk

from rnalysis.utils import io

PICARD_JAR = Path(os.environ.get('PICARDTOOLS_JAR', io.get_data_dir().joinpath('picard.jar')))
JDK_ROOT = io.get_data_dir().joinpath('jdk17')


def get_jdk_path():
    # List all files and directories in the given directory
    items = os.listdir(JDK_ROOT)

    # Filter for directories starting with 'jdk-'
    jdk_directories = sorted([item for item in items if
                              item.startswith('jdk-') and os.path.isdir(os.path.join(JDK_ROOT, item))], reverse=True)
    if len(jdk_directories) == 0:
        raise FileNotFoundError('No JDK directory found')
    return JDK_ROOT.joinpath(f'{jdk_directories[0]}/bin')


def is_jdk_installed():
    try:
        if not JDK_ROOT.exists():
            return False
        # Run the "java -version" command and capture the output
        jdk_path = get_jdk_path()
        output = subprocess.check_output([f'{jdk_path}/java', "-version"], stderr=subprocess.STDOUT, text=True)
        # Check if the output contains "version 17" or "17." (for Java 17)
        if "version 17" in output or " 17." in output:
            return True
    except (subprocess.CalledProcessError, FileNotFoundError):
        # If the "java -version" command returns an error, Java 17 is not installed.
        return False
    return False


def install_jdk():
    if not is_jdk_installed():
        print("Installing Java...")
        if JDK_ROOT.exists():
            shutil.rmtree(JDK_ROOT)
        JDK_ROOT.mkdir(parents=True)
        jdk.install('17', path=JDK_ROOT.as_posix())
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
        print("Picard is installed")
        return
    picard_url = 'https://github.com/broadinstitute/picard/releases/latest/download/picard.jar'
    print(f'downloading picard.jar from {picard_url} to {PICARD_JAR}...')
    with urlopen(picard_url) as (response), open(PICARD_JAR, 'wb') as (f):
        copyfileobj(response, f)
    print('Done')


def install_limma(r_installation_folder: Union[str, Path, Literal['auto']] = 'auto'):
    script_path = Path(__file__).parent.parent.joinpath('data_files/r_templates/limma_install.R')
    try:
        io.run_r_script(script_path, r_installation_folder)
    except AssertionError:
        raise AssertionError("Failed to install limma. "
                             "Please make sure you have write premission to R's library folder, "
                             "or try to install limma manually.")


def install_deseq2(r_installation_folder: Union[str, Path, Literal['auto']] = 'auto'):
    script_path = Path(__file__).parent.parent.joinpath('data_files/r_templates/deseq2_install.R')
    try:
        io.run_r_script(script_path, r_installation_folder)
    except AssertionError:
        raise AssertionError("Failed to install DESeq2. "
                             "Please make sure you have write premission to R's library folder, "
                             "or try to install DESeq2 manually.")
