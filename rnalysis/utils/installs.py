import os.path
import subprocess
from shutil import copyfileobj
from urllib.request import urlopen

import jdk

from rnalysis.utils import io

PICARD_JAR = os.environ.get('PICARDTOOLS_JAR', os.path.join(os.path.dirname(__file__), 'picard.jar'))
JDK_PATH = io.get_data_dir().joinpath('jdk17')


def is_jdk_installed():
    try:
        if not JDK_PATH.exists():
            return False
        # Run the "java -version" command and capture the output
        output = subprocess.check_output([f'{JDK_PATH}/java', "-version"], stderr=subprocess.STDOUT, text=True)
        # Check if the output contains "version 17" or "17." (for Java 17)
        if "version 17" in output or " 17." in output:
            return True
    except subprocess.CalledProcessError:
        # If the "java -version" command returns an error, Java 17 is not installed.
        return False
    return False


def install_jdk():
    if not is_jdk_installed():
        print("Installing Java...")
        if not JDK_PATH.exists():
            JDK_PATH.mkdir(parents=True)
        jdk.install('17', path=JDK_PATH)
        print('Done')


def is_picard_installed():
    try:
        _, stderr = io.run_subprocess([f'{JDK_PATH}/java', '-jar', PICARD_JAR, 'SortVcf', '--version'])
        return len(stderr) >= 1 and 'Version' in stderr[0]
    except FileNotFoundError:
        return False


def install_picard():
    install_jdk()
    if is_picard_installed():
        print("Picard is already installed")
        return
    picard_url = 'https://github.com/broadinstitute/picard/releases/latest/download/picard.jar'
    print(f'downloading picard.jar from {picard_url} to {PICARD_JAR}...')
    with urlopen(picard_url) as (response), open(PICARD_JAR, 'wb') as (f):
        copyfileobj(response, f)
    print('Done')
