import re
import warnings
from pathlib import Path
from typing import Union

import rst_to_myst
from typing import Literal


def get_change_log_for(version: Union[str, Literal['latest']] = 'latest'):
    regex_pattern = r'(\d+\.\d+\.\d+)\s+\(\d{4}-\d{2}-\d{2}\)\s*-+\s*'
    with open(Path(__file__).parent.parent.joinpath('HISTORY.rst')) as hfile:
        text = hfile.read()

    versions = re.findall(regex_pattern, text)
    if version == 'latest':
        version = versions[0]

    if version not in versions:
        warnings.warn(f'changelog not found for version {version}')
        return ''

    start = text.find(version)
    version_ind = versions.index(version)
    if version_ind == len(versions) - 1:
        end = len(text)
    else:
        end = text.find(versions[version_ind + 1])

    change_logs = text[start:end]
    return rst_to_myst.rst_to_myst(change_logs).text


if __name__ == '__main__':
    txt = get_change_log_for('latest')
    with open(Path(__file__).parent.parent.joinpath('latest_changelog.md'), 'w') as f:
        f.write(txt)
