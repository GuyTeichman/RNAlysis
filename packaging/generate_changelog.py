import re
import warnings
from pathlib import Path
from typing import Literal, Union

import rst_to_myst


def get_change_log_for(version: Union[str, Literal['latest']] = 'latest',
                       history_path: Union[str, Path, None] = None):
    regex_pattern = r'(\d+\.\d+\.\d+)\s+\(\d{4}-\d{2}-\d{2}\)\s*-+\s*'
    if history_path is None:
        history_path = Path(__file__).parent.parent.joinpath('HISTORY.rst')
    with open(history_path) as hfile:
        text = hfile.read()

    versions = re.findall(regex_pattern, text)
    if version == 'latest':
        # guard: 'latest' must never silently fall back to the previous version because the
        # newest section is still undated (e.g. "4.3.0 (unreleased)")
        newest = re.search(r'(\d+\.\d+\.\d+)\s+\(([^)\n]*)\)\s*\n-+', text)
        if newest is not None and re.fullmatch(r'\d{4}-\d{2}-\d{2}', newest.group(2)) is None:
            raise ValueError(f"The top HISTORY.rst section ('{newest.group(1)} ({newest.group(2)})') is "
                             f"unreleased/undated -- date it before generating the changelog.")
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
    with open(Path(__file__).parent.parent.joinpath('rnalysis/data_files/latest_changelog.md'), 'w') as f:
        f.write(txt)
