"""Import-regression guards: the heavy scientific stack must stay lazily imported.

``import rnalysis.filtering`` used to cost ~3.5s warm because seaborn, scikit-learn and pandas were
imported at module level (issue #257). Every script, notebook kernel, CI test process and -- worst of
all -- every ``spawn`` multiprocessing worker paid that tax before computing anything. The SPEC 1
lazy-loading conversion (``lazy_loader``) removed it, and these tests exist so it cannot silently
come back.

Each check runs in a **subprocess**, because the pytest process has long since imported pandas &
friends through other test modules. ``lazy_loader.load()`` registers a lazy *stub* under the
package's top-level name in ``sys.modules``, so "``'seaborn' in sys.modules``" is not by itself
evidence of an eager import -- what matters is whether the package body actually ran, which the
helper below detects both by the stub's type and by the presence of any of its submodules.
"""

import json
import subprocess
import sys
import textwrap
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parent.parent

# Packages that dominated the import cost and must not be imported eagerly by the API modules.
# matplotlib/matplotlib.pyplot are deliberately NOT listed: they stay eager (see issue #257) because
# the GUI's backend selection and the many ``plt.Figure`` return annotations depend on them.
LAZY_PACKAGES = ('pandas', 'seaborn', 'sklearn')

_PROBE = textwrap.dedent("""
    import json, sys
    {import_statement}
    packages = {packages!r}
    states = {{}}
    for name in packages:
        module = sys.modules.get(name)
        if module is None:
            states[name] = 'absent'
        elif type(module).__name__ == '_LazyModule':
            # importlib.util.LazyLoader's placeholder -- the package body has not run yet
            states[name] = 'lazy'
        else:
            states[name] = 'imported'
    # a package's body cannot run without pulling in submodules, so any submodule in sys.modules is
    # independent proof of a real (eager) import, whatever the top-level entry looks like
    loaded_submodules = sorted(name for name in sys.modules if name.split('.')[0] in packages
                               and '.' in name)
    print(json.dumps({{'states': states, 'loaded_submodules': loaded_submodules}}))
""")


def _probe_imports(import_statement: str) -> dict:
    """Run ``import_statement`` in a fresh interpreter and report the state of each heavy package."""
    code = _PROBE.format(import_statement=import_statement, packages=LAZY_PACKAGES)
    process = subprocess.run([sys.executable, '-c', code], capture_output=True, text=True, cwd=str(REPO_ROOT))
    assert process.returncode == 0, f'{import_statement!r} failed:\n{process.stderr}'
    return json.loads(process.stdout.splitlines()[-1])


@pytest.mark.parametrize(
    'import_statement',
    [
        'import rnalysis.filtering',
        'from rnalysis import filtering',
        'import rnalysis.enrichment',
        'from rnalysis import enrichment',
    ],
)
def test_importing_the_api_does_not_import_heavy_packages(import_statement):
    result = _probe_imports(import_statement)
    eager = {name: state for name, state in result['states'].items() if state == 'imported'}
    assert eager == {}, f'{import_statement!r} eagerly imported {sorted(eager)}'
    assert result['loaded_submodules'] == [], f'{import_statement!r} eagerly imported {result["loaded_submodules"]}'


@pytest.mark.parametrize('module,attribute', [('filtering', 'CountFilter'), ('enrichment', 'FeatureSet')])
def test_importing_the_api_still_works(module, attribute):
    """The lazy proxies must not break the module at import time (nothing may touch their attributes)."""
    code = f'import rnalysis.{module} as m; print(m.{attribute}.__name__)'
    process = subprocess.run([sys.executable, '-c', code], capture_output=True, text=True, cwd=str(REPO_ROOT))
    assert process.returncode == 0, process.stderr
    assert process.stdout.strip() == attribute
