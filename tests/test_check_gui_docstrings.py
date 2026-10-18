"""Tests for ``packaging/check_gui_docstrings.py``, the ``check-gui-docstrings`` pre-commit hook.

The hook statically (via ``ast``) checks that every parameter the GUI renders for a public API function has
help text, as parsed by the same parser the GUI uses (``rnalysis/utils/docstring_parsing.py``).
It runs in pre-commit's isolated environment, where RNAlysis's dependencies are not installed, so it must
never import the ``rnalysis`` package. ``packaging/`` is not an importable package, so the script is loaded
from its file path (see ``tests/test_bench_equal.py``).
"""

import importlib.util
import subprocess
import sys
import textwrap
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parent.parent
SCRIPT_PATH = REPO_ROOT / 'packaging' / 'check_gui_docstrings.py'
PARSER_PATH = REPO_ROOT / 'rnalysis' / 'utils' / 'docstring_parsing.py'


def _load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


@pytest.fixture(scope='module')
def hook():
    return _load('rnalysis_check_gui_docstrings', SCRIPT_PATH)


def _check(hook, source, module_functions=True, root_classes=('Filter',), path='rnalysis/example.py'):
    return hook.check_source(textwrap.dedent(source), path, module_functions, root_classes)


def _kinds(problems):
    return sorted((p.qualname, p.param, p.kind) for p in problems)


GOOD_FUNC = '''
def good(a: int, b: str = 'x', *args, **kwargs):
    """
    Description.

    :param a: the first parameter; it has {special} characters.
    :type a: int
    :param b: the second parameter
    :type b: str
    :return: nothing
    """
'''


def test_good_function_passes(hook):
    assert _check(hook, GOOD_FUNC) == []


def test_missing_param_help(hook):
    src = '''
    def func(a: int, b: int):
        """
        Description.

        :param a: documented
        :type a: int
        :return: nothing
        """
    '''
    assert _kinds(_check(hook, src)) == [('func', 'b', hook.MISSING)]


def test_empty_param_help(hook):
    src = '''
    def func(a: int):
        """
        Description.

        :param a:
        :type a: int
        """
    '''
    assert _kinds(_check(hook, src)) == [('func', 'a', hook.MISSING)]


def test_unknown_param_typo(hook):
    src = '''
    def func(organism: str):
        """
        Description.

        :param organsim: the organism
        :type organism: str
        """
    '''
    assert _kinds(_check(hook, src)) == [('func', 'organism', hook.MISSING), ('func', 'organsim', hook.UNKNOWN)]


def test_documented_varargs_are_not_typos(hook):
    src = '''
    def func(*others, **function_kwargs):
        """
        Description.

        :param others: other objects
        :param function_kwargs: keyword arguments
        """
    '''
    assert _check(hook, src) == []


@pytest.mark.parametrize('bad_line', ['param a: no leading colon', ':param a no trailing colon', ':param a : spaced'])
def test_malformed_param_field(hook, bad_line):
    src = f'''
    def func(a: int):
        """
        Description.

        {bad_line}
        :type a: int
        """
    '''
    kinds = _kinds(_check(hook, src))
    assert ('func', 'a', hook.MALFORMED) in kinds


def test_missing_blank_line_before_params(hook):
    src = '''
    def func(a: int, b: int):
        """
        Description directly followed by params.
        :param a: first
        :param b: second
        """
    '''
    # one clear problem for the function, instead of one "missing help" per swallowed parameter
    assert _kinds(_check(hook, src)) == [('func', None, hook.NO_BLANK_LINE)]


def test_private_and_ignored_names_are_skipped(hook):
    src = '''
    def _private(a):
        """Description.

        """

    class Filter:
        def _private_method(self, a):
            """Description."""

        def method(self, a, gui_mode: bool = False, legacy_args: bool = False):
            """
            Description.

            :param a: documented
            """

        @classmethod
        def build(cls, a):
            """
            Description.

            :param a: documented
            """

        @staticmethod
        def static(a):
            """
            Description.

            :param a: documented
            """

        @property
        def prop(self):
            return 1
    '''
    assert _check(hook, src) == []


def test_staticmethod_first_param_is_checked(hook):
    src = '''
    class Filter:
        @staticmethod
        def static(a):
            """
            Description.

            :return: nothing
            """
    '''
    assert _kinds(_check(hook, src)) == [('Filter.static', 'a', hook.MISSING)]


def test_only_gui_classes_are_checked(hook):
    src = '''
    class Filter:
        def method(self, a):
            """
            Description.

            :return: nothing
            """

    class CountFilter(Filter):
        def other(self, b):
            """
            Description.

            :return: nothing
            """

    class Pipeline:
        def add(self, c):
            """
            Description.

            :return: nothing
            """
    '''
    assert _kinds(_check(hook, src)) == [('CountFilter.other', 'b', hook.MISSING), ('Filter.method', 'a', hook.MISSING)]


def test_module_functions_can_be_skipped(hook):
    src = '''
    def func(a):
        """
        Description.

        :return: nothing
        """
    '''
    assert _check(hook, src, module_functions=False) == []


def test_inherited_docstring_is_used(hook):
    src = '''
    class Filter:
        def method(self, a):
            """
            Description.

            :param a: documented in the base class
            """

    class CountFilter(Filter):
        def method(self, a):
            return 1

        def method_with_new_param(self, a, b):
            return 1
    '''
    # CountFilter.method inherits Filter.method's docstring (as inspect.getdoc() would);
    # method_with_new_param has no docstring anywhere and no readable_name, so it is skipped
    assert _check(hook, src) == []


def test_inherited_docstring_missing_new_param(hook):
    src = '''
    class Filter:
        def method(self, a):
            """
            Description.

            :param a: documented in the base class
            """

    class CountFilter(Filter):
        def method(self, a, b):
            return 1
    '''
    assert _kinds(_check(hook, src)) == [('CountFilter.method', 'b', hook.MISSING)]


def test_gui_function_without_docstring(hook):
    src = """
    class Filter:
        @readable_name('Do something')
        def method(self, a):
            return 1

        @readable_name('No parameters, no help needed')
        def method_no_params(self):
            return 1
    """
    assert _kinds(_check(hook, src)) == [('Filter.method', None, hook.NO_DOCSTRING)]


def test_docstring_possibly_inherited_from_another_module_is_skipped(hook):
    src = """
    class Filter(generic.GenericBase):
        @readable_name('Do something')
        def method(self, a):
            return 1
    """
    # inspect.getdoc() may find a docstring in generic.GenericBase, which the static check cannot see
    assert _check(hook, src) == []


def test_builtin_base_does_not_hide_undocumented_gui_method(hook):
    src = """
    class FeatureSet(set):
        @readable_name('Do something')
        def new_method(self, a: int = 1):
            pass

        def update(self, other):
            pass
    """
    # `set` has no `new_method`, so nothing can be inherited; `set.update` may supply a docstring
    assert _kinds(_check(hook, src, root_classes=('FeatureSet',))) == [
        ('FeatureSet.new_method', None, hook.NO_DOCSTRING)
    ]


def test_module_functions_only_hide_legacy_args(hook):
    src = '''
    def func(gui_mode: bool = False, legacy_args: bool = False):
        """
        Description.

        :return: nothing
        """
    '''
    # FuncExternalWindow renders gui_mode for a module-level function; only tab-page methods hide it
    assert _kinds(_check(hook, src)) == [('func', 'gui_mode', hook.MISSING)]


def test_prose_line_starting_with_param_is_not_malformed(hook):
    src = '''
    def func(a: int):
        """
        Description.

        :param a: the first parameter. Its
        param values: are listed elsewhere
        """
    '''
    assert _check(hook, src) == []


def test_typed_param_field_names_the_right_param(hook):
    src = '''
    def func(x: int):
        """
        Description.

        :param int x: help
        """
    '''
    assert _kinds(_check(hook, src)) == [('func', 'x', hook.MALFORMED)]


def test_problem_message_is_informative(hook):
    src = '''
    def first():
        pass


    def func(a: int, b: int):
        """
        Description.

        :param a: documented
        :param bb: typo
        """
    '''
    problems = _check(hook, src)
    messages = {p.param: str(p) for p in problems}
    assert messages['b'].startswith("rnalysis/example.py:6: func: parameter 'b'")
    assert messages['bb'].startswith("rnalysis/example.py:11: func: parameter 'bb'")


def test_main_ignores_unrelated_files(hook, capsys):
    assert hook.main(['README.rst', 'rnalysis/utils/io.py']) == 0


def test_main_current_tree_passes(hook, capsys):
    assert hook.main(sorted(hook.CHECKED_MODULES)) == 0, capsys.readouterr().out


def test_main_checks_everything_when_the_hook_changes(hook, capsys, monkeypatch):
    checked = []
    monkeypatch.setattr(hook, 'KNOWN_ISSUES', set())
    monkeypatch.setattr(hook, 'check_file', lambda path: checked.append(path) or [])
    assert hook.main(['packaging/check_gui_docstrings.py']) == 0
    assert sorted(checked) == sorted(hook.CHECKED_MODULES)


def test_main_reports_failures(hook, capsys, monkeypatch):
    problem = hook.Problem('rnalysis/fastq.py', 3, 'func', 'a', hook.MISSING)
    monkeypatch.setattr(hook, 'check_file', lambda path: [problem] if path == 'rnalysis/fastq.py' else [])
    assert hook.main(['rnalysis/fastq.py']) == 1
    assert "rnalysis/fastq.py:3: func: parameter 'a'" in capsys.readouterr().out


def test_known_issues_are_tolerated_and_stale_ones_reported(hook, capsys, monkeypatch):
    known = hook.Problem('rnalysis/fastq.py', 3, 'func', 'a', hook.MISSING)
    monkeypatch.setattr(hook, 'KNOWN_ISSUES', {known.key(), ('rnalysis/fastq.py', 'gone', 'b', hook.MISSING)})
    monkeypatch.setattr(hook, 'check_file', lambda path: [known] if path == 'rnalysis/fastq.py' else [])
    assert hook.main(['rnalysis/fastq.py']) == 1
    out = capsys.readouterr().out
    assert "parameter 'a'" not in out
    assert 'gone' in out and 'KNOWN_ISSUES' in out


def test_known_issues_all_exist(hook):
    # every baseline entry must point at a module the hook checks
    assert {key[0] for key in hook.KNOWN_ISSUES} <= set(hook.CHECKED_MODULES)


def test_hook_runs_without_rnalysis_dependencies(tmp_path):
    # pre-commit runs the hook in an isolated venv with no polars/numpy/PyQt: '-S' drops site-packages,
    # and running from a foreign cwd means `import rnalysis` cannot find the checkout either
    proc = subprocess.run(
        [sys.executable, '-S', str(SCRIPT_PATH), 'rnalysis/filtering.py'],
        cwd=REPO_ROOT,
        capture_output=True,
        text=True,
    )
    assert proc.returncode == 0, proc.stdout + proc.stderr


def test_docstring_parsing_is_stdlib_only():
    code = (
        'import importlib.util, sys\n'
        f'spec = importlib.util.spec_from_file_location("dp", r"{PARSER_PATH}")\n'
        'module = importlib.util.module_from_spec(spec)\n'
        'spec.loader.exec_module(module)\n'
        'assert not any(m.split(".")[0] in {"rnalysis", "polars", "numpy", "pandas"} for m in sys.modules)\n'
        'print(module.parse_docstring("Desc.\\n\\n:param a: text\\n:type a: int"))\n'
    )
    proc = subprocess.run([sys.executable, '-S', '-c', code], capture_output=True, text=True)
    assert proc.returncode == 0, proc.stderr
    assert proc.stdout.strip() == "('Desc.', {'a': ' text\\n'})"
