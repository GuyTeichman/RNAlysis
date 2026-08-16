import pytest

from rnalysis.exceptions import InternalError, InvalidTypeError, InvalidValueError, RNAlysisError, RNAlysisInputError

BUG_REPORT_SUFFIX = (
    'This is likely a bug in RNAlysis - please report it at https://github.com/GuyTeichman/RNAlysis/issues'
)


def test_root_is_an_exception():
    assert issubclass(RNAlysisError, Exception)


@pytest.mark.parametrize('cls', [RNAlysisInputError, InvalidTypeError, InvalidValueError, InternalError])
def test_everything_inherits_the_root(cls):
    assert issubclass(cls, RNAlysisError)


@pytest.mark.parametrize(
    'cls,builtin', [(InvalidTypeError, TypeError), (InvalidValueError, ValueError), (InternalError, RuntimeError)]
)
def test_subclasses_inherit_matching_builtin(cls, builtin):
    assert issubclass(cls, builtin)


@pytest.mark.parametrize('cls', [InvalidTypeError, InvalidValueError])
def test_input_errors_inherit_input_root(cls):
    assert issubclass(cls, RNAlysisInputError)


def test_internal_error_is_not_an_input_error():
    assert not issubclass(InternalError, RNAlysisInputError)


@pytest.mark.parametrize('cls,other_builtin', [(InvalidTypeError, ValueError), (InvalidValueError, TypeError)])
def test_subclasses_do_not_cross_inherit(cls, other_builtin):
    assert not issubclass(cls, other_builtin)


@pytest.mark.parametrize('cls', [RNAlysisError, RNAlysisInputError, InvalidTypeError, InvalidValueError])
def test_message_is_preserved(cls):
    err = cls('some message')
    assert str(err) == 'some message'


@pytest.mark.parametrize('cls', [RNAlysisError, RNAlysisInputError, InvalidTypeError, InvalidValueError, InternalError])
def test_nothing_inherits_assertionerror(cls):
    """The exception-type change is a clean break: no compat shim through AssertionError."""
    assert not issubclass(cls, AssertionError)


def test_internal_error_appends_bug_report_line():
    err = InternalError('the invariant blew up.')
    assert str(err) == f'the invariant blew up. {BUG_REPORT_SUFFIX}'


def test_internal_error_without_message_is_just_the_bug_report_line():
    assert str(InternalError()) == BUG_REPORT_SUFFIX
    assert str(InternalError('')) == BUG_REPORT_SUFFIX


@pytest.mark.parametrize('message', ['the invariant blew up.', '', None])
def test_internal_error_survives_pickle_round_trip_without_doubling_the_suffix(message):
    """joblib pickles worker exceptions back to the parent under the loky/multiprocessing backends,
    so an InternalError raised in a CLICOM worker crosses at least one pickle boundary. If the suffix
    were re-appended on reconstruction it would accumulate one copy per hop."""
    import pickle

    err = InternalError() if message is None else InternalError(message)
    original = str(err)

    round_tripped = pickle.loads(pickle.dumps(err))
    assert str(round_tripped) == original
    assert str(round_tripped).count(BUG_REPORT_SUFFIX) == 1

    # several hops (worker -> parent -> re-raised) must stay stable too
    twice = pickle.loads(pickle.dumps(round_tripped))
    assert str(twice) == original
    assert str(twice).count(BUG_REPORT_SUFFIX) == 1


@pytest.mark.parametrize('copier', ['copy', 'deepcopy'])
def test_internal_error_survives_copy_round_trip_without_doubling_the_suffix(copier):
    import copy as copy_module

    err = InternalError('the invariant blew up.')
    original = str(err)

    copied = getattr(copy_module, copier)(err)
    assert str(copied) == original
    assert str(copied).count(BUG_REPORT_SUFFIX) == 1


def test_internal_error_keeps_the_suffix_in_args():
    """The GUI error dialog renders exception.args rather than str(exception), so the bug-report line
    has to live in args to reach the user."""
    err = InternalError('the invariant blew up.')
    assert BUG_REPORT_SUFFIX in ''.join(str(arg) for arg in err.args)


@pytest.mark.parametrize('cls', [RNAlysisError, RNAlysisInputError, InvalidTypeError, InvalidValueError])
def test_other_exceptions_survive_pickle_round_trip(cls):
    import pickle

    err = cls('some message')
    assert str(pickle.loads(pickle.dumps(err))) == 'some message'


def test_module_imports_nothing_heavy():
    """The exceptions module must stay import-cheap: no third-party imports at all."""
    import ast
    from pathlib import Path

    import rnalysis.exceptions

    tree = ast.parse(Path(rnalysis.exceptions.__file__).read_text(encoding='utf-8'))
    imported = [node for node in ast.walk(tree) if isinstance(node, (ast.Import, ast.ImportFrom))]
    assert imported == []
