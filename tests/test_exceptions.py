import pytest

from rnalysis.exceptions import InvalidTypeError, InvalidValueError, RNAlysisInputError


def test_root_is_an_exception():
    assert issubclass(RNAlysisInputError, Exception)


@pytest.mark.parametrize('cls,builtin', [(InvalidTypeError, TypeError), (InvalidValueError, ValueError)])
def test_subclasses_inherit_root_and_builtin(cls, builtin):
    assert issubclass(cls, RNAlysisInputError)
    assert issubclass(cls, builtin)


@pytest.mark.parametrize('cls', [RNAlysisInputError, InvalidTypeError, InvalidValueError])
def test_message_is_preserved(cls):
    err = cls('some message')
    assert str(err) == 'some message'


@pytest.mark.parametrize('cls,other_builtin', [(InvalidTypeError, ValueError), (InvalidValueError, TypeError)])
def test_subclasses_do_not_cross_inherit(cls, other_builtin):
    assert not issubclass(cls, other_builtin)


def test_module_imports_nothing_heavy():
    """The exceptions module must stay import-cheap: no third-party imports at all."""
    import ast
    from pathlib import Path

    import rnalysis.exceptions

    tree = ast.parse(Path(rnalysis.exceptions.__file__).read_text(encoding='utf-8'))
    imported = [node for node in ast.walk(tree) if isinstance(node, (ast.Import, ast.ImportFrom))]
    assert imported == []
