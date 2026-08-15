"""
Exception types deliberately raised by RNAlysis.

This module imports nothing at all (not even from the standard library), so that it stays cheap to import from
anywhere in the package, including hot paths and modules that are imported at startup.

Every exception raised on purpose by RNAlysis inherits from :class:`RNAlysisError`, and falls into one of two
families:

* **User-input errors** (:class:`RNAlysisInputError`) - the user supplied something RNAlysis cannot work with. \
  These also inherit from the matching Python built-in, so code that already handles those built-ins keeps working:
  a wrong-**type** argument raises :class:`InvalidTypeError` (a :class:`TypeError`), and a bad-**value** argument \
  (out of range, not a legal choice, wrong shape) raises :class:`InvalidValueError` (a :class:`ValueError`)::

      try:
          counts.filter_low_reads(threshold=-5)
      except RNAlysisInputError as err:
          print(err)

  A broken *file* the user asked RNAlysis to open belongs to this family too: \
  :class:`CorruptSessionError` is raised when a session file cannot be loaded because it is \
  corrupt or incomplete.

* **Internal invariants** (:class:`InternalError`) - a condition that cannot be false unless RNAlysis itself has a \
  bug. These used to be bare `assert` statements, which meant they vanished under `python -O`; they are now real \
  exceptions, so a violated invariant always fails loudly instead of silently corrupting an analysis.

None of these inherit from :class:`AssertionError`. Code that used to catch `AssertionError` from RNAlysis must
catch :class:`RNAlysisError` (or one of its subclasses) instead.
"""

#: appended to every :class:`InternalError` message, since reaching one always means RNAlysis has a bug
_BUG_REPORT_SUFFIX = (
    'This is likely a bug in RNAlysis - please report it at https://github.com/GuyTeichman/RNAlysis/issues'
)


class RNAlysisError(Exception):
    """Root for all exceptions deliberately raised by RNAlysis."""


class RNAlysisInputError(RNAlysisError):
    """Root for all user-input validation errors raised by RNAlysis."""


class InvalidTypeError(RNAlysisInputError, TypeError):
    """Raised when a user-supplied argument is of the wrong type."""


class InvalidValueError(RNAlysisInputError, ValueError):
    """Raised when a user-supplied argument has an illegal value, shape, or is not one of the legal choices."""


class CorruptSessionError(RNAlysisInputError, ValueError):
    """
    Raised when an RNAlysis session file cannot be loaded because it is corrupt or incomplete.

    A truncated, hand-edited, or partially-downloaded session file is a broken *input*, not a
    violated internal invariant - loading one must therefore not ask the user to report a bug.
    """


class InternalError(RNAlysisError, RuntimeError):
    """An internal invariant was violated - indicates a bug in RNAlysis itself."""

    def __init__(self, *args):
        # the suffix is baked into args (not just into __str__) because the GUI's error dialog
        # renders exception.args, so this is the only way the bug-report line reaches the user
        self._raw_args = args
        if args:
            first = f'{args[0]} {_BUG_REPORT_SUFFIX}' if str(args[0]) else _BUG_REPORT_SUFFIX
            args = (first,) + args[1:]
        else:
            args = (_BUG_REPORT_SUFFIX,)
        super().__init__(*args)

    def __reduce__(self):
        # BaseException's default __reduce__ returns (cls, self.args), which would feed the
        # already-suffixed message back through __init__ and append a second copy of the suffix.
        # joblib pickles worker exceptions back to the parent under the loky/multiprocessing
        # backends, so without this an error would gain one suffix per hop. Reconstruct from the
        # original arguments instead, which makes round-trips idempotent.
        return self.__class__, self._raw_args
