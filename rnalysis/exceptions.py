"""
Exception types deliberately raised by RNAlysis.

This module imports nothing at all (not even from the standard library), so that it stays cheap to import from
anywhere in the package, including hot paths and modules that are imported at startup.

Every exception raised on purpose by RNAlysis inherits from :class:`RNAlysisError`, and falls into one of three
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

* **External-service failures** (:class:`ExternalServiceError`) - a remote service RNAlysis depends on (UniProt, \
  Ensembl, PANTHER, PhylomeDB, OrthoInspector, KEGG, GO) was slow, changed, or went down. Neither the user's fault \
  nor an RNAlysis bug; code that hits one should degrade gracefully (a partial or empty result plus a clear \
  message), never crash. :class:`IDMappingTimeoutError` is raised when a UniProt idmapping job never becomes ready \
  within the allotted time.

None of these inherit from :class:`AssertionError`. Code that used to catch `AssertionError` from RNAlysis must
catch :class:`RNAlysisError` (or one of its subclasses) instead.
"""

#: appended to every :class:`InternalError` message, since reaching one always means RNAlysis has a bug
_BUG_REPORT_SUFFIX = ('This is likely a bug in RNAlysis - '
                      'please report it at https://github.com/GuyTeichman/RNAlysis/issues')


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


class ExternalServiceError(RNAlysisError):
    """
    Root for failures of a remote service RNAlysis depends on (UniProt, Ensembl, PANTHER, PhylomeDB,
    OrthoInspector, KEGG, GO).

    These are neither the user's fault (:class:`RNAlysisInputError`) nor an RNAlysis bug
    (:class:`InternalError`) - a third party was slow, changed, or went down. Code that raises one
    should have already degraded to a partial or empty result with a clear, user-facing message.
    """


class IDMappingTimeoutError(ExternalServiceError):
    """
    Raised when a UniProt gene-ID mapping job does not become ready within the allotted time.

    UniProt idmapping runs asynchronously and its job latency is entirely server-side; a job can
    occasionally stay in a running state for far longer than any legitimate job takes. Bounding the
    wait lets a wedged job degrade to a partial mapping (retried on the next run) instead of hanging
    the whole analysis indefinitely.
    """


class InternalError(RNAlysisError, RuntimeError):
    """An internal invariant was violated - indicates a bug in RNAlysis itself."""

    def __init__(self, *args):
        if args:
            first = f'{args[0]} {_BUG_REPORT_SUFFIX}' if str(args[0]) else _BUG_REPORT_SUFFIX
            args = (first,) + args[1:]
        else:
            args = (_BUG_REPORT_SUFFIX,)
        super().__init__(*args)
