"""
Exception types raised by RNAlysis when a user supplies an invalid input.

This module deliberately imports nothing at all (not even from the standard library), so that it stays cheap to
import from anywhere in the package, including hot paths and modules that are imported at startup.

All user-input validation in the public API (:mod:`rnalysis.filtering`, :mod:`rnalysis.enrichment`, \
:mod:`rnalysis.fastq`) raises one of the exceptions defined here. Every one of them inherits from \
:class:`RNAlysisInputError`, so API users can catch all RNAlysis input errors with a single `except` clause::

    try:
        counts.filter_low_reads(threshold=-5)
    except RNAlysisInputError as err:
        print(err)

They *also* inherit from the matching Python built-in (:class:`TypeError` / :class:`ValueError`), so code that
already handles those built-ins keeps working:

    * a wrong-**type** argument raises :class:`InvalidTypeError`, which is a :class:`TypeError`.
    * a bad-**value** argument (out of range, not a legal choice, wrong shape) raises :class:`InvalidValueError`, \
      which is a :class:`ValueError`.

Genuine internal invariants - conditions that cannot be false unless RNAlysis itself has a bug - are *not* \
expressed with these exceptions, and remain plain `assert` statements.
"""


class RNAlysisInputError(Exception):
    """Root for all user-input validation errors raised by RNAlysis."""


class InvalidTypeError(RNAlysisInputError, TypeError):
    """Raised when a user-supplied argument is of the wrong type."""


class InvalidValueError(RNAlysisInputError, ValueError):
    """Raised when a user-supplied argument has an illegal value, shape, or is not one of the legal choices."""
