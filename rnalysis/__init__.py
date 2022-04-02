# -*- coding: utf-8 -*-

"""Top-level package for RNA sequencing analysis pipeline."""
import warnings

__all__ = ['general', 'filtering', 'enrichment', '__attr_file_key__', '__biotype_file_key__']
__name__ = "rnalysis"
__author__ = """Guy Teichman"""
__email__ = "guyteichman@gmail.com"
__version__ = "2.0.1"
__license__ = "MIT"
__attr_file_key__ = "attribute_reference_table"""
__biotype_file_key__ = "biotype_reference_table"


def _simple_warning_format(msg, *args, **kwargs):
    # ignore everything except the warning message
    return f"Warning: {msg}\n"


warnings.formatwarning = _simple_warning_format
warnings.simplefilter('always', UserWarning)
