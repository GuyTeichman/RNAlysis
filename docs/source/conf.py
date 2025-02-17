# Configuration file for the Sphinx documentation builder.
#
# This file only contains a selection of the most common options. For a full
# list see the documentation:
# http://www.sphinx-doc.org/en/master/config

# -- Path setup --------------------------------------------------------------

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.
#
import os
import sys

# sys.path.insert(0, os.path.abspath('.'))


sys.path.insert(0, os.path.abspath('../../rnalysis/gui/videos'))
sys.path.insert(0, os.path.abspath('../../rnalysis'))
sys.path.insert(0, os.path.abspath('../..'))

# -- Project information -----------------------------------------------------

project = 'RNAlysis'
copyright = '2024, Guy Teichman'
author = 'Guy Teichman'

# The full version, including alpha/beta/rc tags
release = '4.1.1'

# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = ['sphinx.ext.autodoc', 'sphinx.ext.autosummary', 'sphinx.ext.doctest']

# Add any paths that contain templates here, relative to this directory.
templates_path = ['_templates']

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = []


def hide_non_private(app, what, name, obj, skip, options):
    # if object is a class, doc it anyway
    if what == "method" or what == "function":
        # if private-members is set, show only private members
        if 'private-members' in options and not name.startswith('_'):
            # skip public methods
            return True
        elif 'private-members' in options and name.startswith('__'):
            return True
        else:
            return None
    elif what == "class":
        # do not modify skip - private methods will be shown
        return None


def setup(app):
    app.add_js_file('copybutton.js')
    app.connect('autodoc-skip-member', hide_non_private)


html_favicon = 'favicon.ico'
html_logo = 'logo.png'

# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
html_theme = 'sphinx_rtd_theme'
# html_theme = 'default'
# html_theme = 'alabaster'

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ['_static']

# automodapi options
automodsumm_inherited_members = True
automodapi_inheritance_diagram = False

autosummary_generate = True

html_theme_options = {
    'logo_only': True,
    'display_version': True,
    'prev_next_buttons_location': 'bottom',
}
