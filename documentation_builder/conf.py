# Configuration file for the Sphinx documentation builder.
#
# This file only contains a selection of the most common options. For a full
# list see the documentation:
# http://www.sphinx-doc.org/en/master/config

# -- Path setup ---------------------------------------------------------------

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute.

import sys
from os.path import dirname

from mass import __version__

sys.path.insert(0, dirname(dirname(__file__)))


# -- Project information ------------------------------------------------------

project = 'masspy'
copyright = '2019, Z. Haiman'
author = 'Z. Haiman'

# The full version, including alpha/beta/rc tags
release = __version__
version = '.'.join(release.split('.')[:2])

# -- General configuration ----------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = [
    'sphinx.ext.autodoc',
    'sphinx.ext.intersphinx',
    'sphinx.ext.mathjax',
    'sphinx.ext.viewcode',
    'sphinx.ext.napoleon',
    'sphinx.ext.autosummary',
    'sphinx.ext.autosectionlabel',
    'autoapi.extension',
    'nbsphinx'
]

# Document Python Code
autoapi_type = 'python'
autoapi_dirs = ['../mass']
autoapi_ignore = [
    '.tox', '.pytest_cache', 'scripts', 'benchmarks', "notebooks"]

# Napoleon settings
napoleon_google_docstring = False
napoleon_numpy_docstring = True

# Add any paths that contain templates here, relative to this directory.
templates_path = ['_templates']

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = ['_build', '**.ipynb_checkpoints', 'Thumbs.db', '.DS_Store']

# The master toctree document.
master_doc = 'index'

# If true, `todo` and `todoList` produce output, else they produce nothing.
todo_include_todos = False

pygments_style = 'sphinx'

# In order to build documentation that requires libraries to import
class Mock(object):
    def __init__(self, *args, **kwargs):
        return

    def __call__(self, *args, **kwargs):
        return Mock()

    @classmethod
    def __getattr__(cls, name):
        if name in ('__file__', '__path__'):
            return '/dev/null'
        else:
            return Mock()

# These modules should correspond to the importable Python packages.
MOCK_MODULES = [
    
]
for mod_name in MOCK_MODULES:
    sys.modules[mod_name] = Mock()

# -- Options for HTML output --------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
html_theme = "sphinx_rtd_theme"

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
# html_static_path = ['_static']

# A list of paths that contain extra files not directly related to the
# documentation, such as robots.txt or .htaccess. Relative paths are taken as
# relative to the configuration directory. They are copied to the output
# directory. They will overwrite any existing file of the same name.
html_extra_path = ["robots.txt"]

# -- Options for LaTeX output -------------------------------------------------

# -- Options for manual page output -------------------------------------------

# -- Options for Texinfo output -----------------------------------------------
