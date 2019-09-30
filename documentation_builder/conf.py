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

sys.path.insert(0, dirname(dirname(__file__)))

from mass import __version__

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
    'nbsphinx',
    'sphinx.ext.viewcode',
    'sphinx.ext.mathjax',
    'sphinx.ext.napoleon',
    'sphinx.ext.intersphinx',
    'sphinx.ext.autodoc',
    'sphinx.ext.autosummary',
    'sphinx.ext.autosectionlabel',
    'autoapi.extension',
    'sphinxcontrib.bibtex',
]

# Automated documention of Python Code
autoapi_type = 'python'
autoapi_dirs = ['../mass']
autoapi_ignore = [
    '.tox', '.pytest_cache', 'scripts', 'benchmarks', "notebooks"]

autosectionlabel_prefix_document = True
autodoc_mock_imports = [
    "cobra",
    "depinfo",
    "libroadrunner",
    "matplotlib",
    "numpy",
    "pandas",
    "scipy",
    "sympy",
    "tabulate"
]

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


# Refer to the Python documentation for other libraries.
intersphinx_mapping = {
    "http://docs.python.org/": None,
    "https://cobrapy.readthedocs.io/en/latest/": None,
    "https://libroadrunner.readthedocs.io/en/latest/": None,
    "https://matplotlib.org/": None,
    "http://docs.scipy.org/doc/numpy/": None,
    "https://pandas.pydata.org/pandas-docs/stable/": None,
    "http://docs.scipy.org/doc/scipy/reference": None,
    "https://docs.sympy.org/latest/": None,
    }
intersphinx_cache_limit = 10  # days to keep the cached inventories
