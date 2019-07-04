# Configuration file for the Sphinx documentation builder.
#
# This file only contains a selection of the most common options. For a full
# list see the documentation:
# http://www.sphinx-doc.org/en/master/config

# -- Path setup ---------------------------------------------------------------

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.
#
# import os
# import sys
# sys.path.insert(0, os.path.abspath('.'))


# -- Project information ------------------------------------------------------

project = 'masspy'
copyright = '2019, Z. Haiman'
author = 'Z. Haiman'

# The full version, including alpha/beta/rc tags
version = '0.1.0a38'
release = '0.1.0a38'


# -- General configuration ----------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = [
    'nbsphinx',
    'sphinx.ext.mathjax',
]
# Add any paths that contain templates here, relative to this directory.
templates_path = ['_templates']

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = ['_build', '**.ipynb_checkpoints']

# The master toctree document.
master_doc = 'index'

# -- Options for HTML output --------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
html_theme = 'alabaster'

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
