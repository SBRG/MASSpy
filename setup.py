# -*- coding: utf-8 -*-

from __future__ import absolute_import

from warnings import warn
from sys import argv, version_info

from setuptools import setup, find_packages

if version_info[:2] == (3, 4):
    warn("Support for Python 3.4 was dropped by pandas. Since cobrapy is a "
         "pure Python package you can still install it but will have to "
         "carefully manage your own pandas and numpy versions. We no longer "
         "include it in our automatic testing.")

setup_kwargs = dict()
setup_requirements = []

extras = {
    'sbml': ["python-libsbml", "lxml"]
}
extras["all"] = sorted(list(extras))

try:
    with open('README.rst') as handle:
        readme = handle.read()
    with open('INSTALL.rst') as handle:
        install = handle.read()
    setup_kwargs["long_description"] = readme + "\n\n" + install
except IOError:
    setup_kwargs["long_description"] = ''
    
