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

try:
    with open('README.rst') as handle:
        readme = handle.read()
    with open('INSTALL.rst') as handle:
        install = handle.read()
    setup_kwargs["long_description"] = readme + "\n\n" + install
except IOError:
    setup_kwargs["long_description"] = ''

setup(
    name="masspy",
    version="0.1.0a2",
    description="MASSpy is a package for kinetic modeling and simulation of "
                "biological networks",
    license="LGPL/GPL v2+",
    url="https://github.com/SBRG/MASSpy",
    classifiers=[
        'Development Status :: 3 - Alpha',
        'Environment :: Console',
        'Intended Audience :: Developers',
        'Topic :: Software Development :: Build Tools',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: GNU Lesser General Public License v2'
            ' or later (LGPLv2+)',
        'License :: OSI Approved :: GNU General Public License v2'
            ' or later (GPLv2+)',
        'Operating System :: OS Independent',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 3.6',
        'Topic :: Scientific/Engineering',
        'Topic :: Scientific/Engineering :: Bio-Informatics'
    ],
    keywords=("metabolism biology kinetic modeling simulation programming"),
    packages=find_packages(),
    python_requires='>=3.6',
    setup_requires=setup_requirements,
    install_requires=[
        "six",
        "future",
        "cobra>=0.8.2",
        "pandas>=0.21.0",
        "numpy>=1.13.1",
        "scipy>=0.19.1",
        "sympy>=1.0",
        "matplotlib>=2.1.1",
        "tabulate",
        "python-libsbml", # Move to extras eventually
        "lxml"            # Move to extras eventually
    ],
    tests_require=[],
    **setup_kwargs
)

# To add at a later date:
# download_url
# author
# author_email
# platforms=""
# extras_require
