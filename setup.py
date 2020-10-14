# -*- coding: utf-8 -*-
"""This module contains the setup for the installation of MASSpy."""
from sys import argv

from setuptools import find_packages, setup


setup_kwargs = dict()
setup_requirements = []
# prevent pytest-runner from being installed on every invocation
if {'pytest', 'test', 'ptr'}.intersection(argv):
    setup_requirements.append("pytest-runner")

extras = {
    "escher": "escher",
    "json": "simplejson",
}
extras["all"] = sorted(extras.values())

try:
    with open('README.rst') as handle:
        readme = handle.read()
    with open('INSTALL.rst') as handle:
        install = handle.read()
    setup_kwargs["long_description"] = readme + "\n\n" + install
except IOError:
    setup_kwargs["long_description"] = ''

if __name__ == "__main__":
    setup(
        name="masspy",
        version="0.1.1",
        packages=find_packages(),
        setup_requires=setup_requirements,
        install_requires=[
            "cobra~=0.19.0",
            "depinfo~=1.5",
            "escher~=1.7",       # Optional for network visualiation
            "libroadrunner~=2.0",
            "matplotlib~=3.2",
            "numpy~=1.19",
            "pandas~=1.0",
            "scipy~=1.2",
            "six~=1.11",
            "sympy~=1.0",
            "tabulate~=0.8"
        ],
        # TODO Add when implemented
        # tests_require=[
        #     "pytest",
        #     "pytest-benchmark"
        # ],
        extras_require=extras,
        package_data={
            '': [
                'test/data/*',
            ]
        },
        author="Zachary B. Haiman",
        author_email="zhaiman@eng.ucsd.edu",
        maintainer="Zachary B. Haiman",
        maintainer_email="zhaiman@eng.ucsd.edu",
        description="MASSpy is a package for kinetic modeling and simulation "
                    "of biological networks",
        license="MIT",
        keywords=("metabolism biology kinetic modeling simulation programming"
                  "cobra"),
        url="https://github.com/SBRG/masspy",
        # TODO test_suite="mass.test.suite", Add when implemented
        download_url="https://pypi.python.org/pypi/masspy",
        classifiers=[
            'Development Status :: 5 - Production/Stable',
            'Environment :: Console',
            'Intended Audience :: Developers',
            'Topic :: Software Development :: Build Tools',
            'Intended Audience :: Science/Research',
            'License :: OSI Approved :: GNU Lesser General Public License v3'
                ' or later (LGPLv3+)',
            'License :: OSI Approved :: GNU General Public License v3'
                ' or later (GPLv3+)',
            'Operating System :: OS Independent',
            'Programming Language :: Python :: 3.6',
            'Programming Language :: Python :: 3.7',
            'Programming Language :: Python :: Implementation :: CPython',
            'Topic :: Scientific/Engineering',
            'Topic :: Scientific/Engineering :: Bio-Informatics'
        ],
        python_requires='>=3.6',
        include_package_data=True,
        platforms="GNU/Linux, Mac OS X >= 10.7, Microsoft Windows >= 7",
        **setup_kwargs
    )
