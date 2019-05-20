# -*- coding: utf-8 -*-
"""TODO Module Docstrings."""
from __future__ import absolute_import

import argparse
import os
import sys
from collections import OrderedDict

from nbconvert.preprocessors import CellExecutionError, ExecutePreprocessor

import nbformat

from six import iterkeys


MAIN_DIR = os.path.dirname(os.path.abspath(__file__))
# Map acceptable category arguments to the corresponding notebook directories
CATEGORIES_DICT = OrderedDict({
    "SB2": "SB2",
    "MassModel": "test-models",
    "EnzymeModule": "test-models",
    "test": "test-models",
    "all": "."})
TEST_MODEL_SUFFIXES = ["MassModel", "EnzymeModule"]


def make_dir_path(category_to_get):
    return os.path.join(MAIN_DIR, CATEGORIES_DICT[category_to_get])


def filter_notebooks(category, dir_path):
    return sorted(list(filter(
        lambda x: category in x, os.listdir(dir_path))))


def get_notebook_paths(category_to_get):
    # Initialize Notebook container
    notebook_paths = []

    # Create iterable for notebook categories.
    if category_to_get == "all":
        category_list = [k for k in iterkeys(CATEGORIES_DICT)
                         if k != "test" and k != "all"]
    elif category_to_get == "test":
        category_list = TEST_MODEL_SUFFIXES
    else:
        category_list = [category_to_get]

    dir_list = [make_dir_path(c) for c in category_list]
    for category, dir_path in zip(category_list, dir_list):
        # Get names of notebooks
        notebooks = filter_notebooks(category, dir_path)
        # Make notebook paths and add to list
        notebook_paths += list(map(
            lambda x: os.path.join(dir_path, x), notebooks))

    return notebook_paths


def execute_notebooks(notebook_paths, verbose):
    for nb_path in notebook_paths:
        nb_dir = os.path.realpath(os.path.join(nb_path, ".."))
        nb_name = nb_path.replace(nb_dir + "/", "")
        if verbose:
            print("Executing notebook: {0}".format(nb_name))
        with open(nb_path) as f:
            nb = nbformat.read(f, as_version=4)

        ep = ExecutePreprocessor(timeout=600, kernel_name='python3')

        try:
            ep.preprocess(nb, {'metadata': {'path': nb_dir}})
        except CellExecutionError:
            msg = 'Error executing the notebook "{0}".\n\n'.format(nb_name)
            msg += 'See notebook "{0}" for the traceback.'.format(nb_name)
            print(msg)
            raise
        finally:
            with open(nb_path, mode='w', encoding='utf-8') as f:
                nbformat.write(nb, f)
    if verbose:
        print("Finished executing notebooks")


def main(args):
    # Determine notebook category
    category_to_get = args.category
    # Get relevant notebooks to execute
    notebook_paths = get_notebook_paths(category_to_get)
    # Execute and save notebooks
    execute_notebooks(notebook_paths, args.verbose)

    return 0


if __name__ == "__main__":
    # Parse arguments
    parser = argparse.ArgumentParser(
        description="Execute Jupyter Notebooks, rebuilding models in process")
    parser.add_argument(
        "-c", "--category",
        help="The category of notebooks to execute. Default is 'test'",
        type=str, choices=list(iterkeys(CATEGORIES_DICT)), default="test")
    parser.add_argument("-v", "--verbose", help="increase output verbosity",
                        action="store_true")
    # Run main script
    sys.exit(main(parser.parse_args()))
