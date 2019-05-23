# -*- coding: utf-8 -*-
"""TODO Module Docstrings."""
from __future__ import absolute_import

import argparse
import shutil
import sys
from os import listdir, path


SCRIPTS_DIR = path.dirname(path.abspath(__file__))
NOTEBOOK_MODELS_DIR = path.abspath(path.join(
    SCRIPTS_DIR, "..", "notebooks", "test-models", "models"))
MASS_TEST_MODELS_DIR = path.abspath(path.join(
    SCRIPTS_DIR, "..", "mass", "test", "data", "models", ""))


def _get_model_names():
    """Get the available models in the notebook model directory."""
    return sorted(list(set([
        m.split(".")[0] for m in listdir(NOTEBOOK_MODELS_DIR)])))


def main(args):
    """Copy models created by the notebooks into the masspy testing data."""
    # Get list of models to copy
    if args.name == "all":
        model_files_to_copy = _get_model_names()
    else:
        model_files_to_copy = [args.name]

    # Iterate through directory, copying models from notebook
    for model_file in sorted(listdir(NOTEBOOK_MODELS_DIR)):
        if model_file.split(".")[0] not in model_files_to_copy:
            continue
        if args.verbose:
            print("Copying " + model_file)
        file_to_copy = path.abspath(path.join(NOTEBOOK_MODELS_DIR, model_file))
        shutil.copy(file_to_copy, MASS_TEST_MODELS_DIR)


if __name__ == "__main__":
    # Parse arguments
    parser = argparse.ArgumentParser(
        description="Replace models in /mass/test/data/models directory with "
                    "models stored in /mass/notebooks/test-models/models.")
    parser.add_argument(
        "-n", "--name",
        help="The name of models to copy and move. Default is 'all' for all "
             "models found in the directory.",
        type=str, choices=list(_get_model_names() + ["all"]), default="all")
    parser.add_argument("-v", "--verbose", help="increase output verbosity",
                        action="store_true")
    # Run main script
    sys.exit(main(parser.parse_args()))
