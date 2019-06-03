# -*- coding: utf-8 -*-
"""TODO Module Docstrings."""
# Compatibility with Python 2.7
from __future__ import absolute_import

from os import listdir
from os.path import abspath, dirname, join

try:
    import pytest
    import pytest_benchmark
except ImportError:
    pytest = None

from mass.io.json import load_json_model
from mass.io.sbml import read_sbml_model

FILETYPES = [".xml", ".json"]
MASS_DIR = abspath(join(dirname(abspath(__file__)), ".."))
MASS_LOC = abspath(join(MASS_DIR, ".."))
DATA_DIR = join(MASS_DIR, "test", "data", "")
MODELS_DIR = join(DATA_DIR, "models", "")
MAPS_DIR = join(DATA_DIR, "maps", "")

def create_test_model(model_name, io="sbml"):
    """Return a mass MassModel for testing.

    Parameters
    ----------
    model_name: str
        The name of the test model to load. Valid model names can be printed
        and viewed using the `view_test_models` method.
    io: {'sbml', 'json'}, optional
        A string representing the mass.io module to use to load the model.
        Default is "sbml". Case sensitive.

    Returns
    -------
    MassModel:
        The loaded MassModel object.

    """
    load_dict = {"sbml": (read_sbml_model, ".xml"),
                 "json": (load_json_model, ".json")}

    try:
        load_function, ext = load_dict[io]
    except KeyError as e:
        raise ValueError(
            "Unrecognized value {0} for io. Value must be a str of one of "
            "the following {1}".format(e, str(set(load_dict))))
    if not model_name.endswith(ext):
        model_name = model_name + ext
    filepath = join(MODELS_DIR, model_name)

    return load_function(filepath)


def view_test_models():
    """Print the test models that can be loaded."""
    for model_name in _get_directory_files(MODELS_DIR):
        print(model_name)


def view_test_maps():
    """Print the test models that can be loaded."""
    for map_name in _get_directory_files(MAPS_DIR):
        print(map_name)


def _get_directory_files(directory):
    """Return a list of files in a given directory."""
    all_filenames = []
    for filename in listdir(directory):
        # Do not include invalid/broken models used in testing suite.
        if "invalid" in filename:
            continue
        if any(list(map(lambda x: filename.endswith(x), FILETYPES))):
            all_filenames.append(filename)

    return sorted(all_filenames)

def test_all(args=None):
    """Alias for running all unit-tests on installed mass."""
    if pytest is None:
        raise ImportError('missing package pytest and pytest_benchmark'
                          ' required for testing')

    args = args if args else []
    return pytest.main(
        ['--pyargs', 'mass', '--benchmark-skip', '-v', '-rs'] + args)
