# -*- coding: utf-8 -*-
"""
Module containing functions for testing various :mod:`mass` methods.

There is the :func:`test_all` function to run all tests. Note that the testing
requirements must be installed (e.g. :mod:`pytest`) for this function to work.

There are also functions for viewing and loading pre-built example
:class:`~.MassModel` objects via the :mod:`~mass.io.json` or
:mod:`~mass.io.sbml` submodules, as well as functions to view
`Escher <https://escher.readthedocs.io/en/stable/escher-python.html>`_
maps that correspond to certain pre-defined models.
"""
from os import listdir
from os.path import abspath, dirname, join

from mass.io.json import load_json_model
from mass.io.sbml import read_sbml_model


FILE_EXTENSIONS = [".xml", ".json"]
"""list: list of recognized file extensions."""

DATA_DIR = join(abspath(join(dirname(abspath(__file__)), "data", "")))
"""str: The directory location of the test data model files and maps."""

MODELS_DIR = join(DATA_DIR, "models", "")
"""str: The directory location of the pre-built :class:`MassModel` files."""

MAPS_DIR = join(DATA_DIR, "maps", "")
"""str: The directory location of the pre-made :mod:`escher` maps files."""


def create_test_model(model_name, io="json"):
    """Return a :class:`mass.MassModel <~.MassModel>` for testing.

    Parameters
    ----------
    model_name: str
        The name of the test model to load. Valid model names can be printed
        and viewed using the :func:`view_test_models` function.
    io: str ``{'sbml', 'json'}``
        A string representing the :mod:`mass.io` module to use to load the
        model. Default is ``"sbml"``. Case sensitive.

    Returns
    -------
    MassModel
        The loaded :class:`~.MassModel`

    """
    load_dict = {"sbml": (read_sbml_model, ".xml"), "json": (load_json_model, ".json")}

    try:
        load_function, ext = load_dict[io]
    except KeyError as e:
        raise ValueError(
            "Unrecognized value {0} for io. Value must be a str of one of "
            "the following {1}".format(e, str(set(load_dict)))
        )
    if not model_name.endswith(ext):
        model_name = model_name + ext
    filepath = join(MODELS_DIR, model_name)

    return load_function(filepath)


def view_test_models():
    """Print the test models that can be loaded."""
    return _get_directory_files(MODELS_DIR)


def view_test_maps():
    """Print the test models that can be loaded."""
    return _get_directory_files(MAPS_DIR)


def _get_directory_files(directory):
    """Return a list of files in a given directory."""
    all_filenames = []
    for filename in listdir(directory):
        # Do not include invalid/broken models used in testing suite.
        if "invalid" in filename:
            continue
        if any(list(map(lambda x: filename.endswith(x), FILE_EXTENSIONS))):
            all_filenames.append(filename)

    return sorted(all_filenames)


__all__ = (
    "FILE_EXTENSIONS",
    "DATA_DIR",
    "MODELS_DIR",
    "MAPS_DIR",
    "create_test_model",
    "view_test_models",
    "view_test_maps",
)
