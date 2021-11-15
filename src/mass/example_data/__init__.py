# -*- coding: utf-8 -*-
"""
Module containing functions for loading various examples.

"""
from os import listdir
from os.path import abspath, dirname, join

from mass.io.json import load_json_model
from mass.io.sbml import read_sbml_model


FILE_EXTENSIONS = [".xml", ".json"]
"""list: list of recognized file extensions."""

MODELS_DIR = join(abspath(dirname(abspath(__file__))), "models", "")
"""str: The directory location of the pre-built :class:`MassModel` files."""

MAPS_DIR = join(abspath(dirname(abspath(__file__))), "maps", "")
"""str: The directory location of the pre-made :mod:`escher` maps files."""


def create_example_model(model_name, io="json"):
    """Return an example :class:`mass.MassModel <~.MassModel>`.

    Parameters
    ----------
    model_name: str
        The name of the example model to load. Valid model names can be printed
        and viewed using the :func:`view_example_models` function.
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


def view_example_models():
    """Print the example models that can be loaded."""
    return _get_directory_files(MODELS_DIR)


def view_example_maps():
    """Print the example models that can be loaded."""
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
    "MODELS_DIR",
    "MAPS_DIR",
    "create_example_model",
    "view_example_models",
    "view_example_maps",
)
