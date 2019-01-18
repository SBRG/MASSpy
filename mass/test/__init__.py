# -*- coding: utf-8 -*-

# Compatibility with Python 2.7
from __future__ import absolute_import

# Import necesary packages
import re
from warnings import warn
from os import listdir
from os.path import abspath, dirname, join

from mass.io.json import load_json_model

from mass.exceptions import MassSBMLError as _MassSBMLError

_filetypes = "xml|json|txt"

try:
    from mass.io.sbml import read_sbml_model
except _MassSBMLError:
    _filetypes = "|".join(_filetypes.split("|")[1:])
try:
    import pytest
    import pytest_benchmark
except ImportError:
    pytest = None

file_re = re.compile("(\s+|\S+)\.(%s*)$" % _filetypes)
json_re = re.compile("json")
xml_re = re.compile("xml")
txt_re = re.compile("txt")
mass_dir = abspath(join(dirname(abspath(__file__)), ".."))
mass_loc = abspath(join(mass_dir, ".."))
data_dir = join(mass_dir, "test", "data", "")

def create_test_model(model_name):
    """Returns a mass.MassModel for testing

    Parameters
    ----------
    model_name : str
        The path to a mass.MassModel, or one of the following prebuilt models"
        {"glycolysis", "ppp", "ampsn", "hemoglobin"}

    Note: The priority will always be to try to load SBML before json
        unless otherwise specified or if SBML dependencies not installed.

    Return
    ------
    model : mass.MassModel
        The desired mass.MassModel
    """
    # Gather list of available models
    model_set = _make_model_set()
    # Check inputs
    if not isinstance(model_name, str):
        raise TypeError("model_name must be a string of either the filepath"
                        "or one of the following models: %s" % model_set)
    # Check to see if model_name is a filepath or just a name for a model
    if (file_re.match(model_name) and not re.search("[/]", model_name)) or \
        model_name in model_set:
        # Check filetype to determine loader, use specified type if provided
        if file_re.match(model_name):
            filetypes = [file_re.match(model_name).group(2)]
        # Otherwise use existing models
        else:
            filetypes = _filetypes.split("|")
        for ext in filetypes:
            filepath = model_name
            if json_re.match(ext):
                model_folder = "json_models"
            if xml_re.match(ext):
                model_folder = "sbml_models"
            if txt_re.match(ext):
                model_folder = "txt_models"
            if not file_re.match(filepath):
                filepath = "%s.%s" %(filepath, ext)
            if filepath in listdir(join(data_dir, "models" ,model_folder, "")):
                filepath = join(data_dir, "models", model_folder, filepath)
            else:
                continue
            # Load SBML model if dependencies installed
            if xml_re.match(ext) and re.search(ext, _filetypes):
                return read_sbml_model(filepath)
            # Raise error on attempt to load SBML without required dependencies
            elif xml_re.match(ext) and not re.search(ext, _filetypes):
                raise _MassSBMLError("Cannot import SBML models without "
                                    "SBML dependencies installed")
            # Load json model
            elif json_re.match(ext) and re.search(ext, _filetypes):
                return load_json_model(filepath)

        warn("Wrong extension specified. Cannot load model")
        return None

    elif file_re.match(model_name):
        ext = file_re.match(model_name).group(2)
        if xml_re.match(ext) and re.search(ext, _filetypes):
            return read_sbml_model(model_name)
        # Raise error on attempt to load SBML without required dependencies
        elif xml_re.match(ext) and not re.search(ext, _filetypes):
            raise _MassSBMLError("Cannot import SBML models without "
                                "SBML dependencies installed")
        # Load json model
        elif json_re.match(ext) and re.search(ext, _filetypes):
            return load_json_model(model_name)
    else:
        raise ValueError("Could not recognize file type")

def _make_model_set():
    model_set = set()
    model_dir = join(data_dir, "models", "")
    for model_folder in listdir(model_dir):
        if re.search("models", model_folder):
            for model in listdir(join(model_dir, model_folder, "")):
                if file_re.match(model):
                    model = file_re.match(model).group(1)
                    model_set.add(model)

    return model_set

def view_available_models():
    for model in list(sorted(_make_model_set())):
        print(model)

def test_all(args=None):
    """ alias for running all unit-tests on installed mass
    """
    if pytest:
        args = args if args else []

        return pytest.main(
            ['--pyargs', 'mass', '--benchmark-skip', '-v', '-rs'] + args
        )
    else:
        raise ImportError('missing package pytest and pytest_benchmark'
                          ' required for testing')