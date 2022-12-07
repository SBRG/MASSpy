# -*- coding: utf-8 -*-
"""Create or load models in `JSON <https://www.json.org/>`_ format.

Models created in JSON format can also be viewed using the
`Escher <https://escher.github.io/>`_ network visualization tool. See the
Escher `web-based tool <https://escher.readthedocs.io/en/stable/index.html>`_
or the
`Python package <https://escher.readthedocs.io/en/stable/escher-python.html>`_
documentation for more information on Escher.

To enable faster JSON export, the :mod:`simplejson` package can be installed
during the :mod:`mass` installation process as follows::

    # Installs simplejson package.
    pip install masspy["json"]
    # Or to install all additional packages.
    pip install masspy["all"]

If the :mod:`~mass.visualization` submodule is installed, see the
:mod:`mass.visualiation.escher` documentation for more information on using
:mod:`mass` with :mod:`escher` (Coming soon).
"""
try:
    import simplejson as json
except ImportError:
    import json

from six import string_types

from mass.io.dict import model_from_dict, model_to_dict


JSON_SPEC = "1"


def to_json(mass_model, sort=False, **kwargs):
    """Return the model as a JSON document.

    ``kwargs`` are passed on to ``json.dumps``

    Parameters
    ----------
    mass_model : MassModel or EnzymeModule
        The :mod:`mass` model to represent.
    sort : bool
        Whether to sort the objects in the lists representing attributes, or
        to maintain the order defined in the model. Default is ``False``.

    Returns
    -------
    str
        String representation of the :mod:`mass` model as a JSON document.

    See Also
    --------
    save_json_model
        Write directly to a file.
    json.dumps
        Base function.

    """
    obj = model_to_dict(mass_model, sort=sort)
    obj["version"] = JSON_SPEC
    return json.dumps(obj, allow_nan=False, **kwargs)


def from_json(document):
    """Load a model from a JSON document.

    Parameters
    ----------
    document : str
        The JSON document representation of a :mod:`mass` model.

    Returns
    -------
    MassModel or EnzymeModule
        The :mod:`mass` model as represented in the JSON document.

    See Also
    --------
    load_json_model
        Load directly from a file.

    """
    return model_from_dict(json.loads(document))


def save_json_model(mass_model, filename, sort=False, pretty=False, **kwargs):
    """Write the model to a file in JSON format.

    ``kwargs`` are passed on to ``json.dump``

    Parameters
    ----------
    mass_model : MassModel or EnzymeModule
        The :mod:`mass` model to represent.
    filename : str or file-like
        File path or descriptor the the JSON representation should be
        written to.
    sort : bool
        Whether to sort the objects in the lists representing attributes, or
        to maintain the order defined in the model. Default is ``False``.
    pretty : bool
        Whether to format the JSON more compactly (default), or in a more
        verbose but easier to read fashion. Default is ``False``.
        Can be partially overwritten by the ``kwargs``.

    See Also
    --------
    to_json
        Create a string represenation of the model in JSON format.
    json.dump
        Base function.

    """
    obj = model_to_dict(mass_model, sort=sort)
    obj["version"] = JSON_SPEC

    if pretty:
        dump_opts = {
            "indent": 4,
            "separators": (",", ": "),
            "sort_keys": True,
            "allow_nan": False,
        }
    else:
        dump_opts = {
            "indent": 4,
            "separators": (",", ": "),
            "sort_keys": True,
            "allow_nan": False,
        }
    dump_opts.update(**kwargs)

    if isinstance(filename, string_types):
        with open(filename, "w") as file_handle:
            json.dump(obj, file_handle, **dump_opts)
    else:
        json.dump(obj, filename, **dump_opts)


def load_json_model(filename):
    """Load the model from a file in JSON format.

    Parameters
    ----------
    filename : str or file-like
        File path or descriptor the contains JSON document describing the
        :mod:`mass` model to be loaded.

    Returns
    -------
    MassModel or EnzymeModule
        The :mod:`mass` model as represented in the JSON formatted file.

    See Also
    --------
    from_json
        Load a model from a string representation in JSON format.

    """
    if isinstance(filename, string_types):
        with open(filename, "r") as file_handle:
            return model_from_dict(json.load(file_handle))
    else:
        return model_from_dict(json.load(filename))


JSON_SCHEMA = {
    "$schema": "http://json-schema.org/draft-07/schema#",
    "title": "MASS",
    "description": "JSON representation of MASS model",
    "type": "object",
    "properties": {
        "id": {"type": "string"},
        "name": {"type": "string"},
        "version": {
            "type": "integer",
            "default": 1,
        },
        "reactions": {
            "type": "array",
            "items": {
                "type": "object",
                "properties": {
                    "id": {"type": "string"},
                    "name": {"type": "string"},
                    "reversible": {"type": "boolean"},
                    "metabolites": {
                        "type": "object",
                        "patternProperties": {
                            ".*": {"type": "number"},
                        },
                    },
                    "gene_reaction_rule": {"type": "string"},
                    "lower_bound": {"type": "number"},
                    "upper_bound": {"type": "number"},
                    "subsystem": {"type": "string"},
                    "steady_state_flux": {"type": "number"},
                    "forward_rate_constant": {"type": "number"},
                    "reverse_rate_constant": {"type": "number"},
                    "equilibriun_constant": {"type": "number"},
                    "objective_coefficient": {
                        "type": "number",
                        "default": 0,
                    },
                    "variable_kind": {
                        "type": "string",
                        "pattern": "integer|continuous",
                        "default": "continuous",
                    },
                    "_rate": {
                        "type": "string",
                    },
                    "notes": {"type": "object"},
                    "annotation": {"type": "object"},
                    "enzyme_module_id": {"type": "string"},
                },
            },
            "required": [
                "id",
                "name",
                "reversible",
                "metabolites",
                "lower_bound",
                "upper_bound",
                "gene_reaction_rule",
            ],
            "additionalProperties": False,
        },
        "metabolites": {
            "type": "array",
            "items": {
                "type": "object",
                "properties": {
                    "id": {"type": "string"},
                    "name": {"type": "string"},
                    "formula": {"type": "string"},
                    "charge": {"type": "integer"},
                    "compartment": {
                        "type": "string",
                        "pattern": "[a-z]{1,2}",
                    },
                    "fixed": {"type": "boolean"},
                    "_initial_condition": {
                        "type": "number",
                        "minimum": 0,
                        "exclusiveMinimum": False,
                    },
                    "_constraint_sense": {
                        "type": "string",
                        "default": "E",
                        "pattern": "E|L|G",
                    },
                    "_bound": {
                        "type": "number",
                        "default": 0,
                    },
                    "notes": {"type": "object"},
                    "annotation": {"type": "object"},
                    "_bound_metabolites": {
                        "type": "object",
                        "patternProperties": {
                            ".*": {"type": "number"},
                        },
                    },
                    "enzyme_module_id": {"type": "string"},
                },
                "required": ["id", "name"],
                "additionalProperties": False,
            },
        },
        "genes": {
            "type": "array",
            "items": {
                "type": "object",
                "properties": {
                    "id": {"type": "string"},
                    "name": {"type": "string"},
                    "notes": {"type": "object"},
                    "annotation": {"type": "object"},
                },
                "required": ["id", "name"],
                "additionalProperties": False,
            },
        },
        "enzyme_modules": {
            "type": "array",
            "items": {
                "type": "object",
                "properties": {
                    "id": {"type": "string"},
                    "name": {"type": "string"},
                    "subsystem": {"type": "string"},
                    "enzyme_module_ligands": {
                        "type": "array",
                        "allOf": {"type": "string"},
                    },
                    "enzyme_module_forms": {
                        "type": "array",
                        "allOf": {"type": "string"},
                    },
                    "enzyme_module_reactions": {
                        "type": "array",
                        "allOf": {"type": "string"},
                    },
                    "enzyme_module_ligands_categorized": {
                        "type": "object",
                        "allOf": {"type": "array", "allOf": {"type": "string"}},
                    },
                    "enzyme_module_forms_categorized": {
                        "type": "object",
                        "allOf": {"type": "array", "allOf": {"type": "string"}},
                    },
                    "enzyme_module_reactions_categorized": {
                        "type": "object",
                        "allOf": {"type": "array", "allOf": {"type": "string"}},
                    },
                    "enzyme_concentration_total": {
                        "type": "number",
                        "minimum": 0,
                        "exclusiveMinimum": False,
                    },
                    "enzyme_rate": {"type": "number"},
                    "enzyme_concentration_total_equation": {"type": "string"},
                    "enzyme_rate_equation": {"type": "string"},
                },
            },
            "required": ["id", "name"],
            "additionalProperties": False,
        },
        "units": {
            "type": "array",
            "items": {
                "type": "object",
                "properties": {
                    "kind": {"type": "string"},
                    "exponent": {"type": "number"},
                    "scale": {"type": "number"},
                    "multiplier": {"type": "number"},
                },
            },
        },
        "boundary_conditions": {
            "type": "object",
            "allOf": {
                "type": "number",
                "minimum": 0,
            },
        },
        "custom_rates": {
            "type": "object",
            "patternProperties": {
                ".*": {"type": "string"},
            },
        },
        "custom_parameters": {
            "type": "object",
            "patternProperties": {
                ".*": {"type": "number"},
            },
        },
        "compartments": {
            "type": "object",
            "patternProperties": {
                "[a-z]{1,2}": {"type": "string"},
            },
        },
        "notes": {"type": "object"},
        "annotation": {"type": "object"},
        "enzyme_module_ligands": {"type": "array", "allOf": {"type": "string"}},
        "enzyme_module_forms": {"type": "array", "allOf": {"type": "string"}},
        "enzyme_module_reactions": {"type": "array", "allOf": {"type": "string"}},
        "_enzyme_module_ligands_categorized": {
            "type": "object",
            "allOf": {"type": "array", "allOf": {"type": "string"}},
        },
        "_enzyme_module_forms_categorized": {
            "type": "object",
            "allOf": {"type": "array", "allOf": {"type": "string"}},
        },
        "_enzyme_module_reactions_categorized": {
            "type": "object",
            "allOf": {"type": "array", "allOf": {"type": "string"}},
        },
        "_enzyme_concentration_total": {
            "type": "number",
            "minimum": 0,
            "exclusiveMinimum": False,
        },
        "_enzyme_rate": {"type": "number"},
        "_enzyme_rate_equation": {"type": "string"},
    },
    "required": ["id", "reactions", "metabolites", "genes"],
    "additionalProperties": False,
}
"""dict: The generic JSON schema for representing a model in :mod:`mass`."""

__all__ = ("to_json", "from_json", "save_json_model", "load_json_model", "JSON_SCHEMA")
