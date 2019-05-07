# -*- coding: utf-8 -*-
"""TODO Module Docstrings."""
from __future__ import absolute_import

try:
    import simplejson as json
except ImportError:
    import json

from six import string_types

from mass.io.dict import model_from_dict, model_to_dict

JSON_SPEC = "1"


def to_json(model, sort=False, **kwargs):
    """Return the model as a JSON document.

    `kwargs`` are passed on to ``json.dumps``
    See documentation for ``json.dumps`` for more details.

    Parameters
    ----------
    model: MassModel
        The MassModel object to represent.
    sort: bool, optional
        Whether to sort the metabolites, reactions, and genes or maintain the
        order defined in the model. Default is false.

    Returns
    -------
    str
        String representation of the MassModel as a JSON documentation

    See Also
    --------
    save_json_model: Write directly to a file.
    mass.io.dict.model_to_dict: Make dict representation of a MassModel.
    json.dumps: Base Function.

    """
    obj = model_to_dict(model, sort=sort)
    obj[u'version'] = JSON_SPEC
    return json.dumps(obj, allow_nan=False, **kwargs)


def from_json(document):
    """Load a MassModel as a JSON document.

    `kwargs`` are passed on to ``json.dumps``
    See documentation for ``json.dumps`` for more details.

    Parameters
    ----------
    document: str
        The JSON document representation of a MassModel.

    Returns
    -------
    MassModel:
        The MassModel as represented in the JSON documentation.

    See Also
    --------
    load_json_model: Load MassModel directly from a file.
    mass.io.dict.model_to_dict: Make MassModel from dict representation.

    """
    return model_from_dict(json.loads(document))


def save_json_model(model, filename, sort=False, pretty=False,
                    extension=True, **kwargs):
    """Write the MassModel to a file in JSON format.

    `kwargs`` are passed on to ``json.dumps``
    See documentation for ``json.dumps`` for more details.

    Parameters
    ----------
    model: MassModel
        The MassModel object to represent.
    filename: str or file-like
        File path or descriptor the the JSON representation should be
        written to.
    sort: bool, optional
        Whether to sort the metabolites, reactions, and genes or maintain the
        order defined in the model. Default is false.
    pretty: bool, optional
        Whether to format the JSON more compactly (default), or in a more
        verbose but easier to read fashion. Default is False. Can be partially
        overwritten by the ``kwargs``.
    extension: bool, optional
        Whether to include the file type extension (*.json) at the end of the
        filename if a string is provided. Default is True.

    See Also
    --------
    to_json: Return a string represenation.
    json.dump: Base function.

    """
    obj = model_to_dict(model, sort=sort)
    obj[u'version'] = JSON_SPEC

    if pretty:
        dump_opts = {"indent": 4, "separators": (",", ": "),
                     "sort_keys": True, "allow_nan": False}
    else:
        dump_opts = {"indent": 4, "separators": (",", ": "),
                     "sort_keys": True, "allow_nan": False}
    dump_opts.update(**kwargs)

    if isinstance(filename, string_types):
        if extension and ".json" not in filename[-5:]:
            filename += ".json"
        with open(filename, "w") as file_handle:
            json.dump(obj, file_handle, **dump_opts)
    else:
        json.dump(obj, filename, **dump_opts)


def load_json_model(filename):
    """Load the MassModel from a file in JSON format.

    Parameters
    ----------
    filename: str or file-like
        File path or descriptor the contains JSON document describing the
        MassModel to be loaded.

    Returns
    -------
    MassModel:
        The MassModel as represented in the JSON documentation.

    See Also
    --------
    from_json: Load from a string representation.
    mass.io.dict.model_to_dict: Make MassModel from dict representation.

    """
    if isinstance(filename, string_types):
        with open(filename, "r") as file_handle:
            return model_from_dict(json.load(file_handle))
    else:
        return model_from_dict(json.load(filename))


json_schema = {
    "$schema": "http://json-schema.org/draft-07/schema#",
    "title": "MASS",
    "description": "JSON representation of MASS model",
    "type": "object",
    "properties": {
        "id": {"type": "string"},
        "name": {"type": "string"},
        "description": {"type": "string"},
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
                        "default": "continuous"
                    },
                    "_rtype": {
                        "type": "integer",
                        "enum": [1, 2, 3]
                    },
                    "notes": {"type": "object"},
                    "annotation": {"type": "object"},
                    "enzyme_id": {"type": "string"},
                },
            },
            "required": ["id", "name", "reversible", "metabolites",
                         "lower_bound", "upper_bound", "gene_reaction_rule"],
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
                    "_bound_catalytic": {
                        "type": "object",
                        "patternProperties": {
                            ".*": {"type": "number"},
                        },
                    },
                    "_bound_effectors": {
                        "type": "object",
                        "patternProperties": {
                            ".*": {"type": "number"},
                        },
                    },
                    "enzyme_id": {"type": "string"},

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
                        "allOf": {"type": "string"}
                    },
                    " enzyme_module_forms": {
                        "type": "array",
                        "allOf": {"type": "string"}
                    },
                    "enzyme_module_reactions": {
                        "type": "array",
                        "allOf": {"type": "string"}
                    },
                    "enzyme_module_ligands_categorized": {
                        "type": "object",
                        "allOf": {
                            "type": "array",
                            "allOf": {"type": "string"}
                        },
                    },
                    "enzyme_module_form_categorizeds": {
                        "type": "object",
                        "allOf": {
                            "type": "array",
                            "allOf": {"type": "string"}
                        },
                    },
                    "enzyme_module_reactions_categorized": {
                        "type": "object",
                        "allOf": {
                            "type": "array",
                            "allOf": {"type": "string"}
                        },
                    },
                    "enzyme_concentration_total": {
                        "type": "number",
                        "minimum": 0,
                        "exclusiveMinimum": False,
                    },
                    "enzyme_net_flux": {"type": "number"},
                    "enzyme_net_flux_equation": {"type": "string"},
                },
            },
            "required": ["id", "name"],
            "additionalProperties": False,
        },
        "initial_conditions": {
            "type": "object",
            "allOf": {
                "type": "number",
                "minimum": 0,
            },
        },
        "fixed_concentrations": {
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
        "units": {
            "type": "object",
            "patternProperties": {
                ".*": {"type": "string"},
            },
        },
        "notes": {"type": "object"},
        "annotation": {"type": "object"},
        "enzyme_module_ligands": {
            "type": "array",
            "allOf": {"type": "string"}
        },
        " enzyme_module_forms": {
            "type": "array",
            "allOf": {"type": "string"}
        },
        "enzyme_module_reactions": {
            "type": "array",
            "allOf": {"type": "string"}
        },
        "_enzyme_module_ligands_categorized": {
            "type": "object",
            "allOf": {
                "type": "array",
                "allOf": {"type": "string"}
            },
        },
        "_enzyme_module_form_categorizeds": {
            "type": "object",
            "allOf": {
                "type": "array",
                "allOf": {"type": "string"}
            },
        },
        "_enzyme_module_reactions_categorized": {
            "type": "object",
            "allOf": {
                "type": "array",
                "allOf": {"type": "string"}
            },
        },
        "_enzyme_concentration_total": {
            "type": "number",
            "minimum": 0,
            "exclusiveMinimum": False,
        },
        "_enzyme_net_flux": {"type": "number"},
        "_enzyme_net_flux_equation": {"type": "string"},
    },
    "required": ["id", "reactions", "metabolites", "genes"],
    "additionalProperties": False,
}
