# -*- coding: utf-8 -*-
r"""Module to convert or create :mod:`mass` objects into or from dictionaries.

Converting objects into dictionaries allow for the exportation of
:class:`~.MassModel`\ s in various formats. These formats include:

    * `JSON <https://www.json.org/>`_ format using the functions in
      :mod:`~mass.io.json`.

"""
from collections import OrderedDict
from operator import attrgetter, itemgetter

import numpy as np
import pandas as pd
from cobra.io.dict import gene_from_dict, gene_to_dict
from cobra.util.solver import set_objective
from six import iteritems, iterkeys, string_types
from sympy import Basic, sympify

from mass.core.mass_metabolite import MassMetabolite
from mass.core.mass_model import MassModel
from mass.core.mass_reaction import MassReaction
from mass.core.units import UnitDefinition
from mass.enzyme_modules.enzyme_module import EnzymeModule
from mass.enzyme_modules.enzyme_module_dict import (
    _ORDERED_ENZYMEMODULE_DICT_DEFAULTS,
    EnzymeModuleDict,
)
from mass.enzyme_modules.enzyme_module_form import EnzymeModuleForm
from mass.enzyme_modules.enzyme_module_reaction import EnzymeModuleReaction


# Global
_INF = float("inf")

_REQUIRED_REACTION_ATTRIBUTES = [
    "id",
    "name",
    "_reversible",
    "metabolites",
    "_lower_bound",
    "_upper_bound",
    "gene_reaction_rule",
]
_ORDERED_OPTIONAL_REACTION_KEYS = [
    "subsystem",
    "steady_state_flux",
    "_forward_rate_constant",
    "_reverse_rate_constant",
    "_equilibrium_constant",
    "objective_coefficient",
    "_rate_type",
    "notes",
    "annotation",
]
_OPTIONAL_REACTION_ATTRIBUTES = {
    "subsystem": "",
    "steady_state_flux": None,
    "_forward_rate_constant": None,
    "_equilibrium_constant": None,
    "_reverse_rate_constant": None,
    "objective_coefficient": 0,
    "_rate_type": 1,
    "notes": {},
    "annotation": {},
}
_REQUIRED_ENZYMEMODULEREACTION_ATTRIBUTES = ["enzyme_module_id"]
_ORDERED_OPTIONAL_ENZYMEMODULEREACTION_KEYS = []
_OPTIONAL_ENZYMEMODULEREACTION_ATTRIBUTES = {}

_REQUIRED_METABOLITE_ATTRIBUTES = ["id", "name"]
_ORDERED_OPTIONAL_METABOLITE_KEYS = [
    "formula",
    "charge",
    "compartment",
    "fixed",
    "_initial_condition",
    "_bound",
    "notes",
    "annotation",
]
_OPTIONAL_METABOLITE_ATTRIBUTES = {
    "charge": None,
    "formula": None,
    "compartment": None,
    "fixed": False,
    "_initial_condition": None,
    "_bound": 0,
    "notes": {},
    "annotation": {},
}

_REQUIRED_ENZYMEMODULEFORM_ATTRIBUTES = ["_bound_metabolites", "enzyme_module_id"]
_ORDERED_OPTIONAL_ENZYMEMODULEFORM_KEYS = []
_OPTIONAL_ENZYMEMODULEFORM_ATTRIBUTES = {}

_REQUIRED_ENZYMEMODULE_ATTRIBUTES = [
    "id",
    "name",
    "enzyme_module_ligands",
    "enzyme_module_forms",
    "enzyme_module_reactions",
]
_ORDERED_OPTIONAL_ENZYMEMODULE_KEYS = [
    key
    for key in iterkeys(_ORDERED_ENZYMEMODULE_DICT_DEFAULTS)
    if key not in _REQUIRED_ENZYMEMODULE_ATTRIBUTES + ["S", "model"]
]
_OPTIONAL_ENZYMEMODULE_ATTRIBUTES = OrderedDict(
    {
        key: _ORDERED_ENZYMEMODULE_DICT_DEFAULTS[key]
        for key in _ORDERED_OPTIONAL_ENZYMEMODULE_KEYS
    }
)

_ORDERED_OPTIONAL_MODEL_KEYS = ["name", "compartments", "notes", "annotation"]
_OPTIONAL_MODEL_ATTRIBUTES = {
    "name": None,
    "compartments": {},
    "notes": {},
    "annotation": {},
}


def model_to_dict(model, sort=False):
    """Convert a :class:`~.MassModel` into a serializable dictionary.

    Parameters
    ----------
    model : MassModel or EnzymeModule
        The model to represent as a dictionary.
    sort : bool
        Whether to sort the metabolites, reactions, genes, and enzyme modules
        or maintain the order defined in the model. If the model is an
        :class:`~.EnzymeModule`, the
        :attr:`~.EnzymeModule.enzyme_module_ligands`,
        :attr:`~.EnzymeModule.enzyme_module_forms`, and
        :attr:`~.EnzymeModule.enzyme_module_reactions` attributes are also
        included. Default is ``False``.

    Returns
    -------
    ~collections.OrderedDict
        A dictionary with elements corresponding to the model attributes as
        which are in turn lists containing dictionaries holding all attribute
        information to form the corresponding object.

    See Also
    --------
    model_from_dict

    """
    obj = OrderedDict()
    obj["id"] = model.id
    obj["metabolites"] = list(map(metabolite_to_dict, model.metabolites))
    obj["reactions"] = list(map(reaction_to_dict, model.reactions))
    obj["genes"] = list(map(gene_to_dict, model.genes))
    obj["enzyme_modules"] = list(map(enzyme_to_dict, model.enzyme_modules))
    obj["units"] = list(map(unit_to_dict, model.units))

    if sort:
        get_id = itemgetter("id")
        obj["metabolites"].sort(key=get_id)
        obj["reactions"].sort(key=get_id)
        obj["genes"].sort(key=get_id)
        obj["enzyme_modules"].sort(key=get_id)
        obj["units"].sort(key=get_id)

    for key in ["custom_rates", "custom_parameters", "boundary_conditions"]:
        values = getattr(model, key, {})
        if values:
            values = OrderedDict(
                (getattr(k, "_id", k), _fix_type(v)) for k, v in iteritems(values)
            )
        if sort:
            values = OrderedDict((k, values[k]) for k in sorted(values))
        obj[key] = values

    _update_optional(
        model, obj, _OPTIONAL_MODEL_ATTRIBUTES, _ORDERED_OPTIONAL_MODEL_KEYS
    )
    # Add EnzymeModule attributes if an EnzymeModule is being saved.
    if isinstance(model, EnzymeModule):
        _add_enzyme_module_attributes_into_dict(model, obj)
        if sort:
            obj["enzyme_module_ligands"].sort(key=get_id)
            obj["enzyme_module_forms"].sort(key=get_id)
            obj["enzyme_module_reactions"].sort(key=get_id)
    return obj


def model_from_dict(obj):
    """Create a :class:`~.MassModel` from a dictionary.

    Notes
    -----
    The :attr:`~.EnzymeModule.enzyme_module_ligands`,
    :attr:`~.EnzymeModule.enzyme_module_forms`, and
    :attr:`~.EnzymeModule.enzyme_module_reactions` attributes are used to
    determine whether the model should be initialized as an
    :class:`~.EnzymeModule` or as a :class:`~.MassModel`. At least one of these
    three attributes must be present in order for an :class:`~.EnzymeModule`
    to be created.

    Parameters
    ----------
    obj : dict
        A dictionary with elements corresponding to the model attributes as
        which are in turn lists containing dictionaries holding all attribute
        information to form the corresponding object.

    Returns
    -------
    MassModel or EnzymeModule
        The generated model or enzyme module.

    See Also
    --------
    model_to_dict

    """
    if "reactions" not in obj:
        raise ValueError("Object has no reactions attribute. Cannot load.")
    if any([k in obj for k in _REQUIRED_ENZYMEMODULE_ATTRIBUTES[2:]]):
        model = EnzymeModule(obj["id"])
    else:
        model = MassModel(obj["id"])
    # Add metabolites to the model
    model.add_metabolites(
        [metabolite_from_dict(metabolite) for metabolite in obj["metabolites"]]
    )

    # Add genes to the model
    model.genes.extend([gene_from_dict(gene) for gene in obj["genes"]])

    # Add reactions to the model
    model.add_reactions(
        [reaction_from_dict(reaction, model) for reaction in obj["reactions"]]
    )

    # Add objective coefficients to the model
    set_objective(
        model,
        {
            model.reactions.get_by_id(rxn["id"]): rxn["objective_coefficient"]
            for rxn in obj["reactions"]
            if rxn.get("objective_coefficient", 0) != 0
        },
    )

    # Add units to the model
    model.add_units([unit_from_dict(unit_def) for unit_def in obj["units"]])
    # Add enzyme modules to the model
    if "enzyme_modules" in obj:
        model.enzyme_modules.extend(
            [enzyme_from_dict(enz, model) for enz in obj["enzyme_modules"]]
        )
    # Repair model once all objects are in model.
    model.repair(rebuild_index=True, rebuild_relationships=True)
    # Add boundary conditions to the model if they exist
    if "boundary_conditions" in obj:
        model.add_boundary_conditions(
            {met: bc for met, bc in iteritems(obj["boundary_conditions"])}
        )

    # Get custom parameters if they exist
    if "custom_parameters" in obj:
        model.custom_parameters.update(
            dict(
                (k, float(v)) if v not in ["", None] else (k, None)
                for k, v in iteritems(obj["custom_parameters"])
            )
        )

    # Add custom rates and any custom parameters if they exist
    if "custom_rates" in obj:
        # Add custom rates to the model
        for reaction, custom_rate in iteritems(obj["custom_rates"]):
            model.add_custom_rate(
                model.reactions.get_by_id(reaction),
                custom_rate,
                model.custom_parameters,
            )
    # Update with any opitonal attributes.
    for k, v in iteritems(obj):
        # Set MassModel attributes (and subsystem attribute for EnzymeModules)
        if k in _ORDERED_OPTIONAL_MODEL_KEYS or k == "subsystem":
            setattr(model, k, v)
        elif k == "enzyme_concentration_total_equation":
            continue
        # Update with EnzymeModule attributes if obj represents an EnzymeModule
        elif k.lstrip("_") in _ORDERED_OPTIONAL_ENZYMEMODULE_KEYS:
            model.__class__.__dict__[k.lstrip("_")].fset(model, v)

    return model


def metabolite_to_dict(metabolite):
    """Convert a :class:`~.MassMetabolite` into a serializable dictionary.

    Parameters
    ----------
    metabolite : ~.MassMetabolite
        The metabolite to represent as a dictionary.

    Returns
    -------
    ~collections.OrderedDict
        A dictionary with elements corresponding to metabolite attributes.

    See Also
    --------
    metabolite_from_dict

    """
    # Turn object into an OrderedDict with the required attributes
    new_met = OrderedDict()
    for key in _REQUIRED_METABOLITE_ATTRIBUTES:
        new_met[key] = _fix_type(getattr(metabolite, key))
    # Update with any opitonal attributes that are not their defaults.
    _update_optional(
        metabolite,
        new_met,
        _OPTIONAL_METABOLITE_ATTRIBUTES,
        _ORDERED_OPTIONAL_METABOLITE_KEYS,
    )

    # Add EnzymeModuleForm attributes if metabolite is an
    # EnzymeModuleForm
    if isinstance(metabolite, EnzymeModuleForm):
        _add_enzyme_module_form_attributes_into_dict(metabolite, new_met)

    return new_met


def metabolite_from_dict(metabolite):
    """Create a :class:`~.MassMetabolite` from a dictionary.

    Notes
    -----
    The presence of the :attr:`~.EnzymeModuleForm.enzyme_module_id`
    attribute is used to determine whether the dictionary should be
    initialized as an :class:`~.EnzymeModuleForm` or as a
    :class:`~.MassMetabolite`.

    Parameters
    ----------
    metabolite : dict
        A dictionary with elements corresponding to the metabolite attributes.

    Returns
    -------
    MassMetabolite or EnzymeModuleForm
        The generated metabolite.

    See Also
    --------
    metabolite_to_dict

    """
    # Determine if saved object should be a MassMetabolite or a subclass
    if "enzyme_module_id" in metabolite:
        new_metabolite = EnzymeModuleForm(metabolite["id"])
    else:
        new_metabolite = MassMetabolite(metabolite["id"])

    # Set object attributes
    for k, v in iteritems(metabolite):
        setattr(new_metabolite, k, v)

    return new_metabolite


def reaction_to_dict(reaction):
    """Convert a :class:`~.MassReaction` into a serializable dictionary.

    Parameters
    ----------
    reaction : ~.MassReaction
        The reaction to represent as a dictionary.

    Returns
    -------
    ~collections.OrderedDict
        A dictionary with elements corresponding to reaction attributes.

    See Also
    --------
    reaction_from_dict

    """
    # Turn object into an OrderedDict with the required attributes
    new_reaction = OrderedDict()
    for key in _REQUIRED_REACTION_ATTRIBUTES:
        if key != "metabolites":
            new_reaction[key] = _fix_type(getattr(reaction, key))
            continue
        # Store metabolites objects as their string identifiers
        mets = OrderedDict()
        for met in sorted(reaction.metabolites, key=attrgetter("id")):
            mets[str(met)] = reaction.metabolites[met]
        new_reaction["metabolites"] = mets
    # Update with any opitonal attributes that are not their defaults.
    _update_optional(
        reaction,
        new_reaction,
        _OPTIONAL_REACTION_ATTRIBUTES,
        _ORDERED_OPTIONAL_REACTION_KEYS,
    )
    # Add EnzymeModuleReaction attributes
    # if reaction is an EnzymeModuleReaction
    if isinstance(reaction, EnzymeModuleReaction):
        _add_enzyme_module_reaction_attributes_into_dict(reaction, new_reaction)

    return new_reaction


def reaction_from_dict(reaction, model):
    """Create a :class:`~.MassReaction` from a dictionary.

    Notes
    -----
    The presence of the :attr:`.EnzymeModuleReaction.enzyme_module_id`
    attribute is used to determine whether the dictionary should be initialized
    as an :class:`~.EnzymeModuleReaction` or as a :class:`~.MassReaction`.

    Parameters
    ----------
    reaction : dict
        A dictionary with elements corresponding to the reaction attributes.
    model : MassModel
        The model to assoicate with the reaction.

    Returns
    -------
    MassReaction or EnzymeModuleReaction
        The generated reaction.

    See Also
    --------
    reaction_to_dict

    """
    # Determine if saved object should be a MassReaction or a subclass
    if "enzyme_module_id" in reaction:
        new_reaction = EnzymeModuleReaction(reaction["id"])
    else:
        new_reaction = MassReaction(reaction["id"])

    # Set object attributes
    for k, v in iteritems(reaction):
        # Change infinity type from a string to a float
        if isinstance(v, string_types) and v == "inf":
            v = _INF
        if k in {"objective_coefficient", "reaction"}:
            continue
        elif k == "metabolites":
            new_reaction.add_metabolites(
                OrderedDict(
                    (model.metabolites.get_by_id(str(met)), coeff)
                    for met, coeff in iteritems(v)
                )
            )
        else:
            setattr(new_reaction, k, v)

    return new_reaction


def enzyme_to_dict(enzyme):
    """Convert an :class:`~.EnzymeModuleDict` into a serializable dictionary.

    Parameters
    ----------
    enzyme : ~.EnzymeModuleDict
        The enzyme module to represent as a dictionary.

    Returns
    -------
    ~collections.OrderedDict
        A dictionary with elements corresponding to the enzyme module
        attributes.

    See Also
    --------
    enzyme_from_dict

    """
    # Turn object into an OrderedDict with the required attributes
    new_enzyme = OrderedDict()
    for key in _REQUIRED_ENZYMEMODULE_ATTRIBUTES:
        new_enzyme[key] = _fix_type(getattr(enzyme, key))

    # Update with any opitonal attributes that are not their defaults.
    _update_optional(
        enzyme,
        new_enzyme,
        _OPTIONAL_ENZYMEMODULE_ATTRIBUTES,
        _ORDERED_OPTIONAL_ENZYMEMODULE_KEYS,
    )

    # Store objects and expressions as string for the attributes
    for key in _REQUIRED_ENZYMEMODULE_ATTRIBUTES[2:]:
        new_enzyme[key] = [i.id for i in getattr(enzyme, key)]
        # Repeat for categorized attribute
        key += "_categorized"
        if getattr(enzyme, key) != _OPTIONAL_ENZYMEMODULE_ATTRIBUTES[key]:
            new_enzyme[key] = {
                g.id: [i.id for i in g.members] for g in getattr(enzyme, key)
            }

    for key in ["enzyme_concentration_total_equation", "enzyme_rate_equation"]:
        if key in new_enzyme:
            new_enzyme[key] = str(getattr(enzyme, key))

    return new_enzyme


def enzyme_from_dict(enzyme, model):
    """Create an :class:`~.EnzymeModuleDict` from a dictionary.

    Parameters
    ----------
    enzyme : dict
        A dictionary with elements corresponding to the enzyme module
        dictionary attributes.
    model : MassModel
        The model to assoicate with the enzyme module dictionary.

    Returns
    -------
    EnzymeModuleDict
        The generated enzyme module dictionary.

    See Also
    --------
    enzyme_to_dict

    """
    # Set object attributes
    new_enzyme = EnzymeModuleDict(id_or_enzyme=enzyme)
    # Update model and get objects from the model to populate the DictLists
    new_enzyme["model"] = model
    new_enzyme._update_object_pointers(model)

    # Make the enzyme equations
    new_enzyme.enzyme_concentration_total_equation = sympify(
        new_enzyme.enzyme_concentration_total_equation
    )
    new_enzyme.enzyme_rate_equation = sympify(new_enzyme.enzyme_rate_equation)

    # Make the stoichiometric matrix and clean up the EnzymeModuleDict
    new_enzyme._make_enzyme_stoichiometric_matrix(update=True)
    new_enzyme._set_missing_to_defaults()
    new_enzyme._fix_order()

    return new_enzyme


def unit_to_dict(unit_definition):
    """Convert an :class:`~.UnitDefintion` into a serializable dictionary.

    Parameters
    ----------
    unit_definition : ~.UnitDefintion
        The unit definition to represent as a dictionary.

    Returns
    -------
    ~collections.OrderedDict
        A dictionary with elements corresponding to the unit definition
        attributes.

    See Also
    --------
    unit_from_dict

    """
    # Turn object into an OrderedDict
    new_unit_definition = OrderedDict()
    for key, value in iteritems(unit_definition.__dict__):
        if value and key != "list_of_units":
            new_unit_definition[key] = _fix_type(value)
        if value and key == "list_of_units":
            new_unit_definition[key] = []
            value = sorted(value, key=attrgetter("kind"))
            # Iterate through list of units and write to dict.
            for unit in value:
                new_unit = OrderedDict(
                    (k, _fix_type(v)) for k, v in iteritems(unit.__dict__)
                )
                new_unit_definition[key] += [new_unit]

    return new_unit_definition


def unit_from_dict(unit_definition):
    """Create an :class:`~.UnitDefintion` from a dictionary.

    Parameters
    ----------
    unit_definition : dict
        A dictionary with elements corresponding to the unit definition
        attributes.

    Returns
    -------
    UnitDefintion
        The generated unit definition.

    See Also
    --------
    unit_to_dict

    """
    # Create the new unit definition
    new_unit_definition = UnitDefinition()
    for key, value in iteritems(unit_definition):
        if key == "list_of_units":
            # Create Unit objects for units in the list_of_units attribute
            for unit in value:
                new_unit_definition.create_unit(
                    kind=unit["_kind"],
                    exponent=unit["_exponent"],
                    scale=unit["_scale"],
                    multiplier=unit["_multiplier"],
                )
        else:
            # Set attribute if not list of units
            new_unit_definition.__dict__[key] = value

    return new_unit_definition


# Internal
def _add_enzyme_module_form_attributes_into_dict(enzyme, new_enzyme):
    """Add EnzymeModuleForm attributes to its dict representation.

    Warnings
    --------
    This method is intended for internal use only.

    """
    # Add attributes to enzyme
    for attr in _REQUIRED_ENZYMEMODULEFORM_ATTRIBUTES:
        if attr == "enzyme_module_id":
            new_enzyme[attr] = getattr(enzyme, attr)
        else:
            bound = {str(k): v for k, v in iteritems(getattr(enzyme, attr))}
            new_enzyme[attr] = OrderedDict((k, bound[k]) for k in sorted(bound))
    # Update optional attributes
    _update_optional(
        enzyme,
        new_enzyme,
        _OPTIONAL_ENZYMEMODULEFORM_ATTRIBUTES,
        _ORDERED_OPTIONAL_ENZYMEMODULEFORM_KEYS,
    )


def _add_enzyme_module_reaction_attributes_into_dict(reaction, new_reaction):
    """Add EnzymeModuleReaction attributes to its dict representation.

    Warnings
    --------
    This method is intended for internal use only.

    """
    # Add ttributes to enzyme
    for attr in _REQUIRED_ENZYMEMODULEREACTION_ATTRIBUTES:
        new_reaction[attr] = getattr(reaction, attr)
    # Update optional attributes
    _update_optional(
        reaction,
        new_reaction,
        _OPTIONAL_ENZYMEMODULEREACTION_ATTRIBUTES,
        _ORDERED_OPTIONAL_ENZYMEMODULEREACTION_KEYS,
    )


def _add_enzyme_module_attributes_into_dict(model, obj):
    """Add EnzymeModule attributes to its dict representation.

    Warnings
    --------
    This method is intended for internal use only.

    """
    # Get a list of keys that should be represented as internal variables
    to_fix = ["_categorized", "enzyme_rate", "enzyme_concentration"]
    for key, value in iteritems(enzyme_to_dict(model)):
        if key in ("id", "name"):
            continue
        # Prefix the key so that it is an internal variable
        if any(
            [
                s in key
                for s in to_fix
                if s and key != "enzyme_concentration_total_equation"
            ]
        ):
            key = "_" + key
        obj[key] = value


def _fix_type(value):
    """Fix the type of the value so it can be exported to a file.

    Warnings
    --------
    This method is intended for internal use only.

    """
    if isinstance(value, (string_types, Basic)):
        return str(value)
    if isinstance(value, np.float_):
        return float(value)
    if isinstance(value, np.bool_):
        return bool(value)
    if isinstance(value, np.int_):
        return int(value)
    if isinstance(value, set):
        return sorted(list(value))
    if isinstance(value, dict):
        return OrderedDict((k, _fix_type(value[k])) for k in sorted(value))
    if value is None:
        return ""
    if isinstance(value, float) and abs(value) == _INF:
        return str(value)
    return value


def _update_optional(mass_object, new_dict, optional_attribute_dict, ordered_keys):
    """Update the dict to be saved with the optional attributes of the object.

    Warnings
    --------
    This method is intended for internal use only.

    """
    for key in ordered_keys:
        default = optional_attribute_dict[key]
        value = getattr(mass_object, key)
        if isinstance(value, (pd.DataFrame, np.ndarray)):
            if np.array(value).all() == np.array(default).all():
                continue
        elif value is None or value == default:
            continue
        else:
            pass
        new_dict[key] = _fix_type(value)


__all__ = (
    "model_to_dict",
    "model_from_dict",
    "metabolite_to_dict",
    "metabolite_from_dict",
    "reaction_to_dict",
    "reaction_from_dict",
    "enzyme_to_dict",
    "enzyme_from_dict",
    "unit_to_dict",
    "unit_from_dict",
)
