# -*- coding: utf-8 -*-
"""TODO Module Docstrings."""
from __future__ import absolute_import

from collections import OrderedDict
from operator import attrgetter, itemgetter

import numpy as np

import pandas as pd

from six import iteritems, iterkeys, string_types

from sympy import Basic, Eq, Symbol, sympify

from cobra.core import Gene

from mass.core import MassMetabolite, MassModel, MassReaction, UnitDefinition
from mass.enzyme_modules import (
    EnzymeModule, EnzymeModuleDict, EnzymeModuleForm, EnzymeModuleReaction,
    _ORDERED_ENZYMEMODULE_DICT_DEFAULTS)

# Global
_INF = float("inf")

_REQUIRED_REACTION_ATTRIBUTES = [
    "id", "name", "_reversible", "metabolites", "_lower_bound", "_upper_bound",
    "gene_reaction_rule"]
_ORDERED_OPTIONAL_REACTION_KEYS = [
    "subsystem", "steady_state_flux", "_forward_rate_constant",
    "_reverse_rate_constant", "_equilibrium_constant", "objective_coefficient",
    "variable_kind", "_rtype", "notes", "annotation"]
_OPTIONAL_REACTION_ATTRIBUTES = {
    "subsystem": "",
    "steady_state_flux": None,
    "_forward_rate_constant": None,
    "_equilibrium_constant": None,
    "_reverse_rate_constant": None,
    "objective_coefficient": 0,
    "variable_kind": "continuous",
    "_rtype": 1,
    "notes": {},
    "annotation": {}
}
_REQUIRED_ENZYMEMODULEREACTION_ATTRIBUTES = ["enzyme_module_id"]
_ORDERED_OPTIONAL_ENZYMEMODULEREACTION_KEYS = []
_OPTIONAL_ENZYMEMODULEREACTION_ATTRIBUTES = {}

_REQUIRED_METABOLITE_ATTRIBUTES = ["id", "name"]
_ORDERED_OPTIONAL_METABOLITE_KEYS = [
    "formula", "charge", "compartment", "fixed", "_initial_condition",
    "_constraint_sense", "_bound", "notes", "annotation"]
_OPTIONAL_METABOLITE_ATTRIBUTES = {
    "charge": None,
    "formula": None,
    "compartment": None,
    "fixed": False,
    "_initial_condition": None,
    "_bound": 0,
    "_constraint_sense": "E",
    "notes": {},
    "annotation": {}
}

_REQUIRED_ENZYMEMODULEFORM_ATTRIBUTES = [
    "_bound_catalytic", "_bound_effectors", "enzyme_module_id"]
_ORDERED_OPTIONAL_ENZYMEMODULEFORM_KEYS = []
_OPTIONAL_ENZYMEMODULEFORM_ATTRIBUTES = {}

_REQUIRED_GENE_ATTRIBUTES = ["id", "name"]
_ORDERED_OPTIONAL_GENE_KEYS = ["notes", "annotation"]
_OPTIONAL_GENE_ATTRIBUTES = {
    "notes": {},
    "annotation": {}
}

_REQUIRED_ENZYMEMODULE_ATTRIBUTES = [
    "id", "name", "enzyme_module_ligands", "enzyme_module_forms",
    "enzyme_module_reactions"]
_ORDERED_OPTIONAL_ENZYMEMODULE_KEYS = [
    key for key in iterkeys(_ORDERED_ENZYMEMODULE_DICT_DEFAULTS)
    if key not in _REQUIRED_ENZYMEMODULE_ATTRIBUTES + ["S", "model"]]
_OPTIONAL_ENZYMEMODULE_ATTRIBUTES = OrderedDict({
    key: _ORDERED_ENZYMEMODULE_DICT_DEFAULTS[key]
    for key in _ORDERED_OPTIONAL_ENZYMEMODULE_KEYS
})

_ORDERED_OPTIONAL_MODEL_KEYS = [
    "name", "description", "boundary_conditions", "custom_parameters",
    "compartments", "notes", "annotation"]
_OPTIONAL_MODEL_ATTRIBUTES = {
    "name": None,
    "description": "",
    "boundary_conditions": {},
    "custom_parameters": {},
    "compartments": {},
    "notes": {},
    "annotation": {}
}


def metabolite_to_dict(metabolite):
    """Represent a MassMetabolite object as a dict.

    Parameters
    ----------
    metabolite: MassMetabolite
        The metabolite to represent as a dict.

    """
    # Turn object into an OrderedDict with the required attributes
    new_met = OrderedDict()
    for key in _REQUIRED_METABOLITE_ATTRIBUTES:
        new_met[key] = _fix_type(getattr(metabolite, key))
    # Update with any opitonal attributes that are not their defaults.
    _update_optional(metabolite, new_met, _OPTIONAL_METABOLITE_ATTRIBUTES,
                     _ORDERED_OPTIONAL_METABOLITE_KEYS)

    # Add EnzymeModuleForm attributes if metabolite is an EnzymeModuleForm
    if isinstance(metabolite, EnzymeModuleForm):
        _add_enzyme_module_form_attributes_into_dict(metabolite, new_met)

    return new_met


def metabolite_from_dict(metabolite):
    """Create a MassMetabolite object from its dict representation.

    Parameters
    ----------
    metabolite: dict
        The dict representation of the metabolite to create.

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
    """Represent a MassReaction object as a dict.

    Parameters
    ----------
    reaction: MassReaction
        The reaction to represent as a dict.

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
    _update_optional(reaction, new_reaction, _OPTIONAL_REACTION_ATTRIBUTES,
                     _ORDERED_OPTIONAL_REACTION_KEYS)
    # Add EnzymeModuleReaction attributes
    # if reaction is an EnzymeModuleReaction
    if isinstance(reaction, EnzymeModuleReaction):
        _add_enzyme_module_reaction_attributes_into_dict(
            reaction, new_reaction)

    return new_reaction


def reaction_from_dict(reaction, model):
    """Create a MassReaction object from its dict representation.

    Parameters
    ----------
    reaction: dict
        The dict representation of the reaction to create.
    model: MassModel
        The MassModel to assoicate with the reaction.

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
            new_reaction.add_metabolites(OrderedDict(
                (model.metabolites.get_by_id(str(met)), coeff)
                for met, coeff in iteritems(v)))
        else:
            setattr(new_reaction, k, v)

    return new_reaction


def gene_to_dict(gene):
    """Represent a gene object as a dict.

    Parameters
    ----------
    gene: cobra.Gene
        The gene to represent as a dict.

    """
    # Turn object into an OrderedDict with the required attributes
    new_gene = OrderedDict()
    for key in _REQUIRED_GENE_ATTRIBUTES:
        new_gene[key] = _fix_type(getattr(gene, key))
    # Update with any opitonal attributes that are not their defaults.
    _update_optional(gene, new_gene, _OPTIONAL_GENE_ATTRIBUTES,
                     _ORDERED_OPTIONAL_GENE_KEYS)
    return new_gene


def gene_from_dict(gene):
    """Create a gene object from its dict representation.

    Parameters
    ----------
    gene: dict
        The dict representation of the gene to create.

    """
    new_gene = Gene(gene["id"])
    # Set object attributes
    for k, v in iteritems(gene):
        setattr(new_gene, k, v)
    return new_gene


def enzyme_to_dict(enzyme):
    """Represent an EnzymeModuleDict object as a dict.

    Parameters
    ----------
    enzyme: EnzymeModuleDict
        The enzyme to represent as a dict.

    """
    # Turn object into an OrderedDict with the required attributes
    new_enzyme = OrderedDict()
    for key in _REQUIRED_ENZYMEMODULE_ATTRIBUTES:
        new_enzyme[key] = _fix_type(getattr(enzyme, key))

    # Update with any opitonal attributes that are not their defaults.
    _update_optional(enzyme, new_enzyme, _OPTIONAL_ENZYMEMODULE_ATTRIBUTES,
                     _ORDERED_OPTIONAL_ENZYMEMODULE_KEYS)

    # Store objects and expressions as string for the attributes
    for key in _REQUIRED_ENZYMEMODULE_ATTRIBUTES[2:]:
        new_enzyme[key] = [i.id for i in getattr(enzyme, key)]
        # Repeat for categorized attribute
        key += "_categorized"
        if getattr(enzyme, key) != _OPTIONAL_ENZYMEMODULE_ATTRIBUTES[key]:
            new_enzyme[key] = {
                category: [i.id for i in old_dictlist]
                for category, old_dictlist in iteritems(getattr(enzyme, key))}

    key = "enzyme_net_flux_equation"
    if key in new_enzyme:
        new_enzyme[key] = str(getattr(enzyme, key).rhs)

    return new_enzyme


def enzyme_from_dict(enzyme, model):
    """Create an EnzymeModuleDict from its dict representation.

    Parameters
    ----------
    enzyme: dict
        The dict representation of the gene to create.

    """
    # Set object attributes
    new_enzyme = EnzymeModuleDict(id_or_enzyme=enzyme)
    # Update model and get objects from the model to populate the DictLists
    new_enzyme["model"] = model
    new_enzyme._update_object_pointers(model)

    # Make the enzyme equations
    new_enzyme.enzyme_net_flux_equation = Eq(Symbol(
        "v_" + new_enzyme.id), sympify(new_enzyme.enzyme_net_flux_equation))

    # Make the stoichiometric matrix and clean up the EnzymeModuleDict
    new_enzyme._make_enzyme_stoichiometric_matrix(update=True)
    new_enzyme._set_missing_to_defaults()
    new_enzyme._fix_order()

    return new_enzyme


def unit_to_dict(unit_definition):
    """Represent a UnitDefinition object as a dict.

    Parameters
    ----------
    unit: UnitDefintion
        The UnitDefinition to represent as a dict.

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
                    (k, _fix_type(v)) for k, v in iteritems(unit.__dict__))
                new_unit_definition[key] += [new_unit]

    return new_unit_definition


def unit_from_dict(unit_definition):
    """Create a UnitDefinition object from its dict representation.

    Parameters
    ----------
    unit: dict
        The dict representation of the UnitDefinition to create.

    """
    # Create the new unit definition
    new_unit_definition = UnitDefinition()
    for key, value in iteritems(unit_definition):
        if key == "list_of_units":
            # Create Unit objects for units in the list_of_units attribute
            for unit in value:
                new_unit_definition.create_unit(
                    kind=unit["_kind"], exponent=unit["_exponent"],
                    scale=unit["_scale"], multiplier=unit["_multiplier"])
        else:
            # Set attribute if not list of units
            new_unit_definition.__dict__[key] = value

    return new_unit_definition


def model_to_dict(model, sort=False):
    """Represent a MassModel object as a dict.

    Parameters
    ----------
    model: MassModel, EnzymeModule
        The model to represent as a dict.
    sort: bool, optional
        Whether to sort the metabolites, reactions, genes, and enzyme_modules
        or maintain the order defined in the model. If the model is an
        EnzymeModule, the enzyme_module_ligands, enzyme_module_forms, and
        enzyme_module_reactions attributes are also included. Default is False.

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
    custom_rates = getattr(model, "custom_rates", {})
    if custom_rates:
        custom_rates = OrderedDict(
            (k.id, _fix_type(v)) for k, v in iteritems(custom_rates))
        if sort:
            custom_rates = OrderedDict(
                (k, custom_rates[k]) for k in sorted(custom_rates))
        obj["custom_rates"] = custom_rates

    _update_optional(model, obj, _OPTIONAL_MODEL_ATTRIBUTES,
                     _ORDERED_OPTIONAL_MODEL_KEYS)
    # Add EnzymeModule attributes if an EnzymeModule is being saved.
    if isinstance(model, EnzymeModule):
        _add_enzyme_module_attributes_into_dict(model, obj)
        if sort:
            obj["enzyme_module_ligands"].sort(key=get_id)
            obj["enzyme_module_forms"].sort(key=get_id)
            obj["enzyme_module_reactions"].sort(key=get_id)
    return obj


def model_from_dict(obj):
    """Create a MassModel object from its dict representation.

    Parameters
    ----------
    obj: dict
        The dict representation of the MassModel to create.

    """
    if "reactions" not in obj:
        raise ValueError("Object has no reactions attribute. Cannot load.")
    if all([k in obj for k in _REQUIRED_ENZYMEMODULE_ATTRIBUTES[2:]]):
        model = EnzymeModule(obj["id"])
    else:
        model = MassModel(obj["id"])
    # Add metabolites to the model
    model.add_metabolites([
        metabolite_from_dict(metabolite) for metabolite in obj["metabolites"]])

    # Add genes to the model
    model.genes.extend([gene_from_dict(gene) for gene in obj["genes"]])

    # Add reactions to the model
    model.add_reactions([
        reaction_from_dict(reaction, model) for reaction in obj["reactions"]])

    # Add units to the model
    model.add_units([unit_from_dict(unit_def) for unit_def in obj["units"]])
    # Add enzyme modules to the model
    if "enzyme_modules" in obj:
        model.enzyme_modules.extend([
            enzyme_from_dict(enz, model) for enz in obj["enzyme_modules"]])
    # Repair model once all objects are in model.
    model.repair(rebuild_index=True, rebuild_relationships=True)
    # Add boundary conditions to the model if they exist
    if "boundary_conditions" in obj:
        model.add_boundary_conditions({
            met: bc for met, bc in iteritems(obj["boundary_conditions"])})

    # Get custom parameters if they exist
    if "custom_parameters" in obj:
        model.custom_parameters.update(obj["custom_parameters"])
    # Add custom rates and any custom parameters if they exist
    if "custom_rates" in obj:
        # Add custom rates to the model
        for reaction, custom_rate in iteritems(obj["custom_rates"]):
            model.add_custom_rate(model.reactions.get_by_id(reaction),
                                  custom_rate, model.custom_parameters)
    # Update with any opitonal attributes.
    for k, v in iteritems(obj):
        # Set MassModel attributes (and subsystem attribute for EnzymeModules)
        if k in _ORDERED_OPTIONAL_MODEL_KEYS or k == "subsystem":
            setattr(model, k, v)
        # Update with EnzymeModule attributes if obj represents an EnzymeModule
        elif k.lstrip("_") in _ORDERED_OPTIONAL_ENZYMEMODULE_KEYS:
            model.__class__.__dict__[k.lstrip("_")].fset(model, v)

    return model


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
            new_enzyme[attr] = OrderedDict((k, bound[k])
                                           for k in sorted(bound))
    # Update optional attributes
    _update_optional(enzyme, new_enzyme, _OPTIONAL_ENZYMEMODULEFORM_ATTRIBUTES,
                     _ORDERED_OPTIONAL_ENZYMEMODULEFORM_KEYS)


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
    _update_optional(reaction, new_reaction,
                     _OPTIONAL_ENZYMEMODULEREACTION_ATTRIBUTES,
                     _ORDERED_OPTIONAL_ENZYMEMODULEREACTION_KEYS)


def _add_enzyme_module_attributes_into_dict(model, obj):
    """Add EnzymeModule attributes to its dict representation.

    Warnings
    --------
    This method is intended for internal use only.

    """
    # Get a list of keys that should be represented as internal variables
    to_fix = ["_categorized", "enzyme_net_flux", "enzyme_concentration"]
    for key, value in iteritems(enzyme_to_dict(model)):
        if key in("id", "name"):
            continue
        # Prefix the key so that it is an internal variable
        if any([s in key for s in to_fix]):
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


def _update_optional(mass_object, new_dict, optional_attribute_dict,
                     ordered_keys):
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
