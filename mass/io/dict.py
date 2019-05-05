# -*- coding: utf-8 -*-
"""TODO Module Docstrings."""
from __future__ import absolute_import

from collections import OrderedDict
from operator import attrgetter, itemgetter

import numpy as np

import pandas as pd

from six import iteritems, iterkeys, string_types

from sympy import Eq, Symbol, sympify

from cobra.core import Gene

from mass.core import MassMetabolite, MassModel, MassReaction
from mass.enzyme_modules import (
    EnzymeModule, EnzymeModuleDict, EnzymeModuleForm,
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
    "_reverse_rate_constant": 0,
    "objective_coefficient": 0,
    "variable_kind": "continuous",
    "_rtype": 1,
    "notes": {},
    "annotation": {}
}

_REQUIRED_METABOLITE_ATTRIBUTES = ["id", "name"]
_ORDERED_OPTIONAL_METABOLITE_KEYS = [
    "formula", "charge", "compartment", "_initial_condition",
    "_constraint_sense", "_bound", "notes", "annotation"]
_OPTIONAL_METABOLITE_ATTRIBUTES = {
    "charge": None,
    "formula": None,
    "compartment": None,
    "_initial_condition": None,
    "_bound": 0,
    "_constraint_sense": "E",
    "notes": {},
    "annotation": {}
}

_REQUIRED_ENZYMEMODULEFORM_ATTRIBUTES = [
    "_bound_catalytic", "_bound_effectors"]
_ORDERED_OPTIONAL_ENZYMEMODULEFORM_KEYS = ["enzyme_id", "enzyme_name"]
_OPTIONAL_ENZYMEMODULEFORM_ATTRIBUTES = {
    "enzyme_id": "",
    "enzyme_name": "",
}

_REQUIRED_GENE_ATTRIBUTES = ["id", "name"]
_ORDERED_OPTIONAL_GENE_KEYS = ["notes", "annotation"]
_OPTIONAL_GENE_ATTRIBUTES = {
    "notes": {},
    "annotation": {}
}

_REQUIRED_ENZYME_ATTRIBUTES = [
    "id", "name", "ligands", "enzyme_forms", "enzyme_reactions"]
_ORDERED_OPTIONAL_ENZYME_KEYS = [
    key for key in iterkeys(_ORDERED_ENZYMEMODULE_DICT_DEFAULTS)
    if key not in _REQUIRED_ENZYME_ATTRIBUTES + ["S", "model"]]
_OPTIONAL_ENZYME_ATTRIBUTES = OrderedDict({
    key: _ORDERED_ENZYMEMODULE_DICT_DEFAULTS[key]
    for key in _ORDERED_OPTIONAL_ENZYME_KEYS
})

_ORDERED_OPTIONAL_MODEL_KEYS = [
    "name", "description", "compartments", "modules", "units", "notes",
    "annotation"]
_OPTIONAL_MODEL_ATTRIBUTES = {
    "name": None,
    "description": "",
    "compartments": {},
    "modules": set(),
    "units": {},
    "notes": {},
    "annotation": {}
}


def metabolite_to_dict(metabolite):
    """Represent a MassMetabolite object as a dict.

    Parameters
    ----------
    metabolite: mass.MassMetabolite
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
        _add_enzyme_form_attributes_into_dict(metabolite, new_met)

    return new_met


def metabolite_from_dict(metabolite):
    """Create a MassMetabolite object from its dict representation.

    Parameters
    ----------
    metabolite: dict
        The dict representation of the metabolite to create.

    """
    # Determine if saved object should be a MassMetabolite or a subclass
    if "_bound_catalytic" in metabolite or "_bound_effectors" in metabolite:
        new_metabolite = EnzymeModuleForm(id=metabolite["id"])
    else:
        new_metabolite = MassMetabolite(id=metabolite["id"])

    # Set object attributes
    for k, v in iteritems(metabolite):
        setattr(new_metabolite, k, v)

    return new_metabolite


def reaction_to_dict(reaction):
    """Represent a MassReaction object as a dict.

    Parameters
    ----------
    reaction: mass.MassReaction
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

    return new_reaction


def reaction_from_dict(reaction, model):
    """Create a MassReaction object from its dict representation.

    Parameters
    ----------
    reaction: dict
        The dict representation of the reaction to create.
    model: mass.MassModel
        The MassModel to assoicate with the reaction.

    """
    new_reaction = MassReaction(id=reaction["id"])
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
    new_gene = Gene(id=gene["id"])
    # Set object attributes
    for k, v in iteritems(gene):
        setattr(new_gene, k, v)
    return new_gene


def enzyme_to_dict(enzyme):
    """Represent an EnzymeModuleDict object as a dict.

    Parameters
    ----------
    enzyme: enzyme.EnzymeModuleDict
        The enzyme to represent as a dict.

    """
    # Turn object into an OrderedDict with the required attributes
    new_enzyme = OrderedDict()
    for key in _REQUIRED_ENZYME_ATTRIBUTES:
        new_enzyme[key] = _fix_type(getattr(enzyme, key))

    # Update with any opitonal attributes that are not their defaults.
    _update_optional(enzyme, new_enzyme, _OPTIONAL_ENZYME_ATTRIBUTES,
                     _ORDERED_OPTIONAL_ENZYME_KEYS)

    # Store objects and expressions as string for the attributes
    for key in _REQUIRED_ENZYME_ATTRIBUTES[2:]:
        new_enzyme[key] = [i.id for i in getattr(enzyme, key)]
        # Repeat for categorized attribute
        key = "categorized_" + key
        if getattr(enzyme, key) != _OPTIONAL_ENZYME_ATTRIBUTES[key]:
            new_enzyme[key] = {
                category: [i.id for i in old_dictlist]
                for category, old_dictlist in iteritems(getattr(enzyme, key))}

    key = "enzyme_net_flux_equation"
    if key in new_enzyme:
        new_enzyme[key] = str(new_enzyme[key].rhs)

    return new_enzyme


def enzyme_from_dict(enzyme, model):
    """Create an EnzymeModuleDict from its dict representation.

    Parameters
    ----------
    gene: dict
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


def model_to_dict(model, sort=False):
    """Represent a MassModel object as a dict.

    Parameters
    ----------
    model: mass.MassModel, mass.EnzymeModule
        The model to represent as a dict.
    sort: bool, optional
        Whether to sort the metabolites, reactions, genes, and enzyme_modules
        or maintain the order defined in the model. If the model is an
        EnzymeModule, the ligands, enzyme_forms, and enzyme_reactions
        attributes are also included. Default is False.

    """
    obj = OrderedDict()
    obj["id"] = model.id
    obj["metabolites"] = list(map(metabolite_to_dict, model.metabolites))
    obj["reactions"] = list(map(reaction_to_dict, model.reactions))
    obj["genes"] = list(map(gene_to_dict, model.genes))
    obj["enzyme_modules"] = list(map(enzyme_to_dict, model.enzyme_modules))

    if sort:
        get_id = itemgetter("id")
        obj["metabolites"].sort(key=get_id)
        obj["reactions"].sort(key=get_id)
        obj["genes"].sort(key=get_id)
        obj["enzyme_modules"].sort(key=get_id)

    for key in ["initial_conditions", "fixed_concentrations",
                "custom_rates", "custom_parameters"]:
        values = getattr(model, key, {})
        if values:
            values = OrderedDict(
                (k.id, str(v)) if key == "custom_rates"
                else (str(k), _fix_type(v)) for k, v in iteritems(values))
            if sort:
                values = OrderedDict((k, values[k]) for k in sorted(values))
            obj[key] = values

    _update_optional(model, obj, _OPTIONAL_MODEL_ATTRIBUTES,
                     _ORDERED_OPTIONAL_MODEL_KEYS)
    # Add EnzymeModule attributes if an EnzymeModule is being saved.
    if isinstance(model, EnzymeModule):
        _add_enzyme_module_attributes_into_dict(model, obj)
        if sort:
            obj["ligands"].sort(key=get_id)
            obj["enzyme_forms"].sort(key=get_id)
            obj["enzyme_reactions"].sort(key=get_id)

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
    if all([k in obj for k in _REQUIRED_ENZYME_ATTRIBUTES[2:]]):
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

    # Add enzyme modules to the model
    if "enzyme_modules" in obj:
        model.enzyme_modules.extend([
            enzyme_from_dict(enz, model) for enz in obj["enzyme_modules"]])

    # Add initial conditions to the model if they exist
    if "initial_conditions" in obj:
        model.update_initial_conditions({
            model.metabolites.get_by_id(met): ic
            for met, ic in iteritems(obj["initial_conditions"])})

    # Add fixed concentrations to the model if they exist
    if "fixed_concentrations" in obj:
        model.add_fixed_concentrations(
            dict((model.metabolites.get_by_id(met), ic)
                 if met in model.metabolites else (str(met), ic)
                 for met, ic in iteritems(obj["fixed_concentrations"])))

    # Add custom rates and any custom parameters if they exist
    if "custom_rates" in obj:
        custom_parameters = {}
        # Get custom parameters if they exist
        if "custom_parameters" in obj:
            custom_parameters.update(obj["custom_parameters"])
        # Add custom rates to the model
        for reaction, custom_rate in iteritems(obj["custom_rates"]):
            model.add_custom_rate(model.reactions.get_by_id(reaction),
                                  custom_rate, custom_parameters)
    # Update with any opitonal attributes.
    for k, v in iteritems(obj):
        # Set MassModel attributes (and subsystem attribute for EnzymeModules)
        if k in _ORDERED_OPTIONAL_MODEL_KEYS or k == "subsystem":
            setattr(model, k, v)
        # Update with EnzymeModule attributes if obj represents an EnzymeModule
        elif k.lstrip("_") in _ORDERED_OPTIONAL_ENZYME_KEYS:
            model.__class__.__dict__[k.lstrip("_")].fset(model, v)

    model.modules = set(sorted(model.modules))

    return model


# Internal
def _add_enzyme_form_attributes_into_dict(enzyme, new_enzyme):
    """Add EnzymeModuleForm attributes to its dict representation.

    Warnings
    --------
    This method is intended for internal use only.

    """
    # Add bound_catalytic and bound_effectors attributes to enzyme
    for attr in _REQUIRED_ENZYMEMODULEFORM_ATTRIBUTES:
        bound = {str(k): v for k, v in iteritems(getattr(enzyme, attr))}
        new_enzyme[attr] = OrderedDict((k, bound[k]) for k in sorted(bound))
    # Update optional attributes
    _update_optional(enzyme, new_enzyme, _OPTIONAL_ENZYMEMODULEFORM_ATTRIBUTES,
                     _ORDERED_OPTIONAL_ENZYMEMODULEFORM_KEYS)


def _add_enzyme_module_attributes_into_dict(model, obj):
    """Add EnzymeModule attributes to its dict representation.

    Warnings
    --------
    This method is intended for internal use only.

    """
    # Get a list of keys that should be represented as internal variables
    to_fix = ["categorized_", "enzyme_net_flux", "enzyme_concentration"]
    for key, value in iteritems(enzyme_to_dict(model)):
        if key == "id" or key == "name":
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
    if isinstance(value, string_types):
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
        return OrderedDict((key, value[key]) for key in sorted(value))
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
