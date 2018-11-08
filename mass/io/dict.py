# -*- coding: utf-8 -*-
"""TODO Module Docstrings."""
from __future__ import absolute_import


from collections import OrderedDict
from operator import attrgetter, itemgetter

from cobra.core import Gene

from mass.core import MassMetabolite, MassModel, MassReaction

import numpy as np

from six import iteritems, string_types

# Global
_INF = float("inf")

_REQUIRED_REACTION_ATTRIBUTES = [
    "id", "name", "_reversible", "metabolites", "lower_bound", "upper_bound",
    "gene_reaction_rule"]
_ORDERED_OPTIONAL_REACTION_KEYS = [
    "subsystem", "steady_state_flux", "_forward_rate_constant",
    "_reverse_rate_constant", "_equilibrium_constant", "objective_coefficient",
    "variable_kind", "_rtype",  "notes", "annotation"]
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

_REQUIRED_GENE_ATTRIBUTES = ["id", "name"]
_ORDERED_OPTIONAL_GENE_KEYS = ["notes", "annotation"]
_OPTIONAL_GENE_ATTRIBUTES = {
    "notes": {},
    "annotation": {}
}

_ORDERED_OPTIONAL_MODEL_KEYS = [
    "name", "compartments", "modules", "units", "notes", "annotation"]
_OPTIONAL_MODEL_ATTRIBUTES = {
    "name": None,
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
    new_met = OrderedDict()
    for key in _REQUIRED_METABOLITE_ATTRIBUTES:
        new_met[key] = _fix_type(getattr(metabolite, key))
    _update_optional(metabolite, new_met, _OPTIONAL_METABOLITE_ATTRIBUTES,
                     _ORDERED_OPTIONAL_METABOLITE_KEYS)
    return new_met


def metabolite_from_dict(metabolite):
    """Create a MassMetabolite object from its dict representation.

    Parameters
    ----------
    metabolite: dict
        The dict representation of the metabolite to create.

    """
    new_metabolite = MassMetabolite()
    for k, v in iteritems(metabolite):
        setattr(new_metabolite, k, v)
    return new_metabolite


def gene_to_dict(gene):
    """Represent a gene object as a dict.

    Parameters
    ----------
    gene: cobra.Gene
        The gene to represent as a dict.

    """
    new_gene = OrderedDict()
    for key in _REQUIRED_GENE_ATTRIBUTES:
        new_gene[key] = _fix_type(getattr(gene, key))
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
    for k, v in iteritems(gene):
        setattr(new_gene, k, v)
    return new_gene


def reaction_to_dict(reaction):
    """Represent a MassReaction object as a dict.

    Parameters
    ----------
    reaction: mass.MassReaction
        The reaction to represent as a dict.

    """
    new_reaction = OrderedDict()
    for key in _REQUIRED_REACTION_ATTRIBUTES:
        if key != "metabolites":
            new_reaction[key] = _fix_type(getattr(reaction, key))
            continue

        mets = OrderedDict()
        for met in sorted(reaction.metabolites, key=attrgetter("id")):
            mets[str(met)] = reaction.metabolites[met]
        new_reaction["metabolites"] = mets
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
    for k, v in iteritems(reaction):
        # Change infinity type from a string to a float
        if isinstance(v, string_types) and "inf" in v:
            v = float(v)
        if k in {"objective_coefficient", "reaction"}:
            continue
        elif k == "metabolites":
            new_reaction.add_metabolites(OrderedDict(
                (model.metabolites.get_by_id(str(met)), coeff)
                for met, coeff in iteritems(v)))
        else:
            setattr(new_reaction, k, v)
    return new_reaction


def model_to_dict(model, sort=False):
    """Represent a MassModel object as a dict.

    Parameters
    ----------
    model: mass.MassModel
        The reaction to represent as a dict.
    sort: bool, optional
        Whether to sort the metabolites, reactions, and genes or maintain the
        order defined in the model. Default is false.

    """
    obj = OrderedDict()
    obj["id"] = model.id
    obj["metabolites"] = list(map(metabolite_to_dict, model.metabolites))
    obj["reactions"] = list(map(reaction_to_dict, model.reactions))
    obj["genes"] = list(map(gene_to_dict, model.genes))

    if sort:
        get_id = itemgetter("id")
        obj["metabolites"].sort(key=get_id)
        obj["reactions"].sort(key=get_id)
        obj["genes"].sort(key=get_id)

    for key in ["initial_conditions", "fixed_concentrations",
                "custom_rates", "custom_parameters"]:
        values = getattr(model, key, {})
        if values:
            values = OrderedDict((k.id, str(v))
                                 if key is "custom_rates" else (str(k), v)
                                 for k, v in iteritems(values))
            if sort:
                values = OrderedDict((k, values[k]) for k in sorted(values))
            obj[key] = values

    _update_optional(model, obj, _OPTIONAL_MODEL_ATTRIBUTES,
                     _ORDERED_OPTIONAL_MODEL_KEYS)
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
    model = MassModel()
    # Add metabolites to the model
    model.add_metabolites([
        metabolite_from_dict(metabolite) for metabolite in obj["metabolites"]])

    # Add genes to the model
    model.genes.extend([gene_from_dict(gene) for gene in obj["genes"]])

    # Add reactions to the model
    model.add_reactions([
        reaction_from_dict(reaction, model) for reaction in obj["reactions"]])

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

    for k, v in iteritems(obj):
        if k == "id" or k in _ORDERED_OPTIONAL_MODEL_KEYS:
            setattr(model, k, v)

    return model


# Internal
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
        return list(value)
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
        if value is None or value == default:
            continue
        new_dict[key] = _fix_type(value)
