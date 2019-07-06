# -*- coding: utf-8 -*-
"""TODO Module Docstrings."""
from six import iteritems, string_types

from cobra.core.metabolite import Metabolite
from cobra.core.model import Model
from cobra.core.reaction import Reaction

from mass.core.massmetabolite import MassMetabolite
from mass.core.massmodel import MassModel
from mass.core.massreaction import MassReaction


# Public
def convert_mass_to_cobra(obj, prefix=None):
    """Create a new cobra Object from a given mass Object.

    Parameters
    ----------
    obj: MassMetabolite, MassModel, or MassReaction
        The mass Object to convert into a new cobra Object.
    prefix: str, optional
        If provided, the string is used to prefix the identifier of the
        returned object and any assoicated objects.

    Returns
    -------
    Object of class 'cobra'
        The new 'cobra' object.

    Warnings
    --------
    All similar fields will initialize to be identical to the input object,
    and all other fields will initialize to default values.

    """
    conversion_dict = {MassMetabolite: mass_to_cobra_metabolite,
                       MassReaction: mass_to_cobra_reaction,
                       MassModel: mass_to_cobra_model}
    try:
        conversion_function = conversion_dict[type(obj)]
    except KeyError:
        raise TypeError("object must be a mass object")

    return conversion_function(obj, prefix=prefix)


def convert_cobra_to_mass(obj, prefix=None):
    """Create a new mass Object from a given cobra Object.

    Parameters
    ----------
    object: Metabolite, Model, or Reaction
        The cobra Object to convert into a new mass Object.
    prefix: str, optional
        If provided, the string is used to prefix the identifier of the
        returned object and any assoicated objects.

    Returns
    -------
    new_object: Object of class 'mass'
        The new 'mass' object.

    Warnings
    --------
    All similar fields will initialize to be identical to the input object,
    and all other fields will initialize to default values.

    """
    conversion_dict = {Metabolite: cobra_to_mass_metabolite,
                       Reaction: cobra_to_mass_reaction,
                       Model: cobra_to_mass_model}
    try:
        conversion_function = conversion_dict[type(obj)]
    except KeyError:
        raise TypeError("obj must be a cobra object")

    return conversion_function(obj, prefix=prefix)


def mass_to_cobra_metabolite(metabolite, prefix=None):
    """Create a new cobra.Metabolite from the given mass.MassMetabolite.

    Parameters
    ----------
    metabolite: mass.MassMetabolite
        The mass.MassMetabolite to convert into a new cobra.Metabolite.
    prefix: str, optional
        If provided, the string is used to prefix the identifier of the
        returned object and any assoicated objects.

    Returns
    -------
    cobra.Metabolite
        The new cobra.Metabolite object.

    Warnings
    --------
    All similar fields will initialize to be identical to the input object,
    and all other fields will initialize to default values.

    """
    return _convert_metabolite(metabolite, Metabolite, MassMetabolite, prefix)


def cobra_to_mass_metabolite(metabolite, prefix=None):
    """Create a new mass.MassMetabolite from the given cobra.Metabolite.

    Parameters
    ----------
    metabolite: cobra.Metabolite
        The cobra.Metabolite to convert into a new mass.MassMetabolite.
    prefix: str, optional
        If provided, the string is used to prefix the identifier of the
        returned object and any assoicated objects.

    Returns
    -------
    mass.MassMetabolite
        The new mass.MassMetabolite object.

    Warnings
    --------
    All similar fields will initialize to be identical to the input object,
    and all other fields will initialize to default values.

    """
    return _convert_metabolite(metabolite, MassMetabolite, Metabolite, prefix)


def mass_to_cobra_reaction(reaction, prefix=None, convert_metabolites=True,
                           lower_bound=None, upper_bound=None):
    """Create a new cobra.Reaction from the given mass.MassReaction.

    Parameters
    ----------
    reaction: mass.MassReaction
        The mass.MassReaction to convert into a new cobra.Reaction.
    prefix: str, optional
        If provided, the string is used to prefix the identifier of the
        returned object and any assoicated objects.

    Returns
    -------
    cobra.Reaction
        The new cobra.Reaction object.

    Notes
    -----
    The lower and upper bounds for the cobra.Reaction are inferred from the
        kinetic reversibility of the mass.MassReaction if not provided.

    Warnings
    --------
    All similar fields will initialize to be identical to the input object,
    and all other fields will initialize to default values.

    """
    new_reaction = _convert_reaction(reaction, Reaction, MassReaction, prefix)
    if lower_bound is not None:
        new_reaction.lower_bound = lower_bound
    if upper_bound is not None:
        new_reaction.upper_bound = upper_bound

    # Add converted metabolites
    if convert_metabolites:
        new_metabolites = {mass_to_cobra_metabolite(met, prefix=prefix): c
                           for met, c in iteritems(reaction._metabolites)}
        new_reaction.add_metabolites(new_metabolites)
    return new_reaction


def cobra_to_mass_reaction(reaction, prefix=None, convert_metabolites=True,
                           kinetic_reversibility=None):
    """Create a new mass.MassReaction from the given cobra.Reaction.

    Parameters
    ----------
    reaction: cobra.Reaction
        The cobra.Reaction to convert into a new mass.MassReaction.
    prefix: str, optional
        If provided, the string is used to prefix the identifier of the
        returned object and any assoicated objects.

    Returns
    -------
    mass.MassReaction
        The new mass.MassReaction object.

    Notes
    -----
    Kinetic reversibility for the mass.MassReaction is inferred from the lower
        and upper bounds of the cobra.Reaction if not provided.


    Warnings
    --------
    All similar fields will initialize to be identical to the input object,
    and all other fields will initialize to default values.

    """
    new_reaction = _convert_reaction(reaction, MassReaction, Reaction, prefix)
    if kinetic_reversibility is not None:
        new_reaction.reversible = kinetic_reversibility

    # Add converted metabolites
    if convert_metabolites:
        new_metabolites = {cobra_to_mass_metabolite(met, prefix=prefix): c
                           for met, c in iteritems(reaction._metabolites)}
        new_reaction.add_metabolites(new_metabolites)
    return new_reaction


def mass_to_cobra_model(model, prefix=None):
    """Create a new cobra.Model from the given mass.MassModel.

    Parameters
    ----------
    model: mass.MassModel
        The mass.MassModel to convert into a new cobra.Model.
    prefix: str, optional
        If provided, the string is used to prefix the identifier of the
        returned object and any assoicated objects.

    Returns
    -------
    cobra.Model
        The new cobra.Model object.

    Warnings
    --------
    All similar fields will initialize to be identical to the input object,
    and all other fields will initialize to default values.

    """
    new_model = _convert_model(model, Model, MassModel, prefix)
    # Add converted reactions and their converted metabolites
    new_reactions = [mass_to_cobra_reaction(rxn, prefix=prefix)
                     for rxn in model.reactions]
    new_model.add_reactions(new_reactions)
    # Copy compartments
    new_model.compartments = model.compartments.copy()
    new_model.repair(rebuild_index=True, rebuild_relationships=True)
    return new_model


def cobra_to_mass_model(model, prefix=None, kinetic_reversibility=None):
    """Create a new mass.MassModel from the given cobra.Model.

    Parameters
    ----------
    model: cobra.Model
        The cobra.Reaction to convert into a new mass.MassModel.
    prefix: str, optional
        If provided, the string is used to prefix the identifier of the
        returned object and any assoicated objects.

    Returns
    -------
    mass.MassModel
        The new mass.MassModel object.

    Warnings
    --------
    All similar fields will initialize to be identical to the input object,
    and all other fields will initialize to default values.

    """
    new_model = _convert_model(model, MassModel, Model, prefix)
    # Add converted reactions and their converted metabolites
    new_reactions = [
        cobra_to_mass_reaction(rxn, prefix=prefix,
                               kinetic_reversibility=kinetic_reversibility)
        for rxn in model.reactions]
    new_model.add_reactions(new_reactions)
    # Copy compartments
    new_model.compartments = model.compartments.copy()
    new_model.repair(rebuild_index=True, rebuild_relationships=True)
    return new_model


# Internal
def _convert_metabolite(metabolite, to_class, from_class, prefix):
    """Convert a metabolite of class 'to_class' from class 'from_class'.

    Warnings
    --------
    This method is intended for internal use only. To safely convert between
    cobra and mass objects, use the corresponding conversion method.
    """
    if not isinstance(metabolite, from_class):
        raise TypeError("metabolite must be of {0}"
                        .format(from_class))

    id_str = _prefix_id_str(metabolite, prefix)

    new_metabolite = to_class(id=id_str, name=metabolite.name,
                              formula=metabolite.formula,
                              charge=metabolite.charge,
                              compartment=metabolite.compartment)
    new_metabolite._constraint_sense = metabolite._constraint_sense
    new_metabolite._bound = metabolite._bound

    return new_metabolite


def _convert_reaction(reaction, to_class, from_class, prefix):
    """Convert a reaction of class 'to_class' from class 'from_class'.

    Warnings
    --------
    This method is intended for internal use only. To safely convert between
    cobra and mass objects, use the corresponding conversion method.
    """
    if not isinstance(reaction, from_class):
        raise TypeError("reaction must be of {0}"
                        .format(from_class))

    id_str = _prefix_id_str(reaction, prefix)

    new_reaction = to_class(id=id_str, name=reaction.name,
                            subsystem=reaction.subsystem)

    new_reaction.lower_bound = reaction.lower_bound
    new_reaction.upper_bound = reaction.upper_bound
    new_reaction.variable_kind = reaction.variable_kind

    # Generate and add new genes
    for gene in reaction._genes:
        new_reaction._genes.add(gene.copy())
    # Add the gene reaction rule
    new_reaction._gene_reaction_rule = reaction._gene_reaction_rule

    return new_reaction


def _convert_model(model, to_class, from_class, prefix):
    """Convert a model of class 'to_class' from class 'from_class'.

    Warnings
    --------
    This method is intended for internal use only. To safely convert between
    cobra and mass objects, use the corresponding conversion method.
    """
    if not isinstance(model, from_class):
        raise TypeError("model must be of {0}"
                        .format(from_class))

    id_str = _prefix_id_str(model, prefix)

    new_model = to_class(id_or_model=id_str, name=model.name)
    return new_model


def _prefix_id_str(item, prefix):
    """Add the prefix to the ID string and return the new string.

    Warnings
    --------
    This method is intended for internal use only.
    """
    if prefix is not None:
        if not isinstance(prefix, string_types):
            raise TypeError("prefix must be a str")
        else:
            id_str = "{0}_{1}".format(prefix, item.id)
    else:
        id_str = item.id

    return id_str
