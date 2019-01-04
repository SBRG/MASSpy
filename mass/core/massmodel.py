# -*- coding: utf-8 -*-
"""TODO Module Docstrings."""
from __future__ import absolute_import

import logging
import re
from copy import copy, deepcopy
from functools import partial
from warnings import warn

import numpy as np

from six import integer_types, iteritems, iterkeys, itervalues, string_types

import sympy as sym

from cobra.core.dictlist import DictList
from cobra.core.object import Object
from cobra.util.context import HistoryManager, get_context

from mass.core.massmetabolite import MassMetabolite
from mass.core.massreaction import MassReaction
from mass.util import expressions
from mass.util.util import (_get_matrix_constructor, convert_matrix,
                            ensure_iterable, strip_time)


# Set the logger
LOGGER = logging.getLogger(__name__)
# Global
CHOPNSQ = ['C', 'H', 'O', 'P', 'N', 'S', 'q']
# Pre-compiled regular expressions for building reactions from strings
_RXN_ID_RE = re.compile("^(\w+):")
_MET_ID_RE = re.compile("^s\[(\S+)[,|\]]")
_REVERSIBLE_ARROW_RE = re.compile("<(-+|=+)>")
_FORWARD_ARROW_RE = re.compile("(-+|=+)>")
_REVERSE_ARROW_RE = re.compile("<(-+|=+)")
_NAME_ARG_RE = re.compile("name=(\w+)")
_FORMULA_ARG_RE = re.compile("formula=(\w+)")
_CHARGE_ARG_RE = re.compile("charge=(\w+)")
_COMPARTMENT_RE = re.compile("\](\[[A-Za-z]\])")
_EQUALS_RE = re.compile("=")


class MassModel(Object):
    """Class representation of a model.

    Parameters
    ----------
    id_or_model: str, mass.MassModel
        Either an identifier to associate with the MassModel given as a string,
        or an existing MassModel object. If an existing MassModel object is
        provided, a new MassModel object is instantiated with the same
        properties as the original MassModel.
    name: str, optional
        A human readable name for the model.
    matrix_type: {'dense', 'dok', 'lil', 'DataFrame', 'symbolic'}, optional
        A string identifiying the desired format for the stoichiometric matrix
        of the model. Types can include 'dense' for a standard numpy.array,
        'dok' or 'lil' to obtain the corresponding scipy.sparse matrix,
        'DataFrame' for a pandas.DataFrame, and 'symbolic' for a
        sympy.MutableDenseMatrix. For all matrix types, species (excluding
        genes) are row indicies and reactions are column indicies. If None,
        defaults to "dense".
    dtype: data-type, optional
        The desired array data-type for the stoichiometric matrix. If None,
        defaults to np.float64.

    Attributes
    ----------
    description: str
        A human-readable description of the model.
    reactions: cobra.DictList
        A cobra.DictList where keys are the reaction identifiers and the
        values are the associated MassReaction objects.
    metabolites: cobra.DictList
        A cobra.DictList where the keys are the metabolite identifiers and the
        values are the associated MassMetabolite objects.
    genes: cobra.DictList
        A cobra.DictList where the keys are the gene identifiers and the
        values are the associated cobra.Gene objects.
    initial_conditions: dict
        A dictionary to store the initial conditions of the metabolites,
        where keys are the MassMetabolites and values are initial conditions.
    custom_rates: dict
        A dictionary to store custom rate expressions for specific reactions,
        where keys are the MassReaction objects and values are the custom rate
        expressions given as sympy objects. Custom rate expressions will always
        be prioritized over automatically generated mass action rate laws.
    custom_parameters: dict
        A dictionary to store the custom parameters for the custom rates,
        where key:value pairs are the string identifiers of the parameters and
        their numerical value. Custom rate expressions will always be
        prioritized over automatically generated mass action rate laws.
    fixed_concentrations: dict
        A dictionary to store fixed metabolite concentrations, where keys can
        be MassMetabolite objects and/or string identifiers for 'external
        metabolites' of exchange reactions, and values are the fixed
        concentrations. Fixed concentrations will always be prioritized over
        time dependent metabolite concentrations and are treated as constants
        in other expressions.
    modules: set
        A set containing the identifiers of the MassModels (model.id)
        merged together to create the MassModel (self).
    compartments: dict
        A dictionary to store the compartment shorthands and their full names.
        Keys are the shorthands while values are the full names.
        Example: {'c': 'cytosol'}
    units: dict
        A dictionary to store the units used in the model for referencing.
        Example: {'N': 'Millimoles', 'Vol': 'Liters', 'Time': 'Hours'}

    Warnings
    --------
    MassModels can have different initial conditions from the ones stored in
        the MassMetabolites for various purposes. However, simulations will
        always utilize the initial conditions stored in the model, which
        can be accessed by the MassModel.initial_conditions attribute.
    The MassModel will not automatically track or convert units. Therefore, it
        is up to the user to ensure unit consistency in the model. The
        MassModel.units attribute is provided as a way to inform the user or
        others what units the model is currently using.

    """

    def __init__(self, id_or_model=None, name=None, matrix_type="dense",
                 dtype=np.float64):
        """Initialize the MassModel Object."""
        # Instiantiate a new MassModel object if a MassModel is given.
        if isinstance(id_or_model, MassModel):
            Object.__init__(self, id_or_model, name=name)
            self.__setstate__(id_or_model.__dict__)
            if not hasattr(self, "name"):
                self.name = None
            self.update_S(matrix_type=matrix_type, dtype=dtype)
        else:
            Object.__init__(self, id_or_model, name=name)
            self.description = ''
            # Initialize DictLists for reactions, metabolites, and genes
            self.reactions = DictList()
            self.metabolites = DictList()
            self.genes = DictList()
            # Initialize dictionaries for initial conditions, custom rates,
            # custom parameters, fixed concentrations, compartments, and units.
            self.initial_conditions = {}
            self.fixed_concentrations = {}
            self.custom_rates = {}
            self.custom_parameters = {}
            self._compartments = {}
            self._units = {}
            # Initialize a set to store the modules
            self.modules = set()
            # Store the stoichiometric matrix, its matrix type, and data type
            self._matrix_type = matrix_type
            self._dtype = dtype
            self._S = self._mk_stoich_matrix(matrix_type=self._matrix_type,
                                             dtype=self._dtype,
                                             update_model=True)
        # Store the HistoryManager contexts for context management
        self._contexts = []

    # Public
    @property
    def stoichiometric_matrix(self):
        """Return the stoichiometric matrix of the MassModel."""
        return self.update_S(matrix_type=self._matrix_type, dtype=self._dtype,
                             update_model=False)

    @property
    def rates(self):
        """Return a dict of reaction rate laws as sympy expressions.

        If a reaction has an associated custom rate expression, the custom rate
        will be prioritized and returned in the dictionary instead of the
        automatically generated mass action rate law expression.
        """
        return self.get_rate_laws(self.reactions, rate_type=0, sympy_expr=True,
                                  update_reactions=True)

    @property
    def ordinary_differential_equations(self):
        """Return a dict of ODEs for the metabolites as sympy expressions.

        MassMetabolites with fixed concentrations are considered constant and
        therefore not included in the returned dictionary.
        """
        return {met: met.ode for met in self.metabolites
                if met not in self.fixed_concentrations}

    @property
    def exchanges(self):
        """Return a list of exchange reactions in the model."""
        return [rxn for rxn in self.reactions if rxn.exchange]

    @property
    def external_metabolites(self):
        """Return a sorted list of all 'external' metabolites in the model."""
        return sorted(list(set(rxn.external_metabolite
                               for rxn in self.reactions if rxn.exchange)))

    @property
    def compartments(self):
        """Return a set of all metabolite compartments."""
        return {met.compartment: self._compartments.get(met.compartment, '')
                for met in self.metabolites if met.compartment is not None}

    @compartments.setter
    def compartments(self, value):
        """Get or set the dictionary of current compartment descriptions.

        Assigning a dictionary to this property updates the model's
        dictionary of compartment descriptions with the new values.

        Parameters
        ----------
        value : dict
            Dictionary mapping compartments abbreviations to full names.
            An empty dictionary will reset the compartments.

        """
        if value:
            self._compartments.update(value)
        else:
            setattr(self, "_compartments", {})

    @property
    def units(self):
        """Return a dictionary of stored model units."""
        return getattr(self, "_units")

    @units.setter
    def units(self, value):
        """Set the dictionary of current unit descriptions.

        Assigning a dictionary to this property updates the model's
        dictionary of unit descriptions with the new values.

        Parameters
        ----------
        value : dict
            Dictionary mapping unit abbreviations to full names.
            An empty dictionary will reset the unit.

        """
        if value:
            self._units.update(value)
        else:
            setattr(self, "_units", {})

    @property
    def irreversible_reactions(self):
        """Return a list of all irreversible reactions in the model."""
        return [rxn for rxn in self.reactions if not rxn.reversible]

    @property
    def steady_state_fluxes(self):
        """Return a dict of all reactions' steady state fluxes."""
        return {rxn: rxn.steady_state_flux for rxn in self.reactions
                if rxn.steady_state_flux is not None}

    @property
    def parameters(self):
        """Return all parameters associateed with the MassModel."""
        parameters = {}
        # Sort rate and equilibrium constants into seperate dictionaries,
        # then add dictionaries to returned parameter dictionary.
        for p_type in ["kf", "Keq", "kr"]:
            p_type_dict = {}
            for rxn in self.reactions:
                p_sym = getattr(rxn, p_type + "_str")
                if p_sym in rxn.parameters:
                    p_type_dict.update({p_sym: rxn.parameters[p_sym]})
            parameters.update({p_type: p_type_dict})
        # Add custom parameters and fixed concentrations.
        parameters.update({"Custom": self.custom_parameters})
        parameters.update({"Fixed": self.fixed_concentrations})

        return parameters

    # Shorthands
    @property
    def S(self):
        """Shorthand method to get the stoichiometric matrix."""
        return self.stoichiometric_matrix

    @property
    def ssfluxes(self):
        """Shorthand method to get all reactions' steady state fluxes."""
        return self.steady_state_fluxes

    @property
    def odes(self):
        """Shorthand method to get ODEs for the metabolites."""
        return self.ordinary_differential_equations

    def update_S(self, reaction_list=None, matrix_type=None, dtype=None,
                 update_model=True):
        """Update the stoichiometric matrix of the model.

        Parameters
        ----------
        reaction_list: list of mass.MassReactions, optional
            A list of MassReactions to be add to the stoichiometric matrix.
            Reactions must already exist in the model in order to update.
            If None, the entire stoichiometric matrix is reconstructed.
        matrix_type: {'dense', 'dok', 'lil', 'DataFrame', 'symbolic'}, optional
            A string identifiying the desired format for the stoichiometric
            matrix of the model. Types can include 'dense' for a standard
            numpy.array, 'dok' or 'lil' to obtain the corresponding
            scipy.sparse matrix, 'DataFrame' for a pandas.DataFrame, and
            'symbolic' for a sympy.MutableDenseMatrix. For all matrix types,
            species (excluding  genes) are the row indicies and reactions are
            the column indicies. If None, defaults to "dense".
        dtype: data-type, optional
            The desired array data-type for the stoichiometric matrix. If None,
            defaults to np.float64.
        update_model: bool, optional
            If True, will update the stored stoichiometric matrix, the matrix
            type, and the data-type for the model.

        Returns
        -------
        stoich_mat: matrix of class 'dtype'
            The stoichiometric matrix for the MassModel returned as the given
            matrix_type and with a data-type of 'dtype'.

        Notes
        -----
        reaction_list is assumed to be at the end of self.reactions.

        """
        # Use the model's stored matrix type if the matrix-type is not given.
        if matrix_type is None:
            matrix_type = self._matrix_type
        # Use the model's stored data-type if the data-type is not specified.
        if dtype is None:
            dtype = self._dtype

        # Check input of update_model
        if not isinstance(update_model, bool):
            raise TypeError("update_model must be a bool.")

        # If a matrix has not been constructed yet, or if there are no changes
        # to the reactions, return a newly constructed stoichiometric matrix.
        if self._S is None or reaction_list is None:
            stoich_mat = self._mk_stoich_matrix(matrix_type=matrix_type,
                                                dtype=dtype,
                                                update_model=update_model)
        # Update the current matrix if changes were made to the reactions.
        else:
            stoich_mat = self._update_stoichiometry(reaction_list,
                                                    matrix_type=matrix_type,
                                                    dtype=dtype)
        # Internally update the model if desired
        if update_model:
            self._update_model_s(stoich_mat, matrix_type=matrix_type,
                                 dtype=dtype)
        return stoich_mat

    def add_metabolites(self, metabolite_list, add_initial_conditons=False):
        """Add a list of metabolites to the MassModel.

        The change is reverted upon exit when using the MassModel as a context.

        Parameters
        ----------
        metabolite_list: list of mass.MassMetabolites
            A list of MassMetabolites to add to the MassModel.
        add_initial_conditons: bool, optional
            If True, the initial conditions associated with each metabolite are
            also added to the model. Otherwise, only the metabolites are added
            without their initial conditions.

        """
        # Ensure list is iterable.
        metabolite_list = ensure_iterable(metabolite_list)

        # Check whether a metabolite is a MassMetabolite object, then check if
        # the metabolite already exists in the model, ignoring those that do.
        for met in metabolite_list:
            if not isinstance(met, MassMetabolite):
                warn("Skipping {0}, not a MassMetabolite".format(met))
                metabolite_list.remove(met)

        existing_mets = [met for met in metabolite_list
                         if met.id in self.metabolites]
        for met in existing_mets:
            LOGGER.info("Skipped MassMetabolite %s as it already exists in"
                        " the MassModel.", met.id)
        metabolite_list = [met for met in metabolite_list
                           if met not in self.metabolites]

        # Have metabolites point to the model object.
        for met in metabolite_list:
            met._model = self
        # Add metabolites to the model
        self.metabolites += metabolite_list
        # Add the initial conditions if desired.
        if add_initial_conditons:
            self.set_initial_conditions(metabolite_list)

        context = get_context(self)
        if context:
            context(partial(self.metabolites.__isub__, metabolite_list))
            context(partial(setattr, met, "_model", None)
                    for met in metabolite_list)

    def remove_metabolites(self, metabolite_list, destructive=False):
        """Remove a list of metabolites from the MassModel.

        The metabolite's initial condition will also be removed from the model.

        The change is reverted upon exit when using the MassModel as a context.

        Parameters
        ----------
        metabolite_list: list of mass.MassMetabolites
            A list of MassMetabolites to add to the MassModel.
        destructive: bool, optional
            If False, the MassMetabolite is removed from all associated
            MassReactions. If True, also remove associated MassReactions from
            the MassModel.

        """
        # Ensure list is iterable.
        metabolite_list = ensure_iterable(metabolite_list)

        # Check whether a metabolite already exists in the model,
        # ignoring those that do not.
        metabolite_list = [met for met in metabolite_list
                           if met in self.metabolites]

        # Remove association to model
        for met in metabolite_list:
            met._model = None
            if not destructive:
                # Remove metabolites from reactions
                for rxn in list(met._reaction):
                    coeff = rxn._metabolites[met]
                    rxn.subtract_metabolites({met: coeff})
            else:
                # Remove associated reactions if destructive
                for rxn in list(met._reaction):
                    rxn.remove_from_model()
        # Remove initial conditions and then the metabolites
        self.remove_initial_conditions(metabolite_list)
        self.metabolites -= metabolite_list

        context = get_context(self)
        if context:
            context(partial(self.metabolites.__iadd__, metabolite_list))
            context(partial(setattr, met, "_model", self)
                    for met in metabolite_list)

    def set_initial_conditions(self, metabolite_list=None):
        """Set the initial conditions for a list of metabolites in the model.

        The change is reverted upon exit when using the MassModel as a context.

        Parameters
        ----------
        metabolite_list: list of mass.MassMetabolites
            A list of MassMetabolites to add to the MassModel. If None

        Notes
        -----
        The metabolite(s) must already exist in the model to set the initial
        conditions. Initial conditions for the metabolites are accessed through
        MassMetabolite.initial_condition. If an initial condition for a
        metabolite already exists in the model, it will be replaced.

        """
        # Select all metabolites if None provided.
        if metabolite_list is None:
            metabolite_list = self.metabolites
        metabolite_list = ensure_iterable(metabolite_list)

        # Check whether a metabolite already exists in the model,
        # ignoring those that do not.
        metabolite_list = [met for met in metabolite_list
                           if met in self.metabolites]
        # Add the initial conditions
        self.update_initial_conditions({met: met.ic
                                        for met in metabolite_list})

    def remove_initial_conditions(self, metabolite_list):
        """Remove initial conditions for a list of metabolites in the model.

        The change is reverted upon exit when using the MassModel as a context.

        Parameters
        ----------
        metabolite_list: list
            A list of MassMetabolite objects.

        """
        # Ensure list is iterable.
        metabolite_list = ensure_iterable(metabolite_list)

        # Check whether a metabolite already exists in the model,
        # ignoring those that do not.
        metabolite_list = [met for met in metabolite_list
                           if met in self.metabolites]

        # Keep track of existing initial conditions for context if needed.
        context = get_context(self)
        if context:
            existing_ics = {met: self.initial_conditions[met]
                            for met in metabolite_list
                            if met in self.initial_conditions}
        # Remove the initial conditions
        for met in metabolite_list:
            del self.initial_conditions[met]

        if context:
            context(partial(self.initial_conditions.update, existing_ics))

    def update_initial_conditions(self, initial_conditions,
                                  update_metabolites=False):
        """Update the initial conditions of the MassModel.

        The change is reverted upon exit when using the MassModel as a context.

        Parameters
        ----------
        initial_conditions: dict
            A dictionary where MassMetabolites are the keys and the initial
            conditions are the values.
        update_metabolites: bool, optional
            If True, will update the initial conditions in the MassMetabolite
            objects as well. Otherwise, only update the model.

        Notes
        -----
        The metabolite(s) must already exist in the model to set the initial
        conditions. Initial conditions for the metabolites are accessed through
        MassMetabolite.initial_condition. If an initial condition for a
        metabolite already exists in the model, it will be replaced.

        """
        if not isinstance(initial_conditions, dict):
            raise TypeError("initial_conditions must be a dictionary.")

        # Check whether a metabolite already exists in the model,
        # ignoring those that do not.
        initial_conditions = {met: ic
                              for met, ic in iteritems(initial_conditions)
                              if met in self.metabolites and ic is not None}
        # Keep track of existing initial conditions for context if needed.
        context = get_context(self)
        if context:
            existing_ics = {met: self.initial_conditions[met]
                            for met in initial_conditions
                            if met in self.initial_conditions}
        # Update initial conditions
        self.initial_conditions.update(initial_conditions)
        if update_metabolites:
            # Update the initial condition stored in the MassMetabolite.
            for met, ic in iteritems(initial_conditions):
                met.initial_condition = ic

        if context:
            context(partial(self.initial_conditions.pop, met)
                    for met in iterkeys(initial_conditions)
                    if met not in iterkeys(existing_ics))
            context(partial(self.initial_conditions.update, existing_ics))
            if update_metabolites:
                context(partial(setattr(met, "_initial_condition",
                                        existing_ics[met]))
                        for met in iterkeys(initial_conditions))

    def add_reactions(self, reaction_list):
        """Add MassReactions to the MassModel.

        MassReaction objects with identifiers identical to an existing reaction
        are ignored.

        The change is reverted upon exit when using the MassModel as a context

        Parameters
        ----------
        reaction_list: list of mass.MassReactions
            A list of MassReaction objects.

        """
        def existing_filter(rxn):
            if rxn.id in self.reactions:
                LOGGER.warning(
                    "Ignoring reaction '%s' since it already exists.", rxn.id)
                return False
            return True

        # Ensure list is iterable.
        reaction_list = ensure_iterable(reaction_list)

        # Check whether a metabolite is a MassReaction object, then check if
        # the reaction already exists in the model, ignoring those that do.
        for rxn in reaction_list:
            if not isinstance(rxn, MassReaction):
                warn("Skipping {0}, not a MassReaction object".format(rxn))
                reaction_list.remove(rxn)

        # Check whether reaxctions exist in the model.
        pruned = DictList(filter(existing_filter, reaction_list))
        context = get_context(self)
        # Add reactions, and have reactions point to the model.
        for rxn in pruned:
            rxn._model = self
            # Loop through metabolites associated with a reaction.
            for met in list(rxn._metabolites):
                # If a metabolite doesn't exist in the model,
                # add it with its associated initial condition.
                if met not in self.metabolites:
                    self.add_metabolites(met, add_initial_conditons=True)
                # Otherwise have reactions point to the metabolite in the model
                else:
                    coeff = rxn._metabolites.pop(met)
                    met = self.metabolites.get_by_id(met.id)
                    rxn._metabolites[met] = coeff
                met._reaction.add(rxn)
                if context:
                    context(partial(met._reaction.remove, rxn))

            # Loop through genes associated with a reaction.
            for gene in list(rxn._genes):
                # If the gene is not in the model, add it and point to model.
                if not self.genes.has_id(gene.id):
                    self.genes += [gene]
                    gene._model = self
                    if context:
                        context(partial(self.genes.__isub__, [gene]))
                        context(partial(setattr, gene, "_model", None))
                # Otherwise, make the gene point to the one in the model.
                else:
                    model_gene = self.genes.get_by_id(gene.id)
                    if model_gene is not gene:
                        rxn._dissociate_gene(gene)
                        rxn._associate_gene(model_gene)
        # Add reactions to the model
        self.reactions += pruned

        if context:
            context(partial(self.reactions.__isub__, pruned))
            for rxn in reaction_list:
                context(partial(setattr, rxn, "_model", None))

    def remove_reactions(self, reaction_list, remove_orphans=False):
        """Remove Massreactions from the MassModel.

        The change is reverted upon exit when using the MassModel as a context.

        Parameters
        ----------
        reaction_list: list of mass.MassReactions
            A list of MassReaction objects to be removed from the model.
        remove_orphans: bool, optional
            If True, will also remove orphaned genes and MassMetabolites from
            the MassModel.

        """
        # Ensure list is iterable.
        reaction_list = ensure_iterable(reaction_list)

        # Check whether a reaction already exists in the model,
        # ignoring those that do not.
        reaction_list = [rxn for rxn in reaction_list
                         if rxn.id in self.reactions]

        context = get_context(self)
        # Remove model association
        for rxn in reaction_list:
            rxn._model = None
            # Remove metabolite association
            for met in rxn._metabolites:
                if rxn in met._reaction:
                    met._reaction.remove(rxn)
                    # Remove orphaned metabolites and their initial conditions.
                    if remove_orphans and not met._reaction:
                        self.remove_metabolites(met)
                    if context:
                        context(partial(met._reaction.add, rxn))
            # Remove gene association
            for gene in rxn._genes:
                if rxn in gene._reaction:
                    gene._reaction.remove(rxn)
                    # Remove orphaned genes
                    if remove_orphans and not gene._reaction:
                        self.genes.remove(gene)
                        if context:
                            context(partial(self.genes.add, gene))
                    if context:
                        context(partial(gene._reaction.add, rxn))
            if context:
                context(partial(setattr, rxn, "_model", self))
        # Remove reactions from the model
        self.reactions -= reaction_list

        if context:
            context(partial(self.reactions.__iadd__, reaction_list))
            context(partial(setattr, rxn, "_model", self)
                    for rxn in reaction_list)

    def add_exchange(self, metabolite, exchange_type="exchange",
                     external_concentration=None):
        """Add a pre-defined exchange reaction for a given metabolite.

        Pre-defined exchange types can be "exchange" for reversibly entering or
        exiting the compartment, "source" for irreversibly entering the
        compartment, and "demand" for irreversibly exiting the compartment.

        The change is reverted upon exit when using the MassModel as a context.

        Parameters
        ----------
        metabolite: MassMetabolite
            The metabolite involved in the exchange reaction.
        exchange_type: str {"demand", "exchange", "source"}, optional
            The type of exchange reaction to create.
        external_concentration: float, optional
            The fixed concentration value to set for the the external species.

        Returns
        -------
        exchange_rxn: mass.MassReaction
            The MassReaction object of the new exchange reaction. If the
            reaction already exists, the existing MassReaction object is
            returned.

        """
        # Check whether a metabolite is a MassMetabolite object:
        if not isinstance(metabolite, MassMetabolite):
            raise TypeError("metabolite must be a MassMetabolite object")

        type_dict = {"demand": ["DM", -1, False],
                     "source": ["S", 1, False],
                     "exchange": ["EX", -1, True]}
        # Set the type of exchange
        if exchange_type not in type_dict:
            raise ValueError("Exchange type must be either "
                             "'exchange', 'source',  or 'sink'")

        else:
            values = type_dict[exchange_type]
            rxn_id = metabolite.id
            # Remove leading underscore if necessary
            if re.match("\_", rxn_id[0]):
                rxn_id = rxn_id[1:]
            # Find compartment and replace with "_e" for external/extracellular
            _c = re.search("^\w*\S(?!<\_)(\_\S+)$", metabolite.id)
            if _c is not None and not re.search("\_L$|\_D$", _c.group(1)):
                rxn_id = re.sub(_c.group(1), "_e", rxn_id, count=1)
            rxn_id = "{}_{}".format(values[0], rxn_id)
            if rxn_id in self.reactions:
                warn("Reaction {0} already exists in model".format(rxn_id))
                return self.reactions.get_by_id(rxn_id)
            c = values[1]
            reversible = values[2]

        rxn_name = "{} {}".format(metabolite.name, exchange_type)
        exchange_rxn = MassReaction(id=rxn_id, name=rxn_name,
                                    subsystem="Transport/Exchange",
                                    reversible=reversible)
        exchange_rxn.add_metabolites({metabolite: c})
        self.add_reactions([exchange_rxn])
        self.add_fixed_concentrations({exchange_rxn.external_metabolite:
                                       external_concentration})

        return exchange_rxn

    def get_rate_laws(self, reaction_list=None, rate_type=0, sympy_expr=True,
                      update_reactions=False):
        """Get the rate laws for a list of reactions in the model.

        Parameters
        ----------
        reaction_list: list, optional
            A list of MassReactions to get the rate laws for. Reactions must
            already exist in the model. If None, then return the rates for all
            reactions in the model.
        rate_type: int {1, 2, 3}, optional
            The type of rate law to display. Must be 0, 1, 2, or 3. If 0,
            the type of rate law is determined by what is currently set in the
            reaction. Default is 0.
            Type 1 will utilize the forward rate and equilibrium constants.
            Type 2 will utilize the forward rate and reverse rate constants.
            Type 3 will utilize the equilibrium and reverse rate constants.
        symp_expr: bool, optional
            If True, then return the rate law as a sympy expression, otherwise
            return the rate law as a human readable string.
        update_reactions: bool, optional
            If True, update the MassReaction default rate type in addition to
            returning the rate laws.

        Returns
        -------
        rate_dict: dictionary of reaction rates where keys are the reaction ids
            and values are the rate law expressions.

        """
        # Check the inputs
        if not isinstance(rate_type, (integer_types, float)):
            raise TypeError("rate_type must be an int or float")
        elif not isinstance(sympy_expr, bool):
            raise TypeError("sympy_expr must be a bool")
        else:
            rate_type = int(rate_type)
            if rate_type not in {0, 1, 2, 3}:
                raise ValueError("rate_type must be 0, 1, 2, or 3")

        # Use the MassModel reactions if no reaction list is given
        if reaction_list is None:
            reaction_list = self.reactions
        # Ensure list is iterable.
        reaction_list = ensure_iterable(reaction_list)

        if rate_type == 0:
            rate_dict = {rxn: rxn.get_rate_law(rxn._rtype, sympy_expr,
                                               update_reactions)
                         for rxn in reaction_list}

        else:
            rate_dict = {rxn: rxn.get_rate_law(rate_type, sympy_expr,
                                               update_reactions)
                         for rxn in reaction_list}

        if not self.custom_rates:
            rate_dict.update(self.custom_rates)

        return rate_dict

    def get_mass_action_ratios(self, reaction_list=None, sympy_expr=True):
        """Get the mass action ratios for a list of reactions in the model.

        Parameters
        ----------
        reaction_list: list, optional
            A list of MassReactions to get the mass action ratios for.
            Reactions must already exist in the model. If None, then return the
            rates for all reactions in the model.
        sympy_expr: bool, optional
            If True, then return the mass action ratios as a sympy expression,
            otherwise return the ratio as a human readable string.

        Returns
        -------
        ratio_dict: dictionary of mass action ratios where keys are the
            reaction ids and values are the ratios.

        """
        # Use the MassModel reactions if no reaction list is given
        if reaction_list is None:
            reaction_list = self.reactions
        # Ensure list is iterable.
        reaction_list = ensure_iterable(reaction_list)

        ratio_dict = {rxn: rxn.get_mass_action_ratio(sympy_expr)
                      for rxn in reaction_list}

        return ratio_dict

    def get_disquilibrium_ratios(self, reaction_list=None, sympy_expr=True):
        """Get the disequilibrium ratios for a list of reactions in the model.

        Parameters
        ----------
        reaction_list: list, optional
            A list of MassReactions to get the disequilibrium ratios for.
            Reactions must already exist in the model. If None, then return the
            rates for all reactions in the model.
        sympy_expr: bool, optional
            If True, then return the disequilibrium ratios as a sympy
            expression, otherwise return the ratio as a human readable string.

        Returns
        -------
        ratio_dict: dictionary of disequilibrium ratios where keys are the
            reaction ids and values are the ratios.

        """
        # Use the MassModel reactions if no reaction list is given
        if reaction_list is None:
            reaction_list = self.reactions
        # Ensure list is iterable.
        reaction_list = ensure_iterable(reaction_list)

        ratio_dict = {rxn: rxn.get_disequilibrium_ratio(sympy_expr)
                      for rxn in reaction_list}

        return ratio_dict

    def add_custom_rate(self, reaction, custom_rate, custom_parameters=None):
        """Add a custom rate for a reaction to the model.

        The change is reverted upon exit when using the MassModel as a context.

        Parameters
        ----------
        reaction: mass.MassReaction
            The MassReaction associated with the custom rate.
        custom_rate: str
            The string representation of the custom rate expression. The string
            representation of the custom rate will be used to create a sympy
            expression that represents the custom rate.
        custom_parameters: dict, optional
            A dictionary of custom parameters for the custom rate where the
            key:value pairs are the strings representing the custom parameters
            and their numerical values. The string representation of the custom
            parametes will be used to create the symbols needed for the sympy
            expression of the custom rate. If None, then parameters are assumed
            to be a one of the MassReaction.

        Notes
        -----
        Metabolites must already exist in the MassReaction. However, the
            default parameters of a MassReaction (kf_RID, Keq_RID, kr_RID) are
            automatically taken into account and do not need to be defined as
            an additional custom parameter.

        """
        if custom_parameters is not None:
            custom_parameter_list = list(iterkeys(custom_parameters))
        else:
            custom_parameters = {}
            custom_parameter_list = []
        # Use any existing custom parameters if they are in the rate law.
        existing_customs = self.custom_parameters
        if existing_customs:
            for custom_parameter in iterkeys(existing_customs):
                if re.search(custom_parameter, custom_rate) and \
                   custom_parameter not in custom_parameter_list:
                    custom_parameter_list.append(custom_parameter)

        custom_rate = expressions.create_custom_rate(reaction, custom_rate,
                                                     custom_parameter_list)
        self.custom_rates.update({reaction: custom_rate})
        self.custom_parameters.update(custom_parameters)

        context = get_context(self)
        if context:
            context(partial(self.custom_rates.pop, reaction))
            context(partial((self.custom_parameters.pop, key)
                            for key in custom_parameter_list
                            if key in iterkeys(self.custom_parameters)))
            context(partial(self.custom_parameters.update, existing_customs))

    def remove_custom_rate(self, reaction, remove_orphans=True):
        """Remove the custom rate for a given reaction from the model.

        The change is reverted upon exit when using the MassModel as a context.

        Parameters
        ----------
        reaction: mass.MassReaction
            The MassReaction assoicated with the custom rate to be removed.
        remove_orphans: bool, optional
            If True, then remove any orphaned custom parameters from the model.

        """
        try:
            rate_to_remove = self.custom_rates[reaction]
        except KeyError:
            warn("Did not find a custom custom rate expression associated "
                 "with reaction {0}.".format(reaction.id))
            return
        # Remove the rate
        del self.custom_rates[reaction]

        # Remove orphaned custom parameters if desired.
        symbols = rate_to_remove.atoms(sym.Symbol)

        # Save currently existing parameters for context management if needed.
        existing = {str(sym): self.custom_parameters[str(symbol)]
                    for symbol in symbols}

        if remove_orphans and self.custom_rates:
            # Determine symbols still in use.
            other_symbols = set()
            for custom_rate in itervalues(self.custom_rates):
                other_symbols.update(custom_rate.atoms(sym.Symbol))

            # Remove those that are not being used in any custom rate.
            for symbol in other_symbols:
                if symbol in symbols.copy():
                    symbols.remove(symbol)

            for symbol in symbols:
                del self.custom_parameters[str(sym)]

        context = get_context(self)
        if context:
            context(partial(self.custom_rates.update,
                            {reaction: rate_to_remove}))
            if remove_orphans:
                context(partial(self.custom_parameters.update, existing))

    def reset_custom_rates(self):
        """Reset all custom rate expressions and parameters in a model.

        Warnings
        --------
        Using this method will remove all custom rates and custom rate
            parameters in the MassModel object. To remove a specific rate
            without affecting the other custom rates or custom parameters, use
            the MassModel.remove_custom_rate method instead.

        """
        self.custom_rates = {}
        self.custom_parameters = {}
        print("All custom rate expressions and parameters have been reset")

    def get_elemental_matrix(self, matrix_type=None, dtype=None):
        """Get the elemental matrix for a model.

        Parameters
        ----------
        matrix_type: {'dense', 'dok', 'lil', 'DataFrame', 'symbolic'}, optional
            A string identifiying the desired format for the elemental matrix
            of the model. Types can include 'dense' for a standard numpy.array,
            'dok' or 'lil' to obtain the corresponding scipy.sparse matrix,
            'DataFrame' for a pandas.DataFrame, and 'symbolic' for a
            sympy.MutableDenseMatrix. For all matrix types, the elements are
            row indicies and species (excluding genes) are column indicies.
            Default is 'dense'.
        dtype: data-type, optional
            The desired data-type for the array. Default is np.float64.

        Returns
        -------
        elem_mat: matrix of class 'dtype'
            The elemental matrix for the model.

        """
        # Set up for matrix construction if matrix types are correct.
        (matrix_constructor, matrix_type, dtype) = _get_matrix_constructor(
            matrix_type=matrix_type, dtype=dtype)
        # Build the elemental matrix
        elem_mat = matrix_constructor((len(CHOPNSQ), len(self.metabolites)))
        # Get indices for elements and metabolites
        e_ind = CHOPNSQ.index
        m_ind = self.metabolites.index

        # Fill the elemental matrix
        for met in self.metabolites:
            for element in CHOPNSQ:
                if element in iterkeys(met.elements):
                    amount = met.elements[element]
                elif re.match("q", element) and met.charge is not None:
                    amount = met.charge
                else:
                    amount = 0
                elem_mat[e_ind(element), m_ind(met)] = amount
        # Convert matrix to a dataframe if matrix type is a dataframe
        elem_mat = convert_matrix(elem_mat, matrix_type=matrix_type,
                                  dtype=dtype, row_ids=CHOPNSQ,
                                  col_ids=[m.id for m in self.metabolites])

        return elem_mat

    def get_elemental_charge_balancing(self, matrix_type=None, dtype=None):
        """Get the elemental charge balance as a matrix for a model.

        Parameters
        ----------
        matrix_type: {'dense', 'dok', 'lil', 'DataFrame', 'symbolic'}, optional
            A string identifiying the desired format for the elemental charge
            balancing matrix of the model. Types can include 'dense' for a
            standard numpy.array, 'dok' or 'lil' to obtain the corresponding
            scipy.sparse matrix, 'DataFrame' for a pandas.DataFrame, and
            'symbolic' for a sympy.MutableDenseMatrix. For all matrix types,
            the elements are row indicies and reactions are column indicies.'
            Default is 'dense'.
        dtype: data-type, optional
            The desired data-type for the array. Default is np.float64.

        Returns
        -------
        charge_mat: matrix of class 'dtype'
            The elemental charge balancing as a matrix for the model.

        """
        elem_mat = self.get_elemental_matrix()
        stoich_mat = self.update_S(matrix_type="dense", update_model=False)
        charge_mat = elem_mat.dot(stoich_mat)

        if matrix_type is None:
            matrix_type = "dense"
        if dtype is None:
            dtype = np.float64

        charge_mat = convert_matrix(charge_mat, matrix_type=matrix_type,
                                    dtype=dtype, row_ids=CHOPNSQ,
                                    col_ids=[r.id for r in self.reactions])
        return charge_mat

    def add_fixed_concentrations(self, fixed_concentrations):
        """Add fixed concentrations values for the given metabolites.

        A fixed metabolite concentration will remain at a constant value when
        simulating the MassModel.

        Parameters
        ----------
        fixed_concentrations: dict
            A dictionary of fixed concentrations where metabolites are the keys
            and fixed concentration value. The metabolite must already exist
            in the model, or it must be the string representation of the
            "external" metabolite in an exchange reaction. (e.g. 'MetabID_e')

        Notes
        -----
        Fixed concentrations always have priority over initial conditions.

        See Also
        --------
        MassModel.external_metabolites

        """
        if not isinstance(fixed_concentrations, dict):
            raise TypeError("fixed_concentrations must be a dict.")
        for met, fixed_conc in iteritems(fixed_concentrations):
            if met not in self.external_metabolites and \
               met not in self.metabolites:
                raise ValueError("Did not find {0} in model metabolites or in "
                                 "exchange reactions.".format(met))

            if not isinstance(fixed_conc, (integer_types, float)):
                raise TypeError("Fixed concentrations must be ints or floats.")
            elif fixed_conc < 0.:
                raise ValueError("Fixed concentrations must be non-negative.")
            else:
                fixed_conc = float(fixed_conc)
        # Keep track of existing concentrations for context management.
        context = get_context(self)
        if context:
            existing_concs = {met: self.fixed_concentrations[met]
                              for met in fixed_concentrations
                              if met in self.fixed_concentrations}

        self.fixed_concentrations.update(fixed_concentrations)

        if context:
            context(partial(self.fixed_concentrations.pop, key)
                    for key in iterkeys(fixed_concentrations)
                    if key not in existing_concs)
            context(partial(self.fixed_concentrations.update, existing_concs))

    def remove_fixed_concentrations(self, metabolite_list):
        """Remove a fixed concentration for a list of metabolites.

        Parameters
        ----------
        metabolite_list: list
            A list of metabolites to remove the fixed concentrations for.
            Metabolites must already exist in the model, or it must be the
            string representation of the "external" metabolite in an exchange
            reaction. (e.g. 'MetabID_Xt').

        See Also
        --------
        MassModel.external_metabolites

        """
        metabolite_list = ensure_iterable(metabolite_list)
        # Check whether a metabolite already exists in the model,
        # ignoring those that do not.
        metabolite_list = [met for met in metabolite_list
                           if met in self.fixed_concentrations]

        # Keep track of existing concentrations for context management.
        context = get_context(self)
        if context:
            existing_concs = {met: self.fixed_concentrations[met]
                              for met in metabolite_list
                              if met in self.fixed_concentrations}
        # Remove the initial conditions
        for met in metabolite_list:
            del self.fixed_concentrations[met]

        if context:
            context(partial(self.fixed_concentrations.update, existing_concs))

    def repair(self, rebuild_index=True, rebuild_relationships=True):
        """Update all indicies and pointers in the model.

        Parameters
        ----------
        rebuild_index: bool, optional
            If True, then rebuild the indicies kept in the reactions,
            metabolites, and genes.
        rebuild_relationships: bool, optional
            If True, then reset all associations between the reactions,
            metabolites genes, and the MassModel, and rebuilds them.

        """
        if not isinstance(rebuild_index, bool) or \
           not isinstance(rebuild_relationships, bool):
            raise TypeError("rebuild_index and rebuild_relationships "
                            "must be bools.")
        # Rebuild the DictList indices
        if rebuild_index:
            self.reactions._generate_index()
            self.metabolites._generate_index()
            self.genes._generate_index()
        # Rebuild the relationships between reactions and their associated
        # genes and metabolites
        if rebuild_relationships:
            for met in self.metabolites:
                met._reaction.clear()
            for gene in self.genes:
                gene._reaction.clear()
            for rxn in self.reactions:
                for met in rxn.metabolites:
                    met._reaction.add(rxn)
                for gene in rxn.genes:
                    gene._reaction.add(rxn)
        # Make all objects point to the model.
        for dictlist in (self.reactions, self.metabolites, self.genes):
            for item in dictlist:
                item._model = self

    def copy(self):
        """Create a partial "deepcopy" of the MassModel.

        All of the MassMetabolite, MassReaction, and Gene objects, the initial
        conditions, fixed concentrations, custom_rates, and the stoichiometric
        matrix are created anew, but in a faster fashion than deepcopy.
        """
        # Define a new model
        new_model = self.__class__()
        # Define items that will not be copied by their references
        do_not_copy_by_ref = {"metabolites", "reactions", "genes",
                              "initial_conditions", "_S", "custom_rates",
                              "custom_parameters", "notes", "annotation"}
        for attr in self.__dict__:
            if attr not in do_not_copy_by_ref:
                new_model.__dict__[attr] = self.__dict__[attr]
        new_model.notes = deepcopy(self.notes)
        new_model.annotation = deepcopy(self.annotation)

        # Copy the metabolites
        do_not_copy_by_ref = {"_reaction", "_model"}
        new_model.metabolites = DictList()
        for met in self.metabolites:
            new_met = met.__class__()
            for attr, value in iteritems(met.__dict__):
                if attr in do_not_copy_by_ref:
                    continue
                elif re.match("formula", attr):
                    new_met.__dict__[attr] = copy(value)
                else:
                    new_met.__dict__[attr] = value
            new_met._model = new_model
            new_model.metabolites.append(new_met)
            # Copy the initial condition
            if met in iterkeys(self.initial_conditions):
                ic = self.initial_conditions[met]
                new_model.initial_conditions[new_met] = ic

        # Copy the genes
        new_model.genes = DictList()
        for gene in self.genes:
            new_gene = gene.__class__(None)
            for attr, value in iteritems(gene.__dict__):
                if attr in do_not_copy_by_ref:
                    continue
                elif re.match("formula", attr):
                    new_met.__dict__[attr] = copy(value)
                else:
                    new_met.__dict__[attr] = value
            new_gene._model = new_model
            new_model.genes.append(new_gene)

        # Copy the reactions
        do_not_copy_by_ref = {"_model", "_metabolites", "_genes"}
        new_model.reactions = DictList()
        for rxn in self.reactions:
            new_rxn = rxn.__class__()
            for attr, value in iteritems(rxn.__dict__):
                if attr not in do_not_copy_by_ref:
                    new_rxn.__dict__[attr] = copy(value)
            new_rxn._model = new_model
            new_model.reactions.append(new_rxn)
            # Copy custom rates for the reaction
            if rxn in iterkeys(self.custom_rates):
                new_model.custom_rates[new_rxn] = self.custom_rates[rxn]
            # Update awareness
            for met, stoich in iteritems(rxn._metabolites):
                new_met = new_model.metabolites.get_by_id(met.id)
                new_rxn._metabolites[new_met] = stoich
                new_met._reaction.add(new_rxn)
            for gene in rxn._genes:
                new_gene = new_model.genes.get_by_id(gene.id)
                new_rxn._genes.add(new_gene)
                new_gene._reaction.add(new_rxn)
        # Copy custom parameters if there are custom rates
        if new_model.custom_rates:
            new_model.custom_parameters = copy(self.custom_parameters)

        # Create the new stoichiometric matrix for the model.
        new_model._S = self._mk_stoich_matrix(matrix_type=self._matrix_type,
                                              dtype=self._dtype,
                                              update_model=True)
        return new_model

    def merge(self, second_model, prefix_existing=None, inplace=False,
              new_model_id=None):
        """Merge two MassModels into one MassModel with the objects from both.

        The reactions, metabolites, genes, initial conditions, fixed
        concentrations, custom rate laws, and rate parameter will be added from
        the second model into the first model. In cases where identifiers for
        objects or relevent dictionaries are identical, priority will be given
        to what exists in the first model.

        The change is reverted upon exit when using the MassModel as a context.

        Parameters
        ----------
        second_model: mass.MassModel
            The second MassModel to merge into the first model.
        prefix_existing: str, optional
            If provided, the string is used to prefix the reaction identifier
            of a reaction in the second model if that reaction already exists
            within the first model.
        inplace: bool, optional
            If True, then add the contents of the second model directly into
            the first model. Otherwise, a new MassModel object is created, and
            the first model is left untouched.
        new_model_id: str, optional
            If provided, the string is used as the identifier for the merged
            model. If None and inplace is True, the model ID of the first model
            will be used. If None and inplace is False, a new combined ID
            will be used for the new MassModel object.

        Returns
        -------
        merged_model: mass.MassModel
            A new MassModel object or self representing the merged model.

        """
        # Set the merged model object and the ID of the merged model
        if inplace:
            merged_model = self
            if new_model_id is None:
                new_model_id = self.id
        else:
            merged_model = self.copy()
            if new_model_id is None:
                new_model_id = self.id + "_" + second_model.id
        merged_model.id = new_model_id

        # Add the reactions of the second model to the first model.
        new_reactions = deepcopy(second_model.reactions)
        if prefix_existing is not None:
            existing = new_reactions.query(lambda r: r.id in self.reactions)
            for rxn in existing:
                rxn.id = prefix_existing + "_" + rxn.id
        merged_model.add_reactions(new_reactions)
        merged_model.repair()

        # Add initial conditions
        existing = [m.id for m in iterkeys(merged_model.initial_conditions)]
        merged_model.update_initial_conditions({
            merged_model.metabolites.get_by_id(m.id): ic
            for m, ic in iteritems(second_model.initial_conditions)
            if m.id not in existing})

        # Add fixed concentrations
        existing = [m.id if isinstance(m, MassMetabolite) else m
                    for m in iterkeys(merged_model.fixed_concentrations)]

        merged_model.add_fixed_concentrations({
            m: fc for m, fc in iteritems(second_model.fixed_concentrations)
            if str(m) not in existing})

        # Add custom parameters
        existing = [cp for cp in iterkeys(merged_model.custom_parameters)]
        merged_model.custom_parameters.update(
            {cp: v for cp, v in iteritems(second_model.custom_parameters)
             if cp not in existing})
        # Add custom rates
        existing = [rxn.id for rxn in iterkeys(merged_model.custom_rates)]
        merged_model.custom_rates.update(
            {merged_model.reactions.get_by_id(rxn.id): cr
             for rxn, cr in iteritems(second_model.custom_rates)
             if rxn.id not in existing})
        # Add old models to the module set
        merged_model.modules.add(second_model.id)
        if not inplace:
            merged_model.modules.add(self.id)

        return merged_model

    def compute_steady_state_fluxes(self, pathways, independent_fluxes,
                                    update_reactions=False):
        """Calculate the unique steady state flux for each reaction.

        The unique steady state flux for each reaction in the MassModel is
        calculated using defined pathways, independently defined fluxes, and
        steady state concentrations.

        Parameters
        ----------
        pathways: array or array-like
            An array or array-like object that define the pathways through the
            reaction network of the MassModel. The given pathway vectors must
            be the same length as the number of reactions in the model, with
            indicies of values in the pathway vector corresponding to the
            indicies of reactions in the MassModel.reactions attribute.
        independent_fluxes: dict
            A dictionary of steady state fluxes where MassReactions are keys
            and fluxes are values to utilize in order to calculate all other
            steady state fluxes. Must be the same length as the number of
            specified pathways.
        update_reactions: bool, optional
            If True, then update the MassReaction.steady_state_flux attribute
            with the calculated steady state flux values.

        Return
        ------
        steady_state_fluxes: np.ndarray
            A numpy array of the calculated steady state fluxes. The indicies
            of the values in the pathway vector correspond to the indicies
            of the reactions in the MassModel.reactions attribute.

        Notes
        -----
        The number of individually defined fluxes must be the same as the
            number of pathways in order to determine the solution. For best
            results, the number of pathways to specify must equal the dimension
            of the right nullspace.

        """
        # Check inputs:
        if not isinstance(pathways, (np.ndarray, list)):
            raise TypeError("Pathways must be numpy.ndarrays or array-like, "
                            "such as a list of lists.")
        pathways = np.array(pathways)
        if len(self.reactions) != pathways.shape[1]:
            raise ValueError("Pathways must have the same number of columns as"
                             " the number of reactions in the model.")
        if not isinstance(independent_fluxes, dict):
            raise TypeError("independent_fluxes must be a dict")
        if not isinstance(update_reactions, bool):
            raise TypeError("update_reactions must be a bool")

        coeffs = []
        values = []
        for i, rxn in enumerate(self.reactions):
            if rxn in independent_fluxes:
                values.append(independent_fluxes[rxn])
                coeffs.append([path[i] for path in pathways])
        # Inverse coefficient matrix
        coeffs = np.linalg.inv(coeffs)

        # Obtain the inner product of values and coefficients,
        # then obtain the inner product of the pathways and first inner product
        steady_state_fluxes = np.inner(pathways.T, np.inner(coeffs, values))
        # Update the reactions if desired
        if update_reactions:
            for i, rxn in enumerate(self.reactions):
                rxn.steady_state_flux = steady_state_fluxes[i]

        return steady_state_fluxes

    def calculate_PERCs(self, steady_state_concentrations=None,
                        steady_state_fluxes=None,
                        at_equilibrium_default=100000, update_reactions=False):
        """Calculate pseudo-order rate constants for reactions in the model.

        Pseudo-order rate constants (PERCs) are forward rate constants, and are
        calculated based on the steady state concentrations and fluxes.

        Parameters
        ----------
        steady_state_concentrations: dict, optional
            A dictionary of steady state concentrations where MassMetabolites
            are keys and concentrations are the values. If None, the
            initial conditions and fixed concentrations that exist in the
            MassModel are used.
        steady_state_fluxes: dict, optional
            A dictionary of steady state fluxes where MassReactions are keys
            and fluxes are the values. If None, the steady state flux values
            stored in each of the reactions in the model are used.
        at_equilibrium_default: float, optional
            The value to set the pseudo-order rate constant if the reaction is
            at equilibrium. Default is 100,000.
        update_reactions: bool, optional
            If True, will update the values for the forward rate constants in
            the MassReactions with the calculated pseudo-order rate constants.

        Returns
        -------
        percs_dict: dict
            A dictionary where keys are strings identifers of the pseudo-order
            rate constants (kf_RID) and values are the calculated PERC value

        """
        # TODO Add QCQA check to ensure calculating percs cannot occur if
        # necessary values are not defined.
        if steady_state_concentrations is None:
            steady_state_concentrations = self.fixed_concentrations
            steady_state_concentrations.update(
                {m.id: ic for m, ic in iteritems(self.initial_conditions)})
        if steady_state_fluxes is None:
            steady_state_fluxes = self.steady_state_fluxes

        def calculate_sol(flux, rate_equation, perc, default):
            sol = sym.solveset(sym.Eq(flux, rate_equation),
                               perc, domain=sym.S.Reals)
            if isinstance(sol, type(sym.S.Reals)) or sol.is_EmptySet \
               or next(iter(sol)) <= 0:
                sol = float(default)
            else:
                sol = float(next(iter(sol)))

            return sol

        percs_dict = {}
        # Collect concentration and parameter values
        for rxn, rate in iteritems(strip_time(self.rates)):
            values = {}
            for symbol in list(rate.atoms(sym.Symbol)):
                if str(symbol) in steady_state_concentrations:
                    values[symbol] = steady_state_concentrations[str(symbol)]
                elif str(symbol) in [rxn.Keq_str, rxn.kf_str, rxn.kr_str]:
                    if str(symbol) == rxn.kf_str:
                        perc = symbol
                    else:
                        values[symbol] = rxn.parameters[str(symbol)]

            sol = calculate_sol(steady_state_fluxes[rxn], rate.subs(values),
                                perc, at_equilibrium_default)

            percs_dict.update({str(perc): sol})
            if update_reactions:
                rxn.kf = percs_dict[str(perc)]

        return percs_dict

    def string_to_mass(self, reaction_strings, term_split="+"):
        """Create MassReaction and MassMetabolite objects from strings.

        To correctly parse a string, it must be in the following format:
            "RID: s[ID, **kwargs] + s[ID, **kwargs] <=>  s[ID, **kwargs]
        where kwargs can be the metabolite attributes 'name', 'formula',
        'charge'. For example:
            "v1: s[x1, name=xOne, charge=2] <=> s[x2, formula=X]"

        To add a compartment for a species, add "[c]" where c is a letter
        representing the compartment for the species. For example:
            "v1: s[x1][c] <=> s[x2][c]"

        When creating bound enzyme forms, it is recommended to use '&' in the
        species ID to represent the bound enzyme-metabolite. For example:
            "E1: s[ENZYME][c] + s[metabolite][c] <=> s[ENZYME&metabolite][c]"

        Note that a reaction ID and a metabolite ID are always required.

        Parameters
        ----------
        reaction_strings: str, list of strs
            String or list of strings representing the reaction. Reversibility
            is inferred from the arrow, and metabolites in the model are used
            if they exist or are created if they do not.
        term_split: str, optional
            Term dividing individual metabolite entries.

        """
        if not isinstance(reaction_strings, list):
            reaction_strings = [reaction_strings]

        for rxn_string in reaction_strings:
            if not isinstance(rxn_string, string_types):
                raise TypeError("reaction_strings must be a string or a list "
                                "of strings")

        _metab_args = [_NAME_ARG_RE, _FORMULA_ARG_RE, _CHARGE_ARG_RE]

        for rxn_string in reaction_strings:
            if not _RXN_ID_RE.search(rxn_string):
                raise ValueError("Could not find an ID for "
                                 "'{0}'".format(rxn_string))
            result = _RXN_ID_RE.search(rxn_string)
            rxn_id = result.group(1)
            rxn_string = rxn_string[result.end():]
            # Determine reaction reversibility
            if _REVERSIBLE_ARROW_RE.search(rxn_string):
                arrow_loc = _REVERSIBLE_ARROW_RE.search(rxn_string)
                reversible = True
                # Reactants left of the arrow, products on the right
                reactant_str = rxn_string[:arrow_loc.start()].strip()
                product_str = rxn_string[arrow_loc.end():].strip()
            elif _FORWARD_ARROW_RE.search(rxn_string):
                arrow_loc = _FORWARD_ARROW_RE.search(rxn_string)
                reversible = False
                # Reactants left of the arrow, products on the right
                reactant_str = rxn_string[:arrow_loc.start()].strip()
                product_str = rxn_string[arrow_loc.end():].strip()
            elif _REVERSE_ARROW_RE.search(rxn_string):
                arrow_loc = _REVERSE_ARROW_RE.search(rxn_string)
                reversible = False
                # Reactants right of the arrow, products on the left
                reactant_str = rxn_string[:arrow_loc.end()].strip()
                product_str = rxn_string[arrow_loc.start():].strip()
            else:
                raise ValueError("Unrecognized arrow for "
                                 "'{0}'".format(rxn_string))
            new_reaction = MassReaction(rxn_id, reversible=reversible)

            d_re = re.compile("(\d) ")
            for substr, factor in zip([reactant_str, product_str], [-1, 1]):
                if not substr:
                    continue
                for term in substr.split(term_split):
                    term = term.strip()
                    if re.match("nothing", term.lower()):
                        continue
                    # Find the compartment if it exists
                    if _COMPARTMENT_RE.search(term):
                        compartment = _COMPARTMENT_RE.search(term)
                        compartment = compartment.group(1).strip("[|]")
                        term = _COMPARTMENT_RE.sub("]", term)
                    else:
                        compartment = None
                    if d_re.search(term):
                        num = float(d_re.search(term).group(1)) * factor
                        metab_to_make = term[d_re.search(term).end():]
                    else:
                        num = factor
                        metab_to_make = term
                    # Find the metabolite's ID
                    if not _MET_ID_RE.search(metab_to_make):
                        raise ValueError("Could not locate the metabolite ID")
                    met_id = _MET_ID_RE.search(metab_to_make).group(1)
                    # Use the metabolite in the model if it exists
                    try:
                        metab = self.metabolites.get_by_id(met_id)
                    except KeyError:
                        metab = MassMetabolite(met_id)
                    # Set attributes for the metabolite
                    for arg in _metab_args:
                        if arg.search(metab_to_make):
                            attr = _EQUALS_RE.split(arg.pattern)[0]
                            val = arg.search(metab_to_make).group(1)
                            metab.__dict__[attr] = val
                    new_reaction.add_metabolites({metab: num})
            self.add_reactions(new_reaction)

    # Internal
    def _mk_stoich_matrix(self, matrix_type=None, dtype=None,
                          update_model=True):
        """Return the stoichiometric matrix for a given MassModel.

        The rows represent the chemical species and the columns represent the
        reaction. S[i, j] therefore contains the quantity of species 'i'
        produced (positive) or consumed (negative) by reaction 'j'.

        Parameters
        ----------
        matrix_type: {'dense', 'dok', 'lil', 'DataFrame', 'symbolic'}, optional
            A string identifiying the desired format for the stoichiometric
            matrix of the model. Types can include 'dense' for a standard
            numpy.array, 'dok' or 'lil' to obtain the corresponding
            scipy.sparse matrix, 'DataFrame' for a pandas.DataFrame, and
            'symbolic' for a sympy.MutableDenseMatrix. For all matrix types,
            species (excluding  genes) are the row indicies and reactions are
            the column indicies. If None, defaults to "dense".
        dtype: data-type, optional
            The desired array data-type for the stoichiometric matrix. If None,
            defaults to np.float64.
        update_model: bool, optional
            If True, will update the stored stoichiometric matrix, the matrix
            type, and the data-type for the model.

        Returns
        -------
        stoich_mat: matrix of class 'dtype'
            The stoichiometric matrix for the MassModel returned as the given
            matrix_type and with a data-type of 'dtype'.

        """
        if not isinstance(update_model, bool):
            raise TypeError("update_model must be a bool")

        # Set up for matrix construction if matrix types are correct.
        (matrix_constructor, matrix_type, dtype) = _get_matrix_constructor(
            matrix_type=matrix_type, dtype=dtype,
            matrix_type_default=self._matrix_type, dtype_default=self._dtype)

        stoich_mat = matrix_constructor((len(self.metabolites),
                                         len(self.reactions)))
        # Get the indicies for the species and reactions
        m_ind = self.metabolites.index
        r_ind = self.reactions.index
        # Build the matrix
        for rxn in self.reactions:
            for met, stoich in iteritems(rxn.metabolites):
                stoich_mat[m_ind(met), r_ind(rxn)] = stoich

        # Convert the matrix to the desired type
        stoich_mat = convert_matrix(stoich_mat, matrix_type=matrix_type,
                                    dtype=dtype,
                                    row_ids=[m.id for m in self.metabolites],
                                    col_ids=[r.id for r in self.reactions])
        # Update the stored stoichiometric matrix for the model if True
        if update_model:
            self._update_model_s(stoich_mat, matrix_type, dtype)

        return stoich_mat

    def _update_stoichiometry(self, reaction_list, matrix_type=None,
                              dtype=None):
        """Update the stoichiometric matrix.

        To efficiently update the stoichiometric matrix with additional
        metabolites and reactions, the matrix is first converted into a 'dok'
        matrix, updated with the new objects, and lastly converted back into
        the desired matrix type.

        Parameters
        ----------
        reaction_list: list of mass.MassReactions, optional
            A list of MassReactions to be added to the stoichiometric matrix.
        matrix_type: {'dense', 'dok', 'lil', 'DataFrame', 'symbolic'}
            The desired type after converting the matrix.

        Warnings
        --------
        This method is intended for internal use only. To safely update the
        stoichiometric matrix, use the MassModel.update_S method instead.

        """
        shape = (len(self.metabolites), len(self.reactions))
        # Check the matrix input type to ensure it is valid matrix.
        if matrix_type is None:
            matrix_type = self._matrix_type
        if dtype is None:
            dtype = self._dtype

        # Get the S matrix as a dok matrix
        stoich_mat = convert_matrix(self._S, matrix_type="dok", dtype=dtype)
        # Resize the matrix
        stoich_mat.resize(shape)

        # Update the matrix
        coefficient_dictionary = {}
        for rxn in reaction_list:
            for metab, coeff in rxn._metabolites.items():
                coefficient_dictionary[(self.metabolites.index(metab.id),
                                        self.reactions.index(rxn.id))] = coeff
        stoich_mat.update(coefficient_dictionary)

        # Convert the matrix to the desired type
        stoich_mat = convert_matrix(stoich_mat, matrix_type=matrix_type,
                                    dtype=dtype,
                                    row_ids=[m.id for m in self.metabolites],
                                    col_ids=[r.id for r in self.reactions])
        return stoich_mat

    def _update_model_s(self, stoich_mat, matrix_type, dtype):
        """Update the model stoichiometric matrix and properties.

        Warnings
        --------
        This method is intended for internal use only. To safely update the
        stoichiometric matrix, use the MassModel.update_S method instead.
        """
        self._S = stoich_mat
        self._matrix_type = matrix_type
        self._dtype = dtype

    def _repr_html_(self):
        """HTML representation of the overview for the MassModel."""
        try:
            dim_S = "{0}x{1}".format(self.S.shape[0], self.S.shape[1])
            rank = np.linalg.matrix_rank(self.S)
        except np.linalg.LinAlgError:
            dim_S = "0x0"
            rank = 0
        return """
            <table>
                <tr>
                    <td><strong>Name</strong></td><td>{name}</td>
                </tr><tr>
                    <td><strong>Memory address</strong></td><td>{address}</td>
                </tr><tr>
                    <td><strong>Stoichiometric Matrix</strong></td>
                    <td>{dim_stoich_mat}</td>
                </tr><tr>
                    <td><strong>Matrix Rank</strong></td>
                    <td>{mat_rank}</td>
                </tr><tr>
                    <td><strong>Matrix Type</strong></td>
                    <td>{S_type}</td>
                </tr><tr>
                    <td><strong>Number of Metabolites</strong></td>
                    <td>{num_metabolites}</td>
                </tr><tr>
                    <td><strong>Number of Reactions</strong></td>
                    <td>{num_reactions}</td>
                </tr><tr>
                    <td><strong>Number of Initial Conditions</strong></td>
                    <td>{num_ic}</td>
                </tr><tr>
                    <td><strong>Number of Forward Rate Constants</strong></td>
                    <td>{num_kfs}</td>
                </tr><tr>
                    <td><strong>Number of Equilibrium Constants</strong></td>
                    <td>{num_Keqs}</td>
                </tr><tr>
                    <td><strong>Number of Irreversible Reactions</strong></td>
                    <td>{num_irreversible}</td>
                </tr><tr>
                    <td><strong>Number of Exchanges</strong></td>
                    <td>{num_exchanges}</td>
                </tr><tr>
                    <td><strong>Number of Fixed Concentrations</strong></td>
                    <td>{num_fixed}</td>
                </tr><tr>
                    <td><strong>Number of Custom Rates</strong></td>
                    <td>{num_custom_rates}</td>
                </tr><tr>
                    <td><strong>Number of Genes</strong></td>
                    <td>{num_genes}</td>
                </tr><tr>
                    <td><strong>Modules</strong></td>
                    <td>{modules}</td>
                </tr><tr>
                    <td><strong>Compartments</strong></td>
                    <td>{compartments}</td>
                </tr><tr>
                    <td><strong>Units</strong></td>
                    <td>{units}</td>
                </tr>
            </table>
        """.format(name=self.id, address='0x0%x' % id(self),
                   dim_stoich_mat=dim_S,
                   mat_rank=rank,
                   S_type="{}, {}".format(self._matrix_type,
                                          self._dtype.__name__),
                   num_metabolites=len(self.metabolites),
                   num_reactions=len(self.reactions),
                   num_ic=len(self.initial_conditions),
                   num_kfs=len(self.parameters["kf"]),
                   num_Keqs=len(self.parameters["Keq"]),
                   num_irreversible=len(self.irreversible_reactions),
                   num_exchanges=len(self.exchanges),
                   num_fixed=len(self.fixed_concentrations),
                   num_custom_rates=len(self.custom_rates),
                   num_genes=len(self.genes),
                   modules="<br> ".join([str(m) for m in self.modules
                                         if m is not None]) + "</br>",
                   compartments=", ".join(v if v else k for k, v in
                                          iteritems(self.compartments)),
                   units=", ".join(v if v else k for
                                   k, v in iteritems(self.units)))

    # Dunders
    def __enter__(self):
        """Record all future changes for context management of the MassModel.

        Changes are undone when a call to __exit__ is received.
        """
        # Create a new context and add it to the stack
        try:
            self._contexts.append(HistoryManager())
        except AttributeError:
            self._contexts = [HistoryManager()]
        return self

    def __exit__(self, exception_type, exception_value, traceback):
        """Pop the top context manager and trigger the undo functions."""
        context = self._contexts.pop()
        context.reset()

    def __setstate__(self, state):
        """Ensure all Objects in the MassModel point to the MassModel."""
        self.__dict__.update(state)
        for attr in ['reactions', 'metabolites', 'genes']:
            for x in getattr(self, attr):
                x._model = self
        if not hasattr(self, "name"):
            self.name = ""

    def __getstate__(self):
        """Get the state for serialization.

        Ensures that the context stack is cleared prior to serialization,
        since partial functions cannot be pickled reliably
        """
        odict = self.__dict__.copy()
        odict['_contexts'] = []
        return odict
