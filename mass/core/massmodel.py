# -*- coding: utf-8 -*-
r"""
MassModel is a class for holding information regarding a :mod:`mass` model.

The :class:`MassModel` class inherits and extends the
:class:`~cobra.core.model.Model` class in :mod:`cobra`. It contains additional
information required for simulations and other :mod:`mass` functions and
workflows.
"""
import re
import warnings
from copy import copy, deepcopy
from functools import partial

import numpy as np

from six import integer_types, iteritems, iterkeys, itervalues, string_types

import sympy as sym

from cobra.core.dictlist import DictList
from cobra.core.object import Object
from cobra.medium import (
    find_boundary_types, find_external_compartment, sbo_terms)
from cobra.util.context import HistoryManager, get_context

from mass.core.massmetabolite import MassMetabolite
from mass.core.massreaction import MassReaction
from mass.core.units import UnitDefinition
from mass.util.expressions import create_custom_rate, strip_time
from mass.util.util import (
    _get_matrix_constructor, _make_logger, convert_matrix, ensure_iterable,
    get_object_attributes, get_subclass_specific_attributes)

# Set the logger
LOGGER = _make_logger(__name__)
"""logging.Logger: Logger for :mod:`~mass.core.massmodel` submodule."""

CHOPNSQ = ['C', 'H', 'O', 'P', 'N', 'S', 'q']
"""list: Contains the six most abundant elements and charge for molecules."""


class MassModel(Object):
    r"""Class representation of a model.

    Parameters
    ----------
    id_or_model : str, ~cobra.core.model.Model, MassModel
        A string identifier to associate with the model, or an existing
        :class:`MassModel`. If an existing :class:`MassModel` is provided,
        a new :class:`MassModel` object is instantiated with the same
        properties as the original model.
    name : str
        A human readable name for the model.
    matrix_type : str
        A string identifiying the desired format for the returned matrix.
        Valid matrix types include ``'dense'``, ``'dok'``, ``'lil'``,
        ``'DataFrame'``, and ``'symbolic'`` Default is ``'DataFrame'``.
        See the :mod:`~.linear` module documentation for more information
        on the ``matrix_type``.
    dtype : data-type
        The desired array data-type for the stoichiometric matrix. If ``None``
        then the data-type will default to ``numpy.float64``.

    Attributes
    ----------
    description : str
        A human-readable description of the model.
    reactions : ~cobra.core.dictlist.DictList
        A :class:`~cobra.core.dictlist.DictList` where the keys are reaction
        identifiers and the values are the associated
        :class:`~.MassReaction`\ s.
    metabolites : ~cobra.core.dictlist.DictList
        A :class:`~cobra.core.dictlist.DictList` where the keys are metabolite
        identifiers and the values are the associated
        :class:`~.MassMetabolite`\ s.
    genes : ~cobra.core.dictlist.DictList
        A :class:`~cobra.core.dictlist.DictList` where the keys are gene
        identifiers and the values are the associated
        :class:`~cobra.core.gene.Gene`\ s.
    enzyme_modules : ~cobra.core.dictlist.DictList
        A :class:`~cobra.core.dictlist.DictList` where the keys are enzyme
        module identifiers and the values are the associated
        :class:`~.EnzymeModuleDict`\ s.
    custom_rates : dict
        A dictionary to store custom rate expressions for specific reactions,
        where the keys are :class:`~.MassReaction`\ s and values are the
        custom rate expressions given as :mod:`sympy` expressions. Custom rate
        expressions will always be prioritized over automatically generated
        mass action rates.
    custom_parameters : dict
        A dictionary to store the custom parameters for the custom rates,
        where key:value pairs are the string identifiers for the parameters
        and their corresponding numerical value.
    boundary_conditions : dict
        A dictionary to store boundary conditions, where keys are string
        identifiers for 'boundary metabolites' of boundary reactions, and
        values are the corresponding boundary condition numerical value or
        function of time.
    compartments : dict
        A dictionary to store the compartment shorthands and their full names.
        Keys are the shorthands while values are the full names.
    units : ~cobra.core.dictlist.DictList
        :class:`~cobra.core.dictlist.DictList` of :class:`~.UnitDefinition`\ s
        to store in the model for referencing.

    Warnings
    --------
    Note that the :class:`MassModel` does NOT track units, and it is therefore
    incumbent upon the user to maintain unit consistency the model.

    """

    def __init__(self, id_or_model=None, name=None, matrix_type="DataFrame",
                 dtype=np.float64):
        """Initialize the MassModel."""
        super(MassModel, self).__init__(id_or_model, name)
        if isinstance(id_or_model, MassModel):
            # Instiantiate a new MassModel with state identical to
            # the provided MassModel object.
            self.__setstate__(id_or_model.__dict__)
            if not hasattr(self, "name"):
                self.name = None
        else:
            self.description = ''
            # Initialize DictLists for storing
            # reactions, metabolites, genes, and enzyme_modules
            self.reactions = DictList()
            self.metabolites = DictList()
            self.genes = DictList()
            self.enzyme_modules = DictList()
            # Initialize dictionaries for custom rates, custom parameters,
            # boundary conditions, compartments, and units.
            self.boundary_conditions = {}
            self.custom_rates = {}
            self.custom_parameters = {}
            self._compartments = {}
            self.units = DictList()
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
        """Return the stoichiometric matrix."""
        return self.update_S(matrix_type=self._matrix_type, dtype=self._dtype,
                             update_model=False)

    @property
    def S(self):
        """Shorthand method to get the stoichiometric matrix."""
        return self.stoichiometric_matrix

    @property
    def ordinary_differential_equations(self):
        """Return a dict of ODEs for the metabolites."""
        return {met: met.ode for met in self.metabolites}

    @property
    def odes(self):
        """Shorthand method to get ODEs for the metabolites."""
        return self.ordinary_differential_equations

    @property
    def initial_conditions(self):
        """Return a dict of all metabolite initial conditions."""
        return {met: met.initial_condition for met in self.metabolites
                if met.initial_condition is not None}

    @property
    def ics(self):
        """Shorthand method to get all metabolite initial conditions."""
        return self.initial_conditions

    @property
    def fixed(self):
        """Return a dict of all metabolite fixed conditions."""
        return {met: ic for met, ic in iteritems(self.initial_conditions)
                if met.fixed}

    @property
    def rates(self):
        """Return a dict of reaction rate expressions.

        If a reaction has an associated custom rate expression, the custom rate
        will be prioritized and returned in the dictionary instead of the
        automatically generated mass action rate law expression.
        """
        return self.get_rate_expressions(
            self.reactions, rtype=0, update_reactions=True)

    @property
    def steady_state_fluxes(self):
        """Return a dict of all reaction steady state fluxes."""
        return {rxn: rxn.steady_state_flux for rxn in self.reactions
                if rxn.steady_state_flux is not None}

    @property
    def v(self):
        """Shorthand method to get all reaction steady state fluxes."""
        return self.steady_state_fluxes

    @property
    def boundary(self):
        """Return a list of boundary reactions in the model."""
        return [rxn for rxn in self.reactions if rxn.boundary]

    @property
    def boundary_metabolites(self):
        """Return a sorted list of all 'boundary metabolites' in the model.

        See Also
        --------
        :attr:`.MassReaction.boundary_metabolite`

        """
        return sorted(list(set(rxn.boundary_metabolite
                               for rxn in self.reactions if rxn.boundary)))

    @property
    def exchanges(self):
        """Return exchange reactions in the model.

        Exchange reactions are reactions that exchange mass with the exterior.
        Uses annotations and heuristics to exclude non-exchanges such as sink
        and demand reactions.
        """
        return find_boundary_types(self, "exchange", None)

    @property
    def demands(self):
        """Return demand reactions in the model.

        Demands are irreversible reactions that accumulate or consume a
        metabolite in the inside of the model.
        """
        return find_boundary_types(self, "demand", None)

    @property
    def sinks(self):
        """Return sink reactions in the model.

        Sinks are reversible reactions that accumulate or consume a metabolite
        in the inside of the model.
        """
        return find_boundary_types(self, "sink", None)

    @property
    def irreversible_reactions(self):
        """Return a list of all irreversible reactions in the model."""
        return [rxn for rxn in self.reactions if not rxn.reversible]

    @property
    def parameters(self):
        """Return all parameters associateed with the model."""
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
        # Add fluxes, custom parameters, and fixed concentrations.
        parameters.update({"v": {
            str(rxn.flux_symbol): flux
            for rxn, flux in iteritems(self.steady_state_fluxes)}})
        parameters.update({"Custom": self.custom_parameters})
        parameters.update({"Boundary": self.boundary_conditions})

        return parameters

    @property
    def compartments(self):
        """Get or set a dict of all metabolite compartments.

        Assigning a dictionary to this property updates the model's
        dictionary of compartment descriptions with the new values.

        Parameters
        ----------
        compartment_dict : dict
            Dictionary mapping compartments abbreviations to full names.
            An empty dictionary will reset the compartments.

        """
        return {met.compartment: self._compartments.get(met.compartment, '')
                for met in self.metabolites if met.compartment is not None}

    @compartments.setter
    def compartments(self, compartment_dict):
        """Set the dictionary of current compartment descriptions."""
        if compartment_dict:
            self._compartments.update(compartment_dict)
        else:
            setattr(self, "_compartments", {})

    def print_attributes(self, sep=r"\n", exclude_parent=False):
        r"""Print the attributes and properties of the MassModel.

        Parameters
        ----------
        sep : str
            The string used to seperate different attrubutes. Affects how the
            final string is printed. Default is ``'\n'``.
        exclude_parent : bool
            If ``True``, only display attributes specific to the current class,
            excluding attributes from the parent class.

        """
        if not isinstance(sep, str):
            raise TypeError("sep must be a string")

        if exclude_parent:
            attributes = get_subclass_specific_attributes(self)
        else:
            attributes = get_object_attributes(self)

        print(sep.join(attributes))

    def update_S(self, reaction_list=None, matrix_type=None, dtype=None,
                 update_model=True):
        r"""Update the stoichiometric matrix of the model.

        Notes
        -----
        ``reaction_list`` is assumed to be at the end of :attr:`reactions`.

        Parameters
        ----------
        reaction_list : list
            A list containing :class:`~.MassReaction`\ s to add to
            the stoichiometric matrix. The reactions must already exist in
            the model in order to update. If ``None``, the entire
            stoichiometric matrix is reconstructed.
        matrix_type : str
            A string identifiying the desired format for the returned matrix.
            Valid matrix types include ``'dense'``, ``'dok'``, ``'lil'``,
            ``'DataFrame'``, and ``'symbolic'``
            Default is the current ``matrix_type``. See the :mod:`~.linear`
            module documentation for more information on the ``matrix_type``.
        dtype : data-type
            The desired array data-type for the stoichiometric matrix.
            If ``None`` then the data-type will default to the
            current ``dtype``.
        update_model : bool
            If ``True``, will update the stored stoichiometric matrix,
            the matrix type, and the data-type for the model.

        Returns
        -------
        matrix of type ``matrix_type``
            The stoichiometric matrix for the :class:`~.MassModel` returned
            as the given ``matrix_type`` and with a data-type of ``dtype``.

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

    def add_metabolites(self, metabolite_list):
        r"""Add a list of metabolites to the model.

        The change is reverted upon exit when using the :class:`~.MassModel`
        as a context.

        Parameters
        ----------
        metabolite_list : list
            A list containing :class:`~.MassMetabolite`\ s to add to
            the model.

        """
        # Ensure list is iterable.
        metabolite_list = ensure_iterable(metabolite_list)

        # Check whether a metabolite is a MassMetabolite object, then check if
        # the metabolite already exists in the model, ignoring those that do.
        for met in metabolite_list:
            if not isinstance(met, MassMetabolite):
                warnings.warn("Skipping {0}, not a MassMetabolite".format(met))
                metabolite_list.remove(met)

        existing_mets = [met for met in metabolite_list
                         if met.id in self.metabolites]
        for met in existing_mets:
            LOGGER.info("Skipped adding MassMetabolite %s as it already exists"
                        " in the MassModel.", met.id)
        metabolite_list = [met for met in metabolite_list
                           if met not in self.metabolites]

        # Have metabolites point to the model object.
        for met in metabolite_list:
            met._model = self
        # Add metabolites to the model
        self.metabolites += metabolite_list

        context = get_context(self)
        if context:
            context(partial(self.metabolites.__isub__, metabolite_list))
            context(partial(setattr, met, "_model", None)
                    for met in metabolite_list)

    def remove_metabolites(self, metabolite_list, destructive=False):
        r"""Remove a list of metabolites from the model.

        The change is reverted upon exit when using the :class:`~.MassModel`
        as a context.

        Parameters
        ----------
        metabolite_list : list
            A list containing :class:`~.MassMetabolite`\ s to remove
                from the model.
        destructive : bool
            If ``False``, the metabolites are removed from all associated
            reactions. If ``True``, also remove associated
            :class:`~.MassReaction`\ s from the model.

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
        # Remove the metabolites
        self.metabolites -= metabolite_list

        context = get_context(self)
        if context:
            context(partial(self.metabolites.__iadd__, metabolite_list))
            context(partial(setattr, met, "_model", self)
                    for met in metabolite_list)

    def add_boundary_conditions(self, boundary_conditions):
        """Add boundary conditions values for the given boundary metabolites.

        Boundary condition values can be a numerical value, or they can be a
        string or :mod:`sympy` expression representing a function of time.
        The function must only depend on time.

        Parameters
        ----------
        boundary_conditions : dict
            A dict of boundary conditions containing the 'boundary metabolites'
            and their corresponding value. The string representing the
            'boundary_metabolite' must exist the list returned by
            :attr:`MassModel.boundary_metabolites`.

        See Also
        --------
        :attr:`boundary_metabolites`

        """
        if not isinstance(boundary_conditions, dict):
            raise TypeError("boundary_conditions must be a dict.")

        boundary_conditions_to_set = boundary_conditions.copy()
        for bound_met, bound_cond in iteritems(boundary_conditions_to_set):
            if bound_met not in self.boundary_metabolites and \
               bound_met not in self.metabolites:
                raise ValueError("Did not find {0} in model metabolites or in "
                                 "boundary reactions.".format(bound_met))
            # Boundary condition is a function
            if isinstance(bound_cond, (sym.Basic, string_types)):
                if isinstance(bound_cond, string_types):
                    bound_cond = sym.sympify(bound_cond)
                for arg in list(bound_cond.atoms(sym.Function)):
                    if arg.atoms(sym.Symbol) != {sym.Symbol("t")}:
                        raise ValueError(
                            "Function '{0}' for '{1}' has independent "
                            "variables {2}, expecting only {{t}}".format(
                                str(bound_cond), bound_met,
                                str(arg.atoms(sym.Symbol))))
                boundary_conditions_to_set[bound_met] = bound_cond
            # Boundary condition is an integer or float
            elif isinstance(bound_cond, (integer_types, float)):
                boundary_conditions_to_set[bound_met] = float(bound_cond)
            else:
                raise TypeError(
                    "Invalid boundary value for '{0}'. Boundary conditions can"
                    "only be numerical values or functions of time.".format(
                        bound_met))
        # Keep track of existing concentrations for context management.
        context = get_context(self)
        if context:
            existing_concs = {bound_met: self.boundary_conditions[bound_met]
                              for bound_met in boundary_conditions
                              if bound_met in self.boundary_conditions}

        self.boundary_conditions.update(boundary_conditions_to_set)

        if context:
            context(partial(self.boundary_conditions.pop, key)
                    for key in iterkeys(boundary_conditions)
                    if key not in existing_concs)
            context(partial(self.boundary_conditions.update, existing_concs))

    def remove_boundary_conditions(self, boundary_metabolite_list):
        """Remove the boundary condition for a list of `boundary metabolites`.

        Parameters
        ----------
        metabolite_list : list
            A list of metabolites to remove the boundary conditions for.
            Boundary metabolites must already exist in the model in order
            for them to be removed.

        See Also
        --------
        :attr:`boundary_metabolites`

        """
        boundary_metabolite_list = ensure_iterable(boundary_metabolite_list)
        # Check whether a metabolite already exists in the model,
        # ignoring those that do not.
        boundary_metabolite_list = [
            bound_met for bound_met in boundary_metabolite_list
            if bound_met in self.boundary_conditions]

        # Keep track of existing concentrations for context management.
        context = get_context(self)
        if context:
            existing_concs = {bound_met: self.boundary_conditions[bound_met]
                              for bound_met in boundary_metabolite_list
                              if bound_met in self.boundary_conditions}
        # Remove the initial conditions
        for bound_met in boundary_metabolite_list:
            del self.boundary_conditions[bound_met]

        if context:
            context(partial(self.boundary_conditions.update, existing_concs))

    def add_reactions(self, reaction_list):
        r"""Add reactions to the model.

        :class:`~.MassReaction`\ s with identifiers identical to an
        existing reaction are ignored.

        The change is reverted upon exit when using the :class:`~.MassModel`
        as a context.

        Parameters
        ----------
        reaction_list : list
            A list of :class:`~.MassReaction`\ s to add to the model.

        """
        # Ensure list is iterable.
        reaction_list = ensure_iterable(reaction_list)

        # Check whether a metabolite is a MassReaction object, then check if
        # the reaction already exists in the model, ignoring those that do.
        for rxn in reaction_list:
            if not isinstance(rxn, MassReaction):
                warnings.warn(
                    "Skipping {0}, not a MassReaction object".format(rxn))
                reaction_list.remove(rxn)

        # Check whether reactions exist in the model.
        pruned = self._existing_obj_filter("reactions", reaction_list)
        context = get_context(self)
        # Add reactions, and have reactions point to the model.
        for rxn in pruned:
            rxn._model = self
            # Loop through metabolites associated with a reaction.
            for met in list(rxn._metabolites):
                # If a metabolite doesn't exist in the model,
                # add it with its associated initial condition.
                if met not in self.metabolites:
                    self.add_metabolites(met)
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
        r"""Remove reactions from the model.

        The change is reverted upon exit when using the :class:`~.MassModel`
        as a context.

        Parameters
        ----------
        reaction_list : list
            A list of :class:`~.MassReaction`\ s to be removed
            from the model.
        remove_orphans : bool
            If ``True``, will also remove orphaned genes and metabolites from
            the model.

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
            if rxn in self.custom_rates:
                self.remove_custom_rate(rxn)

            if context:
                context(partial(setattr, rxn, "_model", self))
        # Remove reactions from the model
        self.reactions -= reaction_list

        if context:
            context(partial(self.reactions.__iadd__, reaction_list))
            context(partial(setattr, rxn, "_model", self)
                    for rxn in reaction_list)

    def add_boundary(self, metabolite, boundary_type="exchange",
                     reaction_id=None, subsystem="", boundary_condition=0.,
                     sbo_term=None):
        """Add a boundary reaction for a given metabolite.

        There are three different types of pre-defined boundary reactions:
        exchange, demand, and sink reactions.

            * An exchange reaction is a reversible, unbalanced reaction that
              adds to or removes an extracellular metabolite from the
              extracellular compartment.
            * A demand reaction is an irreversible reaction that consumes an
              intracellular metabolite.
            * A sink is similar to an exchange but specifically for
              intracellular metabolites.

        Notes
        ------
        To set the reaction ``boundary_type`` to something else, the desired
        identifier of the created reaction must be specified. The name will
        be given by the metabolite name and the given ``boundary_type``, and
        the reaction will be set its reversible attribute to ``True``.

        Parameters
        ----------
        metabolite : MassMetabolite
            Any :class:`~.MassMetabolite`, or an identifier of a metabolite
            that already exists in the model. The metabolite compartment is
            not checked but it is encouraged to stick to the definition of
            exchanges, demands, and sinks.
        boundary_type : str
            One of the pre-defined boundary types, or a user-defined type.
            Pre-defined boundary types include ``"exchange"``, ``"demand"``,
            and ``"sink"``. Using one of the pre-defined reaction types is
            easiest. To  create a user-defined kind of boundary reaction
            choose any other string, e.g., 'my-boundary'.
        reaction_id : str
            The ID of the resulting reaction. This takes precedence over the
            auto-generated identifiers but beware that it might make boundary
            reactions harder to identify afterwards when using
            :attr:`~.MassModel.boundary` or specifically
            :attr:`~.MassModel.exchanges` etc.
        subsystem : str
            The subsystem where the reaction is meant to occur.
        boundary_condition : float, str, ~sympy.core.basic.Basic
            The boundary condition value to set. Must be an ``int``,
            ``float``, or a :mod:`sympy` expression dependent only on time.
            Default value is ``0.``
        sbo_term : str
            A correct SBO term is set for the available types. If a custom
            type is chosen, a suitable SBO term should also be set.

        Returns
        -------
        MassReaction
            The :class:`~.MassReaction` of the new boundary reaction.

        """
        # Check whether a metabolite is a MassMetabolite object:
        if isinstance(metabolite, string_types):
            try:
                metabolite = self.metabolites.get_by_id(metabolite)
            except KeyError:
                raise ValueError("metabolite must exist in the model")
        elif not isinstance(metabolite, MassMetabolite):
            raise TypeError("metabolite must be a MassMetabolite object")

        boundary_types = {
            "exchange": ("EX", True, sbo_terms["exchange"]),
            "demand": ("DM", False, sbo_terms["demand"]),
            "sink": ("SK", True, sbo_terms["sink"])}

        if boundary_type == "exchange":
            external = find_external_compartment(self)
            if metabolite.compartment != external:
                raise ValueError("The metabolite is not an external metabolite"
                                 " (compartment is `{0}` but should be `{1}`)."
                                 " Did you mean to add a demand or sink? "
                                 "If not, either change its compartment or "
                                 "rename the model compartments to fix this."
                                 .format(metabolite.compartment, external))

        if boundary_type in boundary_types:
            prefix, reversible, default_term = boundary_types[boundary_type]
            if reaction_id is None:
                reaction_id = "{0}_{1}".format(prefix, metabolite.id)
            if sbo_term is None:
                sbo_term = default_term
        if reaction_id is None:
            raise ValueError(
                "Custom types of boundary reactions require a custom "
                "identifier. Please set the `reaction_id`.")
        # Return the existing reaction if the reaction already exists
        if reaction_id in self.reactions:
            raise ValueError(
                "Boundary reaction '{0}' already exists.".format(reaction_id))
        # Set reaction name, subsystem, and reversible attributes
        name = "{0}_{1}".format(metabolite.name, boundary_type)
        rxn = MassReaction(reaction_id, name=name, subsystem=subsystem,
                           reversible=reversible)
        # Always add metabolite as a reactant
        rxn.add_metabolites({metabolite: -1})
        # Add SBO annotation
        if sbo_term:
            rxn.annotation["sbo"] = sbo_term
        self.add_reactions(rxn)
        self.add_boundary_conditions({
            rxn.boundary_metabolite: boundary_condition})

        return rxn

    def get_rate_expressions(self, reaction_list=None, rtype=0,
                             update_reactions=False):
        r"""Get the rate expressions for a list of reactions in the model.

        Parameters
        ----------
        reaction_list : list
            A list of :class:`~.MassReaction`\ s to get the rate
            expressions for. Reactions must already exist in the model.
            If ``None``, then return the rates for all reactions in the model.
        rtype : int
            The type of rate law to display. Must be 0, 1, 2, or 3.

                * If 0, the currrent default rate law type is used.
                * Type 1 will utilize the :attr:`forward_rate_constant` and the
                  :attr:`equilibrium_constant`.
                * Type 2 will utilize the :attr:`forward_rate_constant` and the
                  :attr:`reverse_rate_constant`.
                * Type 3 will utilize the :attr:`equilibrium_constant` and the
                  :attr:`reverse_rate_constant`.

            Default is `` 0.``
        update_reactions : bool
            If ``True``, update the :class:`~.MassReaction` default rate type
            in addition to returning the rate expressions.

        Returns
        -------
        dict
            Dictionary of reaction rates where keys are the reaction ids
            and values are the rate law expressions.

        """
        # Check the inputs
        if rtype not in {0, 1, 2, 3}:
            raise ValueError("rate_type must be 0, 1, 2, or 3")

        # Use the MassModel reactions if no reaction list is given
        if reaction_list is None:
            reaction_list = self.reactions
        # Ensure list is iterable.
        reaction_list = ensure_iterable(reaction_list)

        if rtype == 0:
            rate_dict = {
                rxn: rxn.get_mass_action_rate(rxn._rtype, update_reactions)
                for rxn in reaction_list}

        else:
            rate_dict = {
                rxn: rxn.get_mass_action_rate(rtype, update_reactions)
                for rxn in reaction_list}

        if self.custom_rates:
            rate_dict.update(self.custom_rates)

        return rate_dict

    def get_mass_action_ratios(self, reaction_list=None, sympy_expr=True):
        r"""Get the mass action ratios for a list of reactions in the model.

        Parameters
        ----------
        reaction_list : list
            A list of :class:`~.MassReaction`\ s to get the mass action
            ratios for. Reactions must already exist in the model.
            If ``None``, then return the ratios for all reactions in the model.
        sympy_expr : bool
            If ``True`` then return the mass action ratios as a :mod:`sympy`
            expression, otherwise return the ratio as a human readable string.

        Returns
        -------
        dict
            Dictionary of mass action ratios where keys are the
            reaction ids and values are the ratios.

        """
        # Use the MassModel reactions if no reaction list is given
        if reaction_list is None:
            reaction_list = self.reactions
        # Ensure list is iterable.
        reaction_list = ensure_iterable(reaction_list)

        ratio_dict = dict((rxn, rxn.get_mass_action_ratio()) if sympy_expr
                          else (rxn, str(rxn.get_disequilibrium_ratio()))
                          for rxn in reaction_list)
        return ratio_dict

    def get_disequilibrium_ratios(self, reaction_list=None, sympy_expr=True):
        r"""Get the disequilibrium ratios for a list of reactions in the model.

        Parameters
        ----------
        reaction_list : list
            A list of :class:`~.MassReaction`\ s to get the disequilibrium
            ratios for. Reactions must already exist in the model.
            If ``None``, then return the ratios for all reactions in the model.
        sympy_expr : bool
            If ``True`` then return the disequilibrium ratios as a :mod:`sympy`
            expression, otherwise return the ratio as a human readable string.

        Returns
        -------
        dict
            Dictionary of mass action ratios where keys are the
            reaction ids and values are the ratios.

        """
        # Use the MassModel reactions if no reaction list is given
        if reaction_list is None:
            reaction_list = self.reactions
        # Ensure list is iterable.
        reaction_list = ensure_iterable(reaction_list)

        ratio_dict = dict((rxn, rxn.get_disequilibrium_ratio()) if sympy_expr
                          else (rxn, str(rxn.get_disequilibrium_ratio()))
                          for rxn in reaction_list)
        return ratio_dict

    def add_custom_rate(self, reaction, custom_rate, custom_parameters=None):
        """Add a custom rate for a reaction to the model.

        The change is reverted upon exit when using the :class:`~.MassModel`
        as a context.

        Notes
        -----
        * Metabolites must already exist in the :class:`MassModel`.
        * Default parameters of a :class:`~.MassReaction` are automatically
          taken into account and do not need to be defined as additional
          custom parameters.

        Parameters
        ----------
        reaction : MassReaction
            The reaction associated with the custom rate.
        custom_rate : str
            The string representation of the custom rate expression. The string
            representation of the custom rate will be used to create a
            :mod:`sympy` expression that represents the custom rate.
        custom_parameters : dict
            A dictionary of custom parameters for the custom rate where the
            key:value pairs are the strings representing the custom parameters
            and their numerical values. The string representation of the custom
            parametes will be used to create the symbols needed for the sympy
            expression of the custom rate. If ``None``, then parameters are
            assumed to already exist in the model.

        See Also
        --------
        ~MassReaction.all_parameter_ids
            List of default reaction parameters automatically accounted for.

        """
        if custom_parameters is not None:
            custom_parameter_list = list(iterkeys(custom_parameters))
        else:
            custom_parameters = {}
            custom_parameter_list = []

        # Ensure custom rate is a string
        if not isinstance(custom_rate, string_types):
            custom_rate = str(custom_rate)

        # Use any existing custom parameters if they are in the rate law.
        existing_customs = self.custom_parameters
        if existing_customs:
            for custom_parameter in iterkeys(existing_customs):
                if re.search(custom_parameter, custom_rate) and \
                   custom_parameter not in custom_parameter_list:
                    custom_parameter_list.append(custom_parameter)
        # Create the custom rate expression
        custom_rate = create_custom_rate(reaction, custom_rate,
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

        The change is reverted upon exit when using the :class:`~.MassModel`
        as a context.

        Parameters
        ----------
        reaction : MassReaction
            The reaction assoicated with the custom rate to be removed.
        remove_orphans : bool
            If ``True``, then remove any orphaned custom parameters from
            the model.

        """
        try:
            rate_to_remove = self.custom_rates[reaction]
        except KeyError:
            warnings.warn("Did not find a custom custom rate expression "
                          "associated with reaction {0}.".format(reaction.id))
            return
        # Remove the rate
        del self.custom_rates[reaction]

        # Remove orphaned custom parameters if desired.
        args = rate_to_remove.atoms(sym.Symbol)
        standards = [
            sym.Symbol(arg) for arg in reaction.all_parameter_ids + ["t"]]
        # Save currently existing parameters for context management if needed.
        existing = {str(arg): self.custom_parameters[str(arg)]
                    for arg in args if arg in self.custom_parameters and
                    arg not in standards}

        if remove_orphans and self.custom_rates:
            # Determine symbols still in use.
            other_args = set()
            for custom_rate in itervalues(self.custom_rates):
                other_args.update(custom_rate.atoms(sym.Symbol))

            # Remove those that are not being used in any custom rate.
            for arg in other_args:
                if arg in args.copy():
                    args.remove(arg)

            for arg in args:
                if arg not in standards and arg in self.custom_parameters:
                    del self.custom_parameters[str(arg)]

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
        parameters in the model. To remove a specific rate without affecting
        the other custom rates or parameters, use :meth:`remove_custom_rate`
        instead.

        """
        self.custom_rates = {}
        self.custom_parameters = {}
        print("All custom rate expressions and parameters have been reset")

    def add_units(self, unit_defs):
        r"""Add a :class:`~.UnitDefinition` to the model :attr:`units`.

        Notes
        -----
        The model will not automatically track or convert units. Therefore,
        it is up to the user to ensure unit consistency in the model.

        Parameters
        ----------
        unit_defs : list
            A list of :class:`~.UnitDefinition`\ s to add to the model.

        """
        # Ensure iterable input and units are valid Unit objects.
        unit_defs = DictList(ensure_iterable(unit_defs))
        for unit in list(unit_defs):
            if not isinstance(unit, UnitDefinition):
                raise ValueError(
                    "'{0}' is not a valid UnitDefinition.".format(str(unit)))
            # Skip existing units.
            if unit.id in self.units.list_attr("id"):
                warnings.warn("Skipping '{0}' for it already exists in the"
                              " model.".format(unit))
                unit_defs.remove(unit)
        # Add new unit definitions to the units attribute
        self.units += unit_defs

    def remove_units(self, unit_defs):
        r"""Remove a :class:`~.UnitDefinition` from the model :attr:`units`.

        Notes
        -----
        The model will not automatically track or convert units. Therefore,
        it is up to the user to ensure unit consistency in the model.

        Parameters
        ----------
        unit_defs : list
            A list of :class:`~.UnitDefinition`\ s or their string identifiers
            to remove from the model.

        """
        unit_defs = ensure_iterable(unit_defs)
        # Iteratre through units, raise ValueError if unit does not exist.
        for unit in unit_defs:
            try:
                unit = self.units.get_by_id(unit)
            except KeyError as e:
                raise ValueError(
                    "'{0}' does not exist in the model.".format(str(e)))
        # Remove unit definitions to the units attribute
        self.units -= unit_defs

    def reset_units(self):
        r"""Reset all unit definitions in a model.

        Warnings
        --------
        Using this method will remove all :class:`~.UnitDefinition`\ s from
        the model. To remove a :class:`~.UnitDefinition` without affecting
        other units, use :meth:`remove_units` instead.

        """
        self.units = DictList()
        print("All unit definitions have been reset")

    def get_elemental_matrix(self, matrix_type=None, dtype=None):
        """Get the elemental matrix for a model.

        Parameters
        ----------
        matrix_type : str
            A string identifiying the desired format for the returned matrix.
            Valid matrix types include ``'dense'``, ``'dok'``, ``'lil'``,
            ``'DataFrame'``, and ``'symbolic'``
            Default is ``'dense'``. See the :mod:`~.linear` module
            documentation for more information on the ``matrix_type``.
        dtype : data-type
            The desired array data-type for the matrix. If ``None`` then
            the data-type will default to ``numpy.float64``.

        Returns
        -------
        matrix of type ``matrix_type``
            The elemntal matrix for the :class:`~.MassModel` returned
            as the given ``matrix_type`` and with a data-type of ``dtype``.

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
        moieties = []
        for met in self.metabolites:
            element_dict = met.elements
            for element in CHOPNSQ:
                if element in iterkeys(element_dict):
                    amount = element_dict.pop(element)
                elif re.match("q", element) and met.charge is not None:
                    amount = met.charge
                else:
                    amount = 0
                elem_mat[e_ind(element), m_ind(met)] = amount
            # Get any additional moieties
            moieties.extend([
                element for element in iterkeys(element_dict)
                if "[" in element and "]" in element
                and element not in moieties])

        row_ids = CHOPNSQ.copy()
        # Add additional moieties to the elemental matrix
        if moieties:
            moiety_mat = matrix_constructor((
                len(moieties), len(self.metabolites)))
            # Create additional matrix for moieties
            for met in self.metabolites:
                element_dict = met.elements
                for moiety in moieties:
                    if moiety in element_dict:
                        amount = 1
                    else:
                        amount = 0
                    moiety_mat[moieties.index(moiety), m_ind(met)] = amount
            # Concatenate matrices
            elem_mat = np.concatenate((elem_mat, moiety_mat), axis=0)
            row_ids.extend(moieties)

        # Convert matrix to a dataframe if matrix type is a dataframe
        elem_mat = convert_matrix(elem_mat, matrix_type=matrix_type,
                                  dtype=dtype, row_ids=row_ids,
                                  col_ids=[m.id for m in self.metabolites])

        return elem_mat

    def get_elemental_charge_balancing(self, matrix_type=None, dtype=None):
        """Get the elemental charge balance as a matrix for a model.

        Parameters
        ----------
        matrix_type : str
            A string identifiying the desired format for the returned matrix.
            Valid matrix types include ``'dense'``, ``'dok'``, ``'lil'``,
            ``'DataFrame'``, and ``'symbolic'``
            Default is ``'dense'``. See the :mod:`~.linear` module
            documentation for more information on the ``matrix_type``.
        dtype : data-type
            The desired array data-type for the matrix. If ``None`` then
            the data-type will default to ``numpy.float64``.

        Returns
        -------
        matrix of type ``matrix_type``
            The charge balancing matrix for the :class:`~.MassModel` returned
            as the given ``matrix_type`` and with a data-type of ``dtype``.

        """
        elem_mat = self.get_elemental_matrix(matrix_type="DataFrame")
        row_ids = elem_mat.index

        stoich_mat = self.update_S(matrix_type="dense", update_model=False)
        charge_mat = np.array(elem_mat).dot(stoich_mat)

        if matrix_type is None:
            matrix_type = "dense"
        if dtype is None:
            dtype = np.float64

        charge_mat = convert_matrix(charge_mat, matrix_type=matrix_type,
                                    dtype=dtype, row_ids=row_ids,
                                    col_ids=[r.id for r in self.reactions])
        return charge_mat

    def repair(self, rebuild_index=True, rebuild_relationships=True):
        """Update all indicies and pointers in the model.

        Parameters
        ----------
        rebuild_index: bool
            If ``True``, then rebuild the indicies kept in the reactions,
            metabolites, and genes.
        rebuild_relationships: bool
            If ``True``, then reset all associations between the reactions,
            metabolites, genes, enzyme_modules, and the :class:`MassModel`,
            and rebuilds them.

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
            self.enzyme_modules._generate_index()
            self.units._generate_index()
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
        for attr in ["reactions", "metabolites", "genes", "enzyme_modules"]:
            for item in getattr(self, attr):
                item._model = self

    def copy(self):
        r"""Create a partial "deepcopy" of the :class:`MassModel`.

        All of the :class:`~.MassMetabolite`\ s, :class:`~.MassReaction`\ s,
        :class:`~cobra.core.gene.Gene`\ s and :class:`~.EnzymeModuleDict`\ s,
        the boundary conditions, custom_rates, custom_parameters, and the
        stoichiometric matrix are created anew, but in a faster fashion
        than deepcopy.

        """
        # Define a new model
        new_model = self.__class__()
        # Define items that will not be copied by their references
        do_not_copy_by_ref = [
            "metabolites", "reactions", "genes", "enzyme_modules", "_S",
            "custom_rates", "units", "boundary_conditions",
            "custom_parameters", "notes", "annotation"]
        for attr in self.__dict__:
            if attr not in do_not_copy_by_ref:
                new_model.__dict__[attr] = self.__dict__[attr]

        for attr in do_not_copy_by_ref[-5:]:
            setattr(new_model, attr, deepcopy(getattr(self, attr)))

        # Copy the metabolites
        new_model = self._copy_model_metabolites(new_model)
        # Copy the genes
        new_model = self._copy_model_genes(new_model)
        # Copy the reactions and rates (including custom rates)
        new_model = self._copy_model_reactions(new_model)
        # Copy any existing enzyme_modules
        new_model = self._copy_model_enzyme_modules(new_model)
        # Create the new stoichiometric matrix for the model.
        new_model._S = self._mk_stoich_matrix(matrix_type=self._matrix_type,
                                              dtype=self._dtype,
                                              update_model=True)
        # Doesn't make sense to retain the context of a copied model so
        # assign a new empty context
        new_model._contexts = []

        return new_model

    def merge(self, right, prefix_existing=None, inplace=False,
              new_model_id=None):
        """Merge two models into one model with the objects from both.

        The reactions, metabolites, genes, enzyme modules, boundary conditions,
        custom rate expressions, rate parameters, compartments, units, notes,
        and annotations from the right model are also copied to left model.
        However, note that in cases where identifiers for objects are identical
        or a dict item has an identical key(s), priority will be given to what
        already exists in the left model.

        Notes
        -----
        When merging an :class:`.~EnzymeModule` into a :class:`MassModel`,
        the enzyme module` is converted to an :class:`.~EnzymeModuleDict` and
        stored in a :class:`~cobra.core.dictlist.DictList` accessible via the
        :attr:`enzyme_modules` attribute. If an :class:`.~EnzymeModuleDict`
        already exists in the model, it will be replaced.

        Parameters
        ----------
        right : MassModel
            The model to merge into the left model.
        prefix_existing : str
            If provided, the string is used to prefix the reaction identifier
            of a reaction in the second model if that reaction already exists
            within the left model. Will also apply prefix to enzyme identifiers
            of enzyme modules in the second model.
        inplace : bool
            If ``True`` then add reactions from right model directly to the
            left model. Otherwise, create a new model leaving the left model
            untouched. When done within the model as context, changes to the
            models are reverted upon exit.
        new_model_id : str
            If provided, the string is used as the identifier for the merged
            model. If ``None`` and inplace is ``True``, the model ID of the
            first model will be used. If ``None`` and inplace is ``False``,
            a new combined ID will be used for the new :class:`MassModel`.

        Returns
        -------
        MassModel
            A new :class:`MassModel` or ``self`` representing the
            merged model.

        """
        # Check whether two MassModels are being merged,
        # or if a MassModel and a MassModel subclass are being merged.
        if right.__class__ is not MassModel \
           and issubclass(right.__class__, MassModel):
            return right._add_self_to_model(self, prefix_existing,
                                            inplace, new_model_id)
        # Set the merged model object and its ID,
        # then add the module attribute of the right model into the left
        if not inplace:
            new_model = self.copy()
        else:
            new_model = self
        if new_model_id is None:
            new_model_id = {True: self.id,
                            False: self.id + "_" + right.id}.get(inplace)
        new_model.id = new_model_id

        # Add the reactions from right to left model.
        new_reactions = deepcopy(right.reactions)
        if prefix_existing is not None:
            existing = new_reactions.query(lambda r: r.id in self.reactions)
            for reaction in existing:
                reaction.id = prefix_existing + "_" + reaction.id
        new_model.add_reactions(new_reactions)
        new_model.repair()

        # Add boundary conditions from right to left model.
        existing = [bc for bc in iterkeys(new_model.boundary_conditions)]
        new_model.add_boundary_conditions({
            met: bc for met, bc in iteritems(right.boundary_conditions)
            if bc not in existing})

        # Add custom parameters from right to left model.
        existing = [cp for cp in iterkeys(new_model.custom_parameters)]
        new_model.custom_parameters.update({
            cp: v for cp, v in iteritems(right.custom_parameters)
            if cp not in existing})

        # Add custom rates from right to left model,
        # prefixing any existing reactions if necessary
        existing = dict((prefix_existing + "_" + rxn.id, rate)
                        if prefix_existing is not None else (rxn.id, rate)
                        for rxn, rate in iteritems(new_model.custom_rates))
        for rate_dict in [right.custom_rates, existing]:
            new_model.custom_rates.update({
                new_model.reactions.get_by_id(getattr(rxn, "_id", rxn)): rate
                for rxn, rate in iteritems(rate_dict)})

        # Add enzyme_modules from right to left model
        if right.enzyme_modules:
            new_enzyme_modules = deepcopy(right.enzyme_modules)
            # Prefix enzyme_modules if necessary
            if prefix_existing is not None:
                existing = new_enzyme_modules.query(
                    lambda e: e.id in self.enzyme_modules)
                for enzyme in existing:
                    enzyme.__dict__["_id"] = prefix_existing + "_" + enzyme.id
            # Check whether reactions exist in the model.
            new_enzyme_modules = self._existing_obj_filter(
                "enzyme_modules", new_enzyme_modules)
            new_model.enzyme_modules += new_enzyme_modules
            for enzyme in new_model.enzyme_modules:
                enzyme.model = new_model

        if right.units:
            new_model.add_units(right.units)

        for attr in ["_compartments", "notes", "annotation"]:
            existing = getattr(new_model, attr).copy()
            setattr(new_model, attr, getattr(right, attr).copy())
            getattr(new_model, attr).update(existing)

        return new_model

    def compute_steady_state_fluxes(self, pathways, independent_fluxes,
                                    update_reactions=False):
        r"""Calculate the unique steady state flux for each reaction.

        The unique steady state flux for each reaction in the
        :class:`MassModel` is calculated using defined pathways, independently
        defined fluxes, and steady state concentrations, where index of values
        in the pathways must correspond to the index of the reaction in
        :attr:`MassModel.reactions`.

        Notes
        -----
        The number of individually defined fluxes must be the same as the
        number of pathways in order to determine the solution. For best
        results, the number of pathways to specify must equal the dimension
        of the right nullspace.

        Parameters
        ----------
        pathways : array-like
            An array-like object that define the pathways through the reaction
            network of the model. The given pathway vectors must be the same
            length as the number of reactions in the model, with indicies of
            values in the pathway vector corresponding to the indicies of
            reactions in the :attr:`reactions` attribute.
        independent_fluxes : dict
            A dict of steady state fluxes where :class:`~.MassReaction`\ s are
            keys and fluxes are values to utilize in order to calculate all
            other steady state fluxes. Must be the same length as the number
            of specified pathways.
        update_reactions : bool
            If True, then update the :attr:`.MassReaction.steady_state_flux`
            with the calculated steady state flux value for each reaction.

        Return
        ------
        dict
            A dict where key:value pairs are the :class:`~.MassReaction`\ s
            with their corresponding calculated steady state fluxes.

        Warnings
        --------
        The indicies of the values in the pathway vector must correspond to the
        indicies of the reactions in the :attr:`reactions` attribute in order
        for the method to work as intended.

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
        flux_vector = np.inner(pathways.T, np.inner(coeffs, values))
        # Update the reactions if desired
        steady_state_fluxes = {}
        for i, rxn in enumerate(self.reactions):
            steady_state_flux = flux_vector[i]
            steady_state_fluxes.update({rxn: steady_state_flux})
            if update_reactions:
                rxn.steady_state_flux = steady_state_flux

        return steady_state_fluxes

    def calculate_PERCs(self, steady_state_fluxes=None,
                        steady_state_concentrations=None,
                        at_equilibrium_default=100000, update_reactions=False):
        r"""Calculate pseudo-order rate constants for reactions in the model.

        Pseudo-order rate constants (PERCs) are considered to be the same as
        :attr:`~.MassReaction.forward_rate_constant` attributes, and are
        calculated based on the steady state concentrations and fluxes.

        Parameters
        ----------
        steady_state_fluxes : dict
            A dict of steady state fluxes where :class:`~MassReaction`\ s
            are keys and fluxes are the values. All reactions provided will
            have their PERCs calculated. If ``None``, PERCs are calculated
            using the current steady state fluxes for all reactions that
            exist in the model.
        steady_state_concentrations : dict
            A dict of all steady state concentrations necessary for the PERC
            calculations, where :class:`~.MassMetabolite`\ s are keys and
            concentrations are the values. If ``None``, the relevant
            concentrations that exist in the model are used. All
            concentrations used in calculations must be provided, including
            relevant boundary conditions.
        at_equilibrium_default : float
            The value to set the pseudo-order rate constant if the reaction is
            at equilibrium. Default is ``100,000``.
        update_reactions : bool
            If ``True`` then will update the values for the
            :attr:`~MassReaction.forward_rate_constant` attributes with the
            calculated PERC values.

        Returns
        -------
        dict
            A dict where keys are strings identifers of the pseudo-order rate
            constants (as given by :attr:`.MassReaction.kf_str`) and values
            are the calculated PERC values.

        """
        # Get the model steady state concentrations if None are provided.
        if steady_state_concentrations is None:
            steady_state_concentrations = self.boundary_conditions.copy()
            steady_state_concentrations.update({
                str(m): ic for m, ic in iteritems(self.initial_conditions)})
        else:
            steady_state_concentrations = {
                str(m): v for m, v in iteritems(steady_state_concentrations)}
        # Get the model reactions and fluxes if None are provided.
        if steady_state_fluxes is None:
            steady_state_fluxes = self.steady_state_fluxes

        # Get defined numerical values
        numerical_values = {
            str(param): value
            for p_type, subdict in iteritems(self.parameters)
            for param, value in iteritems(subdict) if p_type != "kf"}
        numerical_values.update(steady_state_concentrations)

        # Function to calculate the solution
        def calculate_sol(flux, rate_equation, perc):
            sol = sym.solveset(sym.Eq(flux, rate_equation),
                               perc, domain=sym.S.Reals)
            if isinstance(sol, type(sym.S.Reals)) or sol.is_EmptySet \
               or next(iter(sol)) <= 0:
                sol = float(at_equilibrium_default)
            else:
                sol = float(next(iter(sol)))

            return sol

        percs_dict = {}
        for reaction, flux in iteritems(steady_state_fluxes):
            # Ensure inputs are correct
            if not isinstance(reaction, MassReaction)\
               or not isinstance(flux, (float, integer_types)):
                raise TypeError(
                    "steady_state_fluxes must be a dict containing the "
                    "MassReactions and their steady state flux values.")
            rate_eq = strip_time(reaction.rate)
            arguments = list(rate_eq.atoms(sym.Symbol))
            # Check for missing numerical values
            missing_values = [
                str(arg) for arg in arguments
                if str(arg) not in numerical_values
                and str(arg) != reaction.kf_str]

            if missing_values:
                raise ValueError(
                    "Cannot calculate the PERC for reaction '{0}' because "
                    "values for {1} not defined.".format(
                        reaction.id, str(missing_values)))

            # Substitute values into rate equation
            rate_eq = rate_eq.subs(numerical_values)

            # Calculate rate equation and update with soluton for PERC
            sol = calculate_sol(flux, rate_eq, sym.Symbol(reaction.kf_str))
            percs_dict.update({reaction.kf_str: sol})
            if update_reactions:
                reaction.kf = percs_dict[reaction.kf_str]

        return percs_dict

    def build_model_from_string(self, model_str, verbose=True, fwd_arrow=None,
                                rev_arrow=None, reversible_arrow=None,
                                term_split="+", reaction_split=";",
                                reaction_id_split=":"):
        """Create a :class:`MassModel` from strings of reaction equations.

        Takes a string representation of the reactions and uses the
        specifications supplied in the optional arguments to first infer a set
        of reactions and their identifiers, then to infer metabolites,
        metabolite compartments, and stoichiometries for the reactions. It
        also infers the reversibility of the reaction from the reaction arrow.
        For example::

            '''
            ReactionID_1: S + E <=> ES
            ReactionID_2: ES -> E + P;
            ReactionID_3: E + I <=> EI
            '''

        Parameters
        ----------
        model : str
            A string representing the reaction formulas (equation) for the
            model.
        verbose : bool
            Setting the verbosity of the function.
        fwd_arrow: re.compile, None
            For forward irreversible reaction arrows. If ``None``, the
            arrow is expected to be ``'-->'`` or ``'==>'``.
        rev_arrow: re.compile, None
            For backward irreversible reaction arrows. If ``None``, the
            arrow is expected to be ``'<--'`` or ``'<=='``.
        reversible_arrow: re.compile, None
            For reversible reaction arrows. If ``None``, the arrow is expected
            to be ``'<=>'`` or ``'<->'``.
        term_split : str
            Dividing individual metabolite entries. Default is ``"+"``.
        reaction_split : str
            Dividing individual reaction entries. Default is ``";"``.
        reaction_id_split : str
            Dividing individual reaction entries from their identifiers.
            Default is ``":"``.

        See Also
        --------
        :meth:`.MassReaction.build_reaction_from_string`
            Base function for building reactions.

        """
        # Use the reaction split arguments to get the reactions and strip them
        reaction_list = [reaction_str.strip()
                         for reaction_str in model_str.split(reaction_split)
                         if reaction_str.strip()]

        # Iterate through reaction strings
        for orig_reaction_str in reaction_list:
            # Split the reaction ID from the reaction equation
            split = orig_reaction_str.split(reaction_id_split)
            try:
                if len(split) != 2:
                    raise ValueError("Could not parse '{0}' for the reaction "
                                     "ID and formula (equation).".format(
                                         orig_reaction_str))
                reaction_id, reaction_str = (s.strip() for s in split)
                # Cannot build reaction without an ID
                if not reaction_id:
                    raise ValueError("No reaction ID found in '{0}'"
                                     .format(orig_reaction_str))
                # Create a new reaction if one does not already exist.
                try:
                    reaction = self.reactions.get_by_id(reaction_id)
                except KeyError:
                    if verbose:
                        print("New reaction {0} created".format(reaction_id))
                    reaction = MassReaction(reaction_id)
                self.add_reactions(reaction)
                # Build the reaction from the reaction string
                reaction.build_reaction_from_string(
                    reaction_str, verbose=verbose, fwd_arrow=fwd_arrow,
                    rev_arrow=rev_arrow, reversible_arrow=reversible_arrow,
                    term_split=term_split)
            except ValueError as e:
                # Log reactions that could not be built.
                warnings.warn(
                    "Failed to build reaction '{0}' due to the following:\n"
                    "{1}".format(orig_reaction_str, str(e)))
                continue
        # Ensure all pointers are updated.
        self.repair(rebuild_index=True, rebuild_relationships=True)

    def update_parameters(self, parameters, verbose=True):
        """Update the parameters associated with the MassModel.

        Parameters can be the following:

            * :attr:`.MassReaction.forward_rate_constant` (
              :attr:`~.MassReaction.kf`)
            * :attr:`.MassReaction.reverse_rate_constant` (
              :attr:`~.MassReaction.kr`)
            * :attr:`.MassReaction.equilibrium_constant` (
              :attr:`~.MassReaction.Keq`)
            * :attr:`.MassReaction.steady_state_flux` (
              :attr:`~.MassReaction.v`)
            * :attr:`~MassModel.boundary_conditions`
            * :attr:`~MassModel.custom_parameters`

        Notes
        -----
        The reactions must already exist in the model in order to change
        associated parameters. Any identifiers that are not identifiers of
        standard reaction parameter or of any 'boundary metabolites' will be
        set as a custom parameter.

        Parameters
        ----------
        parameters : dict
            A dict containing the parameter identifiers as strings and their
            corresponding values to set in the model.
        verbose : bool
            If ``True`` then display the warnings that may be raised when
            setting reaction parameters. Default is ``True``.

        See Also
        --------
        :attr:`.MassReaction.all_parameter_ids`
            List of default reaction parameter identifiers.
        :attr:`MassModel.boundary_metabolites`
            List of 'boundary metabolites' found in the model.

        """
        if not isinstance(parameters, dict):
            raise TypeError("parameters must be a dict.")

        for key, value in iteritems(parameters):
            if not isinstance(key, string_types):
                raise TypeError(
                    "Keys must be strings. '{0}' not a string.".format(key))

        for key, value in iteritems(parameters):
            # Check the parameter type
            if key in self.boundary_metabolites:
                self.add_boundary_conditions({key: value})
            elif key.split("_", 1)[0] in ["kf", "Keq", "kr", "v"]:
                # See if the reaction exists and if none found, assume
                # parameter is a custom parameter
                try:
                    p_type, reaction = key.split("_", 1)
                    reaction = self.reactions.get_by_id(reaction)
                    with warnings.catch_warnings():
                        if not verbose:
                            warnings.simplefilter("ignore")
                        setattr(reaction, p_type, value)
                except (KeyError, ValueError):
                    self.custom_parameters.update({key: value})
            # If parameter not found, assume parameter is a custom parameter
            else:
                self.custom_parameters.update({key: value})

    def update_initial_conditions(self, initial_conditions):
        """Update the initial conditions of the model.

        Can also be used to update initial conditions of fixed metabolites to
        change the concentration value at which the metabolite is fixed.

        Notes
        -----
        The metabolite(s) must already exist in the model to set the initial
        conditions. Initial conditions for the metabolites are accessed through
        :attr:`.MassMetabolite.initial_condition`. If an initial condition for
        a metabolite already exists, it will be replaced.

        Parameters
        ----------
        initial_conditions : dict
            A dictionary where metabolites are the keys and the initial
            conditions are the values.

        """
        if not isinstance(initial_conditions, dict):
            raise TypeError("initial_conditions must be a dictionary.")
        for metabolite, ic_value in iteritems(initial_conditions):
            # Try getting metabolite object from model
            try:
                metabolite = self.metabolites.get_by_id(str(metabolite))
            except KeyError as e:
                warnings.warn("No metabolite found for {0}".format(str(e)))
                continue
            # Try setting the initial condition
            try:
                metabolite.initial_condition = ic_value
            except (TypeError, ValueError) as e:
                warnings.warn(
                    "Cannot set initial condition for {0} due to the "
                    "following: {1}".format(metabolite.id, str(e)))
                continue

    def update_custom_rates(self, custom_rates, custom_parameters=None):
        r"""Update the custom rates of the model.

        Parameters
        ----------
        custom_rates : dict
            A dictionary where :class:`.MassReaction`\ s or their string
            identifiers are the keys and the rates are the string
            representations of the custom rate expression.
        custom_parameters : dict
            A dictionary of custom parameters for the custom rates, where the
            key:value pairs are the strings representing the custom parameters
            and their numerical values. If a custom parameter already exists in
            the model, it will be updated.

        Notes
        -----
        The reaction(s) must already exist in the model to set the custom rate.

        See Also
        --------
        :meth:`add_custom_rate`

        """
        if custom_parameters is not None:
            if not isinstance(custom_parameters, dict):
                raise TypeError("custom_parameters must be a dict.")
            self.custom_parameters.update(custom_parameters)

        for reaction, custom_rate in iteritems(custom_rates):
            if not isinstance(reaction, MassReaction):
                try:
                    reaction = self.reactions.get_by_id(reaction)
                except KeyError as e:
                    warnings.warn("No reaction found for {0}".format(str(e)))
                    continue
            try:
                self.add_custom_rate(reaction, custom_rate=custom_rate)
            except sym.SympifyError:
                print(custom_rate, "uh oh")
                warnings.warn("Unable to sympify rate equation for "
                              "'{0}'.".format(reaction.id))

    def has_equivalent_odes(self, right, verbose=False):
        """Determine if :attr:`odes` between two models are equivalent.

        Notes
        -----
        The ODEs between two models are compared to determine whether the
        models can be considered equivalent, meaning that the models contain
        the same metabolites, reactions, and rate expressions such that they
        require the same set of parameters and initial conditions for
        simulation.

        Parameters
        ----------
        right : MassModel
            The :class:`MassModel` to compare to the left model (``self``).
        verbose : bool
            If ``True``, display the reason(s) for the differences in the left
            and right models. Default is ``False``.

        Returns
        -------
        bool
            Returns a bool indicating whether the model ODEs are equivalent.

        """
        equivalent = True
        # Determine whether ODE dicts have the same ODEs
        l_odes = {met.id: ode for met, ode in iteritems(self.odes)}
        r_odes = {met.id: ode for met, ode in iteritems(right.odes)}
        if l_odes != r_odes:
            equivalent = False

        if not equivalent and verbose:
            l_rates = {r.id: rate for r, rate in iteritems(self.rates)}
            r_rates = {r.id: rate for r, rate in iteritems(right.rates)}
            for i, (l_dict, r_dict) in enumerate(zip([l_odes, l_rates],
                                                     [r_odes, r_rates])):
                # Determine which metabolites do not exist in both models
                missing = set(l_dict)
                missing.symmetric_difference_update(set(r_dict))
                # Determine which metabolites have different ODEs
                diff_equations = set(
                    l_key for l_key, l_value in iteritems(l_dict)
                    if l_key in r_dict and l_value != r_dict[l_key])

                if i == 0:
                    msgs = ["Metabolites", "ODEs"]
                else:
                    msgs = ["Reactions", "rates"]

                msgs = ["{0} in one model only:".format(msgs[0]),
                        "{0} with different {1}: ".format(*msgs)]
                for item, msg in zip([missing, diff_equations], msgs):
                    if item:
                        warnings.warn(msg + str(sorted(list(item))))

        return equivalent

    # Internal
    def _mk_stoich_matrix(self, matrix_type=None, dtype=None,
                          update_model=True):
        """Return the stoichiometric matrix for a given MassModel.

        The rows represent the chemical species and the columns represent the
        reaction. S[i, j] therefore contains the quantity of species 'i'
        produced (positive) or consumed (negative) by reaction 'j'.

        Warnings
        --------
        This method is intended for internal use only. To safely update the
        stoichiometric matrix, use :meth:`~MassModel.update_S` instead.

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
        stoich_mat = convert_matrix(
            stoich_mat, matrix_type=matrix_type, dtype=dtype,
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

        Warnings
        --------
        This method is intended for internal use only. To safely update the
        stoichiometric matrix, use :meth:`~MassModel.update_S` instead.

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
        """Update the model stoichiometric matrix and its properties.

        Warnings
        --------
        This method is intended for internal use only. To safely update the
        stoichiometric matrix, use :meth:`~MassModel.update_S` instead.

        """
        self._S = stoich_mat
        self._matrix_type = matrix_type
        self._dtype = dtype

    def _get_all_parameters(self):
        """Get a dict containing all of defined model parameters in the model.

        Warnings
        --------
        This method is intended for internal use only.

        """
        return {str(param): value for subdict in itervalues(self.parameters)
                for param, value in iteritems(subdict)}

    def _copy_model_metabolites(self, new_model):
        """Copy the metabolites in creating a partial "deepcopy" of model.

        Warnings
        --------
        This method is intended for internal use only.

        """
        # Initialize DictList and set attributes to not copy by ref.
        new_model.metabolites = DictList()
        do_not_copy_by_ref = {"_reaction", "_model"}
        # Copy the metabolites
        for metabolite in self.metabolites:
            new_metabolite = metabolite.__class__()
            for attr, value in iteritems(metabolite.__dict__):
                if attr not in do_not_copy_by_ref:
                    new_metabolite.__dict__[attr] = copy(
                        value) if attr == "formula" else value
            new_metabolite._model = new_model
            new_model.metabolites.append(new_metabolite)

        return new_model

    def _copy_model_genes(self, new_model):
        """Copy the genes in creating a partial "deepcopy" of model.

        Warnings
        --------
        This method is intended for internal use only.

        """
        # Initialize DictList and set attributes to not copy by ref.
        new_model.genes = DictList()
        do_not_copy_by_ref = {"_reaction", "_model"}
        # Copy the genes
        for gene in self.genes:
            new_gene = gene.__class__(None)
            for attr, value in iteritems(gene.__dict__):
                if attr not in do_not_copy_by_ref:
                    new_gene.__dict__[attr] = copy(
                        value) if attr == "formula" else value
            new_gene._model = new_model
            new_model.genes.append(new_gene)

        return new_model

    def _copy_model_reactions(self, new_model):
        """Copy the reactions in creating a partial "deepcopy" of model.

        Warnings
        --------
        This method is intended for internal use only.

        """
        # Initialize DictList and set attributes to not copy by ref.
        new_model.reactions = DictList()
        do_not_copy_by_ref = {"_model", "_metabolites", "_genes"}
        # Copy the reactions
        for reaction in self.reactions:
            new_reaction = reaction.__class__()
            for attr, value in iteritems(reaction.__dict__):
                if attr not in do_not_copy_by_ref:
                    new_reaction.__dict__[attr] = copy(value)
            new_reaction._model = new_model
            new_model.reactions.append(new_reaction)
            # Update awareness for metabolites and genes
            for metabolite, stoich in iteritems(reaction._metabolites):
                new_metabolite = new_model.metabolites.get_by_id(metabolite.id)
                new_reaction._metabolites[new_metabolite] = stoich
                new_metabolite._reaction.add(new_reaction)
            for gene in reaction._genes:
                new_gene = new_model.genes.get_by_id(gene.id)
                new_reaction._genes.add(new_gene)
                new_gene._reaction.add(new_reaction)
            # Copy the custom rate for the reaction:
            if reaction in iterkeys(self.custom_rates):
                new_model.custom_rates.update({
                    new_reaction: self.custom_rates[reaction]})
            # Copy custom parameters if there are custom rates
        if self.custom_parameters:
            new_model.custom_parameters.update(self.custom_parameters)

        return new_model

    def _copy_model_enzyme_modules(self, new_model):
        """Copy the enzyme_modules in creating a partial "deepcopy" of model.

        Warnings
        --------
        This method is intended for internal use only.

        """
        new_model.enzyme_modules = DictList()
        do_not_copy_by_ref = {
            "ligands", "enzyme_forms", "enzyme_reactions",
            "categorized_ligands", "categorized_enzyme_forms",
            "categorized_enzyme_reactions", "model"}
        # Copy the enzyme_modules
        for enzyme in self.enzyme_modules:
            new_enzyme = enzyme.__class__()
            for attr, value in iteritems(enzyme):
                if attr not in do_not_copy_by_ref:
                    new_enzyme[attr] = copy(value)
            # Update associated model and object pointers for enzyme
            new_enzyme["model"] = new_model
            new_enzyme._update_object_pointers(new_model)
            new_model.enzyme_modules.append(new_enzyme)

        return new_model

    def _existing_obj_filter(self, attr, item_list):
        """Filter for existing items and return a new DictList.

        Warnings
        --------
        This method is intended for internal use only.

        """
        def existing_filter(item):
            if item.id in getattr(self, attr):
                LOGGER.warning("Ignoring %s '%s' since it already exists.",
                               attr[:-1], item.id)
                return False
            return True

        return DictList(filter(existing_filter, item_list))

    # TODO Fix when changes finished
    def _repr_html_(self):
        """HTML representation of the overview for the MassModel.

        Warnings
        --------
        This method is intended for internal use only.

        """
        try:
            dim_S = "{0}x{1}".format(self.S.shape[0], self.S.shape[1])
            rank = np.linalg.matrix_rank(self.S)
        except (np.linalg.LinAlgError, ValueError):
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
                    <td><strong>Number of Initial Conditions</strong></td>
                    <td>{num_ic}</td>
                </tr><tr>
                    <td><strong>Number of Fixed Metabolites</strong></td>
                    <td>{num_fixed}</td>
                </tr><tr>
                    <td><strong>Number of Reactions</strong></td>
                    <td>{num_reactions}</td>
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
                    <td><strong>Number of Boundary Reactions</strong></td>
                    <td>{num_boundary}</td>
                </tr><tr>
                    <td><strong>Number of Boundary Conditions</strong></td>
                    <td>{num_boundary_conditions}</td>
                </tr><tr>
                    <td><strong>Number of Custom Rates</strong></td>
                    <td>{num_custom_rates}</td>
                </tr><tr>
                    <td><strong>Number of Genes</strong></td>
                    <td>{num_genes}</td>
                </tr><tr>
                    <td><strong>Number of Enzymes</strong></td>
                    <td>{num_enzyme_modules}</td>
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
                   S_type="{0}, {1}".format(self._matrix_type,
                                            self._dtype.__name__),
                   num_metabolites=len(self.metabolites),
                   num_ic=len(self.initial_conditions),
                   num_fixed=len(self.fixed),
                   num_reactions=len(self.reactions),
                   num_kfs=len(self.parameters["kf"]),
                   num_Keqs=len(self.parameters["Keq"]),
                   num_irreversible=len(self.irreversible_reactions),
                   num_boundary=len(self.boundary),
                   num_boundary_conditions=len(self.boundary_conditions),
                   num_custom_rates=len(self.custom_rates),
                   num_genes=len(self.genes),
                   num_enzyme_modules=len(self.enzyme_modules),
                   compartments=", ".join(v if v else k for k, v in
                                          iteritems(self.compartments)),
                   units=", ".join([u.id for u in self.units]))

    # Dunders
    def __enter__(self):
        """Record all future changes for context management of the MassModel.

        Changes are undone when a call to __exit__ is received.

        Warnings
        --------
        This method is intended for internal use only.

        """
        # Create a new context and add it to the stack
        try:
            self._contexts.append(HistoryManager())
        except AttributeError:
            self._contexts = [HistoryManager()]
        return self

    def __exit__(self, exception_type, exception_value, traceback):
        """Pop the top context manager and trigger the undo functions.

        Warnings
        --------
        This method is intended for internal use only.

        """
        context = self._contexts.pop()
        context.reset()

    def __setstate__(self, state):
        """Ensure all objects in the model point to the MassModel.

        Warnings
        --------
        This method is intended for internal use only.

        """
        self.__dict__.update(state)
        for attr in ['reactions', 'metabolites', 'genes', 'enzyme_modules']:
            for x in getattr(self, attr):
                x._model = self
        if not hasattr(self, "name"):
            self.name = ""

    def __getstate__(self):
        """Get the state for serialization.

        Ensures that the context stack is cleared prior to serialization,
        since partial functions cannot be pickled reliably

        Warnings
        --------
        This method is intended for internal use only.

        """
        odict = self.__dict__.copy()
        odict['_contexts'] = []
        return odict


__all__ = ("MassModel", "LOGGER", "CHOPNSQ")
