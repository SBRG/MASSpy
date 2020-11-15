# -*- coding: utf-8 -*-
r"""
MassModel is a class for holding information regarding a :mod:`mass` model.

The :class:`MassModel` class inherits and extends the
:class:`~cobra.core.model.Model` class in :mod:`cobra`. It contains additional
information required for simulations and other :mod:`mass` functions and
workflows.

Some key differences between the
:class:`cobra.Model <cobra.core.model.Model>` and the
:class:`mass.MassModel <mass.core.mass_model.MassModel>` are
listed below:

    * When instantiating a :class:`MassModel` from a
      :class:`cobra.Model <cobra.core.model.Model>`, any associated
      :class:`cobra.Metabolite <cobra.core.metabolite.Metabolite>` will also
      be converted into a :class:`.MassMetabolite`, and any associated
      :class:`cobra.Reaction <cobra.core.reaction.Reaction>` will also
      be converted into a :class:`.MassReaction`.

      Additionally, any groups associated with the model containing
      :class:`cobra.Metabolites <cobra.core.metabolite.Metabolite>` or
      :class:`cobra.Reactions <cobra.core.reaction.Reaction>` will be updated
      with the correspondining :mod:`mass` objects.

    * No :mod:`cobra` object can be directly added to the :class:`MassModel`
      with the exception of :class:`~cobra.core.gene.Gene` and
      :class:`~cobra.core.group.Group` objects. Any :mod:`cobra` object will
      first be instantiated into its equivalent :mod:`mass` object before
      beimg added to the :class:`MassModel`.

"""
import re
import warnings
from copy import copy, deepcopy
from functools import partial

import numpy as np
import sympy as sym
from cobra.core.dictlist import DictList
from cobra.core.gene import Gene
from cobra.core.group import Group
from cobra.core.metabolite import Metabolite
from cobra.core.model import Model
from cobra.core.reaction import Reaction
from cobra.exceptions import SolverNotFound
from cobra.util.array import create_stoichiometric_matrix
from cobra.util.context import get_context
from cobra.util.util import format_long_string
from six import integer_types, iteritems, iterkeys, itervalues, string_types

from mass.core.mass_configuration import MassConfiguration
from mass.core.mass_metabolite import MassMetabolite
from mass.core.mass_reaction import MassReaction
from mass.core.units import UnitDefinition
from mass.util.expressions import create_custom_rate, strip_time
from mass.util.matrix import _get_matrix_constructor, convert_matrix, matrix_rank
from mass.util.util import (
    _check_kwargs,
    _make_logger,
    ensure_iterable,
    get_public_attributes_and_methods,
)


# Set the logger
LOGGER = _make_logger(__name__)
"""logging.Logger: Logger for :mod:`~mass.core.mass_model` submodule."""

CHOPNSQ = ["C", "H", "O", "P", "N", "S", "q"]
"""list: Contains the six most abundant elements and charge for molecules."""

MASSCONFIGURATION = MassConfiguration()


class MassModel(Model):
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
    array_type : str
        A string identifiying the desired format for the returned matrix.
        Valid matrix types include ``'dense'``, ``'dok'``, ``'lil'``,
        ``'DataFrame'``, and ``'symbolic'`` Default is ``'DataFrame'``.
        See the :mod:`~.matrix` module documentation for more information
        on the ``array_type``.
    dtype : data-type
        The desired array data-type for the stoichiometric matrix. If ``None``
        then the data-type will default to ``numpy.float64``.

    Attributes
    ----------
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
    groups : ~cobra.core.dictlist.DictList
        A :class:`~cobra.core.dictlist.DictList` where the keys are group
        identifiers and the values are the associated
        :class:`~cobra.core.group.Group`\ s.
    enzyme_modules : ~cobra.core.dictlist.DictList
        A :class:`~cobra.core.dictlist.DictList` where the keys are enzyme
        module identifiers and the values are the associated
        :class:`~.EnzymeModuleDict`\ s.
    custom_rates : dict
        A ``dict`` to store custom rate expressions for specific reactions,
        where the keys are :class:`~.MassReaction`\ s and values are the
        custom rate expressions given as :mod:`sympy` expressions. Custom rate
        expressions will always be prioritized over automatically generated
        mass action rates.
    custom_parameters : dict
        A ``dict`` to store the custom parameters for the custom rates,
        where key:value pairs are the string identifiers for the parameters
        and their corresponding numerical value.
    boundary_conditions : dict
        A ``dict`` to store boundary conditions, where keys are string
        identifiers for 'boundary metabolites' of boundary reactions, and
        values are the corresponding boundary condition numerical value or
        function of time. Note that boundary conditions are treated as
        parameters and NOT as species.
    units : ~cobra.core.dictlist.DictList
        :class:`~cobra.core.dictlist.DictList` of :class:`~.UnitDefinition`\ s
        to store in the model for referencing.

    Warnings
    --------
    * Note that the :class:`MassModel` does NOT track units, and it is
      therefore incumbent upon the user to maintain unit consistency the model.

    * Note that boundary conditions are considered parameters and NOT as
      species in a reaction.

    """

    def __init__(
        self, id_or_model=None, name=None, array_type="dense", dtype=np.float64
    ):
        """Initialize the MassModel."""
        # Instiantiate a new MassModel with state identical to
        # the provided Model or MassModel object via cobra model initialization
        # and its call to the __setstate__ method.
        super(MassModel, self).__init__(id_or_model, name)
        if isinstance(id_or_model, Model) and not isinstance(id_or_model, MassModel):
            # Turn cobra objects into mass objects and fix groups.
            self._cobra_to_mass_repair()

        if not isinstance(id_or_model, MassModel):
            # Initialize DictLists for storing enzyme modules and units.
            # Reactions, metabolites, genes, and groups are initialized with
            # the model
            self.enzyme_modules = DictList()
            self.units = DictList()
            # Initialize dictionaries for custom rates, custom parameters,
            # boundary conditions, compartments, and units.
            self.boundary_conditions = {}
            self.custom_rates = {}
            self.custom_parameters = {}

            # Store the stoichiometric matrix, its matrix type, and data type
            self._array_type = array_type
            self._dtype = dtype
            self._S = self._mk_stoich_matrix(
                array_type=self._array_type, dtype=self._dtype, update_model=True
            )

    # Public
    @property
    def stoichiometric_matrix(self):
        """Return the stoichiometric matrix."""
        return self.update_S(
            array_type=self._array_type, dtype=self._dtype, update_model=False
        )

    @property
    def S(self):
        """Alias for the :attr:`stoichiometric_matrix`."""
        return self.stoichiometric_matrix

    @property
    def ordinary_differential_equations(self):
        """Return a ``dict`` of ODEs for the metabolites."""
        return {met: met.ode for met in self.metabolites}

    @property
    def odes(self):
        """Alias for the :attr:`ordinary_differential_equations`."""
        return self.ordinary_differential_equations

    @property
    def initial_conditions(self):
        r"""Get ``dict`` of :attr:`.MassMetabolite.initial_condition`\ s."""
        return {
            met: met.initial_condition
            for met in self.metabolites
            if met.initial_condition is not None
        }

    @property
    def ics(self):
        """Alias for the :attr:`initial_conditions`."""
        return self.initial_conditions

    @property
    def fixed(self):
        """Return a ``dict`` of all metabolite fixed conditions."""
        return {met: ic for met, ic in iteritems(self.initial_conditions) if met.fixed}

    @property
    def rates(self):
        """Return a ``dict`` of reaction rate expressions.

        If a reaction has an associated custom rate expression, the custom rate
        will be prioritized and returned in the ``dict`` instead of the
        automatically generated rate law expression.
        """
        return self.get_rate_expressions(self.reactions, rate_type=0)

    @property
    def steady_state_fluxes(self):
        """Return a ``dict`` of all reaction steady state fluxes."""
        return {
            rxn: rxn.steady_state_flux
            for rxn in self.reactions
            if rxn.steady_state_flux is not None
        }

    @property
    def v(self):
        """Alias for the :attr:`steady_state_fluxes`."""
        return self.steady_state_fluxes

    @property
    def boundary(self):
        """Return a ``list`` of boundary reactions in the model."""
        return super(MassModel, self).boundary

    @property
    def boundary_metabolites(self):
        """Return a sorted ``list`` of all 'boundary metabolites' in the model.

        See Also
        --------
        :attr:`.MassReaction.boundary_metabolite`

        """
        return sorted(
            list(set(rxn.boundary_metabolite for rxn in self.reactions if rxn.boundary))
        )

    @property
    def exchanges(self):
        """Return exchange reactions in the model.

        Exchange reactions are reactions that exchange mass with the exterior.
        Uses annotations and heuristics to exclude non-exchanges such as sink
        and demand reactions.
        """
        return super(MassModel, self).exchanges

    @property
    def demands(self):
        """Return demand reactions in the model.

        Demands are irreversible reactions that accumulate or consume a
        metabolite in the inside of the model.
        """
        return super(MassModel, self).demands

    @property
    def sinks(self):
        """Return sink reactions in the model.

        Sinks are reversible reactions that accumulate or consume a metabolite
        in the inside of the model.
        """
        return super(MassModel, self).sinks

    @property
    def irreversible_reactions(self):
        """Return a ``list`` of all irreversible reactions in the model."""
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
        parameters.update(
            {
                "v": {
                    rxn.flux_symbol_str: flux
                    for rxn, flux in iteritems(self.steady_state_fluxes)
                }
            }
        )
        parameters.update({"Custom": self.custom_parameters})
        parameters.update({"Boundary": self.boundary_conditions})

        return parameters

    @property
    def compartments(self):
        """Get or set a ``dict`` of all metabolite compartments.

        Assigning a ``dict`` to this property updates the model's
        ``dict`` of compartment descriptions with the new values.

        Notes
        -----
        * Setter extends :meth:`~cobra.core.model.Model.compartments` of
          the :class:`cobra.Model <cobra.core.model.Model>` to enable
          resetting the attribute by setting an empty ``dict``

        Parameters
        ----------
        compartment_dict : dict
            A ``dict`` mapping compartments abbreviations to full names.
            An empty ``dict`` will reset the compartments.

        """
        return super(MassModel, self).compartments

    @compartments.setter
    def compartments(self, compartment_dict):
        """Set the ``dict`` of current compartment descriptions."""
        if compartment_dict:
            super(MassModel, self.__class__).compartments.fset(self, compartment_dict)
        else:
            setattr(self, "_compartments", {})

    @property
    def conc_solver(self):
        """Return the :class:`.ConcSolver` associated with the model."""
        if hasattr(self, "_conc_solver"):
            return getattr(self, "_conc_solver")

        raise SolverNotFound("No ConcSolver is associated with this model.")

    def update_S(self, array_type=None, dtype=None, update_model=True):
        r"""Update the stoichiometric matrix of the model.

        Parameters
        ----------
        array_type : str
            A string identifiying the desired format for the returned matrix.
            Valid matrix types include ``'dense'``, ``'dok'``, ``'lil'``,
            ``'DataFrame'``, and ``'symbolic'``
            Default is the current ``array_type``. See the :mod:`~.matrix`
            module documentation for more information on the ``array_type``.
        dtype : data-type
            The desired array data-type for the stoichiometric matrix.
            If ``None`` then the data-type will default to the
            current ``dtype``.
        update_model : bool
            If ``True``, will update the stored stoichiometric matrix,
            the matrix type, and the data-type for the model.

        Returns
        -------
        matrix of type ``array_type``
            The stoichiometric matrix for the :class:`~.MassModel` returned
            as the given ``array_type`` and with a data-type of ``dtype``.

        """
        # Use the model's stored matrix type if the matrix-type is not given.
        if array_type is None:
            array_type = self._array_type
        # Use the model's stored data-type if the data-type is not specified.
        if dtype is None:
            dtype = self._dtype

        # Check input of update_model
        if not isinstance(update_model, bool):
            raise TypeError("update_model must be a bool.")

        # If a matrix has not been constructed yet, or if there are no changes
        # to the reactions, return a newly constructed stoichiometric matrix.
        stoich_mat = self._mk_stoich_matrix(
            array_type=array_type, dtype=dtype, update_model=update_model
        )
        # Internally update the model if desired
        if update_model:
            self._S = stoich_mat
            self._array_type = array_type
            self._dtype = dtype

        return stoich_mat

    def add_metabolites(self, metabolite_list):
        r"""Add a ``list`` of metabolites to the model.

        The change is reverted upon exit when using the :class:`~.MassModel`
        as a context.

        Parameters
        ----------
        metabolite_list : list
            A ``list`` containing :class:`~.MassMetabolite`\ s to add to
            the model.

        """
        super(MassModel, self).add_metabolites(metabolite_list)

    def remove_metabolites(self, metabolite_list, destructive=False):
        r"""Remove a ``list`` of metabolites from the model.

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
        super(MassModel, self).remove_metabolites(metabolite_list, destructive)

    def add_boundary_conditions(self, boundary_conditions):
        """Add boundary conditions values for the given boundary metabolites.

        Boundary condition values can be a numerical value, or they can be a
        string or :mod:`sympy` expression representing a function of time.
        The function must only depend on time.

        Parameters
        ----------
        boundary_conditions : dict
            A ``dict`` of boundary conditions containing the
            'boundary metabolites' and their corresponding value. The string
            representing the 'boundary_metabolite' must exist in the ``list``
            returned by :attr:`MassModel.boundary_metabolites`.

        See Also
        --------
        :attr:`boundary_metabolites`

        """
        if not isinstance(boundary_conditions, dict):
            raise TypeError("boundary_conditions must be a dict.")

        boundary_conditions_to_set = boundary_conditions.copy()
        for bound_met, bound_cond in iteritems(boundary_conditions_to_set):
            if (
                bound_met not in self.boundary_metabolites
                and bound_met not in self.metabolites
            ):
                raise ValueError(
                    "Did not find {0} in model metabolites or in "
                    "boundary reactions.".format(bound_met)
                )
            if bound_cond in ["", None]:
                boundary_conditions_to_set[bound_met] = None
            # Boundary condition is a function
            elif isinstance(bound_cond, (sym.Basic, string_types)):
                if isinstance(bound_cond, string_types):
                    bound_cond = sym.sympify(bound_cond)
                boundary_conditions_to_set[bound_met] = bound_cond
            # Boundary condition is an integer or float
            elif isinstance(bound_cond, (integer_types, float)):
                boundary_conditions_to_set[bound_met] = float(bound_cond)
            else:
                raise TypeError(
                    "Invalid boundary value for '{0}'. Boundary conditions can"
                    "only be numerical values or functions of time.".format(bound_met)
                )
        # Keep track of existing concentrations for context management.
        context = get_context(self)
        if context:
            existing_concs = {
                bound_met: self.boundary_conditions[bound_met]
                for bound_met in boundary_conditions
                if bound_met in self.boundary_conditions
            }

        self.boundary_conditions.update(boundary_conditions_to_set)

        if context:
            context(
                partial(self.boundary_conditions.pop, key)
                for key in iterkeys(boundary_conditions)
                if key not in existing_concs
            )
            context(partial(self.boundary_conditions.update, existing_concs))

    def remove_boundary_conditions(self, boundary_metabolite_list):
        """Remove the boundary condition for a list of `boundary metabolites`.

        Parameters
        ----------
        metabolite_list : list
            A ``list`` of metabolites to remove the boundary conditions for.
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
            bound_met
            for bound_met in boundary_metabolite_list
            if bound_met in self.boundary_conditions
        ]

        # Keep track of existing concentrations for context management.
        context = get_context(self)
        if context:
            existing_concs = {
                bound_met: self.boundary_conditions[bound_met]
                for bound_met in boundary_metabolite_list
                if bound_met in self.boundary_conditions
            }
        # Remove the initial conditions
        for bound_met in boundary_metabolite_list:
            del self.boundary_conditions[bound_met]

        if context:
            context(partial(self.boundary_conditions.update, existing_concs))

    def add_reactions(self, reaction_list):
        r"""Add reactions to the model.

        The change is reverted upon exit when using the :class:`~.MassModel`
        as a context.

        Notes
        -----
        * :class:`~.MassReaction`\ s with identifiers identical to an
          existing reaction are ignored.

        * Extends :meth:`cobra.core.model.Model.add_reactions` by
          first ensuring that the reactions to be added are
          :class:`.MassReactions`\ s and not
          :class:`cobra.Reactions <cobra.core.reaction.Reaction>`.
          and error message raised reflects the :mod:`mass` object.

        * If a :class:`cobra.Reaction <cobra.core.reaction.Reaction>` is
          provided. a warning is raised and a :class:`.MassReaction`
          will be instantiated using the
          :class:`cobra.Reaction <cobra.core.reaction.Reaction>`.

        Parameters
        ----------
        reaction_list : list
            A ``list`` of :class:`~.MassReaction`\ s to add to the model.

        """
        reaction_list = ensure_iterable(reaction_list)
        for i, reaction in enumerate(reaction_list):
            if isinstance(reaction, MassReaction):
                # No need to change MassReaction objects
                continue
            elif isinstance(reaction, Reaction):
                # Convert reaction to a MassReaction and raise a warning
                warnings.warn(
                    "'{0}' is not a mass.MassReaction, therefore "
                    "converting reaction before adding.".format(reaction.id)
                )
                mass_reaction = MassReaction(reaction)
                reaction_list[i] = mass_reaction
            else:
                # Input not recognized.
                raise TypeError("Unrecognized input {0}".format(str(reaction)))
        super(MassModel, self).add_reactions(reaction_list)

    def remove_reactions(self, reactions, remove_orphans=False):
        r"""Remove reactions from the model.

        The change is reverted upon exit when using the :class:`~.MassModel`
        as a context.

        Notes
        -----
        Extends :meth:`cobra.core.model.Model.remove_reactions` by also
        removing any custom rates along with the reaction (and custom
        parameters if ``remove_orphans=True``). Also removes
        the boundary condition if the reaction is a boundary reaction with
        a defined boundary condition.

        Parameters
        ----------
        reaction_list : list
            A ``list`` of :class:`~.MassReaction`\ s to be removed
            from the model.
        remove_orphans : bool
            If ``True``, will also remove orphaned genes and metabolites from
            the model. If a custom rate is removed, the orphaned custom
            parameters will also be removed.

        """
        super(MassModel, self).remove_reactions(reactions, remove_orphans)
        reactions = ensure_iterable(reactions)
        # Remove reaction from the custom rates and
        # the orphaned custom parameters if set
        for reaction in reactions:
            if reaction in self.custom_rates:
                self.remove_custom_rate(reaction, remove_orphans=remove_orphans)
            if (
                reaction.boundary
                and reaction.boundary_metabolite in self.boundary_conditions
            ):
                self.remove_boundary_conditions([reaction.boundary_metabolite])

    def add_boundary(
        self,
        metabolite,
        boundary_type="exchange",
        reaction_id=None,
        boundary_condition=0.0,
        **kwargs
    ):
        """Add a boundary reaction for a given metabolite.

        Accepted ``kwargs`` are passed to the underlying function for boundary
        reaction creation,
        :func:`cobra.Model.add_boundary <cobra.core.model.Model.add_boundary>`,
        and initialization of the :class:`.MassReaction`.

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
        * Extends :func:`cobra.core.model.Model.add_boundary` by allowing
          for metabolite identifers of existing metabolites in the model,
          boundary conditions and reaction subsystem to be set, and
          utilizes default bounds from the :class:`.MassConfiguration`. for
          creation of custom boundary reaction types.

        * To set the reaction ``boundary_type`` to something else, the desired
          identifier of the created reaction must be specified. The name will
          be given by the metabolite name and the given ``boundary_type``, and
          the reaction will be set its reversible attribute to ``True``.

          Bounds will be set to the defaults specified in the
          :class:`.MassConfiguration`.

        Parameters
        ----------
        metabolite : MassMetabolite or str
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
        boundary_condition : float, str, ~sympy.core.basic.Basic
            The boundary condition value to set. Must be an ``int``,
            ``float``, or a :mod:`sympy` expression dependent only on time.
            Default value is 0.
        **kwargs
            subsystem :
                ``str`` for subsystem where the reaction is meant to occur.
            lb :
                ``float`` for the lower bound of the resulting reaction, or
                ``None`` to use the default specified in the
                :class:`.MassConfiguration`.
            ub :
                ``float`` for the upper bound of the resulting reaction, or
                ``None`` to use the default specified in the
                :class:`.MassConfiguration`.
            sbo_term :
                ``str`` for the SBO term. A correct SBO term is set for the
                available boundary types. If a custom boundary type is chosen,
                a suitable SBO term should also be set.

        Returns
        -------
        MassReaction
            The :class:`~.MassReaction` of the new boundary reaction.

        """
        kwargs = _check_kwargs(
            {
                "subsystem": "",
                "lb": MASSCONFIGURATION.lower_bound,
                "ub": MASSCONFIGURATION.upper_bound,
                "sbo_term": None,
            },
            kwargs,
        )
        # Check for existance of the metabolite
        try:
            metabolite = self.metabolites.get_by_id(str(metabolite))
        except KeyError as e:
            raise ValueError("'{0}' does not exist in the model.".format(str(e)))

        # Create boundary reaction
        with warnings.catch_warnings():
            warnings.filterwarnings(
                "ignore",
                ".*is not a mass.MassReaction, therefore " "converting reaction.*",
            )
            reaction = super(MassModel, self).add_boundary(
                metabolite,
                type=boundary_type,
                reaction_id=reaction_id,
                lb=kwargs.get("lb", None),
                ub=kwargs.get("ub", None),
                sbo_term=kwargs.get("sbo_term", None),
            )
        # Get the added MassReaction
        reaction = self.reactions.get_by_id(reaction.id)
        # Add subsystem or boundary condition
        reaction.subsystem = kwargs.get("subsystem")
        self.add_boundary_conditions({reaction.boundary_metabolite: boundary_condition})

        if boundary_type == "demand":
            reaction.reversible = False

        return reaction

    def get_rate_expressions(
        self, reaction_list=None, rate_type=0, update_reactions=False
    ):
        r"""Get the rate expressions for a ``list`` of reactions in the model.

        Notes
        -----
        If a reaction has a custom rate in the :attr:`MassModel.custom_rates`
        attribute, it will be returned only when the ``rate_type=0``.

        Parameters
        ----------
        reaction_list : list
            A ``list`` of :class:`~.MassReaction`\ s to get the rate
            expressions for. Reactions must already exist in the model.
            If ``None``, then return the rates for all reactions in the model.
        rate_type : int
            The type of rate law to return. Must be 0, 1, 2, or 3.

                * If 0, the currrent rate expression is returned.
                * Type 1 will utilize the :attr:`forward_rate_constant` and the
                  :attr:`equilibrium_constant`.
                * Type 2 will utilize the :attr:`forward_rate_constant` and the
                  :attr:`reverse_rate_constant`.
                * Type 3 will utilize the :attr:`equilibrium_constant` and the
                  :attr:`reverse_rate_constant`.

            Default is ``0.``
        update_reactions : bool
            If ``True``, update the :class:`~.MassReaction` rate in addition
            to returning the rate expressions. Will not remove a custom
            rate.

        Returns
        -------
        dict
            A ``dict`` of reaction rates where keys are the reaction ids
            and values are the rate law expressions.

        """
        # Check the inputs
        if rate_type not in {0, 1, 2, 3}:
            raise ValueError("rate_type must be 0, 1, 2, or 3")

        # Use the MassModel reactions if no reaction list is given
        if reaction_list is None:
            reaction_list = self.reactions
        # Ensure list is iterable.
        reaction_list = ensure_iterable(reaction_list)

        if rate_type == 0:
            rate_dict = {rxn: rxn.rate for rxn in reaction_list}

        else:
            rate_dict = {
                rxn: rxn.get_mass_action_rate(rate_type, update_reactions)
                for rxn in reaction_list
            }

        return rate_dict

    def get_mass_action_ratios(self, reaction_list=None, sympy_expr=True):
        r"""Get mass action ratios for a ``list`` of reactions in the model.

        Parameters
        ----------
        reaction_list : list
            A ``list`` of :class:`~.MassReaction`\ s to get the mass action
            ratios for. Reactions must already exist in the model.
            If ``None``, then return the ratios for all reactions in the model.
        sympy_expr : bool
            If ``True`` then return the mass action ratios as a :mod:`sympy`
            expression, otherwise return the ratio as a human readable string.

        Returns
        -------
        dict
            A ``dict`` of mass action ratios where keys are the
            reaction ids and values are the ratios.

        """
        # Use the MassModel reactions if no reaction list is given
        if reaction_list is None:
            reaction_list = self.reactions
        # Ensure list is iterable.
        reaction_list = ensure_iterable(reaction_list)

        ratio_dict = dict(
            (rxn, rxn.get_mass_action_ratio())
            if sympy_expr
            else (rxn, str(rxn.get_disequilibrium_ratio()))
            for rxn in reaction_list
        )
        return ratio_dict

    def get_disequilibrium_ratios(self, reaction_list=None, sympy_expr=True):
        r"""Get disequilibrium ratios for a ``list`` of reactions in the model.

        Parameters
        ----------
        reaction_list : list
            A ``list`` of :class:`~.MassReaction`\ s to get the disequilibrium
            ratios for. Reactions must already exist in the model.
            If ``None``, then return the ratios for all reactions in the model.
        sympy_expr : bool
            If ``True`` then return the disequilibrium ratios as a :mod:`sympy`
            expression, otherwise return the ratio as a human readable string.

        Returns
        -------
        dict
            A ``dict`` of mass action ratios where keys are the
            reaction ids and values are the ratios.

        """
        # Use the MassModel reactions if no reaction list is given
        if reaction_list is None:
            reaction_list = self.reactions
        # Ensure list is iterable.
        reaction_list = ensure_iterable(reaction_list)

        ratio_dict = dict(
            (rxn, rxn.get_disequilibrium_ratio())
            if sympy_expr
            else (rxn, str(rxn.get_disequilibrium_ratio()))
            for rxn in reaction_list
        )
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
            A ``dict`` of custom parameters for the custom rate where the
            key:value pairs are the strings representing the custom parameters
            and their numerical values. The string representation of the custom
            parametes will be used to create the symbols needed for the sympy
            expression of the custom rate. If ``None``, then parameters are
            assumed to already exist in the model.

        See Also
        --------
        ~MassReaction.all_parameter_ids
            Lists the default reaction parameters automatically accounted for.

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
                if (
                    re.search(custom_parameter, custom_rate)
                    and custom_parameter not in custom_parameter_list
                ):
                    custom_parameter_list.append(custom_parameter)
        # Create the custom rate expression
        custom_rate = create_custom_rate(reaction, custom_rate, custom_parameter_list)
        self.custom_rates.update({reaction: custom_rate})
        self.custom_parameters.update(custom_parameters)

        context = get_context(self)
        if context:
            context(partial(self.custom_rates.pop, reaction))
            context(
                partial(
                    (self.custom_parameters.pop, key)
                    for key in custom_parameter_list
                    if key in iterkeys(self.custom_parameters)
                )
            )
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
            warnings.warn(
                "Did not find a custom custom rate expression "
                "associated with reaction {0}.".format(reaction.id)
            )
            return
        # Remove the rate
        del self.custom_rates[reaction]

        # Remove orphaned custom parameters if desired.
        args = rate_to_remove.atoms(sym.Symbol)
        standards = [sym.Symbol(arg) for arg in reaction.all_parameter_ids + ["t"]]
        # Save currently existing parameters for context management if needed.
        existing = {
            str(arg): self.custom_parameters[str(arg)]
            for arg in args
            if arg in self.custom_parameters and arg not in standards
        }

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
            context(partial(self.custom_rates.update, {reaction: rate_to_remove}))
            if remove_orphans:
                context(partial(self.custom_parameters.update, existing))

    def reset_custom_rates(self):
        """Reset all custom rate expressions and parameters in a model.

        The change is reverted upon exit when using the :class:`~.MassModel`
        as a context.

        Warnings
        --------
        Using this method will remove all custom rates and custom rate
        parameters in the model. To remove a specific rate without affecting
        the other custom rates or parameters, use :meth:`remove_custom_rate`
        instead.

        """
        context = get_context(self)
        if context:
            existing_customs = self.custom_rates.copy()
            existing_parameters = self.parameters.copy()

        self.custom_rates = {}
        self.custom_parameters = {}
        LOGGER.info("All custom rate expressions and parameters have been reset")

        if context:
            context(partial(self.custom_parameters.update, existing_parameters))
            context(partial(self.custom_parameters.update, existing_customs))

    def add_units(self, unit_defs):
        r"""Add a :class:`~.UnitDefinition` to the model :attr:`units`.

        The change is reverted upon exit when using the :class:`~.MassModel`
        as a context.

        Notes
        -----
        The model will not automatically track or convert units. Therefore,
        it is up to the user to ensure unit consistency in the model.

        Parameters
        ----------
        unit_defs : list
            A ``list`` of :class:`~.UnitDefinition`\ s to add to the model.

        """
        # Ensure iterable input and units are valid Unit objects.
        unit_defs = DictList(ensure_iterable(unit_defs))
        for unit in list(unit_defs):
            if not isinstance(unit, UnitDefinition):
                raise ValueError(
                    "'{0}' is not a valid UnitDefinition.".format(str(unit))
                )
            # Skip existing units.
            if unit.id in self.units.list_attr("id"):
                warnings.warn(
                    "Skipping '{0}' for it already exists in the" " model.".format(unit)
                )
                unit_defs.remove(unit)
        # Add new unit definitions to the units attribute
        self.units += unit_defs
        context = get_context(self)
        if context:
            context(partial(self.units.__isub__, unit_defs))

    def remove_units(self, unit_defs):
        r"""Remove a :class:`~.UnitDefinition` from the model :attr:`units`.

        The change is reverted upon exit when using the :class:`~.MassModel`
        as a context.

        Notes
        -----
        The model will not automatically track or convert units. Therefore,
        it is up to the user to ensure unit consistency in the model.

        Parameters
        ----------
        unit_defs : list
            A ``list`` of :class:`~.UnitDefinition`\ s or their string
            identifiers to remove from the model.

        """
        unit_defs = ensure_iterable(unit_defs)
        # Iteratre through units, raise ValueError if unit does not exist.
        for unit in unit_defs:
            try:
                unit = self.units.get_by_id(unit)
            except KeyError as e:
                raise ValueError("'{0}' does not exist in the model.".format(str(e)))
        # Remove unit definitions to the units attribute
        self.units -= unit_defs
        context = get_context(self)
        if context:
            context(partial(self.units.__iadd__, unit_defs))

    def reset_units(self):
        r"""Reset all unit definitions in a model.

        The change is reverted upon exit when using the :class:`~.MassModel`
        as a context.

        Warnings
        --------
        Using this method will remove all :class:`~.UnitDefinition`\ s from
        the model. To remove a :class:`~.UnitDefinition` without affecting
        other units, use :meth:`remove_units` instead.

        """
        context = get_context(self)
        if context:
            existing_units = list(self.units)

        self.units = DictList()
        LOGGER.info("All unit definitions have been reset.")

        if context:
            context(partial(self.units.__iadd__, existing_units))

    def get_elemental_matrix(self, array_type=None, dtype=None):
        """Get the elemental matrix for a model.

        Parameters
        ----------
        array_type : str
            A string identifiying the desired format for the returned matrix.
            Valid matrix types include ``'dense'``, ``'dok'``, ``'lil'``,
            ``'DataFrame'``, and ``'symbolic'``
            Default is ``'dense'``. See the :mod:`~.matrix` module
            documentation for more information on the ``array_type``.
        dtype : data-type
            The desired array data-type for the matrix. If ``None`` then
            the data-type will default to ``numpy.float64``.

        Returns
        -------
        matrix of type ``array_type``
            The elemntal matrix for the :class:`~.MassModel` returned
            as the given ``array_type`` and with a data-type of ``dtype``.

        """
        # Set up for matrix construction if matrix types are correct.
        (matrix_constructor, array_type, dtype) = _get_matrix_constructor(
            array_type=array_type, dtype=dtype
        )

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
            moieties.extend(
                [
                    element
                    for element in iterkeys(element_dict)
                    if "[" in element and "]" in element and element not in moieties
                ]
            )

        row_ids = CHOPNSQ.copy()
        # Add additional moieties to the elemental matrix
        if moieties:
            moiety_mat = matrix_constructor((len(moieties), len(self.metabolites)))
            # Create additional matrix for moieties
            for met in self.metabolites:
                element_dict = met.elements
                for element in moieties:
                    if element in element_dict:
                        amount = 1
                    else:
                        amount = 0
                    moiety_mat[moieties.index(element), m_ind(met)] = amount
            # Concatenate matrices
            elem_mat = np.concatenate((elem_mat, moiety_mat), axis=0)
            row_ids.extend(moieties)

        # Convert matrix to a dataframe if matrix type is a dataframe
        elem_mat = convert_matrix(
            elem_mat,
            array_type=array_type,
            dtype=dtype,
            row_ids=row_ids,
            col_ids=[m.id for m in self.metabolites],
        )

        return elem_mat

    def get_elemental_charge_balancing(self, array_type=None, dtype=None):
        """Get the elemental charge balance as a matrix for a model.

        Parameters
        ----------
        array_type : str
            A string identifiying the desired format for the returned matrix.
            Valid matrix types include ``'dense'``, ``'dok'``, ``'lil'``,
            ``'DataFrame'``, and ``'symbolic'``
            Default is ``'dense'``. See the :mod:`~.matrix` module
            documentation for more information on the ``array_type``.
        dtype : data-type
            The desired array data-type for the matrix. If ``None`` then
            the data-type will default to ``numpy.float64``.

        Returns
        -------
        matrix of type ``array_type``
            The charge balancing matrix for the :class:`~.MassModel` returned
            as the given ``array_type`` and with a data-type of ``dtype``.

        """
        elem_mat = self.get_elemental_matrix(array_type="DataFrame")
        row_ids = elem_mat.index

        stoich_mat = self.update_S(array_type="dense", update_model=False)
        charge_mat = np.array(elem_mat).dot(stoich_mat)

        if array_type is None:
            array_type = "dense"
        if dtype is None:
            dtype = np.float64

        charge_mat = convert_matrix(
            charge_mat,
            array_type=array_type,
            dtype=dtype,
            row_ids=row_ids,
            col_ids=[r.id for r in self.reactions],
        )
        return charge_mat

    def repair(self, rebuild_index=True, rebuild_relationships=True):
        """Update all indicies and pointers in the model.

        Notes
        -----
        Extends :meth:`~cobra.core.model.Model.repair` of the
        :class:`cobra.Model <cobra.core.model.Model>` to include the
        :attr:`MassModel.enzyme_modules` and :attr:`MassModel.units`.

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
        super(MassModel, self).repair(rebuild_index, rebuild_relationships)
        # Rebuild the DictList indices
        if rebuild_index:
            self.enzyme_modules._generate_index()
            self.units._generate_index()

        # Make all objects point to the model.
        for e in self.enzyme_modules:
            e._model = self

    def copy(self):
        r"""Create a partial "deepcopy" of the :class:`MassModel`.

        All of the :class:`~.MassMetabolite`\ s, :class:`~.MassReaction`\ s,
        :class:`~cobra.core.gene.Gene`\ s and :class:`~.EnzymeModuleDict`\ s,
        the boundary conditions, custom_rates, custom_parameters, and the
        stoichiometric matrix are created anew, but in a faster fashion
        than ``deepcopy``.

        Notes
        -----
        Overrides :meth:`~cobra.core.model.Model.copy` of the
        :class:`cobra.Model <cobra.core.model.Model>` so that all objects
        are :mod:`mass` objects and additional attributes of specific to the
        :class:`MassModel`.

        """
        # Define a new model
        new_model = self.__class__()
        # Define items that will not be copied by their references
        do_not_copy_by_ref = [
            "metabolites",
            "reactions",
            "genes",
            "enzyme_modules",
            "groups",
            "_S",
            "custom_rates",
            "units",
            "boundary_conditions",
            "custom_parameters",
            "notes",
            "annotation",
        ]
        for attr in self.__dict__:
            if attr not in do_not_copy_by_ref:
                new_model.__dict__[attr] = self.__dict__[attr]

        for attr in do_not_copy_by_ref[-5:]:
            setattr(new_model, attr, deepcopy(getattr(self, attr)))

        # Copy the metabolites
        new_model.metabolites += self._copy_model_metabolites(new_model)
        # Copy the genes
        new_model.genes += self._copy_model_genes(new_model)
        # Copy the reactions and rates (including custom rates)
        new_model.reactions += self._copy_model_reactions(new_model)
        # Copy the custom rate for the reaction:
        if self.custom_rates:
            new_model.custom_rates.update(
                {
                    new_model.reactions.get_by_id(reaction.id): custom_rate
                    for reaction, custom_rate in iteritems(self.custom_rates)
                }
            )
        # Copy custom parameters
        if self.custom_parameters:
            new_model.custom_parameters.update(self.custom_parameters)

        # Copy any existing groups
        new_model.groups += self._copy_model_groups(new_model)

        # Copy any existing enzyme_modules
        new_model.enzyme_modules += self._copy_model_enzyme_modules(new_model)
        # Create the new stoichiometric matrix for the model.
        new_model._S = self._mk_stoich_matrix(
            array_type=self._array_type, dtype=self._dtype, update_model=True
        )

        solvers_dict = {"_solver": self.solver}
        if hasattr(self, "_conc_solver"):
            solvers_dict["_conc_solver"] = self.conc_solver

        for attr, solver in iteritems(solvers_dict):
            try:
                setattr(new_model, attr, deepcopy(solver))
                # Cplex has an issue with deep copies
            except Exception:
                setattr(new_model, attr, copy(solver))

        # Doesn't make sense to retain the context of a copied model so
        # assign a new empty context
        new_model._contexts = []

        return new_model

    def merge(self, right, prefix_existing=None, inplace=True, objective="left"):
        """Merge two models into one model with the objects from both.

        The reactions, metabolites, genes, enzyme modules, boundary conditions,
        custom rate expressions, rate parameters, compartments, units, notes,
        and annotations from the right model are also copied to left model.
        However, note that in cases where identifiers for objects are identical
        or a ``dict`` item has an identical key(s), priority will be given
        to what already exists in the left model.

        Custom constraints and variables from right models are also copied
        to left model, however note that, constraints and variables are
        assumed to be the same if they have the same name.

        Notes
        -----
        * When merging an :class:`.~EnzymeModule` into a :class:`MassModel`,
          the enzyme module is converted to an :class:`.~EnzymeModuleDict` and
          stored in a :class:`~cobra.core.dictlist.DictList` accessible via the
          :attr:`enzyme_modules` attribute. If an :class:`.~EnzymeModuleDict`
          already exists in the model, it will be replaced.

        * Extends :meth:`~cobra.core.model.Model.merge` of the
          :class:`cobra.Model <cobra.core.model.Model>` by including
          additional :class:`MassModel` attributes.

        Parameters
        ----------
        right : MassModel
            The model to merge into the left model.
        prefix_existing : str
            If provided, the string is used to prefix the reaction identifier
            of a reaction in the second model if that reaction already exists
            within the first model. Will also apply prefix to identifiers
            of enzyme modules in the second model.
        inplace : bool
            If ``True`` then add reactions from second (right) model directly
            to the first (left) model. Otherwise, create a new model leaving
            the left model untouched. When done within the model as context,
            changes to the models are reverted upon exit.
        objective : str
            One of ``"left"``, ``"right"`` or ``"sum"`` for setting the
            objective of the resulting model to that of the corresponding
            model or the sum of both. Default is ``"left"``.

        Returns
        -------
        MassModel
            A new :class:`MassModel` or ``self`` representing the
            merged model.

        """
        # Check whether two MassModels are being merged,
        # or if a MassModel and a MassModel subclass are being merged.
        if right.__class__ is not MassModel and issubclass(right.__class__, MassModel):
            return right._add_self_to_model(self, prefix_existing, inplace, objective)

        # Set the merged model object and its ID,
        # then add the module attribute of the right model into the left
        new_model = super(MassModel, self).merge(
            right, prefix_existing, inplace, objective
        )

        # Add boundary conditions from right to left model.
        existing = [bc for bc in iterkeys(new_model.boundary_conditions)]
        new_model.add_boundary_conditions(
            {
                m: bc
                for m, bc in iteritems(right.boundary_conditions)
                if bc not in existing
            }
        )

        # Add custom parameters from right to left model.
        existing = [cp for cp in iterkeys(new_model.custom_parameters)]
        new_model.custom_parameters.update(
            {
                cp: v
                for cp, v in iteritems(right.custom_parameters)
                if cp not in existing
            }
        )

        def prefix_existing_id(to_prefix):
            """Prefix the ID and return it."""
            return "{0}{1}".format(prefix_existing, to_prefix)

        # Add custom rates from right to left model,
        # prefixing any existing reactions if necessary
        if prefix_existing is not None:
            existing = dict(
                (prefix_existing_id(r.id), rate)
                if prefix_existing_id(r.id) in new_model.reactions
                else (r.id, rate)
                for r, rate in iteritems(right.custom_rates)
            )
        else:
            existing = {}

        existing.update(right.custom_rates)
        new_model.custom_rates.update(
            {
                new_model.reactions.get_by_id(getattr(r, "_id", r)): rate
                for r, rate in iteritems(existing)
            }
        )

        new_items = deepcopy(right.groups)
        if prefix_existing is not None:
            existing = new_items.query(lambda group: group.id in self.groups)
            for group in existing:
                group.id = prefix_existing_id(group.id)
        new_model.add_groups(new_items)

        def existing_enzyme_filter(enzyme):
            """Filter existing EnzymeModules."""
            if enzyme.id in self.enzyme_modules:
                LOGGER.warning(
                    "Ignoring enzyme module '%s' since it already exists.", enzyme.id
                )
                return False
            return True

        # Add enzyme_modules from right to left model
        if right.enzyme_modules:
            new_items = deepcopy(right.enzyme_modules)
            # Prefix enzyme_modules if necessary
            if prefix_existing is not None:
                existing = new_items.query(lambda e: e.id in self.enzyme_modules)
                for enzyme in existing:
                    enzyme.__dict__["_id"] = prefix_existing_id(enzyme.id)

            # Check whether reactions exist in the model.
            new_items = DictList(filter(existing_enzyme_filter, new_items))
            new_model.enzyme_modules += new_items
            for enzyme in new_model.enzyme_modules:
                enzyme.model = new_model

        new_model.add_units([u for u in right.units if u not in new_model.units])

        for attr in ["_compartments", "notes", "annotation"]:
            existing = getattr(new_model, attr).copy()
            setattr(new_model, attr, getattr(right, attr).copy())
            getattr(new_model, attr).update(existing)

        return new_model

    def compute_steady_state_fluxes(
        self, pathways, independent_fluxes, update_reactions=False
    ):
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
            A ``dict`` of steady state fluxes where :class:`~.MassReaction`\ s
            are keys and fluxes are values to utilize in order to calculate
            all other steady state fluxes. Must be the same length as the
            number of specified pathways.
        update_reactions : bool
            If ``True`` then update the :attr:`.MassReaction.steady_state_flux`
            with the calculated steady state flux value for each reaction.

        Return
        ------
        dict
            A ``dict`` where key:value pairs are the :class:`~.MassReaction`\ s
            with their corresponding calculated steady state fluxes.

        Warnings
        --------
        The indicies of the values in the pathway vector must correspond to the
        indicies of the reactions in the :attr:`reactions` attribute in order
        for the method to work as intended.

        """
        # Check inputs:
        if not isinstance(pathways, (np.ndarray, list)):
            raise TypeError(
                "Pathways must be numpy.ndarrays or array-like, "
                "such as a list of lists."
            )
        pathways = np.array(pathways)
        if len(self.reactions) != pathways.shape[1]:
            raise ValueError(
                "Pathways must have the same number of columns as"
                " the number of reactions in the model."
            )
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

    def calculate_PERCs(
        self,
        at_equilibrium_default=100000,
        update_reactions=False,
        verbose=False,
        **kwargs
    ):
        r"""Calculate pseudo-order rate constants for reactions in the model.

        Pseudo-order rate constants (PERCs) are considered to be the same as
        :attr:`~.MassReaction.forward_rate_constant` attributes, and are
        calculated based on the steady state concentrations and fluxes.

        Notes
        -----
        * All fluxes and concentrations used in calculations must be provided,
          including relevant boundary conditions. By default, the relevant
          values are taken from objects associated with the model.
        * To calculate PERCs for a subset of model reactions, use the
          ``fluxes`` kwawrg.

        Parameters
        ----------
        at_equilibrium_default : float
            The value to set the pseudo-order rate constant if the reaction is
            at equilibrium. Default is ``100,000``.
        update_reactions : bool
            If ``True`` then will update the values for the
            :attr:`~MassReaction.forward_rate_constant` attributes with the
            calculated PERC values.
        verbose : bool
            Whether to output more verbose messages for errors and logging.
        **kwargs
            fluxes :
                A ``dict`` of reaction fluxes where :class:`~MassReaction`\ s
                are keys and fluxes are the values. Only reactions provided
                will have their PERCs calculated. If ``None``, PERCs are
                calculated using the current steady state fluxes for all
                reactions in the model.

                Default is ``None``.
            concentrations : dict
                A ``dict`` of concentrations necessary for the PERC
                calculations, where :class:`~.MassMetabolite`\ s are keys and
                concentrations are the values. If ``None``, the relevant
                concentrations that exist in the model are used.

                Default is ``None``.

        Returns
        -------
        dict
            A ``dict`` where keys are strings identifers of the pseudo-order
            rate constants (as given by :attr:`.MassReaction.kf_str`) and
            values are the calculated PERC values.

        """
        kwargs = _check_kwargs({"fluxes": None, "concentrations": None}, kwargs)
        # Get the model steady state concentrations if None are provided.
        if kwargs.get("concentrations") is None:
            concentrations = self.boundary_conditions.copy()
            concentrations.update(
                {str(m): ic for m, ic in iteritems(self.initial_conditions)}
            )
        else:
            concentrations = {
                str(m): v for m, v in iteritems(kwargs.get("concentrations"))
            }

        # Get the model reactions and fluxes if None are provided.
        if kwargs.get("fluxes") is None:
            fluxes = self.steady_state_fluxes
        else:
            fluxes = kwargs.get("fluxes")
            invalid = {
                rxn: flux
                for rxn, flux in iteritems(fluxes)
                if not isinstance(rxn, MassReaction)
                or not isinstance(flux, (float, integer_types))
            }
            if invalid:
                raise TypeError("Invalid ``fluxes``: '{0!r}'".format(invalid))

        # Get defined numerical values
        numerical_values = {
            str(param): value
            for p_type, subdict in iteritems(self.parameters)
            for param, value in iteritems(subdict)
            if p_type not in ["kf", "v"]
        }
        numerical_values.update(concentrations)

        # Function to calculate the solution
        def calculate_sol(flux, rate_equation, perc):
            sol = sym.solveset(sym.Eq(flux, rate_equation), perc, domain=sym.S.Reals)
            if (
                isinstance(sol, type(sym.S.Reals))
                or sol.is_empty
                or next(iter(sol)) == 0
            ):
                sol = float(at_equilibrium_default)
            else:
                sol = float(next(iter(sol)))

            return sol

        # Calculate the PERCs
        percs_dict = {}
        for reaction, flux in iteritems(fluxes):
            rate_eq = strip_time(reaction.rate)
            if reaction.kf_str not in str(rate_eq):
                msg = "Skipping reaction {0}, no PERC in rate".format(reaction.id)
                LOGGER.info(msg)
                if verbose:
                    print(msg)
                continue

            # Get arguments
            args = [
                a for a in list(rate_eq.atoms(sym.Symbol)) if str(a) != reaction.kf_str
            ]

            # Get numerical values for arguments
            vals = {
                a: numerical_values[str(a)] for a in args if str(a) in numerical_values
            }

            if len(args) != len(vals):
                warnings.warn(
                    "Cannot calculate the PERC for reaction '{0}' missing "
                    "values for {1!r}.".format(
                        reaction.id, [str(a) for a in args if a not in vals]
                    )
                )
                continue

            # Substitute values into rate equation
            rate_eq = rate_eq.subs(vals)
            # Calculate rate equation and update with soluton for PERC
            sol = calculate_sol(flux, rate_eq, sym.Symbol(reaction.kf_str))
            percs_dict.update({reaction.kf_str: sol})
            if update_reactions:
                reaction.kf = percs_dict[reaction.kf_str]

        return percs_dict

    def build_model_from_string(
        self,
        model_str,
        verbose=True,
        reaction_split=";",
        reaction_id_split=":",
        **kwargs
    ):
        """Create a :class:`MassModel` from strings of reaction equations.

        Accepted ``kwargs`` are passed to the underlying function for reaction
        creation, :meth:`.MassReaction.build_reaction_from_string`.

        Takes a string representation of the reactions and uses the
        specifications supplied in the optional arguments to first infer a set
        of reactions and their identifiers, then to infer metabolites,
        metabolite compartments, and stoichiometries for the reactions. It
        also infers the reversibility of the reaction from the reaction arrow.
        For example::

            '''
            RID_1: S + E <=> ES;
            RID_2: ES -> E + P;
            RID_3: E + I <=> EI;
            '''

        where ``RID`` represents the identifier to assign the
        :class:`~.MassReaction`.

        Parameters
        ----------
        model : str
            A string representing the reaction formulas (equation) for the
            model.
        verbose : bool
            Setting the verbosity of the function.
        reaction_split : str
            Dividing individual reaction entries. Default is ``";"``.
        reaction_id_split : str
            Dividing individual reaction entries from their identifiers.
            Default is ``":"``.
        **kwargs
            fwd_arrow :
                :func:`re.compile` or ``None`` for forward irreversible
                reaction arrows. If ``None``, the arrow is expected to
                be ``'-->'`` or ``'==>'``.
            rev_arrow :
                :func:`re.compile` or ``None`` for backward irreversible
                reaction arrows. If ``None``, the arrow is expected to
                be ``'<--'`` or ``'<=='``.
            reversible_arrow :
                :func:`re.compile` or ``None`` for reversible reaction arrows.
                If ``None``, the arrow is expected to
                be ``'<=>'`` or ``'<->'``.
            term_split : str
                Dividing individual metabolite entries. Default is ``"+"``.

        See Also
        --------
        :meth:`.MassReaction.build_reaction_from_string`
            Base method for building reactions.

        """
        # Use the reaction split arguments to get the reactions and strip them
        reaction_list = [
            reaction_str.strip()
            for reaction_str in model_str.split(reaction_split)
            if reaction_str.strip()
        ]

        # Iterate through reaction strings
        for orig_reaction_str in reaction_list:
            # Split the reaction ID from the reaction equation
            split = orig_reaction_str.split(reaction_id_split)
            try:
                if len(split) != 2:
                    raise ValueError(
                        "Could not parse '{0}' for the reaction "
                        "ID and formula (equation).".format(orig_reaction_str)
                    )
                reaction_id, reaction_str = (s.strip() for s in split)
                # Cannot build reaction without an ID
                if not reaction_id:
                    raise ValueError(
                        "No reaction ID found in '{0}'".format(orig_reaction_str)
                    )
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
                    reaction_str, verbose=verbose, **kwargs
                )
            except ValueError as e:
                # Raise warnings for reactions that could not be built.
                warnings.warn(
                    "Failed to build reaction '{0}' due to the following:\n"
                    "{1}".format(orig_reaction_str, str(e))
                )
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
            A ``dict`` containing the parameter identifiers as strings and
            their corresponding values to set in the model.
        verbose : bool
            If ``True`` then display the warnings that may be raised when
            setting reaction parameters. Default is ``True``.

        See Also
        --------
        :attr:`.MassReaction.all_parameter_ids`
            Lists the default reaction parameter identifiers.
        :attr:`MassModel.boundary_metabolites`
            Lists the 'boundary metabolites' found in the model.

        """
        if not isinstance(parameters, dict):
            raise TypeError("parameters must be a dict.")

        for key, value in iteritems(parameters):
            if not isinstance(key, string_types):
                raise TypeError("Keys must be strings. '{0}' not a string.".format(key))

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
                            warnings.filterwarnings(
                                "ignore", ".*constant for an irreversible reaction.*"
                            )
                        setattr(reaction, p_type, value)
                except (KeyError, ValueError):
                    self.custom_parameters.update({key: value})
            # If parameter not found, assume parameter is a custom parameter
            else:
                self.custom_parameters.update({key: value})

    def update_initial_conditions(self, initial_conditions, verbose=True):
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
            A ``dict`` where metabolites are the keys and the initial
            conditions are the values.
        verbose : bool
            If ``True`` then display the warnings that may be raised when
            setting metabolite initial conditions. Default is ``True``.

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
                if verbose:
                    warnings.warn(
                        "Cannot set initial condition for {0} due to the "
                        "following: {1}".format(metabolite.id, str(e))
                    )
                continue

    def update_custom_rates(self, custom_rates, custom_parameters=None):
        r"""Update the custom rates of the model.

        Parameters
        ----------
        custom_rates : dict
            A ``dict`` where :class:`.MassReaction`\ s or their string
            identifiers are the keys and the rates are the string
            representations of the custom rate expression.
        custom_parameters : dict
            A ``dict`` of custom parameters for the custom rates, where the
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

        # Iterate through custrom rates
        for reaction, custom_rate in iteritems(custom_rates):
            if not isinstance(reaction, MassReaction):
                # Ensure reaction exists
                try:
                    reaction = self.reactions.get_by_id(reaction)
                except KeyError as e:
                    warnings.warn("No reaction found for {0}".format(str(e)))
                    continue
            # Try to create the rate expression
            try:
                self.add_custom_rate(reaction, custom_rate=custom_rate)
            except sym.SympifyError:
                # Warn if fail
                warnings.warn(
                    "Unable to sympify rate equation for " "'{0}'.".format(reaction.id)
                )

    def has_equivalent_odes(self, right, verbose=False):
        """Determine whether :attr:`odes` between two models are equivalent.

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
        l_odes_and_rates = [
            {met.id: ode for met, ode in iteritems(self.odes)},
            {r.id: rate for r, rate in iteritems(self.rates)},
        ]
        r_odes_and_rates = [
            {met.id: ode for met, ode in iteritems(right.odes)},
            {r.id: rate for r, rate in iteritems(right.rates)},
        ]
        if l_odes_and_rates[0] != r_odes_and_rates[0]:
            equivalent = False

        if not equivalent and verbose:
            msg_out = "{0} vs. {1}:\n".format(self.id, right.id)
            for i, (l_dict, r_dict) in enumerate(
                zip(l_odes_and_rates, r_odes_and_rates)
            ):
                # Determine which metabolites do not exist in both models
                missing = set(l_dict)
                missing.symmetric_difference_update(set(r_dict))
                # Determine which metabolites have different ODEs
                diff_equations = set(
                    l_key
                    for l_key, l_value in iteritems(l_dict)
                    if l_key in r_dict and l_value != r_dict[l_key]
                )

                # Format and return messages for differences
                if i == 0:
                    msgs = ["Metabolites", "ODEs"]
                else:
                    msgs = ["Reactions", "rates"]

                msgs = [
                    "{0} in one model only: ".format(msgs[0]),
                    "{0} with different {1}: ".format(*msgs),
                ]

                for item, msg in zip([missing, diff_equations], msgs):
                    if item:
                        msg_out += msg + str(sorted(list(item))) + "\n"

            print(msg_out.rstrip("\n"))

        return equivalent

    def set_steady_state_fluxes_from_solver(self):
        """Set reaction steady state fluxes based on the state of the solver.

        Only works when reaction is associated with a model that has been
        optimized.

        """
        for reaction in self.reactions:
            reaction.steady_state_flux = reaction.flux

    # Internal
    def _cobra_to_mass_repair(self):
        """Convert associated cobra objects to mass objects for self.

        Warnings
        --------
        This method is intended for internal use only.

        """
        # Convert Metabolites into MassMetabolites
        if self.metabolites:
            self.metabolites = DictList(
                [MassMetabolite(metabolite) for metabolite in self.metabolites]
            )
        # Copy genes
        if self.genes:
            self.genes = DictList([gene.copy() for gene in self.genes])
        # Convert Reactions into MassReactions
        if self.reactions:
            self.reactions = DictList(
                [MassReaction(reaction) for reaction in self.reactions]
            )
        # Convert groups containing metabolites and reactions
        # into groups containing MassMetabolites and MassReactions
        if self.groups:
            for group in self.groups:
                old_members = list(group.members)
                new_members = []
                for member in old_members:
                    if isinstance(member, Metabolite) and not isinstance(
                        member, MassMetabolite
                    ):
                        new_members.append(MassMetabolite(member))
                    if isinstance(member, Reaction) and not isinstance(
                        member, MassReaction
                    ):
                        new_members.append(MassReaction(member))
                    if not isinstance(member, (Metabolite, Reaction)):
                        new_members.append(member)
                group.remove_members(old_members)
                group.add_members(new_members)
        # Ensure all objects in the model point to the MassModel
        self.__setstate__(self.__dict__)

    def _mk_stoich_matrix(self, array_type=None, dtype=None, update_model=True):
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

        # Use current array type and dtype if None provided
        if array_type is None:
            array_type = self._array_type

        if dtype is None:
            dtype = self._dtype

        stoich_mat = create_stoichiometric_matrix(self)

        # Convert the matrix to the desired type
        stoich_mat = convert_matrix(
            stoich_mat,
            array_type=array_type,
            dtype=dtype,
            row_ids=[m.id for m in self.metabolites],
            col_ids=[r.id for r in self.reactions],
        )
        # Update the stored stoichiometric matrix for the model if True
        if update_model:
            self._S = stoich_mat
            self._array_type = array_type
            self._dtype = dtype

        return stoich_mat

    def _get_all_parameters(self):
        """Get a dict containing all of defined model parameters in the model.

        Warnings
        --------
        This method is intended for internal use only.

        """
        return {
            str(param): value
            for subdict in itervalues(self.parameters)
            for param, value in iteritems(subdict)
        }

    def _copy_model_metabolites(self, new_model):
        """Copy the metabolites in creating a partial "deepcopy" of model.

        Warnings
        --------
        This method is intended for internal use only.

        """
        # Initialize DictList and set attributes to not copy by ref.
        new_metabolites = DictList()
        do_not_copy_by_ref = {"_reaction", "_model"}
        # Copy the metabolites
        for metabolite in self.metabolites:
            new_metabolite = metabolite.__class__()
            for attr, value in iteritems(metabolite.__dict__):
                if attr not in do_not_copy_by_ref:
                    new_metabolite.__dict__[attr] = (
                        copy(value) if attr == "formula" else value
                    )
            new_metabolite._model = new_model
            new_metabolites.append(new_metabolite)

        return new_metabolites

    def _copy_model_genes(self, new_model):
        """Copy the genes in creating a partial "deepcopy" of model.

        Warnings
        --------
        This method is intended for internal use only.

        """
        # Initialize DictList and set attributes to not copy by ref.
        new_genes = DictList()
        do_not_copy_by_ref = {"_reaction", "_model"}
        # Copy the genes
        for gene in self.genes:
            new_gene = gene.__class__(None)
            for attr, value in iteritems(gene.__dict__):
                if attr not in do_not_copy_by_ref:
                    new_gene.__dict__[attr] = (
                        copy(value) if attr == "formula" else value
                    )
            new_gene._model = new_model
            new_genes.append(new_gene)

        return new_genes

    def _copy_model_reactions(self, new_model):
        """Copy the reactions in creating a partial "deepcopy" of model.

        Warnings
        --------
        This method is intended for internal use only.

        """
        # Initialize DictList and set attributes to not copy by ref.
        new_reactions = DictList()
        do_not_copy_by_ref = {"_model", "_metabolites", "_genes"}
        # Copy the reactions
        for reaction in self.reactions:
            new_reaction = reaction.__class__()
            for attr, value in iteritems(reaction.__dict__):
                if attr not in do_not_copy_by_ref:
                    new_reaction.__dict__[attr] = copy(value)
            new_reaction._model = new_model
            new_reactions.append(new_reaction)
            # Update awareness for metabolites and genes
            for metabolite, stoich in iteritems(reaction._metabolites):
                new_metabolite = new_model.metabolites.get_by_id(metabolite.id)
                new_reaction._metabolites[new_metabolite] = stoich
                new_metabolite._reaction.add(new_reaction)
            for gene in reaction._genes:
                new_gene = new_model.genes.get_by_id(gene.id)
                new_reaction._genes.add(new_gene)
                new_gene._reaction.add(new_reaction)

        return new_reactions

    def _copy_model_enzyme_modules(self, new_model):
        """Copy the enzyme_modules in creating a partial "deepcopy" of model.

        Warnings
        --------
        This method is intended for internal use only.

        """
        new_enzyme_modules = DictList()
        do_not_copy_by_ref = {"model"}
        # Copy the enzyme_modules
        for enzyme in self.enzyme_modules:
            new_enzyme = enzyme.__class__()
            for attr, value in iteritems(enzyme):
                if attr not in do_not_copy_by_ref:
                    new_enzyme[attr] = copy(value)
            # Update associated model and object pointers for enzyme
            new_enzyme["model"] = new_model
            new_enzyme._update_object_pointers(new_model)
            new_enzyme_modules.append(new_enzyme)

        return new_enzyme_modules

    def _copy_model_groups(self, new_model):
        """Copy the groups in creating a partial "deepcopy" of model.

        Warnings
        --------
        This method is intended for internal use only.

        """
        new_groups = DictList()
        do_not_copy_by_ref = {"_model", "_members"}
        # Groups can be members of other groups. We initialize them first and
        # then update their members.
        for group in self.groups:
            new_group = group.__class__(group.id)
            for attr, value in iteritems(group.__dict__):
                if attr not in do_not_copy_by_ref:
                    new_group.__dict__[attr] = copy(value)
            new_group._model = new_model
            new_groups.append(new_group)
        # Iterate through groups
        for group in self.groups:
            new_group = new_groups.get_by_id(group.id)
            # Update awareness
            new_objects = []
            for member in group.members:
                if isinstance(member, MassMetabolite):
                    new_object = new_model.metabolites.get_by_id(member.id)
                elif isinstance(member, MassReaction):
                    new_object = new_model.reactions.get_by_id(member.id)
                elif isinstance(member, Gene):
                    new_object = new_model.genes.get_by_id(member.id)
                elif isinstance(member, Group):
                    new_object = new_model.genes.get_by_id(member.id)
                elif member.__class__.__name__ == "EnzymeModuleDict":
                    new_object = new_model.enzyme_modules.get_by_id(member.id)
                else:
                    raise TypeError(
                        "The group member {!r} is unexpectedly not a "
                        "metabolite, reaction, gene, enzyme module, nor "
                        "another group.".format(member)
                    )
                new_objects.append(new_object)
            new_group.add_members(new_objects)

        return new_groups

    def _repr_html_(self):
        """HTML representation of the overview for the MassModel.

        Warnings
        --------
        This method is intended for internal use only.

        """
        try:
            dim_S = "{0}x{1}".format(self.S.shape[0], self.S.shape[1])
            rank = matrix_rank(self.S)
        except (np.linalg.LinAlgError, ValueError, IndexError):
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
                    <td><strong>Number of metabolites</strong></td>
                    <td>{num_metabolites}</td>
                </tr><tr>
                    <td><strong>Initial conditions defined</strong></td>
                    <td>{num_ic}/{num_metabolites}</td>
                </tr><tr>
                    <td><strong>Number of reactions</strong></td>
                    <td>{num_reactions}</td>
                </tr><tr>
                    <td><strong>Number of genes</strong></td>
                    <td>{num_genes}</td>
                </tr><tr>
                    <td><strong>Number of enzyme modules</strong></td>
                    <td>{num_enzyme_modules}</td>
                </tr><tr>
                    <td><strong>Number of groups</strong></td>
                    <td>{num_groups}</td>
                </tr><tr>
                    <td><strong>Objective expression</strong></td>
                    <td>{objective}</td>
                </tr><tr>
                    <td><strong>Compartments</strong></td>
                    <td>{compartments}</td>
                </tr>
            </table>
        """.format(
            name=self.id,
            address="0x0%x" % id(self),
            dim_stoich_mat=dim_S,
            mat_rank=rank,
            num_metabolites=len(self.metabolites),
            num_ic=len(self.initial_conditions),
            num_reactions=len(self.reactions),
            num_genes=len(self.genes),
            num_enzyme_modules=len(self.enzyme_modules),
            num_groups=len(self.groups),
            objective=format_long_string(str(self.objective.expression), 100),
            compartments=", ".join(
                v if v else k for k, v in iteritems(self.compartments)
            ),
        )

    def __setstate__(self, state):
        """Ensure all objects in the model point to the MassModel.

        Extends
        :meth:`Model.__setstate__ <cobra.core.model.Model.__setstate__>`
        to include enzyme_modules

        Warnings
        --------
        This method is intended for internal use only.

        """
        super(MassModel, self).__setstate__(state)
        for attr in ["enzyme_modules"]:
            # Ensure attribute exists, needed for preventing inheritance issues
            if hasattr(self, attr):
                for x in getattr(self, attr):
                    x._model = self

    def __getstate__(self):
        """Get the state for serialization.

        Ensures that the context stack is cleared prior to serialization,
        since partial functions cannot be pickled reliably

        Warnings
        --------
        This method is intended for internal use only.

        """
        return super(MassModel, self).__getstate__()

    def __dir__(self):
        """Override default dir() implementation to list only public items.

        Warnings
        --------
        This method is intended for internal use only.

        """
        return get_public_attributes_and_methods(self)


__all__ = ("MassModel", "LOGGER", "CHOPNSQ")
