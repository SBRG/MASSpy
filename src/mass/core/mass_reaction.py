# -*- coding: utf-8 -*-
"""MassReaction is a class for holding information regarding reactions.

The :class:`MassReaction` class inherits and extends the
:class:`~cobra.core.reaction.Reaction` class in :mod:`cobra`. It contains
additional information required for simulations and other :mod:`mass`
functions and workflows.

Some key differences between the
:class:`cobra.Reaction <cobra.core.reaction.Reaction>` and the
:class:`mass.MassReaction <mass.core.mass_reaction.MassReaction>` are
listed below:

    * When instantiating a :class:`MassReaction` from a
      :class:`cobra.Reaction <cobra.core.reaction.Reaction>`, any associated
      :class:`cobra.Metabolite <cobra.core.metabolite.Metabolite>` will also
      be converted into a :class:`.MassMetabolite`.

    * Unlike the :class:`cobra.Reaction <cobra.core.reaction.Reaction>`
      which initializes the :attr:`~cobra.core.reaction.Reaction.lower_bound`
      at a default value of ``0.0``, the :class:`MassReaction` initializes the
      :attr:`~cobra.core.reaction.Reaction.lower_bound` as ``None`` to utilize
      the default value set in the :class:`.MassConfiguration`.

    * The :class:`MassReaction` contains both the
      :attr:`~MassReaction.steady_state_flux` and the inherited
      :attr:`cobra.Reaction.flux <cobra.core.reaction.Reaction.flux>`
      attributes. Note that these attributes **DO NOT** refer to the same flux
      value unless specifically set as such. The
      :attr:`~cobra.core.reaction.Reaction.flux` refers to the flux value in
      the most recent optimization solution, while the
      :attr:`~MassReaction.steady_state_flux` refers to the net steady
      state flux through the reaction.

    * The :class:`MassReaction` contains both the
      :attr:`~MassReaction.reversible`
      and :attr:`~cobra.core.reaction.Reaction.reversibility` attributes.

      Note that these attributes are **NOT** the same. While the
      :attr:`~MassReaction.reversible` refers the the kinetic reversibility of
      the reaction and thus determines how the :attr:`~MassReaction.rate` is
      set, while the inhertied
      :attr:`cobra.Reaction.reversibility <cobra.core.reaction.Reaction.reversibility>`
      refers to the computed reversibility based on the current lower
      and upper bounds.

    * The arrow in the string representation of the reaction is dependent on
      the :attr:`MassReaction.reversible` rather than the inhertied
      :attr:`cobra.Reaction.reversibility <cobra.core.reaction.Reaction.reversibility>`
      attribute.

"""  # noqa: E501
import re
import warnings
from copy import copy, deepcopy
from functools import partial
from operator import attrgetter

from cobra.core.metabolite import Metabolite
from cobra.core.reaction import (
    Reaction,
    _reverse_arrow_finder,
    _reversible_arrow_finder,
)
from cobra.util.context import get_context, resettable
from cobra.util.util import format_long_string
from six import iteritems, itervalues, string_types

from mass.core.mass_configuration import MassConfiguration
from mass.core.mass_metabolite import MassMetabolite
from mass.util.expressions import (
    generate_disequilibrium_ratio,
    generate_forward_mass_action_rate_expression,
    generate_mass_action_rate_expression,
    generate_mass_action_ratio,
    generate_reverse_mass_action_rate_expression,
)
from mass.util.util import (
    _check_kwargs,
    ensure_non_negative_value,
    get_public_attributes_and_methods,
)


MASSCONFIGURATION = MassConfiguration()


class MassReaction(Reaction):
    """Class for holding kinetic information regarding a biochemical reaction.

    Accepted ``kwargs`` are passed onto the initialization method for the
    base class :class:`~cobra.core.reaction.Reaction`.

    Parameters
    ----------
    id_or_reaction : str, ~cobra.core.reaction.Reaction, MassReaction
        A string identifier to associate with the reaction, or an existing
        reaction. If an existing reaction object is provided, a new
        :class:`MassReaction` object is instantiated with the same properties
        as the original reaction.
    name : str
        A human readable name for the reaction.
    subsystem : str
        The subsystem where the reaction is meant to occur.
    reversible : bool
        The kinetic reversibility of the reaction. Irreversible reactions have
        an equilibrium constant and a reverse rate constant as set in the
        :attr:`~.MassBaseConfiguration.irreversible_Keq` and
        :attr:`~.MassBaseConfiguration.irreversible_kr` attributes of the
        :class:`~.MassConfiguration`. Default is ``True``.
    **kwargs
        lower_bound : float or None
            The lower flux bound for optimization. If ``None`` then
            the default bound from the :class:`MassConfiguration` is used.

            Default is ``None``.
        upper_bound : float or None
            The upper flux bound for optimization. If ``None`` then
            the default bound from the :class:`MassConfiguration` is used.

            Default is ``None``.

    """

    def __init__(
        self,
        id_or_reaction=None,
        name="",
        subsystem="",
        reversible=True,
        steady_state_flux=None,
        **kwargs
    ):
        """Initialize the MassReaction."""
        # Check kwargs
        kwargs = _check_kwargs(
            {
                "lower_bound": MASSCONFIGURATION.lower_bound,
                "upper_bound": MASSCONFIGURATION.upper_bound,
            },
            kwargs,
        )
        # Get the identifer and initialize
        super(MassReaction, self).__init__(
            getattr(id_or_reaction, "id", id_or_reaction),
            name=name,
            subsystem=subsystem,
            **kwargs
        )
        if isinstance(id_or_reaction, (Reaction, MassReaction)):
            # Instiantiate a new MassReaction with state identical to
            # the provided MassReaction object.
            self.__dict__.update(id_or_reaction.__dict__)

            # Change associated cobra.Metabolites to MassMetabolites
            if isinstance(id_or_reaction, Reaction) and not isinstance(
                id_or_reaction, MassReaction
            ):
                self._cobra_to_mass_repair()

        # If is not a MassReaction object, initialize additional attributes
        if not isinstance(id_or_reaction, MassReaction):
            self._reversible = reversible
            self._steady_state_flux = steady_state_flux
            # Rate and equilibrium constant parameters for reactions. The
            # reverse and equilibrium constants for irreversible reactions are
            # set here.  Upper and lower bounds are also set.
            self._forward_rate_constant = None
            if self._reversible:
                self._reverse_rate_constant = None
                self._equilibrium_constant = None
            else:
                self._reverse_rate_constant = MASSCONFIGURATION.irreversible_kr
                self._equilibrium_constant = MASSCONFIGURATION.irreversible_Keq

            # Rate type as a sympy expression for simulation.
            self._rate_type = 1

    # Public
    @property
    def reversible(self):
        """Get or set the kinetic reversibility of the reaction.

        When using a `HistoryManager` context, this attribute can be set
        temporarily, reversed when the exiting the context.

        Parameters
        ----------
        reversible : bool
            The kinetic reversibility of the reaction.

        Warnings
        --------
        Changing the :attr:`reversible` attribute will reset the
        :attr:`equilibrium_constant` and the :attr:`reverse_rate_constant`
        to the initialization defaults.

        """
        return getattr(self, "_reversible")

    @reversible.setter
    def reversible(self, reversible):
        """Set the kinetic reversibility of the reaction."""
        if not isinstance(reversible, bool):
            raise TypeError("value must be a bool")

        # No need to do anything if reversible hasn't changed
        if reversible != self.reversible:
            setattr(self, "_reversible", reversible)

            context = get_context(self)
            if context:
                existing = [self._reverse_rate_constant, self._equilibrium_constant]

            if reversible:
                setattr(self, "_reverse_rate_constant", None)
                setattr(self, "_equilibrium_constant", None)
            else:
                setattr(
                    self, "_reverse_rate_constant", MASSCONFIGURATION.irreversible_kr
                )
                setattr(
                    self, "_equilibrium_constant", MASSCONFIGURATION.irreversible_Keq
                )

            if context:
                context(partial(setattr, self, "_reverse_rate_constant", existing[0]))
                context(partial(setattr, self, "_equilibrium_constant", existing[1]))

    @property
    def steady_state_flux(self):
        """Get or set the steady state flux of the reaction.

        When using a `HistoryManager` context, this attribute can be set
        temporarily, reversed when the exiting the context.

        Parameters
        ----------
        flux_value : bool
            The steady state flux value of the reaction.

        """
        return getattr(self, "_steady_state_flux")

    @steady_state_flux.setter
    @resettable
    def steady_state_flux(self, flux_value):
        """Set the steady state flux of the reaction."""
        setattr(self, "_steady_state_flux", flux_value)

    @property
    def forward_rate_constant(self):
        """Get or set the forward rate constant (kf) of the reaction.

        When using a `HistoryManager` context, this attribute can be set
        temporarily, reversed when the exiting the context.

        Notes
        -----
        Forward rate constants cannot be negative.

        Parameters
        ----------
        value : float
            A non-negative number for the forward rate constant (kf) of the
            reaction.

        Raises
        ------
        ValueError
            Occurs when trying to set a negative value.

        """
        return getattr(self, "_forward_rate_constant")

    @forward_rate_constant.setter
    @resettable
    def forward_rate_constant(self, value):
        """Set the forward rate constant (kf) of the reaction."""
        value = ensure_non_negative_value(value)
        setattr(self, "_forward_rate_constant", value)

    @property
    def reverse_rate_constant(self):
        """Get or set the reverse rate constant (kr) of the reaction.

        When using a `HistoryManager` context, this attribute can be set
        temporarily, reversed when the exiting the context.

        Notes
        -----
        * Reverse rate constants cannot be negative.

        * If reaction is not reversible, will warn the user instead of setting
          the parameter value.

        Parameters
        ----------
        value : float
            A non-negative number for the reverse rate constant (kr) of the
            reaction.

        Raises
        ------
        ValueError
            Occurs when trying to set a negative value.

        """
        return getattr(self, "_reverse_rate_constant")

    @reverse_rate_constant.setter
    @resettable
    def reverse_rate_constant(self, value):
        """Set the reverse rate constant (kr) of the reaction."""
        if self.reversible or value in [MASSCONFIGURATION.irreversible_kr, None]:
            value = ensure_non_negative_value(value)
            setattr(self, "_reverse_rate_constant", value)
        else:
            warnings.warn(
                "Cannot set the reverse rate constant for an irreversible "
                "reaction '{0}'".format(self.id)
            )

    @property
    def equilibrium_constant(self):
        """Get or set the equilibrium constant (Keq) of the reaction.

        When using a `HistoryManager` context, this attribute can be set
        temporarily, reversed when the exiting the context.

        Notes
        -----
        * Equilibrium constants cannot be negative.

        * If reaction is not reversible, will warn the user instead of setting
          the parameter value.

        Parameters
        ----------
        value : float
            A non-negative number for the equilibrium constant (Keq)
            of the reaction.

        Raises
        ------
        ValueError
            Occurs when trying to set a negative value.

        """
        return getattr(self, "_equilibrium_constant")

    @equilibrium_constant.setter
    @resettable
    def equilibrium_constant(self, value):
        """Set the equilibrium constant (Keq) of the reaction."""
        if self.reversible or value in [MASSCONFIGURATION.irreversible_Keq, None]:
            value = ensure_non_negative_value(value)
            setattr(self, "_equilibrium_constant", value)
        else:
            warnings.warn(
                "Cannot set the equilibrium constant for an irreversible "
                "reaction '{0}'".format(self.id)
            )

    @property
    def parameters(self):
        """Return a ``dict`` of rate and equilibrium constants.

        Notes
        -----
        The :attr:`reverse_rate_constant` is only included for reversible
        reactions. Additionally, only rate and equilibrium constants are
        accessed here. Steady state fluxes can be accessed through the
        :attr:`steady_state_flux` attribute, and custom parameters can only be
        accessed through the model.

        """
        keys = [self.kf_str, self.Keq_str, self.kr_str]
        attrs = [
            "_forward_rate_constant",
            "_equilibrium_constant",
            "_reverse_rate_constant",
        ]
        # Return reverse rate constants for reversible reactions.
        parameters = {
            key: getattr(self, attr)
            for key, attr in zip(keys, attrs)
            if getattr(self, attr) is not None
        }

        return parameters

    @property
    def metabolites(self):
        """Return the metabolites of a reaction as a read only copy."""
        return getattr(self, "_metabolites").copy()

    @property
    def reactants(self):
        """Return a ``list`` of reactants for the reaction."""
        return super(MassReaction, self).reactants

    @property
    def products(self):
        """Return a ``list`` of products for the reaction."""
        return super(MassReaction, self).products

    @property
    def stoichiometry(self):
        """Return a ``list`` containing the stoichiometry for the reaction."""
        return [c for c in itervalues(self._metabolites)]

    @property
    def rate(self):
        """Return the current rate as a :mod:`sympy` expression.

        If reaction has a custom rate in its associated :class:`~.MassModel`,
        the custom rate will be returned instead.
        """
        if self.model is not None and self in self.model.custom_rates:
            return self.model.custom_rates[self]

        return self.get_mass_action_rate(rate_type=getattr(self, "_rate_type"))

    @property
    def model(self):
        """Return the  :class:`~.MassModel` associated with the reaction."""
        return getattr(self, "_model")

    @property
    def reaction(self):
        """Get or set the reaction as a human readable string.

        Parameters
        ----------
        reaction_str : str
            String representation of the reaction.

        Warnings
        --------
        Care must be taken when setting a reaction using this method.
        See documentation for :meth:`build_reaction_from_string` for more
        information.

        See Also
        --------
        build_reaction_string: Base function for getter method.

        """
        return self.build_reaction_string()

    @reaction.setter
    def reaction(self, reaction_str):
        """Set the reaction using a human readable string."""
        return self.build_reaction_from_string(reaction_str)

    @property
    def compartments(self):
        """Return the set of compartments where the metabolites are located."""
        return super(MassReaction, self).compartments

    @property
    def boundary(self):
        """Determine whether or not the reaction is a boundary reaction.

        Will return ``True`` if the reaction has no products or no reactants
        and only one metabolite.

        Notes
        -----
        These are reactions with a sink or a source term (e.g. 'A --> ')

        """
        return super(MassReaction, self).boundary

    @property
    def boundary_metabolite(self):
        """Return an 'boundary' metabolite for bounary reactions.

        Notes
        -----
        The 'boundary_metabolite' represents the metabolite that corresponds
        to the empty part of a boundary reaction through a string. It's primary
        use is for setting of the :attr:`~.MassModel.boundary_conditions`
        without creating a :class:`~.MassMetabolite` object. Therefore it is
        not counted as a metabolite but instead as a parameter.

        Returns
        -------
        boundary_metabolite : str
            String representation of the boundary metabolite of the reaction,
            or ``None`` if the reaction is not considered a boundary reaction.

        See Also
        --------
        boundary: Method must return ``True`` to get the boundary_metabolite.

        """
        if self.boundary:
            for metabolite in self.metabolites:
                bc_metabolite = metabolite._remove_compartment_from_id_str()
                bc_metabolite += "_" + str(
                    next(iter(MASSCONFIGURATION.boundary_compartment))
                )
        else:
            bc_metabolite = None

        return bc_metabolite

    @property
    def genes(self):
        """Return a ``frozenset`` of the genes associated with the reaction."""
        return frozenset(getattr(self, "_genes"))

    @property
    def gene_reaction_rule(self):
        """Get or set the gene reaction rule for the reaction.

        Parameters
        ----------
        new_rule : str
            String representation of the new reaction rule.

        Notes
        -----
        New genes will be associated with the reaction and old genes will be
        dissociated from the reaction.

        """
        return super(MassReaction, self).gene_reaction_rule

    @gene_reaction_rule.setter
    def gene_reaction_rule(self, new_rule):
        """Set the gene reaction rule of a reaction using a string."""
        super(MassReaction, self.__class__).gene_reaction_rule.fset(self, new_rule)

    @property
    def gene_name_reaction_rule(self):
        """Display gene_reaction_rule with names.

        Warnings
        --------
        Do NOT use this string for computation. It is intended to give a
        representation of the rule using more familiar gene names instead
        of the often cryptic ids.

        """
        return super(MassReaction, self).gene_name_reaction_rule

    @property
    def functional(self):
        """Check if all required enzymes for the reaction are functional.

        Returns
        -------
        bool
            Returns ``True`` if the gene-protein-reaction (GPR) rule is
            fulfilled for the reaction, or if the reaction does not have an
            assoicated :class:`~.MassModel`. Otherwise returns ``False``.

        """
        return super(MassReaction, self).functional

    @property
    def flux_symbol_str(self):
        """Return the string representation for the reaction flux symbol."""
        if self.id is None:
            return None

        return str("v_" + self.id)

    @property
    def all_parameter_ids(self):
        """Return ``list`` of strings representing non-custom parameters."""
        return [self.kf_str, self.Keq_str, self.kr_str, self.flux_symbol_str]

    @property
    def kf_str(self):
        """Return the string representation of the forward rate constant."""
        if self.id is None:
            warnings.warn("No parameter ID. Define the reaction ID first.")
            return None

        return "kf_" + self.id

    @property
    def Keq_str(self):
        """Return the string representation of the equilibrium constant."""
        if self.id is None:
            warnings.warn("No parameter ID. Define the reaction ID first.")
            return None

        return "Keq_" + self.id

    @property
    def kr_str(self):
        """Return the string representation of the reverse rate constant."""
        if self.id is None:
            warnings.warn("No parameter ID. Define the reaction ID first.")
            return None

        return "kr_" + self.id

    @property
    def kf(self):
        """Alias for the :attr:`forward_rate_constant`."""
        return self.forward_rate_constant

    @kf.setter
    def kf(self, value):
        """Alias for the :attr:`forward_rate_constant`."""
        self.forward_rate_constant = value

    @property
    def kr(self):
        """Alias for the :attr:`reverse_rate_constant`."""
        return self.reverse_rate_constant

    @kr.setter
    def kr(self, value):
        """Alias for the :attr:`reverse_rate_constant`."""
        self.reverse_rate_constant = value

    @property
    def Keq(self):
        """Alias for the :attr:`equilibrium_constant`."""
        return self.equilibrium_constant

    @Keq.setter
    def Keq(self, value):
        """Alias for the :attr:`equilibrium_constant`."""
        self.equilibrium_constant = value

    @property
    def S(self):
        """Alias for the :attr:`stoichiometry`."""
        return self.stoichiometry

    @property
    def v(self):
        """Alias for the :attr:`steady_state_flux`."""
        return self.steady_state_flux

    @v.setter
    def v(self, value):
        """Alias for the :attr:`steady_state_flux`."""
        self.steady_state_flux = value

    def reverse_stoichiometry(
        self,
        inplace=False,
        reverse_parameters=False,
        reverse_bounds=True,
        reverse_flux=True,
    ):
        """Reverse the stoichiometry of the reaction.

        Reversing the stoichiometry will turn the products into the reactants
        and the reactants into the products.

        Notes
        -----
        To avoid errors when reversing the reaction equilibrium constant:

            * If ``self.equilibrium_constant=0.`` then
              ``new_reaction.equilibrium_constant=float("inf")``
            * If ``self.equilibrium_constant=float("inf")`` then
              ``new_reaction.equilibrium_constant=0.``

        Parameters
        ----------
        inplace : bool
            If ``True``, modify the reaction directly. Otherwise a new reaction
            is created, modified, and returned.
        reverse_parameters : bool
            If ``True`` then also switch the reaction rate constants and
            inverse the equilibrium constants such that::

                new_reaction.forward_rate_constant = self.reverse_rate_constant
                new_reaction.reverse_rate_constant = self.forward_rate_constant
                new_reaction.equilibrium_constant = 1/self.equilibrium_constant

            Default is ``False``.
        reverse_bounds : bool
            If ``True`` then also switch the lower and upper bounds with one
            another such that::

                new_reaction.bounds = (-self.upper_bound, -self.lower_bound)

            Default is ``True``.
        reverse_flux: bool
             If ``True`` then also switch the direction of the flux such that::

                new_reaction.steady_state_flux = -self.steady_state_flux

            Default is ``True``.

        Returns
        -------
        new_reaction : MassReaction
            Returns the original :class:`MassReaction` if ``inplace=True``.
            Otherwise return a modified copy of the original reaction.

        """
        if inplace:
            new_reaction = self
        else:
            new_reaction = self.copy()

        for met, coeff in iteritems(new_reaction.metabolites):
            new_reaction._metabolites[met] = -1 * coeff

        if reverse_parameters:
            kf, kr, Keq = (self.kf, self.kr, self.Keq)
            new_reaction.kf = kr
            new_reaction.kr = kf
            if Keq == float("inf"):
                Keq = 0
            elif Keq == 0:
                Keq = float("inf")
            else:
                Keq = 1 / Keq
            new_reaction.Keq = round(Keq, MASSCONFIGURATION.decimal_precision)

        if reverse_bounds and self.bounds != (None, None):
            if self.lower_bound is None:
                new_reaction.bounds = (-self.upper_bound, None)
            elif self.upper_bound is None:
                new_reaction.bounds = (None, -self.lower_bound)
            else:
                new_reaction.bounds = (-self.upper_bound, -self.lower_bound)

        if reverse_flux and self.steady_state_flux is not None:
            new_reaction.steady_state_flux = -self.steady_state_flux

        new_reaction.get_mass_action_rate(update_reaction=True)

        return new_reaction

    def get_mass_action_rate(
        self, rate_type=1, update_reaction=False, destructive=False
    ):
        """Get the mass action rate law for the reaction.

        Parameters
        ----------
        rate_type : int
            The type of rate law to return. Must be 1, 2, or 3.

                * Type 1 will utilize the :attr:`forward_rate_constant` and the
                  :attr:`equilibrium_constant`.
                * Type 2 will utilize the :attr:`forward_rate_constant` and the
                  :attr:`reverse_rate_constant`.
                * Type 3 will utilize the :attr:`equilibrium_constant` and the
                  :attr:`reverse_rate_constant`.

            Default is ``1``.
        update_reaction : bool
            Whether to update the :attr:`MassReaction.rate` attribute
            in addition to returning the rate law. Default is ``False``.
        destructive : bool
            If ``True`` and the reaction has a custom rate law in its
            associated model, then setting ``update_reaction=True`` will
            replace the rate law and remove the custom rate law and orphaned
            parameters from the model. Default is ``False``.

        Returns
        -------
        rate_expression : :class:`sympy.core.basic.Basic` or ``None``
            The rate law as a :mod:`sympy` expression. If the reaction has no
            metabolites associated, ``None`` will be returned.

        Warnings
        --------
        Setting ``update_reaction=True`` will not remove any associated
        custom rate laws from the model unless ``destructive=True`` as well.

        """
        rate = generate_mass_action_rate_expression(self, rate_type)
        if update_reaction:
            self._rate_type = rate_type

        if destructive:
            if self.model is not None and self in self.model.custom_rates:
                self.model.remove_custom_rate(self, remove_orphans=True)
            elif self.model is None:
                warnings.warn("No MassModel associated with this reaction.")
            else:
                warnings.warn("No custom rate associated with this reaction.")

        return rate

    def get_forward_mass_action_rate_expression(self, rate_type=None):
        """Get the forward mass action rate expression for the reaction.

        Parameters
        ----------
        rate_type : int, None
            The type of rate law to return. Must be 1, 2, or 3.

                * Type 1 and 2 will utilize the :attr:`forward_rate_constant`.
                * Type 3 will utilize the :attr:`equilibrium_constant` and the
                  :attr:`reverse_rate_constant`.

            If ``None``, the current rate type will be used.
            Default is ``None``.

        Returns
        -------
        fwd_rate : :class:`sympy.core.basic.Basic` or ``None``
            The forward rate as a :mod:`sympy` expression. If the reaction
            has no metabolites associated, ``None`` will be returned.

        """
        if rate_type is None:
            rate_type = getattr(self, "_rate_type")
        return generate_forward_mass_action_rate_expression(self, rate_type)

    def get_reverse_mass_action_rate_expression(self, rate_type=1):
        """Get the reverse mass action rate expression for the reaction.

        Parameters
        ----------
        rate_type : int, None
            The type of rate law to return. Must be 1, 2, or 3.

                * Type 1 will utilize the :attr:`forward_rate_constant` and the
                  :attr:`equilibrium_constant`.
                * Type 2 and 3 will utilize the :attr:`reverse_rate_constant`.

            If ``None``, the current rate type will be used.
            Default is ``None``.

        Returns
        -------
        rev_rate : :class:`sympy.core.basic.Basic` or ``None``
            The reverse rate as a :mod:`sympy` expression. If the reaction
            has no metabolites associated, ``None`` will be returned.

        """
        if rate_type is None:
            rate_type = getattr(self, "_rate_type")
        return generate_reverse_mass_action_rate_expression(self, rate_type)

    def get_mass_action_ratio(self):
        """Get the mass action ratio as a :mod:`sympy` expression.

        Returns
        -------
        :class:`sympy.core.basic.Basic`
            The mass action ratio as a sympy expression.

        """
        return generate_mass_action_ratio(self)

    def get_disequilibrium_ratio(self):
        """Get the disequilibrium ratio as a :mod:`sympy` expression.

        Returns
        -------
        :class:`sympy.core.basic.Basic`
            The disequilibrium ratio as a sympy expression.

        """
        return generate_disequilibrium_ratio(self)

    def copy(self):
        """Copy a reaction.

        The reaction parameters, referenced metabolites, and genes are also
        copied.
        """
        return super(MassReaction, self).copy()

    def get_coefficient(self, metabolite_id):
        """Return the coefficient of a metabolite in the reaction.

        Parameters
        ----------
        metabolite_id : str or MassMetabolite
            The :class:`~.MassMetabolite` or the string identifier of the
            metabolite whose coefficient is desired.

        """
        return super(MassReaction, self).get_coefficient(metabolite_id)

    def get_coefficients(self, metabolite_ids):
        r"""Return coefficients for a ``list`` of metabolites in the reaction.

        Parameters
        ----------
        metabolite_ids : iterable
            Iterable containing the :class:`~.MassMetabolite`\ s or
            their string identifiers.

        """
        return super(MassReaction, self).get_coefficient(metabolite_ids)

    def add_metabolites(self, metabolites_to_add, combine=True, reversibly=True):
        r"""Add metabolites and their coefficients to the reaction.

        If the final coefficient for a metabolite is 0 then it is removed from
        the reaction.

        The change is reverted upon exit when using the :class:`~.MassModel`
        as a context.

        Notes
        -----
        * A final coefficient of < 0 implies a reactant and a final
          coefficient of > 0 implies a product.

        * Extends :meth:`~cobra.core.reaction.Reaction.add_metabolites` of the
          :class:`cobra.Reaction <cobra.core.reaction.Reaction>` by first
          ensuring that the metabolites to be added are
          :class:`.MassMetabolite`\ s and not
          :class:`cobra.Metabolites <cobra.core.metabolite.Metabolite>`.
          and error message raised reflects the :mod:`mass` object.

        * If a :class:`cobra.Metabolite <cobra.core.metabolite.Metabolite>` is
          provided. a warning is raised and a :class:`.MassMetabolite`
          will be instantiated using the
          :class:`cobra.Metabolite <cobra.core.metabolite.Metabolite>`.

        Parameters
        ----------
        metabolites_to_add : dict
            A ``dict`` with :class:`.MassMetabolite`\ s or metabolite
            identifiers as keys and stoichiometric coefficients as values. If
            keys are strings (id of a metabolite), the reaction must already
            be part of a :class:`~.MassModel` and a metabolite with the given
            id must already exist in the :class:`~.MassModel`.
        combine : bool
            Describes the behavior of existing metabolites.
            If ``True``, the metabolite coefficients are combined together.
            If ``False`` the coefficients are replaced.
        reversibly : bool
            Whether to add the change to the context to make the change
            reversible (primarily intended for internal use).

        See Also
        --------
        :meth:`subtract_metabolites`

        """
        # Ensure metabolites to be added are all MassMetabolites
        # or string identifiers of metabolites.
        for met in list(metabolites_to_add):
            if isinstance(met, (MassMetabolite, string_types)):
                # No need to change MassMetabolite objects
                continue
            elif isinstance(met, Metabolite):
                # Convert metabolite to a MassMetabolite and raise a warning
                warnings.warn(
                    "'{0}' is not a mass.MassMetabolite, therefore "
                    "converting metabolite before adding.".format(str(met))
                )
                mass_met = MassMetabolite(met)
                metabolites_to_add[mass_met] = metabolites_to_add.pop(met)
            else:
                # Input not recognized.
                raise TypeError("Unrecognized input {0}".format(str(met)))
        try:
            super(MassReaction, self).add_metabolites(
                metabolites_to_add, combine, reversibly
            )
        except ValueError as e:
            for cobra_obj_str in ["Reaction", "Metabolite"]:
                e = str(e).replace(cobra_obj_str, "Mass" + cobra_obj_str)
            raise ValueError(e)

    def subtract_metabolites(self, metabolites, combine=True, reversibly=True):
        r"""Subtract metabolites and their coefficients from the reaction.

        This function will 'subtract' metabolites from a reaction by adding
        the given metabolites with ``-1 * coeffcient``. If the final
        coefficient for a metabolite is 0, the metabolite is removed from the
        reaction.

        The change is reverted upon exit when using the :class:`~.MassModel`
        as a context.

        Notes
        -----
        * A final coefficient of < 0 implies a reactant and a final
          coefficient of > 0 implies a product.

        * Extends :meth:`cobra.core.reaction.Reaction.subtract_metabolites` of
          the :class:`cobra.Reaction <cobra.core.reaction.Reaction>` by first
          ensuring that the metabolites to be added are
          :class:`.MassMetabolite`\ s and not
          :class:`cobra.Metabolites <cobra.core.metabolite.Metabolite>`.
          and error message raised reflects the :mod:`mass` object.

        * If a :class:`cobra.Metabolite <cobra.core.metabolite.Metabolite>` is
          provided. a warning is raised and a :class:`.MassMetabolite`
          will be instantiated using the
          :class:`cobra.Metabolite <cobra.core.metabolite.Metabolite>`.

        Parameters
        ----------
        metabolites : dict
            A ``dict`` with :class:`~.MassMetabolite`\ s or their identifiers
            as keys and stoichiometric coefficients as values. If keys are
            strings (id of a metabolite), the reaction must already be part of
            a :class:`~.MassModel` and a metabolite with the given id must
            already exist in the :class:`~.MassModel`.
        combine : bool
            Describes the behavior of existing metabolites.
            If ``True``, the metabolite coefficients are combined together.
            If ``False`` the coefficients are replaced.
        reversibly : bool
            Whether to add the change to the context to make the change
            reversible (primarily intended for internal use).

        See Also
        --------
        :meth:`add_metabolites`

        """
        super(MassReaction, self).subtract_metabolites(metabolites, combine, reversibly)

    def build_reaction_string(self, use_metabolite_names=False):
        """Generate a human readable string to represent the reaction.

        Notes
        -----
        Overrides :meth:`~cobra.core.reaction.Reaction.build_reaction_string`
        of the :class:`cobra.Reaction <cobra.core.reaction.Reaction>` so that
        the reaction arrow depends on :attr:`MassReaction.reversible` rather
        than the inherited
        :attr:`cobra.Reaction.reversibility <cobra.core.reaction.Reaction.reversibility>`
        attribute.

        Parameters
        ----------
        use_metabolite_names : bool
            If ``True``, use the metabolite names instead of their identifiers.
            Default is ``False``.

        Returns
        -------
        reaction_string : str
            A string representation of the reaction.

        """  # noqa: E501

        def _format(number):
            return "" if number == 1 else str(number).rstrip(".") + " "

        id_type = "id"
        if use_metabolite_names:
            id_type = "name"
        # Seperate reactants and products
        reactant_bits = []
        product_bits = []
        for metab in sorted(self._metabolites, key=attrgetter("id")):
            coefficient = self._metabolites[metab]
            metab_name = str(getattr(metab, id_type))
            if coefficient >= 0:
                product_bits.append(_format(coefficient) + metab_name)
            else:
                reactant_bits.append(_format(abs(coefficient)) + metab_name)

        # Create reaction string
        reaction_string = " + ".join(reactant_bits)
        if self.reversible:
            reaction_string += " <=> "
        else:
            reaction_string += " --> "
        reaction_string += " + ".join(product_bits)

        return reaction_string

    def check_mass_balance(self):
        """Compute the mass and charge balances for the reaction.

        Returns
        -------
        dict
            Returns a ``dict`` of ``{element: amount}`` for unbalanced
            elements, with the "charge" treated as an element in the ``dict``.
            For a balanced reaction, an empty ``dict`` is returned.

        """
        return super(MassReaction, self).check_mass_balance()

    def build_reaction_from_string(
        self,
        reaction_str,
        verbose=True,
        fwd_arrow=None,
        rev_arrow=None,
        reversible_arrow=None,
        term_split="+",
    ):
        """Build reaction from reaction equation ``reaction_str`` using parser.

        Takes a string representation of the reaction and uses the
        specifications supplied in the optional arguments to infer a set of
        metabolites, metabolite compartments, and stoichiometries for the
        reaction. It also infers the refversibility of the reaction from the
        reaction arrow.

        For example:

            * ``'A + B <=> C'`` for reversible reactions, A & B are reactants.
            * ``'A + B --> C'`` for irreversible reactions, A & B are
              reactants.
            * ``'A + B <-- C'`` for irreversible reactions, A & B are products.

        The change is reverted upon exit when using the :class:`~.MassModel`
        as a context.

        Notes
        -----
        Extends
        :meth:`~cobra.core.reaction.Reaction.build_reaction_from_string`
        of the :class:`cobra.Reaction <cobra.core.reaction.Reaction>`
        in order to change how the irreversible backwards arrow is
        interpreted, affecting the assignment of reactants and products
        rather than how the bounds are set.

        Parameters
        ----------
        reaction_str: str
            A string containing the reaction formula (equation).
        verbose: bool
            Setting the verbosity of the function. Default is ``True``.
        fwd_arrow: re.compile, None
            For forward irreversible reaction arrows. If ``None``, the
            arrow is expected to be ``'-->'`` or ``'==>'``.
        rev_arrow: re.compile, None
            For backward irreversible reaction arrows. If ``None``, the
            arrow is expected to be ``'<--'`` or ``'<=='``.
        reversible_arrow: re.compile, None
            For reversible reaction arrows. If ``None``, the arrow is expected
            to be ``'<=>'`` or ``'<->'``.
        term_split: str
            Dividing individual metabolite entries. Default is ``"+"``.

        """
        is_reversible = self.reversible
        with warnings.catch_warnings():
            warnings.filterwarnings(
                "ignore",
                ".*is not a mass.MassMetabolite, therefore" " converting metabolite.*",
            )
            super(MassReaction, self).build_reaction_from_string(
                reaction_str,
                verbose,
                fwd_arrow,
                rev_arrow,
                reversible_arrow,
                term_split,
            )

        reversible_arrow_finder = (
            _reversible_arrow_finder
            if reversible_arrow is None
            else re.compile(re.escape(reversible_arrow))
        )

        reverse_arrow_finder = (
            _reverse_arrow_finder
            if rev_arrow is None
            else re.compile(re.escape(rev_arrow))
        )

        if not is_reversible and reversible_arrow_finder.search(reaction_str):
            warnings.warn(
                "Reaction '{0}' was previously set as `reversible=False`, but "
                "now has `reversible=True` due to the inferred reaction arrow "
                "in `reaction_str`.".format(self.id)
            )

        if reverse_arrow_finder.search(
            reaction_str
        ) and not reversible_arrow_finder.search(reaction_str):
            # Reverse the stoichiometry of the reaction
            self.reverse_stoichiometry(inplace=True, reverse_bounds=True)

    def knock_out(self):
        """Knockout reaction by setting its bounds to zero."""
        super(MassReaction, self).knock_out()

    # Internal
    def _cobra_to_mass_repair(self):
        """Convert associated cobra.Metabolites to MassMetabolites for self.

        Warnings
        --------
        This method is intended for internal use only.

        """
        metabolites = {}
        genes = set()
        model = self.model
        if self.metabolites:
            for metabolite in self.metabolites:
                # See if there is an associated MassModel with MassMetabolites
                # that already exist to add to the reaction
                if (
                    model.__class__.__name__ == "MassModel"
                    and str(metabolite) in model.metabolites
                ):
                    mass_met = model.metabolites.get_by_id(str(metabolite))
                else:
                    # Otherewise create a new MassMetabolite
                    mass_met = MassMetabolite(metabolite)
                # Set MassMetabolite stoichiometry
                metabolites[mass_met] = self.metabolites[metabolite]
        if self.genes:
            for gene in list(self.genes):
                if model.__class__.__name__ == "MassModel" and str(gene) in model.genes:
                    gene = model.genes.get_by_id(str(gene))
                else:
                    gene = gene.copy()
            genes.add(gene)
        # Remove cobra Metabolites and add MassMetabolites
        setattr(self, "_metabolites", metabolites)
        # Remove old genes and add new ones
        setattr(self, "_genes", genes)
        # Remove model reference
        setattr(self, "_model", None)
        self._update_awareness()

    def _associate_gene(self, cobra_gene):
        """Associates a :class:`~cobra.core.gene.Gene` with the reaction.

        Parameters
        ----------
        cobra_gene: Gene
            :class:`~cobra.core.gene.Gene` to be assoicated with the reaction.

        Warnings
        --------
        This method is intended for internal use only.

        """
        super(MassReaction, self)._associate_gene(cobra_gene)

    def _dissociate_gene(self, cobra_gene):
        """Dissociates a :class:`~cobra.core.gene.Gene` with the reaction.

        Parameters
        ----------
        cobra_gene: Gene
            :class:`~cobra.core.gene.Gene` to be disassociated with the
            reaction.

        Warnings
        --------
        This method is intended for internal use only.

        """
        super(MassReaction, self)._dissociate_gene(cobra_gene)

    def _make_boundary_metabolites(self):
        """Make the boundary metabolite.

        Warnings
        --------
        This method is intended for internal use only.

        """
        bc_metabolites = []
        for metabolite in list(self.metabolites):
            bc_metabolite = metabolite._remove_compartment_from_id_str()
            bc_metabolite += "_" + str(
                next(iter(MASSCONFIGURATION.boundary_compartment))
            )
            bc_metabolites += [bc_metabolite]

        return bc_metabolites

    def _repr_html_(self):
        """HTML representation of the overview for the MassReaction.

        Warnings
        --------
        This method is intended for internal use only.

        """
        return """
            <table>
                <tr>
                    <td><strong>Reaction identifier</strong></td>
                    <td>{id}</td>
                </tr><tr>
                    <td><strong>Name</strong></td>
                    <td>{name}</td>
                </tr><tr>
                    <td><strong>Memory address</strong></td>
                    <td>{address}</td>
                </tr><tr>
                    <td><strong>Subsystem</strong></td>
                    <td>{subsystem}</td>
                </tr><tr>
                    <td><strong>Kinetic Reversibility</strong></td>
                    <td>{reversibility}</td>
                </tr><tr>
                    <td><strong>Stoichiometry</strong></td>
                    <td>
                        <p style='text-align:right'>{stoich_id}</p>
                        <p style='text-align:right'>{stoich_name}</p>
                    </td>
                </tr><tr>
                    <td><strong>GPR</strong></td><td>{gpr}</td>
                </tr><tr>
                    <td><strong>Bounds</strong></td><td>({lb}, {ub})</td>
                </tr>
            </table>
        """.format(
            id=format_long_string(self.id, 100),
            name=format_long_string(self.name, 100),
            address="0x0%x" % id(self),
            subsystem=self.subsystem,
            reversibility=self._reversible,
            stoich_id=format_long_string(self.build_reaction_string(), 200),
            stoich_name=format_long_string(self.build_reaction_string(True), 200),
            gpr=format_long_string(self.gene_reaction_rule, 200),
            lb=self.lower_bound,
            ub=self.upper_bound,
        )

    # Dunders
    def __copy__(self):
        """Create a copy of the MassReaction.

        Warnings
        --------
        This method is intended for internal use only.

        """
        return copy(super(MassReaction, self))

    def __deepcopy__(self, memo):
        """Create a deepcopy of the MassReaction.

        Warnings
        --------
        This method is intended for internal use only.

        """
        return deepcopy(super(MassReaction, self), memo)

    def __str__(self):
        """Create an id string with the stoichiometry.

        Warnings
        --------
        This method is intended for internal use only.

        """
        return "{id}: {stoichiometry}".format(
            id=self.id, stoichiometry=self.build_reaction_string()
        )

    def __dir__(self):
        """Override default dir() implementation to list only public items.

        Warnings
        --------
        This method is intended for internal use only.

        """
        return get_public_attributes_and_methods(self)


__all__ = ("MassReaction",)
