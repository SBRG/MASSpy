# -*- coding: utf-8 -*-
"""
MassReaction is a class for holding information regarding reactions.

The :class:`MassReaction` class inherits and extends the
:class:`~cobra.core.reaction.Reaction` class in :mod:`cobra`. It contains
additional information required for simulations and other :mod:`mass`
functions and workflows.
"""
import re
from collections import defaultdict
from copy import copy, deepcopy
from functools import partial
from operator import attrgetter
from warnings import warn

from cobra.core.gene import Gene, ast2str, eval_gpr, parse_gpr
from cobra.core.object import Object
from cobra.core.reaction import (
    _forward_arrow_finder, _reverse_arrow_finder, _reversible_arrow_finder,
    and_or_search, compartment_finder, gpr_clean, uppercase_AND, uppercase_OR)
from cobra.util.context import get_context, resettable

from six import iteritems, iterkeys, itervalues, string_types

from sympy import Symbol

from mass.core.mass_configuration import MassConfiguration
from mass.core.mass_metabolite import MassMetabolite
from mass.util.expressions import (
    generate_disequilibrium_ratio, generate_foward_mass_action_rate_expression,
    generate_mass_action_rate_expression, generate_mass_action_ratio,
    generate_reverse_mass_action_rate_expression)
from mass.util.util import (
    ensure_non_negative_value, get_object_attributes,
    get_subclass_specific_attributes)

MASSCONFIGURATION = MassConfiguration()


class MassReaction(Object):
    """Class for holding kinetic information regarding a biochemical reaction.

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

    Attributes
    ----------
    steady_state_flux : float
        The stored (typically steady state) flux for the reaction. Stored flux
        values can be accessed for operations such as PERC calculations.

    """

    def __init__(self, id_or_reaction=None, name="", subsystem="",
                 reversible=True, steady_state_flux=None):
        """Initialize the MassReaction."""
        # Get the identifer and initialize
        super(MassReaction, self).__init__(
            getattr(id_or_reaction, "id", id_or_reaction), name)
        if isinstance(id_or_reaction, MassReaction):
            # Instiantiate a new MassReaction with state identical to
            # the provided MassReaction object.
            self.__dict__.update(id_or_reaction.__dict__)
        else:
            self.subsystem = subsystem
            self._reversible = reversible
            self.steady_state_flux = steady_state_flux
            # Rate and equilibrium constant parameters for reactions. The
            # reverse and equilibrium constants for irreversible reactions are
            # set here.  Upper and lower bounds are also set.
            self._forward_rate_constant = None
            if self._reversible:
                self._reverse_rate_constant = None
                self._equilibrium_constant = None
                self._lower_bound = MASSCONFIGURATION.lower_bound
            else:
                self._reverse_rate_constant = MASSCONFIGURATION.irreversible_kr
                self._equilibrium_constant = MASSCONFIGURATION.irreversible_Keq
                self._lower_bound = 0.
            self._upper_bound = MASSCONFIGURATION.upper_bound
            # Set objective coefficient and variable kind
            self.objective_coefficient = 0.
            self.variable_kind = 'continuous'
            # Rate type and law and as a sympy expression for simulation.
            self._rtype = 1
            # A dict of metabolites and their stoichiometric coefficients.
            self._metabolites = {}
            # The compartments where partaking metabolites are located.
            self._compartments = None
            # The associated MassModel.
            self._model = None
            # The genes associated with the reaction.
            self._genes = set()
            self._gene_reaction_rule = ""

    # Public
    @property
    def reversible(self):
        """Get or set the kinetic reversibility of the reaction.

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

        if reversible != self.reversible:
            self._reversible = reversible
            if reversible:
                setattr(self, "_reverse_rate_constant", None)
                setattr(self, "_equilibrium_constant", None)
            else:
                setattr(self, "_reverse_rate_constant",
                        MASSCONFIGURATION.irreversible_kr)
                setattr(self, "_equilibrium_constant",
                        MASSCONFIGURATION.irreversible_Keq)

    @property
    def forward_rate_constant(self):
        """Get or set the forward rate constant (kf) of the reaction.

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
    def forward_rate_constant(self, value):
        """Set the forward rate constant (kf) of the reaction."""
        value = ensure_non_negative_value(value)
        setattr(self, "_forward_rate_constant", value)

    @property
    def reverse_rate_constant(self):
        """Get or set the reverse rate constant (kr) of the reaction.

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
    def reverse_rate_constant(self, value):
        """Set the reverse rate constant (kr) of the reaction."""
        if self.reversible:
            value = ensure_non_negative_value(value)
            setattr(self, "_reverse_rate_constant", value)
        else:
            warn("Cannot set the reverse rate constant for an irreversible "
                 "reaction '{0}'".format(self.id))

    @property
    def equilibrium_constant(self):
        """Get or set the equilibrium constant (Keq) of the reaction.

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
    def equilibrium_constant(self, value):
        """Set the equilibrium constant (Keq) of the reaction."""
        if self._reversible:
            value = ensure_non_negative_value(value)
            setattr(self, "_equilibrium_constant", value)
        else:
            warn("Cannot set the equilibrium constant for an irreversible "
                 "reaction '{0}'".format(self.id))

    @property
    def parameters(self):
        """Return a dict of rate and equilibrium constants.

        Notes
        -----
        The :attr:`reverse_rate_constant` is only included for reversible
        reactions. Additionally, only rate and equilibrium constants are
        accessed here. Steady state fluxes can be accessed through the
        :attr:`steady_state_flux` attribute, and custom parameters can only be
        accessed through the model.

        """
        keys = [self.kf_str, self.Keq_str]
        attrs = ["_forward_rate_constant", "_equilibrium_constant"]
        # Return reverse rate constants for reversible reactions.
        if self.reversible:
            keys += [self.kr_str]
            attrs += ["_reverse_rate_constant"]
        parameters = {key: getattr(self, attr)
                      for key, attr in zip(keys, attrs)
                      if getattr(self, attr) is not None}

        return parameters

    @property
    def metabolites(self):
        """Return the metabolites of a reaction as a read only copy."""
        return getattr(self, "_metabolites").copy()

    @property
    def reactants(self):
        """Return a list of reactants for the reaction."""
        return [m for m, c in iteritems(self._metabolites) if c < 0]

    @property
    def products(self):
        """Return a list of products for the reaction."""
        return [m for m, c in iteritems(self._metabolites) if c >= 0]

    @property
    def stoichiometry(self):
        """Return a list of containing the stoichiometry for the reaction."""
        return [c for c in itervalues(self._metabolites)]

    @property
    def rate(self):
        """Return the current rate as a :mod:`sympy` expression.

        If reaction has a custom rate in its associated :class:`~.MassModel`,
        the custom rate will be returned instead.
        """
        if self.model is not None and self in self.model.custom_rates:
            rate = self._model.custom_rates[self]
        else:
            rate = self.get_mass_action_rate(rtype=self._rtype,
                                             update_reaction=True)
        return rate

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
        build_reaction_from_string: Base function for setter method.

        """
        return self.build_reaction_string()

    @reaction.setter
    def reaction(self, reaction_str):
        """Set the reaction using a human readable string."""
        return self.build_reaction_from_string(reaction_str)

    @property
    def compartments(self):
        """Return the set of compartments where the metabolites are located."""
        setattr(self, "_compartments", set(met.compartment
                                           for met in self._metabolites
                                           if met is not None))
        return getattr(self, "_compartments")

    @property
    def boundary(self):
        """Determine whether or not the reaction is a boundary reaction.

        Will return ``True`` if the reaction has no products or no reactants
        and only one metabolite.

        Notes
        -----
        These are reactions with a sink or a source term (e.g. 'A --> ')

        """
        return (len(self.metabolites) == 1
                and not (self.reactants and self.products))

    @property
    def boundary_metabolite(self):
        """Return an 'boundary' metabolite for bounary reactions.

        Notes
        -----
        The 'boundary_metabolite' represents the metabolite that corresponds
        to the empty part of a boundary reaction through a string. It's primary
        use is for setting of the :attr:`~.MassModel.boundary_conditions`
        without creating a :class:`~.MassMetabolite` object. Therefore it is
        not counted as a metabolite.

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
                bc_metabolite += "_" + str(next(iter(
                    MASSCONFIGURATION.boundary_compartment)))
        else:
            bc_metabolite = None

        return bc_metabolite

    @property
    def genes(self):
        """Return a frozenset of the genes associated with the reaction."""
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
        return getattr(self, "_gene_reaction_rule")

    @gene_reaction_rule.setter
    def gene_reaction_rule(self, new_rule):
        """Set the gene reaction rule of a reaction using a string."""
        if get_context(self):
            warn("Context management not implemented for "
                 "gene reaction rules")

        self._gene_reaction_rule = new_rule.strip()
        try:
            _, gene_names = parse_gpr(self._gene_reaction_rule)
        except (SyntaxError, TypeError):
            if "AND" in new_rule or "OR" in new_rule:
                warn("uppercase AND/OR found in rule '%s' for '%s'" %
                     (new_rule, repr(self)))
                new_rule = uppercase_AND.sub("and", new_rule)
                new_rule = uppercase_OR.sub("or", new_rule)
                self.gene_reaction_rule = new_rule
                return
            warn("malformed gene_reaction_rule '%s' for %s" %
                 (new_rule, repr(self)))
            tmp_str = and_or_search.sub('', self._gene_reaction_rule)
            gene_names = set((gpr_clean.sub(' ', tmp_str).split(' ')))
        if '' in gene_names:
            gene_names.remove('')
        old_genes = self._genes
        if self._model is None:
            self._genes = {Gene(i) for i in gene_names}
        else:
            model_genes = self._model.genes
            self._genes = set()
            for id in gene_names:
                if model_genes.has_id(id):
                    self._genes.add(model_genes.get_by_id(id))
                else:
                    new_gene = Gene(id)
                    new_gene._model = self._model
                    self._genes.add(new_gene)
                    model_genes.append(new_gene)

        # Make the genes aware that it is involved in this reaction
        for g in self._genes:
            g._reaction.add(self)

        # make the old genes aware they are no longer involved in this reaction
        for g in old_genes:
            if g not in self._genes:  # if an old gene is not a new gene
                try:
                    g._reaction.remove(self)
                except KeyError:
                    warn("could not remove old gene %s from reaction %s" %
                         (g.id, self.id))

    @property
    def gene_name_reaction_rule(self):
        """Display gene_reaction_rule with names.

        Warnings
        --------
        Do NOT use this string for computation. It is intended to give a
        representation of the rule using more familiar gene names instead
        of the often cryptic ids.

        """
        names = {i.id: i.name for i in self._genes}
        ast = parse_gpr(self._gene_reaction_rule)[0]

        return ast2str(ast, names=names)

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
        if self._model:
            tree, _ = parse_gpr(self.gene_reaction_rule)
            return eval_gpr(tree, {gene.id for gene in self.genes
                                   if not gene.functional})

        return True

    @property
    def flux_symbol(self):
        """Return the symbol representation for the reaction flux."""
        if self.id is not None:
            return Symbol("v_" + self.id)

        return None

    @property
    def all_parameter_ids(self):
        """Return a list of strings representing all non-custom parameters."""
        return [self.kf_str, self.Keq_str, self.kr_str, str(self.flux_symbol)]

    @property
    def lower_bound(self):
        """Get or set the lower bound of the reaction.

        Notes
        -----
        Infeasible combinations, such as a upper bound lower than the current
        lower bound will update the other bound.

        Parameters
        ----------
        value : float
            The new value for the lower bound.

        """
        return getattr(self, "_lower_bound")

    @lower_bound.setter
    @resettable
    def lower_bound(self, value):
        """Set the lower bound of the reaction."""
        if self.upper_bound < value:
            self.upper_bound = value

        setattr(self, "_lower_bound", value)

    @property
    def upper_bound(self):
        """Get or set the upper bound of the reaction.

        Notes
        -----
        Infeasible combinations, such as a upper bound lower than the current
        lower bound will update the other bound.

        Parameters
        ----------
        value : float
            The new value for the upper bound.

        """
        return getattr(self, "_upper_bound")

    @upper_bound.setter
    @resettable
    def upper_bound(self, value):
        """Set the upper bound of the reaction."""
        if self.lower_bound > value:
            self.lower_bound = value

        setattr(self, "_upper_bound", value)

    @property
    def bounds(self):
        """Get or set the bounds directly from a tuple.

        Convenience method for setting upper and lower bounds in one line
        using a tuple of lower and upper bound.

        When using a `HistoryManager` context, this attribute can be set
        temporarily, reversed when the exiting the context.

        Raises
        ------
        AssertionError
            Occurs when setting invalid bound values.

        """
        return self.lower_bound, self.upper_bound

    @bounds.setter
    @resettable
    def bounds(self, value):
        """Set the bounds directly from a tuple."""
        lower, upper = value
        self.lower_bound = lower
        self.upper_bound = upper

    @property
    def kf_str(self):
        """Return the string representation of the forward rate constant."""
        if self.id is not None:
            return "kf_" + self.id

        return None

    @property
    def Keq_str(self):
        """Return the string representation of the equilibrium constant."""
        if self.id is not None:
            return "Keq_" + self.id

        return None

    @property
    def kr_str(self):
        """Return the string representation of the reverse rate constant."""
        if self.id is not None:
            return "kr_" + self.id

        return None

    # Shorthands
    @property
    def kf(self):
        """Shorthand method to get the forward rate constant (kf)."""
        return self.forward_rate_constant

    @kf.setter
    def kf(self, value):
        """Shorthand method to set the forward rate constant (kf)."""
        self.forward_rate_constant = value

    @property
    def kr(self):
        """Shorthand method to get the reverse rate constant (kr)."""
        return self.reverse_rate_constant

    @kr.setter
    def kr(self, value):
        """Shorthand method to set the reverse rate constant (kr)."""
        self.reverse_rate_constant = value

    @property
    def Keq(self):
        """Shorthand method to get the equilibrium constant (Keq)."""
        return self.equilibrium_constant

    @Keq.setter
    def Keq(self, value):
        """Shorthand method to set the equilibrium constant (Keq)."""
        self.equilibrium_constant = value

    @property
    def S(self):
        """Shorthand method to get the reaction stoichiometry."""
        return self.stoichiometry

    @property
    def v(self):
        """Shorthand method to get the reaction steady state flux."""
        return self.steady_state_flux

    @v.setter
    def v(self, value):
        """Shorthand method to set the reaction steady state flux."""
        self.steady_state_flux = value

    def print_attributes(self, sep=r"\n", exclude_parent=False):
        r"""Print the attributes and properties of the :class:`MassReaction`.

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

    def reverse_stoichiometry(self, inplace=False):
        """Reverse the stoichiometry of the reaction.

        Reversing the stoichiometry will turn the products into the reactants
        and the reactants into the products.

        Notes
        -----
        Only the stoichiometry of the reaction is modified. The reaction
        parameters (e.g. rate and equilibrium constants) are not altered.

        Parameters
        ----------
        inplace : bool
            If ``True``, modify the reaction directly. Otherwise a new reaction
            is created, modified, and returned.

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

        return new_reaction

    def get_mass_action_rate(self, rtype=1, update_reaction=False):
        """Get the mass action rate law for the reaction.

        Parameters
        ----------
        rtype : int
            The type of rate law to display. Must be 1, 2, or 3.

                * Type 1 will utilize the :attr:`forward_rate_constant` and the
                  :attr:`equilibrium_constant`.
                * Type 2 will utilize the :attr:`forward_rate_constant` and the
                  :attr:`reverse_rate_constant`.
                * Type 3 will utilize the :attr:`equilibrium_constant` and the
                  :attr:`reverse_rate_constant`.

            Default is ``1``.
        update_reaction : bool
            Whether to update the reaction in addition to returning the
            rate law. Default is ``False``.

        Returns
        -------
        rate_expression : :class:`sympy.core.basic.Basic` or ``None``
            The rate law as a :mod:`sympy` expression. If the reaction has no
            metabolites associated, ``None`` will be returned.

        """
        return generate_mass_action_rate_expression(self, rtype,
                                                    update_reaction)

    def get_foward_mass_action_rate_expression(self, rtype=1):
        """Get the foward mass action rate expression for the reaction.

        Parameters
        ----------
        rtype : int
            The type of rate law to display. Must be 1, 2, or 3.

                * Type 1 and 2 will utilize the :attr:`forward_rate_constant`.
                * Type 3 will utilize the :attr:`equilibrium_constant` and the
                  :attr:`reverse_rate_constant`.

            Default is `1`.

        Returns
        -------
        fwd_rate : :class:`sympy.core.basic.Basic` or ``None``
            The forward rate as a :mod:`sympy` expression. If the reaction
            has no metabolites associated, ``None`` will be returned.

        """
        return generate_foward_mass_action_rate_expression(self, rtype)

    def get_reverse_mass_action_rate_expression(self, rtype=1):
        """Get the reverse mass action rate expression for the reaction.

        Parameters
        ----------
        rtype : int
            The type of rate law to display. Must be 1, 2, or 3.

                * Type 1 will utilize the :attr:`forward_rate_constant` and the
                  :attr:`equilibrium_constant`.
                * Type 2 and 3 will utilize the :attr:`reverse_rate_constant`.

            Default is `1`.

        Returns
        -------
        rev_rate : :class:`sympy.core.basic.Basic` or ``None``
            The reverse rate as a :mod`sympy` expression. If the reaction
            has no metabolites associated, ``None`` will be returned.

        """
        return generate_reverse_mass_action_rate_expression(self, rtype)

    def get_mass_action_ratio(self):
        """Get the mass action ratio as a :`mod`sympy` expression.

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

    def remove_from_model(self, remove_orphans=False):
        """Remove the reaction from the :class:`~.MassModel`.

        This removes all associations between a reaction, the associated
        model, metabolites, and genes.

        The change is reverted upon exit when using the :class:`~.MassModel`
        as a context.

        Parameters
        ----------
        remove_orphans: bool
            Remove orphaned genes and metabolites from the :class:`~.MassModel`
            as well.

        """
        return self._model.remove_reactions([self], remove_orphans)

    def copy(self):
        """Copy a reaction.

        The reaction parameters, referenced metabolites, and genes are also
        copied.
        """
        # No references to the MassModel when copying the MassReaction
        model = self._model
        setattr(self, "_model", None)
        for i in self._metabolites:
            setattr(i, "_model", None)
        for i in self._genes:
            setattr(i, "_model", None)

        # The reaction can now be copied
        reaction_copy = deepcopy(self)
        # Restore references for the original MassReaction
        setattr(self, "_model", model)
        for i in self._metabolites:
            setattr(i, "_model", model)
        for i in self._genes:
            setattr(i, "_model", model)

        return reaction_copy

    def get_coefficient(self, metabolite_id):
        """Return the coefficient of a metabolite in the reaction.

        Parameters
        ----------
        metabolite_id : str or MassMetabolite
            The :class:`~.MassMetabolite` or the string identifier of the
            metabolite whose coefficient is desired.

        """
        if isinstance(metabolite_id, MassMetabolite):
            return self._metabolites[metabolite_id]

        _id_to_mets = {m.id: m for m in self._metabolites}
        return self._metabolites[_id_to_mets[metabolite_id]]

    def get_coefficients(self, metabolite_ids):
        r"""Return the coefficients for a list of metabolites in the reaction.

        Parameters
        ----------
        metabolite_ids : iterable
            Iterable containing the :class:`~.MassMetabolite`\ s or
            their string identifiers.

        """
        return map(self.get_coefficient, metabolite_ids)

    def get_compartments(self):
        """Return a list of compartments where the metabolites are located."""
        return list(self.compartments)

    def add_metabolites(self, metabolites_to_add, combine=True,
                        reversibly=True):
        r"""Add metabolites and their coefficients to the reaction.

        If the final coefficient for a metabolite is 0 then it is removed from
        the reaction.

        The change is reverted upon exit when using the :class:`~.MassModel`
        as a context.

        Notes
        -----
        A final coefficient of < 0 implies a reactant and a final
        coefficient of > 0 implies a product.

        Parameters
        ----------
        metabolites_to_add : dict
            A ``dict`` with :class:`.MassMetabolite`\ s or metabolite
            identifiers as keys and stoichiometric coefficients as values. If
            keys are strings (id of a metabolite), the reaction must already
            be part of a :class:`~.MassModel` and a metabolite with the given
            id must already exist in the :class:`~.MassModel`.
        combine : bool
            If ``True``, the metabolite coefficients are combined together.
            Otherwise the coefficients are replaced.
        reversibly : bool
            Whether to add the change to the context to make the change
            reversible (primarily intended for internal use).

        Warnings
        --------
        A :class:`cobra.Metabolite <cobra.core.metabolite.Metabolite>` cannot
        be directly added to a :class:`MassReaction`. Instead, the metabolite
        must first be instantiated as a :class:`~.MassMetabolite`.

        See Also
        --------
        :meth:`subtract_metabolites`

        """
        old_coefficients = self.metabolites
        new_metabolites = []
        _id_to_metabolites = {x.id: x for x in self._metabolites}

        for metabolite, coefficient in iteritems(metabolites_to_add):
            # Make sure metabolites being added belong to the same model, or
            # else copy them.
            if isinstance(metabolite, MassMetabolite):
                if metabolite.model is not None and \
                   metabolite.model is not self._model:
                    metabolite = metabolite.copy()

            met_id = str(metabolite)
            # If a metabolite already exists in the reaction,
            # just add the coefficients.
            if met_id in _id_to_metabolites:
                reaction_metabolite = _id_to_metabolites[met_id]
                if combine:
                    self._metabolites[reaction_metabolite] += coefficient
                else:
                    self._metabolites[reaction_metabolite] = coefficient
            else:
                # If the reaction is in a MassModel, ensure a duplicate
                # metabolite is not added.
                if self._model:
                    try:
                        metabolite = self._model.metabolites.get_by_id(met_id)
                    except KeyError as e:
                        if isinstance(metabolite, MassMetabolite):
                            new_metabolites.append(metabolite)
                        else:
                            raise e
                elif isinstance(metabolite, string_types):
                    raise ValueError("Reaction '{0}' does not belong to a "
                                     "MassModel. Either add the reaction to a "
                                     "MassModel or use the MassMetabolite "
                                     "objects as keys instead of strings."
                                     .format(self.id))
                self._metabolites[metabolite] = coefficient
                # Make the metabolite aware of its involvement in the reaction.
                metabolite._reaction.add(self)

        model = self.model
        if model is not None:
            model.add_metabolites(new_metabolites)

        for metabolite, coefficient in list(iteritems(self._metabolites)):
            if coefficient == 0:
                # Make the metabolite aware of it no longer in the reaction.
                metabolite._reaction.remove(self)
                self._metabolites.pop(metabolite)

        context = get_context(self)
        if context and reversibly:
            if combine:
                # Just subtract the previously added metabolites
                context(partial(self.subtract_metabolites, metabolites_to_add,
                                combine=True, reversibly=False))
            else:
                # Reset the metabolites with add_metabolites
                mets_to_reset = {
                    key: old_coefficients[model.metabolites.get_by_any(key)[0]]
                    for key in iterkeys(metabolites_to_add)}

                context(partial(self.add_metabolites, mets_to_reset,
                                combine=False, reversibly=False))

    def subtract_metabolites(self, metabolites_to_subtract, combine=True,
                             reversibly=True):
        r"""Subtract metabolites and their coefficients from the reaction.

        This function will 'subtract' metabolites from a reaction by adding
        the given metabolites with ``-1 * coeffcient``. If the final
        coefficient for a metabolite is 0, the metabolite is removed from the
        reaction.

        The change is reverted upon exit when using the :class:`~.MassModel`
        as a context.

        Notes
        -----
        A final coefficient of < 0 implies a reactant and a final
        coefficient of > 0 implies a product.

        Parameters
        ----------
        metabolites_to_subtract : dict
            A ``dict`` with :class:`~.MassMetabolite`\ s or their identifiers
            as keys and stoichiometric coefficients as values. If keys are
            strings (id of a metabolite), the reaction must already be part of
            a :class:`~.MassModel` and a metabolite with the given id must
            already exist in the :class:`~.MassModel`.
        combine : bool
            If ``True``, the metabolite coefficients are combined together.
            Otherwise the coefficients are replaced.
        reversibly : bool
            Whether to add the change to the context to make the change
            reversible (primarily intended for internal use).

        Warnings
        --------
        A :class:`cobra.Metabolite <cobra.core.metabolite.Metabolite>` cannot
        be directly added to a :class:`MassReaction`. Instead, the metabolite
        must first be instantiated as a :class:`~.MassMetabolite`.

        See Also
        --------
        :meth:`add_metabolites`

        """
        self.add_metabolites({k: -v
                              for k, v in iteritems(metabolites_to_subtract)},
                             combine=combine, reversibly=reversibly)

    def build_reaction_string(self, use_metabolite_names=False):
        """Generate a human readable string to represent the reaction.

        Parameters
        ----------
        use_metabolite_names : bool
            If ``True``, use the metabolite names instead of their identifiers.
            Default is ``False``.

        Returns
        -------
        reaction_string : str
            A string representation of the reaction.

        """
        def format(number):
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
                product_bits.append(format(coefficient) + metab_name)
            else:
                reactant_bits.append(format(abs(coefficient)) + metab_name)

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
            elements, with the "charge" treated as an element in the dict.
            For a balanced reaction, an empty dict is returned.

        """
        reaction_element_dict = defaultdict(int)
        for metabolite, coeff in iteritems(self._metabolites):
            if metabolite.charge is not None:
                reaction_element_dict["charge"] += coeff * metabolite.charge
            if metabolite.elements is None:
                raise ValueError("No elements found in metabolite {0}"
                                 .format(metabolite.id))
            for element, amount in iteritems(metabolite.elements):
                reaction_element_dict[element] += coeff * amount
        # Filter out any 0 values.
        return {k: v for k, v in iteritems(reaction_element_dict) if v != 0}

    def build_reaction_from_string(self, reaction_str, verbose=True,
                                   fwd_arrow=None, rev_arrow=None,
                                   reversible_arrow=None, term_split="+"):
        """Build reaction from reaction equation ``reaction_str`` using parser.

        Takes a string representation of the reaction and uses the
        specifications supplied in the optional arguments to infer a set of
        metabolites, metabolite compartments, and stoichiometries for the
        reaction. It also infers the refversibility of the reaction from the
        reaction arrow.

        For example:

            * 'A + B <=> C' for reversible reactions, A & B are reactants.
            * 'A + B --> C' for irreversible reactions, A & B are reactants.
            * 'A + B <-- C' for irreversible reactions, A & B are products.

        The change is reverted upon exit when using the :class:`~.MassModel`
        as a context.

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
        # Set the arrows
        forward_arrow_finder = _forward_arrow_finder if fwd_arrow is None \
            else re.compile(re.escape(fwd_arrow))
        reverse_arrow_finder = _reverse_arrow_finder if rev_arrow is None \
            else re.compile(re.escape(rev_arrow))
        reversible_arrow_finder = _reversible_arrow_finder \
            if reversible_arrow is None \
            else re.compile(re.escape(reversible_arrow))
        if self._model is None:
            warn("no model found")
            model = None
        else:
            model = self._model
        found_compartments = compartment_finder.findall(reaction_str)
        if len(found_compartments) == 1:
            compartment = found_compartments[0]
            reaction_str = compartment_finder.sub("", reaction_str)
        else:
            compartment = ""

        # Parse reaction string to seperate reactants and products
        def split_reaction_str(arrow_match):
            left_str = reaction_str[:arrow_match.start()].strip()
            right_str = reaction_str[arrow_match.end():].strip()
            return left_str, right_str

        # Check for a reversible reaction
        arrow_match = reversible_arrow_finder.search(reaction_str)
        if arrow_match is not None:
            # Set reversibility, determine the reactants and products.
            self.reversible = True
            # Reactants are on the left, products are on the right.
            reactant_str, product_str = split_reaction_str(arrow_match)
        else:  # Irreversible reaction
            # Try the forward arrow
            arrow_match = forward_arrow_finder.search(reaction_str)
            if arrow_match is not None:
                # Set reversibility, determine the reactants and products.
                self.reversible = False
                # Reactants are on the left, products are on the right.
                reactant_str, product_str = split_reaction_str(arrow_match)
            else:
                # Try the reverse arrow
                arrow_match = reverse_arrow_finder.search(reaction_str)
                if arrow_match is not None:
                    # Set reversibility, determine the reactants and products.
                    self.reversible = False
                    # Reactants are on the right, products are on the left.
                    product_str, reactant_str = split_reaction_str(arrow_match)
                else:
                    raise ValueError("No suitable arrow found in '{0}'"
                                     .format(reaction_str))

        self.subtract_metabolites(self.metabolites, combine=True)

        for substr, factor in ((reactant_str, -1), (product_str, 1)):
            if not substr:
                continue

            for term in substr.split(term_split):
                term = term.strip()
                if term.lower() == "nothing":
                    continue
                if " " in term:
                    num_str, met_id = term.split()
                    num = float(num_str.lstrip("(").rstrip(")")) * factor
                else:
                    met_id = term
                    num = factor
                met_id += compartment
                try:
                    met = model.metabolites.get_by_id(met_id)
                except KeyError:
                    if verbose:
                        print("Unknown metabolite {0} created".format(met_id))
                    met = MassMetabolite(met_id)
                self.add_metabolites({met: num})

    def knock_out(self):
        """Knockout reaction by setting its bounds to zero."""
        self.lower_bound, self.upper_bound = (0, 0)

    # Internal
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
        self._genes.add(cobra_gene)
        cobra_gene._reaction.add(self)
        cobra_gene._model = self._model

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
        self._genes.discard(cobra_gene)
        cobra_gene._reaction.discard(self)

    def _set_id_with_model(self, value):
        """Set the id of the MassReaction to the associated MassModel.

        Warnings
        --------
        This method is intended for internal use only.

        """
        if value in self.model.reactions:
            raise ValueError("The model already contains a reaction with "
                             "the id: ", value)
        self._id = value
        self.model.reactions._generate_index()

    def _update_awareness(self):
        """Make species aware of their involvement with the reaction.

        Warnings
        --------
        This method is intended for internal use only.

        """
        for metab in self._metabolites:
            metab._reaction.add(self)
        for gene in self._genes:
            gene._reaction.add(self)

    def _make_boundary_metabolites(self):
        """Make the boundary metabolite.

        Warnings
        --------
        This method is intended for internal use only.

        """
        bc_metabolites = []
        for metabolite in list(self.metabolites):
            bc_metabolite = metabolite._remove_compartment_from_id_str()
            bc_metabolite += "_" + str(next(iter(
                MASSCONFIGURATION.boundary_compartment)))
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
                    <td><strong>Stoichiometry</strong></td>
                    <td>
                        <p style='text-align:right'>{stoich_id}</p>
                        <p style='text-align:right'>{stoich_name}</p>
                    </td>
                </tr><tr>
                    <td><strong>GPR</strong></td>
                    <td>{gpr}</td>
                </tr><tr>
                    <td><strong>Kinetic Reversibility</strong></td>
                    <td>{reversibility}</td>
                </tr>
            </table>
        """.format(id=self.id, name=self.name,
                   subsystem=self.subsystem, address='0x0%x' % id(self),
                   stoich_id=self.build_reaction_string(),
                   stoich_name=self.build_reaction_string(True),
                   gpr=self.gene_reaction_rule,
                   reversibility=self._reversible)

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

    def __setstate__(self, state):
        """Set state of MassReaction object upon unpickling.

        Probably not necessary to set _model as the MassModel that
        contains self sets the _model attribute for all metabolites and genes
        in the reaction.

        However, to increase performance speed, let the metabolite
        and gene know that they are employed in this reaction

        Warnings
        --------
        This method is intended for internal use only.

        """
        # These are necessary for old pickles which store attributes
        # which have since been superceded by properties.
        if "reaction" in state:
            state.pop("reaction")
        if "gene_reaction_rule" in state:
            state["_gene_reaction_rule"] = state.pop("gene_reaction_rule")
        if "lower_bound" in state:
            state['_lower_bound'] = state.pop('lower_bound')
        if "upper_bound" in state:
            state['_upper_bound'] = state.pop('upper_bound')

        self.__dict__.update(state)
        for x in state['_metabolites']:
            setattr(x, '_model', self._model)
            x._reaction.add(self)
        for x in state['_genes']:
            setattr(x, '_model', self._model)
            x._reaction.add(self)

    def __add__(self, other):
        """Add two reactions.

        The stoichiometry will be the combined stoichiometry of the two
        reactions, and the gene reaction rule will be both rules combined by an
        and. All other attributes (i.e. rate constants) will match those of
        the first reaction, Return a new MassReaction object.

        Warnings
        --------
        This method is intended for internal use only.

        """
        new_reaction = self.copy()
        new_reaction += other
        return new_reaction

    def __iadd__(self, other):
        """Add two reactions.

        The stoichiometry will be the combined stoichiometry of the two
        reactions, and the gene reaction rule will be both rules combined by an
        and. All other attributes (i.e. rate constants) will match those of the
        first reaction. Return the same MassReaction object with its updates.

        Warnings
        --------
        This method is intended for internal use only.

        """
        self.add_metabolites(other._metabolites, combine=True)
        gpr1 = self.gene_reaction_rule.strip()
        gpr2 = other.gene_reaction_rule.strip()
        if gpr1 != "" and gpr2 != "":
            self.gene_reaction_rule = ("(%s) and (%s)" %
                                       (self.gene_reaction_rule,
                                        other.gene_reaction_rule))
        elif gpr1 != "" and gpr2 == "":
            self.gene_reaction_rule = gpr1
        elif gpr1 == "" and gpr2 != "":
            self.gene_reaction_rule = gpr2
        return self

    def __sub__(self, other):
        """Subtract two reactions.

        The stoichiometry will be the combined stoichiometry of the two
        reactions, and the gene reaction rule will be both rules combined by an
        'and'. All other attributes (i.e. rate constants) will match those of
        the first reaction. Return a new MassReaction object.

        Warnings
        --------
        This method is intended for internal use only.

        """
        new_reaction = self.copy()
        new_reaction -= other
        return new_reaction

    def __isub__(self, other):
        """Subtract two reactions.

        The stoichiometry will be the combined stoichiometry of the two
        reactions. Returns the same MassReaction object with its updates.

        Warnings
        --------
        This method is intended for internal use only.

        """
        self.subtract_metabolites(other._metabolites, combine=True)
        return self

    def __mul__(self, coefficient):
        """Scale coefficients in a reaction by a given value.

        Returns a new MassReaction object
        E.g. A -> B becomes 2A -> 2B.

        Warnings
        --------
        This method is intended for internal use only.

        """
        new_reaction = self.copy()
        new_reaction *= coefficient
        return new_reaction

    def __imul__(self, coefficient):
        """Scale coefficients in a reaction by a given value.

        Returns the same MassReaction object with its updates
        E.g. A -> B becomes 2A -> 2B.

        Warnings
        --------
        This method is intended for internal use only.

        """
        self._metabolites = {met: value * coefficient
                             for met, value in iteritems(self._metabolites)}
        return self

    def __str__(self):
        """Create an id string with the stoichiometry.

        Warnings
        --------
        This method is intended for internal use only.

        """
        return "{id}: {stoichiometry}".format(
            id=self.id, stoichiometry=self.build_reaction_string())


__all__ = ("MassReaction",)
