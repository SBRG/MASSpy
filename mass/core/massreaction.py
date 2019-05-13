# -*- coding: utf-8 -*-
"""TODO Module Docstrings."""
from __future__ import absolute_import

import re
from collections import defaultdict
from copy import copy, deepcopy
from functools import partial
from operator import attrgetter
from warnings import warn

from six import iteritems, iterkeys, itervalues, string_types

from sympy import Symbol

from cobra.core.gene import Gene, ast2str, eval_gpr, parse_gpr
from cobra.core.object import Object
from cobra.core.reaction import (
    _forward_arrow_finder, _reverse_arrow_finder, _reversible_arrow_finder,
    and_or_search, compartment_finder, gpr_clean, uppercase_AND, uppercase_OR)
from cobra.util.context import get_context, resettable

from mass.core.massmetabolite import MassMetabolite
from mass.util.expressions import (
    generate_disequilibrium_ratio, generate_mass_action_ratio,
    generate_rate_law)


# Global
_INF = float("inf")
_BOUNDARY_PREFIX = "bc_"


class MassReaction(Object):
    """Class for holding kinetic information regarding a biochemical reaction.

    Parameters
    ----------
    id: str
        The identifier associated with the MassReaction.
    name: str, optional
        A human readable name for the reaction.
    subsystem: str, optional
        The subsystem where the reaction is meant to occur.
    reversible: bool, optional
        The kinetic reversibility of the reaction. Irreversible reactions have
        an equilibrium constant of infinity and a reverse rate constant of 0.
        If not provided, the reaction is assumed to be reversible.

    Attributes
    ----------
    steady_state_flux: float, optional
        The stored (typically steady state) flux for the reaction. Stored flux
        values can be accessed for operations such as PERC calculations.

    """

    def __init__(self, id=None, name="", subsystem="", reversible=True,
                 steady_state_flux=None):
        """Initialize the MassReaction Object."""
        super(MassReaction, self).__init__(id, name)
        self.subsystem = subsystem
        self._reversible = reversible
        self.steady_state_flux = steady_state_flux
        # Rate and equilibrium constant parameters for reactions. The reverse
        # and equilibrium constants for irreversible reactions are set here.
        # For cobra compatibility, lower and upper bounds are also set.
        self._forward_rate_constant = None
        if self._reversible:
            self._reverse_rate_constant = None
            self._equilibrium_constant = None
            self._lower_bound = -1000.
            self._upper_bound = 1000.
        else:
            self._reverse_rate_constant = 0.
            self._equilibrium_constant = _INF
            self._lower_bound = 0.
            self._upper_bound = 1000.

        self.objective_coefficient = 0.
        self.variable_kind = 'continuous'

        # Rate type and law and as a sympy expression for simulation.
        self._rtype = 1

        # A dictionary of metabolites and their stoichiometric coefficients.
        self._metabolites = {}

        # The compartments where partaking metabolites are located.
        self._compartments = None

        # The associated MassModel.
        self._model = None

        # The genes associated with the reaction.
        self._genes = set()
        self._gene_reaction_rule = ""

        # The Gibbs reaction energy associated with the reaction.
        self._gibbs_reaction_energy = None

    # Public
    @property
    def reversible(self):
        """Return the kinetic reversibility of the reaction."""
        return getattr(self, "_reversible")

    @reversible.setter
    def reversible(self, value):
        """Set the kinetic reversibility of the reaction.

        Warnings
        --------
        Changing the reversibility will reset the equilibrium constant and the
            reverse rate constant  to the defaults.

        """
        if not isinstance(value, bool):
            raise TypeError("value must be a bool")

        if value != self.reversible:
            self._reversible = value
            if value:
                setattr(self, "_reverse_rate_constant", None)
                setattr(self, "_equilibrium_constant", None)
            else:
                setattr(self, "_reverse_rate_constant", 0)
                setattr(self, "_equilibrium_constant", _INF)

    @property
    def forward_rate_constant(self):
        """Return the forward rate constant (kf) of the reaction."""
        return getattr(self, "_forward_rate_constant")

    @forward_rate_constant.setter
    def forward_rate_constant(self, value):
        """Set the forward rate constant (kf) of the reaction."""
        setattr(self, "_forward_rate_constant", value)

    @property
    def reverse_rate_constant(self):
        """Return the reverse rate constant (kr) of the reaction."""
        return getattr(self, "_reverse_rate_constant")

    @reverse_rate_constant.setter
    def reverse_rate_constant(self, value):
        """Set the reverse rate constant (kr) of the reaction.

        If reaction is not reversible, will warn the user instead.
        """
        if self.reversible:
            setattr(self, "_reverse_rate_constant", value)
        else:
            warn("Cannot set the reverse rate constant for an irreversible "
                 "reaction")

    @property
    def equilibrium_constant(self):
        """Return the equilibrium constant (Keq) of the reaction."""
        return getattr(self, "_equilibrium_constant")

    @equilibrium_constant.setter
    def equilibrium_constant(self, value):
        """Set the equilibrium constant (Keq) of the reaction.

        If reaction is not reversible, will warn the user instead.
        """
        if self._reversible:
            setattr(self, "_equilibrium_constant", value)
        else:
            warn("Cannot set the equilibrium constant for an irreversible "
                 "reaction")

    @property
    def parameters(self):
        """Return a dictionary of rate and equilibrium constants.

        Notes
        -----
        Reverse rate constants are only included for reversible reactions.
        Only rate and equilibrium constantx are accessed here. Steady state
            fluxes can be accessed through the steady_state_flux attribute,
            and custom parameters can only be accessed through the model.

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
        """Return the current rate for the reaction as a sympy expression.

        If reaction has a custom rate law in its associated MassModel, the
        custom rate law will be returned instead.
        """
        if self.model is not None and self in self.model.custom_rates:
            rate = self._model.custom_rates[self]
        else:
            rate = self.get_rate_law(rate_type=self._rtype, sympy_expr=True,
                                     update_reaction=True)

        return rate

    @property
    def model(self):
        """Return the MassModel associated with the reaction."""
        return getattr(self, "_model")

    @property
    def reaction(self):
        """Return the reaction as a human readable string."""
        return self.build_reaction_string()

    @reaction.setter
    def reaction(self, value):
        """Set the reaction using a human readable string.

         For example:
            'A + B <=> C' for reversible reactions, A & B are reactants.
            'A + B --> C' for irreversible reactions, A & B are reactants.
            'A + B <-- C' for irreversible reactions, A & B are products.


        Parameters
        ----------
        reaction_string: str
            String representation of the reaction.

        Notes
        -----
        The direction of the arrow is used to determine reversibility.

        Warnings
        --------
        Care must be taken when setting a reaction in this manner to ensure the
            reaction identifier matches those in an assoicated model.

        For forward and reversible arrows, reactants are searched on the left
            side of the arrow while products are searched on the right side.
            For reverse arrows, reactants are searched on the right side of the
            arrow while products are searched on the left side.

        """
        return self.build_reaction_from_string(value)

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

        Will return True if the reaction has no products or no reactants.

        Notes
        -----
        These are reactions with a sink or a source term (e.g. 'A --> ')

        """
        return (len(self.metabolites) == 1
                and not (self.reactants and self.products))

    @property
    def boundary_metabolite(self):
        """Return an 'boundary' metabolite for bounary reactions.

        Returns
        -------
        boundary_metabolite: str
            String representation of the boundary metabolite of the reaction,
            or None if the reaction is not considered a boundary reaction.

        See Also
        --------
        MassReaction.boundary

        """
        if self.boundary:
            bc_metabolite = _BOUNDARY_PREFIX + str(list(self.metabolites)[0])
        else:
            bc_metabolite = None

        return bc_metabolite

    @property
    def genes(self):
        """Return a frozenset of the genes associated with the reaction."""
        return frozenset(getattr(self, "_genes"))

    @property
    def gene_reaction_rule(self):
        """Return the gene reaction rule for the reaction as a string."""
        return getattr(self, "_gene_reaction_rule")

    @gene_reaction_rule.setter
    def gene_reaction_rule(self, new_rule):
        """Set the gene reaction rule of a reaction using a string.

        New genes will be associated with the reaction and old genes will be
        dissociated from the reaction.

        Parameters
        ----------
        new_rule: str
            String representation of the new reaction rule.

        """
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
        True if the gene-protein-reaction (GPR) rule is fulfilled for
            the reaction, or if the reaction does not have an assoicated
            MassModel. Otherwise returns False.

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
        else:
            return None

    @property
    def all_parameter_ids(self):
        """Return a list of strings representing all non-custom parameters."""
        return [self.kf_str, self.Keq_str, self.kr_str, str(self.flux_symbol)]

    @property
    def lower_bound(self):
        """Get the lower bound of the reaction."""
        return getattr(self, "_lower_bound")

    @lower_bound.setter
    @resettable
    def lower_bound(self, value):
        """Set the lower bound of the reaction.

        Infeasible combinations, such as a upper bound lower than the current
        lower bound will update the other bound.

        Parameters
        ----------
        value: float
            The new value for the lower bound.

        """
        if self.upper_bound < value:
            self.upper_bound = value

        setattr(self, "_lower_bound", value)

    @property
    def upper_bound(self):
        """Get the lower bound of the reaction."""
        return getattr(self, "_upper_bound")

    @upper_bound.setter
    @resettable
    def upper_bound(self, value):
        if self.lower_bound > value:
            self.lower_bound = value

        setattr(self, "_upper_bound", value)

    @property
    def kf_str(self):
        """Return the string representation of the forward rate constant."""
        return "kf_" + self.id

    @property
    def Keq_str(self):
        """Return the string representation of the equilibrium constant."""
        return "Keq_" + self.id

    @property
    def kr_str(self):
        """Return the string representation of the reverse rate constant."""
        return "kr_" + self.id

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

    def reverse_stoichiometry(self, inplace=False):
        """Reverse the stoichiometry of the reaction.

        Reversing the stoichiometry will turn the products into the reactants
        and the reactants into the products.

        Parameters
        ----------
        inplace: bool, optional
            If True, modify the reaction directly. Otherwise a new MassReaction
            Object is created, and modified.

        Returns
        -------
        new_reaction: MassReaction
            Returns the original MassReaction if inplace=True. Otherwise return
            a modified copy of the original MassReaction.

        Warnings
        --------
        Only the stoichiometry of the reaction is modified. The reaction
        parameters (e.g. rate and equilibrium constants) are not altered.

        """
        if inplace:
            new_reaction = self
        else:
            new_reaction = self.copy()

        for met, coeff in iteritems(new_reaction.metabolites):
            new_reaction._metabolites[met] = -1 * coeff

        return new_reaction

    def get_rate_law(self, rate_type=1, sympy_expr=True,
                     update_reaction=False):
        """Get the rate law for the reaction.

        Parameters
        ----------
        rate type: int {1, 2, 3}, optional
            The type of rate law to display. Must be 1, 2, or 3.
            Type 1 will utilize kf and Keq.
            Type 2 will utilize kf and kr.
            Type 3 will utilize kr and Keq.
        sympy_expr: bool, optional
            If True, output is a sympy expression. Otherwise the output is a
            string.
        update_reaction: bool, optional
            If True, update the MassReaction in addition to returning the rate
            law.

        Returns
        -------
        The rate law expression as a str or sympy expression (sympy.Basic).

        """
        return generate_rate_law(self, rate_type, sympy_expr, update_reaction)

    def get_mass_action_ratio(self):
        """Get the mass action ratio of the reaction as a sympy expression.

        Returns
        -------
        The mass action ratio as a sympy expression (sympy.Basic).

        """
        return generate_mass_action_ratio(self)

    def get_disequilibrium_ratio(self):
        """Get the disequilibrium ratio of the reaction as a sympy expression.

        Returns
        -------
        The disequilibrium ratio as a sympy expression (sympy.Basic).

        """
        return generate_disequilibrium_ratio(self)

    def remove_from_model(self, remove_orphans=False):
        """Remove the reaction from the MassModel.

        This removes all associations between a reaction, the associated
        MassModel, metabolites, and genes.

        The change is reverted upon exit when using the MassModel as a context.

        Parameters
        ----------
        remove_orphans: bool, optional
            Remove orphaned genes and metabolites from the MassModel as well.

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
        metabolite_id: str, MassMetabolite
            The MassMetabolite or the string identifier of the MassMetabolite
            whose coefficient is desired.

        """
        if isinstance(metabolite_id, MassMetabolite):
            return self._metabolites[metabolite_id]

        _id_to_mets = {m.id: m for m in self._metabolites}
        return self._metabolites[_id_to_mets[metabolite_id]]

    def get_coefficients(self, metabolite_ids):
        """Return the coefficients for a list of metabolites in the reaction.

        Parameters
        ----------
        metabolite_ids: iterable
            Iterable of the MassMetabolites or their string identifiers.

        """
        return map(self.get_coefficient, metabolite_ids)

    def get_compartments(self):
        """Return a list of compartments where the metabolites are located."""
        return list(self.compartments)

    def add_metabolites(self, metabolites_to_add, combine=True,
                        reversibly=True):
        """Add metabolites and their coefficients to the reaction.

        If the final coefficient for a metabolite is 0 then it is removed from
        the reaction.

        The change is reverted upon exit when using the MassModel as a context.

        Parameters
        ----------
        metabolites_to_add: dict
            A dictionary with MassMetabolite objects or metabolite identifiers
            as keys and stoichiometric coefficients as values. If keys are
            strings (id of a metabolite), the reaction must already be part of
            a MassModel and a MassMetabolite with the given id must already
            exist in the MassModel.
        combine: bool, optional
            If True, the metabolite coefficients are combined together.
            Otherwise the coefficients are replaced.
        reversibly: bool, optional
            Whether to add the change to the context to make the change
            reversible (primarily intended for internal use).

        Warnings
        --------
        A cobra Metabolite cannot be directly added to a MassReaction. Instead,
            the cobra Metabolite must first be converted to a MassMetabolite
            through the mass.util.conversion class.

        Notes
        -----
        A final coefficient of < 0 implies a reactant and a final
            coefficient of > 0 implies a product.

        See Also
        --------
        MassReaction.subtract_metabolites

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
        """Subtract metabolites and their coefficients from the reaction.

        This function will 'subtract' metabolites from a reaction by adding
        the given metabolites with -1*coeffcient. If the final coefficient for
        a metabolite is 0, the metabolite is removed from the reaction.

        The change is reverted upon exit when using the MassModel as a context.

        Parameters
        ----------
        metabolites_to_subtract: dict
            A dictionary with MassMetabolite objects or metabolite identifiers
            as keys and stoichiometric coefficients as values. If keys are
            strings (id of a metabolite), the reaction must already be part of
            a MassModel and a MassMetabolite with the given id must already
            exist in the MassModel.
        combine: bool, optional
            If True, the metabolite coefficients are combined together.
            Otherwise the coefficients are replaced.
        reversibly: bool, optional
            Whether to add the change to the context to make the change
            reversible (primarily intended for internal use).

        Warnings
        --------
        A cobra Metabolite cannot be directly added to a MassReaction. Instead,
            the cobra Metabolite must first be converted to a MassMetabolite
            through the mass.util.conversion class.

        Notes
        -----
        A final coefficient of < 0 implies a reactant and a final
            coefficient of > 0 implies a product.

        See Also
        --------
        MassReaction.add_metabolites

        """
        self.add_metabolites({k: -v
                              for k, v in iteritems(metabolites_to_subtract)},
                             combine=combine, reversibly=reversibly)

    def build_reaction_string(self, use_metabolite_names=False):
        """Generate a human readable string to represent the reaction.

        Parameters
        ----------
        use_metabolite_names: bool, optional
            If True, use the metabolite names instead of their identifiers.
            Default is false.

        Returns
        -------
        reaction_string: str
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
        """Compute tbhe mass and charge balances for the reaction.

        Returns a dictionary of {element: amount} for unbalanced elements,
        with the "charge" treated as an element in this dictionary.
        For a balanced reaction, an empty dictionary is returned.
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
        """Build a reaction from reaction equation reaction_str using parser.

        Takes a string representation of the reaction and uses the
        specifications supplied in the optional arguments to infer a set of
        metabolites, metabolite compartments, and stoichiometries for the
        reaction. It also infers the refversibility of the reaction from the
        reaction arrow.

        For example:
            'A + B <=> C' for reversible reactions, A & B are reactants.
            'A + B --> C' for irreversible reactions, A & B are reactants.
            'A + B <-- C' for irreversible reactions, A & B are products.

        The change is reverted upon exit when using the MassModel as a context.

        Parameters
        ----------
        reaction_str: str
            A string containing the reaction formula (equation).
        verbose: bool, optional
            Setting the verbosity of the function.
        fwd_arrow: re.compile, optional
            For forward irreversible reaction arrows.
        rev_arrow: re.compile, optional
            For backward irreversible reaction arrows.
        reversible_arrow: re.compile, optional
            For reversible reaction arrows.
        term_split: str, optional
            Dividing individual metabolite entries. Default is "+".

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
        """Associates a cobra.Gene object with the reaction.

        Parameters
        ----------
        cobra_gene: cobra.core.Gene.Gene
            Gene object to be assoicated with the reaction.

        Warnings
        --------
        This method is intended for internal use only.

        """
        self._genes.add(cobra_gene)
        cobra_gene._reaction.add(self)
        cobra_gene._model = self._model

    def _dissociate_gene(self, cobra_gene):
        """Dissociates a cobra.Gene object with the reaction.

        Parameters
        ----------
        cobra_gene: cobra.core.Gene.Gene
            Gene object to be assoicated with the reaction.

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

    def _repr_html_(self):
        """HTML representation of the overview for the MassReaction."""
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
        """Create a copy of the MassReaction."""
        return copy(super(MassReaction, self))

    def __deepcopy__(self, memo):
        """Create a deepcopy of the MassReaction."""
        return deepcopy(super(MassReaction, self), memo)

    def __setstate__(self, state):
        """Set state of MassReaction object upon unpickling.

        Probably not necessary to set _model as the MassModel that
        contains self sets the _model attribute for all metabolites and genes
        in the reaction.

        However, to increase performance speed, let the metabolite
        and gene know that they are employed in this reaction

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

        Similar to the method in cobra.core.reaction
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
        """
        new_reaction = self.copy()
        new_reaction -= other
        return new_reaction

    def __isub__(self, other):
        """Subtract two reactions.

        The stoichiometry will be the combined stoichiometry of the two
        reactions. Returns the same MassReaction object with its updates.
        """
        self.subtract_metabolites(other._metabolites, combine=True)
        return self

    def __mul__(self, coefficient):
        """Scale coefficients in a reaction by a given value.

        Returns a new MassReaction object
        E.g. A -> B becomes 2A -> 2B.
        """
        new_reaction = self.copy()
        new_reaction *= coefficient
        return new_reaction

    def __imul__(self, coefficient):
        """Scale coefficients in a reaction by a given value.

        Returns the same MassReaction object with its updates
        E.g. A -> B becomes 2A -> 2B.
        """
        self._metabolites = {met: value * coefficient
                             for met, value in iteritems(self._metabolites)}
        return self

    def __str__(self):
        """Create an id string with the stoichiometry."""
        return "{id}: {stoichiometry}".format(
            id=self.id, stoichiometry=self.build_reaction_string())
