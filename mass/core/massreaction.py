# -*- coding: utf-8 -*-

# Compatibility with Python 2.7
from __future__ import absolute_import

# Import necesary packages
import re
from copy import copy, deepcopy
from functools import partial
from operator import attrgetter
from warnings import warn
from six import iteritems, iterkeys, string_types


from sympy import sympify, S, var, Add, Mul, Pow, Integer

# from cobra
from cobra.core.object import Object
from cobra.core.metabolite import Metabolite
from cobra.core.gene import Gene, ast2str, parse_gpr, eval_gpr
from cobra.util.context import resettable, get_context

# from mass
from mass.core.massmetabolite import MassMetabolite

# Class begins
## precompiled regular expressions
### Matches and/or in a gene reaction rule
and_or_search = re.compile(r'\(| and| or|\+|\)', re.IGNORECASE)
uppercase_AND = re.compile(r'\bAND\b')
uppercase_OR = re.compile(r'\bOR\b')
gpr_clean = re.compile(' {2,}')
### This regular expression finds any single letter compartment enclosed in
### square brackets at the beginning of the string.
### For example [c] : foo --> bar
compartment_finder = re.compile("^\s*(\[[A-Za-z]\])\s*:*")
### Regular expressions to match the arrows for building reactions from strings
_reversible_arrow_finder = re.compile("<(-+|=+)>")
_forward_arrow_finder = re.compile("(-+|=+)>")
_reverse_arrow_finder = re.compile("<(-+|=+)")

# Class definition
class MassReaction(Object):
    """MassReaction is a class for holding kinetic information regarding a
    biochemical reaction in a mass.MassModel object

    Parameters
    ----------
    id : string
        The identifier to associate with this reaction
    name : string
        A human readable name for the reaction
    subsystem : string
        The subsystem where the reaction is meant to occur
    reversibility : bool
        The kinetic reversibility of the reaction

    Attributes
    ----------
    sym_kf : string
        String representation of the symbol for the forward rate constant.
    sym_kr : string
        String representation of the symbol for the reverse rate constant.
    sym_Keq : string
        String representation of the symbol for the equilibrium rate constant.
    """

    def __init__(self, id=None, name="", subsystem="", reversibility=True):
        """Initialize the MassReaction Object"""
        # Check inputs to ensure they are the correct types
        if not isinstance(name, string_types):
            raise TypeError("name must be a string type")
        elif not isinstance(subsystem, string_types):
            raise TypeError("subsystem must be a string type")
        elif not isinstance(reversibility, bool):
            raise TypeError("reversibility must be a boolean")
        else:
            pass

        Object.__init__(self, id, name)
        self._reversibility = reversibility
        self._subsystem = subsystem
        # The forward, reverse, and equilibrium constants as strings
        # for symbolic expressions
        self.sym_kf = ("kf_%s" % id)
        self.sym_kr = ("kr_%s" % id)
        self.sym_Keq = ("Keq_%s" % id)

        # The forward, reverse, and equilbrium constants for simulation
        # initialized as strings of their symbolic representations
        self._forward_rate_constant = self.sym_kf
        # No reverse rate constant for an irreversible reaction.
        # Therefore initialized to 0.
        if self._reversibility == True:
            self._reverse_rate_constant = self.sym_kr
        else:
            self._reverse_rate_constant = 0.
        self._equilibrium_constant = self.sym_Keq

        # The rate law equation for simulation and the symbolic representation
        self._rate_law = None
        self._rate_law_expr = None

        # A dictionary of metabolites and their stoichiometric
        # coefficients for this kinetic reaction
        self._metabolites = dict()

        # The compartments where partaking metabolites are located
        self._compartments = None

        # The massmodel that the reaction is associated with
        self._model = None

        # The genes associated with the kinetic reaction
        self._genes = set()
        self._gene_reaction_rule = ""

    # Properties
    @property
    def reversibility(self):
        """Returns the kinetic reversibility of the reaction"""
        return self._reversibility

    @reversibility.setter
    def reversibility(self, value):
        """Set the kinetic reversibility of the reaction. Will initialize to
        default values for the reverse rate constant when the kinetic
        reversibility is changed to True, and will set the reverse rate
        constant to 0 if changed to False

        Parameters
        ----------
        reversibility : bool
            True for kinetically reversible reaction, False for irreversible
        """
        if not isinstance(value, bool):
            raise TypeError("Must be a boolean True or False")

        self._reversibility = value
        if value == True:
            self._reverse_rate_constant = self.sym_kr
        else:
            self._reverse_rate_constant = 0.

    @property
    def subsystem(self):
        """Returns the subsystem associated with this reaction"""
        return self._subsystem

    @property
    def forward_rate_constant(self):
        """Returns the forward rate constant associated with this reaction"""
        return self._forward_rate_constant

    @forward_rate_constant.setter
    def forward_rate_constant(self, value):
        """Set the forward rate constant for this reaction"""
        self._forward_rate_constant = value

    @property
    def reverse_rate_constant(self):
        """Returns the reverse rate constant associated with this reaction
        If the reaction is not reversible, warns the user and return a 0."""
        if self._reversibility != True:
            warn("No reverse rate constant for irreversible reactions")
        return self._reverse_rate_constant

    @reverse_rate_constant.setter
    def reverse_rate_constant(self, value):
        """Set the reverse rate constant for this reaction.
        If the reaction is not reversible, warns the user"""
        if self._reversibility != True:
            warn("Cannot set reverse rate constant for irreversible reactions")
        else:
            self._reverse_rate_constant = value

    @property
    def equilibrium_constant(self):
        """Returns the equilibrium constant associated with this reaction"""
        return self._equilibrium_constant

    @equilibrium_constant.setter
    def equilibrium_constant(self, value):
        """Set the equilibrium constant for this reaction"""
        self._equilibrium_constant = value

    @property
    def rate_constants(self):
        """Returns a list containing the rate constants."""
        return [self._forward_rate_constant, self._reverse_rate_constant]

    @property
    def metabolites(self):
        """Returns the metabolites of the reaction as a read-only copy"""
        return self._metabolites.copy()

    @property
    def reactants(self):
        """Returns a list of the reactants for the reaction

        Identical to the method in cobra.core.reaction
        """
        return [m for m, c in iteritems(self._metabolites) if c < 0]

    @property
    def products(self):
        """Returns a list of the products for the reaction

        Identical to the method in cobra.core.reaction
        """
        return [m for m, c in iteritems(self._metabolites) if c >= 0]

    @property
    def stoichiometry(self):
        """Returns a list containing the stoichiometry of the reaction"""
        return [c for m, c in iteritems(self._metabolites)]

    @property
    def forward_rate(self):
        """Returns the forward rate law as a human readable string"""
        return self.generate_forward_rate(num_values=False)

    @property
    def forward_rate_expr(self):
        """Returns the forward rate law as a sympy expression"""
        return self.generate_forward_rate_expr(num_values=False)

    @property
    def reverse_rate(self):
        """Returns the reverse rate law as a human readable string.
        If the reaction is irreversible, warn the user and return float 0.
        """
        return self.generate_reverse_rate(num_values=False)

    @property
    def reverse_rate_expr(self):
        """Returns the reverse rate law as a sympy expression
        If the reaction is irreversible, warn the user and return symbolic 0
        """
        return self.generate_reverse_rate_expr(num_values=False)

    @property
    def rate_law(self):
        """Returns the rate law as a human readable string"""
        if self._rate_law == None:
            self._rate_law = self.generate_rate_law(num_values=False)
        return self._rate_law

    @property
    def rate_law_expr(self):
        """Returns the rate law as a sympy expression"""
        if self._rate_law_expr == None:
            self._rate_law_expr = self.generate_rate_law_expr(num_values=False)
        return self._rate_law_expr

    @property
    def model(self):
        """Returns the massmodel the reaction is associated with"""
        return self._model

    @property
    def reaction(self):
        """Return the reaction as a human readable string

        Similar to the method in cobra.core.reaction
        """
        return self.build_reaction_string()

    @reaction.setter
    def reaction(self, reaction_string):
        """Use a human readable string to set the reaction.
        The direction of the arrow is used to determine reversibility.

        Similar to the method in cobra.core.reaction
        Parameters
        ----------
        reaction_string : string
            String representation of the reaction. For example:
            'A + B <=> C' for reversible reactions
            'A + B --> C' for irreversible reactions where A & B are reactants
            'A + B <-- C' for irreversible reactions where A & B are products

        Warnings
        --------
        Care must be taken when doing setting a reaction in this manner to
            ensure the reaction's id matches those in the model.

        For forward and reversible arrows, reactants are searched on the
            left side of the arrow while products are searched on the right.
            For reverse arrows, reactants are searched on the left side of
            the arrow while products are searched on the right.
        """
        return self.build_reaction_from_string(reaction_string)

    @property
    def compartments(self):
        """Returns a list of compartments that metabolites are in

        Identical to the method in cobra.core.reaction
        """
        if self._compartments is None:
            self._compartments = {met.compartment for met in self._metabolites
                                  if met.compartment is not None}
        return self._compartments

    @property
    def boundary(self):
        """Whether or not this reaction is an exchange reaction
        Returns True if the reaction has either no products or reactants

        Identical to the method in cobra.core.reaction

        .. note:: These are reactions with a sink or source
        """
        return (len(self.metabolites) ==1 and
            not (self.reactants and self.products))

    @property
    def transport(self):
        """Whether or not this reaction is a transport reaction
        Returns True if the reaction has two different compartments involved.

        .. note:: These are reactions where at least one metabolite
                    crosses into a different compartment.
        """
        return (len(self.compartments) != 1)
    @property
    def genes(self):
        """Returns a frozenset of the genes associated with the reaction"""
        return frozenset(self._genes)

    @property
    def gene_reaction_rule(self):
        """Returns the gene reaction rule as a string"""
        return self._gene_reaction_rule

    @gene_reaction_rule.setter
    def gene_reaction_rule(self, new_rule):
        """Set the gene reaction rule using a string, associated the new genes
        and dissociated the old genes for the reaction

        Similar to the method in cobra.core.reaction

        Parameters
        ----------
        new_rule : string
            String representation of the new reaction rule
        """
        if get_context(self):
            warn("Context management not implemented for gene reaction rules")

        self._gene_reaction_rule = new_rule.strip()
        try:
            _, gene_names = parse_gpr(self._gene_reaction_rule)
        except (SyntaxError, TypeError) as e:
            if "AND" in new_rule or "OR" in new_rule:
                warn("uppercase AND/OR found in rule '%s' for %s" %
                    (new_rule, repr(self)))
                tmp_str = and_or_search.sub('', self._gene_reaction_rule)
                gene_names = set((gpr_clean.sub(' ', tmp_str).split(' ')))
            if '' in gene_names:
                gene_names.remove('')
            old_genes = self._genes
            if self._model is None:
                self._genes = {Gene(i) for i in gene_names}
            else:
                massmodel_genes = self._model.genes
                self._genes = set()
                for id in gene_names:
                    if massmodel_genes.has_id(id):
                        self._genes.add(massmodel_genes.get_by_id(id))
                    else:
                        new_gene = Gene(id)
                        # Must be new_gene._model due to inheritance
                        # of cobra gene class
                        new_gene._model = self._model
                        self._genes.add(new_gene)
                        massmodel_genes.append(new_gene)

            # Make the genes aware that it is involved in this reaction
            for g in self._genes:
                g._reaction.add(self)

            # Make the old genes aware that they are no
            # longer involved in this reaction
            for g in old_genes:
                if g not in self._genes: # if an old gene is not a new gene
                    try:
                        g._reaction.remove(self)
                    except:
                        warn("Could not remove old gene %s from reaction %s" %
                            (g.id, self.id))

    @property
    def gene_name_reaction_rule(self):
        """Display gene_reaction_rule with names

        Identical to the method in cobra.core.reaction

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
        """Check if all required enzymes for reaction are functional.

        Returns True if the gene-protein-reaction (GPR) rule is fulfilled for
            this reaction, or if reaction is not associated to a massmodel,
            otherwise returns False.

        Identical to the method in cobra.core.reaction
        """
        if self._model:
            tree, _ = parse_gpr(self.gene_reaction_rule)
            return eval_gpr(tree, {gene.id for gene in self.genes if
                                   not gene.functional})
        return True

    # Methods
    def generate_forward_rate(self, num_values=False):
        """Generates the forward rate law for the reaction and
        returns a human readable string. Returns None if the
        reaction does not have reactants

        Parameters
        ----------
        num_values : bool
            If True, the value of the rate constant is used.
            Otherwise use a symbol for the rate constant.
        """
        if len(self.reactants) == 0:
            warn("Cannot generate a forward rate when there are no"
                                " reactants for this reaction")
            return None

        if num_values != False:
            self._forward_rate = str(self._forward_rate_constant)
        else:
            self._forward_rate = self.sym_kf

        for metab in self.reactants:
            coeff = self.get_coefficient(metab.id)
            if abs(coeff) == 1:
                self._forward_rate += "*%s" % metab.id
            else:
                self._forward_rate += "*%s**%s" % \
                                    (metab.id, coeff)
        return self._forward_rate

    def generate_forward_rate_expr(self, num_values=False):
        """Generates the forward rate law for the reaction and
        returns a sympy expression. Returns None if the
        reaction does not have reactants

        Parameters
        ----------
        num_values : bool
            If True, the value of the rate constant is used.
            Otherwise use a symbol for the rate constant.
        """
        if len(self.reactants) == 0:
            warn("Cannot generate a forward rate when there are no"
                                " reactants for this reaction")
            return None

        if num_values != False:
            self._forward_rate_expr = sympify(self._forward_rate_constant)
        else:
            self._forward_rate_expr = sympify(self.sym_kf)

        for metab in self.reactants:
            coeff = self.get_coefficient(metab.id)
            if abs(coeff) == 1:
                self._forward_rate_expr = Mul(self._forward_rate_expr,
                                                var(metab.id))
            else:
                self._forward_rate_expr = Mul(self._forward_rate_expr,
                                            Pow(var(metab.id), coeff))
        return self._forward_rate_expr

    def generate_reverse_rate(self, num_values=False):
        """Generates the reverse rate law for the reaction and
        returns a human readable string. If the reaction is not reversible,
        returns a 0. Returns None if the reaction does not have products.

        Parameters
        ----------
        num_values : bool
            If True, the value of the rate constant is used.
            Otherwise use a symbol for the rate constant.
        """
        if len(self.products) == 0:
            warn("Cannot generate a reverse rate when there are no"
                                " products for this reaction")
            return None
        if self._reversibility != True:
            return 0.

        if num_values != False:
            self._reverse_rate = str(self._reverse_rate_constant)
        else:
            self._reverse_rate = self.sym_kr

        for metab in self.products:
            coeff = self.get_coefficient(metab.id)
            if abs(coeff) == 1:
                self._reverse_rate += "*%s" % metab.id
            else:
                self._reverse_rate += "*%s**%s" % \
                                    (metab.id, coeff)
        return self._reverse_rate

    def generate_reverse_rate_expr(self, num_values=False):
        """Generates the reverse rate law for the reaction and
        returns a human readable string. If the reaction is not reversible,
        returns a 0. Returns None if the reaction does not have products.

        Parameters
        ----------

        num_values : bool
            If True, the value of the rate constant is used.
            Otherwise use a symbol for the rate constant.
        """
        if len(self.products) == 0:
            warn("Cannot generate a reverse rate when there are no"
                                " products for this reaction")
            return None
        if self._reversibility != True:
            return S.Zero

        if num_values != False:
            self._reverse_rate_expr = sympify(self._reverse_rate_constant)
        else:
            self._reverse_rate_expr = sympify(self.sym_kr)

        for metab in self.products:
            coeff = self.get_coefficient(metab.id)
            if abs(coeff) == 1:
                self._reverse_rate_expr = Mul(self._reverse_rate_expr,
                                                var(metab.id))
            else:
                self._reverse_rate_expr = Mul(self._reverse_rate_expr,
                                            Pow(var(metab.id), coeff))
        return self._reverse_rate_expr

    def generate_rate_law(self, num_values=False):
        """Generates the rate law for the reaction and
        returns a human readable string. If products and reactants are not
        defined, generate a warning and return None.

        Parameters
        ----------
        num_values : bool
            If True, the value of the rate constant is used.
            Otherwise use a symbol for the rate constant.

        Warnings
        --------
        Using the generate_rate_law method will replace a custom rate law.
        """
        self._forward_rate = self.generate_forward_rate(num_values)
        self._reverse_rate = self.generate_reverse_rate(num_values)
        if self._forward_rate == None and self._reverse_rate == None:
            self._rate_law = None
        else:
            self._rate_law = ("%s - %s" % (self._forward_rate, self._reverse_rate))
        return self._rate_law

    def generate_rate_law_expr(self, num_values=False):
        """Generates the rate law for the reaction and
        returns a sympy expression

        Parameters
        ----------
        num_values : bool
            If True, the value of the rate constant is used.
            Otherwise use a symbol for the rate constant.

        Warnings
        --------
        Using the generate_rate_law_expr method will replace a
        custom rate law expression.
        """
        self._forward_rate_expr = self.generate_forward_rate_expr(num_values)
        self._reverse_rate_expr = self.generate_reverse_rate_expr(num_values)
        if self._forward_rate_expr == None and self._reverse_rate_expr == None:
            self._rate_law_expr = None
        else:
            self._rate_law_expr = Add(self._forward_rate_expr,
                                    Mul(Integer(-1), self._reverse_rate_expr))

        return self._rate_law_expr

    def reset_rate_law(self, num_values=False):
        """Reset the custom rate law and custom rate law expression
        to the automatically generated rate law and rate law expression

        Parameters
        ----------
        num_values : bool
            If True, the value of the rate constant is used.
            Otherwise use a symbol for the rate constant.
        """
        self._rate_law = self.generate_rate_law(num_values)
        self._rate_law_expr = self.generate_rate_law_expr(num_values)

    def remove_from_model(self, remove_orphans=False):
        """Removes the reaction from a massmodel.

        This removes all associations between a reaction, the associated
        massmodel, metabolites and genes.

        The change is reverted upon exit when using the massmodel as a context.

        Identical to the method in cobra.core.reaction

        Parameters
        ----------
        remove_orphans : bool
            Remove orphaned genes and metabolites from the massmodel as well
        """
        return self._model.remove_reactions([self],remove_orphans)

    def copy(self):
        """Copy a reaction.

        The rate constants, rate laws, referenced metabolites, and genes are
        also copied
        """
        # No references to massmodel when copying
        massmodel = self._model
        self._model = None
        for i in self._metabolites:
            i._model = None
        for i in self._genes:
            i._model = None

        # The reaction can be copied
        new_massreaction = deepcopy(self)
        # Restore the references
        self._model = massmodel
        for i in self._metabolites:
            i._model = massmodel
        for i in self._genes:
            i._model = massmodel

        return new_massreaction

    def get_coefficient(self, metabolite_id):
        """Return the stoichiometric coefficients of a metabolite
        in the reaction

        Similar to the method in cobra.core.reaction

        Parameters
        ----------
        metabolite_id : string or mass.MassMetabolite object.
        """
        if isinstance(metabolite_id, MassMetabolite):
            return self._metabolites[metabolite_id]

        _id_to_metabolites = {m.id: m for m in self._metabolites}
        return self._metabolites[_id_to_metabolites[metabolite_id]]

    def get_coefficients(self, metabolite_ids):
        """Return the stoichiometric coefficients for a list of
        metabolites in the reaction.

        Identical to the method in cobra.core.reaction

        Parameters
        ----------
        metabolite_ids : iterable
            Containing strings or mass.MassMetabolite objects.
        """
        return map(self.get_coefficient, metabolite_ids)

    def add_metabolites(self, metabolites_to_add, combine=True,
                        reversibly=True):
        """Add metabolites and stoichiometric coefficients to the reaction.
        If the final coefficient for a metabolite is 0 then it is removed
        from the reaction.

        The change is reverted upon exit when using the massmodel as a context.

        Similar to the method in cobra.core.reaction

        Parameters
        ----------
        metabolites_to_add : dict
            Dictionary with MassMetabolite objects or metabolite identifiers as
            keys and coefficients as values. If keys are strings (name of a
            metabolite) the reaction must already be part of a massmodel and a
            metabolite with the given name must exist in the massmodel.

        combine : bool
            Describes behavior a metabolite already exists in the reaction.
            True causes the coefficients to be added.
            False causes the coefficient to be replaced.

        reversibly : bool
            Whether to add the change to the context to make the change
            reversibly or not (primarily intended for internal use).

        Warnings
        --------
        To add a cobra Metabolite object to a MassReaction object, the
            cobra Metabolite must first be converted to a MassMetabolite
            through the from_cobra method in the mass.core.massmetabolite class
        """
        for metabolite, coefficient in iteritems(metabolites_to_add):
            if isinstance(metabolite, Metabolite):
                raise TypeError("Must be a MassMetabolite object and not a "
                            "cobra Metabolite object. Create a MassMetabolite "
                            "by using the from_cobra method on the Metabolite: "
                            "%s" % metabolite.id)

        old_coefficients = self.metabolites
        new_metabolites = []
        _id_to_metabolites = dict([(x.id, x) for x in self._metabolites])

        for metabolite, coefficient in iteritems(metabolites_to_add):
            met_id = str(metabolite)
            # If a metabolite already exists in the reaction then
            # just add them.
            if met_id in _id_to_metabolites:
                reaction_metabolite = _id_to_metabolites[met_id]
                if combine:
                    self._metabolites[reaction_metabolite] += coefficient
                else:
                    self._metabolites[reaction_metabolite] = coefficient
            else:
                # If the reaction is in a massmodel, ensure we aren't using
                # a duplicate metabolite.
                if self._model:
                    try:
                        metabolite = \
                            self._model.metabolites.get_by_id(met_id)
                    except KeyError as e:
                        if isinstance(metabolite, MassMetabolite):
                            new_metabolites.append(metabolite)
                        else:
                            # do we want to handle creation here?
                            raise e
                elif isinstance(metabolite, string_types):
                    # if we want to handle creation, this should be changed
                    raise ValueError("Reaction '%s' does not belong to a "
                                     "massmodel. Either add the reaction to a "
                                     "massmodel or use MassMetabolite objects "
                                     "instead of strings as keys."
                                     % self.id)
                self._metabolites[metabolite] = coefficient
                # make the metabolite aware that it is involved in this
                # reaction
                metabolite._reaction.add(self)

        for metabolite, the_coefficient in list(self._metabolites.items()):
            if the_coefficient == 0:
                # make the metabolite aware that it no longer participates
                # in this reaction
                metabolite._reaction.remove(self)
                self._metabolites.pop(metabolite)

        massmodel = self.model
        context = get_context(self)
        if context and reversibly:
            if combine:
                # Just subtract the metabolites that were added
                context(partial(
                        self.subtract_metabolites, metabolites_to_add,
                        combine=True, reversibly=False))
            else:
                # Reset the metabolites with add_metabolites
                mets_to_reset = {key: old_coefficients[
                                massmodel.metabolites.get_by_any(key)[0]]
                                for key in iterkeys(metabolites_to_add)}

    def subtract_metabolites(self, metabolites_to_subtract, combine=True,
                            reversibly=True):
        """This function will 'subtract' metabolites from a reaction, which
        means add the metabolites with -1*coefficient. If the final coefficient
        for a metabolite is 0 then the metabolite is removed from the reaction.

        The change is reverted upon exit when using the massmodel as a context.

        Similar to the method in cobra.core.reaction

        Parameters
        ----------
        metabolites_to_subtract : dict
            Dictionary with MassMetabolite objects or metabolite identifiers as
            keys and coefficients as values. If keys are strings (name of a
            metabolite) the reaction must already be part of a massmodel and a
            metabolite with the given name must exist in the massmodel.

        combine : bool
            Describes behavior a metabolite already exists in the reaction.
            True causes the coefficients to be added.
            False causes the coefficient to be replaced.

        reversibly : bool
            Whether to add the change to the context to make the change
            reversibly or not (primarily intended for internal use).

        .. note:: A final coefficient < 0 implies a reactant.

        Warnings
        --------
        To add a cobra Metabolite object to a MassReaction object, the
            cobra Metabolite must first be converted to a MassMetabolite through
            the from_cobra method in the mass.core.massmetabolite class
        """
        self.add_metabolites({
            k: -v for k, v in iteritems(metabolites_to_subtract)},
            combine=combine, reversibly=reversibly)

    def _set_id_with_model(self, value):
        """Set the id of the MassReaction object to the associated massmodel.

        Similar to the method in cobra.core.reaction
        """

        if value in self.massmodel.reactions:
            raise ValueError("The massmodel already contains a reaction with "
                                "the id:", value)
        self._id = value
        self.massmodel.reactions._generate_index()

    def _update_awareness(self):
        """Make sure all metabolites and genes that are associated with
        this reaction are aware of it.
        """
        for metab in self._metabolites:
            metab._reaction.add(self)
        for gene in self._genes:
            gene._reaction.add(self)

    def build_reaction_string(self, use_metabolite_names=False):
        """Generate a human readable reaction string

        Similar to the method in cobra.core.reaction
        """
        def format(number):
            return "" if number == 1 else str(number).rstrip(".") + " "

        id_type = "id"
        if use_metabolite_names:
            id_type = "name"
        reactant_bits = []
        product_bits = []
        for metab in sorted(self._metabolites, key=attrgetter("id")):
            coefficient = self._metabolites[metab]
            metab_name = str(getattr(metab, id_type))
            if coefficient >= 0:
                product_bits.append(format(coefficient) + metab_name)
            else:
                reactant_bits.append(format(abs(coefficient)) + metab_name)
        reaction_string = " + ".join(reactant_bits)
        if self._reversibility != True:
            reaction_string += " --> "
        else:
            reaction_string += " <=> "
        reaction_string += " + ".join(product_bits)
        return reaction_string


    def check_mass_balance(self):
        """Compute mass and charge balance for the reaction

        returns a dict of {element: amount} for unbalanced elements.
        "charge" is treated as an element in this dict
        This should be empty for balanced reactions.

        Identical to the method in cobra.core.reaction
        """
        reaction_element_dict = defaultdict(int)
        for metabolite, coefficient in iteritems(self._metabolites):
            if metabolite.charge is not None:
                reaction_element_dict["charge"] += \
                    coefficient * metabolite.charge
            if metabolite.elements is None:
                raise ValueError("No elements found in metabolite %s"
                                 % metabolite.id)
            for element, amount in iteritems(metabolite.elements):
                reaction_element_dict[element] += coefficient * amount
        # filter out 0 values
        return {k: v for k, v in iteritems(reaction_element_dict) if v != 0}

    def get_compartments(self):
        """Return a list of compartments the metabolites are in

        Identical to the method in cobra.core.reaction
        """
        return list(self.compartments)

    def build_reaction_from_string(self, reaction_string, verbose=True,
                                   fwd_arrow=None, rev_arrow=None,
                                   reversible_arrow=None, term_split="+"):
        """Builds reaction from reaction equation reaction_string using parser

        Takes a string and using the specifications supplied in the optional
        arguments infers a set of metabolites, metabolite compartments and
        stoichiometries for the reaction.  It also infers the reversibility
        of the reaction from the reaction arrow.

        Changes to the associated massmodel are reverted upon exit when using
        the massmodel as a context.

        Similar to the method in cobra.core.reaction.

        Parameters
        ----------
        reaction_str : string
            a string containing a reaction formula (equation)
        verbose: bool
            setting verbosity of function
        fwd_arrow : re.compile
            for forward irreversible reaction arrows
        rev_arrow : re.compile
            for backward irreversible reaction arrows
        reversible_arrow : re.compile
            for reversible reaction arrows
        term_split : string
            dividing individual metabolite entries

        """
        # Set the arrows
        forward_arrow_finder = _forward_arrow_finder if fwd_arrow is None \
            else re.compile(re.escape(fwd_arrow))
        reverse_arrow_finder = _reverse_arrow_finder if rev_arrow is None \
            else re.compile(re.escape(rev_arrow))
        reversible_arrow_finder = _reverse_arrow_finder \
            if reversible_arrow is None \
            else re.compile(re.escape(reversible_arrow))
        if self._model is None:
            warn("No massmodel found")
            massmodel = None
        else:
            massmodel = self._model
        found_compartments = compartment_finder.findall(reaction_string)
        if len(found_compartments) == 1:
            compartment = found_compartments[0]
        else:
            compartment = ""

        #Reversible reaction
        arrow_match = reversible_arrow_finder.search(reaction_string)
        if arrow_match is not None:
            self._reversibility = True
            # Reactants left of the arrow, products on the right
            reactant_str = reaction_string[:arrow_match.start()].strip()
            product_str = reaction_string[arrow_match.end():].strip()
        else: # Irreversible reaction
            # Try forward reaction
            arrow_match = forward_arrow_finder.search(reaction_string)
            if arrow_match is not None:
                self._reversibility = False
                # Reactants left of the arrow, products on the right
                reactant_str = reaction_string[:arrow_match.start()].strip()
                product_str = reaction_string[arrow_match.end():].strip()
            else: # Try the reverse arrow
                arrow_match = reverse_arrow_finder.search(reaction_string)
                if arrow_match is None:
                    raise ValueError("No suitable arrow found in '%s'" %
                                        reaction_str)
                else:
                    self._reversibility = False
                    # Reactants right of the arrow, products on the left
                    reactant_str = reaction_string[arrow_match.end():].strip()
                    product_str = reaction_string[:arrow_match.start()].strip()

        self.subtract_metabolites(self.metabolites, combine=True)

        for substr, factor in ((reactant_str, -1), (product_str, 1)):
            if len(substr) == 0:
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
                    met = massmodel.metabolites.get_by_id(met_id)
                except KeyError:
                    if verbose:
                        print("Unknown metabolite '%s' created" % met_id)
                    met = MassMetabolite(met_id)
                self.add_metabolites({met: num})

    def _associate_gene(self, cobra_gene):
        """Associates a cobra.Gene object with a mass.MassReaction.

        Identical to the method in cobra.core.reaction

        Parameters
        ----------
        cobra_gene : cobra.core.Gene.Gene
        """
        self._genes.add(cobra_gene)
        cobra_gene._reaction.add(self)
        cobra_gene._model = self._model

    def _dissociate_gene(self, cobra_gene):
        """Dissociates a cobra.Gene object with a mass.MassReaction.

        Identical to the method in cobra.core.reaction

        Parameters
        ----------
        cobra_gene : cobra.core.Gene.Gene
        """
        self._genes.discard(cobra_gene)
        cobra_gene._reaction.discard(self)

    def knock_out(self):
        """Knockout reaction by setting its rate_constants to 0.

        Similar to the method in cobra.core.reaction"""
        self.forward_rate_constant = 0.
        if self._reversibility == True:
            self.reverse_rate_constant = 0.

    def set_custom_rate_law(self, custom_rate_law):
        """Use a string to set a custom rate law for the reaction

        Parameters
        ----------
        custom_rate_law : string
            String representation of the custom rate law.

        Warnings
        --------
        The custom rate law must be set as a string and will replace any
            previously generated rate law for this reaction. Using the
            generate_rate_law method or the reset_rate_law method will
            replace the custom rate law.
        """
        print("FIXME: Implement")
        return

    def generate_custom_rate_law_expr(self, custom_rate_law):
        """Generate the custom rate law expression from the custom rate law
        string

        Parameters
        ----------
        custom_rate_law : string
            String representation of the custom rate law.

        Warnings
        --------
        The custom rate law must be set as a string and will replace any
            previously generated rate law for this reaction. Using the
            generate_rate_law_expr method or the reset_rate_law method will
            replace the custom rate law expression.
        """
        print("FIXME: Implement")
        return

    def _repr_html_(self):
        return """
            <table>
                <tr>
                    <td><strong>Reaction identifier</strong></td>
                    <td>{id}</td>
                </tr><tr>
                    <td><strong>Name</strong></td>
                    <td>{name}</td>
                </tr><tr>
                    <td><strong>Subsystem</strong></td>
                    <td>{subsystem}</td>
                </tr><tr>
                    <td><strong>Memory address</strong></td>
                    <td>{address}</td>
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
                subsystem=self._subsystem, address='0x0%x' % id(self),
                stoich_id=self.build_reaction_string(),
                stoich_name=self.build_reaction_string(True),
                gpr=self.gene_reaction_rule,
                reversibility=self._reversibility)

    # Shorthands
    @property
    def kf(self):
        """Shorthand for getting the forward rate constant"""
        return self._forward_rate_constant

    @kf.setter
    def kf(self, value):
        """Shorthand for setting the forward rate constant"""
        self.forward_rate_constant = value

    @property
    def kr(self):
        """Shorthand for getting the reverse rate constant"""
        return self.reverse_rate_constant

    @kr.setter
    def kr(self, value):
        """Shorthand for setting the reverse rate constant"""
        self.reverse_rate_constant = value

    @property
    def Keq(self):
        """Shorthand for getting the equilibrium constant"""
        return self._equilibrium_constant

    @Keq.setter
    def Keq(self, value):
        """Shorthand for setting the forward rate constant"""
        self.equilibrium_constant = value

    @property
    def S(self):
        """Shorthand for the reaction stoichiometry"""
        return [c for m, c in iteritems(self._metabolites)]

    # Compatibility functions
    def to_cobra_reaction(self, cobra_id=None,
                            lower_bound=None,
                            upper_bound=None):
        """To convert a MassReaction object into a cobra Reaction object

        If the lower and/or upper bounds are not specified,the reversibility
        will be used to determine the bounds for initializing the reaction.

        For reversible MassReaction objects:
            lower_bound =-1000 and/or upper_bound=1000
        For irreversible MassReaction objects:
            lower_bound =0 and upper_bound=1000

        Warnings
		--------
        All other fields in a cobra Reaction will initialize to their defaults
        """
        try:
            from cobra.core.reaction import Reaction
        except:
            raise ImportError("Failed to import the Reaction Object from "
                    "cobra.core.reaction. Ensure cobra is installed properly")

        if cobra_id == None:
            cobra_id = self._id + "_cobra"

        if lower_bound == None:
            if self._reversibility == True:
                lb = -1000
            else:
                lb = 0.

        if upper_bound == None:
            ub = 1000

        cobra_rxn = Reaction(id=cobra_id, name=self.name,
                            subsystem=self._subsystem, lower_bound=lb,
                            upper_bound=ub, objective_coefficient=0.)

        print("FIXME: Add the current metabolites, model (if any) "
                "and genes to the new object")
        return cobra_rxn

    def from_cobra_reaction(self, CobraReaction=None, mass_id=None,
                            kinetic_reversibility=None):
        """To convert a cobra Reaction object into a MassReaction object

        If kinetic_reversibility is not specifiied, will try to infer
        reversibility from the upper and lower bounds of the cobra object.

        Warnings
		--------
        All other fields in a MassReaction will initialize to their defaults

        A Reaction Object from cobra.core.reaction must be imported into
            the enviroment in order for this method to work properly.
        """
        if CobraReaction == None:
            warn("No cobra Reaction Object was given")
            return None

        if not isinstance(CobraReaction, Reaction):
            raise TypeError("Reaction must be a cobra Reaction Object")

        if mass_id == None:
            mass_id = self._id + "_cobra"

        if kinetic_reversibility == None:
            kinetic_reversibility = CobraReaction.reversibility

        mass_rxn = MassReaction(id=mass_id, name=self.name,
                            subsystem=self._subsystem,
                            reversibility=kinetic_reversibility)

        print("FIXME: Add the current metabolites, model (if any) "
                "and genes to the new object")
        return mass_rxn

    # Module Dunders
    # All dunders are similar or identical to cobra.core.reaction dunders
    def __copy__(self):
        """Create a copy of the mass reaction

        Similar to the method in cobra.core.reaction
        """
        massreaction_copy = copy(super(MassReaction, self))
        massreaction_copy.reset_rate_law()
        return massreaction_copy

    def __deepcopy__(self, memo):
        """Create a deepcopy of the mass reaction

        Similar to the method in cobra.core.reaction
        """
        massreaction_deepcopy = deepcopy(super(MassReaction, self), memo)
        massreaction_deepcopy.reset_rate_law()
        return massreaction_deepcopy

    def __setstate__(self, state):
        """Probably not necessary to set _model as the mass.MassModel that
        contains self sets the _model attribute for all metabolites and
        genes in the reaction.

        However, to increase performance speed we do want to let the metabolite
        and gene know that they are employed in this reaction

        Similar to the method in cobra.core.reaction
        """
        self.__dict__.update(state)
        for x in state["_metabolites"]:
            setattr(x, "_model", self._model)
            x._reaction.add(self)
        for x in state["_genes"]:
            setattr(x, "_model", self._model)
            x._reaction.add(self)

    def __add__(self, other):
        """Add two mass reactions

        The stoichiometry will be the combined stoichiometry of the two
        reactions, and the gene reaction rule will be both rules combined by an
        and. All other attributes (i.e. rate constants) will match those of
        the first reaction, Return a new MassReaction object.

        Similar to the method in cobra.core.reaction
        """
        new_massreaction = self.copy()
        new_massreaction += other
        return new_massreaction

    def __iadd__(self, other):
        """Add two mass reactions

        The stoichiometry will be the combined stoichiometry of the two
        reactions, and the gene reaction rule will be both rules combined by an
        and. All other attributes (i.e. rate constants) will match those of the
        first reaction. Return the same MassReaction object with its updates

        Identical to the method in cobra.core.reaction
        """
        self.add_metabolites(other._metabolites, combine=True)
        gpr1 = self.gene_reaction_rule.strip()
        gpr2 = other.gene_reaction_rule.strip()
        if gpr1 != "" and gpr2 != "":
            self.gene_reaction_rule = ("(%s) and (%s)" % \
                                        (self.gene_reaction_rule,
                                        other.gene_reaction_rule))
        elif gpr1 != "" and gpr2 == "":
            self.gene_reaction_rule = gpr1
        elif gpr1 == "" and gpr2 != "":
            self.gene_reaction_rule = gpr2
        return self

    def __sub__(self, other):
        """Subtract two mass reactions

        The stoichiometry will be the combined stoichiometry of the two
        reactions, and the gene reaction rule will be both rules combined by an
        'and'. All other attributes (i.e. rate constants) will match those of
        the first reaction. Return a new MassReaction object.

        Similar to the method in cobra.core.reaction
        """
        new_massreaction = self.copy()
        new_massreaction -= other
        return new_massreaction

    def __isub__(self, other):
        """Subtract two mass reactions

        The stoichiometry will be the combined stoichiometry of the two
        reactions.

        Similar to the method in cobra.core.reaction
        """
        self.subtract_metabolites(other._metabolites, combine = True)
        return self


    def __mul__(self, coefficient):
        """Scale coefficients in a reaction by a given value.
        Return a new MassReaction object

        E.g. A -> B becomes 2A -> 2B.
        """
        new_massreaction = self.copy()
        new_massreaction *= coefficient
        return new_massreaction

    def __imul__(self, coefficient):
        """Scale coefficients in a reaction by a given value.
        Return the same MassReaction object with its updates

        E.g. A -> B becomes 2A -> 2B.

        Similar to the method in cobra.core.reaction
        """
        self._metabolites = {met: value * coefficient
                            for met, value in iteritems(self._metabolites)}
        return self

    def __str__(self):
        """Create an id string with the stoichiometry

        Identical to the method in cobra.core.reaction
        """
        return "{id}: {stoichiometry}".format(
            id=self.id, stoichiometry=self.build_reaction_string())
