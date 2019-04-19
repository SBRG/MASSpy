# -*- coding: utf-8 -*-
"""TODO Module Docstrings."""
from __future__ import absolute_import

import re
import warnings
from collections import defaultdict
from copy import deepcopy
from functools import partial

import numpy as np

from six import integer_types, iteritems, iterkeys, itervalues, string_types

import sympy as sym

from cobra.core.dictlist import DictList
from cobra.util.context import get_context

from mass.core.massmetabolite import MassMetabolite
from mass.core.massmodel import MassModel
from mass.core.massreaction import MassReaction
from mass.enzymes.enzyme_dict import EnzymeDict
from mass.enzymes.enzyme_form import EnzymeForm
from mass.util.expressions import _mk_met_func
from mass.util.util import _mk_new_dictlist, ensure_iterable, strip_time

_AUTOMATIC_RE = re.compile("^Automatic$")
_UNDEFINED_RE = re.compile("^Undefined$")
_EQUATION_RE = re.compile("^Equation$")


class EnzymeModel(MassModel):
    """Class representation of an EnzymeModel.

    Parameters
    ----------
    id_or_model: str, mass.MassModel
        Either an identifier to associate with the EnzymeModel given as a 
        string, or an existing EnzymeModel object. If an existing EnzymeModel
        is provided, a new EnzymeModel is instantiated with the same
        properties as the original EnzymeModel.
    name: str, optional
        A human readable name for the EnzymeModel.
    subsystem: str, optional
        The subsystem where the enzyme catalyzed net reaction that is
        represented by this EnzymeModel is meant to occur.
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
    enzyme_forms: cobra.DictList
        A cobra.DictList where the keys are the enzyme form identifiers and the
        values are the associated EnzymeForm objects. 
    enzyme_reactions: cobra.DictList
        A cobra.DictList where the keys are the enzyme reaction identifiers and
        the values are the associated MassReactions objects. A reaction is
        considered an enzyme reaction if it involves an EnzymeForm object.
    categorized_ligands: dict
        A dict of user-categorized ligands where keys are categories and
        values are DictLists of corresponding MassMetabolites.
    categorized_enzyme_forms: dict
        A dict of user-categorized enzyme forms where keys are categories and
        values are DictLists of the corresponding EnzymeForms.
    categorized_enzyme_reactions: dict
        A dict of user-categorized reactions involving EnzymeForms where keys
        are categories and values are DictLists of the corresponding 
        MassReactions.
    enzyme_concentration_total_equation: sympy.Basic
        A sympy expression representing the net reaction rate equation for the
        enzyme represented by the EnzymeModel.
    enzyme_net_flux_equation: sympy.Basic
        A sympy expression representing the net reaction rate equation for the
        enzyme represented by the EnzymeModel.
    enzyme_concentration_total: dict
        A dict containing the total enzyme concentration symbol as a sympy
        symbol and the total enzyme concentration value as a float.
    enzyme_net_flux: dict
        A dict containing the enzyme net flux symbol as a sympy
        symbol and the enzyme net flux value as a float.

    """

    def __init__(self, id_or_model=None, name=None, subsystem="",
                 matrix_type="dense", dtype=np.float64):
        """Initialize the EnzymeModel Object."""
        MassModel.__init__(self, id_or_model, name, matrix_type, dtype)
        if not isinstance(subsystem, string_types):
            raise TypeError("subsystem must be a str")
        self.subsystem = subsystem
        # Initialize DictLists for enzyme ligands, forms, and reaction objects
        self.ligands = DictList()
        self.enzyme_forms = DictList()
        self.enzyme_reactions = DictList()

        # Initialize a dict of DictLists for storing categorized objects
        self._categorized_ligands = {"Undefined": DictList()}
        self._categorized_enzyme_forms = {"Undefined": DictList()}
        self._categorized_enzyme_reactions = {"Undefined": DictList()}

        # Initialize EnzymeModel attributes
        self._enzyme_concentration_total = None
        self._enzyme_net_flux = None
        self.enzyme_net_flux_equation = None

    @property
    def enzyme_total_symbol(self):
        """Return the sympy symbol for the total enzyme concentration."""
        if self.id is not None:
            return sym.Symbol(self.id + "_Total")

    @property
    def enzyme_flux_symbol(self):
        """Return the sympy symbol for the net flux through the enzyme."""
        if self.id is not None:
            return sym.Symbol("v_" + self.id)

    @property
    def enzyme_concentration_total(self):
        """Return the total concentration value."""
        return getattr(self, "_enzyme_concentration_total")

    @enzyme_concentration_total.setter
    def enzyme_concentration_total(self, value):
        """Set the expected total enzyme concentration.

        Parameters
        ----------
        value: float
            A non-negative number for the total enzyme concentration.

        Warnings
        --------
        Concentrations cannot be negative.

        """
        if not isinstance(value, (integer_types, float)) and \
           value is not None:
            raise TypeError("Must be an int or float")
        elif value is None:
            pass
        elif value < 0.:
            raise ValueError("Must be a non-negative number")
        setattr(self, "_enzyme_concentration_total", value)

    @property
    def enzyme_net_flux(self):
        """Return a dict with the net flux symbol and value."""
        return getattr(self, "_enzyme_net_flux")

    @enzyme_net_flux.setter
    def enzyme_net_flux(self, value):
        """Set the expected net flux through the enzyme.

        Parameters
        ----------
        value: float
            The numerical value of the net flux through the enzyme.

        """
        if not isinstance(value, (integer_types, float)) and \
           value is not None:
            raise TypeError("Must be an int or float")
        else:
            setattr(self, "_enzyme_net_flux", value)

    @property
    def enzyme_concentration_total_equation(self):
        """Return the total concentration equation as a sympy.Equation."""
        if not self.enzyme_forms:
            return None
        return sym.Eq(
            self.enzyme_total_symbol, self.sum_enzyme_form_concentrations(
                self.enzyme_forms, use_values=False))

    @property
    def enzyme_net_flux_equation(self):
        """Return the net rate equation of the enzyme as a sympy.Equation."""
        return getattr(self, "_enzyme_net_flux_equation", None)

    @enzyme_net_flux_equation.setter
    def enzyme_net_flux_equation(self, value):
        """Set the net rate equation of the enzyme to a new sympy expression.

        The left side of the rate equation will always be the flux symbol of
        the enzyme, accessible self.enzyme_flux_symbol attribute.

        Parameters
        ----------
        value: sympy.Basic, str
            A sympy expression representing the right hand side of the rate
            equation, or a string to be turned into a symbolic expression via
            sympy.sympify.

        """
        if value is not None:
            if not isinstance(value, (sym.Basic, string_types)):
                raise TypeError("value must be a sympy expression.")
            elif isinstance(value, string_types):
                value = sym.sympify(value)
            elif value.lhs == self.enzyme_flux_symbol:
                value = value.rhs
            elif value.rhs == self.enzyme_flux_symbol:
                value = value.lhs
            else:
                pass
            value = sym.Eq(self.enzyme_flux_symbol, value)
        setattr(self, "_enzyme_net_flux_equation", value)

    @property
    def categorized_ligands(self):
        """Return a dict of ligands in user-defined categories."""    
        return self._remove_empty_categories("_categorized_ligands")

    @categorized_ligands.setter
    def categorized_ligands(self, value):
        """Set categories(s) for ligands using a dict.

        Parameters
        ----------
        value: dict
            A dict where keys are the strings representing the categories of
            the ligands, and values are lists of the corresponding ligands. If
            a categories already exists, its current contents will be replaced.
            An empty dict will cause a reset, placing all metabolites into an 
            "Undefined" category.

        Notes
        -----
        A metabolite must already exist in the EnzymeModel in order to set its
        category. Categories with empty lists are removed.

        See Also
        --------
        add_metabolites: 
            Method to add metabolites to the model.
        categorize_ligands: 
            Method to categorize "Undefined" ligands.

        """
        if not isinstance(value, dict):
            raise TypeError("value must be a dict")

        self._set_category_attribute_dict(value, "_categorized_ligands")

    @property
    def categorized_enzyme_forms(self):
        """Return a dict of enzyme forms in user-defined categories."""    
        return self._remove_empty_categories("_categorized_enzyme_forms")

    @categorized_enzyme_forms.setter
    def categorized_enzyme_forms(self, value):
        """Set categories(s) for enzyme forms using a dict.

        Parameters
        ----------
        value: dict
            A dict where keys are the strings representing the categories of
            the enzyme forms, and values are lists of the corresponding enzyme
            form. If a categories already exists, its current contents will be 
            replaced. An empty dict will cause a reset, placing all enzyme 
            forms into an "Undefined" category.

        Notes
        -----
        An enzyme form must already exist in the EnzymeModel in order to set
        its category. Categories with empty lists are removed.

        See Also
        --------
        add_metabolites: 
            Method to add metabolites or EnzymeForms to the model.
        categorize_enzyme_forms: 
            Method to categorize "Undefined" enzyme forms.

        """
        if not isinstance(value, dict):
            raise TypeError("value must be a dict")

        self._set_category_attribute_dict(value, "_categorized_enzyme_forms")

    @property
    def categorized_enzyme_reactions(self):
        """Return the enzyme reactions in user-defined categories."""    
        return self._remove_empty_categories("_categorized_enzyme_reactions")

    @categorized_enzyme_reactions.setter
    def categorized_enzyme_reactions(self, value):
        """Set categories(s) for enzyme binding reactions using a dict.

        Parameters
        ----------
        value: dict
            A dict where keys are the strings representing the categories
            of the binding reactions, and values are lists of the corresponding
            reactions. If a category already exists, its current contents
            will be replaced. An empty dict will cause a reset and place all
            reactions in an "Undefined" category.

        Notes
        -----
        A reaction must already exist in the EnzymeModel in order to set its
        category. Categories with empty lists are removed.

        See Also
        --------
        add_reactions: 
            Method to add reactions to the model.
        categorize_enzyme_reactions: 
            Method to categorize "Undefined" binding reactions.

        """
        if not isinstance(value, dict):
            raise TypeError("value must be a dict")

        self._set_category_attribute_dict(
            value, "_categorized_enzyme_reactions")

    def add_categorized_ligands(self, category, ligands):
        """Add a list of ligands to a new or existing category.

        If a category already exists, the ligands will be added to the
        existing list. If a ligand is categorized as "Undefined", 
        it will be removed from all other existing categories.

        Parameters
        ----------
        category: str
            A string representing the category for the list of ligands
        ligands: iterable of mass.MassMetabolites
            An iterable containing the MassMetabolite objects or their
            identifiers to be categorized.

        """
        self._add_categorized_items("_categorized_ligands", category, ligands)

    def add_categorized_enzyme_forms(self, category, enzyme_forms):
        """Add a list of enzyme forms to a new or existing category.

        If a category already exists, the enzyme forms will be added to the
        existing list. If an enzyme form is categorized as "Undefined", 
        it will be removed from all other existing categories.

        Parameters
        ----------
        category: str
            A string representing the category for the list of enzyme forms
        enzyme_forms: iterable of mass.EnzymeForms
            An iterable containing the EnzymeForm objects or their
            identifiers to be categorized.

        """
        self._add_categorized_items(
            "_categorized_enzyme_forms", category, enzyme_forms)

    def add_categorized_enzyme_reactions(self, category, enzyme_reactions):
        """Add a list of enzyme reactions to a new or existing category.

        If a category already exists, the reactions will be added to the
        existing list. If a reaction is categorized as "Undefined", 
        it will be removed from all other existing categories.

        Parameters
        ----------
        category: str
            A string representing the category for the list of reactions
        enzyme_reactions: iterable of mass.MassReactions
            An iterable containing the MassReactions objects or their
            identifiers to be categorized.

        """
        self._add_categorized_items(
            "_categorized_enzyme_reactions", category, enzyme_reactions)

    def make_enzyme_form(self, id=None, name="automatic", 
                         categories="Undefined", bound_catalytic=None, 
                         bound_effectors=None, compartment=None):
        """Make an EnzymeForm object to add to the EnzymeModel.

        Parameters
        ----------
        id: str
            The identifier associated with the EnzymeForm.
        name: str, optional
            A human readable name for the enzyme form. If name is set to match 
            "Automatic", a name will be generated based on the EnzymeForm and 
            its bound ligands.
        categories: str, iterable of str
            A string representing the category, or an iterable of strings 
            containing several categories for the EnzymeForm.
        bound_catalytic: dict, optional
            A dict representing the metabolites bound to the enzyme's active 
            site(s), with MassMetabolites or their identifiers as keys and the
            number bound as values.
        bound_effectors: dict, optional
            A dict representing the metabolites bound to the enzyme's 
            regulatory site(s), with MassMetabolites or their identifiers as
            keys and the number bound as values.
        compartment: str, optional
            The compartment where the EnzymeForm is located.

        Returns
        -------
        enzyme_form: EnzymeForm
            The newly created EnzymeForm object. 

        Notes
        -----
        If a metabolite does not exist in the EnzymeModel, it will be added to
            the module in addition to the EnzymeForm

        See Also
        --------
        EnzymeForm.generate_enzyme_form_name: 
            Automatic generation of name attribute for EnzymeForm

        """
        # Ensure metabolites for EnzymeForms exist in the EnzymeModel.
        for bound_dict in [bound_catalytic, bound_effectors]:
            if bound_dict is None:
                bound_dict = {}
                continue
            try:
                bound_dict = {self.metabolites.get_by_id(str(met)): num 
                              for met, num in iteritems(bound_dict)}
            except KeyError:
                # Add the metabolites into the module that don't already exist
                self.add_metabolites([met for met in iterkeys(bound_dict) 
                                      if met not in self.metabolites])
                bound_dict = {self.metabolites.get_by_id(str(met)): num 
                              for met, num in iteritems(bound_dict)}
        # Make EnzymeForm object
        enzyme_form = EnzymeForm(
            id=id, name=name, enzyme_id=self.id, enzyme_name=self.name,
            bound_catalytic=bound_catalytic, bound_effectors=bound_effectors, 
            compartment=compartment)
        # Generate a name for the EnzymeForm if the name is set to "automatic"
        if _AUTOMATIC_RE.match(name):
            enzyme_form.generate_enzyme_form_name(
                use_enzyme_name=True, update_enzyme=True)

        # Add the enzyme form to the module and place in respective categories
        self.add_metabolites(enzyme_form)
        categories = ensure_iterable(categories)
        for category in categories:
            self.add_categorized_enzyme_forms(category, enzyme_form)

        return enzyme_form

    def make_enzyme_reaction(self, id=None, name="", subsystem=None, 
                             reversible=True, categories="Undefined", 
                             metabolites_to_add=None):
        """Make an MassReaction object to add to the EnzymeModel.

        Parameters
        ----------
        id: str
            The identifier associated with the MassReaction.
        name: str, optional
            A human readable name for the enzyme reaction. If name is set to 
            match "Automatic", a name will be generated based on the
            EnzymeForms and their bound ligands.
        subsystem: str, optional
            The subsystem where the reaction is meant to occur.
        reversible: bool, optional
            The kinetic reversibility of the reaction. Irreversible reactions 
            havean equilibrium constant of infinity and a reverse rate constant 
            of 0. If not provided, the reaction is assumed to be reversible.
        categories: str, iterable of str
            A string representing the category, or an iterable of strings 
            containing several categories for the reaction.
        metabolites_to_add: dict
            A dictionary with MassMetabolite and EnzymeForm objects or their
            identifiers as keys and stoichiometric coefficients as values. If 
            keys are string identifiers of the objects, the MassMetabolite
            and EnzymeForm objects must already be a part of a MassModel.

        Returns
        -------
        enzyme_reaction: MassReaction
            The newly created MassReaction object. 

        Notes
        -----
        A final coefficient of < 0 implies a reactant and a final coefficient
            of > 0 implies a product.

        See Also
        --------
        generate_enzyme_reaction_name: 
            Automatic generation of name attribute for enzyme reactions.

        """
        # Make MassReaction object
        enzyme_reaction = MassReaction(id, name, subsystem, reversible)

        # Add metabolites
        if metabolites_to_add:
            try:
                metabolites_to_add = dict(
                    (met, c) if isinstance(met, MassMetabolite)
                    else (self.metabolites.get_by_id(met), c)
                    for met, c in iteritems(metabolites_to_add))
            except KeyError as e:
                raise KeyError(str(e) + "not found in model metabolites.")
            enzyme_reaction.add_metabolites(metabolites_to_add)

        # Add reaction to EnzymeModel
        self.add_reactions(enzyme_reaction)
        if enzyme_reaction not in self.enzyme_reactions:
            self.enzyme_reactions.add(enzyme_reaction)

        # Add categories for reaction
        categories = ensure_iterable(categories)
        for category in categories:
            self.add_categorized_enzyme_reactions(category, enzyme_reaction)

        # Set enzyme name if set to Automatic
        if _AUTOMATIC_RE.match(name):
            name = self.generate_enzyme_reaction_name(enzyme_reaction, 
                                                      update_enzyme=True)
        return enzyme_reaction

    def generate_enzyme_reaction_name(self, enzyme_reaction, 
                                      update_enzyme=False):
        """Generate a name for the enzymatic reaction based on bound ligands.

        The name is generated based on, bound_catalytic and bound_effector 
        attributes of the EnzymeForms involved in the enzyme reactions.

        Parameters
        ----------
        update_enzyme: bool, optional
            If True, update the name attribute of the enzyme form in
            addition to returning the automatically generated name of the 
            enzyme form as a str. Default is False.

        Returns
        -------
        name: A str representing the name of the enzyme reaction.

        """
        name = ""
        items = defaultdict(list)
        for met in enzyme_reaction.metabolites:
            key = "Enz" if isinstance(met, EnzymeForm) else "Lig"
            key += " React" if met in enzyme_reaction.reactants else " Prod"
            items[key].append(met)

        for attr in ["bound_catalytic", "bound_effectors"]:
            for enz_r, enz_p in zip(items["Enz React"], items["Enz Prod"]):
                r_dict, p_dict = (getattr(enz_r, attr), getattr(enz_p, attr))
                diff = {}
                for k in list(set(p_dict).union(set(r_dict))):
                    if k in p_dict and k in r_dict:
                        coeff = abs(p_dict[k] - r_dict[k])
                    elif k in p_dict or k in r_dict:
                        coeff = [d[k] for d in [r_dict, p_dict] 
                                 if k in d].pop()
                    if coeff != 0:
                        diff[k] = coeff

                if diff:
                    if list(diff) == list(items["Lig React"]):
                        mets = "-".join([m._remove_compartment_from_id_str() 
                                         for m in [enz_r] + list(diff)])
                        action = " binding"
                    elif list(diff) == list(items["Lig Prod"]):
                        mets = "-".join([m._remove_compartment_from_id_str() 
                                         for m in [enz_r] + list(diff)])
                        action = " release"
                    else:
                        mets = enz_r._remove_compartment_from_id_str()
                        action = " catalyzation"
                    name = mets + action

                if not name:
                    name = "-".join([enz_form._remove_compartment_from_id_str() 
                                     for enz_form in [enz_r, enz_p]])
                    name += " transition"

        # Assign the new name to the name attribute
        if update_enzyme:
            enzyme_reaction.name = name

        return name

    def unify_rate_parameters(self, reaction_list, new_parameter_id, 
                              rate_type=0, enzyme_prefix=False):
        """Unify the parameters in the rate laws for a list of reaction.

        After unification, the new parameters and rate laws are placed into the
        custom_parameters and custom_rates attributes, repsectively.

        Parameters
        ----------
        reaction_list: list of mass.MassReactions
            A list of MassReaction objects or their identifiers. MassReactions
            must already exist in the EnzymeModel.
        new_parameter_id: str
            The new parameter ID to use for the reaction parameters. The
            forward rate, reverse rate, and equilibrium constants in the
            current rate law will have the reaction ID portion of the parameter
            replaced with the new_parameter_id.
        rate_type: int {0, 1, 2, 3}, optional
            The type of rate law to display. Must be 0, 1, 2, or 3. 
                If 0, the currrent default rate law type is used. Default is 0.
                Type 1 will utilize the forward rate and equilibrium constants.
                Type 2 will utilize the forward and reverse rate constants.
                Type 3 will utilize the equilibrium and reverse rate constants.
        enzyme_prefix: bool, optional
            If True, add EnzymeModel identifier as a prefix to the 
            new_parameter_id before using the new_parameter_id in the rate
            paramter unification. Default is False.

        """
        if not isinstance(new_parameter_id, string_types):
            raise TypeError("new_parameter_id must be a str.")
        if not isinstance(enzyme_prefix, bool):
            raise TypeError("enzyme_prefix must be a bool")

        reaction_list = ensure_iterable(reaction_list)

        if enzyme_prefix:
            new_parameter_id = self.id + "_" + new_parameter_id

        try:
            reaction_list = self.reactions.get_by_any(reaction_list)
        except KeyError as e:
            raise KeyError(str(e) + " not found in model reactions")

        for reaction in reaction_list:
            if rate_type == 0:
                rate_type = reaction._rtype
            # Create a string representation of the rate and replace the 
            # reaction id portions of the parameters with new_parameter_id
            custom_rate = str(strip_time(reaction.get_rate_law(rate_type)))
            custom_rate = custom_rate.replace(reaction.id, new_parameter_id)
            self.add_custom_rate(reaction, custom_rate)

    def make_enzyme_net_flux_equation(self, enzyme_reactions, use_rates=False,
                                      update_enzyme=True):
        """Create an equation representing the net flux through the enzyme.

        The left side of the rate equation will always be the flux symbol of
        the enzyme, accessible via enzyme_flux_symbol attribute.

        Parameters
        ----------
        enzyme_reactions: iterable of mass.MassReactions
            An iterable containing the MassReactions objects or their
            identifiers to be combined.
        use_rates: bool, optional
            If True, then the rate laws of the provided reactions are used in
            creating the expression. Otherwise the arguments in the expression
            are left as MassReaction.flux_symbol sympy symbols.
        update_enzyme: bool, optional
            If True, update the enzyme_net_flux_equation attribute and, if 
            needed, the enzyme_reactions attribute of the module in addition to
            returning the generated equation. Otherwise just return the net 
            flux equation without making any updates. Default is True.

        Returns
        -------
        enzyme_net_flux_equation: A sympy expression of the net flux equation.

        """
        # Ensure enzyme_reactions are iterable and exist in model
        enzyme_rxns = ensure_iterable(enzyme_reactions)
        enzyme_rxns = [enz_rxn for enz_rxn in enzyme_rxns 
                       if enz_rxn in self._get_current_enzyme_reactions(
                           update_enzyme=update_enzyme)]

        enzyme_net_flux_equation = self.sum_enzyme_reaction_fluxes(enzyme_rxns)
        if use_rates:
            enzyme_net_flux_equation = sym.simplify(
                enzyme_net_flux_equation.subs({
                    enz_rxn.flux_symbol: enz_rxn.rate 
                    for enz_rxn in enzyme_rxns}))

        enzyme_net_flux_equation = sym.Eq(self.enzyme_flux_symbol, 
                                          enzyme_net_flux_equation)
        if update_enzyme:
            self.enzyme_net_flux_equation = enzyme_net_flux_equation

        return enzyme_net_flux_equation

    def sum_enzyme_form_concentrations(self, enzyme_forms, use_values=False):
        """Sum the enzyme form concentrations for a list of enzyme forms.

        Parameters
        ----------
        enzyme_forms: iterable of mass.EnzymeForms
            An iterable containing the EnzymeForm objects or their
            identifiers to be summed. Must already exist in the module.
        use_values: bool, optional
            If True, then numerical values are substituted into the expression. 
            Otherwise arguments in the expression are left as sympy symbols.

        Returns
        -------
        concentration: float, sympy.Basic, 
            The sum of the concentrations for the given enzyme form as a float 
            if use_values is True and all values are present. Otherwise returns
            a sympy expression representing the sum of the given enzyme forms.

        """
        enzyme_forms = ensure_iterable(enzyme_forms)
        # Get concentration formula
        concentration = self._make_summation_expr(enzyme_forms, EnzymeForm)
        # Substitute numerical values into the formula
        if use_values:
            concentration = self._sub_values_into_expr(
                strip_time(concentration), EnzymeForm)

        return concentration

    def sum_enzyme_reaction_fluxes(self, enzyme_reactions, use_values=False):
        """Sum the enzyme reaction steady state fluxes for a list of reactions.

        Parameters
        ----------
        enzyme_reactions: iterable of mass.MassReactions
            An iterable containing the MassReaction objects or their
            identifiers to be summed. Must already exist in the module and must
            be considered an enzyme reaction.
        use_values: bool, optional
            If True, then numerical values are substituted into the expression. 
            Otherwise arguments in the expression are left as sympy symbols.

        Returns
        -------
        flux: float, sympy.Basic, 
            The sum of the steady state fluxes for the given enzyme reaction as
            a float if use_values is True and all values are present. Otherwise
            returns a sympy expression representing the sum of the enzyme
            reaction fluxes.

        """
        enzyme_reactions = ensure_iterable(enzyme_reactions)
        # Get flux formula
        flux = self._make_summation_expr(enzyme_reactions, MassReaction)
        # Substitute numerical values into the formula
        if use_values:
            flux = self._sub_values_into_expr(
                strip_time(flux), MassReaction)

        return flux

    def enzyme_concentration_total_error(self, use_values=False):
        """Return the error for the total enzyme concentrations.

        The error of the total enzyme concentration is defined to be the 
        difference between the enzyme_concentration_total attribute value and
        the sum of all EnzymeForm initial conditions in the model.

        Parameters
        ----------
        use_values: bool, optional
            If True, then numerical values are substituted into the expression. 
            Otherwise arguments in the expression are left as sympy symbols.

        Returns
        -------
        error: float, sympy.Basic, 
            The error between the set enzyme_concentration_total and the sum of
            the EnzymeForm initial conditions in the model as a float 
            if use_values is True and all values are present. Otherwise returns
            a sympy expression representing the error.

        Notes
        -----
        Positive values indicate the value in the enzyme_concentration_total
        attribute is greater than the value calculated using the expression 
        from the enzyme_concentration_total_equation attribute.

        """
        if self.enzyme_concentration_total_equation is None:
            warnings.warn("No enzyme total concentration equation found. "
                          "Ensure that the model contains EnzymeForms and an e"
                          "quation for the enzyme_concentration_total_equation"
                          "attribute returns an expression")
            return None

        # Make error expression
        error = (self.enzyme_concentration_total_equation.lhs 
                 - self.enzyme_concentration_total_equation.rhs)
        # Try to substitute values into equations
        if use_values:
            error = self._sub_values_into_expr(
                strip_time(error), EnzymeForm, {
                    self.enzyme_total_symbol: self.enzyme_concentration_total})

        return error

    def enzyme_net_flux_error(self, use_values=False):
        """Return the error for the net flux through the enzyme.

        The error of the enzyme net flux is defined to be the 
        difference between the enzyme_net_flux attribute value and
        the calculated value for the enzyme_net_flux_equation attribute.

        Parameters
        ----------
        use_values: bool, optional
            If True, then numerical values are substituted into the expression. 
            Otherwise arguments in the expression are left as sympy symbols.

        Returns
        -------
        error: float, sympy.Basic, 
            The error between the set enzyme_net_flux and the calculated value
            for the enzyme_net_flux_equation attribute as a float if use_values
            is True and all values are present in the model. Otherwise returns
            a sympy expression representing the error.

        Notes
        -----
        Positive values indicate the value in enzyme_net_flux attribute is 
        greater than the value calculated using the expression from the 
        enzyme_net_flux_equation attribute.

        """
        if self.enzyme_net_flux_equation is None:
            warnings.warn("No net flux equation found. Ensure that an equation"
                          " for the enzyme_net_flux_equation attribute has "
                          "been set using the make_enzyme_net_flux_equation "
                          "method")
            return None
        # Make error expression
        error = (self.enzyme_net_flux_equation.lhs 
                 - self.enzyme_net_flux_equation.rhs)
        # Try to substitute values into equations
        if use_values:
            error = self._sub_values_into_expr(
                strip_time(error), MassReaction, {
                    self.enzyme_flux_symbol: self.enzyme_net_flux})

        return error

    def make_enzyme_fraction(self, categorized_attr, top, bottom="Equation",
                             use_values=False):
        """Make the expression for a ratio of categorized enzyme objects.

        Parameters
        ----------
        categorized_attr: str {'forms', 'reactions'}, dict
            A string representing which categorized dict attribute to use or 
            the attribute itself to use in making the enzyme ratio expression.
            Use the string 'forms' for the categorized_enzyme_forms attribute, 
            or 'reactions' for the categorized_enzyme_reactions attribute.
        top: str
            A string representing a category in the dict corresponding to the
            given categorized_attr. The summation expression of the objects in 
            the corresponding list is used as the top (numerator) of the 
            fraction to be made. Alternatively, the string "Equation" can be 
            provided to utilize an equation attribute.
        bottom: str, optional
            A string representing a category in the dict corresponding to the
            given categorized_attr. The summation expression of the objects in 
            the corresponding list is used as the bottom (denominator) of the 
            fraction to be made. Alternatively, the string 'Equation' can be 
            provided to utilize an equation attribute. Default is 'Equation'
        use_values: bool, optional
            If True, then numerical values are substituted into the expression. 
            Otherwise arguments in the expression are left as sympy symbols.

        Returns
        -------
        fraction: float, sympy.Basic, 
            The fraction either calculated and returned as float if use_values
            is True and all values are present in the model, or a sympy 
            expression representing the formula for the fraction.

        Notes
        -----
        The string "Equation" can be passed to either the top or the bottom arg
            to utilize the equation in the corresponding attribute (i.e. 
            EnzymeModel.enzyme_concentration_total_equation for 'forms'
            and EnzymeModel.enzyme_net_flux_equation for 'reactions').

        """
        # Check categorized_attr input, and get corresponding categorized dict
        if isinstance(categorized_attr, dict):
            if categorized_attr == self.categorized_enzyme_forms:
                categorized_attr = "forms"
            elif categorized_attr == self.categorized_enzyme_reactions:
                categorized_attr = "reactions"
            else:
                raise ValueError("Must be the dict accessible through '"
                                 "EnzymeModel.categorized_enzyme_forms' or '"
                                 "EnzymeModel.categorized_enzyme_reactions'.")

        if categorized_attr.lower() in {"forms", "reactions"}:
            categorized_attr = categorized_attr.lower()
            item_dict = self.__class__.__dict__[
                "categorized_enzyme_" + categorized_attr].fget(self)
        else:
            raise ValueError("Must be a string of the following: "
                             "{'forms', 'reactions'}.")

        # Get object_type for summation expression
        object_type = {"forms": EnzymeForm,
                       "reactions": MassReaction}[categorized_attr]
        # Check top & bottom inputs
        expr = sym.S.One
        for category, coeff in zip([top, bottom], [1, -1]):
            if category in item_dict:
                # Get sum of items if category is not "equation"
                summation_expr = self._make_summation_expr(
                    item_dict[category], object_type)
            elif _EQUATION_RE.match(category):
                # Get equation if category is "Equation"
                summation_expr = {
                    "forms": self.enzyme_concentration_total_equation,
                    "reactions": self.enzyme_net_flux_equation
                }[categorized_attr]
                if summation_expr is not None:
                    summation_expr = summation_expr.rhs
                else:
                    raise ValueError(
                        "No equation found for '{0}' attribute".format({
                            "forms": "enzyme_concentration_total_equation",
                            "reactions": "enzyme_net_flux_equation"
                        }[categorized_attr]))
            else:
                raise ValueError(
                    "Unrecognized category: '{0}' provided for '{1}' argument"
                    .format(category, {1: "top", -1: "bottom"}))
            # Get object type and make expression
            expr = sym.Mul(sym.Pow(summation_expr, coeff), expr)

        # Try to substitute values into equations
        if use_values:
            expr = self._sub_values_into_expr(strip_time(expr), object_type)

        return expr

    # Extended Methods
    def add_metabolites(self, metabolite_list, add_initial_conditions=False):
        """Add a list of metabolites and enzyme forms to the EnzymeModel.

        Metabolites are added under the ligand category as "Undefined".

        The change is reverted upon exit when using the EnzymeModel as a 
        context.

        Parameters
        ----------
        metabolite_list: list of mass.MassMetabolites and mass.EnzymeForms
            A list of MassMetabolites and EnzymeForms to add to the
            EnzymeModel.
        add_initial_conditons: bool, optional
            If True, the initial conditions associated with each metabolite and 
            enzyme form are also added to the model. Otherwise, the 
            metabolites and enzyme forms are added without their initial 
            conditions.

        Notes
        -----
        Extends from MassModel.add_metabolites method.

        """
        # Ensure list is iterable.
        metabolite_list = ensure_iterable(metabolite_list)

        # Add metabolites using inherited method
        super(EnzymeModel, self).add_metabolites(
            metabolite_list, add_initial_conditions)

        # Check whether a ligand is not an EnzymeForm object, and check if the
        # ligand already exists in the EnzymeModel, ignoring those that do.
        ligands = [met for met in metabolite_list 
                   if not isinstance(met, EnzymeForm)
                   and met in self.metabolites and met not in self.ligands]
        # Add ligands to the ligand DictList
        if ligands:
            self.ligands += ligands
        # Add ligands to the ligand dict as "Undefined"
        self.add_categorized_ligands("Undefined", ligands)

        # Check whether an enzyme form is a EnzymeForm object, and check if the
        # form already exists in the EnzymeModel, ignoring those that do.
        enzyme_forms = [enzyme_form for enzyme_form in metabolite_list 
                        if isinstance(enzyme_form, EnzymeForm)
                        and enzyme_form not in self.enzyme_forms]

        if enzyme_forms:
            # Add the enzyme forms to the model
            self.enzyme_forms += enzyme_forms
        # Add enzyme forms to the enzyme_form dict as "Undefined".
        self.add_categorized_enzyme_forms("Undefined", enzyme_forms)

        context = get_context(self)
        if context:
            context(partial(self.enzyme_forms.__isub__, ligands))
            context(partial(self.enzyme_forms.__isub__, enzyme_forms))

    def remove_metabolites(self, metabolite_list, destructive=False):
        """Remove a list of metabolites and enzyme forms from the EnzymeModel.

        The species' initial conditions will also be removed from the model.

        The change is reverted upon exit when using the EnzymeModel as a 
        context.

        Parameters
        ----------
        metabolite_list: list of mass.MassMetabolites and mass.EnzymeForms
            A list of species to add to the EnzymeModel.
        destructive: bool, optional
            If False, the metabolites and enzyme forms are removed from all 
            associated MassReactions.If True, also remove associated 
            MassReactions from the EnzymeModel.

        Notes
        -----
        Extends from MassModel.remove_metabolites method.

        """
        # Ensure list is iterable.
        metabolite_list = ensure_iterable(metabolite_list)

        # Set ligands as "Undefined" to reset their categories.
        ligands = [met for met in metabolite_list 
                   if not isinstance(met, EnzymeForm)
                   and met in self.metabolites and met in self.ligands]
        self.add_categorized_ligands("Undefined", ligands)

        # Remove metabolites from model using inherited method
        super(EnzymeModel, self).remove_metabolites(
            metabolite_list, destructive)

        # Remove ligands from the ligand DictList.
        if ligands:
            self.ligands -= ligands
        # Remove ligands from the "Undefined" category
        for met in ligands:
            self._categorized_ligands["Undefined"].remove(met)

        # Check whether an enzyme form is a EnzymeForm object, and check if the
        # form already exists in the EnzymeModel, ignoring those that do not.
        enzyme_forms = [enzyme_form for enzyme_form in metabolite_list
                        if isinstance(enzyme_form, EnzymeForm)
                        and enzyme_form in self.enzyme_forms]
        self.add_categorize_enzyme_forms("Undefined", enzyme_forms)
        # Remove the enzyme forms to the model
        if enzyme_forms:
            self.enzyme_forms -= enzyme_forms

        # Remove enzyme forms from the "Undefined" category
        for enzyme_form in enzyme_forms:
            self._categorized_enzyme_forms["Undefined"].remove(enzyme_form)

        context = get_context(self)
        if context:
            context(partial(self.ligands.__iadd__, enzyme_forms))
            context(partial(self.enzyme_forms.__iadd__, enzyme_forms))

    def add_reactions(self, reaction_list):
        """Add MassReactions to the EnzymeModel.

        MassReaction objects with identifiers identical to an existing reaction
        are ignored.

        The change is reverted upon exit when using the MassModel as a context

        Parameters
        ----------
        reaction_list: list of mass.MassReactions
            A list of MassReaction objects.

        """
        # Ensure list is iterable.
        reaction_list = ensure_iterable(reaction_list)

        # Add reactions using inherited method
        super(EnzymeModel, self).add_reactions(reaction_list)

        # Get the enzyme reactions by checking if an EnzymeForm is involved, 
        # and check whether reaction exists, ignoring those that do.  
        enzyme_reactions = [
            rxn for rxn in reaction_list if [
                met for met in rxn.metabolites if isinstance(met, EnzymeForm)] 
            and rxn in self.reactions and rxn not in self.enzyme_reactions]

        # Add enzyme reactions to the enzyme reaction DictList
        if enzyme_reactions:
            self.enzyme_reactions += enzyme_reactions
        # Add enzyme reactions to the enzyme reactions dict as "Undefined".
        self.add_categorized_enzyme_reactions("Undefined", enzyme_reactions)

    def remove_reactions(self, reaction_list, remove_orphans=False):
        """Remove MassReactions from the EnzymeModel.

        The change is reverted upon exit when using the MassModel as a context.

        Parameters
        ----------
        reaction_list: list of mass.MassReactions
            A list of MassReaction objects to be removed from the model.
        remove_orphans: bool, optional
            If True, will also remove orphaned genes and MassMetabolites from
            the EnzymeModel.

        """
        # Ensure list is iterable.
        reaction_list = ensure_iterable(reaction_list)
        # Get the enzyme reactions by checking if an EnzymeForm is involved, 
        # and check whether reaction exists, ignoring those that do not.
        enzyme_reactions = [
            rxn for rxn in reaction_list if [
                met for met in rxn.metabolites if isinstance(met, EnzymeForm)] 
            and rxn in self.reactions and rxn in self.enzyme_reactions]
        # Set enzyme reactions as "Undefined" to reset their categories.
        self.add_categorized_enzyme_reactions("Undefined", enzyme_reactions)

        # Remove reactions using inherited method
        super(EnzymeModel, self).remove_reactions(reaction_list)

        # Remove enzyme reactions from DictList
        if self.enzyme_reactions:
            self.enzyme_reactions -= enzyme_reactions

        # Remove enzyme reactions from the "Undefined" category
        for reaction in enzyme_reactions:
            self._categorized_enzyme_reactions["Undefined"].remove(reaction)

    def repair(self, rebuild_index=True, rebuild_relationships=True):
        """Update all indicies and pointers in the model.

        In addition to updating indicies and pointers, the enzyme_reactions
        attribute will be updated to ensure it contains all existing reactions
        involving EnzymeForms.

        Parameters
        ----------
        rebuild_index: bool, optional
            If True, then rebuild the indicies kept in the reactions,
            metabolites, and genes.
        rebuild_relationships: bool, optional
            If True, then reset all associations between the reactions,
            metabolites genes, and the MassModel, and rebuilds them.

        Notes
        -----
        Extends from MassModel.repair method.

        """
        # Repair using inherited method
        super(EnzymeModel, self).repair(rebuild_index, rebuild_relationships)
        # Repair enzyme_reactions DictList
        self._get_current_enzyme_reactions(update_enzyme=True)
        self._update_object_pointers()
        # Rebuild DictList indices
        if rebuild_index: 
            for attr in ["ligands", "enzyme_forms", "enzyme_reactions"]:
                getattr(self, attr)._generate_index()
                for value in itervalues(getattr(self, "categorized_" + attr)):
                    value._generate_index()

        for enzyme_form in self.enzymes_forms:
            enzyme_form._repair_bound_pointers()

    # Overridden methods
    def copy(self):
        """Create a partial "deepcopy" of the EnzymeModel.

        All of the MassMetabolite, MassReaction, EnzymeForm, and Gene objects,
        the initial conditions, fixed concentrations, custom_rates, and the 
        stoichiometric matrix are created anew, but in a faster fashion than 
        deepcopy.

        Notes
        -----
        Overrides MassModel.copy method.

        """
        # Define a new model
        new_model = self.__class__()
        # Define items that will not be copied by their references
        do_not_copy_by_ref = {
            "metabolites", "reactions", "genes", "initial_conditions", "_S", 
            "custom_rates", "ligands", "enzyme_forms", "enzyme_reactions", 
            "_categorized_ligands", "_categorized_enzyme_forms", 
            "_categorized_enzyme_reactions", "notes", "annotation"}
        for attr in self.__dict__:
            if attr not in do_not_copy_by_ref:
                new_model.__dict__[attr] = self.__dict__[attr]
        new_model.notes = deepcopy(self.notes)
        new_model.annotation = deepcopy(self.annotation)

        # Copy the metabolites
        new_model = self._copy_model_metabolites(new_model)
        # Copy the genes
        new_model = self._copy_model_genes(new_model)
        # Copy the reactions and rates
        new_model = self._copy_model_reactions(new_model)
        # Create the new stoichiometric matrix for the model.
        new_model._S = self._mk_stoich_matrix(matrix_type=self._matrix_type,
                                              dtype=self._dtype,
                                              update_model=True)
        # Add the newly copied objects to their appropriate DictLists 
        # in the ligands, enzyme_forms and enzyme_reactions attributes
        for metabolite in new_model.metabolites:
            if isinstance(metabolite, EnzymeForm):
                new_model.enzyme_forms.append(metabolite)
            else:
                new_model.ligands.append(metabolite)
        new_model._get_current_enzyme_reactions(update_enzyme=True)
        new_model._update_object_pointers()

        # Doesn't make sense to retain the context of a copied model so
        # assign a new empty context
        new_model._contexts = []
        return new_model

    def merge(self, right, prefix_existing=None, inplace=False,
              new_model_id=None):
        """TODO DOCSTRING."""
        # Merge the two EnzymeModels together if right is an EnzymeModel
        if isinstance(right, EnzymeModel):
            print("TODO: FINISH ENZYMEMODULE MERGE")
            new_model = None
        else:
            # Always merge the EnzymeModel into the MassModel
            new_model = right.merge(self, prefix_existing=prefix_existing, 
                                    inplace=inplace, new_model_id=new_model_id)
        return new_model

    # Internal
    def _update_object_pointers(self):
        """Update objects in the attributes to point to the model.

        Warnings
        --------
        This method is intended for internal use only. 

        """
        for attr in ["ligands", "enzyme_forms", "enzyme_reactions"]:
            model_dictlist = {
                "ligands": self.metabolites, 
                "enzyme_forms": self.metabolites,
                "enzyme_reactions": self.reactions}.get(attr)
            setattr(self, attr, 
                    _mk_new_dictlist(model_dictlist, getattr(self, attr)))
            attr = "_categorized_" + attr
            setattr(self, attr, {
                key: _mk_new_dictlist(model_dictlist, old_dictlist)
                for key, old_dictlist in iteritems(getattr(self, attr))})

    def _get_current_enzyme_reactions(self, update_enzyme=False):
        """Get the enzyme reactions that currently exist in the model.

        Parameters
        ----------
        update_enzyme: bool, optional
            If True, update the enzyme_reactions attribute of the EnzymeModel.

        Warnings
        --------
        This method is intended for internal use only. 

        """
        enzyme_reactions = []
        for enzyme_form in self.enzyme_forms:
            enzyme_reactions.extend([
                rxn for rxn in list(enzyme_form._reaction) 
                if rxn not in enzyme_reactions])

        enzyme_reactions = DictList(enzyme_reactions)
        enzyme_reactions.sort()

        # Update enzyme_reaction attribute if necessary
        if set(enzyme_reactions) ^ set(self.enzyme_reactions) \
           and update_enzyme:
            self.enzyme_reactions = enzyme_reactions

        return enzyme_reactions

    def _make_summation_expr(self, items, object_type):
        """Create a sympy expression of the summation of the given items.

        Warnings
        --------
        This method is intended for internal use only. 

        """
        # Get appropriate list based on object type
        item_list = {
            EnzymeForm: self.enzyme_forms, 
            MassReaction: self.enzyme_reactions}.get(object_type)

        # Ensure object_type are appropriate type and exist in model
        items = [item for item in items if isinstance(item, object_type) 
                 and item in item_list]

        return sum([x.flux_symbol if isinstance(x, MassReaction)
                    else _mk_met_func(x.id) for x in items])

    def _sub_values_into_expr(self, expr, object_type, additional=None):
        """Substitute values into an expression and try to return a float.

        Warnings
        --------
        This method is intended for internal use only. 

        """
        expr = strip_time(expr)
        # Get the dictlist of objects and the values for substitution
        values = {str(k): v for k, v in iteritems(self.initial_conditions)}
        if object_type == MassReaction:
            values.update(self._get_all_parameters())
        values = {str(key): values[str(key)] 
                  for key in list(expr.atoms(sym.Symbol)) 
                  if str(key) in values}

        if additional is not None:
            values.update({str(k): v for k, v in iteritems(additional)})

        # Substitute values
        expr = expr.subs(values)

        # Attempt to convert expression to a float.
        try:
            expr = float(expr)
        except TypeError:
            message_strs = {
                EnzymeForm: ("forms", "initial conditions"),
                MassReaction: ("reactions", "steady state fluxes")
            }.get(object_type)
            warnings.warn("Not all enzyme {0} have {1} defined in model, will "
                          "return a sympy expression".format(*message_strs))

        return expr

    def _add_categorized_items(self, attribute, category, items):
        """Create a category for items in the given attribute.

        Warnings
        --------
        This method is intended for internal use only. 

        """
        if not isinstance(category, string_types):
            raise TypeError("category must be a str")
        items = ensure_iterable(items)

        if category in getattr(self, attribute):
            items = set(getattr(self, attribute)[category]).union(set(items))

        self.__class__.__dict__[attribute[1:]].fset(self, {category: items})

    def _remove_empty_categories(self, attribute):
        """Remove categories with empty lists from an attribute dict.

        Warnings
        --------
        This method is intended for internal use only. 

        """
        # Get the attribute to modify
        attr_dict = getattr(self, attribute)
        # Remove categorizies with empty lists
        to_remove = [c for c, items in iteritems(attr_dict) if not items]
        for category in to_remove:
            del attr_dict[category]

        return attr_dict

    def _set_category_attribute_dict(self, value, attribute):
        """Set the ligand, enzyme_form, or enzyme_reaction attribute_dict.

        Warnings
        --------
        This method is intended for internal use only. 

        """
        if not isinstance(value, dict):
            raise TypeError("value must be a dict")

        # Get the attribute to modify
        attr_dict = getattr(self, attribute)
        # Reset items to the undefined category if given an empty dict
        if not value:
            items = [item for lst in itervalues(attr_dict) for item in lst]
            setattr(self, attribute, {"Undefined": DictList()})
            self._set_category_attribute_dict({"Undefined": items}, attribute)

        # Add items into their category lists
        undefined = []
        for category, items in iteritems(value):
            if not items:
                continue
            items = ensure_iterable(items)
            dictlist_dict = {"_categorized_ligands": "ligands",
                             "_categorized_enzyme_forms": "enzyme_forms",
                             "_categorized_enzyme_reactions": "reactions"}
            dictlist = getattr(self, dictlist_dict[attribute])
            # Try to make a new DictList with correct references.
            try:
                items = _mk_new_dictlist(dictlist, items, ensure_unique=True)
            except KeyError as e:
                raise KeyError(str(e) + " not in model " 
                               + dictlist_dict[attribute])
            items.sort()
            if _UNDEFINED_RE.match(category):
                undefined += items
            attr_dict[category] = items

        # Set the categorized items if not being reset as undefined.
        categorized = []
        for category, items in iteritems(attr_dict):
            if not _UNDEFINED_RE.match(category):
                categorized += [x for x in items if x not in undefined]
                attr_dict[category] = DictList([
                    item for item in attr_dict[category] 
                    if item not in undefined])

        # Add metabolites without a role into the "Undefined" role
        if "Undefined" in attr_dict:
            categorized = list(set(categorized))  # Only unique necessary
            # Set undefined metabolites
            undefined += [x for x in attr_dict["Undefined"] 
                          if x not in categorized]
            undefined = DictList(set(undefined))
            undefined.sort()
            attr_dict["Undefined"] = undefined

        self._remove_empty_categories(attribute)

    def _convert_self_into_enzyme_dict(self):
        """Convert self into an EnzymeDict.

        Primarily used for merging an EnzymeModel into a MassModel while 
        retaining defined EnzymeModel attributes. Also used for checking if
        subclass is an EnzymeModel to avoid circular imports

        Warnings
        --------
        This method is intended for internal use only. 

        """
        return EnzymeDict(self)

    def _add_self_to_model(self, model, prefix_existing=None, inplace=False,
                           new_model_id=None):
        """Add self to the model and return the MassModel object.

        Warnings
        --------
        This method is intended for internal use only. 

        """
        # Create a MassModel instance of self to merge normally
        model = model.merge(MassModel(self), prefix_existing=prefix_existing, 
                            inplace=inplace, new_model_id=new_model_id)
        # Turn EnzymeModel into an EnzymeDict to store in MassModel.enzymes
        enzyme_dict = self._convert_self_into_enzyme_dict()
        # Update EnzymeDict with attributes
        enzyme_dict._update_object_pointers(model)

        if enzyme_dict.id in model.enzymes:
            model.enzymes.remove(enzyme_dict.id)
        model.enzymes.append(enzyme_dict)
        enzyme_dict.model = model

        return model

    def _repr_html_(self):
        """HTML representation of the overview for the EnzymeModel.

        Overrides MassModel._repr_html method.
        """
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
                    <td><strong>Subsystem</strong></td>
                    <td>{subsystem}</td>
                </tr><tr>
                    <td><strong>Number of Ligands</strong></td>
                    <td>{num_ligands}</td>
                </tr><tr>
                    <td><strong>Number of EnzymeForms</strong></td>
                    <td>{num_enz_forms}</td>
                </tr><tr>
                    <td><strong>Number of Enzyme Reactions</strong></td>
                    <td>{num_enz_reactions}</td>
                </tr><tr>
                    <td><strong>Total Enzyme Concentration</strong></td>
                    <td>{enz_conc}</td>
                </tr><tr>
                    <td><strong>Enzyme Net Flux</strong></td>
                    <td>{enz_flux}</td>
                </tr><tr>
                    <td><strong>Number of Initial Conditions</strong></td>
                    <td>{num_ic}</td>
                </tr><tr>
                    <td><strong>Number of Parameters</strong></td>
                    <td>{num_parameters}</td>
                </tr><tr>
                    <td><strong>Number of Irreversible Reactions</strong></td>
                    <td>{num_irreversible}</td>
                </tr><tr>
                    <td><strong>Number of Custom Rates</strong></td>
                    <td>{num_custom_rates}</td>
                </tr><tr>
                    <td><strong>Compartments</strong></td>
                    <td>{compartments}</td>
                </tr><tr>
                    <td><strong>Units</strong></td>
                    <td>{units}</td>
                </tr>
            </table>
        """.format(name=self.id, address='0x0%x' % id(self),
                   dim_stoich_mat=dim_S, mat_rank=rank,
                   S_type="{}, {}".format(self._matrix_type,
                                          self._dtype.__name__),
                   subsystem=self.subsystem,
                   num_ligands=len(self.ligands),
                   num_enz_forms=len(self.enzyme_forms),
                   num_enz_reactions=len(self.enzyme_reactions),
                   enz_conc=self.enzyme_concentration_total, 
                   enz_flux=self.enzyme_net_flux,
                   num_ic=len(self.initial_conditions),
                   num_parameters=len(self._get_all_parameters()),
                   num_irreversible=len(self.irreversible_reactions),
                   num_custom_rates=len(self.custom_rates),
                   compartments=", ".join(v if v else k for k, v in
                                          iteritems(self.compartments)),
                   units=", ".join(v if v else k for
                                   k, v in iteritems(self.units)))
