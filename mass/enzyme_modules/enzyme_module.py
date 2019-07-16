# -*- coding: utf-8 -*-
"""TODO Module Docstrings."""
import re
from collections import defaultdict
from copy import deepcopy
from functools import partial
from warnings import warn

import numpy as np

from six import integer_types, iteritems, iterkeys, itervalues, string_types

import sympy as sym

from cobra.core.dictlist import DictList
from cobra.util.context import get_context

from mass.core.massmodel import MassModel
from mass.enzyme_modules.enzyme_module_dict import EnzymeModuleDict
from mass.enzyme_modules.enzyme_module_form import EnzymeModuleForm
from mass.enzyme_modules.enzyme_module_reaction import EnzymeModuleReaction
from mass.util.expressions import _mk_met_func, strip_time
from mass.util.util import _mk_new_dictlist, ensure_iterable

_AUTOMATIC_RE = re.compile("^Automatic$")
_UNDEFINED_RE = re.compile("^Undefined$")
_EQUATION_RE = re.compile("^Equation$")


class EnzymeModule(MassModel):
    """Class representation of an EnzymeModule.

    Parameters
    ----------
    id_or_model: str, MassModel, EnzymeModule
        Either an identifier to associate with the EnzymeModule given as a
        string, or an existing EnzymeModule object. If an existing EnzymeModule
        or MassModel is provided, a new EnzymeModule is instantiated with the
        same properties as the original object.
    name: str, optional
        A human readable name for the EnzymeModule.
    subsystem: str, optional
        The subsystem where the enzyme catalyzed net reaction that is
        represented by this EnzymeModule is meant to occur.
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
    enzyme_module_ligands: cobra.DictList
        A cobra.DictList where the keys are the metabolite identifiers and the
        values are the associated MassMetabolite objects.
    enzyme_module_forms: cobra.DictList
        A cobra.DictList where the keys are the enzyme form identifiers and the
        values are the associated EnzymeModuleForm objects.
    enzyme_module_reactions: cobra.DictList
        A cobra.DictList where keys are the EnzymeModuleReaction identifiers
        and the values are the associated EnzymeModuleReactions objects. A
        reaction is considered an EnzymeModuleReaction if it involves an
        EnzymeModuleForm object.
    enzyme_module_ligands_categorized: dict
        A dict of user-categorized ligands where keys are categories and
        values are DictLists of corresponding MassMetabolites.
    enzyme_module_forms_categorized: dict
        A dict of user-categorized enzyme forms where keys are categories and
        values are DictLists of the corresponding EnzymeModuleForms.
    enzyme_module_reactions_categorized: dict
        A dict of user-categorized reactions involving EnzymeModuleForms where
        keys are categories and values are DictLists of the corresponding
        EnzymeModuleReaction.
    enzyme_concentration_total_equation: sympy.Basic
        A sympy expression representing the net reaction rate equation for the
        enzyme represented by the EnzymeModule.
    enzyme_net_flux_equation: sympy.Basic
        A sympy expression representing the net reaction rate equation for the
        enzyme represented by the EnzymeModule.
    enzyme_concentration_total: dict
        A dict containing the total enzyme concentration symbol as a sympy
        symbol and the total enzyme concentration value as a float.
    enzyme_net_flux: dict
        A dict containing the enzyme net flux symbol as a sympy
        symbol and the enzyme net flux value as a float.

    """

    # pylint: disable=too-many-instance-attributes
    def __init__(self, id_or_model=None, name=None, subsystem="",
                 matrix_type="dense", dtype=np.float64):
        """Initialize the EnzymeModule Object."""
        super(EnzymeModule, self).__init__(
            id_or_model=id_or_model, name=name, matrix_type=matrix_type,
            dtype=dtype)

        self.subsystem = subsystem
        # Initialize DictLists for enzyme ligands, forms, and reactions
        self.enzyme_module_ligands = DictList()
        self.enzyme_module_forms = DictList()
        self.enzyme_module_reactions = DictList()

        # Initialize a dict of DictLists for storing categorized objects
        self._enzyme_module_ligands_categorized = {"Undefined": DictList()}
        self._enzyme_module_forms_categorized = {"Undefined": DictList()}
        self._enzyme_module_reactions_categorized = {"Undefined": DictList()}

        # Initialize EnzymeModule attributes
        self._enzyme_concentration_total = None
        self._enzyme_net_flux = None
        self.enzyme_net_flux_equation = None

    @property
    def enzyme_total_symbol(self):
        """Return the sympy symbol for the total enzyme concentration."""
        if self.id is None:
            return sym.Symbol("Enzyme_Total")

        return sym.Symbol(self.id + "_Total")

    @property
    def enzyme_flux_symbol(self):
        """Return the sympy symbol for the net flux through the enzyme."""
        if self.id is None:
            return sym.Symbol("v")

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
        if value is None:
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

        setattr(self, "_enzyme_net_flux", value)

    @property
    def enzyme_concentration_total_equation(self):
        """Return the total concentration equation as a sympy.Equation.

        Will first try to sum forms based on the enzyme module forms that have
        enzyme_module_id attributes that match the EnzymeModule id attribute.

        If no enzyme modules forms with matching enzyme_module_id attributes
        are found, then try to sum all enzyme module forms in the the model.
        """
        if not self.enzyme_module_forms:
            return None

        # First try only using enzyme module forms that reference this module
        enzyme_module_forms = [
            form for form in self.enzyme_module_forms
            if form.enzyme_module_id == self.id]
        # If None found, use all enzyme module forms in the EnzymeModule.
        if not enzyme_module_forms:
            enzyme_module_forms = self.enzyme_module_forms

        return sym.Eq(self.enzyme_total_symbol,
                      self.sum_enzyme_module_form_concentrations(
                          enzyme_module_forms, use_values=False))

    @property
    def enzyme_net_flux_equation(self):
        """Return the net rate equation of the enzyme as a sympy.Equation."""
        value = getattr(self, "_enzyme_net_flux_equation", None)
        if value is not None:
            value = sym.Eq(self.enzyme_flux_symbol, value)
        return value

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
            if isinstance(value, string_types):
                value = sym.sympify(value)
            elif hasattr(value, "lhs") and hasattr(value, "rhs"):
                if value.lhs == self.enzyme_flux_symbol:
                    value = value.rhs
                if value.rhs == self.enzyme_flux_symbol:
                    value = value.lhs
            else:
                pass
        setattr(self, "_enzyme_net_flux_equation", value)

    @property
    def enzyme_module_ligands_categorized(self):
        """Return a dict of ligands in user-defined categories."""
        return self._remove_empty_categories(
            "_enzyme_module_ligands_categorized")

    @enzyme_module_ligands_categorized.setter
    def enzyme_module_ligands_categorized(self, value):
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
        A metabolite must already exist in the EnzymeModule in order to set its
        category. Categories with empty lists are removed.

        See Also
        --------
        add_metabolites:
            Method to add metabolites to the model.
        categorize_enzyme_module_ligands:
            Method to categorize "Undefined" ligands.

        """
        if not isinstance(value, dict):
            raise TypeError("value must be a dict")

        self._set_category_attribute_dict(
            value, "_enzyme_module_ligands_categorized")

    @property
    def enzyme_module_forms_categorized(self):
        """Return a dict of enzyme forms in user-defined categories."""
        return self._remove_empty_categories(
            "_enzyme_module_forms_categorized")

    @enzyme_module_forms_categorized.setter
    def enzyme_module_forms_categorized(self, value):
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
        An enzyme form must already exist in the EnzymeModule in order to set
        its category. Categories with empty lists are removed.

        See Also
        --------
        add_metabolites:
            Method to add metabolites or EnzymeModuleForms to the model.
        categorize_enzyme_module_forms:
            Method to categorize "Undefined" enzyme forms.

        """
        if not isinstance(value, dict):
            raise TypeError("value must be a dict")

        self._set_category_attribute_dict(
            value, "_enzyme_module_forms_categorized")

    @property
    def enzyme_module_reactions_categorized(self):
        """Return the enzyme module reactions in user-defined categories."""
        return self._remove_empty_categories(
            "_enzyme_module_reactions_categorized")

    @enzyme_module_reactions_categorized.setter
    def enzyme_module_reactions_categorized(self, value):
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
        A reaction must already exist in the EnzymeModule in order to set its
        category. Categories with empty lists are removed.

        See Also
        --------
        add_reactions:
            Method to add reactions to the model.
        categorize_enzyme_module_reactions:
            Method to categorize "Undefined" binding reactions.

        """
        if not isinstance(value, dict):
            raise TypeError("value must be a dict")

        self._set_category_attribute_dict(
            value, "_enzyme_module_reactions_categorized")

    def set_enzyme_object_category(self, category, object_list):
        """Add a list of EnzymeModule objects to a new or existing category.

        If a category already exists, the objects will be added to the existing
        list. If an object is categorized as "Undefined", it will be removed
        from all other existing categories.

        Parameters
        ----------
        category: str
            A string representing the category for the list of objects to be
            categorized.
        object_list: iterable
            An iterable containing the MassMetabolite objects, the
            EnzymeModuleForm objects, EnzymeModuleReaction objects, or their
            identifiers to be categorized.

        """
        if not isinstance(category, string_types):
            raise TypeError("category must be a str")
        object_list = ensure_iterable(object_list)
        sep_objs = defaultdict(list)
        for obj in object_list:
            if isinstance(obj, EnzymeModuleForm.__base__)\
               and not isinstance(obj, EnzymeModuleForm):
                sep_objs["_enzyme_module_ligands_categorized"].append(obj)
            if isinstance(obj, EnzymeModuleForm):
                sep_objs["_enzyme_module_forms_categorized"].append(obj)
            if isinstance(obj, EnzymeModuleReaction):
                sep_objs["_enzyme_module_reactions_categorized"].append(obj)

        for key, object_list in iteritems(sep_objs):
            if category in getattr(self, key):
                object_list = set(
                    getattr(self, key)[category]).union(set(object_list))

            self.__class__.__dict__[key[1:]].fset(
                self, {category: object_list})

    def make_enzyme_module_form(self, id=None, name="automatic",
                                categories="Undefined", bound_catalytic=None,
                                bound_effectors=None, compartment=None):
        """Make an EnzymeModuleForm object to add to the EnzymeModule.

        Parameters
        ----------
        id: str
            The identifier associated with the EnzymeModuleForm.
        name: str, optional
            A human readable name for the enzyme form. If name is set to match
            "Automatic", a name will be generated based on the EnzymeModuleForm
            and its bound ligands.
        categories: str, iterable of str
            A string representing the category, or an iterable of strings
            containing several categories for the EnzymeModuleForm.
        bound_catalytic: dict, optional
            A dict representing the metabolites bound to the enzyme's active
            site(s), with MassMetabolites or their identifiers as keys and the
            number bound as values.
        bound_effectors: dict, optional
            A dict representing the metabolites bound to the enzyme's
            regulatory site(s), with MassMetabolites or their identifiers as
            keys and the number bound as values.
        compartment: str, optional
            The compartment where the EnzymeModuleForm is located.

        Returns
        -------
        enzyme_module_form: EnzymeModuleForm
            The newly created EnzymeModuleForm object.

        Notes
        -----
        If a metabolite does not exist in the EnzymeModule, it will be added to
            the module in addition to the EnzymeModuleForm

        See Also
        --------
        EnzymeModuleForm.generate_enzyme_module_form_name:
            Automatic generation of name attribute for EnzymeModuleForm

        """
        # Ensure metabolites for EnzymeModuleForms exist in the EnzymeModule.
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
        # Make EnzymeModuleForm object
        enzyme_module_form = EnzymeModuleForm(
            id_or_specie=id, name=name, enzyme_module_id=self.id,
            bound_catalytic=bound_catalytic, bound_effectors=bound_effectors,
            compartment=compartment)
        # Generate a name for the EnzymeModuleForm if name set to "Automatic"
        if _AUTOMATIC_RE.match(name):
            enzyme_module_form.generate_enzyme_module_form_name(True)

        # Add the enzyme form to the module and place in respective categories
        self.add_metabolites(enzyme_module_form)
        categories = ensure_iterable(categories)
        for category in categories:
            self.set_enzyme_object_category(category, enzyme_module_form)

        return enzyme_module_form

    def make_enzyme_module_reaction(self, id=None, name="", subsystem=None,
                                    reversible=True, categories="Undefined",
                                    metabolites_to_add=None):
        """Make an EnzymeModuleReaction object to add to the EnzymeModule.

        Parameters
        ----------
        id: str
            The identifier associated with the EnzymeModuleReaction.
        name: str, optional
            A human readable name for the enzyme reaction. If name is set to
            match "Automatic", a name will be generated based on the
            EnzymeModuleForms and their bound ligands.
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
            A dict with MassMetabolite and EnzymeModuleForm objects or their
            identifiers as keys and stoichiometric coefficients as values. If
            keys are string identifiers of the objects, the MassMetabolite and
            EnzymeModuleForm objects must already be a part of a MassModel.

        Returns
        -------
        enzyme_module_reaction: EnzymeModuleReaction
            The newly created EnzymeModuleReaction object.

        Notes
        -----
        A final coefficient of < 0 implies a reactant and a final coefficient
            of > 0 implies a product.

        See Also
        --------
        EnzymeModuleReaction.generate_enzyme_module_reaction_name:
            Automatic generation of name attribute for enzyme module reactions.

        """
        # Make EnzymeModuleReaction object
        new_reaction = EnzymeModuleReaction(
            id_or_reaction=id, name=name, subsystem=subsystem,
            reversible=reversible, enzyme_module_id=self.id)

        # Add metabolites
        if metabolites_to_add:
            try:
                metabolites_to_add = dict(
                    (met, c) if isinstance(met, EnzymeModuleForm.__base__)
                    else (self.metabolites.get_by_id(met), c)
                    for met, c in iteritems(metabolites_to_add))
            except KeyError as e:
                raise KeyError(str(e) + "not found in model metabolites.")
            new_reaction.add_metabolites(metabolites_to_add)

        # Add reaction to EnzymeModule
        self.add_reactions(new_reaction)

        # Add categories for reaction
        categories = ensure_iterable(categories)
        for category in categories:
            self.set_enzyme_object_category(category, new_reaction)

        # Set enzyme name if set to Automatic
        if _AUTOMATIC_RE.match(name):
            name = new_reaction.generate_enzyme_module_reaction_name(True)

        return new_reaction

    def unify_rate_parameters(self, reaction_list, new_parameter_id,
                              rate_type=0, enzyme_prefix=False):
        """Unify the parameters in the rate laws for a list of reaction.

        After unification, the new parameters and rate laws are placed into the
        custom_parameters and custom_rates attributes, repsectively.

        Parameters
        ----------
        reaction_list: list of EnzymeModuleReactions
            A list of EnzymeModuleReaction objects or their identifiers.
            EnzymeModuleReactions must already exist in the EnzymeModule.
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
            If True, add the EnzymeModule identifier as a prefix to the
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

        for reaction in reaction_list:
            if not isinstance(reaction, EnzymeModuleReaction):
                try:
                    reaction = self.reactions.get_by_id(reaction)
                except KeyError as e:
                    raise KeyError(str(e) + " not found in model reactions")
            if rate_type == 0:
                rate_type = reaction._rtype
            # Create a string representation of the rate and replace the
            # reaction id portions of the parameters with new_parameter_id
            custom_rate = str(strip_time(
                reaction.get_mass_action_rate(rate_type)))
            custom_rate = custom_rate.replace(reaction.id, new_parameter_id)
            self.add_custom_rate(reaction, custom_rate)

    def make_enzyme_net_flux_equation(self, enzyme_module_reactions,
                                      use_rates=False, update_enzyme=False):
        """Create an equation representing the net flux through the enzyme.

        The left side of the rate equation will always be the flux symbol of
        the enzyme, accessible via enzyme_flux_symbol attribute.

        Parameters
        ----------
        enzyme_module_reactions: iterable of EnzymeModuleReactions
            An iterable containing the EnzymeModuleReactions objects or their
            identifiers to be combined.
        use_rates: bool, optional
            If True, then the rate laws of the provided reactions are used in
            creating the expression. Otherwise the arguments in the expression
            are left as EnzymeModuleReaction.flux_symbol sympy symbols.
        update_enzyme: bool, optional
            If True, update the enzyme_net_flux_equation attribute and, if
            needed, the enzyme_module_reactions attribute of the module in
            addition to returning the generated equation. Otherwise just return
            the net flux equation without making any updates. Default is False.

        Returns
        -------
        enzyme_net_flux_equation: A sympy expression of the net flux equation.

        """
        # Ensure enzyme_module_reactions are iterable and exist in model
        enzyme_rxns = ensure_iterable(enzyme_module_reactions)
        enzyme_rxns = [
            enz_rxn for enz_rxn in enzyme_rxns
            if enz_rxn in self._get_current_enzyme_module_objs(
                "reactions", update_enzyme=update_enzyme)]

        enzyme_net_flux_equation = self.sum_enzyme_module_reaction_fluxes(
            enzyme_rxns)
        if use_rates:
            enzyme_net_flux_equation = sym.simplify(
                enzyme_net_flux_equation.subs({
                    enz_rxn.flux_symbol: enz_rxn.rate
                    for enz_rxn in enzyme_rxns}))

        if update_enzyme:
            self.enzyme_net_flux_equation = enzyme_net_flux_equation

        return sym.Eq(self.enzyme_flux_symbol, enzyme_net_flux_equation)

    def sum_enzyme_module_form_concentrations(self, enzyme_module_forms,
                                              use_values=False):
        """Sum the enzyme form concentrations for a list of enzyme forms.

        Parameters
        ----------
        enzyme_module_forms: iterable of EnzymeModuleForms
            An iterable containing the EnzymeModuleForm objects or their
            identifiers to be summed. Must already exist in the module.
        use_values: bool, optional
            If True, then numerical values are substituted into the expression.
            Otherwise arguments in the expression are left as sympy symbols.

        Returns
        -------
        concentration: float, sympy.Basic
            The sum of the concentrations for the given enzyme form as a float
            if use_values is True and all values are present. Otherwise returns
            a sympy expression representing the sum of the given enzyme forms.

        """
        enzyme_module_forms = ensure_iterable(enzyme_module_forms)
        # Get concentration formula
        concentration = self._make_summation_expr(
            enzyme_module_forms, EnzymeModuleForm)
        # Substitute numerical values into the formula
        if use_values:
            concentration = self._sub_values_into_expr(
                strip_time(concentration), EnzymeModuleForm)

        return concentration

    def sum_enzyme_module_reaction_fluxes(self, enzyme_module_reactions,
                                          use_values=False):
        """Sum the enzyme reaction steady state fluxes for a list of reactions.

        Parameters
        ----------
        enzyme_module_reactions: iterable of EnzymeModuleReactions
            An iterable containing the EnzymeModuleReaction objects or their
            identifiers to be summed. Must already exist in the module and must
            be considered an enzyme reaction.
        use_values: bool, optional
            If True, then numerical values are substituted into the expression.
            Otherwise arguments in the expression are left as sympy symbols.

        Returns
        -------
        flux: float, sympy.Basic
            The sum of the steady state fluxes for the given enzyme reaction as
            a float if use_values is True and all values are present. Otherwise
            returns a sympy expression representing the sum of the enzyme
            reaction fluxes.

        """
        enzyme_module_reactions = ensure_iterable(enzyme_module_reactions)
        # Get flux formula
        flux = self._make_summation_expr(
            enzyme_module_reactions, EnzymeModuleReaction)
        # Substitute numerical values into the formula
        if use_values:
            flux = self._sub_values_into_expr(
                strip_time(flux), EnzymeModuleReaction)

        return flux

    def enzyme_concentration_total_error(self, use_values=False):
        """Return the error for the total enzyme concentrations.

        The error of the total enzyme concentration is defined to be the
        difference between the enzyme_concentration_total attribute value and
        the sum of all EnzymeModuleForm initial conditions in the model.

        Parameters
        ----------
        use_values: bool, optional
            If True, then numerical values are substituted into the expression.
            Otherwise arguments in the expression are left as sympy symbols.

        Returns
        -------
        error: float, sympy.Basic
            The error between the set enzyme_concentration_total and the sum of
            the EnzymeModuleForm initial conditions in the model as a float
            if use_values is True and all values are present. Otherwise returns
            a sympy expression representing the error.

        Notes
        -----
        Positive values indicate the value in the enzyme_concentration_total
        attribute is greater than the value calculated using the expression
        from the enzyme_concentration_total_equation attribute.

        """
        if self.enzyme_concentration_total_equation is None:
            warn("No enzyme total concentration equation found. Ensure that "
                 "the model contains EnzymeModuleForms and an equation for the"
                 "enzyme_concentration_total_equation attribute returns an "
                 "expression")
            return None

        # Make error expression
        error = (self.enzyme_concentration_total_equation.lhs
                 - self.enzyme_concentration_total_equation.rhs)
        # Try to substitute values into equations
        if use_values:
            error = self._sub_values_into_expr(
                strip_time(error), EnzymeModuleForm, {
                    self.enzyme_total_symbol: self.enzyme_concentration_total})

        return error

    def enzyme_net_flux_error(self, use_values=False):
        """Return the error for the net flux through the enzyme.

        The error of the enzyme net flux is defined to be the difference
        between the enzyme_net_flux attribute value and the calculated value
        for the enzyme_net_flux_equation attribute.

        Parameters
        ----------
        use_values: bool, optional
            If True, then numerical values are substituted into the expression.
            Otherwise arguments in the expression are left as sympy symbols.

        Returns
        -------
        error: float, sympy.Basic
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
            warn("No net flux equation found. Ensure that an equation for the "
                 "enzyme_net_flux_equation attribute has been set.")
            return None
        # Make error expression
        error = (self.enzyme_net_flux_equation.lhs
                 - self.enzyme_net_flux_equation.rhs)
        # Try to substitute values into equations
        if use_values:
            error = self._sub_values_into_expr(
                strip_time(error), EnzymeModuleReaction, {
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
            Use the string 'forms' for the enzyme_module_forms_categorized
            attribute, or 'reactions' for the
            enzyme_module_reactions_categorized attribute.
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
        fraction: float, sympy.Basic
            The fraction either calculated and returned as float if use_values
            is True and all values are present in the model, or a sympy
            expression representing the formula for the fraction.

        Notes
        -----
        The string "Equation" can be passed to either the top or the bottom arg
            to utilize the equation in the corresponding attribute (i.e.
            EnzymeModule.enzyme_concentration_total_equation for 'forms'
            and EnzymeModule.enzyme_net_flux_equation for 'reactions').

        """
        # Check categorized_attr input, and get corresponding categorized dict
        if isinstance(categorized_attr, dict):
            if categorized_attr == self.enzyme_module_forms_categorized:
                categorized_attr = "forms"
            elif categorized_attr == self.enzyme_module_reactions_categorized:
                categorized_attr = "reactions"
            else:
                raise ValueError(
                    "Must be the dict accessible through '"
                    "EnzymeModule.enzyme_module_forms_categorized' or "
                    "'EnzymeModule.enzyme_module_reactions_categorized'.")

        if categorized_attr.lower() in {"forms", "reactions"}:
            categorized_attr = categorized_attr.lower()
            item_dict = self.__class__.__dict__[
                "enzyme_module_" + categorized_attr + "_categorized"].fget(
                    self)
        else:
            raise ValueError("Must be a string of the following: "
                             "{'forms', 'reactions'}.")

        # Get object_type for summation expression
        object_type = {"forms": EnzymeModuleForm,
                       "reactions": EnzymeModuleReaction}[categorized_attr]
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
    def add_metabolites(self, metabolite_list):
        """Add a list of metabolites and enzyme forms to the EnzymeModule.

        Metabolites are added under the ligand category as "Undefined".

        The change is reverted upon exit when using the EnzymeModule as a
        context.

        Parameters
        ----------
        metabolite_list: list of MassMetabolites and EnzymeModuleForms
            A list of MassMetabolites and EnzymeModuleForms to add to the
            EnzymeModule.
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
        super(EnzymeModule, self).add_metabolites(metabolite_list)

        # Get items that are not EnzymeModuleForm objects, and check if the
        # ligands already exists in the EnzymeModule, ignoring those that do.
        ligands = [met for met in metabolite_list
                   if not isinstance(met, EnzymeModuleForm)
                   and met in self.metabolites
                   and met not in self.enzyme_module_ligands]
        # Add ligands to the ligand DictList
        if ligands:
            self.enzyme_module_ligands += ligands
        # Add ligands to the ligand dict as "Undefined"
        self.set_enzyme_object_category("Undefined", ligands)

        # Get items that are EnzymeModuleForm objects, and check if the enzyme
        # forms already exist in the EnzymeModule, ignoring those that do.
        enzyme_module_forms = [
            enzyme_module_form for enzyme_module_form in metabolite_list
            if isinstance(enzyme_module_form, EnzymeModuleForm)
            and enzyme_module_form not in self.enzyme_module_forms]

        if enzyme_module_forms:
            # Add the enzyme forms to the model
            self.enzyme_module_forms += enzyme_module_forms
        # Add enzyme forms to the enzyme_module_form dict as "Undefined".
        self.set_enzyme_object_category("Undefined", enzyme_module_forms)

        context = get_context(self)
        if context:
            context(partial(self.enzyme_module_forms.__isub__, ligands))
            context(partial(self.enzyme_module_forms.__isub__,
                            enzyme_module_forms))

    def remove_metabolites(self, metabolite_list, destructive=False):
        """Remove a list of metabolites and enzyme forms from the EnzymeModule.

        The species' initial conditions will also be removed from the model.

        The change is reverted upon exit when using the EnzymeModule as a
        context.

        Parameters
        ----------
        metabolite_list: list of MassMetabolites and EnzymeModuleForms
            A list of species to add to the EnzymeModule.
        destructive: bool, optional
            If False, the metabolites and enzyme forms are removed from all
            associated EnzymeModuleReactions.If True, also remove associated
            EnzymeModuleReactions from the EnzymeModule.

        Notes
        -----
        Extends from MassModel.remove_metabolites method.

        """
        # Ensure list is iterable.
        metabolite_list = ensure_iterable(metabolite_list)

        # Set ligands as "Undefined" to reset their categories.
        ligands = [met for met in metabolite_list
                   if not isinstance(met, EnzymeModuleForm)
                   and met in self.metabolites
                   and met in self.enzyme_module_ligands]
        self.set_enzyme_object_category("Undefined", ligands)

        # Remove metabolites from model using inherited method
        super(EnzymeModule, self).remove_metabolites(
            metabolite_list, destructive)

        # Remove ligands from the enzyme_module_ligands DictList.
        if ligands:
            self.enzyme_module_ligands -= ligands
        # Remove ligands from the "Undefined" category
        for met in ligands:
            self._enzyme_module_ligands_categorized["Undefined"].remove(met)

        # Get items that are EnzymeModuleForm objects, and check if the enzyme
        # forms already exists in the EnzymeModule, ignoring those that do not.
        enzyme_module_forms = [
            enzyme_module_form for enzyme_module_form in metabolite_list
            if isinstance(enzyme_module_form, EnzymeModuleForm)
            and enzyme_module_form in self.enzyme_module_forms]

        self.set_enzyme_object_category("Undefined", enzyme_module_forms)
        # Remove the enzyme forms to the model
        if enzyme_module_forms:
            self.enzyme_module_forms -= enzyme_module_forms

        # Remove enzyme forms from the "Undefined" category
        for enzyme_module_form in enzyme_module_forms:
            self._enzyme_module_forms_categorized.get(
                "Undefined").remove(enzyme_module_form)

        context = get_context(self)
        if context:
            context(partial(self.enzyme_module_ligands.__iadd__,
                            enzyme_module_forms))
            context(partial(self.enzyme_module_forms.__iadd__,
                            enzyme_module_forms))

    def add_reactions(self, reaction_list):
        """Add MassReactions to the EnzymeModule.

        MassReaction and EnzymeModuleReaction objects with identifiers
        identical to an existing reaction are ignored.

        The change is reverted upon exit when using the MassModel as a context

        Parameters
        ----------
        reaction_list: list of MassReactions and EnzymeModuleReactions
            A list of MassReactio and EnzymeModuleReaction objects.

        """
        # Ensure list is iterable.
        reaction_list = ensure_iterable(reaction_list)

        # Add reactions using inherited method
        super(EnzymeModule, self).add_reactions(reaction_list)

        # Get the enzyme module reactions by checking if an EnzymeModuleForm(s)
        # are involved, and check whether reaction exists,
        # ignoring those that do.
        enzyme_module_reactions = [
            r for r in reaction_list if isinstance(r, EnzymeModuleReaction)
            and r in self.reactions and r not in self.enzyme_module_reactions]

        # Add enzyme module reactions to the enzyme reaction DictList
        if enzyme_module_reactions:
            self.enzyme_module_reactions += enzyme_module_reactions
        # Add enzyme module reactions to the categorized attr as "Undefined".
        self.set_enzyme_object_category("Undefined", enzyme_module_reactions)

    def remove_reactions(self, reaction_list, remove_orphans=False):
        """Remove MassReactions from the EnzymeModule.

        The change is reverted upon exit when using the MassModel as a context.

        Parameters
        ----------
        reaction_list: list of MassReactions and EnzymeModuleReaction
            A list of MassReaction objects to be removed from the model.
        remove_orphans: bool, optional
            If True, will also remove orphaned genes, MassMetabolites, and
            EnzymeModuleForms from the EnzymeModule.

        """
        # Ensure list is iterable.
        reaction_list = ensure_iterable(reaction_list)
        # Get the enzyme module reactions by checking if EnzymeModuleForm(s)
        # are involved, then check whether reaction exists,
        # ignoring those that do not.
        enzyme_module_reactions = [
            r for r in reaction_list if isinstance(r, EnzymeModuleReaction)
            and r in self.reactions and r in self.enzyme_module_reactions]

        # Set enzyme module reactions as "Undefined" to reset their categories.
        self.set_enzyme_object_category("Undefined", enzyme_module_reactions)

        # Remove reactions using inherited method
        super(EnzymeModule, self).remove_reactions(reaction_list)

        # Remove enzyme module reactions from DictList
        if self.enzyme_module_reactions:
            self.enzyme_module_reactions -= enzyme_module_reactions

        # Remove enzyme module reactions from the "Undefined" category
        for reaction in enzyme_module_reactions:
            self._enzyme_module_reactions_categorized.get(
                "Undefined").remove(reaction)

    def repair(self, rebuild_index=True, rebuild_relationships=True):
        """Update all indicies and pointers in the model.

        In addition to updating indicies and pointers, the
        enzyme_module_reactions attribute will be updated to ensure it contains
        all existing reactions involving EnzymeModuleForms.

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
        super(EnzymeModule, self).repair(rebuild_index, rebuild_relationships)
        # Repair enzyme_module_reactions DictList
        self._get_current_enzyme_module_objs("reactions", update_enzyme=True)
        self._update_object_pointers()
        # Rebuild DictList indices
        if rebuild_index:
            for attr in ["enzyme_module_ligands", "enzyme_module_forms",
                         "enzyme_module_reactions"]:
                getattr(self, attr)._generate_index()
                for value in itervalues(getattr(self, attr + "_categorized")):
                    value._generate_index()

        for enzyme_module_form in self.enzyme_module_forms:
            enzyme_module_form._repair_bound_obj_pointers()

    # Overridden methods
    def copy(self):
        """Create a partial "deepcopy" of the EnzymeModule.

        All of the MassMetabolite, MassReaction, Gene, EnzymeModuleForm,
        EnzymeModuleReaction, and EnzymeModuleDict objects, the boundary
        conditions, custom_rates, custom_parameters, and the stoichiometric
        matrix are created anew, but in a faster fashion than deepcopy.

        Notes
        -----
        Overrides MassModel.copy method.

        """
        # Define a new model
        new_model = self.__class__()
        # Define items that will not be copied by their references
        do_not_copy_by_ref = [
            "metabolites", "reactions", "genes", "enzyme_modules", "_S",
            "enzyme_module_ligands", "enzyme_module_forms",
            "enzyme_module_reactions", "_enzyme_module_ligands_categorized",
            "_enzyme_module_forms_categorized",
            "_enzyme_module_reactions_categorized", "boundary_conditions",
            "custom_rates", "custom_parameters", "notes", "annotation",
            "modules"]
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
        # Copy any existing enzyme_modules
        new_model = self._copy_model_enzyme_modules(new_model)
        # Create the new stoichiometric matrix for the model.
        new_model._S = self._mk_stoich_matrix(matrix_type=self._matrix_type,
                                              dtype=self._dtype,
                                              update_model=True)
        # Add the newly copied objects to their appropriate DictLists
        # in the enzyme_module_ligands, enzyme_module_forms and
        # enzyme_module_reactions attributes
        for metabolite in new_model.metabolites:
            if isinstance(metabolite, EnzymeModuleForm):
                new_model.enzyme_module_forms.append(metabolite)
            else:
                new_model.enzyme_module_ligands.append(metabolite)
        new_model._get_current_enzyme_module_objs(attr="reactions",
                                                  update_enzyme=True)
        new_model._update_object_pointers()

        # Doesn't make sense to retain the context of a copied model so
        # assign a new empty context
        new_model._contexts = []
        return new_model

    def merge(self, right, prefix_existing=None, inplace=False,
              new_model_id=None):
        """Merge two MassModels into one MassModel with the objects from both.

        The reactions, metabolites, genes, enzyme modules, initial conditions,
        fixed concentrations, custom rate laws, rate parameters, compartments,
        units, notes, and annotations from right model are also copied to left
        model. However, note that in cases where identifiers for objects are
        identical or a dict item has an identical key(s), priority will be
        given to what already exists in the left model.

        Parameters
        ----------
        right: MassModel
            The MassModel to merge into the left model.
        prefix_existing: str, optional
            If provided, the string is used to prefix the reaction identifier
            of a reaction in the second model if that reaction already exists
            within the left model. Will also apply prefix to enzyme identifiers
            of an enzyme in the second model.
        inplace : bool
            Add reactions from right directly to left model object.
            Otherwise, create a new model leaving the left model untouched.
        new_model_id: str, optional
            If provided, the string is used as the identifier for the merged
            model. If None and inplace is True, the model ID of the first model
            will be used. If None and inplace is False, a new combined ID
            will be used for the new MassModel object.

        Returns
        -------
        new_model: MassModel
            A new MassModel object or self representing the merged model.

        Notes
        -----
        When merging an EnzymeModule into a MassModel, the EnzymeModule is
            converted to an EnzymeDict and stored in a DictList accessible
            via MassModel.enzyme_modules.
        If an EnzymeModule already exists in the model, it will be replaced.
        When merging an EnzymeModule with another EnzymeModule, a new
            EnzymeModule object will be returned, where the EnzymeModule is
            a copy of the 'left' model (self)  and the 'right' model is
            contained within.

        Extends from MassModel.merge

        """
        if not isinstance(right, EnzymeModule):
            # Always merge the EnzymeModule into the MassModel
            new_model = right.merge(self, prefix_existing=prefix_existing,
                                    inplace=inplace, new_model_id=new_model_id)
        else:
            # Merge the two EnzymeModules together if right is an EnzymeModule
            new_model = MassModel(self).merge(
                right, prefix_existing=prefix_existing, inplace=inplace,
                new_model_id=new_model_id)
            if inplace:
                new_model = self
                enzyme_modules_attrs_to_add = [right]
            else:
                # Reclassify as an EnzymeModule
                new_model = EnzymeModule(new_model)
                enzyme_modules_attrs_to_add = [self, right]
                # Set EnzymeModule attributes in the new model
                for attr in ["subsystem", "_enzyme_concentration_total",
                             "_enzyme_net_flux", "enzyme_net_flux_equation"]:
                    setattr(new_model, attr, getattr(self, attr))

            # Fix enzyme module ligands, forms, and reactions
            for attr in ["ligands", "forms", "reactions"]:
                # Update DictList attributes
                new_model._get_current_enzyme_module_objs(attr=attr,
                                                          update_enzyme=True)
                # Update categorized dict attributes
                attr = "enzyme_module_" + attr + "_categorized"
                new_categorized_dict = {}
                for enzyme_module in enzyme_modules_attrs_to_add:
                    old_categorized_dict = getattr(enzyme_module, attr)
                    # Add the old categorized dict to the new one to set
                    new_categorized_dict.update({
                        enzyme_module.id + " " + cat: values
                        for cat, values in iteritems(old_categorized_dict)})
                # Add new categorized dict attributes
                setattr(new_model, attr, new_categorized_dict)

        return new_model

    # Internal
    def _update_object_pointers(self):
        """Update objects in the attributes to be the objects from the model.

        Warnings
        --------
        This method is intended for internal use only.

        """
        for attr in ["enzyme_module_ligands", "enzyme_module_forms",
                     "enzyme_module_reactions"]:
            model_dictlist = {
                "enzyme_module_ligands": self.metabolites,
                "enzyme_module_forms": self.metabolites,
                "enzyme_module_reactions": self.reactions}.get(attr)
            setattr(self, attr,
                    _mk_new_dictlist(model_dictlist, getattr(self, attr)))
            attr += "_categorized"
            setattr(self, attr, {
                key: _mk_new_dictlist(model_dictlist, old_dictlist)
                for key, old_dictlist in iteritems(getattr(self, attr))})

    def _get_current_enzyme_module_objs(self, attr, update_enzyme=True):
        """Get the enzyme module objects for 'attr' that exist in the model.

        Parameters
        ----------
        attr: str {'ligands', 'forms', 'reactions'}
            A string representing which attribute to update.
        update_enzyme: bool, optional
            If True, update the enzyme_module_reactions attribute of the
            EnzymeModule.

        Warnings
        --------
        This method is intended for internal use only.

        """
        if attr not in {"ligands", "forms", "reactions"}:
            raise ValueError("Unrecognized attribute: '{0}'.".format(attr))
        item_list = []
        if attr == "ligands":
            item_list += list(filter(
                lambda x: not isinstance(x, EnzymeModuleForm),
                self.metabolites))
        if attr == "forms":
            item_list += list(filter(
                lambda x: isinstance(x, EnzymeModuleForm), self.metabolites))
        if attr == "reactions":
            for enzyme_module_form in self.enzyme_module_forms:
                item_list += [
                    rxn for rxn in list(enzyme_module_form.reactions)
                    if rxn not in item_list]

        item_list = DictList(item_list)
        item_list.sort()
        if set(item_list) ^ set(getattr(self, "enzyme_module_" + attr)) \
           and update_enzyme:
            setattr(self, "enzyme_module_" + attr, item_list)

        return item_list

    def _make_summation_expr(self, items, object_type):
        """Create a sympy expression of the summation of the given items.

        Warnings
        --------
        This method is intended for internal use only.

        """
        # Get appropriate list based on object type
        item_list = {
            EnzymeModuleForm: self.enzyme_module_forms,
            EnzymeModuleReaction: self.enzyme_module_reactions
        }.get(object_type)

        # Ensure object_type are appropriate type and exist in model
        items = list(filter(
            lambda x: isinstance(x, object_type) and x in item_list, items))

        return sum([x.flux_symbol if isinstance(x, EnzymeModuleReaction)
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
        if object_type == EnzymeModuleReaction:
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
                EnzymeModuleForm: ("forms", "initial conditions"),
                EnzymeModuleReaction: ("reactions", "steady state fluxes")
            }.get(object_type)
            warn("Not all enzyme {0} have {1} defined in model, will return a "
                 "sympy expression".format(*message_strs))

        return expr

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
        """Set the a categorized attribute dictionary.

        Applies to categorized attributes for enzyme_module_ligands,
        enzyme_module_forms, and enzyme_module_reactions

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
            dictlist_dict = {
                "_enzyme_module_ligands_categorized": "enzyme_module_ligands",
                "_enzyme_module_forms_categorized": "enzyme_module_forms",
                "_enzyme_module_reactions_categorized": "reactions"}
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

    def _convert_self_into_enzyme_module_dict(self):
        """Convert self into an EnzymeModuleDict.

        Primarily used for merging an EnzymeModule into a MassModel while
        retaining defined EnzymeModule attributes. Also used for checking if
        subclass is an EnzymeModule to avoid circular imports

        Warnings
        --------
        This method is intended for internal use only.

        """
        return EnzymeModuleDict(self)

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
        # Turn EnzymeModule into an EnzymeModuleDict
        # to store in MassModel.enzyme_modules
        enzyme_dict = self._convert_self_into_enzyme_module_dict()
        # Update EnzymeModuleDict with attributes
        enzyme_dict._update_object_pointers(model)

        if enzyme_dict.id in model.enzyme_modules:
            model.enzyme_modules.remove(enzyme_dict.id)
        model.enzyme_modules.append(enzyme_dict)
        enzyme_dict.model = model

        return model

    def _repr_html_(self):
        """HTML representation of the overview for the EnzymeModule.

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
                    <td>{num_enzyme_module_ligands}</td>
                </tr><tr>
                    <td><strong>Number of EnzymeModuleForms</strong></td>
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
                   num_enzyme_module_ligands=len(self.enzyme_module_ligands),
                   num_enz_forms=len(self.enzyme_module_forms),
                   num_enz_reactions=len(self.enzyme_module_reactions),
                   enz_conc=self.enzyme_concentration_total,
                   enz_flux=self.enzyme_net_flux,
                   num_ic=len(self.initial_conditions),
                   num_parameters=len(self._get_all_parameters()),
                   num_irreversible=len(self.irreversible_reactions),
                   num_custom_rates=len(self.custom_rates),
                   compartments=", ".join(v if v else k for k, v in
                                          iteritems(self.compartments)),
                   units=", ".join([u.id for u in self.units]))
