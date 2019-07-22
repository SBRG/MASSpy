# -*- coding: utf-8 -*-
r"""EnzymeModule is a class for handling reconstructions of enzymes.

The :class:`EnzymeModule` is a reconstruction an enzyme's mechanism and
behavior in a context of a larger system. To aid in the reconstruction process,
the :class:`EnzymeModule` contains various methods to build and add associated
:class:`~.EnzymeModuleSpecies` and :class:`~.EnzymeModuleReaction`\ s
of the enzyme module (see :meth:`~EnzymeModule.make_enzyme_module_species` and
:meth:`~~EnzymeModule.make_enzyme_module_reaction` methods, respecitvely).

Given the wide variety of enzymes and the various interactions it can have
with ligands (e.g. catalyzation, inhibition, activation, etc.), the enzyme
module has the following "categorized dict" attributes:

    * :attr:`EnzymeModule.enzyme_module_ligands_categorized`
    * :attr:`EnzymeModule.enzyme_module_species_categorized`
    * :attr:`EnzymeModule.enzyme_module_reactions_categorized`

These "categorized dict" attributes allow the user to define categories and
place various objects into one or more categories in the respective
"categorized dict" attribute (ligands a.k.a.
:class:`~.MassMetabolite`\ s, :class:`~.EnzymeModuleSpecies`, and
:class:`~.EnzymeModuleReaction`\ s). Utilizing categories with these attributes
can help with the management of complex enzymes, and are preserved upon
merging an :class:`EnzymeModule` into a :class:`~.MassModel`.

Because the :class:`EnzymeModule` is a subclass of the :class:`~.MassModel`, it
can be merged with a :class:`~.MassModel` representing the larger network in
which the enzyme is a part of.  For best results, an :class:`EnzymeModule`
should always be merged into the model as follows::

    model.merge(enzyme_module, inplace=False)
    # OR
    new_model = model.merge(enzyme_module, inplace=True)

Once merged, the :class:`~.EnzymeModuleSpecies` and
:class:`~.EnzymeModuleReaction`\ s of the :class:`EnzymeModule` are treated
like any other :class:`~.MassMetabolite` or :class:`~.MassReaction`.

Therefore, to prevent the loss of the enzyme specific information that was
stored in the :class:`EnzymeModule`, the enzyme module is converted into an
ordered dictionary known as an :class:`~.EnzymeModuleDict`, which contains
most of the enzyme specific attribute information. Note that all enzyme
specific attribute names start with either ``"enzyme"`` or ``"enzyme_module"``.

During the model merging process, the :class:`~.EnzymeModuleDict` is created,
then stored in the :attr:`.MassModel.enzyme_modules` attribute for access at
a later time. See the :mod:`~.enzyme_module_dict` documentation for more
information about the :class:`~.EnzymeModuleDict`.
"""
import re
from collections import defaultdict
from copy import deepcopy
from functools import partial
from warnings import warn

from cobra.core.dictlist import DictList
from cobra.util.context import get_context

import numpy as np

from six import integer_types, iteritems, iterkeys, itervalues, string_types

import sympy as sym

from mass.core.mass_model import MassModel
from mass.enzyme_modules.enzyme_module_dict import EnzymeModuleDict
from mass.enzyme_modules.enzyme_module_reaction import EnzymeModuleReaction
from mass.enzyme_modules.enzyme_module_species import EnzymeModuleSpecies
from mass.util.expressions import _mk_met_func, strip_time
from mass.util.util import (
    _mk_new_dictlist, ensure_iterable, ensure_non_negative_value)

_AUTOMATIC_RE = re.compile("^Automatic$")
_UNDEFINED_RE = re.compile("^Undefined$")
_EQUATION_RE = re.compile("^Equation$")


class EnzymeModule(MassModel):
    r"""Class representation of an enzyme module reconstruction.

    Parameters
    ----------
    id_or_model : str, MassModel, EnzymeModule
        A string identifier to associate with the :class:`EnzymeModule`, or an
        existing model object. If an existing model object is provided, a new
        :class:`EnzymeModule` object is instantiated with the same properties
        as the original model.
    name : str
        A human readable name for the model.
    subsystem : str
        The subsystem in which the enzyme module is a part of.
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
    enzyme_module_ligands : ~cobra.core.dictlist.DictList
        A :class:`~cobra.core.dictlist.DictList` where the keys are the
        metabolite identifiers and the values are the associated
        :class:`~.MassMetabolite`\ s.
    enzyme_module_species : ~cobra.core.dictlist.DictList
        A :class:`~cobra.core.dictlist.DictList` where the keys are the
        :class:`~.EnzymeModuleSpecies` identifiers and the values are the
        associated :class:`~.EnzymeModuleSpecies`.
    enzyme_module_reactions : ~cobra.core.dictlist.DictList
        A :class:`~cobra.core.dictlist.DictList` where keys are the
        reaction identifiers and the values are the associated
        :class:`~.EnzymeModuleReaction`\ s.
    enzyme_module_ligands_categorized : dict
        A dictionary of user-categorized ligands where keys are categories and
        values are :class:`~cobra.core.dictlist.DictList` of corresponding
        :class:`~.MassMetabolite`\ s.
    enzyme_module_species_categorized : dict
        A dictionary of user-categorized :class:`~.EnzymeModuleSpecies` where
        keys are categories and values are
        :class:`~cobra.core.dictlist.DictList` of the corresponding
        :class:`~.EnzymeModuleSpecies`.
    enzyme_module_reactions_categorized : dict
        A dictionary of user-categorized enzymatic binding reactions where
        keys are categories and values are
        :class:`~cobra.core.dictlist.DictList` of the corresponding
        :class:`~.EnzymeModuleReaction`\ s.
    enzyme_concentration_total : float
        A non-negative value representing the total enzyme concentration.
    enzyme_net_flux : float
        The value of the net flux through the enzyme.
    enzyme_concentration_total_equation : ~sympy.core.basic.Basic
        A :mod:`sympy` expression representing the net reaction rate equation
        for the enzyme represented by the module.
    enzyme_net_flux_equation :  ~sympy.core.basic.Basic
        A :mod:`sympy` expression representing the net reaction rate equation
        for the enzyme represented by the module.

    """

    # pylint: disable=too-many-instance-attributes
    def __init__(self, id_or_model=None, name=None, subsystem="",
                 matrix_type="dense", dtype=np.float64):
        """Initialize the EnzymeModule."""
        super(EnzymeModule, self).__init__(
            id_or_model=id_or_model, name=name, matrix_type=matrix_type,
            dtype=dtype)

        self.subsystem = subsystem
        # Initialize DictLists for enzyme ligands, species, and reactions
        self.enzyme_module_ligands = DictList()
        self.enzyme_module_species = DictList()
        self.enzyme_module_reactions = DictList()

        # Initialize a dict of DictLists for storing categorized objects
        self._enzyme_module_ligands_categorized = {"Undefined": DictList()}
        self._enzyme_module_species_categorized = {"Undefined": DictList()}
        self._enzyme_module_reactions_categorized = {"Undefined": DictList()}

        # Initialize EnzymeModule attributes
        self._enzyme_concentration_total = None
        self._enzyme_net_flux = None
        self.enzyme_net_flux_equation = None

    @property
    def enzyme_total_symbol(self):
        """Get the symbol for the total enzyme concentration."""
        if self.id is None:
            return sym.Symbol("Enzyme_Total")

        return sym.Symbol(self.id + "_Total")

    @property
    def enzyme_flux_symbol(self):
        """Get the symbol for the net flux through the enzyme."""
        if self.id is None:
            return sym.Symbol("v")

        return sym.Symbol("v_" + self.id)

    @property
    def enzyme_concentration_total(self):
        """Get or set the total concentration value.

        Notes
        -----
        The total concentration of the enzyme cannot be negative.

        Parameters
        ----------
        concentration : float
            A non-negative number for the concentration of the enzyme.

        Raises
        ------
        ValueError
            Occurs when trying to set a negative value.

        """
        return getattr(self, "_enzyme_concentration_total")

    @enzyme_concentration_total.setter
    def enzyme_concentration_total(self, concentration):
        """Set the expected total enzyme concentration."""
        concentration = ensure_non_negative_value(concentration)
        setattr(self, "_enzyme_concentration_total", concentration)

    @property
    def enzyme_net_flux(self):
        """Get or set the net flux through the enzyme.

        Parameters
        ----------
        flux_value : float
            The value of the net flux through the enzyme.

        """
        return getattr(self, "_enzyme_net_flux")

    @enzyme_net_flux.setter
    def enzyme_net_flux(self, flux_value):
        """Set the expected net flux through the enzyme."""
        if not isinstance(flux_value, (integer_types, float)) and \
           flux_value is not None:
            raise TypeError("Must be an int or float")

        setattr(self, "_enzyme_net_flux", flux_value)

    @property
    def enzyme_concentration_total_equation(self):
        r"""Return the total concentration equation.

        Notes
        -----
        * Will sum the :class:`~.EnzymeModuleSpecies` that have their
          :attr:`~.EnzymeModuleSpecies.enzyme_module_id` match the
          :attr:`EnzymeModule.id`.
        * If no :class:`~.EnzymeModuleSpecies` are found to have an
          :attr:`~.EnzymeModuleSpecies.enzyme_module_id` that matches the
          :attr:`EnzymeModule.id`, all :class:`~.EnzymeModuleSpecies`
          in the the model will be used.

        Returns
        -------
        ~sympy.core.relational.Equality
            A :mod:`sympy` equation with the left-hand side containing the
            :attr:`EnzymeModule.enzyme_total_symbol` and the right-hand side
            containing the sum of the :class:`~.EnzymeModuleSpecies`.

        """
        if not self.enzyme_module_species:
            return None

        # First try only using enzyme module species that reference this module
        enzyme_module_species = [
            species for species in self.enzyme_module_species
            if species.enzyme_module_id == self.id]
        # If None found, use all enzyme module species in the EnzymeModule.
        if not enzyme_module_species:
            enzyme_module_species = self.enzyme_module_species

        return sym.Eq(self.enzyme_total_symbol,
                      self.sum_enzyme_module_species_concentrations(
                          enzyme_module_species, use_values=False))

    @property
    def enzyme_net_flux_equation(self):
        """Get or set the net flux rate equation of the enzyme.

        Parameters
        ----------
        equation_rhs : str, ~sympy.core.basic.Basic
            Either a string representing the right-hand side of the equation
            that will be sympified via the :func:`~sympy.core.sympify.sympify`
            function., or a :mod:`sympy` expression representing the
            right-hand side of the equation.

        Returns
        -------
        ~sympy.core.relational.Equality
            A :mod:`sympy` equation with the left-hand side containing the
            :attr:`EnzymeModule.enzyme_flux_symbol` and the right-hand side
            containing the equation repredenting the net flux through the
            enzyme.

        """
        value = getattr(self, "_enzyme_net_flux_equation", None)
        if value is not None:
            value = sym.Eq(self.enzyme_flux_symbol, value)
        return value

    @enzyme_net_flux_equation.setter
    def enzyme_net_flux_equation(self, equation_rhs):
        """Set the net rate equation of the enzyme."""
        if equation_rhs is not None:
            if not isinstance(equation_rhs, (sym.Basic, string_types)):
                raise TypeError("equation_rhs must be a sympy expression.")
            if isinstance(equation_rhs, string_types):
                equation_rhs = sym.sympify(equation_rhs)
            elif hasattr(equation_rhs, "lhs") and hasattr(equation_rhs, "rhs"):
                if equation_rhs.lhs == self.enzyme_flux_symbol:
                    equation_rhs = equation_rhs.rhs
                if equation_rhs.rhs == self.enzyme_flux_symbol:
                    equation_rhs = equation_rhs.lhs
            else:
                pass
        setattr(self, "_enzyme_net_flux_equation", equation_rhs)

    @property
    def enzyme_module_ligands_categorized(self):
        r"""Get or set categories for ligands.

        Notes
        -----
        * A ligand must already exist in the :class:`EnzymeModule` as a
          :class:`~.MassMetabolite` in order to set its category.
        * If categories already exists, their existing contents are replaced.
        * Categories with empty :class:`~cobra.core.dictlist.DictList`\ s
          are removed.

        Parameters
        ----------
        value : dict
            A dictionary where keys are the strings representing the
            categories of the ligands, and values are ``list``\ s containing
            the corresponding :class:`~.MassMetabolite`\ s.
            An empty ``dict`` will cause a reset, placing all ligands into an
            ``"Undefined"`` category.

        See Also
        --------
        set_enzyme_object_category
            Method to categorize ligands as :class:`~.MassMetabolite`\ s
            without replacing existing category contents.

        """
        return self._remove_empty_categories(
            "_enzyme_module_ligands_categorized")

    @enzyme_module_ligands_categorized.setter
    def enzyme_module_ligands_categorized(self, value):
        """Set categories(s) for ligands using a ``dict``."""
        if not isinstance(value, dict):
            raise TypeError("value must be a dict")

        self._set_category_attribute_dict(
            value, "_enzyme_module_ligands_categorized")

    @property
    def enzyme_module_species_categorized(self):
        r"""Get or set categories for :class:`~.EnzymeModuleSpecies`.

        Notes
        -----
        * An enzymatic species must already exist in the :class:`EnzymeModule`
          as a :class:`~.EnzymeModuleSpecies` in order to set its category.
        * If categories already exists, their existing contents are replaced.
        * Categories with empty :class:`~cobra.core.dictlist.DictList`\ s
          are removed.

        Parameters
        ----------
        value : dict
            A dictionary where keys are the strings representing the
            categories of the enzymatic species, and values are ``list``\ s
            containing the corresponding :class:`~.EnzymeModuleSpecies`\ s.
            An empty ``dict`` will cause a reset, placing all enzymatic species
            into an ``"Undefined"`` category.

        See Also
        --------
        set_enzyme_object_category
            Method to categorize :class:`~.EnzymeModuleSpecies`\ s without
            replacing existing category contents.

        """
        return self._remove_empty_categories(
            "_enzyme_module_species_categorized")

    @enzyme_module_species_categorized.setter
    def enzyme_module_species_categorized(self, value):
        """Set categories(s) for enzymatic species using a ``dict``."""
        if not isinstance(value, dict):
            raise TypeError("value must be a dict")

        self._set_category_attribute_dict(
            value, "_enzyme_module_species_categorized")

    @property
    def enzyme_module_reactions_categorized(self):
        r"""Get or set categories for :class:`~.EnzymeModuleReaction`\ s.

        Notes
        -----
        * An enzymatic binding reaction must already exist in the
          :class:`EnzymeModule` as a :class:`~.EnzymeModuleReaction` in order
          to set its category.
        * If categories already exists, their existing contents are replaced.
        * Categories with empty :class:`~cobra.core.dictlist.DictList`\ s
          are removed.

        Parameters
        ----------
        value : dict
            A dictionary where keys are the strings representing the
            categories of the enzymatic binding reactions, and values are
            ``list``\ s containing the corresponding
            :class:`~.EnzymeModuleReaction`\ s. An empty ``dict`` will cause
            a reset, placing all enzymatic binding reactions into
            an ``"Undefined"`` category.

        See Also
        --------
        set_enzyme_object_category
            Method to categorize :class:`~.EnzymeModuleReaction`\ s without
            replacing existing category contents.

        """
        return self._remove_empty_categories(
            "_enzyme_module_reactions_categorized")

    @enzyme_module_reactions_categorized.setter
    def enzyme_module_reactions_categorized(self, value):
        """Set categories(s) for enzyme module reactions using a ``dict``."""
        if not isinstance(value, dict):
            raise TypeError("value must be a dict")

        self._set_category_attribute_dict(
            value, "_enzyme_module_reactions_categorized")

    def set_enzyme_object_category(self, category, object_list):
        r"""Add a list of objects to a new or existing category.

        Notes
        -----
        * If a category already exists, the objects will be added to the
          existing list. If an object is categorized as ``"Undefined"``, it
          will be removed from all other existing categories.
        * The objects to be categorized must already exist in the
          :class:`EnzymeModule`.

        Parameters
        ----------
        category : str
            A string representing the category for the list of objects to be
            categorized.
        object_list : list
            A ``list`` containing the objects or their identifiers to be
            categorized. The ``list`` must contain ONLY one of the following
            :mod:`mass` object types:

                * :class:`~.MassMetabolite`\ s representing enzyme ligands.
                * :class:`~.EnzymeModuleSpecies`\ s representing enzymatic
                  species.
                * :class:`~.EnzymeModuleReaction`\ s representing enzymatic
                  binding reactions ligands.

        """
        if not isinstance(category, string_types):
            raise TypeError("category must be a str")
        object_list = ensure_iterable(object_list)
        sep_objs = defaultdict(list)
        for obj in object_list:
            if isinstance(obj, EnzymeModuleSpecies.__base__)\
               and not isinstance(obj, EnzymeModuleSpecies):
                sep_objs["_enzyme_module_ligands_categorized"].append(obj)
            if isinstance(obj, EnzymeModuleSpecies):
                sep_objs["_enzyme_module_species_categorized"].append(obj)
            if isinstance(obj, EnzymeModuleReaction):
                sep_objs["_enzyme_module_reactions_categorized"].append(obj)

        for key, object_list in iteritems(sep_objs):
            if category in getattr(self, key):
                object_list = set(
                    getattr(self, key)[category]).union(set(object_list))

            self.__class__.__dict__[key[1:]].fset(
                self, {category: object_list})

    def make_enzyme_module_species(self, id=None, name="automatic",
                                   categories="Undefined",
                                   bound_catalytic=None, bound_effectors=None,
                                   compartment=None):
        r"""Create and add an :class:`~.EnzymeModuleSpecies` to the module.

        Notes
        -----
        * Upon creation, the :class:`~.EnzymeModuleSpecies` is added to the
          :class:`EnzymeModule`.
        * If :class:`~.MassMetabolite`\ s in the ``bound_catalytic`` or the
          ``bound_effectors`` arguments do not already exist in the
          :class:`EnzymeModule`, they will also be added.

        Parameters
        ----------
        id : str
            A string identifier to associate with the enzymatic species.
        name : str
            Either a human readable name for the enzyme module species, or the
            string ``"automatic"``. If set to ``"automatic"``, a name will be
            generated based on the identifier of the enzyme module species and
            its bound ligands.
        categories : str or list
            A string representing the category, or a list of strings
            containing several categories for the enzyme module species.
        bound_catalytic : dict
            A ``dict`` representing the ligands bound to the enzyme's active
            site(s), with :class:`~.MassMetabolite`\ s or their identifiers as
            keys and the number bound as values.
        bound_effectors : dict
            A ``dict`` representing the ligands bound to the enzyme's
            regulatory site(s), with :class:`~.MassMetabolite`\ s or their
            identifiers as keys and the number bound as values.
        compartment : str
            The compartment where the enzyme module species is located.

        Returns
        -------
        EnzymeModuleSpecies
            The newly created :class:`~.EnzymeModuleSpecies`.

        See Also
        --------
        .EnzymeModuleSpecies.generate_enzyme_module_species_name
            Automatic generation of the ``name`` for an
            :class:`~.EnzymeModuleSpecies`.

        """
        # Ensure metabolites for EnzymeModuleSpecies exist in the EnzymeModule.
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
        # Make EnzymeModuleSpecies object
        enzyme_module_species = EnzymeModuleSpecies(
            id_or_specie=id, name=name, enzyme_module_id=self.id,
            bound_catalytic=bound_catalytic, bound_effectors=bound_effectors,
            compartment=compartment)
        # Generate name for the EnzymeModuleSpecies if name set to "Automatic"
        if _AUTOMATIC_RE.match(name):
            enzyme_module_species.generate_enzyme_module_species_name(True)

        # Add the enzyme species to the module and place in respective
        # categories
        self.add_metabolites(enzyme_module_species)
        categories = ensure_iterable(categories)
        for category in categories:
            self.set_enzyme_object_category(category, enzyme_module_species)

        return enzyme_module_species

    def make_enzyme_module_reaction(self, id=None, name="", subsystem=None,
                                    reversible=True, categories="Undefined",
                                    metabolites_to_add=None):
        r"""Create and add an :class:`~.EnzymeModuleReaction` to the module.

        Notes
        -----
        * When adding metabolites, a final coefficient of < 0 implies a
          reactant and a final coefficient of > 0 implies a product.

        Parameters
        ----------
        id : str
            The identifier associated with the enzyme module reaction.
        name : str
            A human readable name for the enzyme module reaction. If name
            is set to match ``"automatic"``, a name will be generated based
            on the :class:`~.EnzymeModuleSpecies` and their bound ligands.
        subsystem : str
            The subsystem where the reaction is meant to occur.
        reversible : bool
            The kinetic reversibility of the reaction. Irreversible reactions
            have an equilibrium constant and a reverse rate constant as set
            in the :attr:`~.MassBaseConfiguration.irreversible_Keq` and
            :attr:`~.MassBaseConfiguration.irreversible_kr` attributes of the
            :class:`~.MassConfiguration`. Default is ``True``.
        categories : str or list
            A string representing the category, or a list of strings
            containing several categories for the enzyme module reactions.
        metabolites_to_add : dict
            A ``dict`` with :class:`~.MassMetabolite`\ s and
            :class:`~.EnzymeModuleSpecies` or their identifiers as keys and
            stoichiometric coefficients as values. If keys are string
            identifiers then the :class:`~.MassMetabolite`\ s and
            :class:`~.EnzymeModuleSpecies` must already be a part of model.

        Returns
        -------
        EnzymeModuleReaction
            The newly created :class:`~.EnzymeModuleReaction`.

        See Also
        --------
        .EnzymeModuleReaction.generate_enzyme_module_reaction_name
            Automatic generation of the ``name`` for an
            :class:`~.EnzymeModuleReaction`.


        """
        # Make EnzymeModuleReaction object
        new_reaction = EnzymeModuleReaction(
            id_or_reaction=id, name=name, subsystem=subsystem,
            reversible=reversible, enzyme_module_id=self.id)

        # Add metabolites
        if metabolites_to_add:
            try:
                metabolites_to_add = dict(
                    (met, c) if isinstance(met, EnzymeModuleSpecies.__base__)
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
                              rtype=1, enzyme_prefix=False):
        r"""Unify rate law parameters for a list of enzyme module reactions.

        After unification, the new parameters and rate laws are placed into the
        ``custom_parameters`` and ``custom_rates`` attributes, repsectively.

        Parameters
        ----------
        reaction_list : list
            A ``list`` containing :class:`~.EnzymeModuleReaction`\ s or
            their string identifiers. :class:`~.EnzymeModuleReaction`\ s must
            already exist in the :class:`~.EnzymeModule`.
        new_parameter_id : str
            The new parameter ID to use for the unified reaction parameters.
            The forward rate, reverse rate, and/or equilibrium constants in
            the current rate law will have the reaction ID component replaced
            with the ``new_parameter_id`` in the parameter ID.
        rtype : int
            The type of rate law to display. Must be 1, 2, or 3.

                * Type 1 will utilize the :attr:`forward_rate_constant` and the
                  :attr:`equilibrium_constant`.
                * Type 2 will utilize the :attr:`forward_rate_constant` and the
                  :attr:`reverse_rate_constant`.
                * Type 3 will utilize the :attr:`equilibrium_constant` and the
                  :attr:`reverse_rate_constant`.

            Default is ``1``.
        enzyme_prefix : bool
            If ``True``, add the :class:`EnzymeModule` ID as a prefix to the
            ``new_parameter_id`` before using the ``new_parameter_id`` in
            the rate parameter unification. Default is ``False``.

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
            # Create a string representation of the rate and replace the
            # reaction id portions of the parameters with new_parameter_id
            custom_rate = str(
                strip_time(reaction.get_mass_action_rate(rtype)))
            custom_rate = custom_rate.replace(reaction.id, new_parameter_id)
            self.add_custom_rate(reaction, custom_rate)

    def make_enzyme_net_flux_equation(self, enzyme_module_reactions,
                                      use_rates=False, update_enzyme=False):
        r"""Create an equation representing the net flux through the enzyme.

        The left side of the rate equation will always be the flux symbol of
        the enzyme, accessible via :attr:`enzyme_flux_symbol`.

        Parameters
        ----------
        enzyme_module_reactions : list
            A list containing the :class:`~.EnzymeModuleReaction`\ s or their
            identifiers to be summed for the equation.
        use_rates : bool
            If ``True``, then the rates of the provided reactions are used in
            creating the expression. Otherwise the arguments in the expression
            are left as the :attr:`.EnzymeModuleReaction.flux_symbol`\ s.
        update_enzyme : bool
            If ``True``, update the :attr:`enzyme_net_flux_equation` attribute
            and, if necessary, the :attr:`enzyme_module_reactions` attribute
            of the module in addition to returning the generated equation.
            Otherwise just return the net flux equation without making any
            updates. Default is ``False``.

        Returns
        -------
        ~sympy.core.basic.Basic
            A :mod:`sympy` expression of the net flux equation.

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

    def sum_enzyme_module_species_concentrations(self, enzyme_module_species,
                                                 use_values=False):
        """Sum the species concentrations for a list of enzyme module species.

        Parameters
        ----------
        enzyme_module_species : list
            A ``list`` containing the :class:`~.EnzymeModuleSpecies` or their
            identifiers to be summed. Species must already exist in the
            :class:`~.EnzymeModule`.
        use_values : bool
            If ``True``, then numerical values are substituted into the
            expression. Otherwise arguments in the expression are left as
            :mod:`sympy` symbols.

        Returns
        -------
        float or ~sympy.core.basic.Basic
            The sum of the concentrations for the given enzyme module species
            as a ``float`` if ``use_values`` is ``True`` and all values are
            present. Otherwise returns a ``sympy`` expression representing
            the sum of the given enzyme module species.

        """
        enzyme_module_species = ensure_iterable(enzyme_module_species)
        # Get concentration formula
        concentration = self._make_summation_expr(
            enzyme_module_species, EnzymeModuleSpecies)
        # Substitute numerical values into the formula
        if use_values:
            concentration = self._sub_values_into_expr(
                strip_time(concentration), EnzymeModuleSpecies)

        return concentration

    def sum_enzyme_module_reaction_fluxes(self, enzyme_module_reactions,
                                          use_values=False):
        """Sum the enzyme reaction steady state fluxes for a list of reactions.

        Parameters
        ----------
        enzyme_module_reactions : list
            A ``list`` a containing the :class:`~.EnzymeModuleReaction` or
            their identifiers to be summed. Reactions must already exist in
            the module and must be considered an enzyme module reaction.
        use_values : bool
            If ``True``, then numerical values are substituted into the
            expression. Otherwise arguments in the expression are left as
            :mod:`sympy` symbols.

        Returns
        -------
        float or ~sympy.core.basic.Basic
            The sum of the steady state fluxes for the given enzyme reaction
            as a ``float`` if ``use_values`` is ``True`` and all values are
            present. Otherwise returns a :mod:`sympy` expression representing
            the sum of the enzyme module reaction fluxes.

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
        difference between the :attr:`enzyme_concentration_total` value and
        the sum of all :class:`.EnzymeModuleSpecies` initial conditions
        in the model.

        Notes
        -----
        Positive values indicate the value in the
        :attr:`enzyme_concentration_total` attribute is greater than the value
        calculated using the expression from the
        :attr:`enzyme_concentration_total_equation` attribute.

        Parameters
        ----------
        use_values : bool
            If ``True``, then numerical values are substituted into the
            expression. Otherwise arguments in the expression are left as
            :mod:`sympy` symbols.

        Returns
        -------
        float or ~sympy.core.basic.Basic
            The error between the set :attr:`enzyme_concentration_total` and
            the sum of the :class:`~.EnzymeModuleSpecies` initial condition
            values in the model as a ``float`` if ``use_values`` is ``True``
            and all values are present. Otherwise returns a :mod:`sympy`
            expression representing the error.


        """
        if self.enzyme_concentration_total_equation is None:
            warn("No enzyme total concentration equation found. Ensure that "
                 "the model contains EnzymeModuleSpecies and an equation for "
                 "the enzyme_concentration_total_equation attribute returns an"
                 " expression")
            return None

        # Make error expression
        error = (self.enzyme_concentration_total_equation.lhs
                 - self.enzyme_concentration_total_equation.rhs)
        # Try to substitute values into equations
        if use_values:
            error = self._sub_values_into_expr(
                strip_time(error), EnzymeModuleSpecies, {
                    self.enzyme_total_symbol: self.enzyme_concentration_total})

        return error

    def enzyme_net_flux_error(self, use_values=False):
        """Return the error for the net flux through the enzyme.

        The error of the enzyme net flux is defined to be the difference
        between the :attr:`enzyme_net_flux` value and the
        calculated value for the :attr:`enzyme_net_flux_equation`.

        Notes
        -----
        Positive values indicate the value in :attr:`enzyme_net_flux`
        attribute is greater than the value calculated using the expression
        from the :attr:`enzyme_net_flux_equation` attribute.

        Parameters
        ----------
        use_values : bool
            If ``True``, then numerical values are substituted into the
            expression. Otherwise arguments in the expression are left as
            :mod:`sympy` symbols.

        Returns
        -------
        float or ~sympy.core.basic.Basic
            The error between the set :attr:`enzyme_net_flux` and the
            calculated value for the :attr:`enzyme_net_flux_equation`
            attribute as a ``float`` if ``use_values`` is ``True`` and all
            values are present in the model. Otherwise returns
            a :mod:`sympy` expression representing the error.

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

        Notes
        -----
        The string ``"Equation"`` can be passed to either ``top`` or
        ``bottom`` to utilize the equation in the corresponding attribute
        (i.e. :attr:`enzyme_concentration_total_equation` for ``'species'``
        and :attr:`enzyme_net_flux_equation` for ``'reactions'``).

        Parameters
        ----------
        categorized_attr: str or dict
            Either a string representing which categorized ``dict`` attribute
            to use or the attribute ``dict`` itself to use in making the
            enzyme ratio expression. Use the string ``'species'`` for
            :attr:`enzyme_module_species_categorized`, or ``'reactions'`` for
            :attr:`enzyme_module_reactions_categorized`.
        top : str
            A string representing a category in the ``dict`` corresponding
            to the given categorized attribute. The summation expression of
            the objects in the corresponding list is used as the top
            (numerator) of the fraction to be made. Alternatively, the
            string ``"Equation"`` can be provided to utilize an
            equation attribute.
        bottom : str
            A string representing a category in the ``dict`` corresponding
            to the given categorized attribute. The summation expression of
            the objects in the corresponding list is used as the top
            (numerator) of the fraction to be made. Alternatively, the
            string ``"Equation"`` can be provided to utilize an
            equation attribute. Default is ``'Equation'``.
        use_values : bool
            If ``True``, then numerical values are substituted into the
            expression. Otherwise arguments in the expression are left as
            :mod:`sympy` symbols.

        Returns
        -------
        float or ~sympy.core.basic.Basic
            The fraction either calculated and returned as float if
            ``use_values`` is ``True`` and all values are present in the
            model, or a :mod:`sympy` expression representing the formula
            for the fraction.

        """
        # Check categorized_attr input, and get corresponding categorized dict
        if isinstance(categorized_attr, dict):
            if categorized_attr == self.enzyme_module_species_categorized:
                categorized_attr = "species"
            elif categorized_attr == self.enzyme_module_reactions_categorized:
                categorized_attr = "reactions"
            else:
                raise ValueError(
                    "Must be the dict accessible through '"
                    "EnzymeModule.enzyme_module_species_categorized' or "
                    "'EnzymeModule.enzyme_module_reactions_categorized'.")

        if categorized_attr.lower() in {"species", "reactions"}:
            categorized_attr = categorized_attr.lower()
            item_dict = self.__class__.__dict__[
                "enzyme_module_" + categorized_attr + "_categorized"].fget(
                    self)
        else:
            raise ValueError("Must be a string of the following: "
                             "{'species', 'reactions'}.")

        # Get object_type for summation expression
        object_type = {"species": EnzymeModuleSpecies,
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
                    "species": self.enzyme_concentration_total_equation,
                    "reactions": self.enzyme_net_flux_equation
                }[categorized_attr]
                if summation_expr is not None:
                    summation_expr = summation_expr.rhs
                else:
                    raise ValueError(
                        "No equation found for '{0}' attribute".format({
                            "species": "enzyme_concentration_total_equation",
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

    def add_metabolites(self, metabolite_list):
        r"""Add a list of metabolites and enzyme species to the module.

        Species are added as the category as ``"Undefined"``.

        The change is reverted upon exit when using the :class:`EnzymeModule`
        as a context.

        Notes
        -----
        Extends from :meth:`.MassModel.add_metabolites`.

        Parameters
        ----------
        metabolite_list : list
            A list of :class:`~.MassMetabolite`\ s and
            :class:`EnzymeModuleSpecies` to add to the :class:`EnzymeModule`.
        add_initial_conditons : bool
            If ``True``, the initial conditions associated with the species
            are also added to the model. Otherwise, the species are added
            without their initial conditions.

        """
        # Ensure list is iterable.
        metabolite_list = ensure_iterable(metabolite_list)

        # Add metabolites using inherited method
        super(EnzymeModule, self).add_metabolites(metabolite_list)

        # Get items that are not EnzymeModuleSpecies objects, and check if the
        # ligands already exists in the EnzymeModule, ignoring those that do.
        ligands = [met for met in metabolite_list
                   if not isinstance(met, EnzymeModuleSpecies)
                   and met in self.metabolites
                   and met not in self.enzyme_module_ligands]
        # Add ligands to the ligand DictList
        if ligands:
            self.enzyme_module_ligands += ligands
        # Add ligands to the ligand dict as "Undefined"
        self.set_enzyme_object_category("Undefined", ligands)

        # Get items that are EnzymeModuleSpecies objects, and check if the
        # enzyme species already exist in the EnzymeModule, ignoring those that
        # do.
        enzyme_module_species = [
            enzyme_module_species for enzyme_module_species in metabolite_list
            if isinstance(enzyme_module_species, EnzymeModuleSpecies)
            and enzyme_module_species not in self.enzyme_module_species]

        if enzyme_module_species:
            # Add the enzyme species to the model
            self.enzyme_module_species += enzyme_module_species
        # Add enzyme species to the enzyme_module_species dict as "Undefined".
        self.set_enzyme_object_category("Undefined", enzyme_module_species)

        context = get_context(self)
        if context:
            context(partial(self.enzyme_module_species.__isub__, ligands))
            context(partial(self.enzyme_module_species.__isub__,
                            enzyme_module_species))

    def remove_metabolites(self, metabolite_list, destructive=False):
        r"""Remove a list of metabolites and enzyme species from the module.

        The species' initial conditions will also be removed from the model.

        The change is reverted upon exit when using the :class:`EnzymeModule`
        as a context.

        Notes
        -----
        Extends from :meth:`.MassModel.remove_metabolites`.

        Parameters
        ----------
        metabolite_list : list
            A list of :class:`~.MassMetabolite`\ s and
            :class:`~.EnzymeModuleSpecies` to remove from the
            :class:`EnzymeModule`.
        destructive : bool
            If ``False``, the species are removed from all associated
            :class:`~.EnzymeModuleReaction`\ s . If ``True``, also remove
            associated :class:`~.EnzymeModuleReaction`\ s from the
            :class:`EnzymeModule`.

        """
        # Ensure list is iterable.
        metabolite_list = ensure_iterable(metabolite_list)

        # Set ligands as "Undefined" to reset their categories.
        ligands = [met for met in metabolite_list
                   if not isinstance(met, EnzymeModuleSpecies)
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

        # Get items that are EnzymeModuleSpecies objects, and check if the
        # enzyme species already exists in the EnzymeModule, ignoring those
        # that do not.
        enzyme_module_species = [
            enzyme_module_species for enzyme_module_species in metabolite_list
            if isinstance(enzyme_module_species, EnzymeModuleSpecies)
            and enzyme_module_species in self.enzyme_module_species]

        self.set_enzyme_object_category("Undefined", enzyme_module_species)
        # Remove the enzyme species to the model
        if enzyme_module_species:
            self.enzyme_module_species -= enzyme_module_species

        # Remove enzyme species from the "Undefined" category
        for enzyme_module_specie in enzyme_module_species:
            self._enzyme_module_species_categorized.get(
                "Undefined").remove(enzyme_module_specie)

        context = get_context(self)
        if context:
            context(partial(self.enzyme_module_ligands.__iadd__,
                            enzyme_module_species))
            context(partial(self.enzyme_module_species.__iadd__,
                            enzyme_module_species))

    def add_reactions(self, reaction_list):
        r"""Add a list of reactions to the :class:`EnzymeModule`.

        :class:`~.MassReaction`\ s and :class:`~.EnzymeModuleReaction`\ s
        with identifiers identical to an existing reaction are ignored.

        The change is reverted upon exit when using the :class:`EnzymeModule`
        as a context.

        Notes
        -----
        Extends from :meth:`.MassModel.add_reactions`.

        Parameters
        ----------
        reaction_list : list
            A list of :class:`~.MassReaction` and
            :class:`~.EnzymeModuleReaction` to add.

        """
        # Ensure list is iterable.
        reaction_list = ensure_iterable(reaction_list)

        # Add reactions using inherited method
        super(EnzymeModule, self).add_reactions(reaction_list)

        # Get the enzyme module reactions by checking if an
        # EnzymeModuleSpecies(s) are involved, and check whether reaction
        # exists, ignoring those that do.
        enzyme_module_reactions = [
            r for r in reaction_list if isinstance(r, EnzymeModuleReaction)
            and r in self.reactions and r not in self.enzyme_module_reactions]

        # Add enzyme module reactions to the enzyme reaction DictList
        if enzyme_module_reactions:
            self.enzyme_module_reactions += enzyme_module_reactions
        # Add enzyme module reactions to the categorized attr as "Undefined".
        self.set_enzyme_object_category("Undefined", enzyme_module_reactions)

    def remove_reactions(self, reaction_list, remove_orphans=False):
        r"""Remove reactions from the :class:`EnzymeModule`.

        The change is reverted upon exit when using the :class:`EnzymeModule`
        as a context.

        Notes
        -----
        Extends from :meth:`.MassModel.remove_reactions`.

        Parameters
        ----------
        reaction_list : list
            A list of :class:`~.MassReaction` and
            :class:`~.EnzymeModuleReaction` to remove from the
            :class:`EnzymeModule`.
        remove_orphans : bool
            If ``True``, will also remove orphaned genes,
            :class:`~.MassMetabolite`\ s, and :class:`~.EnzymeModuleSpecies`
            from the :class:`EnzymeModule`.

        """
        # Ensure list is iterable.
        reaction_list = ensure_iterable(reaction_list)
        # Get the enzyme module reactions by checking if EnzymeModuleSpecies(s)
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
        :attr:`enzyme_module_reactions` attribute will be updated to
        ensure it contains all existing reactions involving
        :class:`~.EnzymeModuleSpecies`.

        Notes
        -----
        Extends from :meth:`.MassModel.repair`.

        Parameters
        ----------
        rebuild_index : bool
            If ``True``, then rebuild the indicies kept in the reactions,
            metabolites, and genes.
        rebuild_relationships: bool
            If ``True``, then reset all associations between the reactions,
            metabolites, genes, and the model, and rebuilds them.

        """
        # Repair using inherited method
        super(EnzymeModule, self).repair(rebuild_index, rebuild_relationships)
        # Repair enzyme_module_reactions DictList
        self._get_current_enzyme_module_objs("reactions", update_enzyme=True)
        self._update_object_pointers()
        # Rebuild DictList indices
        if rebuild_index:
            for attr in ["enzyme_module_ligands", "enzyme_module_species",
                         "enzyme_module_reactions"]:
                getattr(self, attr)._generate_index()
                for value in itervalues(getattr(self, attr + "_categorized")):
                    value._generate_index()

        for enzyme_module_specie in self.enzyme_module_species:
            enzyme_module_specie._repair_bound_obj_pointers()

    # Overridden methods
    def copy(self):
        r"""Create a partial "deepcopy" of the EnzymeModule.

        All of the :class:`~.MassMetabolite`\ s, :class:`~.MassReaction`\ s,
        :class:`~cobra.core.gene.Gene`\ s, :class:`~.EnzymeModuleSpecies`,
        :class:`~.EnzymeModuleReaction`\ s, and :class:`~.EnzymeModuleDict`\ s,
        the boundary conditions, custom rates, custom parameters, and the
        stoichiometric matrix are created anew, but in a faster fashion than
        ``deepcopy``.

        Notes
        -----
        Overrides :meth:`.MassModel.copy`.

        """
        # Define a new model
        new_model = self.__class__()
        # Define items that will not be copied by their references
        do_not_copy_by_ref = [
            "metabolites", "reactions", "genes", "enzyme_modules", "_S",
            "enzyme_module_ligands", "enzyme_module_species",
            "enzyme_module_reactions", "_enzyme_module_ligands_categorized",
            "_enzyme_module_species_categorized",
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
        # in the enzyme_module_ligands, enzyme_module_species and
        # enzyme_module_reactions attributes
        for metabolite in new_model.metabolites:
            if isinstance(metabolite, EnzymeModuleSpecies):
                new_model.enzyme_module_species.append(metabolite)
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
        """Merge two models into one model with the objects from both.

        The reactions, metabolites, genes, enzyme modules, boundary conditions,
        custom rate expressions, rate parameters, compartments, units, notes,
        and annotations from the right model are also copied to left model.
        However, note that in cases where identifiers for objects are identical
        or a dict item has an identical key(s), priority will be given to what
        already exists in the left model.

        Notes
        -----
        * When merging an :class:`.~EnzymeModule` into a :class:`.MassModel`,
          the enzyme module is converted to an :class:`.~EnzymeModuleDict` and
          stored in a :class:`~cobra.core.dictlist.DictList` accessible via the
          :attr:`enzyme_modules` attribute. If an :class:`.~EnzymeModuleDict`
          already exists in the model, it will be replaced.
        * If an :class:`EnzymeModule` already exists in the model, it will
          be replaced.
        * When merging an :class:`EnzymeModule` with another
          :class:`EnzymeModule`, a new :class:`EnzymeModule` will be returned,
          where the EnzymeModule is a copy of the 'left' model (``self``)
          with the ``'right'`` model is contained within.
        * Overrides :meth:`.MassModel.merge`.

        Parameters
        ----------
        right : MassModel
            The model to merge into the left model.
        prefix_existing : str
            If provided, the string is used to prefix the reaction identifier
            of a reaction in the second model if that reaction already exists
            within the left model. Will also apply prefix to enzyme identifiers
            of an enzyme in the second model.
        inplace : bool
            Add reactions from ``right`` directly to left model. Otherwise,
            create a new model leaving the left model untouched.
        new_model_id : str
            If provided, the string is used as the identifier for the merged
            model. If ``None`` and ``inplace=True``, the model ID of the
            first model will be used. If ``None`` and ``inplace=False``, a
            new combined ID will be used for the new model.

        Returns
        -------
        MassModel or EnzymeModule
            A new :class:`~.MassModel` or :class:`EnzymeModule`
            representing the merged model.

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

            # Fix enzyme module ligands, species, and reactions
            for attr in ["ligands", "species", "reactions"]:
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
        for attr in ["enzyme_module_ligands", "enzyme_module_species",
                     "enzyme_module_reactions"]:
            model_dictlist = {
                "enzyme_module_ligands": self.metabolites,
                "enzyme_module_species": self.metabolites,
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
        attr: str {'ligands', 'species', 'reactions'}
            A string representing which attribute to update.
        update_enzyme : bool
            If True, update the enzyme_module_reactions attribute of the
            EnzymeModule.

        Warnings
        --------
        This method is intended for internal use only.

        """
        if attr not in {"ligands", "species", "reactions"}:
            raise ValueError("Unrecognized attribute: '{0}'.".format(attr))
        item_list = []
        if attr == "ligands":
            item_list += list(filter(
                lambda x: not isinstance(x, EnzymeModuleSpecies),
                self.metabolites))
        if attr == "species":
            item_list += list(filter(
                lambda x: isinstance(x, EnzymeModuleSpecies),
                self.metabolites))
        if attr == "reactions":
            for enzyme_module_specie in self.enzyme_module_species:
                item_list += [
                    rxn for rxn in list(enzyme_module_specie.reactions)
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
            EnzymeModuleSpecies: self.enzyme_module_species,
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
                EnzymeModuleSpecies: ("species", "initial conditions"),
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
        enzyme_module_species, and enzyme_module_reactions

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
                "_enzyme_module_species_categorized": "enzyme_module_species",
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
        enzyme_dict = EnzymeModuleDict(self)
        # Update EnzymeModuleDict with attributes
        enzyme_dict._update_object_pointers(model)

        if enzyme_dict.id in model.enzyme_modules:
            model.enzyme_modules.remove(enzyme_dict.id)
        model.enzyme_modules.append(enzyme_dict)
        enzyme_dict.model = model

        return model

    def _repr_html_(self):
        """HTML representation of the overview for the EnzymeModule.

        Warnings
        --------
        This method is intended for internal use only.

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
                    <td><strong>Number of EnzymeModuleSpecies</strong></td>
                    <td>{num_enz_species}</td>
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
                   num_enz_species=len(self.enzyme_module_species),
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


__all__ = ("EnzymeModule",)
