# -*- coding: utf-8 -*-
r"""EnzymeModule is a class for handling reconstructions of enzymes.

The :class:`EnzymeModule` is a reconstruction an enzyme's mechanism and
behavior in a context of a larger system. To aid in the reconstruction process,
the :class:`EnzymeModule` contains various methods to build and add associated
:class:`~.EnzymeModuleForm` and :class:`~.EnzymeModuleReaction`\ s
of the enzyme module (see :meth:`~EnzymeModule.make_enzyme_module_forms` and
:meth:`~~EnzymeModule.make_enzyme_module_reaction` methods, respecitvely).

Given the wide variety of enzymes and the various interactions it can have
with ligands (e.g. catalyzation, inhibition, activation, etc.), the enzyme
module has the following "categorized dict" attributes:

    * :attr:`EnzymeModule.enzyme_module_ligands_categorized`
    * :attr:`EnzymeModule.enzyme_module_forms_categorized`
    * :attr:`EnzymeModule.enzyme_module_reactions_categorized`

These "categorized dict" attributes allow the user to define categories and
place various objects into one or more categories in the respective
"categorized dict" attribute (ligands a.k.a.
:class:`~.MassMetabolite`\ s, :class:`~.EnzymeModuleForm`, and
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

Once merged, the :class:`~.EnzymeModuleForm` and
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
from copy import copy, deepcopy
from functools import partial
from warnings import warn

import numpy as np
import sympy as sym
from cobra.core.dictlist import DictList
from cobra.core.group import Group
from cobra.util.context import get_context
from six import integer_types, iteritems, iterkeys, string_types

from mass.core.mass_model import MassModel
from mass.enzyme_modules.enzyme_module_dict import EnzymeModuleDict
from mass.enzyme_modules.enzyme_module_form import EnzymeModuleForm
from mass.enzyme_modules.enzyme_module_reaction import EnzymeModuleReaction
from mass.util.expressions import _mk_met_func, strip_time
from mass.util.matrix import matrix_rank
from mass.util.util import _mk_new_dictlist, ensure_iterable, ensure_non_negative_value


_AUTOMATIC_RE = re.compile("^automatic$")
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
    enzyme_module_ligands : ~cobra.core.dictlist.DictList
        A :class:`~cobra.core.dictlist.DictList` where the keys are the
        metabolite identifiers and the values are the associated
        :class:`~.MassMetabolite`\ s.
    enzyme_module_forms : ~cobra.core.dictlist.DictList
        A :class:`~cobra.core.dictlist.DictList` where the keys are the
        :class:`~.EnzymeModuleForm` identifiers and the values are the
        associated :class:`~.EnzymeModuleForm`.
    enzyme_module_reactions : ~cobra.core.dictlist.DictList
        A :class:`~cobra.core.dictlist.DictList` where keys are the
        reaction identifiers and the values are the associated
        :class:`~.EnzymeModuleReaction`\ s.
    """

    def __init__(
        self,
        id_or_model=None,
        name=None,
        subsystem="",
        array_type="dense",
        dtype=np.float64,
    ):
        """Initialize the EnzymeModule."""
        super(EnzymeModule, self).__init__(
            id_or_model=id_or_model, name=name, array_type=array_type, dtype=dtype
        )

        self.subsystem = subsystem
        # Initialize DictLists for enzyme ligands, forms, and reactions
        self.enzyme_module_ligands = DictList()
        self.enzyme_module_forms = DictList()
        self.enzyme_module_reactions = DictList()

        # Initialize a dict of DictLists for storing categorized objects
        self._enzyme_module_ligands_categorized = DictList()
        self._enzyme_module_forms_categorized = DictList()
        self._enzyme_module_reactions_categorized = DictList()

        # Initialize EnzymeModule attributes
        self._enzyme_concentration_total = None
        self._enzyme_rate = None
        self._enzyme_rate_equation = None

    @property
    def enzyme_total_symbol_str(self):
        """Get the symbol as a string for the total enzyme concentration."""
        if self.id is None:
            return None

        return str(self.id + "_Total")

    @property
    def enzyme_flux_symbol_str(self):
        """Get the symbol as a string for the net flux through the enzyme."""
        if self.id is None:
            return None

        return str("v_" + self.id)

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
    def enzyme_rate(self):
        """Get or set the flux through the enzyme.

        Parameters
        ----------
        value : float
            The value of the net flux through the enzyme.

        """
        return getattr(self, "_enzyme_rate")

    @enzyme_rate.setter
    def enzyme_rate(self, value):
        """Set the expected flux through the enzyme."""
        if not isinstance(value, (integer_types, float)) and value is not None:
            raise TypeError("Must be an int or float")

        setattr(self, "_enzyme_rate", value)

    @property
    def enzyme_concentration_total_equation(self):
        r"""Return the total concentration equation.

        Notes
        -----
        * Will sum the :class:`~.EnzymeModuleForm` that have their
          :attr:`~.EnzymeModuleForm.enzyme_module_id` match the
          :attr:`EnzymeModule.id`.
        * If no :class:`~.EnzymeModuleForm` are found to have an
          :attr:`~.EnzymeModuleForm.enzyme_module_id` that matches the
          :attr:`EnzymeModule.id`, all :class:`~.EnzymeModuleForm`
          in the the model will be used.

        Returns
        -------
        ~sympy.core.basic.Basic
            A :mod:`sympy` expression of the sum of the
            :class:`~.EnzymeModuleForm`.

        """
        if not self.enzyme_module_forms:
            warn("No EnzymeModuleForm found in EnzymeModule.")
            return None

        # First try only using enzyme module forms that reference this module
        enzyme_module_forms = [
            forms
            for forms in self.enzyme_module_forms
            if forms.enzyme_module_id == self.id
        ]
        # If None found, use all enzyme module forms in the EnzymeModule.
        if not enzyme_module_forms:
            enzyme_module_forms = self.enzyme_module_forms

        return self.sum_enzyme_module_form_concentrations(
            enzyme_module_forms, use_values=False
        )

    @property
    def enzyme_rate_equation(self):
        """Get or set the net flux rate equation of the enzyme.

        Parameters
        ----------
        equation : str, ~sympy.core.basic.Basic
            Either a string representing  the equationcthat will be sympified
            via the :func:`~sympy.core.sympify.sympify` function., or a
            :mod:`sympy` expression representing the of the expression.

        Returns
        -------
        ~sympy.core.basic.Basic
            A :mod:`sympy` expression representing the net flux through the
            enzyme.

        """
        return getattr(self, "_enzyme_rate_equation", None)

    @enzyme_rate_equation.setter
    def enzyme_rate_equation(self, equation):
        """Set the net rate equation of the enzyme."""
        # Set the equation
        if not isinstance(equation, (sym.Basic, string_types)):
            raise TypeError("`equation` must be a str or a sympy expression.")
        if isinstance(equation, string_types):
            equation = sym.sympify(equation)
        setattr(self, "_enzyme_rate_equation", equation)

    @property
    def enzyme_module_ligands_categorized(self):
        r"""Get or set categories for ligands.

        Notes
        -----
        * A ligand must already exist in the :class:`EnzymeModule` as a
          :class:`~.MassMetabolite` in order to set its category.
        * If categories already exists, their existing contents are updated.
        * Setting an empty ``list`` for a category in the dict will cause that
          particular category group to be removed completely from the model.
        * Setting an empty ``dict`` will cause ALL category groups to be
          removed completely from the model.

        Parameters
        ----------
        value : ~cobra.core.group.Group or dict
            Either a :class:`cobra.Group <cobra.core.group.Group>` to add
            to the categorized ligands, or a ``dict`` where keys are strings
            representing categories for the ligands, and values are ``lists``
            containing the corresponding :class:`~.MassMetabolite`\ s or
            their identifiers.

            An empty ``list`` will remove the corresponding category from the
            model and attribute.

            An empty ``dict`` will remove all categories from the attribute
            and the model.

        """
        return getattr(self, "_enzyme_module_ligands_categorized")

    @enzyme_module_ligands_categorized.setter
    def enzyme_module_ligands_categorized(self, value):
        """Create category group(s) for ligands."""
        self._set_category_attribute(
            value, attr="enzyme_module_ligands", to_filter="MassMetabilite"
        )

    @property
    def enzyme_module_forms_categorized(self):
        r"""Get or set categories for enzyme module forms.

        Notes
        -----
        * An enzyme module form must already exist in the
          :class:`EnzymeModule` as an :class:`~.EnzymeModuleForm` in order
          to set its category.
        * If categories already exists, their existing contents are updated.
        * Setting an empty ``list`` for a category in the dict will cause that
          particular category group to be removed completely from the model.
        * Setting an empty ``dict`` will cause ALL category groups to be
          removed completely from the model.

        Parameters
        ----------
        value : ~cobra.core.group.Group or dict
            Either a :class:`cobra.Group <cobra.core.group.Group>` to add
            to the categorized enzyme module forms, or a ``dict`` where keys
            are strings representing categories for the enzyme module forms,
            and values are ``lists`` containing the corresponding
            :class:`~.EnzymeModuleForm`\ s or their identifiers.

            An empty ``list`` will remove the corresponding category from the
            model and attribute.

            An empty ``dict`` will remove all categories from the attribute
            and the model.

        """
        return getattr(self, "_enzyme_module_forms_categorized")

    @enzyme_module_forms_categorized.setter
    def enzyme_module_forms_categorized(self, value):
        """Create category group(s) for enzyme module forms."""
        self._set_category_attribute(
            value, attr="enzyme_module_forms", to_filter="EnzymeModuleForm"
        )

    @property
    def enzyme_module_reactions_categorized(self):
        r"""Get or set categories for enzyme module reactions.

        Notes
        -----
        * An enzyme module reaction must already exist in the
          :class:`EnzymeModule` as an :class:`~.EnzymeModuleReaction` in order
          to set its category.
        * If categories already exists, their existing contents are updated.
        * Setting an empty ``list`` for a category in the dict will cause that
          particular category group to be removed completely from the model.
        * Setting an empty ``dict`` will cause ALL category groups to be
          removed completely from the model.

        Parameters
        ----------
        value : ~cobra.core.group.Group or dict
            Either a :class:`cobra.Group <cobra.core.group.Group>` to add
            to the categorized enzyme module reaction, or a ``dict`` where keys
            are strings representing categories for the enzyme module
            reactions, and values are ``lists`` containing the corresponding
            :class:`~.EnzymeModuleReactions`\ s or their identifiers.

            An empty ``list`` will remove the corresponding category from the
            model and attribute.

            An empty ``dict`` will remove all categories from the attribute
            and the model.

        """
        return getattr(self, "_enzyme_module_reactions_categorized")

    @enzyme_module_reactions_categorized.setter
    def enzyme_module_reactions_categorized(self, value):
        """Create category group(s) for enzyme module reactions."""
        self._set_category_attribute(
            value, attr="enzyme_module_reactions", to_filter="EnzymeModuleReactions"
        )

    def make_enzyme_module_form(
        self,
        id=None,
        name="automatic",
        categories=None,
        bound_metabolites=None,
        compartment=None,
    ):
        r"""Create and add an :class:`~.EnzymeModuleForm` to the module.

        Notes
        -----
        * Upon creation, the :class:`~.EnzymeModuleForm` is added to the
          :class:`EnzymeModule`.
        * If :class:`~.MassMetabolite`\ s in the ``bound_metabolites``
          argument do not already exist in the :class:`EnzymeModule`,
          they will also be added.

        Parameters
        ----------
        id : str
            A string identifier to associate with the enzymatic forms.
        name : str
            Either a human readable name for the enzyme module forms, or the
            string ``"automatic"``. If set to ``"automatic"``, a name will be
            generated based on the identifier of the enzyme module forms and
            its bound ligands.
        categories : str or list
            A string representing the category, or a list of strings
            containing several categories for the enzyme module forms.
        bound_metabolites : dict
            A ``dict`` representing the ligands bound to the enzyme,
            with :class:`~.MassMetabolite`\ s or their identifiers as
            keys and the number bound as values.
        compartment : str
            The compartment where the enzyme module forms is located.

        Returns
        -------
        EnzymeModuleForm
            The newly created :class:`~.EnzymeModuleForm`.

        See Also
        --------
        .EnzymeModuleForm.generate_enzyme_module_forms_name
            Automatic generation of the ``name`` for an
            :class:`~.EnzymeModuleForm`.

        """
        # Ensure metabolites for EnzymeModuleForm exist in the EnzymeModule.
        if bound_metabolites is None:
            bound_metabolites = {}
        else:
            try:
                bound_metabolites = {
                    self.metabolites.get_by_id(str(met)): num
                    for met, num in iteritems(bound_metabolites)
                }
            except KeyError:
                # Add the metabolites into the module that don't already exist
                self.add_metabolites(
                    [
                        met
                        for met in iterkeys(bound_metabolites)
                        if met not in self.metabolites
                    ]
                )
                bound_metabolites = {
                    self.metabolites.get_by_id(str(met)): num
                    for met, num in iteritems(bound_metabolites)
                }

        # Make EnzymeModuleForm object
        enzyme_module_forms = EnzymeModuleForm(
            id_or_specie=id,
            name=name,
            enzyme_module_id=self.id,
            bound_metabolites=bound_metabolites,
            compartment=compartment,
        )
        # Generate name for the EnzymeModuleForm if name set to "automatic"
        if _AUTOMATIC_RE.match(name.lower()):
            enzyme_module_forms.generate_enzyme_module_form_name(True)

        # Add the enzyme forms to the module and place in respective
        # categories if desired.
        self.add_metabolites([enzyme_module_forms])
        if categories is not None:
            categories = ensure_iterable(categories)
            for category in categories:
                self.enzyme_module_forms_categorized = {category: enzyme_module_forms}

        return enzyme_module_forms

    def make_enzyme_module_reaction(
        self,
        id=None,
        name="",
        subsystem=None,
        reversible=True,
        categories=None,
        metabolites_to_add=None,
    ):
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
            on the :class:`~.EnzymeModuleForm` and their bound ligands.
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
            :class:`~.EnzymeModuleForm` or their identifiers as keys and
            stoichiometric coefficients as values. If keys are string
            identifiers then the :class:`~.MassMetabolite`\ s and
            :class:`~.EnzymeModuleForm` must already be a part of model.

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
            id_or_reaction=id,
            name=name,
            subsystem=subsystem,
            reversible=reversible,
            enzyme_module_id=self.id,
        )

        # Add metabolites
        if metabolites_to_add:
            try:
                metabolites_to_add = dict(
                    (met, c)
                    if isinstance(met, EnzymeModuleForm.__base__)
                    else (self.metabolites.get_by_id(met), c)
                    for met, c in iteritems(metabolites_to_add)
                )
            except KeyError as e:
                raise KeyError(str(e) + "not found in model metabolites.")
            new_reaction.add_metabolites(metabolites_to_add)

        # Add reaction to EnzymeModule
        self.add_reactions([new_reaction])

        if categories is not None:
            categories = ensure_iterable(categories)
            for category in categories:
                self.enzyme_module_reactions_categorized = {category: new_reaction}

        # Set enzyme name if set to automatic
        if _AUTOMATIC_RE.match(name.lower()):
            name = new_reaction.generate_enzyme_module_reaction_name(True)

        return new_reaction

    def unify_rate_parameters(
        self, reaction_list, new_parameter_id, rate_type=1, enzyme_prefix=False
    ):
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
        rate_type : int
            The type of rate law to utilize in unification. Must be 1, 2, or 3.

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
            custom_rate = str(strip_time(reaction.get_mass_action_rate(rate_type)))
            custom_rate = custom_rate.replace(reaction.id, new_parameter_id)
            self.add_custom_rate(reaction, custom_rate)

    def make_enzyme_rate_equation(
        self, enzyme_module_reactions, use_rates=False, update_enzyme=False
    ):
        r"""Create an equation representing the net flux through the enzyme.

        The left side of the rate equation will always be the flux symbol of
        the enzyme, accessible via :attr:`enzyme_flux_symbol_str`.

        Parameters
        ----------
        enzyme_module_reactions : list
            A list containing the :class:`~.EnzymeModuleReaction`\ s or their
            identifiers to be summed for the equation.
        use_rates : bool
            If ``True``, then the rates of the provided reactions are used in
            creating the expression. Otherwise the arguments in the expression
            are left as the :attr:`.EnzymeModuleReaction.flux_symbol_str`\ s.
        update_enzyme : bool
            If ``True``, update the :attr:`enzyme_rate_equation` attribute
            and, if necessary, the :attr:`enzyme_module_reactions` attribute
            of the module in addition to returning the generated equation.
            Otherwise just return the net flux equation without making any
            updates. Default is ``False``.

        Returns
        -------
        ~sympy.core.basic.Basic
            A :mod:`sympy` expression of the net flux equation.

        """
        if self.enzyme_flux_symbol_str is None:
            warn("No enzyme flux symbol. Define the EnzymeModule ID first.")
            return None

        # Ensure enzyme_module_reactions are iterable and exist in model
        enzyme_rxns = ensure_iterable(enzyme_module_reactions)
        enzyme_rxns = [
            enz_rxn
            for enz_rxn in enzyme_rxns
            if enz_rxn
            in self._get_current_enzyme_module_objs(
                "reactions", update_enzyme=update_enzyme
            )
        ]

        enzyme_rate_equation = self.sum_enzyme_module_reaction_fluxes(enzyme_rxns)
        if use_rates:
            enzyme_rate_equation = sym.simplify(
                enzyme_rate_equation.subs(
                    {enz_rxn.flux_symbol_str: enz_rxn.rate for enz_rxn in enzyme_rxns}
                )
            )

        if update_enzyme:
            self.enzyme_rate_equation = enzyme_rate_equation

        return enzyme_rate_equation

    def sum_enzyme_module_form_concentrations(
        self, enzyme_module_forms, use_values=False
    ):
        """Sum the forms concentrations for a list of enzyme module forms.

        Parameters
        ----------
        enzyme_module_forms : list
            A ``list`` containing the :class:`~.EnzymeModuleForm` or their
            identifiers to be summed. Forms must already exist in the
            :class:`~.EnzymeModule`.
        use_values : bool
            If ``True``, then numerical values are substituted into the
            expression. Otherwise arguments in the expression are left as
            :mod:`sympy` symbols.

        Returns
        -------
        float or ~sympy.core.basic.Basic
            The sum of the concentrations for the given enzyme module forms
            as a ``float`` if ``use_values`` is ``True`` and all values are
            present. Otherwise returns a ``sympy`` expression representing
            the sum of the given enzyme module forms.

        """
        enzyme_module_forms = ensure_iterable(enzyme_module_forms)
        # Get concentration formula
        concentration = self._make_summation_expr(enzyme_module_forms, EnzymeModuleForm)
        # Substitute numerical values into the formula
        if use_values:
            concentration = self._sub_values_into_expr(
                strip_time(concentration), EnzymeModuleForm
            )

        return concentration

    def sum_enzyme_module_reaction_fluxes(
        self, enzyme_module_reactions, use_values=False
    ):
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
        flux = self._make_summation_expr(enzyme_module_reactions, EnzymeModuleReaction)
        # Substitute numerical values into the formula
        if use_values:
            flux = self._sub_values_into_expr(strip_time(flux), EnzymeModuleReaction)

        return flux

    def enzyme_concentration_total_error(self, use_values=False):
        """Return the error for the total enzyme concentrations.

        The error of the total enzyme concentration is defined to be the
        difference between the :attr:`enzyme_concentration_total` value and
        the sum of all :class:`.EnzymeModuleForm` initial conditions
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
            the sum of the :class:`~.EnzymeModuleForm` initial condition
            values in the model as a ``float`` if ``use_values`` is ``True``
            and all values are present. Otherwise returns a :mod:`sympy`
            expression representing the error.


        """
        if self.enzyme_total_symbol_str is None:
            warn("No enzyme total symbol. Define the EnzymeModule ID first.")
            return None

        if self.enzyme_concentration_total_equation is None:
            warn(
                "No enzyme total concentration equation found. Ensure that "
                "the model contains EnzymeModuleForm and an equation for "
                "the enzyme_concentration_total_equation attribute returns an"
                " expression"
            )
            return None

        # Make error expression
        error = (
            sym.Symbol(self.enzyme_total_symbol_str)
            - self.enzyme_concentration_total_equation
        )
        # Try to substitute values into equations
        if use_values:
            values = {self.enzyme_total_symbol_str: self.enzyme_concentration_total}
            error = self._sub_values_into_expr(
                strip_time(error), EnzymeModuleForm, values
            )

        return error

    def enzyme_rate_error(self, use_values=False):
        """Return the error for the net flux through the enzyme.

        The error of the enzyme net flux is defined to be the difference
        between the :attr:`enzyme_rate` value and the
        calculated value for the :attr:`enzyme_rate_equation`.

        Notes
        -----
        Positive values indicate the value in :attr:`enzyme_rate`
        attribute is greater than the value calculated using the expression
        from the :attr:`enzyme_rate_equation` attribute.

        Parameters
        ----------
        use_values : bool
            If ``True``, then numerical values are substituted into the
            expression. Otherwise arguments in the expression are left as
            :mod:`sympy` symbols.

        Returns
        -------
        float or ~sympy.core.basic.Basic
            The error between the set :attr:`enzyme_rate` and the
            calculated value for the :attr:`enzyme_rate_equation`
            attribute as a ``float`` if ``use_values`` is ``True`` and all
            values are present in the model. Otherwise returns
            a :mod:`sympy` expression representing the error.

        """
        if self.enzyme_flux_symbol_str is None:
            warn("No enzyme flux symbol. Define the EnzymeModule ID first.")
            return None

        if self.enzyme_rate_equation is None:
            warn(
                "No net flux equation found. Ensure that an equation for the "
                "enzyme_rate_equation attribute has been set."
            )
            return None

        # Make error expression
        error = sym.Symbol(self.enzyme_flux_symbol_str) - self.enzyme_rate_equation
        # Try to substitute values into equations
        if use_values:
            error = self._sub_values_into_expr(
                strip_time(error),
                EnzymeModuleReaction,
                {self.enzyme_flux_symbol_str: self.enzyme_rate},
            )

        return error

    def make_enzyme_fraction(self, categorized_attr, top, bottom, use_values=False):
        """Make the expression for a ratio of categorized enzyme objects.

        Notes
        -----
        The string ``"Equation"`` can be passed to either ``top`` or
        ``bottom`` to utilize the equation in the corresponding attribute
        (i.e. :attr:`enzyme_concentration_total_equation` for ``'forms'``
        and :attr:`enzyme_rate_equation` for ``'reactions'``).

        Parameters
        ----------
        categorized_attr: str
            Either a string representing which categorized attribute
            to use or the attribute  itself to use in making the
            enzyme ratio expression. Use the string ``'forms'`` for
            :attr:`enzyme_module_forms_categorized`, or ``'reactions'`` for
            :attr:`enzyme_module_reactions_categorized`.
        top : str
            A string representing a category in the categorized attribute.
            The summation expression of
            the objects in the corresponding list is used as the top
            (numerator) of the fraction to be made. Alternatively, the
            string ``"Equation"`` can be provided to utilize an
            equation attribute.
        bottom : str
            A string representing a category in the categorized attribute.
            The summation expression of the objects in the corresponding list
            is used as the bottom (denominator) of the fraction to be made.
            Alternatively, the string ``"Equation"`` can be provided to
            utilize an equation attribute.
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
        if isinstance(categorized_attr, DictList):
            if categorized_attr == self.enzyme_module_forms_categorized:
                categorized_attr = "forms"
            elif categorized_attr == self.enzyme_module_reactions_categorized:
                categorized_attr = "reactions"
            else:
                raise ValueError(
                    "Must be the attribute accessible through '"
                    "EnzymeModule.enzyme_module_forms_categorized' or "
                    "'EnzymeModule.enzyme_module_reactions_categorized'."
                )

        if not categorized_attr.lower() in {"forms", "reactions"}:
            raise ValueError(
                "Must be a string of the following: " "{'forms', 'reactions'}."
            )

        categorized_attr = categorized_attr.lower()
        # Get object_type for summation expression
        object_type = {"forms": EnzymeModuleForm, "reactions": EnzymeModuleReaction}[
            categorized_attr
        ]

        attr_dictlist = self.__class__.__dict__[
            "enzyme_module_" + categorized_attr + "_categorized"
        ].fget(self)

        # Check top & bottom inputs
        expr = sym.S.One
        for category, coeff in zip([top, bottom], [1, -1]):
            if category in attr_dictlist:
                # Get sum of items if category is not "equation"
                summation_expr = self._make_summation_expr(
                    attr_dictlist.get_by_id(category).members, object_type
                )
            elif _EQUATION_RE.match(category):
                # Get equation if category is "Equation"
                summation_expr = {
                    "forms": self.enzyme_concentration_total_equation,
                    "reactions": self.enzyme_rate_equation,
                }[categorized_attr]
                if summation_expr is None:
                    raise ValueError(
                        "No equation found for '{0}' attribute".format(
                            {
                                "forms": "enzyme_concentration_total_equation",
                                "reactions": "enzyme_rate_equation",
                            }.get(categorized_attr)
                        )
                    )
            else:
                raise ValueError(
                    "Unrecognized category: '{0}' provided for "
                    "'{1}' argument.".format(
                        category, {1: "top", -1: "bottom"}.get(coeff)
                    )
                )
            # Get object type and make expression
            expr = sym.Mul(sym.Pow(summation_expr, coeff), expr)

        # Try to substitute values into equations
        if use_values:
            expr = self._sub_values_into_expr(strip_time(expr), object_type)

        return expr

    def add_metabolites(self, metabolite_list):
        r"""Add a list of metabolites and enzyme forms to the module.

        The change is reverted upon exit when using the :class:`EnzymeModule`
        as a context.

        Notes
        -----
        Extends from :meth:`.MassModel.add_metabolites`.

        Parameters
        ----------
        metabolite_list : list
            A list of :class:`~.MassMetabolite`\ s and
            :class:`EnzymeModuleForm` to add to the :class:`EnzymeModule`.
        add_initial_conditons : bool
            If ``True``, the initial conditions associated with the species
            are also added to the model. Otherwise, the species are added
            without their initial conditions.

        """
        # Ensure list is iterable.
        metabolite_list = ensure_iterable(metabolite_list)

        # Add metabolites using inherited method
        super(EnzymeModule, self).add_metabolites(metabolite_list)

        # Get items that are not EnzymeModuleForm objects, and check if the
        # ligands already exists in the EnzymeModule, ignoring those that do.
        ligands = [
            met
            for met in metabolite_list
            if not isinstance(met, EnzymeModuleForm)
            and met in self.metabolites
            and met not in self.enzyme_module_ligands
        ]
        # Add ligands to the ligand DictList
        if ligands:
            self.enzyme_module_ligands += ligands

        # Get items that are EnzymeModuleForm, and check if the enzyme
        # forms already exist in the EnzymeModule, ignoring those that do.
        enzyme_module_forms = [
            enzyme_module_forms
            for enzyme_module_forms in metabolite_list
            if isinstance(enzyme_module_forms, EnzymeModuleForm)
            and enzyme_module_forms not in self.enzyme_module_forms
        ]

        if enzyme_module_forms:
            # Add the enzyme forms to the model
            self.enzyme_module_forms += enzyme_module_forms

        # Context manager
        context = get_context(self)
        if context:
            context(partial(self.enzyme_module_ligands.__isub__, ligands))
            context(partial(self.enzyme_module_forms.__isub__, enzyme_module_forms))

    def remove_metabolites(self, metabolite_list, destructive=False):
        r"""Remove a list of metabolites and enzyme forms from the module.

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
            :class:`~.EnzymeModuleForm` to remove from the
            :class:`EnzymeModule`.
        destructive : bool
            If ``False``, the species are removed from all associated
            :class:`~.EnzymeModuleReaction`\ s . If ``True``, also remove
            associated :class:`~.EnzymeModuleReaction`\ s from the
            :class:`EnzymeModule`.

        """
        # Ensure list is iterable.
        metabolite_list = ensure_iterable(metabolite_list)

        ligands = [
            met
            for met in metabolite_list
            if not isinstance(met, EnzymeModuleForm)
            and met in self.metabolites
            and met in self.enzyme_module_ligands
        ]

        # Remove metabolites from model using inherited method
        super(EnzymeModule, self).remove_metabolites(metabolite_list, destructive)

        # Remove ligands from the enzyme_module_ligands DictList.
        if ligands:
            self.enzyme_module_ligands -= ligands

        # Get items that are EnzymeModuleForm and check if the
        # enzyme forms already exists in the EnzymeModule, ignoring those
        # that do not.
        enzyme_module_forms = [
            enzyme_module_forms
            for enzyme_module_forms in metabolite_list
            if isinstance(enzyme_module_forms, EnzymeModuleForm)
            and enzyme_module_forms in self.enzyme_module_forms
        ]

        # Remove the enzyme forms to the model
        if enzyme_module_forms:
            self.enzyme_module_forms -= enzyme_module_forms

        # Context manager
        context = get_context(self)
        if context:
            context(partial(self.enzyme_module_ligands.__iadd__, ligands))
            context(partial(self.enzyme_module_forms.__iadd__, enzyme_module_forms))

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

        # Get the enzyme module reactions and check whether reaction
        # exists, ignoring those that do.
        enzyme_module_reactions = [
            r
            for r in reaction_list
            if isinstance(r, EnzymeModuleReaction)
            and r in self.reactions
            and r not in self.enzyme_module_reactions
        ]

        # Add enzyme module reactions to the enzyme reaction DictList
        if enzyme_module_reactions:
            self.enzyme_module_reactions += enzyme_module_reactions

        # Context manager
        context = get_context(self)
        if context:
            context(
                partial(self.enzyme_module_reactions.__isub__, enzyme_module_reactions)
            )

    def remove_reactions(self, reactions, remove_orphans=False):
        r"""Remove reactions from the :class:`EnzymeModule`.

        The change is reverted upon exit when using the :class:`EnzymeModule`
        as a context.

        Notes
        -----
        Extends from :meth:`.MassModel.remove_reactions`.

        Parameters
        ----------
        reactions : list
            A list of :class:`~.MassReaction` and
            :class:`~.EnzymeModuleReaction` to remove from the
            :class:`EnzymeModule`.
        remove_orphans : bool
            If ``True``, will also remove orphaned genes,
            :class:`~.MassMetabolite`\ s, and :class:`~.EnzymeModuleForm`
            from the :class:`EnzymeModule`.

        """
        # Ensure list is iterable.
        reactions = ensure_iterable(reactions)
        # Get the enzyme module reactions and then check whether reaction
        # exists, ignoring those that do not.
        enzyme_module_reactions = [
            r
            for r in reactions
            if isinstance(r, EnzymeModuleReaction)
            and r in self.reactions
            and r in self.enzyme_module_reactions
        ]

        # Remove reactions using inherited method
        super(EnzymeModule, self).remove_reactions(reactions)

        # Remove enzyme module reactions from DictList
        if self.enzyme_module_reactions:
            self.enzyme_module_reactions -= enzyme_module_reactions

        # Context manager
        context = get_context(self)
        if context:
            context(
                partial(self.enzyme_module_reactions.__iadd__, enzyme_module_reactions)
            )

    def repair(self, rebuild_index=True, rebuild_relationships=True):
        """Update all indicies and pointers in the model.

        In addition to updating indicies and pointers, the
        :attr:`enzyme_module_reactions` attribute will be updated to
        ensure it contains all existing reactions involving
        :class:`~.EnzymeModuleForm`.

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
        # # Repair enzyme_module_reactions DictList
        self._get_current_enzyme_module_objs("reactions", update_enzyme=True)
        self._update_object_pointers()
        # # Rebuild DictList indices
        if rebuild_index:
            for attr in [
                "enzyme_module_ligands",
                "enzyme_module_forms",
                "enzyme_module_reactions",
            ]:
                getattr(self, attr)._generate_index()
                getattr(self, attr + "_categorized")._generate_index()

        for enzyme_module_form in self.enzyme_module_forms:
            enzyme_module_form._repair_bound_obj_pointers()

    def copy(self):
        r"""Create a partial "deepcopy" of the EnzymeModule.

        All of the :class:`~.MassMetabolite`\ s, :class:`~.MassReaction`\ s,
        :class:`~cobra.core.gene.Gene`\ s, :class:`~.EnzymeModuleForm`,
        :class:`~.EnzymeModuleReaction`\ s, and :class:`~.EnzymeModuleDict`\ s,
        the boundary conditions, custom rates, custom parameters, and the
        stoichiometric matrix are created anew, but in a faster fashion than
        ``deepcopy``.

        Notes
        -----
        * Overrides :meth:`.MassModel.copy` in order to exclude more items
          to not copy by ref.

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
            "enzyme_module_ligands",
            "enzyme_module_forms",
            "enzyme_module_reactions",
            "_enzyme_module_ligands_categorized",
            "_enzyme_module_forms_categorized",
            "_enzyme_module_reactions_categorized",
            "boundary_conditions",
            "custom_rates",
            "custom_parameters",
            "notes",
            "annotation",
        ]
        for attr in self.__dict__:
            if attr not in do_not_copy_by_ref:
                new_model.__dict__[attr] = self.__dict__[attr]
        new_model.notes = deepcopy(self.notes)
        new_model.annotation = deepcopy(self.annotation)

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

        # Add the newly copied objects to their appropriate DictLists
        # in the enzyme_module_ligands, enzyme_module_forms and
        # enzyme_module_reactions attributes
        for attr in ["ligands", "forms", "reactions"]:
            new_model._get_current_enzyme_module_objs(attr, update_enzyme=True)
            # Update categorized dict attributes
            attr = "enzyme_module_" + attr + "_categorized"
            new_model_categorized_attr = getattr(new_model, attr)
            new_model_categorized_attr += new_model.groups.get_by_any(
                [g.id for g in getattr(self, attr)]
            )

        # Create the new stoichiometric matrix for the model.
        new_model._S = self._mk_stoich_matrix(
            array_type=self._array_type, dtype=self._dtype, update_model=True
        )
        try:
            new_model._solver = deepcopy(self.solver)
            # Cplex has an issue with deep copies
        except Exception:
            new_model._solver = copy(self.solver)

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
            The model to merge into the left model. If a :class:`.MassModel`
            then the first model refers to the ``right`` model and the second
            model refers to the ``left`` model. Otherwise the first model
            refers to the ``left`` model and the second model refers to the
            ``right`` model.
        prefix_existing : str
            If provided, the string is used to prefix the reaction identifier
            of a reaction in the second model if that reaction already exists
            within the first model. Will also apply prefix to identifiers
            of enzyme modules in the second model.
        inplace : bool
            If ``True`` then add reactions from second model directly to the
            first model. Otherwise, create a new model leaving the first model
            untouched. When done within the model as context, changes to the
            models are reverted upon exit.
        objective : str
            One of ``"left"``, ``"right"`` or ``"sum"`` for setting the
            objective of the resulting model to that of the corresponding
            model or the sum of both. Default is ``"left"``. Note that when
            merging a :class:`.MassModel` with an :class:`EnzymeModule`,
            ``"left"`` will refer to the :class:`.MassModel`.

        Returns
        -------
        MassModel or EnzymeModule
            A new :class:`~.MassModel` or :class:`EnzymeModule`
            representing the merged model.

        """
        if not isinstance(right, EnzymeModule):
            # Always merge the EnzymeModule into the MassModel
            return right.merge(self, prefix_existing, inplace, objective)

        # Merge the two EnzymeModules together if right is an EnzymeModule
        new_model = MassModel(self).merge(right, prefix_existing, inplace, objective)
        if inplace:
            new_model = self
        else:
            # Reclassify as an EnzymeModule
            new_model = EnzymeModule(new_model)
            # Set EnzymeModule attributes in the new model
            for attr in [
                "subsystem",
                "_enzyme_concentration_total",
                "_enzyme_rate",
                "enzyme_rate_equation",
            ]:
                setattr(new_model, attr, getattr(self, attr))

            # Fix enzyme module ligands, forms, and reactions
            for attr in ["ligands", "forms", "reactions"]:
                # Update DictList attributes
                new_model._get_current_enzyme_module_objs(attr=attr, update_enzyme=True)
                # Update categorized dict attributes
                attr = "_enzyme_module_" + attr + "_categorized"
                new_categorized_attr_ids = [g.id for g in getattr(new_model, attr)]
                for g in getattr(right, attr):
                    if prefix_existing is not None:
                        gid = "{0}{1}".format(prefix_existing, g.id)
                    else:
                        gid = g.id
                    if gid in new_model.groups and gid not in new_categorized_attr_ids:
                        new_categorized_attr_ids += [gid]
                setattr(
                    new_model,
                    attr,
                    _mk_new_dictlist(new_model.groups, getattr(new_model, attr)),
                )

        return new_model

    # Internal
    def _set_category_attribute(self, item, attr, to_filter):
        """Set the categorized attribute after ensuring it is valid.

        Warnings
        --------
        This method is intended for internal use only.

        """
        categorized_attr = getattr(self, "_" + attr + "_categorized")
        # If a group is provided, ensure all members are of the correct type
        # before adding to the attribute and model if not already existing.
        if isinstance(item, Group):
            filter_func = _make_category_filter(to_filter)
            if list(filter(filter_func, item.members)) != len(item.members):
                raise ValueError(
                    "Not all Group members are {0}s. Cannot add the"
                    " Group to {1}_categorized attribute.".format(to_filter, attr)
                )
            # Add to the model and the attribute
            if item.id not in self.groups:
                self.add_groups([item])
            categorized_attr += [item]

        # If a dict is provided, either create new groups or remove categories
        elif isinstance(item, dict):
            # An empty dict means to remove ALL categories from the attribute.
            if not item:
                item = {k.id: [] for k in categorized_attr}
            else:
                item = {k: ensure_iterable(v) for k, v in iteritems(item)}

            for key, value in iteritems(item):
                # Ensure objects exist before adding them.
                for i, v in enumerate(value):
                    invalid = []
                    try:
                        v = getattr(self, attr).get_by_id(getattr(v, "_id", v))
                    except KeyError as e:
                        invalid += [str(e)]
                    else:
                        value[i] = v
                    if invalid:
                        raise ValueError(
                            "Could not find the following in the "
                            "model '{0}': {1}".format(attr, str(invalid))
                        )
                self._set_enzyme_object_category("_" + attr, key, value)

        else:
            raise TypeError(
                "Unrecognized input value. Must be a cobra.Group or a dict."
            )

    def _set_enzyme_object_category(self, attr, category, object_list):
        r"""Add a list of objects to a new or existing category.

        Notes
        -----
        * If a category already exists, the objects will be added to the
          corresponding :class:`cobra.Group <cobra.core.group.Group>`.
        * The objects to be categorized must already exist in the
          :class:`EnzymeModule`.
        * An empty ``object_list`` will cause the group representing the
          category to be removed.

        Parameters
        ----------
        category : str
            A string representing the category for the list of objects to be
            categorized.
        object_list : list
            A ``list`` containing the objects to be categorized.
            The ``list`` must contain ONLY one of the following :mod:`mass`
            object types:

                * :class:`~.MassMetabolite`\ s representing enzyme ligands.
                * :class:`~.EnzymeModuleForm`\ s representing enzymatic
                  forms.
                * :class:`~.EnzymeModuleReaction`\ s representing enzymatic
                  binding reactions ligands.

        """
        if not isinstance(category, string_types):
            raise TypeError("category must be a str")

        # Ensure object list is iterable
        object_list = ensure_iterable(object_list)
        # Seperate objects and ensure all are of the same type
        filter_for = {
            "_enzyme_module_ligands": "MassMetabolite",
            "_enzyme_module_forms": "EnzymeModuleForm",
            "_enzyme_module_reactions": "EnzymeModuleReaction",
        }[attr]
        pruned = list(filter(_make_category_filter(filter_for), object_list))
        # Raise error if more than one tyoe found.
        if len(object_list) != len(pruned):
            raise TypeError(
                "Objects of different types found. Only one object" " type is allowed."
            )

        categorized_attr = getattr(self, attr + "_categorized")
        if category in categorized_attr:
            # Get existing Group and add additional members
            group = categorized_attr.get_by_id(category)
            if object_list:
                group.add_members(object_list)
            else:
                # Remove group if category value is empty
                categorized_attr.remove(group)
                self.remove_groups([group])

        else:
            # Otherwise create a new group with new members.
            if object_list:
                group = Group(id=category, members=object_list)
                self.add_groups([group])
                categorized_attr += [group]
            else:
                warn("No existing group '{0}' to be removed.".format(category))

    def _update_object_pointers(self):
        """Update objects in the attributes to be the objects from the model.

        Warnings
        --------
        This method is intended for internal use only.

        """
        for attr in [
            "enzyme_module_ligands",
            "enzyme_module_forms",
            "enzyme_module_reactions",
        ]:
            model_dictlist = {
                "enzyme_module_ligands": self.metabolites,
                "enzyme_module_forms": self.metabolites,
                "enzyme_module_reactions": self.reactions,
            }.get(attr)
            setattr(self, attr, _mk_new_dictlist(model_dictlist, getattr(self, attr)))
            attr = "_" + attr + "_categorized"
            setattr(self, attr, _mk_new_dictlist(self.groups, getattr(self, attr)))

    def _get_current_enzyme_module_objs(self, attr, update_enzyme=True):
        """Get the enzyme module objects for 'attr' that exist in the model.

        Parameters
        ----------
        attr: str {'ligands', 'forms', 'reactions'}
            A string representing which attribute to update.
        update_enzyme : bool
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
            item_list += list(
                filter(lambda x: not isinstance(x, EnzymeModuleForm), self.metabolites)
            )
        if attr == "forms":
            item_list += list(
                filter(lambda x: isinstance(x, EnzymeModuleForm), self.metabolites)
            )
        if attr == "reactions":
            for enzyme_module_form in self.enzyme_module_forms:
                item_list += [
                    rxn
                    for rxn in list(enzyme_module_form.reactions)
                    if rxn not in item_list
                ]

        item_list = DictList(item_list)
        item_list.sort()
        if (
            set(item_list) ^ set(getattr(self, "enzyme_module_" + attr))
            and update_enzyme
        ):
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
            EnzymeModuleReaction: self.enzyme_module_reactions,
        }.get(object_type)

        # Ensure object_type are appropriate type and exist in model
        items = list(
            filter(lambda x: isinstance(x, object_type) and x in item_list, items)
        )

        return sum(
            [
                sym.Symbol(x.flux_symbol_str)
                if isinstance(x, EnzymeModuleReaction)
                else _mk_met_func(x.id)
                for x in items
            ]
        )

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
        values = {
            str(key): values[str(key)]
            for key in list(expr.atoms(sym.Symbol))
            if str(key) in values
        }

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
                EnzymeModuleReaction: ("reactions", "steady state fluxes"),
            }.get(object_type)
            warn(
                "Not all enzyme {0} have {1} defined in model, will return a "
                "sympy expression".format(*message_strs)
            )

        return expr

    def _add_self_to_model(self, model, prefix_existing, inplace, objective):
        """Add self to the model and return the MassModel object.

        Warnings
        --------
        This method is intended for internal use only.

        """
        # Switch the objective to match the switch in the merge order.
        objective = {"left": "right", "right": "left", "sum": "sum"}[objective]

        # Create a MassModel instance of self to merge normally
        model = model.merge(
            MassModel(self),
            prefix_existing=prefix_existing,
            inplace=inplace,
            objective=objective,
        )

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
            rank = matrix_rank(self.S)
        except (np.linalg.LinAlgError, ValueError, IndexError):
            dim_S = "{0}x{1}".format(len(self.metabolites), len(self.reactions))
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
                    <td><strong>Subsystem</strong></td>
                    <td>{subsystem}</td>
                </tr><tr>
                    <td><strong>Number of ligands</strong></td>
                    <td>{num_enzyme_module_ligands}</td>
                </tr><tr>
                    <td><strong>Number of enzyme module forms</strong></td>
                    <td>{num_enz_forms}</td>
                </tr><tr>
                    <td><strong>Initial conditions defined</strong></td>
                    <td>{num_ic}/{num_metabolites}</td>
                </tr><tr>
                    <td><strong>Number of enzyme module reactions</strong></td>
                    <td>{num_enz_reactions}</td>
                </tr><tr>
                    <td><strong>Total enzyme concentration</strong></td>
                    <td>{enz_conc}</td>
                </tr><tr>
                    <td><strong>Enzyme rate</strong></td>
                    <td>{enz_flux}</td>
                </tr><tr>
                    <td><strong>Number of groups</strong></td>
                    <td>{num_groups}</td>
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
            subsystem=self.subsystem,
            num_enzyme_module_ligands=len(self.enzyme_module_ligands),
            num_enz_forms=len(self.enzyme_module_forms),
            num_ic=len(self.initial_conditions),
            num_metabolites=len(self.metabolites),
            num_enz_reactions=len(self.enzyme_module_reactions),
            enz_conc=self.enzyme_concentration_total,
            enz_flux=self.enzyme_rate,
            num_groups=len(self.groups),
            compartments=", ".join(
                v if v else k for k, v in iteritems(self.compartments)
            ),
        )


def _make_category_filter(filter_for):
    """Make a filter function for one of the categorized attributes.

    Warnings
    --------
    This method is intended for internal use only.

    """
    if filter_for in "MassMetabolite":

        def ligand_filter_function(obj):
            """Make type filter for MassMetabolites."""
            if isinstance(obj, EnzymeModuleForm.__base__) and not isinstance(
                obj, EnzymeModuleForm
            ):
                return True

            return False

        type_filter_function = ligand_filter_function

    if filter_for == "EnzymeModuleForm":

        def enzyme_forms_filter_function(obj):
            """Make type filter for EnzymeModuleForm."""
            if isinstance(obj, EnzymeModuleForm):
                return True

            return False

        type_filter_function = enzyme_forms_filter_function

    if filter_for == "EnzymeModuleReaction":

        def enzyme_reaction_filter_function(obj):
            """Make type filter for EnzymeModuleReaction."""
            if isinstance(obj, EnzymeModuleReaction):
                return True

            return False

        type_filter_function = enzyme_reaction_filter_function

    return type_filter_function


__all__ = ("EnzymeModule",)
