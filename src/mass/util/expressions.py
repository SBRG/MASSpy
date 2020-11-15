# -*- coding: utf-8 -*-
"""Handles generation and manipulation of :mod:`sympy` expressions."""
import re
from warnings import warn

import sympy as sym
from six import iteritems, iterkeys, string_types
from sympy.physics.vector import dynamicsymbols

from mass.core.mass_configuration import MassConfiguration


MASSCONFIGURATION = MassConfiguration()


# Public
def Keq2k(sympy_expr, simplify=False):
    """Replace ``'Keq'`` symbols with ``'kf/kr'`` in :mod:`sympy` expressions.

    Parameters
    ----------
    sympy_expr : ~sympy.core.basic.Basic, dict, or list
        A :mod:`sympy` expression, a list of :mod:`sympy` expressions, or a
        dictionary with :mod:`sympy` expressions as the values.
    simplify : bool
        If ``True`` then try to simplify the expression after making the
        substitution. Otherwise leave the expression as is.

    Returns
    -------
    ~sympy.core.basic.Basic, dict, or list
        The :mod:`sympy` expression(s) with the substitution made, returned as
        the same type as the original input.

    """

    def _replace_Keq(expr, simplify):
        """Replace the Keq symbol with kf/kr."""
        identifiers = [
            str(symbol).split("_", 1)[1]
            for symbol in list(expr.atoms(sym.Symbol))
            if str(symbol).startswith("Keq_")
        ]
        # Return the expression if no Keq found.
        if not identifiers:
            return expr

        substituion_dict = {}
        # Create substitution dict
        for pid in identifiers:
            kf, kr, Keq = (
                sym.Symbol(param_type + "_" + str(pid))
                for param_type in ["kf", "kr", "Keq"]
            )
            substituion_dict[Keq] = kf / kr
        # Substitute Keq for kf/kr
        new_expr = expr.subs(substituion_dict)
        if len(identifiers) == 1:
            new_expr = sym.collect(new_expr, kf)

        # Simplify if desired
        if simplify:
            new_expr = sym.simplify(new_expr)

        return new_expr

    new_expr = _apply_func_to_expressions(sympy_expr, _replace_Keq, [simplify])

    return new_expr


def k2Keq(sympy_expr, simplify=False):
    """Replace ``'kr'`` symbols with ``'kf/Keq'`` in :mod:`sympy` expressions.

    Parameters
    ----------
    sympy_expr : ~sympy.core.basic.Basic, dict, or list
        A :mod:`sympy` expression, a list of :mod:`sympy` expressions, or a
        dictionary with :mod:`sympy` expressions as the values.
    simplify : bool
        If ``True`` then try to simplify the expression after making the
        substitution. Otherwise leave the expression as is.

    Returns
    -------
    ~sympy.core.basic.Basic, dict, or list
        The :mod:`sympy` expression(s) with the substitution made, returned as
        the same type as the original input.

    """

    def _replace_kr(expr, simplify):
        """Replace the Keq symbol with kf/kr."""
        if not isinstance(expr, sym.Basic):
            raise TypeError("{0} is not a sympy expression".format(str(expr)))
        identifiers = [
            str(symbol).split("_", 1)[1]
            for symbol in list(expr.atoms(sym.Symbol))
            if str(symbol).startswith("kr_")
        ]
        # Return the expression if no kr found.
        if not identifiers:
            return expr

        substituion_dict = {}
        # Create substitution dict
        for pid in identifiers:
            kf, kr, Keq = (
                sym.Symbol(param_type + "_" + str(pid))
                for param_type in ["kf", "kr", "Keq"]
            )
            substituion_dict[kr] = kf / Keq
        # Substitute kr for kf/Keq
        new_expr = expr.subs(substituion_dict)

        if len(identifiers) == 1:
            new_expr = sym.collect(new_expr, kf)

        # Simplify if desired
        if simplify:
            new_expr = sym.simplify(new_expr)

        return new_expr

    new_expr = _apply_func_to_expressions(sympy_expr, _replace_kr, [simplify])

    return new_expr


def strip_time(sympy_expr):
    """Strip the time dependency in :mod:`sympy` expressions.

    Parameters
    ----------
    sympy_expr : ~sympy.core.basic.Basic, dict, or list
        A :mod:`sympy` expression, a list of :mod:`sympy` expressions, or a
        dictionary with :mod:`sympy` expressions as the values.

    Returns
    -------
    ~sympy.core.basic.Basic, dict, or list
        The :mod:`sympy` expression(s) with the time dependency removed,
        returned as the same type as the original input.

    """
    # Helper function to strip a single expression
    def _strip_single_expr(expr):
        if not isinstance(expr, sym.Basic):
            raise TypeError("{0} is not a sympy expression".format(str(expr)))
        # Get the functions of only time.
        subs_dict = {}
        funcs = list(expr.atoms(sym.Function))
        for func in funcs:
            if len(func.atoms(sym.Function)) == 1 and func.atoms(
                sym.Symbol
            ).pop() == sym.Symbol("t"):
                # Make symbol to replace function
                subs_dict[func] = sym.Symbol(str(func)[:-3])
        # Substitute functions for symbols
        new_expr = expr.subs(subs_dict)

        return new_expr

    new_expr = _apply_func_to_expressions(sympy_expr, _strip_single_expr)

    return new_expr


def generate_mass_action_rate_expression(reaction, rate_type=1):
    """Generate the mass action rate law for the reaction.

    Parameters
    ----------
    reaction : MassReaction
        The reaction to generate the rate expression for.
    rate_type : int
        The type of rate law to return. Must be 1, 2, or 3.

            * Type 1 will utilize the
              :attr:`~.MassReaction.forward_rate_constant` and the
              :attr:`~.MassReaction.equilibrium_constant`.
            * Type 2 will utilize the
              :attr:`~.MassReaction.forward_rate_constant` and the
              :attr:`~.MassReaction.reverse_rate_constant`.
            * Type 3 will utilize the
              :attr:`~.MassReaction.equilibrium_constant` and the
              :attr:`~.MassReaction.reverse_rate_constant`.

        Default is ``1``.

    Returns
    -------
    ~sympy.core.basic.Basic or None
        The rate law as a :mod:`sympy` expression. If the reaction has no
        metabolites associated, ``None`` will be returned.

    """
    if not reaction.metabolites:
        warn("No metabolites exist in reaction '{0}'.".format(reaction.id))
        return None

    # Generate forward rate expression
    fwd_rate = generate_forward_mass_action_rate_expression(reaction, rate_type)

    # Ignore reverse rate if it is mathematically equal to 0, or if
    # the equilibrium and rate constants are None and reaction is irreversible
    if not reaction.reversible:
        rate_expression = fwd_rate
    else:
        # Generate reverse rate expression
        rev_rate = generate_reverse_mass_action_rate_expression(reaction, rate_type)
        rate_expression = sym.Add(fwd_rate, sym.Mul(-sym.S.One, rev_rate))

    # Try to group the forward rate constants
    if rate_type == 1:
        rate_expression = sym.collect(rate_expression, reaction.kf_str)

    # Try to group compartments in the rate
    if not MASSCONFIGURATION.exclude_compartment_volumes_in_rates:
        for c in list(reaction.compartments):
            rate_expression = sym.collect(rate_expression, "volume_" + c)

    return rate_expression


def generate_forward_mass_action_rate_expression(reaction, rate_type=1):
    """Generate the forward mass action rate expression for the reaction.

    Parameters
    ----------
    reaction : MassReaction
        The reaction to generate the rate expression for.
    rate_type : int
        The type of rate law to return. Must be 1, 2, or 3.

            * Type 1 and 2 will utilize the
              :attr:`~.MassReaction.forward_rate_constant`.
            * Type 3 will utilize the
              :attr:`~.MassReaction.equilibrium_constant` and the
              :attr:`~.MassReaction.reverse_rate_constant`.

        Default is `1`.

    Returns
    -------
    ~sympy.core.basic.Basic or None
        The forward rate as a :mod:`sympy` expression. If the reaction
        has no metabolites associated, ``None`` will be returned.

    """
    if not reaction.metabolites:
        warn("No metabolites exist in reaction '{0}'.".format(reaction.id))
        return None

    if MASSCONFIGURATION.exclude_metabolites_from_rates and not reaction.boundary:
        reaction = _remove_metabolites_from_rate(reaction)

    fwd_rate = _format_metabs_sym(sym.S.One, reaction, reaction.reactants)
    if rate_type == 3:
        fwd_rate = sym.Mul(
            sym.Mul(sym.var(reaction.kr_str), sym.var(reaction.Keq_str)), fwd_rate
        )
    else:
        fwd_rate = sym.Mul(sym.var(reaction.kf_str), fwd_rate)

    # Remove time dependency from fixed metabolites
    fwd_rate = _set_fixed_metabolites_in_rate(reaction, fwd_rate)

    # Add compartments
    if not MASSCONFIGURATION.exclude_compartment_volumes_in_rates:
        compartments = set(
            met.compartment for met in reaction.reactants if met is not None
        )
        for c in list(compartments):
            fwd_rate = sym.Mul(fwd_rate, sym.Symbol("volume_" + c))

    return fwd_rate


def generate_reverse_mass_action_rate_expression(reaction, rate_type=1):
    """Generate the reverse mass action rate expression for the reaction.

    Parameters
    ----------
    reaction : MassReaction
        The reaction to generate the rate expression for.
    rate_type : int
        The type of rate law to return. Must be 1, 2, or 3.

            * Type 1 will utilize the
              :attr:`~.MassReaction.forward_rate_constant` and the
              :attr:`~.MassReaction.equilibrium_constant`.
            * Type 2 and 3 will utilize the
              :attr:`~.MassReaction.reverse_rate_constant`.

        Default is `1`.

    Returns
    -------
    ~sympy.core.basic.Basic or None
        The reverse rate as a :mod:`sympy` expression. If the reaction
        has no metabolites associated, ``None`` will be returned.

    """
    if not reaction.metabolites:
        warn("No metabolites exist in reaction '{0}'.".format(reaction.id))
        return None

    if MASSCONFIGURATION.exclude_metabolites_from_rates and not reaction.boundary:
        reaction = _remove_metabolites_from_rate(reaction)

    rev_rate = _format_metabs_sym(sym.S.One, reaction, reaction.products)
    if rate_type == 1:
        rev_rate = sym.Mul(
            sym.Mul(sym.var(reaction.kf_str), sym.Pow(sym.var(reaction.Keq_str), -1)),
            rev_rate,
        )
    else:
        rev_rate = sym.Mul(sym.var(reaction.kr_str), rev_rate)

    # Remove time dependency from fixed metabolites
    rev_rate = _set_fixed_metabolites_in_rate(reaction, rev_rate)

    # Add compartments
    if not MASSCONFIGURATION.exclude_compartment_volumes_in_rates:
        compartments = set(
            met.compartment for met in reaction.products if met is not None
        )
        for c in list(compartments):
            rev_rate = sym.Mul(rev_rate, sym.Symbol("volume_" + c))

    return rev_rate


def generate_mass_action_ratio(reaction):
    """Generate the mass action ratio for a given reaction.

    Parameters
    ----------
    reaction : MassReaction
        The reaction to generate the mass action ratio for.

    Returns
    -------
    ~sympy.core.basic.Basic
        The mass action ratio as a :mod:`sympy` expression.

    """
    if MASSCONFIGURATION.exclude_metabolites_from_rates and not reaction.boundary:
        reaction = _remove_metabolites_from_rate(reaction)

    # Handle reactants
    r_bits = _format_metabs_sym(sym.S.One, reaction, reaction.reactants)
    # Handle products
    p_bits = _format_metabs_sym(sym.S.One, reaction, reaction.products)
    # Combine to make the mass action ratio
    ma_ratio = sym.Mul(p_bits, sym.Pow(r_bits, -1))

    return ma_ratio


def generate_disequilibrium_ratio(reaction):
    """Generate the disequilibrium ratio for a given reaction.

    Parameters
    ----------
    reaction: MassReaction
        The reaction to generate the disequilibrium ratio for.

    Returns
    -------
    ~sympy.core.basic.Basic
        The disequilibrium ratio as a :mod:`sympy` expression.

    """
    diseq_ratio = sym.Mul(
        generate_mass_action_ratio(reaction), sym.Pow(sym.var(reaction.Keq_str), -1)
    )

    return diseq_ratio


def create_custom_rate(reaction, custom_rate, custom_parameters=None):
    """Create a :mod:`sympy` expression for a given custom rate law.

    Notes
    -----
    * Metabolites must already exist in the :class:`~.MassModel` or
      :class:`~.MassReaction`.
    * Default parameters of a :class:`~.MassReaction` are automatically
      taken into account and do not need to be defined as additional
      custom parameters.

    Parameters
    ----------
    reaction : MassReaction
        The reaction associated with the custom rate.
    custom_rate : str
        The custom rate law as a str. The string representation of the
        custom rate law will be used to create the expression through the
        :func:`~sympy.core.sympify.sympify` function.
    custom_parameters : list of str
        The custom parameter(s) of the custom rate law as a list of strings.
        The string representation of the custom parameters will be used for
        creation and recognition of the custom parameter symbols in the
        :mod:`sympy` expression. If ``None`` then parameters are assumed to be
        one or more of the reaction rate or equilibrium constants.

    Returns
    -------
    ~sympy.core.basic.Basic or None
        A :mod:`sympy` expression of the custom rate. If no metabolites are
        assoicated with the reaction, ``None`` will be returned.

    See Also
    --------
    :attr:`.MassReaction.all_parameter_ids`
        List of default reaction parameters automatically accounted for.

    """
    # Check inputs
    if not reaction._metabolites:
        warn("No metabolites exist in reaction '{0}'.".format(reaction.id))

    model = reaction.model

    if not isinstance(custom_rate, string_types):
        raise TypeError("custom_rate must be a string")

    if custom_parameters:
        if not hasattr(custom_parameters, "__iter__"):
            custom_parameters = [custom_parameters]
        for custom_param in custom_parameters:
            if not isinstance(custom_param, string_types):
                raise TypeError(
                    "custom_parameters must be a string or " "a list of strings"
                )
    else:
        custom_parameters = []

    custom_rate_expr = custom_rate.replace("(t)", "")

    # Get metabolites as symbols if they are in the custom rate law
    obj_iter = iterkeys(reaction.metabolites) if model is None else model.metabolites
    met_syms = {
        str(met): _mk_met_func(met)
        for met in obj_iter
        if re.search(str(met), custom_rate_expr)
    }

    # Get fixed concentrations as symbols if they are in the custom rate law
    fix_syms = {}
    if reaction._model is not None:
        for attr in ["fixed", "boundary_metabolites"]:
            fix_syms = {
                str(met): sym.Symbol(str(met))
                for met in getattr(reaction._model, attr)
                if re.search("[" + str(met) + "]", custom_rate_expr)
            }

    # Get rate parameters as symbols if they are in the custom rate law
    rate_syms = {
        getattr(reaction, p): sym.Symbol(getattr(reaction, p))
        for p in ["kf_str", "Keq_str", "kr_str"]
        if re.search(str(getattr(reaction, p)), custom_rate_expr)
    }
    # Get custom parameters as symbols
    custom_syms = {custom: sym.Symbol(custom) for custom in custom_parameters}

    # Create custom rate expression
    symbol_dict = {}
    for dictionary in [met_syms, fix_syms, rate_syms, custom_syms]:
        symbol_dict.update(dictionary)
    custom_rate_expr = sym.sympify(custom_rate_expr, locals=symbol_dict)
    custom_rate_expr = _set_fixed_metabolites_in_rate(reaction, custom_rate_expr)
    return custom_rate_expr


def generate_ode(metabolite):
    """Generate the ODE for a given metabolite as a :mod:`sympy` expression.

    Parameters
    ----------
    metabolite : MassMetabolite
        The metabolite to generate the ODE for.

    Returns
    -------
    ode : ~sympy.core.basic.Basic or None
        A :mod:`sympy` expression of the metabolite ODE. If the metabolite
        is not associated with any reactions, then ``None`` will be returned.

    """
    if metabolite._reaction:
        ode = sym.S.Zero
        if not metabolite.fixed:
            for rxn in metabolite._reaction:
                ode = sym.Add(
                    ode, sym.Mul(rxn.get_coefficient(metabolite.id), rxn.rate)
                )
    else:
        ode = None

    return ode


def _remove_metabolites_from_rate(reaction):
    """Remove metabolites from a copy of the reaction before creating the rate.

    Warnings
    --------
    This method is intended for internal use only.

    """
    rxn = reaction.copy()
    # Get exclusion criteria and reaction metabolites
    exclusion_criteria_dict = MASSCONFIGURATION.exclude_metabolites_from_rates
    metabolites_to_exclude = []
    # Iterate through attributes and exclusion values
    for attr, exclusion_values in iteritems(exclusion_criteria_dict):
        exclusion_values = [
            getattr(value, attr) if hasattr(value, attr) else value
            for value in exclusion_values
        ]
        # Iterate through reaction metabolites
        for met in list(rxn.metabolites):
            met_value = getattr(met, attr)
            # Add metabolite to be excluded if it matches the criteria
            if met_value in exclusion_values:
                metabolites_to_exclude += [str(met)]

    # Remove metabolites from reaction copy
    rxn.subtract_metabolites(
        {
            met: coeff
            for met, coeff in iteritems(rxn.metabolites)
            if str(met) in metabolites_to_exclude
        }
    )

    # If all metabolites were removed, leave the reaction as is
    if not rxn.metabolites:
        rxn = reaction.copy()

    return rxn


# Internal
def _mk_met_func(met):
    """Make an undefined sympy.Function of time.

    Warnings
    --------
    This method is intended for internal use only.

    """
    return dynamicsymbols(str(met))


def _set_fixed_metabolites_in_rate(reaction, rate):
    """Strip time dependency of fixed metabolites in the rate expression.

    Warnings
    --------
    This method is intended for internal use only.

    """
    to_strip = [
        str(metabolite) for metabolite in list(reaction.metabolites) if metabolite.fixed
    ]

    if reaction.model is not None and reaction.model.boundary_conditions:
        to_strip += [
            met
            for met, value in iteritems(reaction.model.boundary_conditions)
            if not isinstance(value, sym.Basic)
        ]
    if to_strip:
        to_sub = {
            _mk_met_func(met): sym.Symbol(met)
            for met in to_strip
            if _mk_met_func(met) in list(rate.atoms(sym.Function))
        }
        rate = rate.subs(to_sub)

    return rate


def _format_metabs_sym(expr, rxn, mets):
    """Format the metabolites for a rate law or ratio sympy expression."""
    # For boundary reactions, generate an "boundary" metabolite for boundary
    if rxn.boundary and not mets:
        expr = sym.Mul(expr, sym.Symbol(rxn.boundary_metabolite))
    # For all other reactions
    else:
        for met in mets:
            met_ode = _mk_met_func(met)
            coeff = abs(rxn.get_coefficient(met.id))
            if coeff == 1:
                expr = sym.Mul(expr, met_ode)
            else:
                expr = sym.Mul(expr, sym.Pow(met_ode, coeff))
    return expr


def _apply_func_to_expressions(sympy_expr, function, args=None):
    """Apply the given function to alter each sympy expression provided.

    Warnings
    --------
    This method is intended for internal use only.

    """
    if args is None:

        def func(expr):
            return function(expr)

    else:

        def func(expr):
            return function(expr, *args)

    if isinstance(sympy_expr, dict):
        new_expr = dict((k, func(expr)) for k, expr in iteritems(sympy_expr))
    elif hasattr(sympy_expr, "__iter__"):
        new_expr = list(func(expr) for expr in sympy_expr)
    else:
        new_expr = func(sympy_expr)

    return new_expr


__all__ = (
    "Keq2k",
    "k2Keq",
    "strip_time",
    "generate_mass_action_rate_expression",
    "generate_forward_mass_action_rate_expression",
    "generate_reverse_mass_action_rate_expression",
    "generate_mass_action_ratio",
    "generate_disequilibrium_ratio",
    "create_custom_rate",
    "generate_ode",
)
