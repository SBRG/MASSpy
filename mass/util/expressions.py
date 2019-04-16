# -*- coding: utf-8 -*-
"""TODO Module Docstrings."""
from __future__ import absolute_import

import re
from warnings import warn

from six import integer_types, iteritems, iterkeys, string_types

import sympy as sym
from sympy.physics.vector import dynamicsymbols


# Public
def generate_rate_law(reaction, rate_type=1, sympy_expr=True,
                      update_reaction=False):
    """Generate the rate law for the reaction.

    Parameters
    ----------
    reaction: mass.MassReaction
        The MassReaction object to generate the rate expression for.
    rate type: int {1, 2, 3}, optional
        The type of rate law to display. Must be 1, 2, or 3.
            Type 1 will utilize the forward rate and equilibrium constants.
            Type 2 will utilize the forward rate and reverse rate constants.
            Type 3 will utilize the equilibrium and reverse rate constants.
    sympy_expr: bool, optional
        If True, output is a sympy expression. Otherwise output is a string.
    update_reaction: bool, optional
        If True, update the MassReaction in addition to returning the rate law.

    Returns
    -------
    rate_law: str or sympy expression
        The rate law expression. Will return a blank string if no metabolites
        are associated with the reaction.

    """
    # Check inputs
    if not isinstance(rate_type, (integer_types, float)):
        raise TypeError("rate_type must be an int or float")
    elif not isinstance(sympy_expr, bool):
        raise TypeError("sympy_expr must be a bool")
    elif not isinstance(update_reaction, bool):
        raise TypeError("update_reaction must be a bool")
    else:
        rate_type = int(rate_type)

    if not reaction.metabolites:
        return ""

    # Remove H+ and H2O from the reaction if necessary
    rxn = _ignore_h_and_h2o(reaction)

    # Construct the rate law
    rate_constructor = {
        1: [_generate_rate_sym_1, _generate_rate_str_1],
        2: [_generate_rate_sym_2, _generate_rate_str_2],
        3: [_generate_rate_sym_3, _generate_rate_str_3]}

    try:
        rate_expr = rate_constructor[rate_type][0 if sympy_expr else 1](rxn)
    except KeyError:
        raise ValueError("rate_type must be 1, 2, or 3")

    if update_reaction:
        reaction._rtype = rate_type

    if reaction.model is not None:
        rate_expr = _strip_time_for_fixed_mets(reaction, rate_expr)

    return rate_expr


def generate_mass_action_ratio(reaction):
    """Generate the mass action ratio for a given reaction.

    Parameters
    ----------
    reaction: mass.MassReaction
        The MassReaction object to generate the mass action ratio for.

    Returns
    -------
    ma_ratio: sympy.Basic
        The mass action ratio as a sympy expression (sympy.Basic).

    """
    reaction = _ignore_h_and_h2o(reaction)

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
    reaction: mass.MassReaction
        The MassReaction object to generate the disequilibrium ratio for.

    Returns
    -------
    diseq_ratio: sympy.Basic
        The disequilibrium ratio as a sympy expression (sympy.Basic).

    """
    diseq_ratio = sym.Mul(generate_mass_action_ratio(reaction),
                          sym.Pow(sym.var(reaction.Keq_str), -1))

    return diseq_ratio


def generate_ode(metabolite):
    """Generate the ODE for a given metabolite as a sympy expression.

    Parameters
    ----------
    metabolite: mass.MassMetabolite
        The MassMetabolite object to generate the ODE for.
    update_metabolite: bool, optional
        If True, update the MassMetabolite in addition to returning the ODE.

    Returns
    -------
    ode: sympy.Basic
        The metabolite ODE as a sympy expression.

    """
    if metabolite._reaction:
        ode = sym.S.Zero
        for rxn in metabolite._reaction:
            ode = sym.Add(ode, sym.Mul(rxn.get_coefficient(metabolite.id),
                                       rxn.rate))
    else:
        ode = None

    return ode


def create_custom_rate(reaction, custom_rate, custom_parameters=None):
    """Create a sympy expression for a given custom rate law.

    Parameters
    ----------
    reaction: mass.MassReaction
        The MassReaction object to generate the custom rate ratio for.
    custom_rate: str
        The custom rate law as a string. The string representation of the
        custom rate law will be used to create the sympy expression through the
        sympy.sympify method.
    custom_parameters: str, list of str, optional
        The custom parameters of the custom rate law as strings. The string
        representation of the custom parameters will be used for creation and
        recognition of the custom parameter symbols in the sympy expression.
        If None, parameters are assumed to be one or more of the MassReaction
        rate or equilibrium constants.

    Returns
    -------
    custom_rate_expression: sympy.Basic, None
        A sympy expression of the custom rate. If no metabolites are assoicated
        with the reaction, None will be returned.

    Warnings
    --------
    Any metabolite in a custom rate expression must already exist as a
        MassMetabolite in the MassModel or MassReaction object.

    Notes
    -----
    The default MassReaction parameters (kf_RID, Keq_RID, and kr_RID), the
        reaction metabolites, and the fixed concentration parameters from the
        associated MassModel are automatically recognized and should not be
        included the list of custom parameters. RID: MassReaction identifier.

    """
    # Check inputs
    if not reaction._metabolites:
        warn("No metabolites associated with this reaction")
        return None
    model = reaction.model

    if not isinstance(custom_rate, string_types):
        raise TypeError("custom_rate must be a string")

    if custom_parameters:
        if not hasattr(custom_parameters, "__iter__"):
            custom_parameters = [custom_parameters]
        for custom_param in custom_parameters:
            if not isinstance(custom_param, string_types):
                raise TypeError("custom_parameters must be a string or "
                                "a list of strings")
    else:
        custom_parameters = []

    custom_rate_expr = custom_rate.replace("(t)", "")

    # Get metabolites as symbols if they are in the custom rate law
    obj_iter = (iterkeys(reaction.metabolites) 
                if model is None else model.metabolites)
    met_syms = {str(met): _mk_met_func(met)
                for met in obj_iter if re.search(str(met), custom_rate_expr)}

    # Get fixed concentrations as symbols if they are in the custom rate law
    if reaction._model is not None:
        fix_syms = {str(met): sym.Symbol(str(met))
                    for met in reaction._model.fixed_concentrations
                    if re.search("[" + str(met) + "]", custom_rate_expr)}

    else:
        fix_syms = {}

    # Get rate parameters as symbols if they are in the custom rate law
    rate_syms = {getattr(reaction, p): sym.Symbol(getattr(reaction, p))
                 for p in ["kf_str", "Keq_str", "kr_str"]
                 if re.search(str(getattr(reaction, p)), custom_rate_expr)}
    # Get custom parameters as symbols
    custom_syms = {custom: sym.Symbol(custom) for custom in custom_parameters}

    # Create custom rate expression
    symbol_dict = {}
    for dictionary in [met_syms, fix_syms, rate_syms, custom_syms]:
        symbol_dict.update(dictionary)
    custom_rate_expr = sym.sympify(custom_rate_expr, locals=symbol_dict)

    return custom_rate_expr


# Internal
def _ignore_h_and_h2o(reaction):
    """Remove hydrogen and water from reactions.

    Designed for internal use. Remove hydrogen and water from reactions to
    prevent their inclusion in simulation. Does not effect water and hydrogen
    exchange reactions.
    """
    reaction = reaction.copy()
    for met, coeff in iteritems(reaction.metabolites):
        # Must be water or hydrogen.
        if met.elements == {"H": 2, "O": 1} or met.elements == {"H": 1}:
            # Must not be an exchange reaction.
            if not reaction.exchange:
                reaction.subtract_metabolites({met: coeff})

    return reaction


def _format_metabs_str(expr, rxn, mets, left_sign):
    """Format the metabolites for a rate law string."""
    if left_sign:
        l, r = "*", ""
    else:
        l, r = "", "*"
    # For exchange reactions
    if rxn.exchange and not mets:
        expr += l + rxn.external_metabolite + "(t)" + r
    # For all other reactions
    else:
        for met in mets:
            coeff = abs(rxn.get_coefficient(met.id))
            if coeff == 1:
                expr += l + met.id + "(t)" + r
            else:
                expr += l + met.id + "(t)**" + coeff + r
    return expr


def _format_metabs_sym(expr, rxn, mets):
    """Format the metabolites for a rate law or ratio sympy expression."""
    # For exchange reactions, generate an "external" metabolite for exchange
    if rxn.exchange and not mets:
        expr = sym.Mul(expr, sym.Symbol(rxn.external_metabolite))
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


def _generate_rate_str_1(reaction):
    """Generate the type 1 rate law as a human readable string.

    Designed for internal use. Generate the rate law for reaction as a human
    readable string. To safely generate a rate law, use generate_rate_law.
    """
    # Generate forward rate
    rate_law = _format_metabs_str("", reaction, reaction.reactants, True)

    # Return rate if reaction is irreversible
    if not reaction.reversible:
        return reaction.kf_str + rate_law

    # Generate reverse rate
    rate_law = reaction.kf_str + "*(" + rate_law.lstrip("*")
    rate_law = _format_metabs_str(rate_law, reaction, reaction.products, False)
    rate_law = rate_law.rstrip("*") + "/" + reaction.Keq_str + ")"

    return rate_law


def _generate_rate_str_2(reaction):
    """Generate the type 2 rate law as a human readable string.

    Designed for internal use. Generate the rate law for reaction as a human
    readable string. To safely generate a rate law, use generate_rate_law.
    """
    # Generate forward rate
    rate_law = _format_metabs_str("", reaction, reaction.reactants, True)

    # Return rate if reaction is irreversible
    if not reaction.reversible:
        return reaction.kf_str + rate_law

    # Generate reverse rate
    rate_law = reaction.kf_str + rate_law + " - " + reaction.kr_str
    rate_law = _format_metabs_str(rate_law, reaction, reaction.products, True)

    return rate_law


def _generate_rate_str_3(reaction):
    """Generate the type 3 rate law as a human readable string.

    Designed for internal use. Generate the rate law for reaction as a human
    readable string. To safely generate a rate law, use generate_rate_law.
    """
    # Generate forward rate
    rate_law = _format_metabs_str("", reaction, reaction.reactants, True)

    # Return rate if reaction is irreversible
    if not reaction.reversible:
        rate_law = reaction.kr_str + "*" + reaction.Keq_str + rate_law
        return rate_law

    # Generate reverse rate
    rate_law = '{0}*({1}{2} - '.format(reaction.kr_str, reaction.Keq_str,
                                       rate_law)

    rate_law = reaction.kr_str + "*(" + reaction.Keq_str + rate_law + " - "
    rate_law = _format_metabs_str(rate_law, reaction, reaction.products, False)
    rate_law = rate_law.rstrip("*") + ')'

    return rate_law


def _generate_rate_sym_1(reaction):
    """Generate the type 1 rate law as a sympy expression.

    Designed for internal use. Generate the rate law for reaction as a sympy
    expression. To safely generate a rate law, use generate_rate_law.
    """
    # Generate forward rate
    rate_law_f = sym.S.One
    rate_law_f = _format_metabs_sym(rate_law_f, reaction, reaction.reactants)

    # Return rate if reaction is irreversible
    if not reaction.reversible:
        return sym.Mul(sym.var(reaction.kf_str), rate_law_f)

    # Generate reverse rate
    rate_law_r = sym.Pow(sym.var(reaction.Keq_str), -1)
    rate_law_r = _format_metabs_sym(rate_law_r, reaction, reaction.products)

    # Combine forward and reverse rates
    rate_law = sym.Mul(sym.var(reaction.kf_str),
                       sym.Add(rate_law_f, sym.Mul(-1, rate_law_r)))
    return rate_law


def _generate_rate_sym_2(reaction):
    """Generate the type 2 rate law as a sympy expression.

    Designed for internal use. Generate the rate law for reaction as a sympy
    expression. To safely generate a rate law, use generate_rate_law.
    """
    # Generate forward rate
    rate_law_f = sym.var(reaction.kf_str)
    rate_law_f = _format_metabs_sym(rate_law_f, reaction, reaction.reactants)

    # Return rate if reaction is irreversible
    if not reaction.reversible:
        return rate_law_f

    # Generate reverse rate
    rate_law_r = sym.var(reaction.kr_str)
    rate_law_r = _format_metabs_sym(rate_law_r, reaction, reaction.products)

    # Combine forward and reverse rates
    rate_law = sym.Add(rate_law_f, sym.Mul(-1, rate_law_r))

    return rate_law


def _generate_rate_sym_3(reaction):
    """Generate the type 3 rate law as a sympy expression.

    Designed for internal use. Generate the rate law for reaction as a sympy
    expression. To safely generate a rate law, use generate_rate_law.
    """
    # Generate forward rate
    rate_law_f = sym.var(reaction.Keq_str)
    rate_law_f = _format_metabs_sym(rate_law_f, reaction, reaction.reactants)

    # Return rate if reaction is irreversible
    if not reaction.reversible:
        return sym.Mul(sym.var(reaction.kr_str), rate_law_f)

    # Generate reverse rate
    rate_law_r = sym.S.One
    rate_law_r = _format_metabs_sym(rate_law_r, reaction, reaction.products)

    # Combine forward and reverse rates
    rate_law = sym.Mul(sym.var(reaction.kr_str),
                       sym.Add(rate_law_f, sym.Mul(-1, rate_law_r)))

    return rate_law


def _determine_reaction_rtype(reaction):
    """Determine rate to use based on a reaction's available parameters.

    Will return the rate law type based on
    numerically defined parameters, with priority given to rate type 1.

    Warnings
    --------
    This method is intended for internal use only.
    """
    # Get available parameters for a reaction.
    p_types = [param[:3] for param in iterkeys(reaction.parameters)]
    # Ensure reversible reactions have at least two parameters, and define
    # rate type accordingly. Defaults to type 1.
    if reaction.reversible:
        if "kf_" in p_types and "Keq" in p_types:
            rt = 1
        elif "kf_" in p_types and "kr_" in p_types:
            rt = 2
        elif "Keq" in p_types and "kr_" in p_types:
            rt = 3
        else:
            rt = 1
    # Irreversible reactions only require a kf and default to type 1
    else:
        rt = 1

    return rt


def _strip_time_for_fixed_mets(reaction, rate_expr):
    """Strip the time dependency of a sympy Function for fixed metabolites.

    Warnings
    --------
    This method is intended for internal use only.
    """
    mets_to_strip = [met for met in reaction._metabolites
                     if met in reaction.model.fixed_concentrations]
    if mets_to_strip:
        sub_mets_stripped = {_mk_met_func(m): sym.Symbol(str(m))
                             for m in mets_to_strip}
        rate_expr = rate_expr.subs(sub_mets_stripped)
    return rate_expr


def _mk_met_func(met):
    """Make an undefined sympy.Function of time.

    Warnings
    --------
    This method is intended for internal use only.
    """
    return dynamicsymbols(str(met))
