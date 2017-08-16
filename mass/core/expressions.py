# -*- coding: utf-8 -*-

# Compatibility with Python 2.7
from __future__ import absolute_import

# Import necesary packages
from six import iterkeys, integer_types
from sympy import Function, Symbol, S, Add, Mul, Pow, var

# from mass
from mass.core import massmetabolite
from mass.core import massreaction
from mass.core import massmodel

# Class begins
## global symbol for time
t = Symbol('t')

# Public
def generate_rate_law(reaction, rate_type=1, sympy_expr=False,
					update_reaction=False):
	"""Generates the rate law for the reaction as a human readable string.
	or as a sympy expression for simulation.

	The type determines which rate law format to return.
	For example: A <=> B

	type=1: kf*(A - B/Keq)
	type=2: kf*A - kr*B
	type=3: kr*(Keq*A - B)

	Parameters
	----------
	rate_type : int {1, 2, 3}
		The type of rate law to display. Must be 1, 2, of 3.
		type 1 will utilize kf and Keq,
		type 2 will utilize kf and kr,
		type 3 will utilize kr and Keq.
	sympy_expr : bool
		If True, will output a sympy expression, otherwise
		will output a human readable string.
	update_reaction : bool
		If True, will update the MassReaction in addition to returning the
		rate law. Otherwise just return the rate law.

	Returns
	-------
	string representation or sympy expression of the rate law
	"""
	# Check inputs
	if not isinstance(reaction, massreaction.MassReaction):
		raise TypeError("reaction must be a MassReaction")
	elif not isinstance(rate_type, integer_types) and \
		not isinstance(rate_type, float):
		raise TypeError("rate_type must be an int or float")
	elif not isinstance(sympy_expr, bool):
		raise TypeError("sympy_expr must be a bool")
	elif not isinstance(update_reaction, bool):
		raise TypeError("update_reaction must be a bool")
	else:
		rate_type = int(rate_type)

	if len(reaction.metabolites) == 0:
		return None

	rate_constructor = {
		1 : [_generate_rate_type_1, _generate_rate_expr_type_1],
		2 : [_generate_rate_type_2, _generate_rate_expr_type_2],
		3 : [_generate_rate_type_3, _generate_rate_expr_type_3]}

	if rate_type not in iterkeys(rate_constructor):
		raise ValueError("rate_type must be 1, 2, or 3")

	# Construct the rate law
	rate_law = rate_constructor[rate_type][0](reaction)
	rate_law_expr = rate_constructor[rate_type][1](reaction)

	if update_reaction:
		reaction._rate_law_expr = rate_law_expr
		reaction._rate_law = rate_law
		reaction._rtype = rate_type

	if sympy_expr:
		return rate_law_expr
	else:
		return rate_law

def get_mass_action_ratio(reaction, sympy_expr=False):
    """Generate the mass action ratio for the reaction as
    a human readable string or a sympy expression for simulation.

    Parameters
    ----------
    sympy_expr : bool
        If True, will output a sympy expression, otherwise
        will output a human readable string.

    Returns
    -------
    string representation or sympy expression of the mass action ratio
    """
    if not isinstance(reaction, massreaction.MassReaction):
        raise TypeError("reaction must be a MassReaction")
    if not isinstance(sympy_expr, bool):
        raise TypeError("sympy_expr must be a bool")

    mass_action_ratio = _get_mass_action_ratio_expr(reaction)
    if sympy_expr:
        return mass_action_ratio
    else:
        mass_action_ratio = re.sub("[(t)]" ,"", str(mass_action_ratio))
        r = re.search("[/]", mass_action_ratio)
        return "(%s)/(%s)" % (mass_action_ratio[:r.start()],
                        mass_action_ratio[r.end():])


def get_disequilibrium_ratio(reaction, sympy_expr=False):
    """Generate the disequilibrium ratio for the reaction as
    a human readable string or a sympy expression for simulation.

    Parameters
    ----------
    sympy_expr : bool
        If True, will output a sympy expression, otherwise
        will output a human readable string.

    Returns
    -------
    string representation or sympy expression of the disequilibrium ratio
    """
    if not isinstance(reaction, massreaction.MassReaction):
        raise TypeError("reaction must be a MassReaction")
    if not isinstance(sympy_expr, bool):
        raise TypeError("sympy_expr must be a bool")

    disequilibrium_ratio = Mul(_get_mass_action_ratio_expr(reaction),
                            Pow(var(reaction._sym_Keq), -1))
    if sympy_expr:
        return disequilibrium_ratio
    else:
        disequilibrium_ratio = re.sub("[(t)]" ,"",
                                str(disequilibrium_ratio))
        r = re.search("[/]", disequilibrium_ratio)
        return "(%s)/(%s)" % (disequilibrium_ratio[:r.start()],
                        disequilibrium_ratio[r.end():])

def generate_ode(metabolite):
	if not isinstance(metabolite, massmetabolite.MassMetabolite):
		raise TypeError("metabolite must be a MassMetabolite")

	if len(metabolite._reaction) == 0:
		return None

	metabolite._ode = S.Zero
	for rxn in metabolite._reaction:
		if rxn._model is not None and rxn in rxn._model.custom_rates:
			print("FIXME: IMPLEMENT CUSTOM RATES")
			return None
		else:
			rate_law_expr = rxn.rate_law_expr

		if metabolite in rxn.reactants:
			rate_law_expr = Mul(-1, rate_law_expr)
		metabolite._ode = Add(metabolite._ode, rate_law_expr)
	return metabolite._ode

## Internal
def _generate_rate_type_1(reaction):
    """Internal use. Generates the type 1 rate law for the reaction as
    a human readable string.

    To safely generate a rate law, use the generate_rate_law method.
    """
    # Generate forward rate
    rate_law = ""
    # For exchange reactions
    if reaction.exchange and len(reaction.reactants) == 0:
        rate_law += reaction.get_exchange_metabolite
    # For all other reactions
    else:
        for metab in reaction.reactants:
            coeff = reaction.get_coefficient(metab.id)
            if abs(coeff) == 1:
                rate_law += "*%s" % metab.id
            else:
                rate_law += "*%s**%s" % (metab.id, coeff)

    # Return rate if reaction is irreversible
    if not reaction._reversible:
        return reaction._sym_kf + rate_law

    # Generate reverse rate
    rate_law = "%s*(%s - " % (reaction._sym_kf, rate_law.lstrip("*"))
    # For exchange reactions
    if reaction.exchange and len(reaction.products) == 0:
        rate_law += reaction.get_exchange_metabolite
    # For all other reactions
    else:
        for metab in reaction.products:
            coeff = reaction.get_coefficient(metab.id)
            if abs(coeff) == 1:
                rate_law += "%s*" % metab.id
            else:
                rate_law += "%s**%s*" % (metab.id, coeff)

    rate_law = "%s / %s)" % (rate_law.rstrip("*"), reaction._sym_Keq)
    return rate_law

def _generate_rate_type_2(reaction):
    """Internal use. Generates the type 2 rate law for the reaction as
    a human readable string.

    To safely generate a rate law, use the generate_rate_law method.
    """
    # Generate forward rate
    rate_law = ""
    # For exchange reactions
    if reaction.exchange and len(reaction.reactants) == 0:
        # Generate an "external" metabolite for exchanges
        rate_law += "*%s" % reaction.get_exchange_metabolite
    # For all other reactions
    else:
        for metab in reaction.reactants:
            coeff = reaction.get_coefficient(metab.id)
            if abs(coeff) == 1:
                rate_law += "*%s" % metab.id
            else:
                rate_law += "*%s**%s" % (metab.id, coeff)
    # Return rate if reaction is irreversible
    if not reaction._reversible:
        return reaction._sym_kf + rate_law

    # Generate reverse rate
    rate_law = "%s%s - %s" % (reaction._sym_kf, rate_law, reaction._sym_kr)
    if reaction.exchange and len(reaction.products) == 0:
        # Generate an "external" metabolite for exchanges
        rate_law += "*%s" % reaction.get_exchange_metabolite
    # For all other reactions
    else:
        for metab in reaction.products:
            coeff = reaction.get_coefficient(metab.id)
            if abs(coeff) == 1:
                rate_law += "*%s" % metab.id
            else:
                rate_law += "*%s**%s" % (metab.id, coeff)

    return rate_law

def _generate_rate_type_3(reaction):
    """Internal use. Generates the type 3 rate law for the reaction as
    a human readable string.

    To safely generate a rate law, use the generate_rate_law method.
    """
    # Generate forward rate
    rate_law = ""
    # For exchange reactions
    if reaction.exchange and len(reaction.reactants) == 0:
        # Generate an "external" metabolite for exchanges
        rate_law += "*%s" % reaction.get_exchange_metabolite
    # For all other reactions
    else:
        for metab in reaction.reactants:
            coeff = reaction.get_coefficient(metab.id)
            if abs(coeff) == 1:
                rate_law += "*%s" % metab.id
            else:
                rate_law += "*%s**%s" % (metab.id, coeff)
    # Return rate if reaction is irreversible
    if not reaction._reversible:
        return "%s*%s%s" % (reaction._sym_kr, reaction._sym_Keq, rate_law)

    # Generate reverse rate
    rate_law = '%s*(%s%s - ' % (reaction._sym_kr, reaction._sym_Keq, rate_law)
    # For exchange reactions
    if reaction.exchange and len(reaction.products) == 0:
        # Generate an "external" metabolite for exchanges
        rate_law += "%s*" % reaction.get_exchange_metabolite
    # For all other reactions
    else:
        for metab in reaction.products:
            coeff = reaction.get_coefficient(metab.id)
            if abs(coeff) == 1:
                rate_law += "%s*" % metab.id
            else:
                rate_law += "%s**%s*" % (metab.id, coeff)
    rate_law = rate_law.rstrip("*") + ')'
    return rate_law

def _generate_rate_expr_type_1(reaction):
    """Internal use. Generates the type 1 rate law for the reaction as
    a sympy expression for simulation.

    To safely generate a rate law, use the generate_rate_law method.
    """
    # Generate forward rate
    rate_law_f = S.One
    # For exchange reactions
    if reaction.exchange and len(reaction.reactants) == 0:
        # Generate an "external" metabolite for exchanges
        metab_ode = Function(reaction.get_exchange_metabolite)(t)
        rate_law_f = Mul(rate_law_f, metab_ode)
    # For all other reactions
    else:
        for metab in reaction.reactants:
            metab_ode = Function(metab.id)(t)
            coeff = reaction.get_coefficient(metab.id)
            if abs(coeff) == 1:
                rate_law_f = Mul(rate_law_f, metab_ode)
            else:
                rate_law_f = Mul(rate_law_f, Pow(metab_ode, coeff))

    # Return rate if reaction is irreversible
    if not reaction._reversible:
        return Mul(var(reaction._sym_kf), rate_law_f)

    # Generate reverse rate
    rate_law_r = Pow(var(reaction._sym_Keq), -1)
    # For exchange reactions
    if reaction.exchange and len(reaction.products) == 0:
        metab_ode = Function(reaction.get_exchange_metabolite)(t)
        rate_law_r = Mul(rate_law_r, metab_ode)
    # For all other reactions
    else:
        for metab in reaction.products:
            metab_ode = Function(metab.id)(t)
            coeff = reaction.get_coefficient(metab.id)
            if abs(coeff) == 1:
                rate_law_r = Mul(rate_law_r, metab_ode)
            else:
                rate_law_r = Mul(rate_law_r, Pow(metab_ode, coeff))

    # Combine forward and reverse rates, and return rate law
    return Mul(var(reaction._sym_kf), Add(rate_law_f, Mul(-1, rate_law_r)))


def _generate_rate_expr_type_2(reaction):
    """Internal use. Generates the type 2 rate law for the reaction as
    a sympy expression for simulation.

    To safely generate a rate law, use the generate_rate_law method.
    """
    # Generate forward rate
    rate_law_f = var(reaction._sym_kf)
    # For exchange reactions
    if reaction.exchange and len(reaction.reactants) == 0:
        # Generate an "external" metabolite for exchanges
        metab_ode = Function(reaction.get_exchange_metabolite)(t)
        rate_law_f = Mul(rate_law_f, metab_ode)
    # For all other reactions
    else:
        for metab in reaction.reactants:
            metab_ode = Function(metab.id)(t)
            coeff = reaction.get_coefficient(metab.id)
            if abs(coeff) == 1:
                rate_law_f = Mul(rate_law_f, metab_ode)
            else:
                rate_law_f = Mul(rate_law_f, Pow(metab_ode, coeff))

    # Return rate if reaction is irreversible
    if not reaction._reversible:
        return Mul(var(reaction._sym_kf), rate_law_f)

    # Generate reverse rate
    rate_law_r = var(reaction._sym_kr)
    # For exchange reactions
    if reaction.exchange and len(reaction.products) == 0:
        # Generate an "external" metabolite for exchanges
        metab_ode = Function(reaction.get_exchange_metabolite)(t)
        rate_law_r = Mul(rate_law_r, metab_ode)
    # For all other reactions
    else:
        for metab in reaction.products:
            metab_ode = Function(metab.id)(t)
            coeff = reaction.get_coefficient(metab.id)
            if abs(coeff) == 1:
                rate_law_r = Mul(rate_law_r, metab_ode)
            else:
                rate_law_r = Mul(rate_law_r, Pow(metab_ode, coeff))

    # Combine forward and reverse rates, and return rate law
    return Add(rate_law_f, Mul(-1, rate_law_r))

def _generate_rate_expr_type_3(reaction):
    """Internal use. Generates the type 3 rate law for the reaction as
    a sympy expression for simulation.

    To safely generate a rate law, use the generate_rate_law method.
    """
    # Generate forward rate
    rate_law_f = var(reaction._sym_Keq)
    # For exchange reactions
    if reaction.exchange and len(reaction.reactants) == 0:
        # Generate an "external" metabolite for exchanges
        metab_ode = Function(reaction.get_exchange_metabolite)(t)
        rate_law_f = Mul(rate_law_f, metab_ode)
    # For all other reactions
    else:
        for metab in reaction.reactants:
            metab_ode = Function(metab.id)(t)
            coeff = reaction.get_coefficient(metab.id)
            if abs(coeff) == 1:
                rate_law_f = Mul(rate_law_f, metab_ode)
            else:
                rate_law_f = Mul(rate_law_f, Pow(metab_ode, coeff))

    # Return rate if reaction is irreversible
    if not reaction._reversible:
        return Mul(var(reaction._sym_kr), rate_law_f)

    # Generate reverse rate
    rate_law_r = S.One
    # For exchange reactions
    if reaction.exchange and len(reaction.products) == 0:
        # Generate an "external" metabolite for exchanges
        metab_ode = Function(reaction.get_exchange_metabolite)(t)
        rate_law_r = Mul(rate_law_r, metab_ode)
    # For all other reactions
    else:
        for metab in reaction.products:
            metab_ode = Function(metab.id)(t)
            coeff = reaction.get_coefficient(metab.id)
            if abs(coeff) == 1:
                rate_law_r = Mul(rate_law_r, metab_ode)
            else:
                rate_law_r = Mul(rate_law_r, Pow(metab_ode, coeff))

    # Combine forward and reverse rates, and return rate law
    return Mul(var(reaction._sym_kr), Add(rate_law_f, Mul(-1, rate_law_r)))

def _get_mass_action_ratio_expr(reaction):
    """Internal use. Generates the mass action ratio for the reaction as
    a human readable string.

    To safely generate the mass action ratio, use the
    get_mass_action_ratio method.
    """
    # For the reactants
    reactant_bits = S.One
    if reaction.exchange and len(reaction.reactants) == 0:
        # Generate an "external" metabolite for exchanges
        metab_ode = Function(reaction.get_exchange_metabolite)(t)
        reactant_bits = Mul(reactant_bits, metab_ode)
    # For all other reactions
    else:
        for metab in reaction.reactants:
            metab_ode = Function(metab.id)(t)
            coeff = reaction.get_coefficient(metab.id)
            if abs(coeff) == 1:
                reactant_bits = Mul(reactant_bits, metab_ode)
            else:
                reactant_bits = Mul(reactant_bits, Pow(metab_ode, coeff))

    # For the products
    product_bits = S.One
    if reaction.exchange and len(reaction.products) == 0:
        # Generate an "external" metabolite for exchanges
        metab_ode = Function(reaction.get_exchange_metabolite)(t)
        product_bits = Mul(product_bits, metab_ode)
    # For all other reactions
    else:
        for metab in reaction.products:
            metab_ode = Function(metab.id)(t)
            coeff = reaction.get_coefficient(metab.id)
            if abs(coeff) == 1:
                product_bits = Mul(product_bits, metab_ode)
            else:
                product_bits = Mul(product_bits, Pow(metab_ode, coeff))

    # Combine to make the mass action ratio
    return Mul(product_bits, Pow(reactant_bits, -1))
