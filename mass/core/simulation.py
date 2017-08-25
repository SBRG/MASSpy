# -*- coding: utf-8 -*-

# Compatibility with Python 2.7
from __future__ import absolute_import

import re
import numpy as np
import sympy as sp
from six import iteritems, iterkeys, itervalues
from scipy.integrate import ode
from mass.core import expressions
from mass.util import qcqa

# Class begins
## Global symbol for time
t = sp.Symbol("t")
## precompiled re for 'external' metabolites
ext_metab_re = re.compile("\_Xt")
# Possible Perturbation Types
kf_re = re.compile("kf|forward_rate_constant")
Keq_re = re.compile("Keq|equilibrium_constant")
kr_re = re.compile("kr|reverse_rate_constant")
ic_re = re.compile("ic|initial_condition")
fixed_re = re.compile("fix|fixed|fixed_concentration")

def simulate(model, time_vector, perturbations=None, solver="vode"):
	# Collect sympy symbols and make dictionariess for odes and rates
	odes, rates, symbols = expressions._sort_symbols(model)
	# Perturb the system if perturbations exist
	if perturbations is not None:
		odes, rates, symbols, perturbations = _perturb(model, odes, rates,
													symbols, perturbations)
	# Get values to substitute into ODEs and the metabolite initial conditions
	values, ics = _get_values(model, perturbations, symbols)

	metab_syms = symbols[0]
	fixed_syms = symbols[2]

	# Make lambda functions of the odes and rates
	[lam_odes, lam_jacb] = _make_lambda_odes(model,metab_syms , odes, values)
	lam_rates = _make_lambda_rates(model, metab_syms, rates, values)

	# Integrate the odes to obtain the concentration solutions
	c = _integrate_odes(time_vector, lam_odes, lam_jacb, ics, solver)

	# Map metbaolite ids to their concentration solutions
	c_profile = dict()
	for i, sym in enumerate(metab_syms):
		c_profile.update({model.metabolites.get_by_id(str(sym)[:-3]) : c[:,i]})
	for i, sym in enumerate(fixed_syms):
		if not ext_metab_re.search(str(sym)):
			metab = model.metabolites.get_by_id(str(sym))
			c_profile.update({metab: values[sym]})

	# Get use the concentrations to get the flux the flux solutions
	# Map reactiom ids to their flux solutions
	f_profile = dict()
	for rxn, lambda_func_and_args in iteritems(lam_rates):
		f = np.zeros(time_vector.shape)
		lambda_func = lambda_func_and_args[1]
		concs = np.array([c_profile[model.metabolites.get_by_id(str(arg))]
							for arg in lambda_func_and_args[0]]).T
		for i in range(0, len(f)):
			f[i] = lambda_func(*concs[i,:])
		f_profile[rxn] = f

	return [c_profile, f_profile]


def _perturb(model, ode_dict, rate_dict, symbol_list, perturbations):
	metabolites = symbol_list[0]
	rate_params = symbol_list[1]
	fixed_concs = symbol_list[2]
	custom_params = symbol_list[3]

	rate_perturbs = dict()
	conc_perturbs = dict()
	# Substitute perturbation values system if there are perturbations
	for pert, value in iteritems(perturbations):
		# Handle Custom Rate parameter perturbations
		if pert in iterkeys(model.custom_parameters):
			print("FIXME: IMPLEMENT CUSTOM PARAMETER PERTURBATIONS")
			continue
		[item_id, to_perturb] = re.split("\.", pert)
		 # Handle rate and equilibrium constant perturbations
		if kf_re.match(to_perturb):
			rate_perturbs.update({"kf_%s" % item_id : value})
			continue
		if Keq_re.match(to_perturb):
			rate_perturbs.update({"Keq_%s" % item_id : value})
			continue
		if kr_re.match(to_perturb):
			rate_perturbs.update({"kr_%s" % item_id : value})
			continue
		# Handle fixed concentration perturbations
		if fixed_re.match(to_perturb):
			# 'External' metabolites
			if ext_metab_re.search(item_id):
				conc_perturbs.update({item_id: value})
				continue
			else:
				# Other metabolites
				metab = model.metabolites.get_by_id(item_id)
				# Use a copy to enable removal of metabolite function
				# from the metabolite set in the iteration
				for metab_func in metabolites.copy():
					if(re.match(metab.id, str(metab_func)[:-3])):
						metab_sym = sp.Symbol(metab.id, nonnegative=True)
						# Remove metabolite from functions and
						# add to fixed_concs set.
						metabolites.remove(metab_func)
						fixed_concs.add(metab_sym)
						# Remove from ODE dictionary
						del ode_dict[metab]
						# Sub fixed concentration into all other ODEs
						for m_key, expression in iteritems(ode_dict):
							ode_dict[m_key] = expression.subs(
												{metab_func: metab_sym})
							for rxn, rate in iteritems(rate_dict):
								rate_dict[rxn] = rate.subs({metab_func:
															metab_sym})
				conc_perturbs.update({metab: value})
				continue
		# Handle initial condition perturbations
		if ic_re.match(to_perturb):
			conc_perturbs.update(
							{model.metabolites.get_by_id(item_id): value})
			continue
	symbol_list = [metabolites, rate_params, fixed_concs, custom_params]
	perturb_list = [conc_perturbs, rate_perturbs]
	return ode_dict, rate_dict, symbol_list, perturb_list

def _get_values(model, perturbations, symbol_list):
	if perturbations is None:
		conc_perturbs = dict()
		rate_perturbs = dict()
	else:
		conc_perturbs = perturbations[0]
		rate_perturbs = perturbations[1]
	metabolites = symbol_list[0]
	rate_params = symbol_list[1]
	fixed_concs = symbol_list[2]
	custom_params = symbol_list[3]

	# For rate parameters
	values = dict()
	for param_sym in rate_params:
		[p_type, rid] = re.split("\_", str(param_sym), maxsplit=1)
		# Use the perturbation value if it exists
		if "%s_%s" % (p_type, rid) in iterkeys(rate_perturbs):
			values.update({param_sym: rate_perturbs["%s_%s" % (p_type, rid)]})
		# Otherwise get the numerical value for the parameter
		else:
			reaction = model.reactions.get_by_id(rid)
			prop_f = reaction.__class__.__dict__[p_type]
			values.update({param_sym: prop_f.fget(reaction)})

	# For fixed concentrations
	for metab_sym in fixed_concs:
		metab = str(metab_sym)
		if ext_metab_re.search(metab):
			if metab in iterkeys(conc_perturbs):
				values.update({metab_sym : conc_perturbs[metab]})
			else:
				values.update({metab_sym : model.fixed_concentrations[metab]})
		else:
			metab = model.metabolites.get_by_id(metab)
			if metab in iterkeys(conc_perturbs):
				values.update({metab_sym : conc_perturbs[metab]})
			else:
				values.update({metab_sym :model.fixed_concentrations[metab]})

	# For custom_parameters
	for c_param in custom_params:
		print("FIXME: IMPLEMENT CUSTOM PARAMETERS")

	# Get initial conditions
	initial_conditions = list()
	for m in metabolites:
		metab = model.metabolites.get_by_id(str(m)[:-3])
		if metab in iterkeys(conc_perturbs):
			initial_conditions.append(conc_perturbs[metab])
		else:
			initial_conditions.append(model.initial_conditions[metab])

	return values, initial_conditions

def _make_lambda_odes(model, metabolites, ode_dict, values):
	# Get the metabolite matrix
	metabolites = tuple(metab for metab in list(metabolites))
	metab_matrix = sp.Matrix(metabolites)
	# Get the ODE equations
	eqs = {}
	for metab in metab_matrix:
		eqs[metab.diff(t)] = ode_dict[model.metabolites.get_by_id(
															str(metab)[:-3])]

	# Get ode equations and Jacobian matrix
	odes = metab_matrix.diff(t).subs(eqs)
	jacb = sp.Matrix([[fj.diff(m) for m in metab_matrix] for fj in odes])
	# Build lambda functions for odes and Jacobian
	lambda_odes = sp.lambdify((t, metab_matrix), odes.subs(values), "numpy")
	lambda_jacb = sp.lambdify((t, metab_matrix), jacb.subs(values), "numpy")

	return [lambda_odes, lambda_jacb]

def _integrate_odes(t_vector, lam_odes, lam_jacb, ics, solver):
	# Fix types from tuples to lists
	def f(t, y):
		res = lam_odes(t, y)
		return list(res)

	def j(t, y):
		res = lam_jacb(t, y)
		return list(res)

	# Set up integrator
	integrator = ode(f, j).set_initial_value(ics, t_vector[0])
	integrator.set_integrator(solver, nsteps=int(1000000))

	# Set up solutions array
	y = np.zeros((len(t_vector), len(ics)))
	dt = t_vector[1] - t_vector[0]
	idx = 0
	while integrator.successful() and integrator.t < t_vector[-1]:
		y[idx,:] = integrator.y
		integrator.integrate(integrator.t + dt)
		idx += 1

	return y

def _make_lambda_rates(model, metabolites, rate_dict, values):
	# Turn metabolite functions into metabolite symbols
	metab_syms = list(sp.Symbol(model.metabolites.get_by_id(
					str(m_func)[:-3]).id) for m_func in metabolites)
	metab_func_to_sym = {}

	metab_func_to_sym = dict((metab_func, metab_syms[i])
							for i, metab_func in enumerate(list(metabolites)))
	for rxn, rate in iteritems(rate_dict):
		rate = rate.subs(metab_func_to_sym).subs(values)
		args = tuple(sp.Symbol(m.id) for m in iterkeys(rxn.metabolites))
		rate_dict[rxn] = [args, sp.lambdify(args, rate)]
	return rate_dict
