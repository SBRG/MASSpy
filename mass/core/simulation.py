# -*- coding: utf-8 -*-

# Compatibility with Python 2.7
from __future__ import absolute_import

import re
import numpy as np
import sympy as sp
from warnings import warn
from scipy.integrate import odeint
from scipy.optimize import root
from scipy.interpolate import interp1d
from six import iteritems, iterkeys, itervalues, integer_types

# from mass
from mass.util import qcqa
from mass.core import expressions
from mass.core.massmodel import MassModel

# Class begins
## Global symbol for time
t = sp.Symbol("t")
## Precompiled re for 'external' metabolites
ext_metab_re = re.compile("\_e")
# Possible Perturbation Types
kf_re = re.compile("kf|forward_rate_constant")
Keq_re = re.compile("Keq|equilibrium_constant")
kr_re = re.compile("kr|reverse_rate_constant")
ic_re = re.compile("ic|initial_condition")
fixed_re = re.compile("fix|fixed")
function_re = re.compile("func|function")


# Public
def simulate(model, time, perturbations=None, numpoints=1000, nsteps=500,
				first_step=0., min_step=0., max_step=0.):
	"""Simulate a MassModel by integrating the ODEs  using the specified solver
	at the given time points and for given perturbation(s) to the model to
	obtain solutions for the metabolite concentrations and reaction fluxes.

	Perturbations can be one of the following types:
		'rxn.kf' or 'rxn.forward_rate_constant'
		'rxn.Keq' or 'rxn.equilibrium_constant'
		'rxn.kr' or 'rxn.reverse_rate_constant'
		'metab.ic' or 'metab.initial_condition'
		'metab.fix' or 'metab.fixed'

	Parameters
	----------
	model : mass.MassModel
		The MassModel object to simulate:
	time : tuple, list, or numpy.ndarray
		A tuple containing the start and end time points for the
		simulation, or a list of numerical values to treat as the time points
		for integration of the ODEs.
	numpoints :  int, optional
		The number of points to plot if the given time is a tuple.
		Default is 1000.
	perturbations : dict, optional
		A dictionary of events to incorporate into the simulation, where keys
		are the event to incorporate, and values are new parameter or initial
		condition. Can be changes to the rate and equilibrium constants,
		a change to an initial condition, or fixing a concentration.
	nsteps : int, optional
		Maximum number of (internally defined) steps allowed during
		one call to the solver.
	first_step :  float, optional
		The value for the first step used by the integrator
	min_step : float, optional
		The minimum allowable step size used by the integrator.
		Does not apply to dopri5 and dop853 solvers
	max_step : float, optional
		 Limits for the step sizes used by the integrator.

	Returns
	-------
	c_profile : dict
		A dictionary containing the concentration solutions, where key:value
		pairs are metabolites: scipy interpolating functions of solutions
	f_profile : dict
		A dictionary containing the flux solutions, where key:value
		pairs are reactions: scipy interpolating functions of solutions
	"""
	# Check inputs
	if not isinstance(model, MassModel):
		raise TypeError("model must be a MassModel")

	if isinstance(time, tuple):
		if not isinstance(numpoints, (float, integer_types)):
			raise TypeError("numpoints must an integer")
		if abs(time[0]) < 1e-9:
			time = (1e-6, time[1])
		time = np.geomspace(time[0], time[1], int(numpoints), endpoint=True)

	if not isinstance(time, (np.ndarray, list)):
		raise TypeError("time must be a list or numpy.ndarray of time "
						"points, or a tuple of form (start_point, end_point).")
	if perturbations is None:
		perturbations = {}
	elif not isinstance(perturbations, dict):
		raise TypeError("Perturbations must be in a dictionary")
	else:
		full_pert_check = re.compile("|".join([pert_type.pattern
					for pert_type in [kf_re, Keq_re, kr_re,
										ic_re, fixed_re, function_re]]))
		for perturb, value in iteritems(perturbations):
			if not full_pert_check.search(perturb):
				raise TypeError("Perturbation not recognized")

	if not isinstance(nsteps, int):
		raise TypeError("nsteps must be an integer")
	if not isinstance(first_step, float):
		raise TypeError("first_step must be an integer")
	if not isinstance(min_step, float):
		raise TypeError("min_step must be an integer")
	if not isinstance(max_step, float):
		raise TypeError("max_step must be an integer")

	sim_check = qcqa.can_simulate(model)
	if not sim_check:
		warn("Unable to simulate")
		qcqa.qcqa_model(model, initial_conditions=True, parameters=True,
							simulation=True)
		return [None, None]
	model.repair()
	# Collect sympy symbols and make dictionaries for odes and rates
	odes, rates, symbols = expressions._sort_symbols(model)
	# Perturb the system if perturbations exist
	if len(perturbations) != 0:
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
	time, c = _integrate_odes(time, lam_odes, lam_jacb, ics,
								nsteps, first_step, min_step, max_step)
	# Map metbaolite ids to their concentration solutions
	c_profile = dict()
	for i, sym in enumerate(metab_syms):
		c_profile.update({model.metabolites.get_by_id(str(sym)[:-3]) : c[:,i]})
	for i, sym in enumerate(fixed_syms):
		if not ext_metab_re.search(str(sym)):
			metab = model.metabolites.get_by_id(str(sym))
			c_profile.update({metab: values[sym]})

	# Use the concentrations to get the flux the flux solutions
	# Map reactiom ids to their flux solutions
	f_profile = dict()
	for rxn, lambda_func_and_args in iteritems(lam_rates):
		f = np.zeros(time.shape)
		lambda_func = lambda_func_and_args[1]
		args = lambda_func_and_args[0]
		concs = np.array([c_profile[model.metabolites.get_by_id(str(arg))]
							for arg in args if arg is not sp.Symbol("t")]).T
		for i in range(0, len(f)):
			if len(args) != 0 and args[-1] != sp.Symbol("t"):
				f[i] = lambda_func(*concs[i,:])
			elif len(args) != 0 and args[-1] == sp.Symbol("t"):
				f[i] = lambda_func(*concs[i,:], time[i])
			else:
				f[i] = lambda_func()
		f_profile[rxn] = f

	# Create interpolating functions for the concentrations and flux profiles
	for profile in [c_profile, f_profile]:
		for key, sol in iteritems(profile):
			if abs(sol[-1]) <= 1e-9:
				sol = sol[:-1]
			if time.shape != sol.shape:
				t = time[:-1]
			else:
				t = time
			profile[key] = interp1d(t, sol, kind='cubic',
										fill_value='extrapolate')

	# Reorder dictionary to match model metabolites and reactions.
	new_profile = {}
	for metab in model.metabolites:
		new_profile[metab] = c_profile[metab]
	c_profile = new_profile
	new_profile = {}
	for rxn in model.reactions:
		new_profile[rxn] = f_profile[rxn]
	f_profile = new_profile

	return [c_profile, f_profile]

def find_steady_state(model, strategy="simulate", perturbations=None,
			update_reactions=False, update_initial_conditions=False):
	"""Find the steady state solution of a model using a given strategy

	Parameters
	----------
	model : mass.MassModel
		The MassModel object to find a steady state for.
	strategy : 'simulate' or 'find_roots'
		The strategy to use to solve for the steady state.
	perturbations : dict, optional
		A dictionary of events to incorporate into the simulation, where keys
		are the event to incorporate, and values are new parameter or initial
		condition. Can be changes to the rate and equilibrium constants,
		a change to an initial condition, or fixing a concentration.
	update_reactions : bool, optional
		If True, update the steady state fluxes (ssflux) in each reaction.
	update_initial_conditions : bool, optional
		If True, update initial conditions in the model for each metabolite.

	Returns
	-------
	list of dictionaries representing the concentration and flux solutions:
		The first item in the list is a dictionary containing the concentration
		solutions where the key:value pairs are metabolite objects and their
		corresponding steady state solutions.
		The second item in the list is a dictionary containing the flux
		solutions where the key:value pairs are reactionobjects and their
		corresponding steady state solutions
	"""
	# Check input
	if not isinstance(model, MassModel):
		raise TypeError("model must be a MassModel")

	sim_check = qcqa.can_simulate(model)
	if not sim_check:
		warn("Unable to find steady state due to missing values")
		qcqa.qcqa_model(model, initial_conditions=True, parameters=True,
							simulation=True)
		return [None, None]

	if perturbations is None:
		perturbations = {}
	elif not isinstance(perturbations, dict):
		raise TypeError("Perturbations must be in a dictionary")
	else:
		full_pert_check = re.compile("|".join([pert_type.pattern
					for pert_type in [kf_re, Keq_re, kr_re,
										ic_re, fixed_re, function_re]]))
		for perturb, value in iteritems(perturbations):
			if not full_pert_check.search(perturb):
				raise TypeError("Perturbation not recognized")

	options = {"simulate": simulate, "find_roots": None}
	# Perform the simulate strategy
	if strategy is "simulate":
		# Start with final time point at 10^3, quit after trying 10^6
		power = 1
		fail_power = 6
		while power <= fail_power:
			retry = False
			time = np.geomspace(1e-6, 10**power, int(10**power), endpoint=True)
			[c_profile, f_profile] = options[strategy](
									model, time,
									perturbations=perturbations)
			for metab, profile in iteritems(c_profile):
				conc = profile(time)
				if abs(conc[-1] - conc[-2]) <= 10**-9:
					continue
				else:
					retry = True
			if retry:
				power += 1
			else:
				break
		if power > fail_power:
			warn("Unable to find a steady state using strategy %s" % strategy)
			return [None, None]
		# Return steady state solutions
		for metab, profile in iteritems(c_profile):
			conc = float(profile(time[-1]))
			c_profile[metab] = round(conc, 6)
			# Update model initial conditions if specified
			if update_initial_conditions:
				model.initial_conditions[metab] = round(conc, 6)
		for reaction, profile in iteritems(f_profile):
			flux = float(profile(time[-1]))
			f_profile[reaction] = round(flux, 6)
			# Update reaction steady state flux if specified
			if update_reactions:
				reaction.ssflux = round(flux, 6)

		return [c_profile, f_profile]

	# Perform the find_roots strategy
	elif strategy is "find_roots":
		# Collect sympy symbols and make dictionariess for odes and rates
		odes, rates, symbols = expressions._sort_symbols(model)
		if len(perturbations) != 0:
			odes, rates, symbols, perturbations = _perturb(model, odes, rates,
														symbols, perturbations)
		# Get values to substitute into ODEs and metabolite initial conditions
		values, ics = _get_values(model, dict(), symbols)

		metab_syms = tuple(sp.Symbol(str(metab)[:-3])
							for metab in list(symbols[0]))
		# Strip time dependency from odes
		odes = expressions.strip_time(odes)
		eqs = dict()
		for metab in metab_syms:
			ode = odes[model.metabolites.get_by_id(str(metab))]
			eqs[metab] = ode.subs(values)
		# Create the lambda function for root finding
		lam_f = sp.lambdify(list(iterkeys(eqs)), list(itervalues(eqs)),"numpy")

		# Create a callable function for root interface
		def f(ics):
			res = lam_f(*ics)
			return res
		# Get roots
		sol = root(f, ics, method='krylov')
		# Return a warning if no roots are found
		if not sol.success:
			warn("Unable to find a steady state using strategy %s" % strategy)
			return [None, None]

		# Map concentration
		c_profile = dict()
		for i, metab in enumerate(iterkeys(eqs)):
			metab = model.metabolites.get_by_id(str(metab))
			c_profile[metab] = round(sol.x[i], 6)
			# Update model initial conditions if specified
			if update_initial_conditions:
				model.initial_conditions[metab] = round(sol.x[i], 6)

		# Make lambda functions for rates
		lam_rates = _make_lambda_rates(model, symbols[0], rates, values)
		# Use the concentrations to get the flux the flux solutions
		# Map reactiom ids to their flux solutions
		f_profile = dict()
		for rxn, lambda_func_and_args in iteritems(lam_rates):
			lambda_func = lambda_func_and_args[1]
			concs = np.array([c_profile[model.metabolites.get_by_id(str(arg))]
								for arg in lambda_func_and_args[0]]).T
			f_profile[rxn] = lambda_func(*concs)
			# Update reaction steady state flux if specified
			if update_reactions:
				rxn.ssflux = f_profile[rxn]

		return [c_profile, f_profile]

	else:
		raise ValueError("Unrecognized strategy. Strategy must be "
						"'simulate' or 'find_roots'.")

# Internal
def _perturb(model, ode_dict, rate_dict, symbol_list, perturbations):
	"""Internal use. Apply the perturbations to the ODEs and rate laws"""
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
			rate_perturbs.update({pert : value})
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
			if not ext_metab_re.search(item_id):
				# Other metabolites
				metab = model.metabolites.get_by_id(item_id)
				# Use a copy to enable removal of metabolite function
				# from the metabolite set in the iteration
				for metab_func in metabolites.copy():
					if(re.match(metab.id, str(metab_func)[:-3])):
						metab_sym = sp.Symbol(metab.id)
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
			else:
				try:
					metab = model.metabolites.get_by_id(item_id)
				except KeyError:
					metab = item_id
				finally:
					conc_perturbs.update({metab: value})
		# Handle initial condition perturbations
		if ic_re.match(to_perturb):
			if item_id not in model.metabolites:
				raise KeyError("%s is not a MassMetabolite object found "
								"in the model.metabolites DictList." % item_id)
			conc_perturbs.update(
							{model.metabolites.get_by_id(item_id): value})
			continue
		if function_re.match(to_perturb):
			try:
				metab = model.metabolites.get_by_id(item_id)
			except KeyError:
				metab = item_id
			finally:
				conc_perturbs.update({metab: sp.sympify(value)})
	symbol_list = [metabolites, rate_params, fixed_concs, custom_params]
	perturb_list = [conc_perturbs, rate_perturbs]
	return ode_dict, rate_dict, symbol_list, perturb_list

def _get_values(model, perturbations, symbol_list):
	"""Internal use. Obtain the numerical values for sympy symbols and
	create a dictionary for of those values for subsitution"""
	if len(perturbations) == 0:
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
		if str(c_param) in iterkeys(rate_perturbs):
			values.update({c_param: rate_perturbs[str(c_param)]})
		else:
			values.update({c_param: model.custom_parameters[str(c_param)]})

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
	"""Internal use. Make the system of ODEs and its jacobian into a
	lambda function for integration"""
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
	lambda_odes = sp.lambdify((metab_matrix, t), odes.subs(values), "numpy")
	lambda_jacb = sp.lambdify((metab_matrix, t), jacb.subs(values), "numpy")

	return [lambda_odes, lambda_jacb]

def _integrate_odes(t_vector, lam_odes, lam_jacb, ics, nsteps, first_step,
					min_step, max_step):
	"""Internal use. Integrate the ODEs using lambda functions that represent
	the system of ODEs and the jacobian"""
	# Fix types from tuples to lists

	def f(y, t):
		res = np.array(lam_odes(y, t)).T[0]
		return list(res)

	def j(y, t):
		res = lam_jacb(y, t)
		return list(res)

	y = odeint(f, ics, t_vector, Dfun=j, h0=first_step, mxstep=nsteps,
											hmax=max_step, hmin=min_step)
	return t_vector, y

def _make_lambda_rates(model, metabolites, rate_dict, values):
	"""Make a lambda function that uses the concentration solutions to
	calculate the flux values"""
	# Turn metabolite functions into metabolite symbols
	metab_objects = list(model.metabolites.get_by_id(str(m_func)[:-3])
						for m_func in metabolites)
	metab_syms = list(sp.Symbol(m.id) for m in metab_objects)

	metab_func_to_sym = dict((metab_func, metab_syms[i])
							for i, metab_func in enumerate(list(metabolites)))
	for rxn, rate in iteritems(rate_dict):
		rate = rate.subs(metab_func_to_sym).subs(values)
		args = tuple(sp.Symbol(m.id) for m in metab_objects
					if sp.Symbol(m.id) in rate.atoms(sp.Symbol))
		if t in rate.atoms(sp.Symbol):
			args += tuple([t])
		rate_dict[rxn] = [args, sp.lambdify(args, rate, "numpy")]
	return rate_dict
