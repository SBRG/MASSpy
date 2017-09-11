# -*- coding: utf-8 -*-

# Compatibility with Python 2.7
from __future__ import absolute_import

# Import necesary packages
import re
import sympy as sp
from math import inf
from tabulate import tabulate
from six import iteritems, iterkeys

# from mass
from mass.core import massmetabolite
from mass.core import massmodel
from mass.core import massreaction

# Class begins
## Global symbol for time
t = sp.Symbol("t")
## Public
def qcqa_model(model, initial_conditions=False, parameters=False,
			simulation=False, superfluous=False, unconserved_metabolites=False,
			param_consistency=False, stoichiometry=False, elemental=False,
			thermodynamics=False):
	"""Run a series of quality control and assessment tests on a massmodel and
	return a summary of the test results

	Parameters
	----------
	model : mass.massmodel
		The MassModel to inspect
	initial_conditions : bool
		Check for missing initial_conditions in the model
	parameters : bool
		Check for missing parameters in the model and ensure the model
		parameters are consistent if there are superfluous parameters
	simulation : bool
		Check to see if the model can be simulated
	superfluous : bool
		Check for superfluous parameters in the model
	unconserved_metabolites : bool
		Check for unconserved metabolites in the model
	param_consistency : bool
		Check for parameter consistency in the model
	stoichiometry : bool
		Check for stoichiometric consistency in the model
	elemental : bool
		Check for elemental consistency in the model. Ignores the exchanges.
	thermodynamics : bool
		Check for thermodynamic consistency in the model.
	"""
	# List of bools indicating what QCQA functions to perform
	check_list = [initial_conditions, parameters, simulation, superfluous,
				unconserved_metabolites, param_consistency, stoichiometry,
				elemental, thermodynamics]

	# Names of the inputs
	name_list = ["initial_conditions", "parameters", "simulation",
				"superfluous","unconserved_metabolites", "param_consistency",
				"stoichiometry","elemental", "thermodynamics"]
	# The functions to perform for each check, and associated arguments if any.
	function_and_args =[[get_missing_initial_conditions, None],
						[get_missing_parameters, [True]*5],
						[can_simulate, [[1,2,3]], False],
						[get_superfluous_parameters, None],
						[get_unconserved_metabolites, None],
						[parameter_consistency, None],
						[stoichiometric_consistency, None],
						[elemental_consistency, None],
						[thermodynamic_consistency, None]]

	# Check inputs
	if not isinstance(model, massmodel.MassModel):
		raise TypeError("model must be a mass.MassModel")
	for i, check in enumerate(check_list):
		if not isinstance(check, bool):
			raise TypeError("%s must be a bool" % name_list[i])

	# Create a list of results to display in the QCQA report
	to_display = []
	for i, check in enumerate(check_list):
		# If the check is set to True, perform the function and
		# pass any existing arguments.
		if check:
			function = function_and_args[i][0]
			args = function_and_args[i][1]
			if args is not None:
				to_display.append(function(model, *args))
			else:
				to_display.append(function(model))
		# Set result to None if function check is not enabled.
		else:
			to_display.append(None)
	# Print report
	reports = _qcqa_summary(to_display)
	for report in reports:
		print("%s\n" % report)

def get_missing_initial_conditions(model):
	"""Get the initial conditions that are missing from the MassModel for the
	metabolites that exist in the Massmodel.

	Parameters
	----------
	model_or_metabolite_list : mass.massmodel
		The MassModel to inspect
	"""
	# Check inputs
	if not isinstance(model, massmodel.MassModel):
		raise TypeError("model must be a mass.MassModel")

	return [metab for metab in model.metabolites
			if metab not in iterkeys(model.initial_conditions)]

def get_missing_parameters(model, kf=False, Keq=False,
							kr=False, ssflux=False, custom_parameters=False):
	"""Get the parameters that are missing from the reactions that exist
	in the Massmodel.

	Parameters
	----------
	model : mass.Massmodel
		The MassModel or list of reactions to inspect.
	kf : bool
		If True, check MassReactions for missing forward rate constants.
	Keq : bool
		If True, check MassReactions for missing equilibrium rate constants.
	kr : bool
		If True, check MassReactions for missing reverse rate constants.
	ssflux : bool
		If True, check MassReactions for missing steady state fluxes.

	Returns
	-------
	dictonary where keys are reactions and values are a list of
		the missing parameters for the reaction.
	"""
	param_checks = [kf, Keq, kr, ssflux]
	param_keys = ["kf", "Keq", "kr", "ssflux"]
	attr_list = ["_forward_rate_constant", "_equilibrium_constant",
				"_reverse_rate_constant", "ssflux"]

	# Check inputs
	if not isinstance(model, massmodel.MassModel):
		raise TypeError("model must be a mass.MassModel")
	for i, param in enumerate(param_checks):
		if not isinstance(param, bool):
			raise TypeError("%s must be a bool" % param_keys[i])
		#If input type is correct, set parameters as their symbolic attributes
		if param_keys[i] != "ssflux":
			param_keys[i] = ("_sym_%s" % param_keys[i])

	missing_param_dict = dict()
	for rxn in model.reactions:
		missing_params = list()
		# If the reaction has custom rates and check is set to True
		if rxn in iterkeys(model.custom_rates) and custom_parameters:
			symbols = model.custom_rates[rxn].atoms(sp.Symbol)
			for sym in symbols:
				# Ignore time symbol
				if sym is t:
					continue
				# If the symbol is in the custom parameters, check the value
				elif str(sym) in iterkeys(model.custom_parameters):
					if model.custom_parameters[str(sym)] is not None:
						continue
					else:
						# Add to missing parameters if value is None
						missing_params.append(str(sym))
				elif re.search("kf|Keq|kr", str(sym)) and \
						re.split("\_", str(sym), maxsplit=1):
						p_type= re.split("\_", str(sym), maxsplit=1)[0]
						prop_f = rxn.__class__.__dict__[p_type]
						if prop_f.fget(rxn) is not None or \
							str(sym) in missing_params:
							continue
						# Add if value is none and not already added
						else:
							missing_params.append(str(sym))
				# Add to missing parameters if not found anywhere
				else:
					missing_params.append(str(sym))
		# If no custom rate, search for parameters from parameter checks
		else:
			for i, param_check in enumerate(param_checks):
				# Move on to next parameter if set to False
				if not param_check:
					# Next parameter if set to False
					continue
				param = rxn.__dict__[attr_list[i]]
				param_name = rxn.__dict__[param_keys[i]]
				if param is not None:
					continue
				elif param_keys[i] != "ssflux":
					missing_params.append(param_name)
				else:
					missing_params.append(param_keys[i])

		if len(missing_params) == 0:
			continue
		else:
			missing_param_dict[rxn] = missing_params
	return missing_param_dict

def can_simulate(model, rate_type=None):
	"""Check to see if the model has the required initial conditions and
	parameters in order to be simulated with the given rate law type(s)

	Parameters
	----------
	model : mass.massmodel
		The MassModel to inspect
	rate_type :  1,2,3, list of types, or None
		What rate type(s) to check for the ability to simulate.
		If None, will use the model's current rate type.

	Returns
	-------
	A dictionary where keys are the rate types and values are bools indicating
		if the model has met conditions for simulation.
	"""
	# Check inputs
	if not isinstance(model, massmodel.MassModel):
		raise TypeError("model must be a mass.MassModel")
	if rate_type is None:
		rate_type = [model._rtype]
	if not isinstance(rate_type, list):
		rate_type = [rate_type]

	for i, rt in enumerate(rate_type):
		rt = int(rt)
		if rt not in {1,2,3}:
			raise TypeError("Rate type must be integers 1,2, or 3")
		else:
			rate_type[i] = rt

	# Check for missing initial conditions
	missing_ics = get_missing_initial_conditions(model)

	# Check for missing parameters
	simulate_checks = {}
	for rt in rate_type:
		if rt == 1:
			missing_params = get_missing_parameters(model, kf=True, Keq=True,
								kr=False, ssflux=False, custom_parameters=True)
		elif rt == 2:
			missing_params = get_missing_parameters(model, kf=True, Keq=False,
								kr=True, ssflux=False, custom_parameters=True)

		else:
			missing_params = get_missing_parameters(model, kf=False, Keq=True,
								kr=True, ssflux=False, custom_parameters=True)

		# Set check results based on missing initial conditions and parameters
		if len(missing_params) != 0 or len(missing_ics) != 0:
			simulate_checks.update({rt: False})
		else:
			simulate_checks.update({rt: True})

	return simulate_checks

def get_superfluous_parameters(model):
	"""Get extra parameters required for massmodel simulation. superfluous
	parameters are extra parameters that are not necessarily required for
	simulating the model.

	Primarily a concern when model.kr != (model.Keq / model.kf).

	Parameters
	----------
	model : mass.massmodel
		The MassModel to inspect
	"""
	# Check inputs
	if not isinstance(model, massmodel.MassModel):
		raise TypeError("model must be a mass.MassModel")

	# Get superfluous parameters
	superfluous_parameters = {}
	for rxn in model.reactions:
		if len(rxn.parameters) == 3:
			superfluous_parameters[rxn] = [rxn._sym_kr]

	return superfluous_parameters

def get_unconserved_metabolites(model):
	"""Get the unconserved metabolites in a massmodel

	Parameters
	----------
	model : mass.massmodel
		The MassModel to inspect
	"""
	# Check inputs
	print("FIXME: IMPLEMENT UNCONSERVED METABS")
	return

def parameter_consistency(model):
	"""Performs a consistency check on the rate constants and the equilibrium
	constant if there are superfluous parameters. If there are no missing or
	superfluous parameters, parameters are considered consistent.

	Parameters
	----------
	model : mass.massmodel
		The MassModel to inspect
	"""
	# Check inputs
	if not isinstance(model, massmodel.MassModel):
		raise TypeError("model must be a mass.MassModel")

	# Get superfluous parameters
	param_consistency = {}
	superfluous_parameters = get_superfluous_parameters(model)
	for rxn, superfluous in iteritems(superfluous_parameters):
		check = (rxn.kr == rxn.kf/rxn.Keq)
		if check:
			param_consistency[rxn] = ("%s: kf/Keq == kr" % check)
		else:
			param_consistency[rxn] = ("%s: kf/Keq != kr" % check)

	return param_consistency

def elemental_consistency(model):
	"""Performs a consistency check on the reactions in the model to ensure
	they are mass balanced and charged balanced.

	Exchange reactions are ignored because they are not typically balanced

	Parameters
	----------
	model : mass.massmodel
		The MassModel to inspect
	"""
	# Check inputs
	if not isinstance(model, massmodel.MassModel):
		raise TypeError("model must be a mass.MassModel")

	# Check for elemental consistency
	elem_consistency = {}
	for rxn in model.reactions:
		if not rxn.exchange and rxn.check_mass_balance() != {}:
			unbalanced = ""
			for elem, amount in iteritems(rxn.check_mass_balance()):
				unbalanced += "%s: %.1f;" % (elem, amount)
			elem_consistency[rxn] = unbalanced.rstrip("; ") + " unbalanced"

	return elem_consistency

def stoichiometric_consistency(model):
	"""Performs a consistency check on the stoichiometry of the model.

	Parameters
	----------
	model : mass.massmodel
		The MassModel to inspect
	"""
	# Check inputs
	print("FIXME: IMPLEMENT S CONSISTENTCY")
	return

def thermodynamic_consistency(model):
	"""Performs a consistency check on the thermodynamics of the model.

	Parameters
	----------
	model : mass.massmodel
		The MassModel to inspect
	"""
	# Check inputs
	print("FIXME: IMPLEMENT THERMO CONSISTENTCY")
	return

## Internal
def _qcqa_summary(to_display):
	"""Internal use. Create reports to print out based on the QCQA results.
	Returns a list of tabulated reports for printing.

	Parameters
	----------
	to_display: list
		A list of the QCQA results to report
	"""
	name_list = ["Missing Initial Conditions", "Missing Parameters",
			"Can Simulate","Superfluous Parameters","Unconserved Metabolites",
			"Parameter Consistency", "Stoichiometric Consistency",
			"Elemental Consistency", "Thermodynamic Consistency"]

	headers = []
	item_list = []
	for i, display_item in enumerate(to_display):
		if display_item is not None:
			headers.append(name_list[i])
			item_list.append(display_item)

	missing_dict = {}
	sim_checks = {}
	consistencies = {}
	reports = []
	for i, item in enumerate(item_list):
		header = headers[i]
		# Set up printout for missing initial conditions
		if re.match("Missing Initial Conditions", header):
			if len(item) != 0:
				missing_dict[header] = item
			continue
		# Set up printout for missing parameters
		if re.match("Missing Parameters", header) or \
			re.match("Superfluous Parameters", header):
			table_list = []
			for rxn, missing in iteritems(item):
				missing_params = ": "
				for param in missing:
					if re.split("\_", param, maxsplit=1):
						missing_params += re.split("\_", param)[0] + "; "
					else:
						missing_params += "; "
				missing_params = "%s%s" % (rxn.id, missing_params.rstrip("; "))
				table_list.append(missing_params)

			if len(table_list) != 0:
				missing_dict[header] = table_list
			continue

		if re.match("Can Simulate", header):
			table_list = []
			for rate_type, sim_check in iteritems(item):
				table_list += ["Rate Type %s: %s" % (rate_type, sim_check)]
			sim_checks[header] = table_list
			continue

		if re.match("Parameter Consistency", header) or \
			re.match("Elemental Consistency", header):
			table_list = []
			if len(item) != 0:
				for rxn, consistency in iteritems(item):
					table_list += ["%s: %s" % (rxn.id, consistency)]
				consistencies[header] = table_list
			continue

	reports.append("QCQA REPORT\n" + "="*79)
	if sim_checks != {}:
		reports.append(tabulate(sim_checks, headers="keys",stralign="center"))
	if missing_dict != {}:
		reports.append(tabulate(missing_dict, headers="keys",stralign="left"))
	if consistencies != {}:
		reports.append(tabulate(consistencies, headers="keys",stralign="left"))

	return reports
