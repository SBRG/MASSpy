# -*- coding: utf-8 -*-

# Compatibility with Python 2.7
from __future__ import absolute_import

# Import necesary packages
import re
import sympy as sp
from six import iteritems, iterkeys
# from cobra
from cobra.core.dictlist import DictList
# from mass
from mass.core import massmetabolite
from mass.core import massmodel
from mass.core import massreaction

# Class begins
## Global symbol for time
t = sp.Symbol("t")
## Public
def qcqa_model(model, missing_params=True, missing_ics=True,
		simulation=True, superflous=True,
		unconserved_metabolites=True, param_consistency=True,
		stoichiometry=True, elemental=True, thermodynamics=True):
	"""Run a series of quality control and assessment tests on a massmodel and
	return a summary of the test results

	Parameters
	----------
	model : mass.massmodel
		The MassModel to inspect
	missing_params : bool
		Check for missing parameters in the model and ensure the model
		parameters are consistent if there are superflous parameters
	missing_ics : bool
		Check for missing initial_conditions in the model
	simulation : bool
		Check to see if the model can be simulated
	superflous : bool
		Check for superflous parameters in the model
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
	check_list = [missing_params, missing_ics, simulation, superflous,
					unconserved_metabolites, param_consistency, stoichiometry,
					elemental, thermodynamics]

	name_list = ["missing_parameters", "missing_initial_conditions",
				"can_simulate","superflous","unconserved_metabolites"
				"param_consistency", "stoichiometry","elemental",
				"thermodynamics"]

	function_list =[get_missing_parameters,
					get_missing_initial_conditions,
					can_simulate,
					get_superflous_parameters,
					get_unconserved_metabolites,
					parameter_consistency,
					stoichiometric_consistency,
					elemental_consistency,
					thermodynamic_consistency]

	# Check inputs
	if not isinstance(model, massmodel.MassModel):
		raise TypeError("model must be a mass.MassModel")
	for i, check in enumerate(check_list):
		if not isinstance(check, bool):
			raise TypeError("%s must be a bool" % name_list[i])

	for i, check in enumerate(check_list):
		if check:
			check_results = function_list[i](model)


	print("\nFIXME: IMPLEMENT QCQA")
	return

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

		if len(missing_params) == 0:
			continue
		elif len(missing_params) == 1:
			missing_param_dict[rxn] = missing_params.pop()
		else:
			missing_param_dict[rxn] = missing_params

	return missing_param_dict

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

def can_simulate(model, show_missing=False):
	"""Check to see if the model has the required initial conditions and
	parameters in order to be simulated.

	Parameters
	----------
	model : mass.massmodel
		The MassModel to inspect
	show_missing : bool
		If True, will printout the missing parameters and conditions required
		for simulation.

	Returns
	-------
	True if the model has met conditions for simulation. Otherwise will return
	False, and if show_missing is True, will print the why False was returned.
	"""
	# Check inputs
	print("FIXME: IMPLEMENT CAN SIMULATE QCQA")
	return True

def get_superflous_parameters(model):
	"""Get extra parameters required for massmodel simulation. Superflous
	parameters are extra parameters that are not necessarily required for
	simulating the model. Will

	Primarily a concern when model.kr != (model.Keq / model.kf).

	Parameters
	----------
	model : mass.massmodel
		The MassModel to inspect
	"""
	# Check inputs
	print("FIXME: IMPLEMENT SUPERFLOUS")
	return

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
	constant if there are superflous parameters. If there are no missing or
	superflous parameters, parameters are considered consistent.

	Parameters
	----------
	model : mass.massmodel
		The MassModel to inspect
	"""
	# Check inputs
	print("FIXME: IMPLEMENT PARAM CONSISTENTCY")
	return

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
	print("FIXME: IMPLEMENT ELEM CONSISTENTCY")
	return

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
