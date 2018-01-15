# -*- coding: utf-8 -*-

# Compatibility with Python 2.7
from __future__ import absolute_import

import numpy as np
import sympy as sp
from scipy import interpolate
from six import string_types, integer_types, iterkeys, iteritems

# from cobra
from cobra import DictList

# Class begins
## Public
def pools_from_string(concentration_profile, time_range, pools,
						parameters=None, pool_ids=None, numpoints=5000):
	"""Create a dictionary of interpolating functions for a list of pools
	defined by string input.

	Example: For the reaction v1: x1 <=> x2 with Keq = 2,  a conservation pool
	and a disequilibrium pool can be defined by providing the following input
	for pools and parameters:

	pools = ['x1 + x2', 'x1 - x2/Keq_v1']
	parameters = {'Keq_v1' : 2}

	Parameters
	----------
	concentration_profile : dict
		A dictionary containing the concentration or flux solutions, where
		keys are MassMetabolite or MassReaction objects, and values are the
		scipy interpolating functions of solutions
	time_range : tuple or list
		Either a tuple containing the start and end time points for the
		simulation, or a list of numerical values to treat as the time points.
		Extrapolation will be utilized for when the interpolating functions
		recieves a value outside this time range after its creation.
	pools : string or list
		Either a string or a list of strings defining the pooling formula. All
		metabolites to be pooled must exist in the concentration profile
		dictionary.
	parameters : dictionary, optional
		A dictionary of aditional parameters to be used in the pools, where the
		key:value pairs are the parameter identifiers and its numerical value.
	pool_ids : string or list, optional
		String identifiers to use for the pools. The number of identifiers must
		match the number of pools. If None, will use default identifiers of
		'p1', 'p2', etc.
	numpoints :  int, optional
		The number of time points to use to create the interpolating function
		for the pools.
		Default is 500.

	Returns
	-------
	pool_dict : dict
		A dictionary containing the pooling interpolating functions, where
		key:value pairs are pool identiiers: scipy interpolating functions
		of the pools
	"""
	pool_dict = _setup_interpol_func(concentration_profile, time_range, pools,
								parameters, pool_ids, numpoints, "pool")
	return pool_dict

def net_fluxes_from_strings(flux_profile, time_range, net_fluxes,
				parameters=None, net_flux_ids=None, numpoints=5000):
	"""Create a dictionary of interpolating functions for a list of net_fluxes
	defined by string input.

	Example: To sum the fluxes of v1 and v2:
	net_fluxes = ['v1 + v2']

	Parameters
	----------
	flux_profile : dict
		A dictionary containing the flux or flux solutions, where
		keys are MassMetabolite or MassReaction objects, and values are the
		scipy interpolating functions of solutions
	time_range : tuple or list
		Either a tuple containing the start and end time points for the
		simulation, or a list of numerical values to treat as the time points.
		Extrapolation will be utilized for when the interpolating functions
		recieves a value outside this time range after its creation.
	net_fluxes : string or list
		Either a string or a list of strings defining the net flux formula. All
		fluxes to be grouped must exist in the flux profile
		dictionary.
	parameters : dictionary, optional
		A dictionary of aditional parameters to be used in the net_fluxes,
		where the key:value pairs are the parameter identifiers and its
		numerical value.
	net_flux_ids : string or list, optional
		String identifiers to use for the net_fluxes. The number of identifiers
		mustmatch the number of net_fluxes. If None, will use default
		identifiers of 'v_net1', 'v_net2', etc.
	numpoints :  int, optional
		The number of time points to use to create the interpolating function
		for the net_fluxes.
		Default is 500.

	Returns
	-------
	net_flux_dict: dict
		A dictionary containing the net flux interpolating functions, where
		key:value pairs are net flux identiiers: scipy interpolating functions
		of the net fluxes
	"""
	net_flux_dict = _setup_interpol_func(flux_profile, time_range, net_fluxes,
							parameters, net_flux_ids, numpoints, "net_flux")
	return net_flux_dict

def pairwise_angles(modal_matrix, correlation_cutoff=0.85):
	"""Obtain pooling matrices for each mode by calculating the angle
	between columns for each column of the modal matrix. The algorithm
	calculates of the columns for the matrix as a function of index k, which
	runs from 1 to n timescales, where the the first k row(s) of the modal
	matrix are zeroed out.

	Angles between two zero vectors are considered as undefined and angles
	between a zero vector and a vector with at least one non-zero element is
	defined at 90 degrees.

	For additional information on how the algorithm works:
	https://doi.org/10.1371/journal.pcbi.1000177

	Parameters
	----------
	modal_matrix :  numpy.array
		A numpy.array representing the modal matrix rank ordered from
		fastest mode to slowest mode.
	correlation_cutoff : float, optional
		A value between 0 and 1. that determines the cutoff for correlation.
		Any cos(angle) < correlation_cutoff will be set at zero.

	Returns
	-------
	p_matrices : numpy.array of length n
		A numpy.array of length n, where the matrix in p_matrices[k]
		corresponds to the kth mode of the modal matrix (modal_matrix[k])
		Each matrix entry represents the value given by cos(angle).
	"""

	# Make a copy of the modal matrix to prevent edits on the original
	M = modal_matrix.copy()
	# Get dimensions of the modal matrix
	n_ts, n_met = M.shape

	p_matrices = []
	for k in range(n_ts):
		if k > 0:
			for i in range(0, k):
				M[i] = np.zeros(M[i].shape)
		angle_mat = np.zeros((n_met, n_met))
		for i in range(n_met):
			for j in range(n_met):
				angle_mat[i, j] = abs(np.dot(M[:, i], M[:, j]) / (
										np.linalg.norm(M[:, i]) *
										np.linalg.norm(M[:, j])))
		F1 = np.zeros(angle_mat.shape)
		F2 = np.zeros(angle_mat.shape)
		for i, row in enumerate(angle_mat):
			for j, val in enumerate(row):
				if i == j and i <= k:
					F1[i, j] = 1
				if val >=correlation_cutoff and i != j:
					F2[i, j] = 1
				else:
					F2[i, j] = 0
		p_matrices.append(np.dot(F1, angle_mat, F2))
	return p_matrices

# Internal
def _setup_interpol_func(solution_profile, time_range, to_create,
						parameters, new_ids, numpoints, func_type):
	# Check inputs
	if not isinstance(solution_profile, dict):
		raise TypeError("concentration_profile must be a dictionary")
	elif len(solution_profile) == 0:
		warn ("No solutions associated with this profile")
		return None

	if isinstance(time_range, tuple):
		if not isinstance(numpoints, (float, integer_types)):
			raise TypeError("numpoints must an integer")
		if abs(time_range[0]) < 1e-9:
			time_range = (1e-6, time_range[1])
		time_range = np.geomspace(time_range[0], time_range[1], int(numpoints))

	if not hasattr(time_range, '__iter__'):
		raise TypeError("time_range must an iterable list of time points or "
						" a tuple of containing start and end points")

	if isinstance(to_create, string_types):
		to_create = [to_create]
	if not isinstance(to_create, list):
		raise TypeError("%s must be a string or list of strings" % func_type)

	if parameters is None:
		parameters = {}
	if not isinstance(parameters, dict):
		raise TypeError("parameters must be a dictionary")

	if isinstance(new_ids, string_types):
		new_ids = [new_ids]
	if not isinstance(new_ids, (list, type(None))):
		raise TypeError("%s_ids must be a string or a list of strings with "
					"the same length as %ss" % (func_type, func_type + "e"))


	else:
		# Generate pool IDs if not provided
		if new_ids is None and func_type is "pool":
			new_ids = ['p%s' % str(i+1) for i in range(0, len(to_create))]
		elif new_ids is None and func_type is "net_flux":
			new_ids = ['v_net%s' % str(i+1) for i in range(0, len(to_create))]
		# Otherwise check length of IDs if provided
		elif len(new_ids) == len(to_create):
			pass
		else:
			raise ValueError("Number of %s IDs does not match the number of "
							"defined %ss" % (func_type, func_type + "e"))

	# Create pools expressions using Sympy
	to_create_dict = dict()
	symbol_dict = dict()
	lookup = DictList()
	# Get metabolite objects in a DictList and local symbols in a dictionary
	for key in iterkeys(solution_profile):
		lookup.append(key)
		symbol_dict.update({key.id: sp.Symbol(key.id)})
	# Add parameters into local symbol dictionary, if any
	for key in iterkeys(parameters):
		symbol_dict.update({key: sp.Symbol(key)})
	for i, id_str in enumerate(to_create):
		# Use sympy.sympify to create expression
		try:
			expr = sp.sympify(id_str, locals=symbol_dict)
			expr = expr.subs(dict((sp.Symbol(k), v)
							for k, v in iteritems(parameters)))
		except:
			raise ValueError("Could not convert (%s) into a %s."
							% (id_str, func_type))
		# Set up substitution for concentration values
		sub_dict = dict((sym, solution_profile[
						lookup.get_by_id(str(sym))])
						for sym in expr.atoms(sp.Symbol))
		# Obtain solutions by substituting values for each time point
		sol = list()
		for t in time_range:
			val_at_t = dict((sym, sub_dict[sym](t))
							for sym, func in iteritems(sub_dict))
			sol.extend([expr.subs(val_at_t)])
		# Remove the last point if needed
		if sol[-1] == 0:
			sol = sol[:-1]
			time_vec = time_range[:-1]
		else:
			time_vec = time_range
		# Add the new pool interpolating funtion into pool_dict
		to_create_dict.update({new_ids[i]: interpolate.interp1d(time_vec,
							sol, kind='cubic', fill_value='extrapolate')})
	return to_create_dict
