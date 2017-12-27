# -*- coding: utf-8 -*-

# Compatibility with Python 2.7
from __future__ import absolute_import

import numpy as np
import sympy as sp
import scipy as sc
from six import string_types, integer_types, iterkeys, iteritems
from cobra import DictList

def pools_from_string(concentration_profile, time_range, pools,
						parameters=None, pool_ids=None, numpoints=5000):
	"""Create a dictionary of interpolating functions for a list of pools
	defined by string input.

	Example: For the reaction v1: x1 <=> x2 with Keq = 2,  a conservation pool and a
	disequilibrium pool can be defined by providing the following input for
	pools and parameters:

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
	pooling_dictionary : dict
		A dictionary containing the pooling interpolating functions, where
		key:value pairs are pool identiiers: scipy interpolating functions
		of the pools
	"""
	# Check inputs
	if not isinstance(concentration_profile, dict):
		raise TypeError("concentration_profile must be a dictionary")
	elif len(concentration_profile) == 0:
		warn ("No concentration solutions associated with this profile")
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

	if isinstance(pools, string_types):
		pools = [pools]
	if not isinstance(pools, list):
		raise TypeError("pools must be a string or list of strings")

	if parameters is None:
		parameters = {}
	if not isinstance(parameters, dict):
		raise TypeError("parameters must be a dictionary")

	if isinstance(pool_ids, string_types):
		pool_ids = [pool_ids]
	if not isinstance(pool_ids, (list, type(None))):
		raise TypeError("pool_ids must be a string or a list of strings with "
						"the same length as pools")
	else:
		# Generate pool IDs if not provided
		if pool_ids is None:
			pool_ids = ['p%s' % str(i+1) for i in range(0, len(pools))]
		# Otherwise check length of pool IDs if provided
		elif len(pool_ids) == len(pools):
			pass
		else:
			raise ValueError("Number of pool IDs does not match the number of "
							"defined pools")
	# Create pools expressions using Sympy
	pool_dict = dict()
	symbol_dict = dict()
	metab_lookup = DictList()
	# Get metabolite objects in a DictList and local symbols in a dictionary
	for key in iterkeys(concentration_profile):
		metab_lookup.append(key)
		symbol_dict.update({key.id: sp.Symbol(key.id)})
	# Add parameters into local symbol dictionary, if any
	for key in iterkeys(parameters):
		symbol_dict.update({key: sp.Symbol(key)})
	for i, p_str in enumerate(pools):
		# Use sympy.sympify to create pool expression
		try:
			pool_expr = sp.sympify(p_str, locals=symbol_dict)
			pool_expr = pool_expr.subs(dict((sp.Symbol(k), v)
											for k, v in iteritems(parameters)))
		except:
			raise ValueError("Could not convert (%s) into a pool." % p_str)
		# Set up substitution for concentration values
		sub_dict = dict((sym, concentration_profile[
						metab_lookup.get_by_id(str(sym))])
						for sym in pool_expr.atoms(sp.Symbol))
		# Obtain pooling solution by substituting concentration values for each
		# time point
		pool_sol = list()
		for t in time_range:
			val_at_t = dict((sym, sub_dict[sym](t))
							for sym, func in iteritems(sub_dict))
			pool_sol.extend([pool_expr.subs(val_at_t)])
		# Remove the last point if needed
		if pool_sol[-1] == 0:
			pool_sol = pool_sol[:-1]
			time_vec = time_range[:-1]
		else:
			time_vec = time_range
		# Add the new pool interpolating funtion into pool_dict
		pool_dict.update({pool_ids[i]: sc.interpolate.interp1d(time_vec,
							pool_sol, kind='cubic', fill_value='extrapolate')})
	return pool_dict

def net_fluxes_from_strings(flux_profile, time_range, net_fluxes,
				parameters=None, net_flux_ids=None, numpoints=500):
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
	# Check inputs
	if not isinstance(flux_profile, dict):
		raise TypeError("flux_profile must be a dictionary")
	elif len(flux_profile) == 0:
		warn ("No flux solutions associated with this profile")
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

	if isinstance(net_fluxes, string_types):
		net_fluxes = [net_fluxes]
	if not isinstance(net_fluxes, list):
		raise TypeError("net_fluxes must be a string or list of strings")

	if parameters is None:
		parameters = {}
	if not isinstance(parameters, dict):
		raise TypeError("parameters must be a dictionary")

	if isinstance(net_flux_ids, string_types):
		net_flux_ids = [net_flux_ids]
	if not isinstance(net_flux_ids, (list, type(None))):
		raise TypeError("net_flux_ids must be a string or a list of strings with "
						"the same length as net_fluxes")
	else:
		# Generate net flux IDs if not provided
		if net_flux_ids is None:
			net_flux_ids = ['v_net%s' % str(i+1) for i in range(0, len(net_fluxes))]
		# Otherwise check length of net_flux IDs if provided
		elif len(net_flux_ids) == len(net_fluxes):
			pass
		else:
			raise ValueError("Number of net_flux IDs does not match the number"
							" of defined net_fluxes")
	# Create net flux expressions using Sympy
	net_flux_dict = dict()
	symbol_dict = dict()
	reaction_lookup = DictList()
	# Get metabolite objects in a DictList and local symbols in a dictionary
	for key in iterkeys(flux_profile):
		reaction_lookup.append(key)
		symbol_dict.update({key.id: sp.Symbol(key.id)})
	# Add parameters into local symbol dictionary, if any
	for key in iterkeys(parameters):
		symbol_dict.update({key: sp.Symbol(key)})
	for i, f_str in enumerate(net_fluxes):
		# Use sympy.sympify to create net_flux expression
		try:
			net_flux_expr = sp.sympify(f_str, locals=symbol_dict)
			net_flux_expr = net_flux_expr.subs(dict((sp.Symbol(k), v)
											for k, v in iteritems(parameters)))
		except:
			raise ValueError("Could not convert (%s) into a net_flux." % f_str)
		# Set up substitution for flux values
		sub_dict = dict((sym, flux_profile[
						reaction_lookup.get_by_id(str(sym))])
						for sym in net_flux_expr.atoms(sp.Symbol))
		# Obtain net flux solution by substituting flux values for each
		# time point
		net_flux_sol = list()
		for t in time_range:
			val_at_t = dict((sym, sub_dict[sym](t))
							for sym, func in iteritems(sub_dict))
			net_flux_sol.extend([net_flux_expr.subs(val_at_t)])
		# Remove the last point if needed
		if net_flux_sol[-1] == 0:
			net_flux_sol = net_flux_sol[:-1]
			time_vec = time_range[:-1]
		else:
			time_vec = time_range
		# Add the new net_flux interpolating funtion into net_flux_dict
		net_flux_dict.update({net_flux_ids[i]: sc.interpolate.interp1d(
								time_vec, net_flux_sol,
								kind='cubic', fill_value='extrapolate')})
	return net_flux_dict

# Class begins
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
