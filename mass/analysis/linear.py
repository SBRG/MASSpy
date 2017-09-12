# -*- coding: utf-8 -*-

# Compatibility with Python 2.7
from __future__ import absolute_import

# Import necesary packages
import re
import scipy as sc
import numpy as np
import sympy as sp
import pandas as pd
from scipy.sparse import dok_matrix, lil_matrix
from six import integer_types, string_types, iteritems, iterkeys

# from mass
from mass.core import expressions
from mass.core import massmodel

# Class begins
## Global symbol for time
t = sp.Symbol("t")
def gradient(model, strip_time=False, sub_parameters=False,
				sub_concentrations=False, matrix_type='symbolic'):
	"""Create the gradient matrix for a given massmodel

	Matrix types can include 'dense' for a standard  numpy.array,
	'DataFrame' for a pandas.DataFrame, and 'symbolic' for a sympy.Matrix.

	Parameters
	----------
	model : mass.MassModel
		The MassModel object to construct the matrix for
	strip_time : bool
		If True, turn metabolite functions into symbols by stripping the time
		dependency.(Metab(t) -> Metab)
	subs_parameters : bool
		If True, will output the gradient matrix with derivatives of the rates
		represented with numerical values substituted for kinetic parameters.
		Otherwise output derivatives in as symbolic expressions without
		substituting numerical values.
	sub_concentrations :  bool
		If True, will output the gradient matrix with derivatives of the rates
		represented with the numerical values of the initial conditions
		substituted for metabolite concentrations. Otherwise output derivatives
		in as symbolic expressions without substituting numerical values.
	matrix_type: {'dense', 'dataframe', 'symbolic'},
		Matrix return type.

	Returns
	-------
	numpy.array, pandas.DataFrame or sympy.Matrix
		The gradient matrix for the given MassModel
	"""
	if not isinstance(matrix_type, string_types):
		raise TypeError("matrix_type must be a string")
	matrix_type = matrix_type.lower()
	if matrix_type not in {'dense', 'dataframe','symbolic'}:
		raise ValueError("matrix_type must be one of the following types: "
						"'dense', 'DataFrame', or 'symbolic'}")

	# Set up for matrix construction if matrix types are correct.
	n_metabolites = len(model.metabolites)
	n_reactions = len(model.reactions)

	# No need to construct a matrix if there are no metabolites or species
	if n_metabolites == 0 or n_reactions == 0:
		return None

	# Construct the gradient matrix
	g_matrix = sp.Matrix(np.zeros((n_reactions, n_metabolites)))
	# Get index for metabolites and reactions
	r_ind = model.reactions.index
	m_ind = model.metabolites.index

	# Get rates and symbols
	odes, rates, symbols = expressions._sort_symbols(model)

	# Strip time if desired
	if strip_time:
		rates = expressions.strip_time(rates)
		metab_list = [sp.Symbol(m.id) for m in model.metabolites]
	else:
		metab_list = [sp.Symbol(m.id)(t) for m in model.metabolites]

	# Get values for substitution if necessary
	values = dict()
	if sub_concentrations:
		values.update(dict((metab_list[i], model.initial_conditions[m])
					for i, m in enumerate(model.metabolites)
					if m in iterkeys(model.initial_conditions)))

	if sub_parameters:
		# For rate parameters
		for param_sym in symbols[1]:
			[p_type, rid] = re.split("\_", str(param_sym), maxsplit=1)
			reaction = model.reactions.get_by_id(rid)
			prop_f = reaction.__class__.__dict__[p_type]
			values.update({param_sym: prop_f.fget(reaction)})
		# For fixed concentrations
		for metab_sym in symbols[2]:
			metab = str(metab_sym)
			if re.search("\_Xt",metab):
				values.update({metab_sym :
								model.fixed_concentrations[metab]})
			else:
				metab = model.metabolites.get_by_id(metab)
				values.update({metab_sym :
								model.fixed_concentrations[metab]})
		# For custom_parameters
		for c_param in symbols[3]:
			values.update({c_param: model.custom_parameters[str(c_param)]})

	# Create gradient matrix
	for rxn, rate in iteritems(rates):
		for i, metab in enumerate(model.metabolites):
			g_matrix[r_ind(rxn), m_ind(metab)] = rate.diff(metab_list[i])

	# Substitution if necessary
	if sub_parameters or sub_concentrations:
		g_matrix = g_matrix.subs(values)

	# Convert to a numpy.array or pandas.DataFrame if specified
	if matrix_type == 'dense':
		g_matrix = np.array(g_matrix)
		# Attempt to change the datatype from object to float64 if values
		# were substituted, otherwise leave as is.
		if sub_parameters and sub_concentrations:
			try:
				g_matrix = g_matrix.astype(np.float64)
			except:
				pass

	if matrix_type == 'dataframe':
		reaction_ids = [rxn.id for rxn in model.reactions]
		metabolite_ids =[metab.id for metab in model.metabolites]
		g_matrix = pd.DataFrame(np.array(g_matrix), index = reaction_ids,
												columns=metabolite_ids)

	return g_matrix

def kappa(model, strip_time=False, sub_parameters=False,
				sub_concentrations=False, matrix_type='symbolic'):
	"""Create the kappa matrix for a given massmodel. The kappa matrix is the
	the diagnolization of the norms for the row vectors in the gradient matrix.

	Matrix types can include 'dense' for a standard  numpy.array,
	'DataFrame' for a pandas.DataFrame, and 'symbolic' for a sympy.Matrix.

	Parameters
	----------
	model : mass.MassModel
		The MassModel object to construct the matrix for
	strip_time : bool
		If True, turn metabolite functions into symbols by stripping the time
		dependency.(Metab(t) -> Metab)
	subs_parameters : bool
		If True, will output the gradient matrix with derivatives of the rates
		represented with numerical values substituted for kinetic parameters.
		Otherwise output derivatives in as symbolic expressions without
		substituting numerical values.
	sub_concentrations :  bool
		If True, will output the gradient matrix with derivatives of the rates
		represented with the numerical values of the initial conditions
		substituted for metabolite concentrations. Otherwise output derivatives
		in as symbolic expressions without substituting numerical values.
	matrix_type: {'dense', 'dataframe', 'symbolic'},
		Matrix return type.

	Returns
	-------
	numpy.array, pandas.DataFrame or sympy.Matrix
		The kappa matrix for the given MassModel
	"""
	g_matrix = gradient(model, strip_time, sub_parameters,
						sub_concentrations, matrix_type='symbolic')
	kappa = sp.diag(*[g_matrix[r,:].norm()
					for r in range(g_matrix.rows)]).subs({sp.nan: sp.S.Zero})

	if matrix_type == 'dense':
		kappa = np.array(kappa)
		# Attempt to change the datatype from object to float64 if values
		# were substituted, otherwise leave as is.
		if sub_parameters and sub_concentrations:
			try:
				kappa = kappa.astype(np.float64)
			except:
				pass

	if matrix_type == 'dataframe':
		reaction_ids = [rxn.id for rxn in model.reactions]
		kappa = pd.DataFrame(np.array(kappa), index=reaction_ids,
											columns=reaction_ids)
	return kappa

def gamma(model, strip_time=False, sub_parameters=False,
				sub_concentrations=False, matrix_type='symbolic'):
	"""Create the gamma matrix for a given massmodel. The gamma matrix is
	the 1-norms of the gradient matrix.

	Matrix types can include 'dense' for a standard  numpy.array,
	'DataFrame' for a pandas.DataFrame, and 'symbolic' for a sympy.Matrix.

	Parameters
	----------
	model : mass.MassModel
		The MassModel object to construct the matrix for
	strip_time : bool
		If True, turn metabolite functions into symbols by stripping the time
		dependency.(Metab(t) -> Metab)
	subs_parameters : bool
		If True, will output the gradient matrix with derivatives of the rates
		represented with numerical values substituted for kinetic parameters.
		Otherwise output derivatives in as symbolic expressions without
		substituting numerical values.
	sub_concentrations :  bool
		If True, will output the gradient matrix with derivatives of the rates
		represented with the numerical values of the initial conditions
		substituted for metabolite concentrations. Otherwise output derivatives
		in as symbolic expressions without substituting numerical values.
	matrix_type: {'dense', 'dataframe', 'symbolic'},
		Matrix return type.

	Returns
	-------
	numpy.array, pandas.DataFrame or sympy.Matrix
		The gamma matrix for the given MassModel
	"""
	g_matrix = gradient(model, strip_time, sub_parameters,
						sub_concentrations, matrix_type='symbolic')
	gamma = sp.Matrix([g_matrix[r, :].normalized()
				for r in range(g_matrix.rows)]).subs({sp.nan: sp.S.Zero})


	if matrix_type == 'dense':
		gamma = np.array(gamma)
		# Attempt to change the datatype from object to float64 if values
		# were substituted, otherwise leave as is.
		if sub_parameters and sub_concentrations:
			try:
				gamma = gamma.astype(np.float64)
			except:
				pass

	if matrix_type == 'dataframe':
		metabolite_ids =[metab.id for metab in model.metabolites]
		reaction_ids = [rxn.id for rxn in model.reactions]
		gamma = pd.DataFrame(np.array(gamma), index=reaction_ids,
											columns=metabolite_ids)
	return gamma

def jacobian(model, jacobian_type="metabolite", strip_time=False,
			sub_parameters=False, sub_concentrations=False,
			matrix_type='symbolic'):
	"""Get the jacobian matrix for a given massmodel

	Matrix types can include 'dense' for a standard  numpy.array,
	'DataFrame' for a pandas.DataFrame, and 'symbolic' for a sympy.Matrix.

	Parameters
	----------
	model : mass.MassModel
		The MassModel object to construct the matrix for
	jacobian_type : {'metabolite', 'reaction'}
		Whether to obtain the jacobian with respect to the metabolites (Jx)
		or to obbtain the jacobian with respect to the reactions (Jv)
	strip_time : bool
		If True, turn metabolite functions into symbols by stripping the time
		dependency.(Metab(t) -> Metab)
	subs_parameters : bool
		If True, will output the gradient matrix with derivatives of the rates
		represented with numerical values substituted for kinetic parameters.
		Otherwise output derivatives in as symbolic expressions without
		substituting numerical values.
	sub_concentrations :  bool
		If True, will output the gradient matrix with derivatives of the rates
		represented with the numerical values of the initial conditions
		substituted for metabolite concentrations. Otherwise output derivatives
		in as symbolic expressions without substituting numerical values.
	matrix_type: {'dense', 'symbolic'},
	   Construct the J matrix with the specified matrix type.

	Returns
	-------
	numpy.array or sympy.Matrix
		The gradient matrix for the given MassModel
	"""
	if jacobian_type not in {'metabolite', 'reaction'}:
		raise ValueError("jacobian_type must be either 'metabolite' "
							"or 'reaction'}")
	if not isinstance(matrix_type, string_types):
		raise TypeError("matrix_type must be a string")
	matrix_type = matrix_type.lower()
	if matrix_type not in {'dense', 'dataframe','symbolic'}:
		raise ValueError("matrix_type must be one of the following types: "
						"'dense', 'DataFrame', or 'symbolic'}")

	# Set up for matrix construction if matrix types are correct.
	n_metabolites = len(model.metabolites)
	n_reactions = len(model.reactions)

	# No need to construct a matrix if there are no metabolites or species
	if n_metabolites == 0 or n_reactions == 0:
		return None

	# Get the stoichiometric and gradient matrices
	g_matrix = gradient(model, strip_time, sub_parameters,
						sub_concentrations, matrix_type='symbolic')
	s_matrix = model._create_stoichiometric_matrix(matrix_type='symbolic',
													update_model=False)
	# Multiply matrices to get the jacobian
	if jacobian_type is "metabolite":
		j_matrix = s_matrix * g_matrix
	else:
		j_matrix = g_matrix * s_matrix
	# Convert to a numpy.array if specified
	if matrix_type == 'dense':
		j_matrix = np.array(j_matrix)
		# Attempt to change the datatype from object to float64 if values
		# were substituted, otherwise leave as is.
		if sub_parameters and sub_concentrations:
			try:
				j_matrix = j_matrix.astype(np.float64)
			except:
				pass

	if matrix_type == 'dataframe':
		if jacobian_type is "metabolite":
			metabolite_ids = [metab.id for metab in model.metabolites]
			j_matrix = pd.DataFrame(np.array(j_matrix), index=metabolite_ids,
														columns=metabolite_ids)
		else:
			reaction_ids = [rxn.id for rxn in model.reactions]
			j_matrix = pd.DataFrame(np.array(j_matrix), index=reaction_ids,
														columns=reaction_ids)
	return j_matrix

def nullspace(A, integers=False):
	"""Compute an approximate basis for the nullspace of A.

	The sympy.Matrix.nullspace method is used to calculate the nullspace.
	The algorithm uses Gaussian elimination with back substitution.

	Parameters
	----------
	A : numpy.ndarray, scipy.sparse dok matrix or lil matrix, pandas.DataFrame
		or sympy.Matrix

		Note: 'A' should be at most 2-D.  A 1-D array with length k will be
		treated as a 2-D with shape (1, k)
	integers : bool
		If true, will find the least common denominator in
		each vector and turn the entries into integers
		Will be ignored if method is not 'gaussian'.

	Returns
	-------
	numpy.ndarray
		If `A` is an array with shape (m, k), then `ns` will be an array
		with shape (k, n), where n is the estimated	dimension of the nullspace
		of `A`. The columns of ns are a basis for nullspace; each element
		in numpy.dot(A, ns) will be approximately zero.

	"""
	if isinstance(A, np.ndarray) or isinstance(A, pd.DataFrame):
		A = np.atleast_2d(A)
		A = sp.Matrix(A)
	elif isinstance(A, dok_matrix) or isinstance(A, lil_matrix):
		A = A.toarray()
		A = np.atleast_2d(A)
		A = sp.Matrix(A)
	elif isinstance(A, sp.Matrix):
		pass
	else:
		raise TypeError("Matrix must be one of the following formats: "
				"numpy.ndarray, scipy.dok_matrix, scipy.lil_matrix, "
				"pandas.DataFrame, or sympy.Matrix.")

	ns = np.array(A.nullspace()).astype(np.float64)
	# Make integers if True
	if integers:
		for i in range(ns.shape[0]):
			# Find abs min value with index
			ma = np.ma.masked_values(np.absolute(ns[i]), 0.0, copy=False)
			# Rationalize by dividing by abs min value
			ns[i] = np.divide(ns[i], ma.min())

	ns = ns.T
	return ns

def left_nullspace(A, integers=False):
	"""Compute an approximate basis for the left nullspace of A.

	The sympy.Matrix.nullspace method is used to calculate the left nullspace.
	The algorithm uses Gaussian	elimination with back substitution.

	Parameters
	----------
	A : numpy.ndarray, scipy.sparse dok matrix or lil matrix, pandas.DataFrame
		or sympy.Matrix

		Note: 'A' should be at most 2-D.  A 1-D array with length k will be
		treated as a 2-D with shape (1, k)
	integers : bool
		If true, will find the least common denominator in
		each vector and turn the entries into integers

	Returns
	-------
	numpy.ndarray
		If `A` is an array with shape (m, k), then `lns` will be an array
		with shape (m, n), where n is the estimated dimension of the left
		nullspace of `A`. The columns of lns are a basis for the left
		nullspace; each element in numpy.dot(lns.T, A) will be approximately
		zero.
	"""
	if isinstance(A, np.ndarray) or isinstance(A, pd.DataFrame):
		A = np.atleast_2d(A)
		A = sp.Matrix(A)
	elif isinstance(A, dok_matrix) or isinstance(A, lil_matrix):
		A = A.toarray()
		A = np.atleast_2d(A)
		A = sp.Matrix(A)
	elif isinstance(A, sp.Matrix):
		pass
	else:
		raise TypeError("Matrix must be one of the following formats: "
				"numpy.ndarray, scipy.dok_matrix, scipy.lil_matrix, "
				"pandas.DataFrame, or sympy.Matrix.")
	lns = nullspace(A.transpose(), integers)
	return lns

def rowspace(A):
	"""Compute an approximate basis for the rowspace of A.

	The sympy.Matrix.rowspace method is used to calculate the rowspace.
	The algorithm uses Gaussian elimination with back substitution.

	Parameters
	----------
	A : numpy.ndarray, scipy.sparse dok matrix or lil matrix, pandas.DataFrame
		or sympy.Matrix

		Note: 'A' should be at most 2-D.  A 1-D array with length k will be
		treated as a 2-D with shape (1, k)

	Returns
	-------
	numpy.ndarray
		If `A` is an array with shape (m, k), then `rs` will be an array
		with shape (k, n), where n is the estimated dimension of the left
		nullspace of `A`. The columns of rs are a basis for the rowspace.
	"""
	if isinstance(A, np.ndarray) or isinstance(A, pd.DataFrame):
		A = np.atleast_2d(A)
		A = sp.Matrix(A)
	elif isinstance(A, dok_matrix) or isinstance(A, lil_matrix):
		A = A.toarray()
		A = np.atleast_2d(A)
		A = sp.Matrix(A)
	elif isinstance(A, sp.Matrix):
		pass
	else:
		raise TypeError("Matrix must be one of the following formats: "
				"numpy.ndarray, scipy.dok_matrix, scipy.lil_matrix, "
				"pandas.DataFrame, or sympy.Matrix.")

	rs = np.array(A.rowspace()).astype(np.float64)
	rs = rs.T
	return rs

def columnspace(A):
	"""Compute an approximate basis for the columnspace of A.

	The sympy.Matrix.columnspace method is used to calculate the columnspace.
	The algorithm uses Gaussian elimination with back substitution.

	Parameters
	----------
	A : numpy.ndarray, scipy.sparse dok matrix or lil matrix, pandas.DataFrame
		or sympy.Matrix

		Note: 'A' should be at most 2-D.  A 1-D array with length k will be
		treated as a 2-D with shape (1, k)

	Returns
	-------
	numpy.ndarray
		If `A` is an array with shape (m, k), then `cs` will be an array
		with shape (m, n), where n is the estimated dimension of the left
		nullspace of `A`. The columns of cs are a basis for the rowspace.
	"""
	if isinstance(A, np.ndarray) or isinstance(A, pd.DataFrame):
		A = np.atleast_2d(A)
		A = sp.Matrix(A)
	elif isinstance(A, dok_matrix) or isinstance(A, lil_matrix):
		A = A.toarray()
		A = np.atleast_2d(A)
		A = sp.Matrix(A)
	elif isinstance(A, sp.Matrix):
		pass
	else:
		raise TypeError("Matrix must be one of the following formats: "
				"numpy.ndarray, scipy.dok_matrix, scipy.lil_matrix, "
				"pandas.DataFrame, or sympy.Matrix.")

	cs = np.array(A.columnspace()).astype(np.float64)
	cs = cs.T
	return cs

def matrix_rank(A):
	"""Get the rank of a matrix.

	The sympy.Matrix.rank method is used to calculate the matrix rank.
	The algorithm uses Gaussian elimination with back substitution.

	Parameters
	----------
	A : numpy.ndarray, scipy.sparse dok matrix or lil matrix, pandas.DataFrame
		or sympy.Matrix

		Note: 'A' should be at most 2-D.  A 1-D array with length k will be
		treated as a 2-D with shape (1, k)
	"""
	if isinstance(A, np.ndarray) or isinstance(A, pd.DataFrame):
		A = np.atleast_2d(A)
		A = sp.Matrix(A)
	elif isinstance(A, dok_matrix) or isinstance(A, lil_matrix):
		A = A.toarray()
		A = np.atleast_2d(A)
		A = sp.Matrix(A)
	elif isinstance(A, sp.Matrix):
		pass
	else:
		raise TypeError("Matrix must be one of the following formats: "
				"numpy.ndarray, scipy.dok_matrix, scipy.lil_matrix, "
				"pandas.DataFrame, or sympy.Matrix.")

	return A.rank()

def svd(A, **kwargs):
	"""Get the singular value decomposition of 'A'

	This function utilizes the scipy.linalg.svd method to obtain the singular
	value decompostion (svd) of matrix 'A'.

	Other kwargs are the same as those for scipy.linalg.svd. For more details:
	https://docs.scipy.org/doc/scipy/reference/generated/scipy.linalg.svd.html

	Parameters
	----------
	A : numpy.ndarray, scipy.sparse dok matrix or lil matrix, pandas.DataFrame
		or sympy.Matrix
	Returns
	-------
	matrix of the same type as 'A'

	See Also
	--------
	scipy.linalg.svd
		svd and its arguments are the same as this method. The only difference
		is that matrices of various formats are converted in order to ensure
		the correct input for scipy.linalg.svd.
	"""
	if isinstance(A, np.ndarray) or isinstance(A, pd.DataFrame):
		pass
	elif isinstance(A, dok_matrix) or isinstance(A, lil_matrix):
		A = A.toarray()
	elif isinstance(A, sp.Matrix):
		try:
			A = np.array(A).astype(np.float64)
		except:
			raise ValueError("Cannot have sympy symbols in the matrix. Try "
							"substituting numerical values in first")
	else:
		raise TypeError("Matrix must be one of the following formats: "
				"numpy.ndarray, scipy.dok_matrix, scipy.lil_matrix, "
				"pandas.DataFrame, or sympy.Matrix.")

	return sc.linalg.svd(A, **kwargs)

def eigenvalues(A, left_eigenvec=False, right_eigenvec=False, **kwargs):
	"""Get the eigenvalues of 'A'

	This function utilizes the scipy.linalg.eig method to obtain the
	eigenvalues and of matrix 'A'.

	Other kwargs are the same as those for scipy.linalg.svd. For more details:
	https://docs.scipy.org/doc/scipy/reference/generated/scipy.linalg.eig.html

	Parameters
	----------
	A : numpy.ndarray, scipy.sparse dok matrix or lil matrix, pandas.DataFrame
		or sympy.Matrix
	left_eigenvec : bool
		Whether to calculate and return left eigenvectors. Default is False.
	right_eigenvec : bool
		Whether to calculate and return right eigenvectors. Default is True.

	Returns
	-------
	w : (M,) or (2, M) double or complex ndarray
		The eigenvalues, each repeated according to its multiplicity.
		The shape is (M,) unless homogeneous_eigvals=True.
	vl : (M, M) double or complex ndarray
		The normalized left eigenvector corresponding to the eigenvalue w[i]
		is the column vl[:,i]. Only returned if left=True.
	vr : (M, M) double or complex ndarray
		The normalized right eigenvector corresponding to the eigenvalue w[i]
		is the column vr[:,i]. Only returned if right=True.

	See Also
	--------
	scipy.linalg.eig
		svd and its arguments are the same as this method. The only difference
		is that matrices of various formats are converted in order to ensure
		the correct input for scipy.linalg.eig.
	"""
	if isinstance(A, np.ndarray) or isinstance(A, pd.DataFrame):
		pass
	elif isinstance(A, dok_matrix) or isinstance(A, lil_matrix):
		A = A.toarray()
	elif isinstance(A, sp.Matrix):
		try:
			A = np.array(A).astype(np.float64)
		except:
			raise ValueError("Cannot have sympy symbols in the matrix. Try "
							"substituting numerical values in first")
	else:
		raise TypeError("Matrix must be one of the following formats: "
				"numpy.ndarray, scipy.dok_matrix, scipy.lil_matrix, "
				"pandas.DataFrame, or sympy.Matrix with no symbolics.")

	return sc.linalg.eig(A, left=left_eigenvec, right=right_eigenvec, **kwargs)
