# -*- coding: utf-8 -*-

# Compatibility with Python 2.7
from __future__ import absolute_import

# Import necesary packages
import re
import numpy as np
import pandas as pd
import sympy as sp
from scipy.sparse import dok_matrix, lil_matrix
from six import string_types, iteritems

# from mass
from mass.core import expressions
from mass.core import massmodel

# Class begins
## Public
def create_stoichiometric_matrix(model, matrix_type=None, dtype=None,
								update_model=True):
	"""Return the stoichiometrix matrix representation for a given massmodel

	The rows represent the chemical species and the columns represent the
	reactions. S[i, j] therefore contain the quantity of species 'i' produced
	(positive) or consumed (negative) by reaction 'j'.

	Matrix types can include 'dense' for a standard  numpy.array, 'dok' or
	'lil' to obtain the scipy matrix of the corresponding type, DataFrame for
	a pandas 'Dataframe' where species (excluding genes) are row indicies and
	reactions are column indicices, and 'symbolic' for a sympy.Matrix.

	Parameters
	----------
	model : mass.MassModel
		The MassModel object to construct the matrix for
	matrix_type: {'dense', 'dok', 'lil', 'dataframe', 'symbolic'}, or None
	   Construct the S matrix with the specified matrix type. If None, will
	   utilize the matrix type in the massmodel. If massmodel does not have a
	   specified matrix type, will default to 'dense'
	   Not case sensitive
	dtype : data-type
		Construct the S matrix with the specified data type. If None, will
		utilize  the data type in the massmodel. If massmodel does not have a
		specified data type, will default to float64
	Returns
	-------
	matrix of class 'dtype'
		The stoichiometric matrix for the given MassModel
	"""
	# Check input of update model
	if not isinstance(update_model, bool):
		raise TypeError("update_model must be a bool")

	# Set up for matrix construction if matrix types are correct.
	(matrix_constructor, dtype) = _setup_matrix_constructor(model,
										matrix_type=matrix_type, dtype=dtype)
	n_metabolites = len(model.metabolites)
	n_reactions = len(model.reactions)

	# No need to construct a matrix if there are no metabolites or species
	if n_metabolites == 0 or n_reactions == 0:
		s_matrix = None

	else:
		# Construct the stoichiometric matrix
		s_matrix = matrix_constructor((n_metabolites, n_reactions),
										dtype=dtype)
		# Get index for metabolites and reactions
		m_ind = model.metabolites.index
		r_ind = model.reactions.index

		# Build matrix
		for rxn in model.reactions:
			for metab, stoic in iteritems(rxn.metabolites):
				s_matrix[m_ind(metab), r_ind(rxn)] = stoic

		# Convert matrix to dataframe if matrix type is a dataframe
		if matrix_type == 'dataframe':
			metabolite_ids =[metab.id for metab in model.metabolites]
			reaction_ids = [rxn.id for rxn in model.reactions]
			s_matrix = pd.DataFrame(s_matrix, index = metabolite_ids,
											columns = reaction_ids)
		if matrix_type == 'symbolic':
			s_matrix = sp.Matrix(s_matrix)

		# Update the model's stored matrix data if True
	if update_model:
		_update_model_s(model, s_matrix, matrix_type, dtype)

	return s_matrix

def create_gradient_matrix(model, symbolic=True, matrix_type='dense'):
	"""Return the gradient matrix representation for a given massmodel

	Matrix types can include 'dense' for a standard  numpy.array and
	'symbolic' for a sympy.Matrix.

	Parameters
	----------
	model : mass.MassModel
		The MassModel object to construct the matrix for
	symbolic : bool
		If True, will output the gradient matrix with derivatives of the rates
		represented in their symbolic forms. Otherwise numerical values are
		substituted for kinetic parameters
	matrix_type: {'dense', 'symbolic'},
	   Construct the G matrix with the specified matrix type.

	Returns
	-------
	sympy.Matrix
		The gradient matrix for the given MassModel
	"""
	if not isinstance(matrix_type, string_types):
		raise TypeError("matrix_type must be a string")
	matrix_type = matrix_type.lower()
	if matrix_type not in {'dense','symbolic'}:
		raise ValueError("matrix_type must be either 'dense', 'symbolic'}")

	# Set up for matrix construction if matrix types are correct.
	n_metabolites = len(model.metabolites)
	n_reactions = len(model.reactions)

	# No need to construct a matrix if there are no metabolites or species
	if n_metabolites == 0 or n_reactions == 0:
		g_matrix = None

	else:
		# Construct the stoichiometric matrix
		g_matrix = sp.Matrix(np.zeros((n_reactions, n_metabolites)))
		# Get index for metabolites and reactions
		r_ind = model.reactions.index
		m_ind = model.metabolites.index

		# Get rates and symbols
		odes, rates, symbols = expressions._sort_symbols(model)

		# Get values for substitution if necessary
		if not symbolic:
			values = dict()
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
			# Subsitute values into rates
			for rxn, rate in iteritems(rates):
				rates[rxn] = rate.subs(values)

		# Strip time in order to take derivatives
		rates = expressions.strip_time(rates)
		# Create gradient matrix
		for rxn, rate in iteritems(rates):
			for metab in model.metabolites:
				sym = sp.Symbol(metab.id)
				g_matrix[r_ind(rxn), m_ind(metab)] = rate.diff(sym)

	return g_matrix

## Internal
def _setup_matrix_constructor(model, matrix_type=None, dtype=None):
	"""Internal use. Check inputs and create a constructor for the specified
	matrix type.

	Parameters
	----------
	model : mass.MassModel
		The MassModel object to construct the matrix for
	matrix_type: {'dense', 'dok', 'lil', 'dataframe', 'symbolic'}, or None
	   Construct the matrix with the specified matrix type. If None, will
	   utilize  the matrix type in the massmodel. If massmodel does not have a
	   specified matrix type, will default to 'dense'
	   Not case sensitive
	dtype : data-type
		Construct the S matrix with the specified data type. If None, will
		utilize  the data type in the massmodel. If massmodel does not have a
		specified data type, will default to float64
	Returns
	-------
	matrix of class 'dtype'
	"""
	# Dictionary for constructing the matrix
	matrix_constructor = {'dense': np.zeros, 'dok': dok_matrix,
				'lil': lil_matrix, 'dataframe': np.zeros, 'symbolic': np.zeros}

	if not isinstance(model, massmodel.MassModel):
		raise TypeError("model must be a MassModel")

	# Check matrix type input if it exists
	if matrix_type is not None:
		if not isinstance(matrix_type, string_types):
			raise TypeError("matrix_type must be a string")
		# Remove case sensitivity
		matrix_type = matrix_type.lower()
	else:
		# Use the models stored matrix type if None is specified
		if model._matrix_type is not None:
			matrix_type = model._matrix_type.lower()
		# Otherwise use the default type, 'dense'
		else:
			matrix_type = 'dense'
			model._matrix_type = 'dense'

	# Check to see if matrix type is one of the defined types
	if matrix_type not in matrix_constructor:
		raise ValueError("matrix_type must be a string of one of the following"
					" types: {'dense', 'dok', 'lil', 'dataframe', 'symbolic'}")

	# Set the data-type if it is none
	if dtype is None:
		# Use the models stored data type if available
		if model._dtype is not None:
			dtype = model._dtype
		# Otherwise use the default type, np.float64
		else:
			dtype = np.float64

	constructor = matrix_constructor[matrix_type]
	return (constructor, dtype)

def _update_S(model, reaction_list=None, matrix_type=None, dtype=None,
			update_model=True):
	"""For internal use only. Update the S matrix of the model.

	NOTE: reaction_list is assumed to be at the end of self.reactions.

	Warnings
	--------
	This method is intended for internal use only. To safely convert a matrix
	to another type of matrix, use the massmodel.update_S method instead

	Parameters
	----------
	model : mass.MassModel
		The MassModel object to construct the matrix for
	reaction_list : list of MassReactions or None
		The list of MassReactions to add to the current stoichiometric matrix.
		Reactions must already exist in the model in order to update.
		If None, the entire stoichiometric matrix is reconstructed
	matrix_type: string {'dense', 'dok', 'lil', 'DataFrame', 'symbolic'}
		The type of matrix
	dtype : data-type
		The desired data-type for the array. If None, defaults to float

	Returns
	-------
	matrix of class 'dtype'
		The stoichiometric matrix for the given MassModel
	"""
	if not isinstance(model, massmodel.MassModel):
		raise TypeError("model must be a MassModel")

	# Check matrix type input if it exists to ensure its a valid matrix type
	if matrix_type is not None:
		if not isinstance(matrix_type, string_types):
			raise TypeError("matrix_type must be a string")
		# Remove case sensitivity
		matrix_type = matrix_type.lower()
		if matrix_type not in {'dense', 'dok', 'lil', 'dataframe', 'symbolic'}:
			raise ValueError("matrix_type must be of one of the following"
					" types: {'dense', 'dok', 'lil', 'dataframe', 'symbolic'}")
	else:
		matrix_type = model._matrix_type

	# Use the model's stored datatype if the datatype is not specified
	if dtype is None:
		dtype = model._dtype

	# Check input of update model
	if not isinstance(update_model, bool):
		raise TypeError("update_model must be a bool")

	# If there is no change to the reactions, just reconstruct the model
	if model._S is None or reaction_list is None:
		s_matrix = create_stoichiometric_matrix(model,
						matrix_type=matrix_type,
						dtype=dtype,
						update_model=update_model)
	else:
		s_matrix = _update_stoichiometry(model, reaction_list,
										matrix_type=matrix_type)

	if update_model:
		_update_model_s(model, s_matrix, matrix_type, dtype)

	return s_matrix

def _update_stoichiometry(model, reaction_list, matrix_type=None):
	"""For internal uses only. To update the stoichometric matrix with
	additional reactions and metabolites efficiently by converting to
	a dok matrix, updating the dok matrix, and converting back to the
	desired type

	Parameters
	----------
	massmodel : mass.MassModel
		The massmodel to update
	reaction_list: list of MassReactions
		The reactions to add to the matrix
	matrix_type: string {'dense', 'dok', 'lil', 'DataFrame', 'symbolic'}
		The type of matrix

	Warnings
	--------
	This method is intended for internal use only. To safely update a matrix
	use the massmodel.update_S method instead
	"""
	# Set defaults
	shape = (len(model.metabolites), len(model.reactions))
	if matrix_type is None:
		matrix_type = 'dense'

	# Get the S matrix as a dok matrix
	s_matrix = _convert_S(model._S, 'dok')
	# Resize the matrix
	s_matrix.resize(shape)

	# Update the matrix
	coefficient_dictionary = {}
	for rxn in reaction_list:
		rxn_index = model.reactions.index(rxn.id)
		for metab, coeff in rxn._metabolites.items():
			coefficient_dictionary[(model.metabolites.index(metab.id),
									rxn_index)] = coeff
	s_matrix.update(coefficient_dictionary)

	# Convert the matrix to the desired type
	s_matrix = _convert_S(s_matrix, matrix_type)
	if matrix_type == 'dataframe':
		metabolite_ids =[metab.id for metab in model.metabolites]
		reaction_ids = [rxn.id for rxn in model.reactions]
		s_matrix = pd.DataFrame(s_matrix, index = metabolite_ids,
										columns = reaction_ids)
	if matrix_type == 'symbolic':
		s_matrix = sp.Matrix(s_matrix)

	return s_matrix

def _convert_S(s_matrix, matrix_type):
	"""For internal uses only. To convert a matrix to a different type.

	Parameters
	----------
	s_matrix : matrix of class "dtype"
		The S matrix for conversion
	matrix_type: string {'dense', 'lil', 'dok', 'DataFrame', 'symbolic'}
		The type of matrix to convert to

	Warnings
	--------
	This method is intended for internal use only. To safely convert a matrix
	to another type of matrix, use the massmodel.update_S method instead
	"""
	def _to_dense(s_mat=s_matrix):
		if isinstance(s_mat, np.ndarray):
			return s_mat
		elif isinstance(s_mat, pd.DataFrame):
			return s_mat.as_matrix()
		elif isinstance(s_mat, sp.Matrix):
			return np.array(s_mat)
		else:
			return s_mat.toarray()

	def _to_lil(s_mat=s_matrix):
		if isinstance(s_mat, sp.Matrix):
			s_mat = np.array(s_mat)
		return lil_matrix(s_mat)

	def _to_dok(s_mat=s_matrix):
		if isinstance(s_mat, sp.Matrix):
			s_mat = np.array(s_mat)
		return dok_matrix(s_mat)

	matrix_conversion = {'dense': _to_dense,
						'lil' : _to_lil,
						'dok' : _to_dok,
						'dataframe' : _to_dense,
						'symbolic' : to_dense}

	s_matrix = matrix_conversion[matrix_type](s_mat=s_matrix)
	return s_matrix

def _update_model_s(model, s_matrix, matrix_type, dtype):
	"""For internal use only. Update the model storage of the s matrix,
	matrix type, and data type

	Warnings
	--------
	This method is intended for internal use only. To safely convert a matrix
	to another type of matrix, use the massmodel.update_S method instead
	"""
	model._S = s_matrix
	model._matrix_type = matrix_type
	model._dtype = dtype
