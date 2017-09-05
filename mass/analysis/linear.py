# -*- coding: utf-8 -*-

# Compatibility with Python 2.7
from __future__ import absolute_import

# Import necesary packages
import numpy as np
import pandas as pd
from sympy import Matrix
from scipy.sparse import dok_matrix, lil_matrix
from six import string_types, iteritems

# from mass
from mass.core import massmodel
# Class begins

def nullspace(A, orient="col", integers=False):
	"""Compute an approximate basis for the nullspace of A.

	The sympy.Matrix.nullspace method is utilized to calculate the ns.
	The algorithm utilizes Gaussian elimination with back substitution

	Parameters
	----------
	A : numpy.ndarray, scipy.sparse dok matrix or lil matrix, pandas
		A should be at most 2-D.  A 1-D array with length k will be treated
		as a 2-D with shape (1, k)
	orient : 'row' or 'col'
		Whether to output the nullspace as row vectors or column vectors
	integers : bool
		If true, will find the least common denominator in
		each vector and turn the entries into integers

	Returns
	-------
	numpy.ndarray of shape
		If `A` is an array with shape (m, k), then `ns` will be an array
		with shape (k, n) for orientation 'col',where n is the estimated
		dimension of the nullspace of `A`. The columns of ns are a basis for
		nullspace; each element in numpy.dot(A, ns) will be approximately zero.
	"""
	if isinstance(A, np.ndarray) or isinstance(A, pd.DataFrame):
		pass
	elif isinstance(A, dok_matrix) or isinstance(A, lil_matrix):
		A = A.toarray()
	else:
		raise TypeError("Matrix must be one of the following formats: "
				"numpy.ndarray, dok_matrix, lil_matrix, or pandas 'DataFrame'")

	if orient not in ("row" "col"):
		raise TypeError("Orientation must be either 'row' or 'col'")

	A = np.atleast_2d(A)
	ns = np.array(Matrix(A).nullspace()).astype(np.float64)
	# Make integers if True
	if integers:
		for i in range(ns.shape[0]):
			# Find abs min value with index
			ma = np.ma.masked_values(np.absolute(ns[i]), 0.0, copy=False)
			# Rationalize by dividing by abs min value
			ns[i] = np.divide(ns[i], ma.min())

	# Set orientation
	if orient is "col":
		ns = ns.T

	return ns


def left_nullspace(A, orient="row", integers=False):
	"""Compute an approximate basis for the left nullspace of A.

	The sympy.Matrix.nullspace method is utilized to calculate the lns
	The algorithm utilizes Gaussian elimination with back substitution

	Parameters
	----------
	A : numpy.ndarray, scipy.sparse dok matrix or lil matrix, pandas
		A should be at most 2-D.  A 1-D array with length k will be treated
		as a 2-D with shape (1, k)
	orient : 'row' or 'col'
		Whether to output the nullspace as row vectors or column vectors
	integers : bool
		If true, will find the least common denominator for fractions in
		each vector and turn the entries into integers

	Returns
	-------
	numpy.ndarray of shape
		If `A` is an array with shape (m, k), then `lns` will be an array
		with shape (n, k) for orientation "row" where n is the estimated
		dimension of the left nullspace of `A`.
		The rows of lns are a basis for the left nullspace; each element
		in numpy.dot(lns, A) will be approximately zero.
	"""
	if isinstance(A, np.ndarray) or isinstance(A, pd.DataFrame):
		pass
	elif isinstance(A, dok_matrix) or isinstance(A, lil_matrix):
		A = A.toarray()
	else:
		raise TypeError("Matrix must be one of the following formats: "
				"numpy.ndarray, dok_matrix, lil_matrix, or pandas 'DataFrame'")
	A = np.atleast_2d(A)
	lns = nullspace(np.transpose(A), orient, integers)
	return lns

def matrix_rank(A):
	"""Get the rank of a matrix.

	Utilizes sympy.Matrix.rank method

	Parameters
	----------
	A : numpy.ndarray, dok_matrix, lil_matrix, or pandas 'DataFrame'
		A should be at most 2-D.  A 1-D array with length k will be treated
		as a 2-D with shape (1, k)
	"""
	if isinstance(A, np.ndarray) or isinstance(A, pd.DataFrame):
		pass
	elif isinstance(A, dok_matrix) or isinstance(A, lil_matrix):
		A = A.toarray()
	else:
		raise TypeError("Matrix must be one of the following formats: "
				"numpy.ndarray, dok_matrix, lil_matrix, or pandas 'DataFrame'")

	return Matrix(A).rank()
