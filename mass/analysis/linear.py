# -*- coding: utf-8 -*-

# Compatibility with Python 2.7
from __future__ import absolute_import

# Import necesary packages
import numpy as np
import pandas as pd
from sympy import Matrix
from scipy.sparse import dok_matrix, lil_matrix
from six import integer_types, string_types, iteritems

# from cobra
from cobra.util import array
# from mass
from mass.core import massmodel
# Class begins
def nullspace(A, orient="col", method='gaussian', integers=False,
				atol=1e-13, rtol=0.):
	"""Compute an approximate basis for the nullspace of A.

	If the method is set to 'gaussian', the sympy.Matrix.nullspace method
	is utilized to calculate the nullspace.The algorithm uses Gaussian
	elimination with back substitution.

	If the method is set to 'svd', the numpy.linalg.svd method is utilized to
	calculate the nullspace. The algorithm is from the numpy cookbook.

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
		Will be ignored if method is not 'gaussian'.
	method : 'gaussian', 'svd'
		Whether to obtain the nullspace through Gaussian elimination or through
		singular value decomposition (svd).
	atol : float (Only for method SVD)
        The absolute tolerance for a zero singular value.  Singular values
        smaller than `atol` are considered to be zero.
		Will be ignored if method is not 'svd'
    rtol : float (Only for method SVD)
        The relative tolerance.  Singular values less than rtol*smax are
        considered to be zero, where smax is the largest singular value.
		Will be ignored if method is not 'svd'.

    If both `atol` and `rtol` are positive, the combined tolerance is the
    maximum of the two; that is::
    tol = max(atol, rtol * smax)
    Singular values smaller than `tol` are considered to be zero.

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
				"numpy.ndarray, scipy.dok_matrix, scipy.lil_matrix, "
				"or pandas.DataFrame.")

	if orient not in ("row" "col"):
		raise TypeError("Orientation must be either 'row' or 'col'.")

	if method not in {'gaussian', 'svd'}:
		raise TypeError("method must be either 'gaussian' or 'svd'.")

	if method is 'gaussian':
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

	else:
		if not isinstance(atol, (integer_types, float)):
			raise TypeError("atol must be a float or integer type")
		if not isinstance(rtol, (integer_types, float)):
			raise TypeError("atol must be a float or integer type")
		ns = array.nullspace(A, atol, rtol)
		if orient is "row":
			ns = ns.T
	return ns

def left_nullspace(A, orient="row", method='gaussian',integers=False,
				atol=1e-13, rtol=0.):
	"""Compute an approximate basis for the left nullspace of A.

	If the method is set to 'gaussian', the sympy.Matrix.nullspace method
	is utilized to calculate the left nullspace.The algorithm uses Gaussian
	elimination with back substitution.

	If the method is set to 'svd', the numpy.linalg.svd method is utilized to
	calculate the left nullspace. The algorithm is from the numpy cookbook.

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
	method : 'gaussian', 'svd'
		Whether to obtain the nullspace through Gaussian elimination or through
		singular value decomposition (svd).
	atol : float (Only for method SVD)
        The absolute tolerance for a zero singular value.  Singular values
        smaller than `atol` are considered to be zero.
		Will be ignored if method is not 'svd'
    rtol : float (Only for method SVD)
        The relative tolerance.  Singular values less than rtol*smax are
        considered to be zero, where smax is the largest singular value.
		Will be ignored if method is not 'svd'.

    If both `atol` and `rtol` are positive, the combined tolerance is the
    maximum of the two; that is::
    tol = max(atol, rtol * smax)
    Singular values smaller than `tol` are considered to be zero.

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
	lns = nullspace(np.transpose(A), orient, method, integers, atol, rtol)

	return lns

def matrix_rank(A, method='gaussian', tol=None):
	"""Get the rank of a matrix.

	If the method is set to 'gaussian', the sympy.Matrix.rank method
	is utilized to calculate the matrix rank.The algorithm uses Gaussian
	elimination with back substitution.

	If the method is set to 'svd', the numpy.linalg.matrix_rank method is
	utilized to calculate the matrix rank. The algorithm uses SVD.

	Utilizes

	Parameters
	----------
	A : numpy.ndarray, dok_matrix, lil_matrix, or pandas 'DataFrame'
		A should be at most 2-D.  A 1-D array with length k will be treated
		as a 2-D with shape (1, k)
	tol : float or None
		Threshold below which SVD values are considered 0. If tol is None,
		and S is an array with singular values for A, and eps is the epsilon
		value for datatype of S, then tol is S.max() * max(A.shape) * eps.
		Will be ignored if method is not 'svd'.
	"""
	if isinstance(A, np.ndarray) or isinstance(A, pd.DataFrame):
		pass
	elif isinstance(A, dok_matrix) or isinstance(A, lil_matrix):
		A = A.toarray()
	else:
		raise TypeError("Matrix must be one of the following formats: "
				"numpy.ndarray, dok_matrix, lil_matrix, or pandas 'DataFrame'")

	if method not in {'gaussian', 'svd'}:
		raise TypeError("method must be either 'gaussian' or 'svd'.")

	if method is 'gaussian':
		return Matrix(A).rank()
	else:
		return np.linalg.matrix_rank(A, tol)
