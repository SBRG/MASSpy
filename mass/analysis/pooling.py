# -*- coding: utf-8 -*-

# Compatibility with Python 2.7
from __future__ import absolute_import

import numpy as np
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
		correlation_cutoff : float
		A value between 0 and 1. that determines the cutoff for correlation.
		Any cos(angle) < correlation_cutoff will be set at zero.

	Returns
	-------
	p_matrices : np.array of length n
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
