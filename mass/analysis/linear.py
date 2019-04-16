# -*- coding: utf-8 -*-
"""TODO Module Docstrings."""
from __future__ import absolute_import

import numpy as np

import pandas as pd

from scipy import linalg
from scipy.sparse import dok_matrix, lil_matrix

from six import iteritems

import sympy as sym

from mass.util.expressions import _mk_met_func
from mass.util.util import convert_matrix


# Public
def gradient(model, use_parameter_values=True, use_concentration_values=True,
             matrix_type="dense"):
    """Create the gradient matrix for a given model.

    Parameters
    ----------
    model: mass.MassModel
        The MassModel object to construct the matrix for.
    use_parameter_values: bool, optional
        Whether to substitute the numerical values for parameters into the
        matrix. If True, then numerical values of the kinetic parameters are
        substituted into the matrix. Otherwise parametrs in the matrix are
        left as symbols.
    use_concentration_values: bool, optional
        Whether to substitute the numerical values for concentrations into the
        matrix. If True, then numerical values of the initial conditions are
        substituted into the matrix. Otherwise species concentrations in the
        matrix are left as symbols.
    matrix_type: {'dense', 'dok', 'lil', 'DataFrame', 'symbolic'}, optional
        A string identifiying the desired format for the returned matrix.
        Types can include 'dense' for a standard numpy.array, 'dok' or 'lil' to
        obtain the corresponding scipy.sparse matrix, 'DataFrame' for a
        pandas.DataFrame, and 'symbolic' for a sympy.MutableDenseMatrix. For
        all matrix types, species (excluding  genes) are the row indicies and
        reactions are the column indicies. Defaults to "dense".

    Returns
    -------
    gradient_mat: matrix of type 'matrix_type'
        The gradient matrix for the model returned as the given matrix_type.

    """
    # Get model rates and metabolites dependent on time.
    rates = model.rates

    # Construct the base gradient matrix
    gradient_mat = sym.Matrix(np.zeros((len(rates), len(model.metabolites))))

    # Get index for metabolites and reactions
    r_ind = model.reactions.index
    m_ind = model.metabolites.index

    # Create the gradient matrix
    for rxn, rate in iteritems(rates):
        for met in model.metabolites:
            gradient_mat[r_ind(rxn), m_ind(met)] = rate.diff(_mk_met_func(met))

    # Get values for substitution
    if use_concentration_values or use_parameter_values:
        values = {}
        if use_parameter_values:
            values.update(model._get_all_parameters())

        if use_concentration_values:
            values.update({_mk_met_func(k): v
                           for k, v in iteritems(model.initial_conditions)})

        # Substitute values into the matrix
        gradient_mat = gradient_mat.subs(values)

    # Try to cast matrix as a float if all values could be computed
    if gradient_mat.is_symbolic():
        dtype = object
    else:
        dtype = np.float64

    gradient_mat = convert_matrix(gradient_mat, matrix_type=matrix_type,
                                  dtype=dtype,
                                  row_ids=[r.id for r in model.reactions],
                                  col_ids=[m.id for m in model.metabolites])
    return gradient_mat


def kappa(model, use_parameter_values=True, use_concentration_values=True,
          matrix_type="dense"):
    """Create the kappa matrix for a given model.

    The kappa matrix is the diagnolization of the norms for the rows in the
    gradient matrix.

    Parameters
    ----------
    model: mass.MassModel
        The MassModel object to construct the matrix for.
    use_parameter_values: bool, optional
        Whether to substitute the numerical values for parameters into the
        matrix. If True, then numerical values of the kinetic parameters are
        substituted into the matrix. Otherwise parametrs in the matrix are
        left as symbols.
    use_concentration_values: bool, optional
        Whether to substitute the numerical values for concentrations into the
        matrix. If True, then numerical values of the initial conditions are
        substituted into the matrix. Otherwise species concentrations in the
        matrix are left as symbols.
    matrix_type: {'dense', 'dok', 'lil', 'DataFrame', 'symbolic'}, optional
        A string identifiying the desired format for the returned matrix.
        Types can include 'dense' for a standard numpy.array, 'dok' or 'lil' to
        obtain the corresponding scipy.sparse matrix, 'DataFrame' for a
        pandas.DataFrame, and 'symbolic' for a sympy.MutableDenseMatrix. For
        all matrix types, species (excluding  genes) are the row indicies and
        reactions are the column indicies. If None, defaults to "dense".

    Returns
    -------
    kappa_matrix: matrix of type 'matrix_type'
        The kappa matrix for the model returned as the given matrix_type.

    """
    gradient_matrix = gradient(model, use_parameter_values,
                               use_concentration_values, matrix_type)
    kappa_matrix = sym.diag(*[gradient_matrix[row, :].norm()
                              for row in range(gradient_matrix.rows)])
    kappa_matrix = kappa_matrix.subs({sym.nan: sym.S.Zero})
    kappa_matrix = convert_matrix(kappa_matrix, matrix_type=matrix_type,
                                  dtype=np.float64,
                                  row_ids=[r.id for r in model.reactions],
                                  col_ids=[r.id for r in model.reactions])
    return kappa_matrix


def gamma(model, use_parameter_values=True, use_concentration_values=True,
          matrix_type="dense"):
    """Create the gamma matrix for a given model.

    The gamma matrix are the 1-norms of the gradient matrix.

    Parameters
    ----------
    model: mass.MassModel
        The MassModel object to construct the matrix for.
    strip_time: bool, optional
        If True, will strip the time dependency on concentration solutions in
        the returned matrix. (e.g. MetabID(t) -> MetabID)
    use_parameter_values: bool, optional
        Whether to substitute the numerical values for parameters into the
        matrix. If True, then numerical values of the kinetic parameters are
        substituted into the matrix. Otherwise parametrs in the matrix are
        left as symbols.
    use_concentration_values: bool, optional
        Whether to substitute the numerical values for concentrations into the
        matrix. If True, then numerical values of the initial conditions are
        substituted into the matrix. Otherwise species concentrations in the
        matrix are left as symbols.
    matrix_type: {'dense', 'dok', 'lil', 'DataFrame', 'symbolic'}, optional
        A string identifiying the desired format for the returned matrix.
        Types can include 'dense' for a standard numpy.array, 'dok' or 'lil' to
        obtain the corresponding scipy.sparse matrix, 'DataFrame' for a
        pandas.DataFrame, and 'symbolic' for a sympy.MutableDenseMatrix. For
        all matrix types, species (excluding  genes) are the row indicies and
        reactions are the column indicies. If None, defaults to "dense".

    Returns
    -------
    gamma_matrix: matrix of type 'matrix_type'
        The gamma matrix for the model returned as the given matrix_type.

    """
    gradient_matrix = gradient(model, use_parameter_values,
                               use_concentration_values, matrix_type)
    gamma_matrix = sym.Matrix([gradient_matrix[row, :].normalized()
                               for row in range(gradient_matrix.rows)])
    gamma_matrix = gamma_matrix.subs({sym.nan: sym.S.Zero})
    gamma_matrix = convert_matrix(gamma_matrix, matrix_type=matrix_type,
                                  dtype=np.float64,
                                  row_ids=[r.id for r in model.reactions],
                                  col_ids=[m.id for m in model.metabolites])
    return gamma_matrix


def jacobian(model, jacobian_type="Jx", use_parameter_values=True,
             use_concentration_values=True, matrix_type="dense"):
    """Get the jacobian matrix for a given model.

    Parameters
    ----------
    model: mass.MassModel
        The MassModel object to construct the matrix for.
    jacobian_type: {'Jx', 'Jv'}
        Whether to obtain the jacobian matrix with respect to metabolites (Jx),
        or to obtain the jacobian matrix with respect to the reactions (Jv).
        Default is 'Jx'.
    strip_time: bool, optional
        If True, will strip the time dependency on concentration solutions in
        the returned matrix. (e.g. MetabID(t) -> MetabID)
    use_parameter_values: bool, optional
        Whether to substitute the numerical values for parameters into the
        matrix. If True, then numerical values of the kinetic parameters are
        substituted into the matrix. Otherwise parametrs in the matrix are
        left as symbols.
    use_concentration_values: bool, optional
        Whether to substitute the numerical values for concentrations into the
        matrix. If True, then numerical values of the initial conditions are
        substituted into the matrix. Otherwise species concentrations in the
        matrix are left as symbols.
    matrix_type: {'dense', 'dok', 'lil', 'DataFrame', 'symbolic'}, optional
        A string identifiying the desired format for the returned matrix.
        Types can include 'dense' for a standard numpy.array, 'dok' or 'lil' to
        obtain the corresponding scipy.sparse matrix, 'DataFrame' for a
        pandas.DataFrame, and 'symbolic' for a sympy.MutableDenseMatrix. For
        all matrix types, species (excluding  genes) are the row indicies and
        reactions are the column indicies. If None, defaults to "dense".

    Returns
    -------
    jacobian_matrix: matrix of type 'matrix_type'
        The jacobian matrix for the model returned as the given matrix_type.

    """
    if jacobian_type not in {"Jx", "Jv"}:
        raise ValueError("jacobian_type must be either 'Jx' or Jv'")

    gradient_matrix = gradient(model, use_parameter_values,
                               use_concentration_values,
                               matrix_type="symbolic")
    stoich_matrix = model._mk_stoich_matrix(matrix_type="symbolic",
                                            update_model=False)
    if jacobian_type == "Jx":
        jacobian_matrix = stoich_matrix * gradient_matrix
        identifiers = [m.id for m in model.metabolites]
    else:
        jacobian_matrix = gradient_matrix * stoich_matrix
        identifiers = [r.id for r in model.reactions]

    jacobian_matrix = convert_matrix(jacobian_matrix, matrix_type=matrix_type,
                                     dtype=np.float64,
                                     row_ids=identifiers,
                                     col_ids=identifiers)
    return jacobian_matrix


def nullspace(matrix, atol=1e-13, rtol=0):
    """Compute an approximate basis for the nullspace of `matrix`.

    The algorithm used by this function is based on the singular value
    decomposition of `matrix`.

    Parameters
    ----------
    matrix: numpy.ndarray, scipy.sparse dok matrix or lil matrix,
        pandas.DataFrame, or sympy.Matrix
        The matrix to decompose. Note: `matrix` should be at most 2-D.  A 1-D
        array with length k will be treated as a 2-D with shape (1, k).
    atol: float, optional
        The absolute tolerance for a zero singular value.  Singular values
        smaller than `atol` are considered to be zero.
    rtol: float, optional
        The relative tolerance.  Singular values less than rtol*smax are
        considered to be zero, where smax is the largest singular value.

    Returns
    -------
    ns: ndarray
        If `matrix` is an array with shape (m, k), then `ns` will be an array
        with shape (k, n), where n is the estimated dimension of the
        nullspace of `matrix`.  The columns of `ns` are a basis for the
        nullspace; each element in numpy.dot(A, ns) will be approximately
        zero.

    Notes
    -----
    If both `atol` and `rtol` are positive, the combined tolerance is the
    maximum of the two; that is::
        tol = max(atol, rtol * smax)
    Singular values smaller than `tol` are considered to be zero.

    """
    matrix = np.atleast_2d(_ensure_dense_matrix(matrix))
    s, vh = linalg.svd(matrix)[1:]
    tol = max(atol, rtol * s[0])
    nnz = (s >= tol).sum()
    ns = vh[nnz:].conj().T

    # Apply zero singular value tolerance
    for i, row in enumerate(ns):
        for j, val in enumerate(row):
            if abs(val) <= tol:
                ns[i, j] = 0.
    return ns


def left_nullspace(matrix, atol=1e-13, rtol=0):
    """Compute an approximate basis for the left nullspace of `matrix`.

    The algorithm used by this function is based on the singular value
    decomposition of `matrix`.

    Parameters
    ----------
    matrix: numpy.ndarray, scipy.sparse dok matrix or lil matrix,
        pandas.DataFrame, or sympy.Matrix
        The matrix to decompose. Note: `matrix` should be at most 2-D.  A 1-D
        array with length k will be treated as a 2-D with shape (1, k).
    atol: float, optional
        The absolute tolerance for a zero singular value.  Singular values
        smaller than `atol` are considered to be zero.
    rtol: float, optional
        The relative tolerance.  Singular values less than rtol*smax are
        considered to be zero, where smax is the largest singular value.

    Returns
    -------
    lns: ndarray
        If `matrix` is an array with shape (m, k), then `lns` will be an array
        with shape (n, m), where n is the estimated dimension of the
        left nullspace of `matrix`.  The rows of `lns` are a basis for the
        left nullspace; each element in numpy.dot(lns A) will be
        approximately zero.

    See Also
    --------
    nullspace: Base function.

    Notes
    -----
    If both `atol` and `rtol` are positive, the combined tolerance is the
    maximum of the two; that is::
        tol = max(atol, rtol * smax)
    Singular values smaller than `tol` are considered to be zero.

    """
    lns = nullspace(matrix.T, atol, rtol).T
    return lns


def columnspace(matrix, atol=1e-13, rtol=0):
    """Compute an approximate basis for the columnspace of `matrix`.

    This function utilizes the scipy.linalg.qr method to obtain an orthogonal
    basis for the columnspace of `matrix`.

    Parameters
    ----------
    matrix: numpy.ndarray, scipy.sparse dok matrix or lil matrix,
        pandas.DataFrame, or sympy.Matrix
        The matrix to decompose. Note: `matrix` should be at most 2-D.  A 1-D
        array with length k will be treated as a 2-D with shape (1, k).
    atol: float, optional
        The absolute tolerance for a zero singular value.  Singular values
        smaller than `atol` are considered to be zero.
    rtol: float, optional
        The relative tolerance.  Singular values less than rtol*smax are
        considered to be zero, where smax is the largest singular value.

    Returns
    -------
    cs: numpy.ndarray
        If `matrix` is an array with shape (m, k), then `cs` will be an array
        with shape (m, n), where n is the estimated dimension of the
        columnspace of `matrix`. The columns of cs are a basis for the
        columnspace.

    Notes
    -----
    If both `atol` and `rtol` are positive, the combined tolerance is the
    maximum of the two; that is::
        tol = max(atol, rtol * smax)
    Singular values smaller than `tol` are considered to be zero.

    """
    matrix = _ensure_dense_matrix(matrix)
    q = linalg.qr(matrix)[0]
    cs = q[:, :matrix_rank(matrix, atol, rtol)]

    # Apply zero singular value tolerance
    s = linalg.svd(matrix, compute_uv=False)
    tol = max(atol, rtol * s[0])
    for i, row in enumerate(cs):
        for j, val in enumerate(row):
            if abs(val) <= tol:
                cs[i, j] = 0.

    return cs


def rowspace(matrix, atol=1e-13, rtol=0):
    """Compute an approximate basis for the columnspace of `matrix`.

    This function utilizes the scipy.linalg.qr method to obtain an orthogonal
    basis for the columnspace of `matrix`.

    Parameters
    ----------
    matrix: numpy.ndarray, scipy.sparse dok matrix or lil matrix,
        pandas.DataFrame, or sympy.Matrix
        The matrix to decompose. Note: `matrix` should be at most 2-D.  A 1-D
        array with length k will be treated as a 2-D with shape (1, k).
    atol: float, optional
        The absolute tolerance for a zero singular value.  Singular values
        smaller than `atol` are considered to be zero.
    rtol: float, optional
        The relative tolerance.  Singular values less than rtol*smax are
        considered to be zero, where smax is the largest singular value.

    Returns
    -------
    rs: numpy.ndarray
        If `matrix` is an array with shape (m, k), then `rs` will be an array
        with shape (m, n), where n is the estimated dimension of the rowspace
        of `matrix`. The columns of rs are a basis for the rowspace.

    See Also
    --------
    columnspace: Base function.

    Notes
    -----
    If both `atol` and `rtol` are positive, the combined tolerance is the
    maximum of the two; that is::
        tol = max(atol, rtol * smax)
    Singular values smaller than `tol` are considered to be zero.

    """
    rs = columnspace(matrix.T, atol, rtol)
    return rs


def matrix_rank(matrix, atol=1e-13, rtol=0):
    """Estimate the rank (i.e. the dimension of the nullspace) of a matrix.

    The algorithm used by this function is based on the singular value
    decomposition of `matrix`.

    Parameters
    ----------
    matrix: numpy.ndarray, scipy.sparse dok matrix or lil matrix,
        pandas.DataFrame, or sympy.Matrix
        The matrix to obtain the rank of. Note: `matrix` should be at most 2-D.
        A 1-D array with length k will be treated as a 2-D with shape (1, k).
    atol: float, optional
        The absolute tolerance for a zero singular value.  Singular values
        smaller than `atol` are considered to be zero.
    rtol: float, optional
        The relative tolerance.  Singular values less than rtol*smax are
        considered to be zero, where smax is the largest singular value.

    Returns
    -------
    rank: int
        The estimated rank of the matrix.

    Notes
    -----
    If both `atol` and `rtol` are positive, the combined tolerance is the
    maximum of the two; that is::
        tol = max(atol, rtol * smax)
    Singular values smaller than `tol` are considered to be zero.
    Taken from the scipy cookbook.

    See Also
    --------
    numpy.linalg.matrix_rank
        matrix_rank is basically the same as this function, but it does not
        provide the option of the absolute tolerance.

    """
    matrix = np.atleast_2d(_ensure_dense_matrix(matrix))
    s = linalg.svd(matrix, compute_uv=False)
    tol = max(atol, rtol * s[0])
    rank = int((s >= tol).sum())
    return rank


def svd(matrix, **kwargs):
    """Get the singular value decomposition of `matrix`.

    `kwargs`` are passed on to ``scipy.linalg.svd``
    See documentation for ``scipy.linalg.svd`` for more details.

    Parameters
    ----------
    matrix: numpy.ndarray, scipy.sparse dok matrix or lil matrix,
        pandas.DataFrame, or sympy.Matrix
        The matrix to decompose. Note: `matrix` should be at most 2-D.  A 1-D
        array with length k will be treated as a 2-D with shape (1, k).


    Returns
    -------
    matrix of the same type as `matrix`.

    See Also
    --------
    scipy.linalg.svd: Base function.
        svd and its arguments are the same as this method. The only difference
        is that matrices of various formats are converted in order to ensure
        the correct input for scipy.linalg.svd.

    """
    matrix = _ensure_dense_matrix(matrix)
    return linalg.svd(matrix, **kwargs)


def eig(matrix, left=False, right=False, **kwargs):
    """Get the eigenvalues of `matrix`.

    `kwargs`` are passed on to ``scipy.linalg.eig``
    See documentation for ``scipy.linalg.eig`` for more details.

    Parameters
    ----------
    matrix: numpy.ndarray, scipy.sparse dok matrix or lil matrix,
        pandas.DataFrame, or sympy.Matrix
    left: bool, optional
        Whether to calculate and return left eigenvectors. Default is False.
    right: bool, optional
        Whether to calculate and return right eigenvectors. Default is True.

    Returns
    -------
    w: (M,) or (2, M) double or complex ndarray
        The eigenvalues, each repeated according to its multiplicity.
        The shape is (M,) unless homogeneous_eigvals=True.
    vl: (M, M) double or complex ndarray
        The normalized left eigenvector corresponding to the eigenvalue w[i]
        is the column vl[:,i]. Only returned if left=True.
    vr: (M, M) double or complex ndarray
        The normalized right eigenvector corresponding to the eigenvalue w[i]
        is the column vr[:,i]. Only returned if right=True.

    See Also
    --------
    scipy.linalg.eig: Base function.
        svd and its arguments are the same as this method. The only difference
        is that matrices of various formats are converted in order to ensure
        the correct input for scipy.linalg.eig.

    """
    matrix = _ensure_dense_matrix(matrix)
    return linalg.eig(matrix, left=left, right=right, **kwargs)


def _ensure_dense_matrix(matrix):
    """Ensure matrix is dense before performing linear algebra operations.

    Warnings
    --------
    This method is intended for internal use only.
    """
    if isinstance(matrix, (np.ndarray, pd.DataFrame)):
        pass
    elif isinstance(matrix, (dok_matrix, lil_matrix)):
        matrix = matrix.toarray()
    elif isinstance(matrix, sym.Matrix):
        try:
            matrix = np.array(matrix).astype(np.float64)
        except TypeError:
            raise ValueError("Cannot have sympy symbols in the matrix. Try "
                             "substituting numerical values in first")
    else:
        raise TypeError("Matrix must be one of the following formats: "
                        "numpy.ndarray, scipy.dok_matrix, scipy.lil_matrix, "
                        "pandas.DataFrame, or sympy.Matrix.")
    return matrix
