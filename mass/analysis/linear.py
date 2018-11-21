# -*- coding: utf-8 -*-
"""TODO Module Docstrings."""
from __future__ import absolute_import

import re
import warnings

from mass.util.util import convert_matrix

import numpy as np

import pandas as pd

from scipy import linalg
from scipy.sparse import dok_matrix, lil_matrix

from six import iteritems

import sympy as sym

# Global
_T_SYM = sym.Symbol("t")
# Precompiled re for 'external' metabolites
ext_metab_re = re.compile("\_e")


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
            met_func = sym.Function(str(met))(_T_SYM)
            gradient_mat[r_ind(rxn), m_ind(met)] = rate.diff(met_func)

    # Get values for substitution
    if use_concentration_values or use_parameter_values:
        values = {}
        if use_parameter_values:
            model_parameters = model.parameters
            # fixed_concs = model_parameters.pop("Fixed")

            for key, dictionary in iteritems(model_parameters):
                values.update({sym.Symbol(str(k)): v
                              for k, v in iteritems(dictionary)})

        if use_concentration_values:
            values.update({sym.Function(str(k))(_T_SYM): v
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
    kappa_mat: matrix of type 'matrix_type'
        The kappa matrix for the model returned as the given matrix_type.

    """
    gradient_mat = gradient(model, use_parameter_values,
                            use_concentration_values, matrix_type)
    kappa_mat = sym.diag(*[gradient_mat[row, :].norm()
                           for row in range(gradient_mat.rows)])
    kappa_mat = kappa_mat.subs({sym.nan: sym.S.Zero})
    kappa_mat = convert_matrix(kappa_mat, matrix_type=matrix_type,
                               dtype=np.float64,
                               row_ids=[r.id for r in model.reactions],
                               col_ids=[r.id for r in model.reactions])
    return kappa_mat


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
    gamma_mat: matrix of type 'matrix_type'
        The gamma matrix for the model returned as the given matrix_type.

    """
    gradient_mat = gradient(model, use_parameter_values,
                            use_concentration_values, matrix_type)
    gamma_mat = sym.Matrix([gradient_mat[row, :].normalized()
                            for row in range(gradient_mat.rows)])
    gamma_mat = gamma_mat.subs({sym.nan: sym.S.Zero})
    gamma_mat = convert_matrix(gamma_mat, matrix_type=matrix_type,
                               dtype=np.float64,
                               row_ids=[r.id for r in model.reactions],
                               col_ids=[m.id for m in model.metabolites])
    return gamma_mat


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
    jacobian_mat: matrix of type 'matrix_type'
        The jacobian matrix for the model returned as the given matrix_type.

    """
    if jacobian_type not in {"Jx", "Jv"}:
        raise ValueError("jacobian_type must be either 'Jx' or Jv'")

    gradient_mat = gradient(model, use_parameter_values,
                            use_concentration_values, matrix_type="symbolic")
    stoich_mat = model._mk_stoich_matrix(matrix_type="symbolic",
                                         update_model=False)
    if "Jx" == jacobian_type:
        jacobian_mat = stoich_mat * gradient_mat
        identifiers = [m.id for m in model.metabolites]
    else:
        jacobian_mat = gradient_mat * stoich_mat
        identifiers = [r.id for r in model.reactions]

    jacobian_mat = convert_matrix(jacobian_mat, matrix_type=matrix_type,
                                  dtype=np.float64,
                                  row_ids=identifiers,
                                  col_ids=identifiers)
    return jacobian_mat


def nullspace(A, atol=1e-13, rtol=0):
    """Compute an approximate basis for the nullspace of A.

    The algorithm used by this function is based on the singular value
    decomposition of `A`.

    Parameters
    ----------
    A: numpy.ndarray, scipy.sparse dok matrix or lil matrix, pandas.DataFrame
        or sympy.Matrix
        Note: 'A' should be at most 2-D.  A 1-D array with length k will be
        treated as a 2-D with shape (1, k)
    atol: float, optional
        The absolute tolerance for a zero singular value.  Singular values
        smaller than `atol` are considered to be zero.
    rtol: float, optional
        The relative tolerance.  Singular values less than rtol*smax are
        considered to be zero, where smax is the largest singular value.

    Returns
    -------
    ns: ndarray
        If `A` is an array with shape (m, k), then `ns` will be an array
        with shape (k, n), where n is the estimated dimension of the
        nullspace of `A`.  The columns of `ns` are a basis for the
        nullspace; each element in numpy.dot(A, ns) will be approximately
        zero.

    Notes
    -----
    If both `atol` and `rtol` are positive, the combined tolerance is the
    maximum of the two; that is::
        tol = max(atol, rtol * smax)
    Singular values smaller than `tol` are considered to be zero.

    """
    A = np.atleast_2d(_ensure_dense_mat(A))
    s, vh = linalg.svd(A)[1:]
    tol = max(atol, rtol * s[0])
    nnz = (s >= tol).sum()
    ns = vh[nnz:].conj().T

    # Apply zero singular value tolerance
    for i, row in enumerate(ns):
        for j, val in enumerate(row):
            if abs(val) <= tol:
                ns[i, j] = 0.
    return ns


def left_nullspace(A, atol=1e-13, rtol=0):
    """Compute an approximate basis for the left nullspace of A.

    The algorithm used by this function is based on the singular value
    decomposition of `A`.

    Parameters
    ----------
    A: numpy.ndarray, scipy.sparse dok matrix or lil matrix, pandas.DataFrame
        or sympy.Matrix
        Note: 'A' should be at most 2-D.  A 1-D array with length k will be
        treated as a 2-D with shape (1, k)
    atol: float, optional
        The absolute tolerance for a zero singular value.  Singular values
        smaller than `atol` are considered to be zero.
    rtol: float, optional
        The relative tolerance.  Singular values less than rtol*smax are
        considered to be zero, where smax is the largest singular value.

    Returns
    -------
    lns: ndarray
        If `A` is an array with shape (m, k), then `lns` will be an array
        with shape (n, m), where n is the estimated dimension of the
        left nullspace of `A`.  The rows of `lns` are a basis for the
        left nullspace; each element in numpy.dot(lns A) will be
        approximately zero.

    See ALso
    --------
    nullspace: Base function.

    Notes
    -----
    If both `atol` and `rtol` are positive, the combined tolerance is the
    maximum of the two; that is::
        tol = max(atol, rtol * smax)
    Singular values smaller than `tol` are considered to be zero.

    """
    lns = nullspace(A.T, atol, rtol).T
    return lns


def columnspace(A, atol=1e-13, rtol=0):
    """Compute an approximate basis for the columnspace of A.

    This function utilizes the scipy.linalg.qr method to obtain an orthogonal
    basis for the columnspace of A.

    Parameters
    ----------
    A: numpy.ndarray, scipy.sparse dok matrix or lil matrix, pandas.DataFrame
        or sympy.Matrix
        Note: 'A' should be at most 2-D.  A 1-D array with length k will be
        treated as a 2-D with shape (1, k)
    atol: float, optional
        The absolute tolerance for a zero singular value.  Singular values
        smaller than `atol` are considered to be zero.
    rtol: float, optional
        The relative tolerance.  Singular values less than rtol*smax are
        considered to be zero, where smax is the largest singular value.

    Returns
    -------
    cs: numpy.ndarray
        If `A` is an array with shape (m, k), then `cs` will be an array
        with shape (m, n), where n is the estimated dimension of the
        columnspace of `A`. The columns of cs are a basis for the columnspace.

    Notes
    -----
    If both `atol` and `rtol` are positive, the combined tolerance is the
    maximum of the two; that is::
        tol = max(atol, rtol * smax)
    Singular values smaller than `tol` are considered to be zero.

    """
    A = _ensure_dense_mat(A)
    q = linalg.qr(A)[0]
    cs = q[:, :matrix_rank(A, atol, rtol)]

    # Apply zero singular value tolerance
    s = linalg.svd(A, compute_uv=False)
    tol = max(atol, rtol * s[0])
    for i, row in enumerate(cs):
        for j, val in enumerate(row):
            if abs(val) <= tol:
                cs[i, j] = 0.

    return cs


def rowspace(A, atol=1e-13, rtol=0):
    """Compute an approximate basis for the columnspace of A.

    This function utilizes the scipy.linalg.qr method to obtain an orthogonal
    basis for the columnspace of A.

    Parameters
    ----------
    A: numpy.ndarray, scipy.sparse dok matrix or lil matrix, pandas.DataFrame
        or sympy.Matrix
        Note: 'A' should be at most 2-D.  A 1-D array with length k will be
        treated as a 2-D with shape (1, k)
    atol: float, optional
        The absolute tolerance for a zero singular value.  Singular values
        smaller than `atol` are considered to be zero.
    rtol: float, optional
        The relative tolerance.  Singular values less than rtol*smax are
        considered to be zero, where smax is the largest singular value.

    Returns
    -------
    rs: numpy.ndarray
        If `A` is an array with shape (m, k), then `rs` will be an array
        with shape (m, n), where n is the estimated dimension of the rowspace
        of `A`. The columns of rs are a basis for the rowspace.

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
    rs = columnspace(A.T, atol, rtol)
    return rs


def matrix_rank(A, atol=1e-13, rtol=0):
    """Estimate the rank (i.e. the dimension of the nullspace) of a matrix.

    The algorithm used by this function is based on the singular value
    decomposition of `A`.

    Parameters
    ----------
    A: ndarray
        A should be at most 2-D.  A 1-D array with length n will be treated
        as a 2-D with shape (1, n)
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
    A = np.atleast_2d(_ensure_dense_mat(A))
    s = linalg.svd(A, compute_uv=False)
    tol = max(atol, rtol * s[0])
    rank = int((s >= tol).sum())
    return rank


def svd(A, **kwargs):
    """Get the singular value decomposition of 'A'.

    `kwargs`` are passed on to ``scipy.linalg.svd``
    See documentation for ``scipy.linalg.svd`` for more details.

    Parameters
    ----------
    A: numpy.ndarray, scipy.sparse dok matrix or lil matrix, pandas.DataFrame
        or sympy.Matrix

    Returns
    -------
    matrix of the same type as 'A'

    See Also
    --------
    scipy.linalg.svd: Base function.
        svd and its arguments are the same as this method. The only difference
        is that matrices of various formats are converted in order to ensure
        the correct input for scipy.linalg.svd.

    """
    A = _ensure_dense_mat(A)
    return linalg.svd(A, **kwargs)


def eig(A, left=False, right=False, **kwargs):
    """Get the eigenvalues of 'A'.

    `kwargs`` are passed on to ``scipy.linalg.eig``
    See documentation for ``scipy.linalg.eig`` for more details.

    Parameters
    ----------
    A: numpy.ndarray, scipy.sparse dok matrix or lil matrix, pandas.DataFrame
        or sympy.Matrix
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
    A = _ensure_dense_mat(A)
    return linalg.eig(A, left=left, right=right, **kwargs)


def temporal_decomposition(model, jacobian_type='metabolite', eigtol=1e-10,
                           zerotol=1e-10, as_percents=True, remove_imag=True,
                           mode_equations=True, dynamic_invariants=True):
    """Perform a temporal decomposition of the jacobian matrix of a model.

    The timescales (ts) and the modal matrix (M) are returned, where ts[i] is
    the time constant for M[i]. Will also return the equations for each mode as
    a sympy expression where ts[i] is the time constant for m[i] if specified.

    Parameters
    ----------
    model: mass.MassModel
        The MassModel object to construct the matrix for
    jacobian_type: {'Jx', 'Jv'}
        Whether to obtain the jacobian matrix with respect to metabolites (Jx),
        or to obtain the jacobian matrix with respect to the reactions (Jv).
        Default is 'Jx'.
    eigtol: float, optional
        The absolute tolerance for a zero singular value of an eigenvalue.
        Singular values smaller than `eigtol` are considered to be zero.
    zerotol: float, optional
        The absolute tolerance for a zero singular value in the modal matrix M.
        Singular values smaller than `zerotol` are considered to be zero.
    as_percents: bool, optional
        If True, will normalize rows of the modal matrix such that the
        sum(abs(M[i])) is equal to 1.
    remove_imag: bool, optional
        If True, will remove the complex part of the the values.
    mode_equations: bool, optional
        If True, will return the equations for each mode as a sympy expression
    dynamic_invariants: bool, optional
        If True, will include the dynamically invariant pools in the returned
        modal matrix and mode equations .

    Returns
    -------
    ts: numpy.array
        A numpy array  where ts[i] is time constant for M[i].
        Time invariants will have a time constant of np.inf.
    M: numpy.array
        A numpy array representing the modal matrix where M[i] is the row of
        the matrix that corresponds to the time constant ts[i].
    m: numpy.array
            A numpy array representing the mode equations where m[i] is the
            equation that corresponds to the time constant ts[i].

    """
    # TODO Refactor function in a new timescale_decomposition class later
    J = jacobian(model, jacobian_type=jacobian_type, strip_time=True,
                 use_parameter_values=True, use_concentration_values=True,
                 matrix_type="dense")
    # get the eigenvalues and matrix of eigenrows
    w, vr = eig(J, left=False, right=True)
    rank = matrix_rank(J, atol=eigtol)
    # Define eigenvalues with real parts below the tolerance as zero.
    for i, eigenvalue in enumerate(w):
        if abs(eigenvalue.real) <= eigtol:
            w[i] = 0
    # Get indicies of eigenvalues and sort timescales and correspodning rows
    # of the model matrix from fastest to slowest
    ts = -1/np.sort(w)[:rank]
    if dynamic_invariants:
        indices = np.argsort(w)
        ts = np.concatenate((ts, [np.inf]*(len(w)-rank)))
    else:
        indices = np.argsort(w)[:rank]
    M = np.array([linalg.inv(vr)[index] for index in indices])
    # Normalize rows by largest weight
    for i, row in enumerate(M):
        row = row//max(abs(row))
        # Convert to percents if specificed
        if as_percents:
            row = row/sum(abs(row))
        for j, val in enumerate(row):
            val = np.real_if_close(val, tol=1/zerotol)
            if abs(val.real) <= zerotol:
                val = 0.
            row[j] = val
        M[i] = row
    if remove_imag:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            M = M.astype(float)
            ts = ts.astype(float)
    # Get mode equations if specified
    if mode_equations:
        metab_funcs = [sym.Symbol(m.id)(_T_SYM) for m in model.metabolites]
        m = [0]*len(ts)
        e = abs(np.floor(np.log10(np.abs(zerotol))).astype(int))
        for i, row in enumerate(M):
            m[i] = sum([round(val, e)*metab_funcs[i]
                        for i, val in enumerate(row)])
        return ts, M, np.array(m)

    return ts, M


def _ensure_dense_mat(A):
    """Ensure matrix is dense before performing linear algebra operations.

    Warnings
    --------
    This method is intended for internal use only.
    """
    if isinstance(A, (np.ndarray, pd.DataFrame)):
        pass
    elif isinstance(A, (dok_matrix, lil_matrix)):
        A = A.toarray()
    elif isinstance(A, sym.Matrix):
        try:
            A = np.array(A).astype(np.float64)
        except TypeError:
            raise ValueError("Cannot have sympy symbols in the matrix. Try "
                             "substituting numerical values in first")
    else:
        raise TypeError("Matrix must be one of the following formats: "
                        "numpy.ndarray, scipy.dok_matrix, scipy.lil_matrix, "
                        "pandas.DataFrame, or sympy.Matrix.")
    return A
