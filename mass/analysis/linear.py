# -*- coding: utf-8 -*-
"""
Containins basic matrix operations that can be applied to a :mod:`mass` model.

To assist with various matrix operations using different packages, the
following values can be provided to the ``matrix_type`` argument to set the
return type of the output. Valid matrix types include:

* ``'dense'`` for a :class:`numpy.ndarray`
* ``'dok'`` for a :class:`scipy.sparse.dok_matrix`
* ``'lil'`` for a :class:`scipy.sparse.lil_matrix`
* ``'DataFrame'`` for a :class:`pandas.DataFrame`
* ``'symbolic'`` for a
  :class:`sympy.MutableDenseMatrix <sympy.matrices.dense.MutableDenseMatrix>`

For all matrix types, species (excluding genes) are the row indicies and
reactions are the column indicies.

There are also several methods that are nearly identical to :mod:`scipy.linalg`
methods, with the main exception being that matrix conversions are performed
beforehand to ensure that valid input is passed to the :mod:`scipy` method.
These methods include:

* :func:`~mass.analysis.linear.svd`
* :func:`~mass.analysis.linear.eig`

"""
import numpy as np

import pandas as pd

from scipy import linalg
from scipy.sparse import dok_matrix, lil_matrix

from six import iteritems

import sympy as sym

from mass.core.mass_configuration import MassConfiguration
from mass.util.expressions import _mk_met_func
from mass.util.util import convert_matrix

MASSCONFIGURATION = MassConfiguration()


# Public
def gradient(model, use_parameter_values=True, use_concentration_values=True,
             matrix_type="dense"):
    """Create the gradient matrix for a given model.

    Parameters
    ----------
    model : MassModel
        The :class:`~.MassModel` to construct the matrix for.
    use_parameter_values : bool
        Whether to substitute the numerical values for parameters into the
        matrix. If ``True`` then numerical values of the kinetic parameters
        are substituted into the matrix. Otherwise parameters in the matrix
        are left as symbols.
    use_concentration_values : bool
        Whether to substitute the numerical values for concentrations into the
        matrix. If ``True`` then numerical values of the initial conditions
        are substituted into the matrix. Otherwise species concentrations in
        the matrix are left as symbols.
    matrix_type : str
        A string identifiying the desired format for the returned matrix.
        Default is ``'dense'``.  See the :mod:`~.linear` module documentation
        for more information on the ``matrix_type``

    Returns
    -------
    matrix of type ``matrix_type``
        The gradient matrix for the model.

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

    Notes
    -----
    The kappa matrix is the diagnolization of the norms for the rows in the
    gradient matrix.

    Parameters
    ----------
    model : MassModel
        The :class:`~.MassModel` to construct the matrix for.
    use_parameter_values : bool
        Whether to substitute the numerical values for parameters into the
        matrix. If ``True`` then numerical values of the kinetic parameters
        are substituted into the matrix. Otherwise parameters in the matrix
        are left as symbols.
    use_concentration_values : bool
        Whether to substitute the numerical values for concentrations into the
        matrix. If ``True`` then numerical values of the initial conditions
        are substituted into the matrix. Otherwise species concentrations in
        the matrix are left as symbols.
    matrix_type : str
        A string identifiying the desired format for the returned matrix.
        Default is ``'dense'``.  See the :mod:`~.linear` module documentation
        for more information on the ``matrix_type``.

    Returns
    -------
    matrix of type ``matrix_type``
        The kappa matrix for the model.

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

    Notes
    -----
    The gamma matrix is composed of the 1-norms of the gradient matrix.

    Parameters
    ----------
    model : MassModel
        The :class:`~.MassModel` to construct the matrix for.
    use_parameter_values : bool
        Whether to substitute the numerical values for parameters into the
        matrix. If ``True`` then numerical values of the kinetic parameters
        are substituted into the matrix. Otherwise parameters in the matrix
        are left as symbols.
    use_concentration_values : bool
        Whether to substitute the numerical values for concentrations into
        the matrix. If ``True`` then numerical values of the initial
        conditions are substituted into the matrix. Otherwise species
        concentrations in the matrix are left as symbols.
    matrix_type : str
        A string identifiying the desired format for the returned matrix.
        Default is ``'dense'``.  See the :mod:`~.linear` module documentation
        for more information on the ``matrix_type``.

    Returns
    -------
    matrix of type ``matrix_type``
        The gamma matrix for the model.

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
    model : MassModel
        The :class:`~.MassModel` to construct the matrix for.
    jacobian_type: str
        Either the string ``'Jx'`` to obtain the jacobian matrix with respect
        to species, or the string ``'Jv'`` to obtain the jacobian matrix
        with respect to the reactions. Default is ``'Jx'``.
    use_parameter_values : bool
        Whether to substitute the numerical values for parameters into the
        matrix. If ``True`` then numerical values of the kinetic parameters
        are substituted into the matrix. Otherwise parameters in the matrix
        are left as symbols.
    use_concentration_values : bool
        Whether to substitute the numerical values for concentrations into the
        matrix. If ``True`` then numerical values of the initial conditions
        are substituted into the matrix. Otherwise species concentrations in
        the matrix are left as symbols.
    matrix_type : str
        A string identifiying the desired format for the returned matrix.
        Default is ``'dense'``.  See the :mod:`~.linear` module documentation
        for more information on the ``matrix_type``.

    Returns
    -------
    matrix of type ``matrix_type``
        The jacobian matrix for the model.

    """
    if jacobian_type not in {"Jx", "Jv"}:
        raise ValueError("jacobian_type must be either 'Jx' or 'Jv'")

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


def nullspace(matrix, atol=1e-13, rtol=0, decimal_precision=False):
    """Compute an approximate basis for the nullspace of a matrix.

    The algorithm used by this function is based on singular value
    decomposition.

    Notes
    -----
    * If both ``atol`` and ``rtol`` are positive, the combined tolerance
      is the maximum of the two; that is::

        tol = max(atol, rtol * smax)

      Singular values smaller than ``tol`` are considered to be zero.

    * Similar to the :class:`cobra.util.array.nullspace` function, but includes
      utlization of the  :attr:`~.MassBaseConfiguration.decimal_precision` in
      the :class:`.MassConfiguration` and sets values below the tolerance to 0.

    * Taken from the numpy cookbook and extended.

    Parameters
    ----------
    matrix : array-like
        The matrix to decompose. The matrix should be at most 2-D. A 1-D
        array with length ``k`` will be treated as a 2-D with shape ``(1, k)``.
    atol : float
        The absolute tolerance for a zero singular value.  Singular values
        smaller than ``atol`` are considered to be zero.
    rtol : float
        The relative tolerance.  Singular values less than ``rtol * smax`` are
        considered to be zero, where ``smax`` is the largest singular value.
    decimal_precision : bool
        Whether to apply the :attr:`~.MassBaseConfiguration.decimal_precision`
        set in the :class:`.MassConfiguration` to the nullspace values
        before comparing to the tolerance. Default is ``False``.

    Returns
    -------
    ns : numpy.ndarray
        If ``matrix`` is an array with shape ``(m, k)``, then ``ns`` will be
        an array with shape ``(k, n)``, where ``n`` is the estimated dimension
        of the nullspace of ``matrix``.  The columns of ``ns`` are a basis for
        the nullspace; each element in the dot product of the matrix
        and the nullspace will be approximately 0.

    """
    matrix = np.atleast_2d(_ensure_dense_matrix(matrix))
    s, vh = linalg.svd(matrix)[1:]
    tol = max(atol, rtol * s[0])
    nnz = (s >= tol).sum()
    ns = vh[nnz:].conj().T

    # Apply zero singular value tolerance
    for i, row in enumerate(ns):
        for j, val in enumerate(row):
            if decimal_precision\
               and MASSCONFIGURATION.decimal_precision is not None:
                val = round(val, MASSCONFIGURATION.decimal_precision)
            if abs(val) <= tol:
                ns[i, j] = 0.
    return ns


def left_nullspace(matrix, atol=1e-13, rtol=0, decimal_precision=False):
    """Compute an approximate basis for the left nullspace of a matrix.

    The algorithm used by this function is based on singular value
    decomposition.

    Notes
    -----
    If both ``atol`` and ``rtol`` are positive, the combined tolerance is the
    maximum of the two; that is::

        tol = max(atol, rtol * smax)

    Singular values smaller than ``tol`` are considered to be zero.

    Parameters
    ----------
    matrix : array-like
        The matrix to decompose. The matrix should be at most 2-D. A 1-D
        array with length ``k`` will be treated as a 2-D with shape ``(1, k)``.
    atol : float
        The absolute tolerance for a zero singular value.  Singular values
        smaller than ``atol`` are considered to be zero.
    rtol : float
        The relative tolerance.  Singular values less than ``rtol * smax`` are
        considered to be zero, where ``smax`` is the largest singular value.
    decimal_precision : bool
        Whether to apply the :attr:`~.MassBaseConfiguration.decimal_precision`
        set in the :class:`.MassConfiguration` to the left nullspace values
        before comparing to the tolerance. Default is ``False``.

    Returns
    -------
    lns : numpy.ndarray
        If ``matrix`` is an array with shape ``(m, k)``, then ``lns`` will be
        an array with shape ``(n, m)``, where ``n`` is the estimated dimension
        of the left nullspace of ``matrix``.  The rows of ``lns`` are a basis
        for the left nullspace; each element in the dot product of the matrix
        and the left nullspace will be approximately 0.

    See Also
    --------
    :func:`nullspace` : Base function.

    """
    lns = nullspace(matrix.T, atol, rtol, decimal_precision).T
    return lns


def columnspace(matrix, atol=1e-13, rtol=0, decimal_precision=False):
    """Compute an approximate basis for the columnspace of a matrix.

    This function utilizes the :func:`scipy.linalg.qr` function to obtain an
    orthogonal basis for the columnspace of the matrix.

    Notes
    -----
    If both ``atol`` and ``rtol`` are positive, the combined tolerance is the
    maximum of the two; that is::

        tol = max(atol, rtol * smax)

    Singular values smaller than ``tol`` are considered to be zero.

    Parameters
    ----------
    matrix : array-like
        The matrix to decompose. The matrix should be at most 2-D. A 1-D
        array with length ``k`` will be treated as a 2-D with shape ``(1, k)``.
    atol : float
        The absolute tolerance for a zero singular value.  Singular values
        smaller than ``atol`` are considered to be zero.
    rtol : float
        The relative tolerance.  Singular values less than ``rtol * smax`` are
        considered to be zero, where ``smax`` is the largest singular value.
    decimal_precision : bool
        Whether to apply the :attr:`~.MassBaseConfiguration.decimal_precision`
        set in the :class:`.MassConfiguration` to the columnspace values
        before comparing to the tolerance. Default is ``False``.

    Returns
    -------
    cs : numpy.ndarray
        If ``matrix`` is an array with shape ``(m, k)``, then ``cs`` will be
        an array with shape ``(m, n)``, where ``n`` is the estimated dimension
        of the columnspace of ``matrix``. The columns of ``cs`` are a basis
        for the columnspace.

    """
    matrix = _ensure_dense_matrix(matrix)
    q = linalg.qr(matrix)[0]
    cs = q[:, :matrix_rank(matrix, atol, rtol)]

    # Apply zero singular value tolerance
    s = linalg.svd(matrix, compute_uv=False)
    tol = max(atol, rtol * s[0])

    # Apply zero singular value tolerance
    for i, row in enumerate(cs):
        for j, val in enumerate(row):
            if decimal_precision\
               and MASSCONFIGURATION.decimal_precision is not None:
                val = round(val, MASSCONFIGURATION.decimal_precision)
            if abs(val) <= tol:
                cs[i, j] = 0.

    return cs


def rowspace(matrix, atol=1e-13, rtol=0, decimal_precision=False):
    """Compute an approximate basis for the rowspace of a matrix.

    This function utilizes the :func:`scipy.linalg.qr` function to obtain an
    orthogonal basis for the rowspace of the matrix.

    Notes
    -----
    If both ``atol`` and ``rtol`` are positive, the combined tolerance is the
    maximum of the two; that is::

        tol = max(atol, rtol * smax)

    Singular values smaller than ``tol`` are considered to be zero.

    Parameters
    ----------
    matrix : array-like
        The matrix to decompose. The matrix should be at most 2-D. A 1-D
        array with length ``k`` will be treated as a 2-D with shape ``(1, k)``.
    atol : float
        The absolute tolerance for a zero singular value.  Singular values
        smaller than ``atol`` are considered to be zero.
    rtol : float
        The relative tolerance.  Singular values less than ``rtol * smax`` are
        considered to be zero, where ``smax`` is the largest singular value.
    decimal_precision : bool
        Whether to apply the :attr:`~.MassBaseConfiguration.decimal_precision`
        set in the :class:`.MassConfiguration` to the rowspace values
        before comparing to the tolerance. Default is ``False``.

    Returns
    -------
    rs : numpy.ndarray
        If ``matrix`` is an array with shape ``(m, k)``, then ``rs`` will be
        an array with shape ``(n, k)``, where ``n`` is the estimated dimension
        of the rowspace of ``matrix``. The columns of ``rs`` are a basis for
        the rowspace.

    See Also
    --------
    :func:`columnspace` : Base function.

    """
    rs = columnspace(matrix.T, atol, rtol, decimal_precision)
    return rs


def matrix_rank(matrix, atol=1e-13, rtol=0):
    """Estimate the rank (i.e. the dimension of the nullspace) of a matrix.

    The algorithm used by this function is based on singular value
    decomposition. Taken from the :mod:`scipy` cookbook.

    Notes
    -----
    If both ``atol`` and ``rtol`` are positive, the combined tolerance is the
    maximum of the two; that is::

        tol = max(atol, rtol * smax)

    Singular values smaller than ``tol`` are considered to be zero.

    Parameters
    ----------
    matrix : array-like
        The matrix to obtain the rank for. The matrix should be at most
        2-D. A 1-D array with length ``k`` will be treated as a 2-D with
        shape ``(1, k)``.
    atol : float
        The absolute tolerance for a zero singular value.  Singular values
        smaller than ``atol`` are considered to be zero.
    rtol : float
        The relative tolerance.  Singular values less than ``rtol * smax`` are
        considered to be zero, where ``smax`` is the largest singular value.

    Returns
    -------
    rank : int
        The estimated rank of the matrix.

    See Also
    --------
    :func:`numpy.linalg.matrix_rank`
        :func:`mass.analysis.linear.matrix_rank` is nearly identical to this
        function, but it does not provide the option of the absolute tolerance.

    """
    matrix = np.atleast_2d(_ensure_dense_matrix(matrix))
    s = linalg.svd(matrix, compute_uv=False)
    tol = max(atol, rtol * s[0])
    rank = int((s >= tol).sum())
    return rank


def svd(matrix, **kwargs):
    """Get the singular value decomposition of a matrix.

    ``kwargs`` are passed on to :func:`scipy.linalg.svd`.

    Parameters
    ----------
    matrix : array-like
        The matrix to decompose. The matrix should be at most 2-D. A 1-D
        array with length ``k`` will be treated as a 2-D with shape ``(1, k)``.

    Returns
    -------
    U : ndarray
        Unitary matrix having left singular vectors as columns.
        Of shape ``(M, M)`` or ``(M, K)``, depending on `full_matrices`.
    s : ndarray
        The singular values, sorted in non-increasing order.
        Of shape ``(K, )``, with ``K = min(M, N)``.
    Vh : ndarray
        Unitary matrix having right singular vectors as rows.
        Of shape ``(N, N)`` or ``(K, N)`` depending on `full_matrices`.
    For ``compute_uv=False``, only ``s`` is returned.

    See Also
    --------
    :func:`scipy.linalg.svd`: Base function.

    """
    matrix = _ensure_dense_matrix(matrix)
    return linalg.svd(matrix, **kwargs)


def eig(matrix, left=False, right=False, **kwargs):
    """Get the eigenvalues of a matrix.

    ``kwargs`` are passed on to :func:`scipy.linalg.eig`

    Parameters
    ----------
    matrix : array-like
        The matrix to decompose. The matrix should be at most 2-D. A 1-D
        array with length ``k`` will be treated as a 2-D with shape ``(1, k)``.
    left : bool
        Whether to calculate and return left eigenvectors.
        Default is ``False``.
    right : bool
        Whether to calculate and return right eigenvectors.
        Default is ``True``.

    Returns
    -------
    w : (M, ) double or complex ndarray
        The eigenvalues, each repeated according to its multiplicity.
    vl : (M, M) double or complex ndarray
        The normalized left eigenvector corresponding to the eigenvalue
        ``w[i]`` is the column v[:,i]. Only returned if ``left=True``.
    vr : (M, M) double or complex ndarray
        The normalized right eigenvector corresponding to the eigenvalue
        ``w[i]`` is the column ``vr[:,i]``.  Only returned if ``right=True``.

    See Also
    --------
    :func:`scipy.linalg.eig`: Base function.

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


__all__ = (
    "columnspace", "eig", "gamma", "gradient", "jacobian", "kappa",
    "left_nullspace", "matrix_rank", "nullspace", "rowspace", "svd")
