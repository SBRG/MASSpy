# -*- coding: utf-8 -*-
"""Contains utility functions to assist in various :mod:`mass` functions."""
import logging
import re
import warnings

from cobra import DictList

from depinfo import print_dependencies

import numpy as np

import pandas as pd

from scipy.sparse import dok_matrix, lil_matrix

from six import integer_types, iteritems, string_types

import sympy as sym

_MATRIX_TYPES = ["dense", "dok", "lil", "DataFrame", "symbolic"]


# Public
def show_versions():
    """Print dependency information."""
    print_dependencies("masspy")


def ensure_iterable(item):
    """Ensure the given item is an returned as an iterable.

    Parameters
    ----------
    item : object
        The item to ensure is returned as an iterable.

    """
    # Make list iterable if necessary
    if item is None:
        item = list()
    if not hasattr(item, "__iter__") or \
       isinstance(item, string_types):
        item = [item]

    item = list(item)
    return item


def ensure_non_negative_value(value):
    """Ensure provided value is a non-negative value, or ``None``.

    Parameters
    ----------
    value : float
        The value to ensure is non-negative

    Raises
    ------
    ValueError
        Occurs if the value is negative.

    """
    if value is None:
        pass
    elif not isinstance(value, (integer_types, float)):
        raise TypeError("Must be an int or float")
    elif value < 0.:
        raise ValueError("Must be a non-negative number")
    return value


def convert_matrix(matrix, matrix_type, dtype, row_ids=None, col_ids=None):
    """Convert a matrix to a different type.

    Parameters
    ----------
    matrix : array-like
        The matrix to convert.
    matrix_type : str
        A string identifiying the desired format for the returned matrix.
        Valid matrix types include ``'dense'``, ``'dok'``, ``'lil'``,
        ``'DataFrame'``, and ``'symbolic'`` See the :mod:`~.linear` module
        documentation for more information on the ``matrix_type``.
    dtype : data-type
        The desired array data-type for the matrix.
    row_ids : array-like
        The idenfifiers for each row. Only used if type is ``'DataFrame'``.
    col_ids : array-like
        The idenfifiers for each column. Only used if type is ``'DataFrame'``.

    Warnings
    --------
    This method is NOT the safest way to convert a matrix to another type.
    To safely convert a matrix into another type, use the ``'matrix_type'``
    argument in the method that returns the desired matrix.

    """
    if matrix_type not in _MATRIX_TYPES:
        raise ValueError("Unrecognized matrix_type.")

    # Convert the matrix type
    conversion_method_dict = dict(zip(
        _MATRIX_TYPES, [_to_dense, _to_dok, _to_lil, _to_dense, _to_dense]))

    try:
        matrix = conversion_method_dict[matrix_type](matrix)
        # Convert the dtype
        if not re.match("symbolic", matrix_type):
            if re.match("DataFrame", matrix_type):
                matrix = pd.DataFrame(matrix, index=row_ids, columns=col_ids)
            try:
                matrix = matrix.astype(dtype)
            except TypeError:
                warnings.warn("Could not cast matrix as the given dtype")
        else:
            matrix = sym.Matrix(matrix)
    except TypeError:
        warnings.warn("Could not cast matrix as the given matrix_type")

    return matrix


def get_public_attributes_and_methods(obj, exclude_parent=False):
    """Return a list of public attributes and methods for an object.

    Parameters
    ----------
    exclude_parent : bool
        If ``True``, only display public attributes and methods specific to
        the current class, excluding those inherited from the parent class.
        Overridden and extended methods are not excluded. 

    """
    all_public = [i.strip("_") for i in obj.__dict__]
    all_public += [i for i in obj.__class__.__dict__
                   if i not in all_public and not i.startswith("_")]
    if not exclude_parent:
        parent_public = get_public_attributes_and_methods(
            obj.__class__.__base__(), exclude_parent=True)
        all_public += [i for i in parent_public if i not in all_public]

    return sorted(all_public, key=str.lower)


# Internal
def _check_kwargs(default_kwargs, kwargs):
    """Check the provided kwargs against the default values for kwargs."""
    if kwargs is not None:
        for key, value in iteritems(default_kwargs):
            if key in kwargs:
                # Check the value type against the default.
                if value is None:
                    continue
                type_ = type(value)
                if not isinstance(kwargs[key], type_):
                    raise TypeError(
                        "'{0}' must be of type: {1}.".format(
                            key, str(type_)))
            else:
                # Set the kwarg as the default
                kwargs[key] = value
        if len(kwargs) != len(default_kwargs):
            warnings.warn("Unrecognized kwargs: {0}".format(
                str([key for key in kwargs if key not in default_kwargs])))
    else:
        # Set the kwargs as the defaults
        kwargs = default_kwargs

    return kwargs


def _mk_new_dictlist(ref_dictlist, old_dictlist, ensure_unique=False):
    """Return a new DictList with object references updated."""
    items = ref_dictlist.get_by_any([i.id if hasattr(i, "id") else str(i)
                                     for i in old_dictlist])
    if ensure_unique:
        items = set(items)
    return DictList(items)


def _get_matrix_constructor(matrix_type, dtype, matrix_type_default="dense",
                            dtype_default=np.float64):
    """Create a matrix constructor for the specified matrix type.

    Parameters
    ----------
    matrix_type: {'dense', 'dok', 'lil', 'DataFrame', 'symbolic'}, optional
        The desired type after for the matrix. If None, defaults to "dense".
    dtype: data-type, optional
        The desired array data-type for the stoichiometric matrix. If None,
        defaults to np.float64.

    Returns
    -------
    matrix: matrix of class 'dtype'
        The matrix for the MassModel returned as the given matrix_type
        and with a data-type of 'dtype'.

    Warnings
    --------
    This method is intended for internal use only. To safely create a
    matrix, use the appropriate MassModel method instead.

    """
    if matrix_type in _MATRIX_TYPES:
        pass
    elif matrix_type is None:
        matrix_type = matrix_type_default
    else:
        raise ValueError("Unrecognized matrix_type.")

    # Use the model's stored data-type if the data-type is not specified.
    if dtype is None:
        dtype = dtype_default

    # Dictionary of options for constructing the matrix
    matrix_constructor = dict(zip(_MATRIX_TYPES,
                                  [np.zeros, dok_matrix, lil_matrix,
                                   np.zeros, np.zeros]))
    constructor = matrix_constructor[matrix_type]
    return (constructor, matrix_type, dtype)


# Define small conversion functions based on the original matrix type.
def _to_dense(matrix):
    """Convert matrix to a numpy array."""
    if isinstance(matrix, np.ndarray):
        pass
    elif isinstance(matrix, pd.DataFrame):
        matrix = matrix.as_matrix()
    elif isinstance(matrix, sym.Matrix):
        matrix = np.array(matrix)
    else:
        matrix = matrix.toarray()

    return matrix


def _to_lil(matrix):
    """Convert matrix to a scipy lil matrix."""
    if isinstance(matrix, sym.Matrix):
        matrix = sym.matrix2numpy(matrix, dtype=float)
    return lil_matrix(matrix)


def _to_dok(matrix):
    """Convert matrix to a scipy dok matrix."""
    if isinstance(matrix, sym.Matrix):
        matrix = sym.matrix2numpy(matrix, dtype=float)
    return dok_matrix(matrix)


def _make_logger(name):
    """Make the logger instance and set the default format."""
    name = name.split(".")[-1]
    logging.basicConfig(format="%(name)s %(levelname)s: %(message)s")
    logger = logging.getLogger(name)
    return logger


__all__ = (
    "show_versions", "ensure_iterable", "ensure_non_negative_value",
    "convert_matrix")
