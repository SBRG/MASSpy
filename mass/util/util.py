# -*- coding: utf-8 -*-
"""TODO Module Docstrings."""
from __future__ import absolute_import, print_function

import inspect
import logging
import re
import warnings

from depinfo import print_dependencies

import numpy as np

import pandas as pd

from scipy.sparse import dok_matrix, lil_matrix

from six import integer_types, string_types

import sympy as sym

from cobra import DictList

_MATRIX_TYPES = ["dense", "dok", "lil", "DataFrame", "symbolic"]


# Public
def show_versions():
    """Print dependency information."""
    print_dependencies("masspy")


def get_object_attributes(obj):
    """Get the public attributes for an object. Includes properties."""
    def filter_function(item):
        """Filter for non-routine and non-class methods."""
        return not inspect.isroutine(item) and not inspect.ismethod(item)

    # Sometimes raises logger warnings, therefore temporary disabled
    logging.disable(logging.WARNING)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        attributes = sorted([
            k for k, _ in inspect.getmembers(obj, filter_function)
            if not k.startswith("_")], key=str.lower)

    # Reenable logger warnings
    logging.disable(logging.NOTSET)
    return attributes


def get_subclass_specific_attributes(obj):
    """Get the public attributes for an object. Includes properties."""
    attributes = get_object_attributes(obj)
    attributes = [
        attr for attr in attributes
        if attr not in get_object_attributes(obj.__class__.__base__())]

    return attributes


def ensure_iterable(list_to_check):
    """Ensure the given list is an iterable.

    Parameters
    ----------
    list_to_check: list
        The list to ensure is iterable.

    """
    # Make list iterable if necessary
    if list_to_check is None:
        list_to_check = list()
    if not hasattr(list_to_check, "__iter__") or \
       isinstance(list_to_check, string_types):
        list_to_check = [list_to_check]

    list_to_check = list(list_to_check)
    return list_to_check


def ensure_non_negative_value(value):
    """Ensure provided value is a non-negative value, or None.

    Will raise a ValueError if the provided value is negative.

    Parameters
    ----------
    value: float
        The value to ensure is non-negative

    """
    if value is None:
        pass
    elif not isinstance(value, (integer_types, float)):
        raise TypeError("Must be an int or float")
    elif value < 0.:
        raise ValueError("Must be a non-negative number")


def convert_matrix(matrix, matrix_type, dtype, row_ids=None, col_ids=None):
    """Convert a matrix to a different type.

    Parameters
    ----------
    matrix: array or array-like
        The matrix to convert.
    matrix_type: {'dense', 'dok', 'lil', 'DataFrame', 'symbolic'}
        The desired type after converting the matrix.
    dtype: data-type
        The desired data-type for the array.
    row_ids: array
        The idenfifiers for each row. Only necessary if type is "dataframe",
        otherwise is ignored.
    col_ids: array
        The idenfifiers for each column. Only necessary if type is "dataframe",
        otherwise is ignored.

    Warnings
    --------
    This method is not the safest way to convert a matrix. To safely convert
    a matrix into another type, use the "matrix_type" argument in the method
    that returns the desired matrix.

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


# Internal
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
    if isinstance(matrix, sym.Matrix):
        matrix = sym.matrix2numpy(matrix, dtype=float)
    return lil_matrix(matrix)


def _to_dok(matrix):
    if isinstance(matrix, sym.Matrix):
        matrix = sym.matrix2numpy(matrix, dtype=float)
    return dok_matrix(matrix)
