# -*- coding: utf-8 -*-
"""Contains utility functions to assist in various :mod:`mass` functions."""
import logging
from operator import lt, le
import warnings

from cobra import DictList

from depinfo import print_dependencies

from six import integer_types, iteritems, string_types


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


def ensure_non_negative_value(value, exclude_zero=False):
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
    if value is not None:
        if not isinstance(value, (integer_types, float)):
            raise TypeError("Must be an int or float")
        if exclude_zero:
            comparision = le
            msg = "Must be a postive number"
        else:
            comparision = lt
            msg = "Must be a non-negative number"

        if comparision(value, 0.):
            raise ValueError(msg)

    return value


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


def _make_logger(name):
    """Make the logger instance and set the default format."""
    name = name.split(".")[-1]
    logging.basicConfig(format="%(name)s %(levelname)s: %(message)s")
    logger = logging.getLogger(name)
    return logger


__all__ = (
    "show_versions", "ensure_iterable", "ensure_non_negative_value",)
