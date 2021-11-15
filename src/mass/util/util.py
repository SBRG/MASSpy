# -*- coding: utf-8 -*-
"""Contains utility functions to assist in various :mod:`mass` functions."""
import logging
import warnings
from copy import copy
from operator import le, lt

import numpy as np
from cobra import DictList
from depinfo import print_dependencies
from six import integer_types, iteritems, string_types


LOG_COLORS = {
    logging.CRITICAL: "\x1b[90m",
    logging.ERROR: "\x1b[91m",
    logging.WARNING: "\x1b[93m",
    logging.INFO: "\x1b[94m",
    logging.DEBUG: "\x1b[92m",
    -1: "\x1b[0m",
}
"""dict: Contains logger levels and corresponding color codes."""


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
    if not hasattr(item, "__iter__") or isinstance(item, string_types):
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

        if comparision(value, 0.0):
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
    try:
        all_public = [i for i in obj.__dict__ if not i.startswith("_")]
    except AttributeError:
        all_public = []
    all_public += [
        i
        for i in obj.__class__.__dict__
        if i not in all_public and not i.startswith("_")
    ]
    if not exclude_parent:
        parent_public = get_public_attributes_and_methods(
            obj.__class__.__base__(), exclude_parent=True
        )
        all_public += [i for i in parent_public if i not in all_public]

    return sorted(all_public, key=str.lower)


def apply_decimal_precision(value, decimal_precision):
    """Apply the decimal precision to the value by through rounding.

    Parameters
    ----------
    value : float
        The value to be rounded.
    decimal_precision : int
        The decimal place right of the decimal to round the value to.

    """
    if decimal_precision is not None and value is not None:
        if hasattr(value, "__iter__"):
            new = [round(v, decimal_precision) for v in value]
            if not isinstance(value, (list, np.ndarray)):
                try:
                    value = value.__class__(new)
                except Exception:
                    pass
            elif isinstance(value, np.ndarray):
                value = np.array(new)

        else:
            value = round(value, decimal_precision)

    return value


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
                if type_ == float:
                    type_ = (float, integer_types)
                if type_ == list:
                    type_ = (list, np.ndarray)
                if not isinstance(kwargs[key], type_):
                    raise TypeError(
                        "'{0}' must be of type: {1}.".format(key, str(type_))
                    )
            else:
                # Set the kwarg as the default
                kwargs[key] = value
        if len(kwargs) != len(default_kwargs):
            warnings.warn(
                "Unrecognized kwargs: {0}".format(
                    str([key for key in kwargs if key not in default_kwargs])
                )
            )
    else:
        # Set the kwargs as the defaults
        kwargs = default_kwargs

    return kwargs


def _mk_new_dictlist(ref_dictlist, old_dictlist, ensure_unique=False):
    """Return a new DictList with object references updated."""
    items = ref_dictlist.get_by_any(
        [i.id if hasattr(i, "id") else str(i) for i in old_dictlist]
    )
    if ensure_unique:
        items = set(items)
    return DictList(items)


class ColorFormatter(logging.Formatter):
    """Colored Formatter for logging output.

    Based on
    http://uran198.github.io/en/python/2016/07/12/colorful-python-logging.html

    """

    def format(self, record, *args, **kwargs):
        """Set logger format."""
        # Copy old record
        new_record = copy(record)
        # Ensure level in log color dict
        if new_record.levelno in LOG_COLORS:
            # Get the color
            color, reset = LOG_COLORS[new_record.levelno], LOG_COLORS[-1]
            # Set the levelname color
            new_record.levelname = "{color}{level}:{reset}".format(
                color=color, level=new_record.levelname, reset=reset
            )

            # Set the message color
            new_record.msg = "{color}{msg}{reset}".format(
                color=color, msg=new_record.msg, reset=reset
            )

        return super(ColorFormatter, self).format(new_record, *args, **kwargs)


def _make_logger(name):
    """Make the logger instance and set the default format."""
    # Create colored formatter
    formatter = ColorFormatter("%(levelname)s %(message)s")
    handler = logging.StreamHandler()
    handler.setFormatter(formatter)

    # Get logger
    logger = logging.getLogger(name)
    # Add handler and return
    logger.addHandler(handler)
    return logger


def _log_msg(logger, level, verbose, msg, *args):
    """Log a message in the logger and print if desired.

    Warnings
    --------
    This method is intended for internal use only.

    """
    logger.log(level, msg, *args)
    if verbose:
        print(msg % args)


__all__ = (
    "show_versions",
    "ensure_iterable",
    "ensure_non_negative_value",
    "LOG_COLORS",
    "ColorFormatter",
)
