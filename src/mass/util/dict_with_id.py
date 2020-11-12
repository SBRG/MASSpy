# -*- coding: utf-8 -*-
"""
DictWithID and OrderedDictWithID are dictionaries with identifers attributes.

The :mod:`~.dict_with_id` submodule utilizes the built-in :class:`type` to
dynamically create the :class:`DictWithID` and :class:`OrderedDictWithID`
classes based on whether the parent class should be a standard :class:`dict`
or an :class:`~.collections.OrderedDict`, therefore allowing each object to
inherit its parent class methods and generable behavior.

The :class:`DictWithID` and :class:`OrderedDictWithID` classes are primarily
used for in order to use dictionaries with the speciailized
:class:`~cobra.core.dictlist.DictList` container for both performance gains
and user convenience.
"""
from collections import OrderedDict
from copy import copy

from six import string_types


def __constructor(self, id=None, data_dict=None):
    """Initialize a DictWithID instance.

    Parameters
    ----------
    id: str, None
        The identifier to associate with the object.
    data_dict: dict, optional
        If provided, the new DictWithID will contain the key:value pairs from
        "data_dict". Otherwise the new DictWithID is initialized as empty.

    Warnings
    --------
    This object is primarily intended for internal use only.

    """
    super(type(self), self).__setattr__("_id", id)
    if isinstance(data_dict, dict):
        self.update(data_dict)
    elif data_dict is None:
        pass
    else:
        raise TypeError("data_dict must be a dict")


@property
def __id(self):
    """Return identifier of the dictionary."""
    return getattr(self, "_id", None)


@__id.setter
def __id(self, value):
    """Set the identifier of the dictionary."""
    if value == self.id:
        pass
    elif not isinstance(value, string_types):
        raise TypeError("ID must be a string")
    elif getattr(self, "_model", None) is not None:
        self._set_id_with_model(value)
    else:
        self._id = value


def __set_id_with_model(self, value):
    """Set the id of the dictionary to the associated MassModel.

    Warnings
    --------
    This method is intended for internal use only.

    """
    self._id = value


def __copy(self):
    """Copy a DictWithID object."""
    return copy(self)


def __copy__(self):
    """Override default copy() implementation."""
    return copy(super(type(self), self))


def __repr__(self):
    """Override of default repr() implementation."""
    return "<%s %s at 0x%x>" % (self.__class__.__name__, self.id, id(self))


def __str__(self):
    """Override of default str() implementation."""
    return dict(self).__str__()


def __doc__(self):
    """Return the docstring for the class."""
    cls_name = self.__class__.__name__
    supercls_str, all_subclasses = {
        "DictWithID": ["dict", ["MassSolution"]],
        "OrderedDictWithID": ["OrderedDict", ["EnzymeDict"]],
    }[cls_name]

    return """A specialized {1} with an additional ID attribute.

    The {0} class is essentially a subclass of an {1} with a string identifier.
    The {2} objects are all subclasses of {0}.

    Parameters
    ----------
    id: str, None
        The identifier to associate with the {0}.
    data_dict: dict, optional
        If provided, the new {0} will contain the key:value pairs from
        "data_dict". Otherwise the new {0} is initialized as empty.

    """.format(
        cls_name, supercls_str, str(tuple(all_subclasses))
    )


def __make_class_constructor(ordered=False):
    """Create a class constructor for a specialized dictionary class.

    Parameters
    ----------
    ordered: bool
        If True, the constructor will make a subclass of an OrderedDict.
        Otherwise,  the constructor will make a subclass of a regular dict.

    Returns:
    --------
    constructor: class constructor
        The class constructor for OrderedDictWithID or DictWithID

    Warnings
    --------
    This object is primarily intended for internal use only.

    """
    if ordered:
        cls_name, obj = ("OrderedDictWithID", OrderedDict)
    else:
        cls_name, obj = ("DictWithID", dict)

    constructor = type(
        cls_name,
        (obj,),
        {
            "__init__": __constructor,
            "id": __id,
            "copy": __copy,
            "_set_id_with_model": __set_id_with_model,
            "__copy__": __copy__,
            "__repr__": __repr__,
            "__str__": __str__,
            "__doc__": __doc__,
        },
    )

    return constructor


DictWithID = __make_class_constructor(False)
"""
class: Has ID attribute, inherits :class:`dict` methods"""

OrderedDictWithID = __make_class_constructor(True)
"""
class: Has ID attribute, inherits :class:`~collections.OrderedDict` methods."""

__all__ = ("DictWithID", "OrderedDictWithID")
