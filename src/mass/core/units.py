# -*- coding: utf-8 -*-
r"""Unit and UnitDefinition implementation based on SBML specifications.

The :mod:`~.units` module contains the :class:`Unit` and
:class:`UnitDefinition` classes based on the implementation of units in
`SBML <https://sbml.org/>`_.

Note that :mod:`mass` does not support automatic unit tracking to ensure
unit consistency. Therefore, it is incumbent upon the user to maintain unit
consistency as they use the various :mod:`mass` modules and functions.

To view valid units, use the :func:`print_defined_unit_values` function.
Please send a PR if you want to add something to the pre-built
:class:`Unit` \s.
"""
from warnings import warn

import libsbml
from cobra.core.object import Object
from six import integer_types, iteritems, iterkeys, itervalues, string_types
from tabulate import tabulate

from mass.util.dict_with_id import DictWithID


SBML_BASE_UNIT_KINDS_DICT = DictWithID(
    id="SBML Base Unit Kinds",
    data_dict={
        "ampere": libsbml.UNIT_KIND_AMPERE,
        "avogadro": libsbml.UNIT_KIND_AVOGADRO,
        "becquerel": libsbml.UNIT_KIND_BECQUEREL,
        "candela": libsbml.UNIT_KIND_CANDELA,
        "coulomb": libsbml.UNIT_KIND_COULOMB,
        "dimensionless": libsbml.UNIT_KIND_DIMENSIONLESS,
        "farad": libsbml.UNIT_KIND_FARAD,
        "gram": libsbml.UNIT_KIND_GRAM,
        "gray": libsbml.UNIT_KIND_GRAY,
        "henry": libsbml.UNIT_KIND_HENRY,
        "hertz": libsbml.UNIT_KIND_HERTZ,
        "item": libsbml.UNIT_KIND_ITEM,
        "joule": libsbml.UNIT_KIND_JOULE,
        "katal": libsbml.UNIT_KIND_KATAL,
        "kelvin": libsbml.UNIT_KIND_KELVIN,
        "kilogram": libsbml.UNIT_KIND_KILOGRAM,
        "liter": libsbml.UNIT_KIND_LITER,
        "litre": libsbml.UNIT_KIND_LITRE,
        "lumen": libsbml.UNIT_KIND_LUMEN,
        "lux": libsbml.UNIT_KIND_LUX,
        "meter": libsbml.UNIT_KIND_METER,
        "metre": libsbml.UNIT_KIND_METRE,
        "mole": libsbml.UNIT_KIND_MOLE,
        "newton": libsbml.UNIT_KIND_NEWTON,
        "ohm": libsbml.UNIT_KIND_OHM,
        "pascal": libsbml.UNIT_KIND_PASCAL,
        "radian": libsbml.UNIT_KIND_RADIAN,
        "second": libsbml.UNIT_KIND_SECOND,
        "siemens": libsbml.UNIT_KIND_SIEMENS,
        "sievert": libsbml.UNIT_KIND_SIEVERT,
        "steradian": libsbml.UNIT_KIND_STERADIAN,
        "tesla": libsbml.UNIT_KIND_TESLA,
        "volt": libsbml.UNIT_KIND_VOLT,
        "watt": libsbml.UNIT_KIND_WATT,
        "weber": libsbml.UNIT_KIND_WEBER,
        "invalid": libsbml.UNIT_KIND_INVALID,
    },
)
""":class:`~.DictWithID`: Contains SBML base units and their ``int`` values."""

SI_PREFIXES_DICT = DictWithID(
    id="SI Unit Scale Prefixes",
    data_dict={
        "atto": -18,
        "femto": -15,
        "pico": -12,
        "nano": -9,
        "micro": -6,
        "milli": -3,
        "centi": -2,
        "deci": -1,
        "deca": 1,
        "hecto": 2,
        "kilo": 3,
        "mega": 6,
        "giga": 9,
        "tera": 12,
        "peta": 15,
        "exa": 18,
    },
)
""":class:`~.DictWithID`: Contains SI unit prefixes and scale values."""


class Unit:
    """Manage units via this implementation of the SBML Unit specifications.

    Parameters
    ----------
    kind : str or int
        A string representing the SBML Level 3 recognized base unit or its
        corresponding SBML integer value as defined in
        :const:`SBML_BASE_UNIT_KINDS_DICT`.
    exponent : int
        The unit exponent.
    scale : int or str
        An integer representing the scale of the unit, or a string for one of
        the pre-defined SI scales in :const:`SI_PREFIXES_DICT`.
    multiplier : float
        A number used to multiply the unit by a real-numbered factor, enabling
        units that are not necessarily a power-of-ten multiple.

    """

    def __init__(self, kind, exponent, scale, multiplier):
        """Initialize the Unit."""
        object.__init__(self)
        # Initialize Unit kind, exponent, scale, and multiplier
        self._kind = None
        self._exponent = None
        self._scale = None
        self._multiplier = None

        for name, value in zip(
            ["kind", "exponent", "scale", "multiplier"],
            [kind, exponent, scale, multiplier],
        ):
            setattr(self, name, value)

    @property
    def kind(self):
        """Return the unit kind of the :class:`Unit`.

        Parameters
        ----------
        kind : str
            An SBML recognized unit kind as a string.

        """
        return getattr(self, "_kind")

    @kind.setter
    def kind(self, kind):
        """Set the unit kind of the :class:`Unit`."""
        # Ensure input is SBML compliant, remove invalid kinds.
        valid_keys = list(iterkeys(SBML_BASE_UNIT_KINDS_DICT))
        valid_keys.remove("invalid")

        valid_values = list(itervalues(SBML_BASE_UNIT_KINDS_DICT))
        valid_values.remove(36)

        if isinstance(kind, string_types) and kind in valid_keys:
            pass
        elif isinstance(kind, integer_types) and kind in valid_values:
            # Get corresponding base unit kind string for the value
            kind = valid_keys[valid_values.index(kind)]
        else:
            # Raise an error if not SBML compliant
            raise ValueError(
                "Invalid SBML Base Unit Kind '{0}'. Allowable values can be "
                "viewed by passing the string 'BaseUnitKinds' to the function "
                "'print_defined_unit_values' from the mass.core.units "
                "submodule.".format(kind)
            )
        setattr(self, "_kind", kind)

    @property
    def exponent(self):
        """Return the exponent of the :class:`Unit`.

        Parameters
        ----------
        exponent : int
            The exponent of the unit as an integer.

        """
        return getattr(self, "_exponent")

    @exponent.setter
    def exponent(self, exponent):
        """Set the exponent of the :class:`Unit`."""
        if (
            isinstance(exponent, (integer_types, float))
            and float(exponent).is_integer()
        ):
            setattr(self, "_exponent", int(exponent))
        else:
            raise TypeError("exponent must be an integer")

    @property
    def scale(self):
        """Return the scale of the :class:`Unit`.

        Parameters
        ----------
        scale : int or str
            An integer representing the scale of the unit, or a
            string from the pre-defined SI prefixes. Not case sensitive.

        """
        return getattr(self, "_scale")

    @scale.setter
    def scale(self, scale):
        """Set the scale of the :class:`Unit`."""
        if isinstance(scale, string_types):
            if scale.lower() not in SI_PREFIXES_DICT:
                raise ValueError(
                    "Invalid SI Scale Prefix '{0}'. Allowable values can be "
                    "viewed by passing the string 'Scales' to the function "
                    "'print_defined_unit_values' from the mass.core.units "
                    "submodule.".format(scale)
                )
            scale = SI_PREFIXES_DICT[scale.lower()]

        if isinstance(scale, (integer_types, float)) and float(scale).is_integer():
            setattr(self, "_scale", int(scale))
        else:
            raise TypeError("scale must be an integer")

    @property
    def multiplier(self):
        """Get or set the multiplier of the :class:`Unit`.

        Parameters
        ----------
        multiplier : float
            A numerical value representing a multiplier for the unit.

        """
        return getattr(self, "_multiplier")

    @multiplier.setter
    def multiplier(self, multiplier):
        """Set the multiplier of the :class:`Unit`."""
        if not isinstance(multiplier, (integer_types, float)):
            raise TypeError("multiplier must be an int.")

        self._multiplier = multiplier

    def __str__(self):
        """Override of default :class:`str` implementation.

        Warnings
        --------
        This method is intended for internal use only.

        """
        return "kind: %s; exponent: %s; scale: %s; multiplier: %s" % (
            self.kind,
            self.exponent,
            self.scale,
            self.multiplier,
        )

    def __repr__(self):
        """Override of default :func:`repr` implementation.

        Warnings
        --------
        This method is intended for internal use only.

        """
        return "<%s at 0x%x %s>" % (self.__class__.__name__, id(self), str(self))


PREDEFINED_UNITS_DICT = DictWithID(
    id="Pre-defined Units",
    data_dict={
        "mole": Unit(kind="mole", exponent=1, scale=0, multiplier=1),
        "millimole": Unit(kind="mole", exponent=1, scale=-3, multiplier=1),
        "litre": Unit(kind="litre", exponent=1, scale=0, multiplier=1),
        "per_litre": Unit(kind="litre", exponent=-1, scale=0, multiplier=1),
        "second": Unit(kind="second", exponent=1, scale=0, multiplier=1),
        "per_second": Unit(kind="second", exponent=-1, scale=0, multiplier=1),
        "hour": Unit(kind="second", exponent=1, scale=0, multiplier=3600),
        "per_hour": Unit(kind="second", exponent=-1, scale=0, multiplier=3600),
        "per_gDW": Unit(kind="gram", exponent=-1, scale=0, multiplier=1),
    },
)
r""":class:`~.DictWithID`: Contains pre-built :class:`Unit`\ s."""


class UnitDefinition(Object):
    r"""Manage units via implementation of SBML UnitDefinition specifications.

    Parameters
    ----------
    id : str
        The identifier to associate with the unit definition
    name : str
        A human readable name for the unit definition.

    Attributes
    ----------
    list_of_units : list
        A list containing :class:`Unit`\ s that are needed to define the
        :class:`UnitDefinition`, or a string that corresponds with the
        pre-defined units. Invalid units are ignored.

    """

    def __init__(self, id=None, name="", list_of_units=None):
        """Initialize a UnitDefinition object."""
        super(UnitDefinition, self).__init__(id=id, name=name)
        if list_of_units is None:
            self.list_of_units = []
        elif not isinstance(list_of_units, list):
            raise TypeError("list_of_units must be a list.")
        else:
            self.list_of_units = []
            self.add_units(list_of_units)

    def create_unit(self, kind, exponent=1, scale=0, multiplier=1):
        """Create a :class:`Unit` and add it to the :class:`UnitDefinition`.

        Parameters
        ----------
        kind : str
            A string representing the SBML Level 3 recognized base unit.
        exponent : int
            The exponent on the unit. Default is ``1``.
        scale : int or str
            An integer representing the scale of the unit, or
            a string for one of the pre-defined scales. Default is ``0.``
        multiplier : float
            A number used to multiply the unit by a real-numbered factor,
            enabling units that are not necessarily a power-of-ten multiple.
            Default is ``1.``

        """
        unit = Unit(kind=kind, exponent=exponent, scale=scale, multiplier=multiplier)
        self.add_units([unit])

    def add_units(self, new_units):
        r"""Add :class:`Unit`\ s to the :attr:`list_of_units`.

        Parameters
        ----------
        new_units : list
            A list of :class:`Unit`\ s and the string identifiers of pre-built
            units to add to the :attr:`list_of_units`

        """
        # Get set of units to add
        to_add = self._units_to_alter(new_units)
        # Get current units in UnitDefinition
        current = set(self.list_of_units)
        # Remove units from UnitDefinition
        current.update(to_add)
        # Update attribute
        self.list_of_units = list(current)

    def remove_units(self, units_to_remove):
        r"""Remove :class:`Unit`\ s from the :attr:`list_of_units`.

        Parameters
        ----------
        units_to_remove : list
            A list of :class:`Unit`\ s and/or the string corresponding to the
            unit :attr:`.Unit.kind` to remove from the :attr:`list_of_units`.

        """
        # Get set of units to remove
        to_remove = self._units_to_alter(units_to_remove)
        # Get current units in UnitDefinition
        current = set(self.list_of_units)
        # Remove units from UnitDefinition
        current.difference_update(to_remove)
        # Update attribute
        self.list_of_units = list(current)

    def _units_to_alter(self, units):
        """Create a set of units to alter in the unit definition.

        Warnings
        --------
        This method is intended for internal use only.

        """
        # Ensure list input
        if isinstance(units, (string_types, Unit)):
            warn("needs to be passed as a list.")
            units = [units]
        if not isinstance(units, list):
            raise TypeError("must be a list.")

        list_of_units = self.list_of_units
        u_kinds = [u.kind for u in list_of_units]

        # Create a set of units to be altered and return the set.
        to_alter = set()
        for unit in units:
            # Add a predefined unit to alter if a string.
            if isinstance(unit, str) and unit in PREDEFINED_UNITS_DICT:
                to_alter.add(PREDEFINED_UNITS_DICT[unit])
            elif isinstance(unit, str) and unit in u_kinds:
                unit = [u for u in list_of_units if u.kind == unit].pop()
                to_alter.add(unit)
            # Add user-defined unit to alter if a Unit object.
            elif isinstance(unit, Unit):
                to_alter.add(unit)
            # Otherwise raise an error and skip.
            else:
                warn("Skipping unrecognized unit: '{0}'.".format(str(unit)))
                continue

        return to_alter

    def __repr__(self):
        """Override of default :func:`repr` implementation.

        Warnings
        --------
        This method is intended for internal use only.

        """
        if self.name:
            name = self.name + ' "' + self.id + '"'
        else:
            name = self.id
        return "<%s %s at 0x%x>" % (self.__class__.__name__, name, id(self))

    def __iter__(self):
        """Override of default :func:`iter` implementation.

        Warnings
        --------
        This method is intended for internal use only.

        """
        for unit in self.list_of_units:
            yield unit


def print_defined_unit_values(value="Units"):
    r"""Print the pre-defined unit quantities in the :mod:`.units` submodule.

    Parameters
    ----------
    value: str
        A string representing which pre-defined values to display.
        Must be one of the following:

            * ``"Scales"``
            * ``"BaseUnitKinds"``
            * ``"Units"``
            * ``"all"``

        Default is ``"Units"`` to display all pre-defined :class:`Unit`\ s.

    """
    predefined_items = ["Units", "Scales", "BaseUnitKinds", "all"]
    # Check input
    if value not in predefined_items:
        raise ValueError(
            "'{0}' not recognized, must be one of the following '{1}'.".format(
                str(value), str(predefined_items)
            )
        )

    # Make a list of the pre-defined items to print.
    if value.lower() == "all":
        predefined_items.remove(value)
    else:
        predefined_items = [value]

    # Print relevant items as tables
    for item in predefined_items:
        # Get dictionary to display
        value_dict, headers = {
            "Units": (PREDEFINED_UNITS_DICT, ["Unit", "Definition"]),
            "Scales": (SI_PREFIXES_DICT, ["Prefix", "Scale Value"]),
            "BaseUnitKinds": (SBML_BASE_UNIT_KINDS_DICT, ["Base Unit", "SBML Value"]),
        }.get(item)
        # Format dictionary items into table
        content = [[k, str(v)] for k, v in iteritems(value_dict)]
        # Add to final tables
        table = [tabulate(content, headers=headers, tablefmt="Simple")]
        print(tabulate([table], headers=[value_dict.id], tablefmt="fancy_grid"))


__all__ = (
    "SBML_BASE_UNIT_KINDS_DICT",
    "SI_PREFIXES_DICT",
    "Unit",
    "PREDEFINED_UNITS_DICT",
    "UnitDefinition",
    "print_defined_unit_values",
)
