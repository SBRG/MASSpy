# -*- coding: utf-8 -*-
"""TODO Module Docstrings."""
from __future__ import absolute_import

from warnings import warn

import libsbml

from six import integer_types, iteritems, string_types

from tabulate import tabulate

from cobra.core.object import Object

from mass.util.DictWithID import DictWithID

# Units
_SI_PREFIXES_DICT = DictWithID(
    id="SI Unit Scale Prefixes", dictionary={
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
    }
)

_SBML_BASE_UNIT_KINDS_DICT = DictWithID(
    id="SBML Base Unit Kinds", dictionary={
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
    }
)


class Unit(object):
    """Manage units via this implementation of the SBML Unit specifications.

    Parameters
    ----------
    kind: str
        A string representing the SBML L3 recognized base unit.
    exponent: int
        The exponent on the Unit.
    scale: int, str
        An integer representing the scale of the unit, or string for one of the
        pre-defined scales.
    multiplier: float, int
        A number used to multiply the unit by a real-numbered factor, enabling
        units that are not necessarily a power-of-ten multiple.

    """

    def __init__(self, kind, exponent, scale, multiplier):
        """Initialize the Unit Object."""
        object.__init__(self)
        # Set Unit kind
        self._kind = None
        self.kind = kind

        # Set Unit exponent
        self._exponent = None
        self.exponent = exponent

        # Set Unit scale
        self._scale = None
        self.scale = scale

        # Set Unit multiplier
        self._multiplier = None
        self.multiplier = multiplier

    @property
    def kind(self):
        """Return the unit kind of the Unit."""
        return list(getattr(self, "_kind"))[0]

    @kind.setter
    def kind(self, kind):
        """Set the unit kind of the Unit.

        Parameters
        ----------
        kind: str
            An SBML recognized unit kind identifier as a string.

        """
        if kind not in _SBML_BASE_UNIT_KINDS_DICT:
            raise ValueError(
                "Invalid SBML Base Unit Kind '{0}'. Allowable values can be "
                "viewed by passing the string 'BaseUnitKinds' to the function "
                "'print_defined_unit_values' from the mass.core.units "
                "submodule.".format(kind))

        setattr(self, "_kind", {kind: _SBML_BASE_UNIT_KINDS_DICT[kind]})

    @property
    def exponent(self):
        """Return the exponent of the Unit."""
        return getattr(self, "_exponent")

    @exponent.setter
    def exponent(self, exponent):
        """Set the exponent of the Unit.

        Parameters
        ----------
        exponent: int
            An integer representing the exponent of the Unit.

        """
        if isinstance(exponent, (integer_types, float))\
           and float(exponent).is_integer():
            setattr(self, "_exponent", int(exponent))
        else:
            raise TypeError("exponent must be an integer")

    @property
    def scale(self):
        """Return the scale of the Unit."""
        return getattr(self, "_scale")

    @scale.setter
    def scale(self, scale):
        """Set the scale of the Unit.

        Parameters
        ----------
        scale: int, str
            An integer representing the scale of the Unit, or a string
            from the predefined SI prefixes. Not case sensitive.

        """
        if isinstance(scale, string_types):
            if scale.lower() not in _SI_PREFIXES_DICT:
                raise ValueError(
                    "Invalid SI Scale Prefix '{0}'. Allowable values can be "
                    "viewed by passing the string 'Scales' to the function "
                    "'print_defined_unit_values' from the mass.core.units "
                    "submodule.".format(scale))
            scale = _SI_PREFIXES_DICT[scale.lower()]

        if isinstance(scale, (integer_types, float))\
           and float(scale).is_integer():
            setattr(self, "_scale", int(scale))
        else:
            raise TypeError("scale must be an integer")

    @property
    def multiplier(self):
        """Return the multiplier of the Unit."""
        return getattr(self, "_multiplier")

    @multiplier.setter
    def multiplier(self, multiplier):
        """Set the multiplier of the Unit.

        Parameters
        ----------
        multiplier: float
            A numerical value representing the multiplier of the Unit.

        """
        if not isinstance(multiplier, (integer_types, float)):
            raise TypeError("multiplier must be an int.")

        self._multiplier = multiplier

    def __str__(self):
        """Override of default str() implementation."""
        return "kind: %s; exponent: %s; scale: %s; multiplier: %s" % (
            self.kind, self.exponent, self.scale, self.multiplier)

    def __repr__(self):
        """Override of default repr() implementation."""
        return "<%s at 0x%x %s>" % (
            self.__class__.__name__, id(self), str(self))


_PREDEFINED_UNITS_DICT = DictWithID(
    id="Pre-defined Units", dictionary={
        "mole": Unit(
            kind="mole", exponent=1, scale=0, multiplier=1),
        "millimole": Unit(
            kind="mole", exponent=1, scale=-3, multiplier=1),
        "litre": Unit(
            kind="litre", exponent=1, scale=0, multiplier=1),
        "per_litre": Unit(
            kind="litre", exponent=-1, scale=0, multiplier=1),
        "second": Unit(
            kind="second", exponent=1, scale=0, multiplier=1),
        "per_second": Unit(
            kind="second", exponent=-1, scale=0, multiplier=1),
        "hour": Unit(
            kind="second", exponent=1, scale=0, multiplier=3600),
        "per_hour": Unit(
            kind="second", exponent=-1, scale=0, multiplier=3600),
    }
)


class UnitDefinition(Object):
    """Manage units via implementation of SBML UnitDefinition specifications.

    Parameters
    ----------
    id: str
        The identifier to associate with this unit definition
    name : str, optional
        A human readable name for this unit definition.

    Attributes
    ----------
    list_of_units: list, optional
        A list iterable containing mass.Unit objects that are needed to
        define the UnitDefinition, or a string that corresponds with the
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
        """Create a Unit object ahd add it to the UnitDefinition.

        Parameters
        ----------
        kind: str
            A string representing the SBML L3 recognized base unit.
        exponent: int
            The exponent on the Unit. Default is 1.
        scale: int, str
            An integer representing the scale of the unit, or string for one of
            the pre-defined scales. Default is 0.
        multiplier: float, int
            A number used to multiply the unit by a real-numbered factor,
            enabling units that are not necessarily a power-of-ten multiple.
            Default is 1.

        Returns
        -------
        unit: mass.Unit
            The newly created mass.Unit object.

        """
        unit = Unit(
            kind=kind, exponent=exponent, scale=scale, multiplier=multiplier)
        self.add_units([unit])

        return unit

    def add_units(self, new_units):
        """Add Unit objects to the UnitDefinition's list of units.

        Parameters
        ----------
        new_units: list
            A list of mass.Units objects and/or the string identifiers of
            pre-defined units to add to the list_of_units attribute of the
            UnitDefinition object.

        """
        # Get set of units to add
        to_add = _units_to_alter(new_units)
        # Get current units in UnitDefinition
        current = set(self.list_of_units)
        # Remove units from UnitDefinition
        current.update(to_add)
        # Update attribute
        self.list_of_units = list(current)

    def remove_units(self, units_to_remove):
        """Remove Unit objects from the UnitDefinition's list of units.

        Parameters
        ----------
        units_to_remove: list
            A list of mass.Units objects and/or the string identifiers of
            pre-defined units to remove to the list_of_units attribute of the
            UnitDefinition object.

        """
        # Get set of units to remove
        to_remove = _units_to_alter(units_to_remove)
        # Get current units in UnitDefinition
        current = set(self.list_of_units)
        # Remove units from UnitDefinition
        current.difference_update(to_remove)
        # Update attribute
        self.list_of_units = list(current)

    def __repr__(self):
        """Override of default repr() implementation."""
        if self.name:
            name = self.name
        else:
            name = self.id
        return "<%s %s at 0x%x>" % (self.__class__.__name__, name, id(self))


def _units_to_alter(units):
    """Create a set of units to alter in the unit definition.

    Warnings
    --------
    This method is intended for internal use only.

    """
    # Ensure list input
    if isinstance(units, (string_types, Unit)):
        warn("needs to be pass in a list.")
        units = [units]
    if not isinstance(units, list):
        raise TypeError("must be a list.")

    # Create a set of units to be altered and return the set.
    to_alter = set()
    for unit in units:
        # Add a predefined unit to the UnitDefinition if a string.
        if isinstance(unit, str) and unit in _PREDEFINED_UNITS_DICT:
            to_alter.add(_PREDEFINED_UNITS_DICT[unit])
        # Add user-defined unit to the UnitDefinition if a Unit object.
        elif isinstance(unit, Unit):
            to_alter.add(unit)
        # Otherwise raise an error and skip.
        else:
            warn("Skipping unrecognized unit: '{0}'.".format(str(unit)))
            continue

    return to_alter


def print_defined_unit_values(value="Units"):
    """Print the pre-defined unit quantities in the mass.core.units submodule.

    Parameters
    ----------
    value: {"Scales", "BaseUnitKinds", "Units", "all"}
        A string representing which pre-defined values to display.
        Default is "Units" to display all pre-defined Unit objects.

    """
    predefined_items = ["Units", "Scales", "BaseUnitKinds", "all"]
    # Check input
    if value not in predefined_items:
        raise ValueError(
            "'{0}' not recognized, must be one of the following '{1}'.".format(
                str(value), str(predefined_items)))

    # Make a list of the pre-defined items to print.
    if value.lower() == "all":
        predefined_items.remove(value)
    else:
        predefined_items = [value]

    # Print relevant items as tables
    for item in predefined_items:
        # Get dictionary to display
        value_dict, headers = {
            "Units": (_PREDEFINED_UNITS_DICT, ["Unit", "Definition"]),
            "Scales": (_SI_PREFIXES_DICT, ["Prefix", "Scale Value"]),
            "BaseUnitKinds": (
                _SBML_BASE_UNIT_KINDS_DICT, ["Base Unit", "SBML Value"])
        }.get(item)
        # Format dictionary items into table
        content = [[k, str(v)] for k, v in iteritems(value_dict)]
        # Add to final tables
        table = [tabulate(content, headers=headers, tablefmt="Simple")]
        print(tabulate(
            [table], headers=[value_dict.id], tablefmt="fancy_grid"))
