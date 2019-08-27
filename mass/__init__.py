# -*- coding: utf-8 -*-
import warnings as _warnings
from os import name as _name
from os.path import abspath as _abspath
from os.path import dirname as _dirname

from mass import (io, enzyme_modules)
from mass.core import (
    MassConfiguration, MassMetabolite, MassModel, MassReaction, Simulation,
    MassSolution, Unit, UnitDefinition)
from mass.util import show_versions, qcqa

__version__ = "0.1.0a41"

# set the warning format to be prettier and fit on one line
_mass_path = _dirname(_abspath(__file__))
if _name == "posix":
    _warning_base = "%s:%s \x1b[1;31m%s\x1b[0m: %s\n"  # colors
else:
    _warning_base = "%s:%s %s: %s\n"


def _warn_format(message, category, filename, lineno, file=None, line=None):
    """Set the warning format to be prettier and fit on one line."""
    shortname = filename.replace(_mass_path, "mass", 1)
    return _warning_base % (shortname, lineno, category.__name__, message)


_warnings.formatwarning = _warn_format

__all__ = ("_warn_format")
