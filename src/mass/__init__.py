# -*- coding: utf-8 -*-
import warnings as _warnings
from os import name as _name
from os.path import abspath as _abspath
from os.path import dirname as _dirname

from mass import enzyme_modules, io
from mass.core import (
    MassConfiguration,
    MassMetabolite,
    MassModel,
    MassReaction,
    MassSolution,
    Unit,
    UnitDefinition,
)
from mass.simulation import Simulation
from mass.util import qcqa_model, show_versions, strip_time


__version__ = "0.1.7"

# set the warning format to be prettier and fit on one line
_MASS_PATH = _dirname(_abspath(__file__))
if _name == "posix":
    _WARNING_BASE = "%s:%s \x1b[1;31m%s\x1b[0m: %s\n"  # colors
else:
    _WARNING_BASE = "%s:%s %s: %s\n"


def _warn_format(message, category, filename, lineno, file=None, line=None):
    """Set the warning format to be prettier and fit on one line."""
    shortname = filename.replace(_MASS_PATH, "mass", 1)
    return _WARNING_BASE % (shortname, lineno, category.__name__, message)


_warnings.formatwarning = _warn_format

__all__ = "_warn_format"
