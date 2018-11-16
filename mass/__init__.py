# -*- coding: utf-8 -*-

from __future__ import absolute_import

# set the warning format to be on a single line
import warnings as _warnings
from os import name as _name
from os.path import abspath as _abspath
from os.path import dirname as _dirname

from mass import (analysis, io)
from mass.core import (
    MassMetabolite, MassModel, MassReaction, Simulation, get_defaults,
    plot_phase_portrait, plot_simulation, plot_tiled_phase_portrait,
    set_defaults)

from mass.util.qcqa import qcqa_model, qcqa_simulation
from mass.util.util import show_versions, strip_time

__version__ = "0.1.0a25"

# set the warning format to be prettier and fit on one line
_mass_path = _dirname(_abspath(__file__))
if _name == "posix":
    _warning_base = "%s:%s \x1b[1;31m%s\x1b[0m: %s\n"  # colors
else:
    _warning_base = "%s:%s %s: %s\n"


def _warn_format(message, category, filename, lineno, file=None, line=None):
    shortname = filename.replace(_mass_path, "mass", 1)
    return _warning_base % (shortname, lineno, category.__name__, message)

_warnings.formatwarning = _warn_format
