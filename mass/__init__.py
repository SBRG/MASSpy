# -*- coding: utf-8 -*-

from __future__ import absolute_import

# set the warning format to be on a single line
import warnings as _warnings
from os import name as _name
from os.path import abspath as _abspath
from os.path import dirname as _dirname

from mass.core import (
	MassMetabolite, MassModel, MassReaction, simulate, find_steady_state,
	plot_simulation, plot_phase_portrait, plot_tiled_phase_portrait,
	get_plot_defaults, get_tiled_defaults)

from mass.analysis import (
	pools_from_string, net_fluxes_from_strings)
from mass import io
from mass import util
from mass.util.version_info import show_versions

__version__ = "0.1.0a4"

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
