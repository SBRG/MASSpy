# -*- coding: utf-8 -*-

# Compatibility with Python 2.7
from __future__ import absolute_import

# Import necesary packages
import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.patches as patches

from pandas.compat import range, lrange, zip
from pandas.io.formats.printing import pprint_thing
from pandas.plotting._style import _get_standard_colors
from pandas.plotting._tools import _subplots, _set_ticks_props

from scipy.interpolate import interp1d
from math import inf
from six import iterkeys, itervalues
from cycler import cycler

# from cobra
from cobra import DictList

# from mass
from mass.core.massmetabolite import MassMetabolite
from mass.core.massreaction import MassReaction
from mass.core.massmodel import MassModel