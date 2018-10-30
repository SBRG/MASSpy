# -*- coding: utf-8 -*-
"""This module contains the exceptions specific to masspy."""
from __future__ import absolute_import


class MassSBMLError(Exception):
    """Class to reclassify any masspy errors specific to SBML issues."""

    pass


class MassSimulationError(Exception):
    """Class to reclassify any masspy errors specific to Simulation issues."""

    pass
