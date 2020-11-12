# -*- coding: utf-8 -*-
"""This module contains Exceptions specific to :mod:`mass` module."""


class MassSBMLError(Exception):
    """SBML error class."""


class MassSimulationError(Exception):
    """Simulation error class."""


class MassEnsembleError(Exception):
    """Simulation error class."""


__all__ = (
    "MassSBMLError",
    "MassSimulationError",
    "MassEnsembleError",
)
