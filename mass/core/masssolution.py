# -*- coding: utf-8 -*-
"""TODO Module Docstrings."""
from __future__ import absolute_import

from warnings import warn

import matplotlib.pyplot as plt

import pandas as pd

from scipy.interpolate import interp1d

from six import iteritems, iterkeys, string_types

from mass.core.massmodel import MassModel
from mass.core.visualization import plot_simulation, plot_tiled_phase_portrait
from mass.util.DictWithID import DictWithID
from mass.util.util import ensure_iterable

# Strings of valid solution types
_CONC_STR = "Conc"
_FLUX_STR = "Flux"


class MassSolution(DictWithID):
    """Container to store the solutions for the simulation of a MassModel.

    MassSolution containers are given an ID of the following form:
        MassSolution.id = MassSolution_{id_or_model}_{solution_type}

    All solution MassSolutuion can be accessed via attribute accessors
    in addition to accessing via standard dict class methods.

    The MassSolution class is essentially a subclass of a dict with additional
    attributes and properties.

    Parameters
    ----------
    id_or_model: str, MassModel
        A string identifier or a MassModel to associate with the solution. If
        a MassModel, the MassModel identifier will be used, and a reference to
        the model is stored.
    solution_type: str
        The type of solution being stored. Can be one of the following:
        {"concentration", "flux"}
    solution_dictionary: dict, optional
        A dict containing the solutions to store in the Solution object. If
        None provided then the Solution object will be initialized with no
        solutions. Solutions can be added or changed later using various dict
        methods (i.e. Solution.update(solution_dictionary)).
    time: array-like, optional
        An array-like object containing the time points used in calculating the
        solutions to be stored.
    interpolate: bool, optional
        A bool indicating whether solutions should be stored as interpolating
        functions. If True, solutions are stored as scipy.interpolate.interp1d
        objects. If False, solutions are stored as numpy arrays. Default value
        is False.

    Notes
    -----
    MassSolution objects can behave as booleans, with empty MassSolution
        objects returning as False, and those with solutions returning as True.

    """

    def __init__(self, id_or_model, data_dict=None, solution_type=None,
                 time=None, interpolate=False):
        """Initialize MassSolution."""
        if solution_type not in {_CONC_STR, _FLUX_STR}:
            raise ValueError(
                "'{0}' is not a valid solution type.".format(solution_type))

        if isinstance(id_or_model, MassModel):
            id_or_model = "{0}_{1}Sol".format(str(id_or_model), solution_type)
        if not isinstance(id_or_model, string_types):
            raise TypeError(
                "'id_or_model' must be a MassModel instance or a string")

        super(MassSolution, self).__init__(id=id_or_model, data_dict=data_dict)
        self.solution_type = solution_type
        self._simulation = None
        self._time = time
        self._interpolate = None
        self.interpolate = interpolate

    @property
    def simulation(self):
        """Return the Simulation associated with the MassSolution."""
        return getattr(self, "_simulation")

    @property
    def df(self):
        """Return the MassSolution object as a pandas.DataFrame."""
        sols = dict((k, v(self.time)) if self.interpolate
                    else (k, v) for k, v in iteritems(self))
        df = pd.DataFrame.from_dict(sols)
        df.index = pd.Series(self.time, name="Time")
        return df

    @property
    def time(self):
        """Return time points stored in the MassSolution."""
        return getattr(self, "_time", None)

    @time.setter
    def time(self, value):
        """Set the time points that are stored in the MassSolution.

        Parameters
        ----------
        value: array-like
            An array-like object containing the time points used in calculating
            the solutions to be stored.

        Notes
        -----
        If the MassSolution is stored as numerical arrays and not as
            interpolating functions, the numerical arrays of the solutions will
            be recomputed to correspond to the new time points.

        """
        value = ensure_iterable(value)
        if self.interpolate:
            setattr(self, "_time", value)
        else:
            self.interpolate = True
            setattr(self, "_time", value)
            self.interpolate = False

    @property
    def t(self):
        """Shorthand to return the time points stored in the MassSolution."""
        return self.time

    @t.setter
    def t(self, value):
        """Shorthand to set the time points stored in the MassSolution.

        Parameters
        ----------
        value: array-like
            An array-like object containing the time points used in calculating
            the solutions to be stored.

        Notes
        -----
        If the MassSolution is stored as numerical arrays and not as
            interpolating functions, the numerical arrays of the solutions will
            be recomputed to correspond to the new time points.

        """
        self.time = value

    @property
    def interpolate(self):
        """Return whether solutions are given as interpolating functions."""
        return getattr(self, "_interpolate", None)

    @interpolate.setter
    def interpolate(self, value):
        """Set whether solutions are stored as interpolating functions.

        Parameters
        ----------
        value: bool
            If True, solutions are stored in the MassSolution object as
            scipy interpolating functions. Otherwise store solutions as arrays.

        """
        if not isinstance(value, bool):
            raise ValueError("value must be a bool")
        if self.time is None:
            warn("No time points associated with MassSolution. Cannot convert "
                 "between numerical arrays and interpolating functions.")
            return

        if value != self.interpolate:
            for key, sol in iteritems(self):
                if value and not isinstance(sol, interp1d):
                    self[key] = interp1d(self.time, sol, kind='cubic',
                                         fill_value='extrapolate')
                if not value and isinstance(sol, interp1d):
                    self[key] = sol(self.time)

        setattr(self, "_interpolate", value)

    @property
    def preview_time_profile(self):
        """Generate a preview of the time profile for the solution.

        Notes
        -----
        Will clear and use the current axis (accessible via plt.gca()).

        """
        ax = plt.gca()
        options = {"plot_function": "semilogx",
                   "grid": ("major", "x"),
                   "title": "Time Profile for " + self.id,
                   "xlabel": "Time",
                   "ylabel": (self.solution_type + "es"
                              if self.solution_type[-1] == "x"
                              else self.solution_type + "s")}
        ax.cla()
        ax = plot_simulation(self, ax=ax, legend="right outside", **options)
        ax.get_figure().set_size_inches((6, 4))

    @property
    def preview_phase_portraits(self):
        """Generate a preview of the phase portraits for the solution.

        Notes
        -----
        Will clear and use the current axis and (accessible via plt.gca()).

        """
        ax = plt.gca()
        ax.cla()
        ax = plot_tiled_phase_portrait(self, ax=ax, empty_tiles=None,
                                       title="Phase portraits for " + self.id)
        ax.get_figure().set_size_inches((7, 7))

    def __getattribute__(self, name):
        """Override of default __getattr__."""
        if name in self:
            return self[name]
        else:
            return super(MassSolution, self).__getattribute__(name)

    def __dir__(self):
        """Override of default __dir__."""
        return list(iterkeys(self)) + super(MassSolution, self).__dir__()
