# -*- coding: utf-8 -*-
r"""
MassSolution is a class for storing the simulation results.

After a :class:`~.Simulation` is used to simulate a :mod:`mass` model,
:class:`MassSolution`\ s are created to store the results computed over the
time interval specified when simulating. These results are divided into two
categories:

    * Concentration solutions (abbreviated as ConcSols)
    * Reaction Flux solutions (abbreviated as FluxSols)

:class:`MassSolution`\ s are therefore given an identifier of the
following format: ``{id_or_model}_{solution_type}Sols`` where ``id_or_model``
and ``solution_type`` correspond to the simulated model and resulting solution
category, respectively.

All solutions in a :class:`MassSolution` can be accessed via attribute
accessors. A :class:`MassSolution` also contains standard ``dict`` methods.

All functions from the :mod:`mass.core.visualization` submodule are designed
to work seamlessly with :class:`MassSolution`\ s, provided they are properly
created.
"""
from warnings import warn

import matplotlib.pyplot as plt

import pandas as pd

from scipy.interpolate import interp1d

from six import iteritems, iterkeys, string_types

from mass.core.mass_model import MassModel
from mass.core.visualization import (
    plot_simulation, plot_tiled_phase_portrait)
from mass.util.dict_with_id import DictWithID
from mass.util.util import ensure_iterable

# Strings of valid solution types
_CONC_STR = "Conc"
_FLUX_STR = "Flux"


class MassSolution(DictWithID):
    r"""Container to store the solutions for the simulation of a model.

    Parameters
    ----------
    id_or_model : str, MassModel
        A :class:`MassModel` or a string identifier to associate with the
        stored solutions. If a :class:`MassModel` is provided, then the model
        identifier will be used.
    solution_type : str
        The type of solution being stored.  Must be  ``'Conc'`` or ``'flux'``.
    data_dict : dict
        A dict containing the solutions to store. If ``None`` provided then the
        :class:`MassSolution` will be initialized with no solutions.
        Solutions can be added or changed later using various ``dict``
        methods (e.g. :meth:`~dict.update`).
    time : array-like
        An array-like object containing the time points used in obtaining the
        solutions.
    interpolate : bool
        If ``True`` then all solutions are converted and stored as
        :class:`scipy.interpolate.interp1d` objects.
        If ``False``, solutions are converted and stored as
        :class:`numpy.ndarray`\ s.

        Default value is ``False``.

    """

    def __init__(self, id_or_model, solution_type=None, data_dict=None,
                 time=None, interpolate=False):
        """Initialize MassSolution."""
        if solution_type not in {_CONC_STR, _FLUX_STR}:
            raise ValueError(
                "'{0}' is not a valid solution type.".format(solution_type))

        if isinstance(id_or_model, MassModel):
            id_or_model = "{0}_{1}Sols".format(str(id_or_model), solution_type)
        if not isinstance(id_or_model, string_types):
            raise TypeError(
                "'id_or_model' must be a MassModel instance or a string")

        super(MassSolution, self).__init__(id=id_or_model, data_dict=data_dict)
        self.solution_type = solution_type
        self._simulation = None
        self._time = time
        self._interpolate = interpolate
        if time is not None:
            self.interpolate = interpolate

    @property
    def simulation(self):
        """Return the associated :class:`~.Simulation`."""
        return getattr(self, "_simulation")

    @property
    def df(self):
        """Return the stored solutions as a :class:`pandas.DataFrame`."""
        sols = dict((k, v(self.time)) if self.interpolate
                    else (k, v) for k, v in iteritems(self))
        df = pd.DataFrame.from_dict(sols)
        df.index = pd.Series(self.time, name="Time")
        return df

    @property
    def time(self):
        r"""Get or set the time points stored in the :class:`MassSolution`.

        Notes
        -----
        If the solutions stored in the :class:`MassSolution` are
        :class:`numpy.ndarray`\ s and the numerical arrays of the solutions
        will be recomputed to correspond to the new time points using
        :class:`scipy.interpolate.interp1d` interpolating functions

        Parameters
        ----------
        value : array-like
            An array-like object containing the time points used in calculating
            the solutions to be stored.

        """
        return getattr(self, "_time", None)

    @time.setter
    def time(self, value):
        """Set the time points that are stored in the :class:`MassSolution`."""
        value = ensure_iterable(value)
        if self.interpolate:
            setattr(self, "_time", value)
        else:
            self.interpolate = True
            setattr(self, "_time", value)
            self.interpolate = False

    @property
    def t(self):
        r"""Shorthand method to get or set the stored time points.

        Notes
        -----
        If the solutions stored in the :class:`MassSolution` are
        :class:`numpy.ndarray`\ s and the numerical arrays of the solutions
        will be recomputed to correspond to the new time points using
        :class:`~scipy.interpolate.interp1d` interpolating functions

        Parameters
        ----------
        value : array-like
            An array-like object containing the time points used in calculating
            the solutions to be stored.

        """
        return self.time

    @t.setter
    def t(self, value):
        """Shorthand to set the time points stored in the MassSolution."""
        self.time = value

    @property
    def interpolate(self):
        """Get or set whether solutions are stored as interpolating functions.

        Parameters
        ----------
        value : bool
            If ``True``, solutions are stored in the :class:`MassSolution` as
            :class:`~scipy.interpolate.interp1d` objects. Otherwise solutions
            are stored as arrays.

        """
        return getattr(self, "_interpolate", None)

    @interpolate.setter
    def interpolate(self, value):
        """Set whether solutions are stored as interpolating functions."""
        if not isinstance(value, bool):
            raise ValueError("value must be a bool")
        if self.time is None:
            warn("No time points associated with MassSolution. Cannot convert "
                 "between numerical arrays and interpolating functions.")
            return

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

        In order to use this method, the :mod:`mass.core.visualization` module
        must be installed.

        See :mod:`~mass.core.visualization` documentation for more information.

        Notes
        -----
        Will clear and use the current axis (accessible via
        :func:`matplotlib.pyplot.gca`).

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

        In order to use this method, the :mod:`mass.core.visualization` module
        must be installed.

        See :mod:`~mass.core.visualization` documentation for more information.

        Notes
        -----
        Will clear and use the current axis (accessible via
        :func:`matplotlib.pyplot.gca`).

        """
        ax = plt.gca()
        ax.cla()
        ax = plot_tiled_phase_portrait(self, ax=ax, empty_tiles=None,
                                       title="Phase portraits for " + self.id)
        ax.get_figure().set_size_inches((7, 7))

    def __getattribute__(self, name):
        """Override of default :func:`getattr` to enable attribute accessors.

        Warnings
        --------
        This method is intended for internal use only.

        """
        if name in self:
            return self[name]
        else:
            return super(MassSolution, self).__getattribute__(name)

    def __dir__(self):
        """Override of default :func:`dir` to include solution accessors.

        Warnings
        --------
        This method is intended for internal use only.

        """
        return list(iterkeys(self)) + super(MassSolution, self).__dir__()


__all__ = ("MassSolution",)