# -*- coding: utf-8 -*-
r"""MassSolution is a class for storing the simulation results.

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

All functions from the :mod:`mass.visualization` submodule are designed
to work seamlessly with :class:`MassSolution`\ s, provided they are properly
created.
"""
from warnings import warn


try:
    import matplotlib.pyplot as plt
except ImportError:
    plt = None
import pandas as pd
from numpy import array
from scipy.interpolate import interp1d
from six import iteritems, iterkeys, itervalues, string_types
from sympy import Symbol, lambdify, sympify

from mass.core.mass_configuration import MassConfiguration
from mass.core.mass_model import MassModel
from mass.util.dict_with_id import DictWithID
from mass.util.util import apply_decimal_precision, ensure_iterable
from mass.visualization import plot_tiled_phase_portraits, plot_time_profile
from mass.visualization.visualization_util import _validate_visualization_packages


# Strings of valid solution types
_CONC_STR = "Conc"
_FLUX_STR = "Flux"

MASSCONFIGURATION = MassConfiguration()


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
        A dict containing the solutions to store. If ``None`` provided then
        the :class:`MassSolution` will be initialized with no solutions.
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

    def __init__(
        self,
        id_or_model,
        solution_type="",
        data_dict=None,
        time=None,
        interpolate=False,
        initial_values=None,
    ):
        """Initialize MassSolution."""
        if isinstance(id_or_model, MassModel):
            id_or_model = "{0}_{1}Sols".format(str(id_or_model), solution_type)
        if not isinstance(id_or_model, string_types):
            raise TypeError("'id_or_model' must be a MassModel instance or a string")

        super(MassSolution, self).__init__(id=id_or_model, data_dict=data_dict)
        self.solution_type = solution_type
        self._simulation = None
        self._time = time
        self._interpolate = interpolate
        if time is not None:
            self.interpolate = interpolate
        self._initial_values = {}
        self.initial_values = initial_values

    @property
    def simulation(self):
        """Return the associated :class:`~.Simulation`."""
        return getattr(self, "_simulation")

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
        """Shorthand method to set the stored time points."""
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
            warn(
                "No time points associated with MassSolution. Cannot convert "
                "between numerical arrays and interpolating functions."
            )
            return

        for key, sol in iteritems(self):
            if value and not isinstance(sol, interp1d):
                self[key] = interp1d(
                    self.time, sol, kind="cubic", fill_value="extrapolate"
                )
            if not value and isinstance(sol, interp1d):
                self[key] = sol(self.time)

        setattr(self, "_interpolate", value)

    @property
    def initial_values(self):
        """Get or set a ``dict`` of the initial values for solution variables.

        Notes
        -----
        Primarily for storing the intial values used in a simulation and for
        calculating deviations from the initial state.

        Parameters
        ----------
        initial_values : dict
            A ``dict`` containining the variables stored in the
            :class:`MassSolution` and their initial values.

        """
        return getattr(self, "_initial_values")

    @initial_values.setter
    def initial_values(self, initial_values):
        """Set the initial values for solution variables."""
        new_values = {}
        if initial_values is not None:
            for key, value in iteritems(initial_values):
                if key not in self:
                    continue
                new_values[key] = apply_decimal_precision(
                    value, MASSCONFIGURATION.decimal_precision
                )
        setattr(self, "_initial_values", new_values)

    def view_time_profile(self, deviation=False, plot_function="loglog"):
        """Generate a quick view of the time profile for the solution.

        See :mod:`~mass.visualization` documentation for more information.

        Parameters
        ----------
        deviation : bool
            Whether to plot time profiles as a deviation from their
            initial value.
        plot_function : str
            The plotting function to use. Accepted values are the following:

              * ``"plot"`` for a linear x-axis and a linear y-axis
                via :meth:`Axes.plot() <matplotlib.axes.Axes.plot>`
              * ``"loglog"`` for a logarithmic x-axis and a logarithmic y-axis
                via :meth:`Axes.loglog() <matplotlib.axes.Axes.loglog>`
              * ``"semilogx``" for a logarithmic x-axis and a linear y-axis
                via :meth:`Axes.semilogx() <matplotlib.axes.Axes.semilogx>`
              * ``"semilogy"`` for a linear x-axis and a logarithmic y-axis
                via :meth:`Axes.semilogy() <matplotlib.axes.Axes.semilogy>`

        Notes
        -----
        Will clear and use the current axis (accessible via
        :func:`matplotlib.pyplot.gca`).

        """
        _validate_visualization_packages("matplotlib")
        # Get solution type for title
        solution_type = ""
        if self.solution_type == _CONC_STR:
            solution_type += " Concentrations"

        if self.solution_type == _FLUX_STR:
            solution_type += " Fluxes"

        # Set plot options
        options = {
            "plot_function": plot_function,
            "grid": ("major", "x"),
            "title": "Time Profile for {0}{1}".format(self.id, solution_type),
            "xlabel": "Time",
            "ylabel": solution_type,
            "deviation": deviation,
        }
        # Get current axis and clear it
        ax = plt.gca()
        ax.cla()
        # Plot time profile
        ax = plot_time_profile(self, ax=ax, legend="right outside", **options)
        # Set figure size
        ax.get_figure().set_size_inches((6, 4))

    def view_tiled_phase_portraits(self):
        """Generate a preview of the phase portraits for the solution.

        See :mod:`~mass.visualization` documentation for more information.

        Notes
        -----
        Will clear and use the current axis (accessible via
        :func:`matplotlib.pyplot.gca`).

        """
        _validate_visualization_packages("matplotlib")
        # Get current axis
        ax = plt.gca()
        ax.cla()
        # Set options
        options = {
            "title": ("Tiled Phase portraits for " + self.id, {"size": "large"}),
            "annotate_time_points": "endpoints",
            "annotate_time_points_color": ["r", "b"],
            "annotate_time_points_marker": ["o", "D"],
            "tile_xlabel_fontdict": {"size": "large"},
            "tile_ylabel_fontdict": {"size": "large"},
            "annotate_time_points_legend": "right outside",
        }
        # Plot with kwargs
        ax = plot_tiled_phase_portraits(self, ax=ax, **options)
        # Set figure size
        ax.get_figure().set_size_inches((7, 7))

    def to_frame(self):
        """Return the stored solutions as a :class:`pandas.DataFrame`."""
        # Ensure solutions are data points
        sols = dict(
            (k, v(self.time)) if self.interpolate else (k, v)
            for k, v in iteritems(self)
        )
        # Make dataframe and set time as index
        try:
            df = pd.DataFrame.from_dict(sols)
        except ValueError:
            df = pd.DataFrame(itervalues(sols), index=iterkeys(sols))
            df.index.name = "ID"
        else:
            df.index = pd.Series(self.time, name="Time")

        return df

    def make_aggregate_solution(
        self, aggregate_id, equation, variables=None, parameters=None, update=True
    ):
        """Make a new aggregate variable and its solution from an equation.

        Parameters
        ----------
        aggregate_id : str
            An identifier for the solution to be made.
        equation : str
            A string representing the equation of the new solution.
        variables : iterable or None
            Either an iterable of object identifiers or the objects themselves
            representing keys in the :class:`MassSolution` that are used as
            variables in equation. If ``None``, then all keys of the solution
            object are checked as variables, potentially resulting in lower
            performance time.
        parameters : dict or None
            A ``dict`` of additional parameters to use, where key:value pairs
            are the parameter identifiers and their numerical values.
            If ``None`` then it is assumed that there are no additional
            parameters in the equation.
        update : bool
            Whether to add the new solution into the :class:`MassSolution`.
            via the :meth:`~dict.update` method. Default is ``True``.

        Returns
        -------
        solution : dict
            A ``dict`` containing where the key is the ``aggregate_id`` and the
            value is the newly created solution as the same type as the
            variable solutions

        Raises
        ------
        SympifyError
            Raised if the ``equation_str`` could not be interpreted.

        """
        if variables is None:
            variables = list(iterkeys(self))
        else:
            variables = sorted([getattr(var, "_id", var) for var in variables])
            invalid = [var for var in variables if var not in self]
            if invalid:
                raise ValueError("'{0!r}' not found in MassSolution".format(invalid))
        local_syms = list(variables)
        if parameters is None:
            parameters = {}
        else:
            local_syms += list(parameters)

        equation = sympify(equation, locals={k: Symbol(k) for k in local_syms})
        equation = lambdify(args=variables, expr=equation.subs(parameters))

        if self.time is not None:
            values = array(
                [self[k](self.time) if self.interpolate else self[k] for k in variables]
            )
        else:
            values = array([self[k] for k in variables])

        solution = equation(*values)
        # # Return solution as type in the MassSolution
        if self.interpolate and self.time is not None:
            solution = interp1d(
                self.time, solution, kind="cubic", fill_value="extrapolate"
            )
        elif solution.size == 1:
            solution = solution.item()

        solution = {aggregate_id: solution}

        # Update MassSolution if d
        if update:
            self.update(solution)
            try:
                self.initial_values[aggregate_id] = equation(
                    *array([self.initial_values[k] for k in variables])
                )
            except KeyError:
                # Not all initial values have been set
                pass

        return solution

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
