# -*- coding: utf-8 -*-
r"""Contains function for visualizing time profiles of simulation results.

This module contains the following functions for visualization of
time-dependent solutions returned in :class:`.MassSolution`\ s after
simulation of models.

    * :func:`~.time_profiles.plot_time_profile`

Note that to use the :mod:`mass.visualization` functions, additional packages
for visualziation included during the :mod:`mass` installation process as
follows::

    # Installs plotting visualization packages.
    pip install mass["plotting"]
    # Or to install all additional packages.
    pip install mass["all"]


In general, it is recommended to pass an :class:`matplotlib.axes.Axes` instance
into the visualization function in order to guaruntee more control over the
generated plot and ensure that the plot be placed onto a figure as expected.
If an :class:`~matplotlib.axes.Axes` is not passed into the function, then the
currrent axes instance will be used, accessed via::

    matplotlib.pyplot.gca()

All functions will return the axes instance utilzed in generating the plot.

The following are optional ``kwargs`` that can be passed to the functions
from this module.

"""
from six import iteritems

from mass.util.util import _check_kwargs
from mass.visualization import visualization_util as v_util


def plot_time_profile(mass_solution, ax=None, observable=None, legend=None,
                      **kwargs):
    """Plot solution time profile(s) in a given :class:`.MassSolution`.

    Accepted ``kwargs`` are passed onto various matplotlib methods in utilized
    in the function. See the :mod:`.visualization.time_profiles` module
    documentation for more detailed information about the possible ``kwargs``.

    Notes
    -----
    * To prevent any changes to the original :class:`.MassSolution`, a copy of
      the :class:`.MassSolution` will be created and used.

    Parameters
    ----------
    mass_solution : MassSolution
        The :class:`.MassSolution` containing the time-dependent solutions
        to be plotted.
    ax : matplotlib.axes.Axes, None
        An :class:`~matplotlib.axes.Axes` instance to plot the data on.
        If ``None`` then the current axes instance is used.
    observable : list
        A list containing string identifiers of the :mod:`mass` objects or
        the objects themselves that correspond to the keys for the desired
        solutions in the  :class:`.MassSolution`.
    legend : iterable, str
        There are three possible inputs for the legend:

            1. An iterable of legend entries as strings.
            2. A str representing the location of the legend.
            3. An iterable of the format ``(labels, loc)`` to set both
               the legend entries and location, where ``entries`` and ``loc``
               follows the format specified in **1** and  **2**.

        If only legend labels are provided, the default legend location
        ``"best"`` will be used.

        If only the legend location is provided, the corresponding keys in
        ``mass_solution`` will be used as default legend labels.
    **kwargs

    Returns
    -------
    ax :  matplotlib.axes.Axes
        The :class:`~matplotlib.axes.Axes` instance containing the newly
        created plot.

    """
    # Validate whether necessary packages are installed.
    v_util._validate_visualization_packages("matplotlib")
    # Check kwargs
    kwargs = _check_kwargs(
        get_time_profile_default_kwargs("plot_time_profile"),
        kwargs)

    # Get the axies instance
    ax = v_util._validate_axes_instance(ax)

    # Make a copy of the MassSolution and validate time input (if any).
    mass_solution = v_util._validate_mass_solution(mass_solution)
    if not mass_solution:
        return ax

    time_points = v_util._validate_time_vector(kwargs.get("time_points"),
                                               mass_solution.time)

    # Get the solutions to be observed.
    observable = v_util._validate_plot_observables(mass_solution, time_points,
                                                   observable)

    # Get the plotting function and plot the observable data
    plot_function = v_util._get_plotting_function(
        ax, plot_function_str=kwargs.get("plot_function"),
        valid={"plot", "semilogx", "semilogy", "loglog"})

    legend_labels, legend_kwargs = v_util._get_legend_args(ax, legend,
                                                           observable,
                                                           **kwargs)
    observable = v_util._map_labels_to_solutions(observable, legend_labels)

    # Set line colors and styles using a custom cycler
    v_util._set_line_properties(ax, n_new=len(observable), **kwargs)

    # Plot lines onto axes using legend entries as labels (if legend valid).
    for label, sol in iteritems(observable):
        plot_function(time_points, sol, label=label)

    # Set the axes options including axis labels, limits, and gridlines.
    v_util._set_axes_labels(ax, **kwargs)
    v_util._set_axes_limits(ax, **kwargs)
    v_util._set_axes_gridlines(ax, **kwargs)
    ax.legend(*ax.get_legend_handles_labels(), **legend_kwargs)
    # Reset default prop_cycle
    ax.set_prop_cycle(v_util._get_default_cycler())

    return ax


def get_time_profile_default_kwargs(function):
    """Get default kwargs for a plotting function in :mod:`time_profiles`.

    Parameters
    ----------
    function_name : str
        The name of the plotting function to get the kwargs for. Valid values
        include the following:

            * plot_time_profile

    Returns
    -------
    dict
        Default kwarg values for the given ``function_name``.

    """
    if function not in __all__[:-1]:
        raise ValueError(
            "Invalid 'function_name'. Valid values include the following: "
            + str(__all__[:-1]))

    default_kwargs = {
        "time_points": None,
        "plot_function": "plot",
        "title": None,
        "xlabel": None,
        "ylabel": None,
        "xlim": None,
        "ylim": None,
        "grid": None,
        "grid_color": None,
        "grid_linestyle": None,
        "grid_linewidth": None,
        "prop_cycle": None,
        "color": None,
        "linestyle": None,
        "linewidth": None,
        "marker": None,
        "markersize": None,
        "legend_ncol": None,
    }

    return default_kwargs


__all__ = ("plot_time_profile", "get_time_profile_default_kwargs")
