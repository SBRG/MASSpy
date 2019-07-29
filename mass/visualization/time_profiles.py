# -*- coding: utf-8 -*-
r"""Contains function for visualizing time profiles of simulation results.

See the  :mod:`mass.visualization` documentation for general information
on :mod:`mass.visualization` functions.

This module contains the following functions for visualization of
time-dependent solutions returned in :class:`.MassSolution`\ s after
simulation of models.

    * :func:`~.time_profiles.plot_time_profile`

The following are optional ``kwargs`` that can be passed to the functions
from this module.

    time_vector
    plot_function
    title
    xlabel
    ylabel
    xlim
    ylim
    grid
    grid_color
    grid_linestyle
    grid_linewidth
    prop_cycle
    color
    linestyle
    linewidth
    marker
    markersize
    legend_ncol
    annotate_time_points
    annotate_time_points_color
    annotate_time_points_marker
    annotate_time_points_markersize
    annotate_time_points_legend_loc

"""
from six import iteritems

from mass.util.util import _check_kwargs
from mass.visualization import visualization_util as v_util


def plot_time_profile(mass_solution, observable=None, ax=None, legend=None,
                      **kwargs):
    """Plot time profiles of solutions in a given :class:`.MassSolution`.

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
    observable : iterable
        An iterable containing string identifiers of the :mod:`mass` objects
        or the objects themselves that correspond to the keys for the desired
        solutions in the :class:`.MassSolution`.
    ax : matplotlib.axes.Axes, None
        An :class:`~matplotlib.axes.Axes` instance to plot the data on.
        If ``None`` then the current axes instance is used.
    legend : iterable, str, int
        There are three possible input formats for the legend:

            1. An iterable of legend labels as strings.
            2. A ``str`` representing the location of the legend, or an
               ``int`` between 0 and 10 (inclusive) corresponding to the
               legend location.
            3. An iterable of the format ``(labels, loc)`` to set both
               the legend labels and location, where ``labels`` and ``loc``
               follows the labels specified in **1** and  **2**.

        See the :mod:`~mass.visualization` documentation for more information
        about legend and valid legend locations.
    **kwargs
        * time_vector
        * plot_function
        * title
        * xlabel
        * ylabel
        * xlim
        * ylim
        * grid
        * grid_color
        * grid_linestyle
        * grid_linewidth
        * prop_cycle
        * color
        * linestyle
        * linewidth
        * marker
        * markersize
        * legend_ncol
        * annotate_time_points
        * annotate_time_points_color
        * annotate_time_points_marker
        * annotate_time_points_markersize
        * annotate_time_points_legend_loc

        See :mod:`~mass.visualization.time_profiles` documentation
        for more information on optional ``kwargs``.

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

    # Validate the MassSolution input, ensure it is not empty.
    mass_solution = v_util._validate_mass_solution(mass_solution)
    if not mass_solution:
        return ax

    # Get the solutions to be observed and validate time vector.
    observable = v_util._validate_plot_observables(mass_solution, observable,
                                                   kwargs.get("time_vector"))

    # Get the plotting function or raise an error if invalid.
    plot_function = v_util._get_plotting_function(
        ax, plot_function_str=kwargs.get("plot_function"),
        valid={"plot", "semilogx", "semilogy", "loglog"})

    # Get the legend arguments if desired.
    if legend is not None:
        legend_labels, legend_kwargs = v_util._get_legend_args(ax, legend,
                                                               observable,
                                                               **kwargs)
        observable = v_util._map_labels_to_solutions(observable, legend_labels)
    else:
        legend_kwargs = None

    # Set line colors and styles using a custom cycler
    prop_cycler = v_util._get_line_property_cycler(
        n_current=len(v_util._get_ax_current(ax)), n_new=len(observable),
        **kwargs)
    if prop_cycler:
        ax.set_prop_cycle(prop_cycler)

    # Plot lines onto axes using legend entries as labels (if legend valid).
    for label, sol in iteritems(observable):
        plot_function(observable.time, sol, label=label)

    # Set the axes options including axis labels, limits, and gridlines.
    v_util._set_axes_labels(ax, **kwargs)
    v_util._set_axes_limits(ax, **kwargs)
    v_util._set_axes_gridlines(ax, **kwargs)

    # Set the legend if desired.
    if legend is not None:
        lines_and_labels = v_util._get_handles_and_labels(ax, False)
        legend = ax.legend(*lines_and_labels, **legend_kwargs)

    if kwargs.get("annotate_time_points", None):
        ax = v_util._set_annotated_time_points(
            ax, observable=observable, type_of_plot="time_profile",
            first_legend=(legend, legend_kwargs), **kwargs)
    # Reset default prop_cycle
    ax.set_prop_cycle(v_util._get_default_cycler())

    return ax


def get_time_profile_default_kwargs(function_name):
    """Get default ``kwargs`` for plotting functions in :mod:`time_profiles`.

    Parameters
    ----------
    function_name : str
        The name of the plotting function to get the ``kwargs`` for.
        Valid values include the following:

            * ``"plot_time_profile"``

    Returns
    -------
    dict
        Default ``kwarg`` values for the given ``function_name``.

    """
    if function_name not in __all__[:-1]:
        raise ValueError(
            "Invalid 'function_name'. Valid values include the following: "
            + str(__all__[:-1]))

    default_kwargs = {
        "time_vector": None,
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
        "annotate_time_points": None,
        "annotate_time_points_color": None,
        "annotate_time_points_marker": None,
        "annotate_time_points_markersize": None,
        "annotate_time_points_legend_loc": None,
    }

    return default_kwargs


__all__ = ("plot_time_profile", "get_time_profile_default_kwargs")
