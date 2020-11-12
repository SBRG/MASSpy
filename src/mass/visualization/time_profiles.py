# -*- coding: utf-8 -*-
r"""Contains function for visualizing time profiles of simulation results.

See the  :mod:`mass.visualization` documentation for general information
on :mod:`mass.visualization` functions.

This module contains the following functions for visualization of
time-dependent solutions returned in :class:`~.MassSolution`\ s after
simulation of models.

    * :func:`~.time_profiles.plot_time_profile`
    * :func:`~.time_profiles.plot_ensemble_time_profile`

"""
from warnings import warn

from six import iteritems

from mass.util.util import _check_kwargs
from mass.visualization import visualization_util as v_util


def plot_time_profile(mass_solution, observable=None, ax=None, legend=None, **kwargs):
    """Plot time profiles of solutions in a given :class:`~.MassSolution`.

    Accepted ``kwargs`` are passed onto various :mod:`matplotlib` methods
    utilized in the function. See the :mod:`~mass.visualization` module
    documentation for more detailed information about the possible ``kwargs``.

    Notes
    -----
    * To prevent any changes to the original :class:`~.MassSolution`, a copy of
      the :class:`~.MassSolution` is created and used.

    Parameters
    ----------
    mass_solution : MassSolution
        The :class:`~.MassSolution` containing the time-dependent solutions
        to be plotted.
    observable : iterable, None
        An iterable containing string identifiers of the :mod:`mass` objects
        or the objects themselves that correspond to the keys for the desired
        solutions in the :class:`~.MassSolution`. If ``None`` then all solutions
        are plotted.
    ax : matplotlib.axes.Axes, None
        An :class:`~matplotlib.axes.Axes` instance to plot the data on.
        If ``None`` then the current axes instance is used.
    legend : iterable, str, int
        There are three possible input formats for the legend:

            1. An iterable of legend labels as strings.
            2. A ``str`` representing the location of the legend, or an
               ``int`` between 0 and 14 (inclusive) corresponding to the
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
        * annotate_time_points_labels
        * annotate_time_points_legend
        * annotate_time_points_zorder
        * deviation
        * deviation_zero_centered
        * deviation_normalization

        See :mod:`~mass.visualization` documentation for more information on
        optional ``kwargs``.

    Returns
    -------
    ax :  matplotlib.axes.Axes
        The :class:`~matplotlib.axes.Axes` instance containing the newly
        created plot.

    """
    # Validate whether necessary packages are installed.
    v_util._validate_visualization_packages("matplotlib")
    # Check kwargs
    kwargs = _check_kwargs(get_time_profile_default_kwargs("plot_time_profile"), kwargs)
    # Get the axies instance
    ax = v_util._validate_axes_instance(ax)

    # Validate the MassSolution input, ensure it is not empty.
    mass_solution = v_util._validate_mass_solution(mass_solution)
    if not mass_solution:
        return ax

    # Get the solutions to be observed and validate time vector.
    observable = v_util._validate_plot_observables(mass_solution, observable, **kwargs)

    # Get the plotting function or raise an error if invalid.
    plot_function = v_util._get_plotting_function(
        ax,
        plot_function_str=kwargs.get("plot_function"),
        valid={"plot", "semilogx", "semilogy", "loglog"},
    )

    # Get the legend arguments if desired.
    if legend is not None:
        legend_labels, legend_kwargs = v_util._get_legend_args(
            ax, legend, observable, **kwargs
        )
        observable = v_util._map_labels_to_solutions(observable, legend_labels)
    else:
        legend_kwargs = None

    # Set line colors and styles using a custom cycler
    prop_cycler = v_util._get_line_property_cycler(
        n_current=len(v_util._get_ax_current(ax)), n_new=len(observable), **kwargs
    )
    if prop_cycler:
        ax.set_prop_cycle(prop_cycler)

    # Plot lines onto axes using legend entries as labels (if legend valid).
    for label, sol in iteritems(observable):
        plot_function(observable.time, sol, label=label, zorder=kwargs.get("zorder"))

    # Set the axes options including axis labels, limits, and gridlines.
    v_util._set_axes_labels(ax, **kwargs)
    v_util._set_axes_limits(ax, **kwargs)
    v_util._set_axes_margins(ax, **kwargs)
    v_util._set_axes_gridlines(ax, **kwargs)

    # Set the legend if desired.
    if legend is not None:
        lines_and_labels = v_util._get_handles_and_labels(ax)
        legend = ax.legend(*lines_and_labels, **legend_kwargs)

    # Annotate time points if desired
    if kwargs.get("annotate_time_points", None) is not None:
        ax = v_util._set_annotated_time_points(
            ax,
            observable=observable,
            type_of_plot="time_profile",
            first_legend=(legend, legend_kwargs),
            time_range=observable.time,
            **kwargs
        )
    # Reset default prop_cycle
    ax.set_prop_cycle(v_util._get_default_cycler())

    return ax


def plot_ensemble_time_profile(
    mass_solution_list, observable, ax=None, legend=None, interval_type=None, **kwargs
):
    """Plot time profiles for an ensemble of class:`~.MassSolution` objects.

    The plotted lines represent the mean for the values of a particular
    solution specified in ``observable``.

    Accepted ``kwargs`` are passed onto various :mod:`matplotlib` methods
    utilized in the function. See the :mod:`~mass.visualization` module
    documentation for more detailed information about the possible ``kwargs``.

    Notes
    -----
    * To prevent any changes to the original :class:`~.MassSolution`, copies of
      :class:`~.MassSolution` objects are created and used.

    Parameters
    ----------
    mass_solution_list : iterable
        An iterable of :class:`~.MassSolution` objects containing the
        time-dependent solutions to be plotted.
    observable : iterable
        An iterable containing string identifiers of the :mod:`mass` objects
        or the objects themselves that correspond to the keys for the desired
        solutions in the :class:`~.MassSolution`.
    ax : matplotlib.axes.Axes, None
        An :class:`~matplotlib.axes.Axes` instance to plot the data on.
        If ``None`` then the current axes instance is used.
    legend : iterable, str, int
        There are three possible input formats for the legend:

            1. An iterable of legend labels as strings.
            2. A ``str`` representing the location of the legend, or an
               ``int`` between 0 and 14 (inclusive) corresponding to the
               legend location.
            3. An iterable of the format ``(labels, loc)`` to set both
               the legend labels and location, where ``labels`` and ``loc``
               follows the labels specified in **1** and  **2**.

        See the :mod:`~mass.visualization` documentation for more information
        about legend and valid legend locations.
    interval_type : str, None
        The type of interval to display with the plotted mean of the solution.
        Can be one of the following:

            * ``"range"``: Interval shading occurs from the minimum to the
              maximum value for each time point.
            * ``"CI="``: Interval shading occurs for a confidence interval.
              (e.g. confidence interval of 95% is specified as ``"CI=95.0"``.)
            * ``None`` to prevent interval shading from occurring.

        Default is ``None``
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
        * annotate_time_points_labels
        * annotate_time_points_legend
        * deviation
        * deviation_zero_centered
        * deviation_normalization
        * mean_line_alpha
        * interval_fill_alpha
        * interval_border_alpha
        * CI_distribution

        See :mod:`~mass.visualization` documentation for more information on
        optional ``kwargs``.

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
        get_time_profile_default_kwargs("plot_ensemble_time_profile"), kwargs
    )
    # Get the axes instance and validate the interval type
    ax = v_util._validate_axes_instance(ax)
    interval_type = v_util._validate_interval_type(interval_type)

    # Validate MassSolutions
    mass_solution_list = [
        sol for sol in mass_solution_list if v_util._validate_mass_solution(sol)
    ]
    if not mass_solution_list:
        warn("No valid MassSolution objects given")
        return ax

    # Get the solutions to be observed and validate time vector.
    observable, time_vector = v_util._validate_ensemble_plot_observables(
        mass_solution_list, observable, **kwargs
    )

    # Get the legend arguments if desired.
    if legend is not None:
        legend_labels, legend_kwargs = v_util._get_legend_args(
            ax, legend, observable, **kwargs
        )
        observable = v_util._map_labels_to_solutions(observable, legend_labels)
    else:
        legend_kwargs = None

    # Set line colors and styles using a custom cycler
    prop_cycler = v_util._get_line_property_cycler(
        n_current=len(v_util._get_ax_current(ax)), n_new=len(observable), **kwargs
    )
    if prop_cycler:
        ax.set_prop_cycle(prop_cycler)

    # Plot ensemble results
    _plot_ensemble_lines(ax, observable, interval_type, **kwargs)

    # Set the axes options including axis labels, limits, and gridlines.
    v_util._set_axes_labels(ax, **kwargs)
    v_util._set_axes_limits(ax, **kwargs)
    v_util._set_axes_margins(ax, **kwargs)
    v_util._set_axes_gridlines(ax, **kwargs)

    # Set the legend if desired.
    if legend is not None:
        lines_and_labels = v_util._get_handles_and_labels(ax)
        legend = ax.legend(*lines_and_labels, **legend_kwargs)

    # Annotate time points if desired
    if kwargs.get("annotate_time_points", None):
        ax = v_util._set_annotated_time_points(
            ax,
            observable=observable,
            type_of_plot="time_profile",
            first_legend=(legend, legend_kwargs),
            time_range=time_vector,
            **kwargs
        )

    # Reset default prop_cycle
    ax.set_prop_cycle(v_util._get_default_cycler())

    return ax


def _plot_ensemble_lines(ax, observable, interval_type, **kwargs):
    """Plot the mean and interval of the ensemble solutions.

    Warnings
    --------
    This method is intended for internal use only.

    """

    # Get the plotting function or raise an error if invalid.
    plot_function = v_util._get_plotting_function(
        ax,
        plot_function_str=kwargs.get("plot_function"),
        valid={"plot", "semilogx", "semilogy", "loglog"},
    )

    # Plot lines onto axes using legend entries as labels (if legend valid).
    for label, sol_df in iteritems(observable):
        mean = sol_df.mean(axis=0)

        # Only plot the mean if no interval type specified
        if interval_type is None:
            plot_function(
                sol_df.columns, mean, label=label, alpha=kwargs.get("mean_line_alpha")
            )
            continue

        if interval_type == "range":
            # Get arrays of smallest and largest values per time point
            lb = sol_df.min(axis=0)
            ub = sol_df.max(axis=0)
        else:
            # Compute the confidence interval
            lb, ub, mean = v_util._calculate_confidence_interval(
                interval_type, sol_df, kwargs.get("CI_distribution").lower()
            )

        # Plot mean and get line color
        lines = plot_function(
            sol_df.columns, mean, label=label, alpha=kwargs.get("mean_line_alpha")
        )
        color = lines[0].get_color()

        # Plot lower and upper bound lines
        plot_function(
            sol_df.columns,
            lb,
            label=label + "_lb",
            color=color,
            alpha=kwargs.get("interval_border_alpha"),
        )
        plot_function(
            sol_df.columns,
            ub,
            label=label + "_ub",
            color=color,
            alpha=kwargs.get("interval_border_alpha"),
        )
        # Shade between lower and upper bound lines
        ax.fill_between(
            sol_df.columns,
            lb,
            ub,
            where=lb <= ub,
            facecolor=color,
            alpha=kwargs.get("interval_fill_alpha"),
            interpolate=True,
        )


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
            + str(__all__[:-1])
        )

    default_kwargs = {
        "time_vector": None,
        "plot_function": "plot",
        "title": None,
        "xlabel": None,
        "ylabel": None,
        "xlim": None,
        "ylim": None,
        "xmargin": None,
        "ymargin": None,
        "color": None,
        "linestyle": None,
        "linewidth": None,
        "marker": None,
        "markersize": None,
        "grid": None,
        "grid_color": None,
        "grid_linestyle": None,
        "grid_linewidth": None,
        "legend_ncol": None,
        "annotate_time_points": None,
        "annotate_time_points_color": None,
        "annotate_time_points_marker": None,
        "annotate_time_points_markersize": None,
        "annotate_time_points_labels": False,
        "annotate_time_points_legend": None,
        "annotate_time_points_zorder": None,
        "prop_cycle": None,
        "deviation": False,
        "deviation_zero_centered": False,
        "deviation_normalization": "initial value",
        "zorder": None,
    }

    if function_name == "plot_ensemble_time_profile":
        default_kwargs.update(
            {
                "mean_line_alpha": 1.0,
                "interval_fill_alpha": 0.5,
                "interval_border_alpha": 0.5,
                "CI_distribution": "t",
            }
        )

    return default_kwargs


__all__ = (
    "plot_time_profile",
    "plot_ensemble_time_profile",
    "get_time_profile_default_kwargs",
)
