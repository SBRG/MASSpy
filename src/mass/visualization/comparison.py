# -*- coding: utf-8 -*-
"""Contains function for visually comparing values in various objects.

See the  :mod:`mass.visualization` documentation for general information
on :mod:`mass.visualization` functions.

This module contains the following functions for visually comparing a set of
values in one object against a similar set of valeus in another object.

    * :func:`~.comparison.plot_comparison`

"""
import math

from mass.util.util import _check_kwargs
from mass.visualization import visualization_util as v_util


def plot_comparison(
    x, y, compare=None, observable=None, ax=None, legend=None, **kwargs
):
    """Plot values of two objects for comparision.

    This function can take two :class:`.MassModel`, :class:`.ConcSolution`,
    :class:`cobra.Solution <cobra.core.solution.Solution>`, or
    :class:`pandas.Series` objects and plot them against one another in a
    calibration plot.

    Accepted ``kwargs`` are passed onto various :mod:`matplotlib` methods
    utilized in the function. See the :mod:`~mass.visualization` module
    documentation for more detailed information about the possible ``kwargs``.

    Notes
    -----
    * If a :class:`pandas.Series`, the index must correspond to the identifier
      of the assoicated object. (e.g. a metabolite identifier for
      ``compare="concentrations"``, or a reaction identifier for
      ``compare="Keqs"``)

    Parameters
    ----------
    x : MassModel, ConcSolution, ~cobra.core.solution.Solution, ~pandas.Series
        The object to access for x-axis values.
    y : MassModel, ConcSolution, ~cobra.core.solution.Solution, ~pandas.Series
        The object to access for y-axis values.
    compare : str
        The values to be compared. Must be one of the following:

            * ``"concentrations"`` for :class:`.MassModel` and
              :class:`.ConcSolution` objects.
            * ``"Keqs"`` for :class:`.MassModel` and
              :class:`.ConcSolution` objects.
            * ``"fluxes"`` for :class:`.MassModel` and
              :class:`cobra.Solution <cobra.core.solution.Solution>` objects.
            * ``"kfs"`` for :class:`.MassModel` objects.

        Not required if both ``x`` and ``y`` are :class:`pandas.Series`.
    observable : iterable
        An iterable containing string identifiers of :mod:`mass` objects
        or the objects themselves corresponding to the object or index where
        the value is located.
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
        * marker
        * markersize
        * legend_ncol
        * xy_line
        * xy_linecolor
        * xy_linewidth
        * xy_linestyle
        * xy_legend

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
    kwargs = _check_kwargs(get_comparison_default_kwargs("plot_comparison"), kwargs)
    kwargs["linestyle"] = " "
    # Get the axies instance
    ax = v_util._validate_axes_instance(ax)

    x = v_util._get_values_as_series(x, compare, name="x")
    y = v_util._get_values_as_series(y, compare, name="y")
    observable = v_util._get_dataframe_of_observables(x, y, compare, observable)

    # Get the plotting function or raise an error if invalid.
    plot_function = v_util._get_plotting_function(
        ax, plot_function_str=kwargs.get("plot_function"), valid={"plot", "loglog"}
    )

    # Get the legend arguments if desired.
    if legend is not None:
        legend_labels, legend_kwargs = v_util._get_legend_args(
            ax, legend, list(observable.index), **kwargs
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
    for sol in observable.itertuples():
        plot_function(sol.x, sol.y, label=sol.Index)

    # Set the axes options including axis labels, limits, and gridlines.
    v_util._set_axes_labels(ax, **kwargs)
    v_util._set_axes_limits(ax, **kwargs)
    v_util._set_axes_margins(ax, **kwargs)
    v_util._set_axes_gridlines(ax, **kwargs)

    # Set the legend if desired.
    if legend is not None:
        lines_and_labels = v_util._get_handles_and_labels(ax, False)
        legend = ax.legend(*lines_and_labels, **legend_kwargs)

    if kwargs.get("xy_line"):
        # Get min and max value, make line extend slightly past those values.
        if kwargs.get("plot_function") == "plot":
            limits = (
                int(math.floor(observable.values.min()) / 1.05 - 1),
                int(math.ceil(observable.values.max()) * 1.05 + 1),
            )
        else:
            limits = (
                int(math.floor(observable.values.min()) / 10),
                int(math.ceil(observable.values.max()) * 10),
            )
        ax = _plot_xy_line(ax, limits, first_legend=(legend, legend_kwargs), **kwargs)
    # Reset default prop_cycle
    ax.set_prop_cycle(v_util._get_default_cycler())

    return ax


def _plot_xy_line(ax, limits, first_legend=None, **kwargs):
    """Plot a a line for ``y=x`` on the comparison plot.

    Warnings
    --------
    This method is intended for internal use only.

    """
    # Validate color
    linecolor = kwargs.get("xy_linecolor")
    if linecolor is None:
        linecolor = "grey"

    # Validate linestyle
    linestyle = kwargs.get("xy_linestyle")
    if linestyle is None:
        linestyle = "--"

    # Validate linewidth
    linewidth = kwargs.get("xy_linewidth")

    plot_function = v_util._get_plotting_function(
        ax, plot_function_str=kwargs.get("plot_function"), valid={"plot", "loglog"}
    )

    # Plot the line using the set kwarg options, set zorder to 1.9 so that
    # it is below original points (default is 2.)
    line = plot_function(
        limits,
        limits,
        label="y=x",
        color=linecolor,
        linestyle=linestyle,
        linewidth=linewidth,
        marker="",
        zorder=1.9,
    )

    if kwargs.get("xy_legend") is not None:
        desired, taken = v_util._check_second_legend_location(
            kwargs.get("xy_legend"), first_legend[1]
        )
        # Set default desired location
        if desired is None:
            desired = "best" if taken != "best" else "right outside"
        # Get kwargs for legend location
        anch = None
        if desired in v_util.OUTSIDE_LEGEND_LOCATION_AND_ANCHORS:
            desired, anch = v_util.OUTSIDE_LEGEND_LOCATION_AND_ANCHORS[desired]

        legend_args = (line, ["y=x"], {"loc": desired, "bbox_to_anchor": anch})
        ax = v_util._set_additional_legend_box(
            ax, legend_args, first_legend=first_legend[0]
        )

    return ax


def get_comparison_default_kwargs(function_name):
    """Get default ``kwargs`` for plotting functions in :mod:`comparison`.

    Parameters
    ----------
    function_name : str
        The name of the plotting function to get the ``kwargs`` for.
        Valid values include the following:

            * ``"plot_comparison"``

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
        "plot_function": "plot",
        "title": None,
        "xlabel": None,
        "ylabel": None,
        "xlim": None,
        "ylim": None,
        "xmargin": None,
        "ymargin": None,
        "color": None,
        "marker": "o",
        "markersize": None,
        "grid": None,
        "grid_color": None,
        "grid_linestyle": None,
        "grid_linewidth": None,
        "legend_ncol": None,
        "prop_cycle": None,
        "xy_line": False,
        "xy_linecolor": None,
        "xy_linewidth": None,
        "xy_linestyle": None,
        "xy_legend": None,
    }

    return default_kwargs


__all__ = ("plot_comparison", "get_comparison_default_kwargs")
