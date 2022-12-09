# -*- coding: utf-8 -*-
r"""Contains function for visualizing phase portraits of simulation results.

See the  :mod:`mass.visualization` documentation for general information
on :mod:`mass.visualization` functions.

This module contains the following functions for visualization of
time-dependent solutions returned in :class:`~.MassSolution`\ s after
simulation of models.

    * :func:`~.phase_portraits.plot_phase_portrait`
    * :func:`~.phase_portraits.plot_ensemble_phase_portrait`
    * :func:`~.phase_portraits.plot_tiled_phase_portrait`

"""
from warnings import warn

import numpy as np
from six import iteritems, iterkeys, itervalues

from mass.util.util import _check_kwargs
from mass.visualization import visualization_util as v_util


def plot_phase_portrait(mass_solution, x, y, ax=None, legend=None, **kwargs):
    """Plot phase portraits of solutions in a given :class:`~.MassSolution`.

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
    x : :mod:`mass` object or its string identifier
        The string identifier of a :mod:`mass` object or the object itself
        that corresponds to the key for the desired solution in the
        :class:`~.MassSolution` for the x-axis of the phase portrait.
    y : :mod:`mass` object or its string identifier
        The string identifier of a :mod:`mass` object or the object itself
        that corresponds to the key for the desired solution in the
        :class:`~.MassSolution` for the y-axis of the phase portrait.
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
               follows the format specified in **1** and  **2**.

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
    kwargs = _check_kwargs(
        get_phase_portrait_default_kwargs("plot_phase_portrait"), kwargs
    )
    # Get the axies instance
    ax = v_util._validate_axes_instance(ax)

    # Validate the MassSolution input, ensure it is not empty.
    mass_solution = v_util._validate_mass_solution(mass_solution)
    if not mass_solution:
        return ax

    # Get the solutions to be observed and validate time vector.
    xy = v_util._validate_plot_observables(mass_solution, (x, y), **kwargs)

    # Get the plotting function or raise an error if invalid.
    plot_function = v_util._get_plotting_function(
        ax,
        plot_function_str=kwargs.get("plot_function"),
        valid={"plot", "semilogx", "semilogy", "loglog"},
    )

    label = "{0} vs. {1}".format(*v_util._group_xy_items(xy, iterkeys))
    sols = v_util._group_xy_items(xy, itervalues)
    observable = {label: sols}

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
    for label, sols in iteritems(observable):
        plot_function(*sols, label=label, zorder=kwargs.get("zorder"))

    # Set the axes options including axis labels, limits, and gridlines.
    v_util._set_axes_labels(ax, **kwargs)
    v_util._set_axes_limits(ax, **kwargs)
    v_util._set_axes_margins(ax, **kwargs)
    v_util._set_axes_gridlines(ax, **kwargs)

    # Set the legend if desired.
    if legend is not None:
        lines_and_labels = v_util._get_handles_and_labels(ax, False)
        legend = ax.legend(*lines_and_labels, **legend_kwargs)

    if kwargs.get("annotate_time_points", None) is not None:
        ax = v_util._set_annotated_time_points(
            ax,
            observable=xy,
            type_of_plot="phase_portrait",
            first_legend=(legend, legend_kwargs),
            time_range=xy.time,
            **kwargs
        )

    # Reset default prop_cycle
    ax.set_prop_cycle(v_util._get_default_cycler())

    return ax


def plot_ensemble_phase_portrait(
    mass_solution_list, x, y, ax=None, legend=None, **kwargs
):
    """Plot a phase portrait for an ensemble of class:`~.MassSolution` objects.

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
    x : :mod:`mass` object or its string identifier
        The string identifier of a :mod:`mass` object or the object itself
        that corresponds to the key for the desired solution in the
        :class:`~.MassSolution` for the x-axis of the phase portrait.
    y : :mod:`mass` object or its string identifier
        The string identifier of a :mod:`mass` object or the object itself
        that corresponds to the key for the desired solution in the
        :class:`~.MassSolution` for the y-axis of the phase portrait.
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
               follows the format specified in **1** and  **2**.

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
    kwargs = _check_kwargs(
        get_phase_portrait_default_kwargs("plot_ensemble_phase_portrait"), kwargs
    )
    # Get the axes instance
    ax = v_util._validate_axes_instance(ax)

    # Validate MassSolutions
    mass_solution_list = [
        sol for sol in mass_solution_list if v_util._validate_mass_solution(sol)
    ]
    if not mass_solution_list:
        warn("No valid MassSolution objects given")
        return ax

    # Get the solutions to be observed and validate time vector.
    xy, time_vector = v_util._validate_ensemble_plot_observables(
        mass_solution_list, (x, y), **kwargs
    )

    label = "{0} vs. {1}".format(*v_util._group_xy_items(xy, iterkeys))
    sols = v_util._group_xy_items(xy, itervalues)
    observable = {label: sols}

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

    _plot_ensemble_lines(ax, observable, **kwargs)

    # Set the axes options including axis labels, limits, and gridlines.
    v_util._set_axes_labels(ax, **kwargs)
    v_util._set_axes_limits(ax, **kwargs)
    v_util._set_axes_margins(ax, **kwargs)
    v_util._set_axes_gridlines(ax, **kwargs)

    # Set the legend if desired.
    if legend is not None:
        lines_and_labels = v_util._get_handles_and_labels(ax, False)
        legend = ax.legend(*lines_and_labels, **legend_kwargs)

    if kwargs.get("annotate_time_points", None):
        ax = v_util._set_annotated_time_points(
            ax,
            observable=xy,
            type_of_plot="phase_portrait",
            first_legend=(legend, legend_kwargs),
            time_range=time_vector,
            **kwargs
        )

    # Reset default prop_cycle
    ax.set_prop_cycle(v_util._get_default_cycler())

    return ax


def _plot_ensemble_lines(ax, observable, **kwargs):
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

    for label, (x_sol_df, y_sol_df) in iteritems(observable):
        plot_function(x_sol_df.mean(axis=0), y_sol_df.mean(axis=0), label=label)


def plot_tiled_phase_portraits(
    mass_solution,
    observable=None,
    ax=None,
    plot_tile_placement="all",
    additional_data=None,
    **kwargs
):
    """Plot phase portraits of solutions in a given :class:`~.MassSolution`.

    Accepted ``kwargs`` are passed onto various matplotlib methods in utilized
    in the function. See the :mod:`~mass.visualization` module documentation
    for more detailed information about the possible ``kwargs``.

    Notes
    -----
    * To prevent any changes to the original :class:`~.MassSolution`, a copy of
      the :class:`~.MassSolution` will be created and used.
    * ``i`` and ``j`` represent the number of rows and columns, respectively.

    Parameters
    ----------
    mass_solution : MassSolution
        The :class:`~.MassSolution` containing the time-dependent solutions
        to be plotted.
    observable : iterable
        An iterable containing string identifiers of the :mod:`mass` objects
        or the objects themselves that correspond to the keys for the desired
        solutions in the :class:`~.MassSolution`.
    ax : matplotlib.axes.Axes, None
        An :class:`~matplotlib.axes.Axes` instance to plot the data on.
        If ``None`` then the current axes instance is used.
    plot_tile_placement : str
        A string representing the location to place the tiles containing
        phase portrait plots. Must be one of the following:

            * ``"lower"`` to place plot tiles on the lower left triangular
              section ``(i < j)`` on the figure tiles.
            * ``"upper"`` to place plot tiles on the upper right triangular
              section ``(i > j)`` on the figure tiles.
            * ``all`` to place plot tiles on the lower left triangular
              section ``(i < j)`` AND on the upper right triangular section
              ``(i > j)`` on the figure tiles.
    additional_data : array_like, None
        A matrix of shape ``(N, N)`` where ``N_obs`` is the number
        of observables provided, or the number of keys in the
        :class:`~.MassSolution` if ``observable=None``. The value at ``(i, j)``
        of the matrix must correspond to the empty tile that the data should
        be displayed on. All other values are ignored. If ``None`` then no
        data will be displayed and tiles will be left empty.
    **kwargs
        * time_vector
        * plot_function
        * title
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
        * annotate_time_points
        * annotate_time_points_color
        * annotate_time_points_marker
        * annotate_time_points_markersize
        * annotate_time_points_legend
        * annotate_time_points_zorder
        * tile_ticks_on
        * tile_xlabel_fontdict
        * tile_ylabel_fontdict
        * data_tile_fontsize
        * data_tile_color
        * diag_tile_color
        * empty_tile_color

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
        get_phase_portrait_default_kwargs("plot_tiled_phase_portraits"), kwargs
    )
    # Get the axes instance
    ax = v_util._validate_axes_instance(ax)

    # Validate the MassSolution input, ensure it is not empty.
    mass_solution = v_util._validate_mass_solution(mass_solution)
    if not mass_solution:
        return ax

    # Get the solutions to be observed and validate time vector.
    observable = v_util._validate_plot_observables(mass_solution, observable, **kwargs)
    if ax.child_axes is not None and len(ax.child_axes) == len(observable) ** 2:
        subaxes = np.reshape(
            np.array(ax.child_axes), (len(observable), len(observable))
        )
    else:
        subaxes = None
    plot_tile_placement = v_util._validate_tile_placement(
        plot_tile_placement, prefix="plot"
    )
    # Split the tiled kwargs and the phase portrait kwargs into seperate dicts
    tile_kwargs, pp_kwargs = _sep_kwargs_for_tiled_phase_portraits(**kwargs)

    # Remove axis lines
    ax.axis("off")
    # Create N x N subplots where N is the number of observable solutions
    # Fraction of larger figure to be utilized by the subplot (Inverse of N)
    sub_ax_placement_vals = [1 / len(observable)] * 2
    # Get width and height. Alter if ticks will be included
    if tile_kwargs.get("tile_ticks_on"):
        sub_ax_placement_vals[1] = sub_ax_placement_vals[1] * 0.75

    for j, y in enumerate(observable):
        for i, x in enumerate(observable):
            # [x0, y0, width, height] from lower left corner of inset axes
            if subaxes is not None:
                sub_ax = subaxes[j, i]
            else:
                sub_ax = ax.inset_axes(
                    bounds=[
                        i * sub_ax_placement_vals[0],
                        1 - sub_ax_placement_vals[0] * (j + 1),
                        sub_ax_placement_vals[1],
                        sub_ax_placement_vals[1],
                    ]
                )
            # Create tile (either phase_portrait, data, or empty)
            sub_ax = _create_tiled_phase_portraits_tile(
                sub_ax,
                observable,
                i,
                j,
                x,
                y,
                plot_tile_placement,
                additional_data,
                tile_kwargs,
                pp_kwargs,
            )
            # If tile_ticks is not True, remove them
            if not tile_kwargs.get("tile_ticks_on"):
                sub_ax.set_xticks([])
                sub_ax.set_yticks([])
            # Set xlabels only on the final row
            if j == len(observable) - 1:
                sub_ax.set_xlabel(x, kwargs.get("tile_xlabel_fontdict"))
            # Set ylabels only on the first column
            if i == 0:
                sub_ax.set_ylabel(y, kwargs.get("tile_ylabel_fontdict"))

    if kwargs.get("annotate_time_points_legend"):
        for sub_ax in ax.get_children():
            if sub_ax.__class__.__name__ == "Axes" and v_util._get_ax_current(
                sub_ax, time_points=True
            ):
                leg_args = v_util._get_annotated_time_points_legend_args(
                    sub_ax, kwargs.get("annotate_time_points_legend")
                )
                break

        ax = v_util._set_additional_legend_box(ax, leg_args, first_legend=None)
    # Set the axes title.
    v_util._set_axes_labels(ax, **kwargs)
    # Reset default prop_cycle
    ax.set_prop_cycle(v_util._get_default_cycler())

    return ax


def get_phase_portrait_default_kwargs(function_name):
    """Get default ``kwargs`` for plotting functions in :mod:`phase_portraits`.

    Parameters
    ----------
    function_name : str
        The name of the plotting function to get the ``kwargs`` for.
        Valid values include the following:

            * ``"plot_phase_portrait"``
            * ``"plot_tiled_phase_portraits"``

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

    if function_name == "plot_phase_portrait":
        default_kwargs.update(
            {
                "xlabel": None,
                "ylabel": None,
                "legend_ncol": None,
            }
        )

    if function_name == "plot_ensemble_phase_portrait":
        default_kwargs.update(
            {
                "xlabel": None,
                "ylabel": None,
                "legend_ncol": None,
            }
        )

    if function_name == "plot_tiled_phase_portraits":
        default_kwargs.update(
            {
                "tile_ticks_on": False,
                "tile_xlabel_fontdict": None,
                "tile_ylabel_fontdict": None,
                "data_tile_fontsize": "large",
                "data_tile_color": None,
                "diag_tile_color": None,
                "empty_tile_color": None,
            }
        )

    return default_kwargs


def _sep_kwargs_for_tiled_phase_portraits(**kwargs):
    """Seperate kwargs for tile properties from kwargs for phase portraits.

    Warnings
    --------
    This method is intended for internal use only.

    """
    # Create dicts for each kwarg
    tile_kwargs = {}
    pp_kwargs = {}
    for key, value in iteritems(kwargs):
        # Get tile kwarg
        if "tile" in key:
            # Add to the tile kwargs
            tile_kwargs[key] = value
        # Get tile kwarge that is not exclusive to tiled phase portraits
        elif key in ["title", "annotate_time_points_legend"]:
            # Add to the tile kwargs
            tile_kwargs[key] = value
        else:
            # Otherwise kwarg belongs to the phase portrait and will be
            # validated in phase portrait function.
            if "margin" in key and value is None:
                value = 0.15
            pp_kwargs[key] = value

    # Iterate through tile colors, setting defaults if no color provided.
    for key, default_color in zip(
        ["data", "diag", "empty"], ["lightgray", "black", "white"]
    ):
        color = tile_kwargs.get(key + "_tile_color")
        if color is None:
            color = default_color
        tile_kwargs[key + "_tile_color"] = color

    # Make endpoints default for tiled phase portrait plots.
    if pp_kwargs.get("annotate_time_points") is None:
        pp_kwargs["annotate_time_points"] = "endpoints"
        if pp_kwargs.get("annotate_time_points_color") is None:
            pp_kwargs["annotate_time_points_color"] = ["red", "blue"]

    return tile_kwargs, pp_kwargs


def _create_tiled_phase_portraits_tile(ax, observable, *args):
    """Create a tile for the tiled phase portrait figure.

    Warnings
    --------
    This method is intended for internal use only.

    """

    def get_plot_tile_bool(i, j, plot_tile_placement):
        """Get a bool indicating if a plot should be made."""
        return {"all": i < j or i > j, "lower": i < j, "upper": i > j}.get(
            plot_tile_placement
        )

    i, j, x, y, plot_tile_placement, data_matrix, tile_kwargs, pp_kwargs = args
    plot_tile_bool = get_plot_tile_bool(i, j, plot_tile_placement)

    # Validate fontsize and set default data tile fontsize as large if needed.
    if not tile_kwargs.get("data_tile_fontsize"):
        tile_kwargs["data_tile_fontsize"] = "large"

    # Set diagonal tile color
    if i == j:
        ax.set_facecolor(tile_kwargs.get("diag_tile_color"))
    # Create a phase portrait for the tile
    elif plot_tile_bool:
        ax = plot_phase_portrait(observable, x=x, y=y, ax=ax, legend=None, **pp_kwargs)
    # Create the data tile
    elif data_matrix is not None and not plot_tile_bool:
        # Create the data tile only if there is information,
        # otherwise set the facecolor as an empty tile
        if data_matrix[j][i] == 0 or not data_matrix[j][i]:
            ax.set_facecolor(tile_kwargs.get("empty_tile_color"))
        else:
            # Place data onto tile and set facecolor
            ax.annotate(
                str(data_matrix[j][i]),
                xy=(0.5, 0.5),
                xycoords="axes fraction",
                va="center",
                ha="center",
                fontsize=tile_kwargs.get("data_tile_fontsize"),
            )
            ax.set_facecolor(tile_kwargs.get("data_tile_color"))
    else:
        ax.set_facecolor(tile_kwargs.get("empty_tile_color"))

    return ax


__all__ = (
    "plot_phase_portrait",
    "plot_ensemble_phase_portrait",
    "plot_tiled_phase_portraits",
    "get_phase_portrait_default_kwargs",
)
