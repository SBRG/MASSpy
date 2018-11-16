# -*- coding: utf-8 -*-
"""TODO Module Docstrings."""
from __future__ import absolute_import

import math
from collections import Hashable, Iterable, OrderedDict
from warnings import warn

from cycler import Cycler

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.axes import Axes

import numpy as np

from six import integer_types, iteritems, iterkeys, string_types

_ZERO_TOL = 1e-6
_FONTSIZES = [
    'xx-small', 'x-small', 'small', 'medium', 'large', 'x-large', 'xx-large']
_LEGEND_LOCS = [
    "best", "upper right", "upper left", "lower left", "lower right", "right",
    "center left", "center right", "lower center", "upper center", "center",
    "left outside", "right outside", "lower outside", "upper outside"]

_custom_default_options = {"plot": {}, "tiled": {}}


def plot_simulation(solution, observable=None, time=None, ax=None, legend=None,
                    **kwargs):
    """Generate a plot of the time profile for the given Solution object.

    Accepted ``kwargs`` are passed on to various matplotlib methods.
    Default ``kwargs`` values can be viewed via get_defaults(plot_type="plot").
    See set_defaults() for a full description of the possible ``kwargs``.

    Parameters
    ----------
    solution: mass.Solution
        The mass.Solution object containing the solutions to be plotted.
    observable: str, iterable of str, optional
        A string or an iterable of strings corresponding to items to display in
        the given 'solution'. Objects can also be provided as long as their ids
        exist as keys in the given 'solution'. If None provided then all items
        in the given 'solution' are displayed.
    time: iterable, optional
        Either an iterable of two values containing the initial and final time
        points, or an iterable of values to treat as time points for the
        solutions. If None provided, the time points used in creating the
        Solution object will be used (accessed via solution.t).
    ax: matplotlib.axes.Axes, optional
        A matplotlib.pyplot.axes instance to plot the data on. If None,
        the current axis is used (accessed via plt.gca()).
    legend: tuple of length 3, str, iterable of str, optional
        If a tuple of length 3 provided, it should be of the following format:
        (list of labels, str for legend loc, str for legend fontsize), allowing
        for full control over the displayed legend.
        If providing a string, can be a single label (when plotting one line),
        the legend location, or the fontsize for the legend.
        If an iterable of strings are provided, they are assumed to be line
        labels unless they correspond to a legend location or legend fontsize.
        If legend is not None but no labels provided, default labels
        representing the keys in the solution object will be used.
        Accepted locations are standard matplotlib legend locations plus
        {"left outside", "right outside", "lower outside", "upper outside"}
        Accepted fontsize values are standard matplotlib legend fontsizes.
        Examples:
            Only legend labels: legend=["A", "B", "C"]
            Only legend properties: legend=("best", "medium")
            Labels and properties: legend=(["A", "B", "C"], "best", "medium")
    **kwargs

    Returns
    -------
    ax: matplotlib.axes.Axes
        A reference to the matplotlib.axes.Axes object used for plotting.

    See Also
    --------
    get_defaults: Default values for options
    set_defaults: Descriptions of accepted kwargs and input formats.

    """
    # Copy the solution object and validate the time input
    solution, time = _format_solution_and_time_input(solution, time)
    # Get observable solutions
    observable_solutions = _set_plot_observables(solution, time, observable)
    # Get current axis if none provided, otherwise check its type
    if ax is None:
        ax = plt.gca()
    if not isinstance(ax, Axes):
        raise TypeError("ax must be a matplotlib.axes Axes object.")

    # Update plot options with provided kwargs
    options = _update_kwargs("plot", **kwargs)
    # Set the plotting function.
    plot_function = _set_plot_function(ax, options["plot_function"])

    # Use the property cycler for the axis if provided
    if options["prop_cycle"] is not None:
        ax.set_prop_cycle(options["prop_cycle"])

    # Use the plotting function to plot the observable data.
    default_labels = []
    for label, sol in iteritems(observable_solutions):
        plot_function(time, sol, label=label)
        default_labels += [label]

    # Set axis options which include labels, limits, and gridlines.
    _set_axis_options(ax, options)
    # Set line colors and styles if no cycler provided.
    if options["prop_cycle"] is None:
        _set_line_properties(ax, options, len(default_labels))

    # Set the legend
    if legend is not None:
        _set_axis_legend(ax, legend, default_labels, options)

    return ax


def plot_phase_portrait(solution, x, y, time=None, ax=None, legend=None,
                        **kwargs):
    """Generate a plot of the time profile for the given Solution object.

    Accepted ``kwargs`` are passed on to various matplotlib methods.
    Default ``kwargs`` values can be viewed via get_defaults(plot_type="plot").
    See get_defaults() for a full description of the possible ``kwargs``.

    Parameters
    ----------
    solution: mass.Solution
        The mass.Solution object containing the solutions to be plotted.
    x: str, iterable of str
        A string or an iterable of strings corresponding to items in
        the given 'solution' to plot as the x-axis. Objects can also be
        provided as long as their ids exist as keys in the given 'solution'.
    y: str, iterable of str
        A string or an iterable of strings corresponding to items in
        the given 'solution' to plot as the y-axis. Objects can also be
        provided as long as their ids exist as keys in the given 'solution'.
    time: iterable, optional
        Either an iterable of two values containing the initial and final time
        points, or an iterable of values to treat as time points for the
        solutions. If None provided, the time points used in creating the
        Solution object will be used (accessed via solution.t).
    ax: matplotlib.axes.Axes, optional
        A matplotlib.pyplot.axes instance to plot the data on. If None,
        the current axis is used (accessed via plt.gca()).
    legend: tuple of length 3, str, iterable of str, optional
        If a tuple of length 3 provided, it should be of the following format:
        (list of labels, str for legend loc, str for legend fontsize), allowing
        for full control over the displayed legend.
        If providing a string, can be a single label (when plotting one line),
        the legend location, or the fontsize for the legend.
        If an iterable of strings are provided, they are assumed to be line
        labels unless they correspond to a legend location or legend fontsize.
        If legend is not None but no labels provided, default labels
        representing the keys in the solution object will be used.
        Accepted locations are standard matplotlib legend locations plus
        {"left outside", "right outside", "lower outside", "upper outside"}
        Accepted fontsize values are standard matplotlib legend fontsizes.
        Examples:
            Only legend labels: legend=["A", "B", "C"]
            Only legend properties: legend=("best", "medium")
            Labels and properties: legend=(["A", "B", "C"], "best", "medium")
    **kwargs

    Returns
    -------
    ax: matplotlib.axes.Axes
        A reference to the matplotlib.axes.Axes object used for plotting.

    See Also
    --------
    get_defaults: Default values for options
    set_defaults: Descriptions of accepted kwargs and input formats.

    """
    # Copy the solution object and validate the time input
    solution, time = _format_solution_and_time_input(solution, time)
    # Get observable solutions
    x_sols = _set_plot_observables(solution, time, x)
    y_sols = _set_plot_observables(solution, time, y)
    # Get current axis if none provided, otherwise check its type
    if ax is None:
        ax = plt.gca()
    if not isinstance(ax, Axes):
        raise TypeError("ax must be a matplotlib.axes Axes object.")

    # Update plot options with provided kwargs
    options = _update_kwargs("plot", **kwargs)
    # Set the plotting function
    plot_function = _set_plot_function(ax, options["plot_function"])

    # Use the property cycler for the axis if provided
    if options["prop_cycle"] is not None:
        ax.set_prop_cycle(options["prop_cycle"])

    # Use the plotting function to plot the observable data
    default_labels = []
    for x_label, x_sol in iteritems(x_sols):
        for y_label, y_sol in iteritems(y_sols):
            label = "{0} vs. {1}".format(x_label, y_label)
            plot_function(x_sol, y_sol, label=label)
            default_labels += [label]

    # Set axis options which include labels, limits, and gridlines.
    _set_axis_options(ax, options)
    # Set line colors and styles if no cycler provided.
    if options["prop_cycle"] is None:
        _set_line_properties(ax, options, len(default_labels))

    # Set the legend
    if legend is not None:
        _set_axis_legend(ax, legend, default_labels, options)

    return ax


def plot_tiled_phase_portrait(solution, observable=None, time=None, ax=None,
                              display_data=None, **kwargs):
    """TODO DOCSTRING."""
    solution, time = _format_solution_and_time_input(solution, time)
    observable_solutions = _set_plot_observables(solution, time, observable)
    n_sols = len(observable_solutions)

    tile_options = _update_kwargs("tiled", **kwargs)
    print(tile_options)
    # Create figure with N x N subplots where N is the number of osolutions
    fig, axes = plt.subplots(nrows=n_sols, ncols=n_sols,
                             figsize=tile_options["figsize"])

    for i, y in enumerate(observable_solutions):
        for j, x in enumerate(observable_solutions):
            # Get the corresponding axis
            ax = axes[i][j]
            # Set common axis labels
            if i == n_sols - 1:
                ax.set_xlabel(x, fontdict={"size": tile_options["fontsize"]})
            if j == 0:
                ax.set_ylabel(y, fontdict={"size": tile_options["fontsize"]})

            if i == j:
                ax.set_facecolor(tile_options["diag_color"])

            elif tile_options["placement"] == "upper":
                if i < j:
                    ax = plot_phase_portrait(solution, x, y, time, ax=ax)
                elif display_data is not None and i > j:
                    ax.annotate(str(display_data[i][j]), xy=(0.5, 0.5),
                                xycoords="axes fraction", va="center",
                                ha="center", fontsize=tile_options["fontsize"])
                    ax.set_facecolor(tile_options["data_color"])
                else:
                    ax.set_facecolor(tile_options["none_color"])

            elif tile_options["placement"] == "lower":
                if i > j:
                    ax = plot_phase_portrait(solution, x, y, time, ax=ax)
                elif display_data is not None and i < j:
                    ax.annotate(str(display_data[i][j]), xy=(0.5, 0.5),
                                xycoords="axes fraction", va="center",
                                ha="center", fontsize=tile_options["fontsize"])
                    ax.set_facecolor(tile_options["data_color"])
                else:
                    ax.set_facecolor(tile_options["none_color"])
            else:
                ax = plot_phase_portrait(solution, x, y, time, ax=ax)

            ax.set_xticklabels([])
            ax.set_yticklabels([])

    return axes


def get_default_colors(n_items, start=0):
    """Return a list of the colors for a given number of items.

    This function is primarily for internal use.

    Parameters
    ----------
    n_items: int
        Number of items that need a color.
    start: int
        The index of where to start the cmap.

    Returns
    -------
    list:
        A list of of (n_items - start) colors.

    """
    for parameter, value in zip(["n_items", "start"], [n_items, start]):
        if not isinstance(value, integer_types):
            raise TypeError(parameter + " must be an integer.")

    if n_items > 10 and n_items <= 20:
        cmap = mpl.cm.tab20.colors
    elif n_items > 20 and n_items <= 60:
        cmap = (mpl.cm.tab20.colors + mpl.cm.tab20b.colors
                + mpl.cm.tab20c.colors)
    elif n_items > 60:
        cmap = mpl.cm.nipy_spectral(np.linspace(0, 1, n_items))
    else:
        cmap = mpl.cm.tab10.colors

    return list(cmap)[start:n_items]


def get_defaults(plot_type):
    """Get the default values for options of a plotting function.

    Parameters
    ----------
    plot_type: str
        The type of plotting function used. Accepted 'plot_type' values are
        {"plot", "tiled"}, with "plot" corresponding to the 'plot_simulation'
        and 'plot_phase_portrait' functions, and with "tiled" corresponding to
        the 'plot_tiled_phase_portrait' function.

    Returns
    -------
    options: dict
        A dict containing key:value pairs of the possible options and their
        default values for a given "plot_type".

    See Also
    --------
    set_defaults: Descriptions of accepted kwargs and input formats.

    """
    if plot_type not in {"plot", "tiled"}:
        raise ValueError("Unrecognized 'plot_type' value.")

    options = {
        "plot_function": "plot",
        "grid": None,
        "prop_cycle": None,
        "linecolor": None,
        "linestyle": None,
    }
    # Update dict for single plot specific options
    if _custom_default_options[plot_type]:
        options = _custom_default_options[plot_type]
    elif plot_type is "plot":
        options.update({
            "title": (None, {"size": "medium"}),
            "xlabel": (None, {"size": "medium"}),
            "ylabel": (None, {"size": "medium"}),
            "xlim": (None, None),
            "ylim": (None, None),
            "default_legend_loc": "best",
            "default_legend_fontsize": "medium",
            "default_legend_ncol": None,
        })
    # Update dict for tiled specific options
    elif plot_type is "tiled":
        options.update({
            "figsize": (8., 8.),
            "placement": "both",
            "fontsize": "medium",
            "diag_color": "black",
            "data_color": "gray",
            "none_color": "white",
        })
    else:
        pass

    return options


def set_defaults(plot_type, **kwargs):
    """Get the default values for options of a plotting function.

    Parameters
    ----------
    plot_type: str
        The type of plotting function used. Accepted 'plot_type' values are
        {"plot", "tiled"}, with "plot" corresponding to the 'plot_simulation'
        and 'plot_phase_portrait' functions, and with "tiled" corresponding to
        the 'plot_tiled_phase_portrait' function.
    kwargs:
        plot_function: str
            A string representing the plotting function to use.
            Accepted values are {"plot", "semilogx", "semilogy", "loglog"}.
                "plot" uses an linear x-axis and a linear y-axis
                "semilogx" uses an logarithmic x-axis and a linear y-axis
                "semilogy" uses an linear x-axis and a logarithmic y-axis
                "loglog" uses an logarithmic x-axis and a logarithmic y-axis
            Default value is "plot".
        grid: tuple, bool
            If providing a tuple, the expected format is ("which", "axis"),
            where the values for "which" and "axis" are passed to the
            corresponding argument passed to axes.grid. Accepted values include
            {"major", "minor", "both"} for "which" and {"x", "y, "both"} for
            "axis". If True is provided, gridlines are set with using default
            values for "which" and "axis", and gridlines are removed if False.
        prop_cycle: cylcer.Cyler
            Set a property cycler to control colors and other style properties
            for multiline plots. The property cycler is validated through the
            matplotlib.rcsetup.cycler function and set through the
            axes.set_prop_cycle function. Providing a property cycler will
            cause "linecolor" and "linestyle" kwargs to be ignored.
        linecolor: str, iterable of str, None
            A string or an iterable of strings representing matplotlib.colors
            to use as line colors. If a single string is provided, the color
            will be applied to all solutions to be plotted. If an iterable of
            strings is provided, it must be equal to the number of lines
            to be newly plotted. All colors are validated using the
            matplotlib.colors.is_color_like function. If None provided, will
            use colors obtained from the get_default_colors function.
            Will be ignored if a property cycler is provided via prop_cycle.
        linestyle: str, iterable of str, None
            A string or an iterable of strings representing matplotlib
            recognized linestyles. If a single string is provided, the
            linestyle will be applied to all solutions to be plotted. If an
            iterable of strings is provided, it must be equal to the number of
            lines to be newly plotted. If None provided, will use solid lines.
            Will be ignored if a property cycler is provided via prop_cycle.
        title: string, tuple
            Either a string to use as a title, or a tuple where the first value
            is the title as a string, and the second value is a dict of font
            options to pass to axes.set_title function. Keys must be one of the
            six matplotlib.font_manager font properties. When setting new
            defaults, the new value must be the tuple.
        xlabel: string, tuple
            Either a string to use as the x-axis label, or a tuple where the
            first value is the label as a string, and the second value is a
            dict of font options to pass to axes.set_xlabel function. Keys must
            be one of the six matplotlib.font_manager font properties. When
            setting new defaults, the new value must be the tuple.
            Only valid for plot_type="plot".
        ylabel: string, tuple
            Either a string to use as the y-axis label, or a tuple where the
            first value is the label as a string, and the second value is a
            dict of font options to pass to axes.set_ylabel function. Keys must
            be one of the six matplotlib.font_manager font properties. When
            setting new defaults, the new value must be the tuple.
            Only valid for plot_type="plot".
        xlim: tuple
            A tuple of integers or floats of form (xmin, xmax) specifying the
            limits of the x-axis. Args are passed to axes.set_xlim function.
            Only valid for plot_type="plot".
        ylim: tuple
            A tuple of integers or floats of form (ymin, ymax) specifying the
            limits of the y-axis. Args are passed to axes.set_ylim function.
            Only valid for plot_type="plot".
        default_legend_loc: str
            The default legend location if no location was provided to the
            legend arg of the plotting function. Ideally, this kwarg should not
            be used directly in the plotting function, but as a kwarg for
            setting a new default location for the set_defaults function.
            Only valid for plot_type="plot".
        default_legend_fontsize: float, str
            The default legend fontsize if no fontsize was provided to the
            legend arg of the plotting function. Ideally, this kwarg should not
            be used directly in the plotting function, but as a kwarg for
            setting a new default fontsize for the set_defaults function.
            Valid values include
            Only valid for plot_type="plot".
        legend_ncol: int, optional
            An int representing the number of columns for the legend. If
            None, then ncols=int(math.ceil(math.sqrt(n_items)/3)).
            Only valid for plot_type="plot".
        figsize:
            Only valid for plot_type="tiled".
        placement:
            Only valid for plot_type="tiled".
        fontsize:
            Only valid for plot_type="tiled".
        diag_color:
            Only valid for plot_type="tiled".
        none_color:
            Only valid for plot_type="tiled".

    Returns
    -------
    dict:
        A dict containing key:value pairs of the possible options and their
        new default values for a given "plot_type".

    See Also
    --------
    get_defaults: Default values for options

    """
    if plot_type not in {"plot", "tiled"}:
        raise ValueError("Unrecognized 'plot_type' value.")

    if kwargs:
        options = _update_kwargs(plot_type, **kwargs)
        _custom_default_options[plot_type].update(options)
    else:
        _custom_default_options[plot_type] = {}

    return _custom_default_options[plot_type]


# Internal
def _format_solution_and_time_input(solution, time):
    """Copy the solution object and format the time input (Helper function).

    Warnings
    --------
    This method is intended for internal use only.

    """
    # Copy the solution object to prevent modifications to the original
    solution = solution.copy()

    # Use time stored in Solution if None provided
    if time is None:
        time = solution.t
    # Create an array of time points to use if time bounds provided.
    elif len(time) == 2:
        # TODO Find a better method of determining numpoints to use.
        if abs(time[0]) < _ZERO_TOL:
            time = (_ZERO_TOL, time[-1])
        time = np.geomspace(time[0], time[1], solution.numpoints)
    # Use the array of time points provided
    elif isinstance(time, Iterable) and not isinstance(time, string_types):
        time = np.array(sorted(time))

    # If time input is not recognized, then raise an error
    else:
        raise TypeError("Unrecognized 'time' input")

    return solution, time


def _set_plot_observables(solution, time, observable):
    """Set plot observables and return the solution values (Helper function).

    Warnings
    --------
    This method is intended for internal use only.

    """
    # Return all items in the solution if no observables are provided.
    if observable is None:
        observable = solution.solutions

    # If a single observable is provided, make it iterable.
    if not hasattr(observable, '__iter__') or \
       isinstance(observable, string_types):
        observable = [observable]

    # Replace objects providfed in observable with their identifiers.
    observable = [getattr(x, "id", x) for x in observable]

    # Check to ensure specified observables are in the Solution object
    if not set(observable).issubset(set(iterkeys(solution))):
        raise ValueError("observable must keys from the solution")

    # Turn solutions into interpolating functions if the timepoints provided
    # are not identical to those used in the simulation.
    if solution.t.all() != time.all():
        solution.interpolate = True

    observable = OrderedDict((x, solution[x])
                             if isinstance(solution[x], np.ndarray)
                             else (x, solution[x](time)) for x in observable)

    return observable


def _set_plot_function(ax, option):
    """Set the plotting function to be used (Helper function).

    Warnings
    --------
    This method is intended for internal use only.

    """
    plotting_functions_dict = {"plot": ax.plot,
                               "semilogx": ax.semilogx,
                               "semilogy": ax.semilogy,
                               "loglog": ax.loglog}

    return plotting_functions_dict[option]


def _set_axis_options(ax, options):
    """Set axis labels, titles, limits, and gridlines (Helper function).

    Warnings
    --------
    This method is intended for internal use only.

    """
    # Set the xlabel, ylabel, and the title if provided.
    for key, setter in zip(["xlabel", "ylabel", "title"],
                           [ax.set_xlabel, ax.set_ylabel, ax.set_title]):
        if options[key][0] is not None:
            setter(options[key][0], fontdict=options[key][1])

    # Set the xlim and the ylim if provided.
    for key, setter in zip(["xlim", "ylim"], [ax.set_xlim, ax.set_ylim]):
        if options[key] != (None, None):
            limits = list(options[key])
            if 0. in limits:
                plot_func = options["plot_function"]
                if key == "xlim" and plot_func in ["loglog", "semilogx"] \
                   or key == "ylim" and plot_func in ["loglog", "semilogy"]:
                    limits[limits.index(0.)] = _ZERO_TOL
            setter(tuple(limits))

    # Set the gridlines if provided.
    if options["grid"] is not None:
        if isinstance(options["grid"], bool):
            ax.grid(options["grid"])
        else:
            ax.grid(True, which=options["grid"][0], axis=options["grid"][1])


def _set_line_properties(ax, options, n_new):
    """Set line colors and styles (Helper function).

    Warnings
    --------
    This method is intended for internal use only.

    """
    # Get the current lines on the plot.
    lines = [line for line in ax.get_lines()]
    # Get default values for linestyles and linecolors.
    default_dict = {
        "linecolor": get_default_colors(len(lines), start=len(lines) - n_new),
        "linestyle": ["-"]*n_new
        }

    for key in ["linecolor", "linestyle"]:
        # Set the property if given
        if options[key] is not None:
            to_apply = options[key].copy()
        # Skip linestyles if not given
        elif options[key] is None and key is "linestyle":
            continue
        else:
            to_apply = default_dict[key]

        # Apply to all newly plotted items if only one value provided.
        if len(to_apply) == 1 and n_new != 1:
            to_apply = to_apply*n_new
        # Use default values if the input does not match the number of newly
        # plotted lines.
        elif len(to_apply) != n_new:
            warn("Wrong number of {0}s provided, will use default {0}s "
                 "instead.".format(key))
            to_apply = default_dict[key]
        else:
            pass
        # Set from the last line to the first line to
        # preserve previously plotted lines.
        to_apply.reverse()
        for i, value in enumerate(to_apply):
            if key is "linestyle":
                lines[len(lines)-1-i].set_linestyle(value)
            else:
                lines[len(lines)-1-i].set_color(value)


def _set_axis_legend(ax, legend, default_labels, options):
    """Set axis legend and associated properties (Helper function).

    Warnings
    --------
    This method is intended for internal use only.

    """
    # Initialize default values for fontsize, legend location
    font = options["default_legend_fontsize"]
    loc = options["default_legend_loc"]
    ncol = options["default_legend_ncol"]
    legend = list(v if isinstance(v, string_types) else list(v)
                  for v in legend)
    # Parse through legend input if provided as a tuple.
    if len(legend) <= 3:
        while legend:
            brk = False
            font, brk = _get_legend_values(legend, _FONTSIZES, font, brk)
            loc, brk = _get_legend_values(legend, _LEGEND_LOCS, loc, brk)
            if legend and isinstance(legend[0], string_types):
                legend = [legend]
            items, brk = _get_legend_values(legend, [], default_labels, brk)
            if legend and not brk:
                if isinstance(legend[-1], string_types):
                    msg = ("Could not interpret legend parameter '{0}'"
                           .format(legend[-1]))
                else:
                    msg = ("'{0}' required != '{1}' provided".format(
                           len(default_labels), len(legend[-1])))
                msg += ", will attempt to use default legend values instead."
                warn(msg)
                del legend[-1]
            if not legend:
                break
    # Set the legend if only labels are provided.
    elif (all(isinstance(entry, string_types) for entry in legend) and
          len(legend) == len(default_labels)):
        items = legend
    # Use default values if input could not be interpreted.
    else:
        warn("Could not interpret the legend input, will use the default "
             "entries and parameters instead.")
        items = default_labels
    # Get the legend properties
    property_dict = _get_legend_properties(options, loc, ncol, len(items))
    property_dict.update({"fontsize": font})

    # Apply new labels to the lines and make/update the legend.
    lines = [line for line in ax.get_lines()]
    for new_label, line in zip(items, lines[len(lines) - len(items):]):
        line.set_label(new_label)

    labels = [line.get_label() for line in lines]

    ax.legend(lines, labels, **property_dict)


def _update_kwargs(plot_type, **kwargs):
    """Validate kwargs and update their corresponding values (Helper function).

    Warnings
    --------
    This method is intended for internal use only.

    """
    options = get_defaults(plot_type)
    # Return default options if no kwargs provided.
    if kwargs:
        update_functions_dict = {
            "plot_function": _update_plot_function,
            "xlabel": _update_labels,
            "ylabel": _update_labels,
            "title": _update_labels,
            "xlim": _update_limits,
            "ylim": _update_limits,
            "grid": _update_grid,
            "linecolor": _update_lines,
            "linestyle": _update_lines,
            "prop_cycle": _update_lines,
            "default_legend_loc": _update_legend_properties,
            "default_legend_fontsize": _update_legend_properties,
            "default_legend_ncol": _update_legend_properties,
            "placement": _update_tiles,
            "diag_color": _update_tiles,
            "data_color": _update_tiles,
            "none_color": _update_tiles,
        }
        # Iterate through kwargs and update options.
        for key, value in iteritems(kwargs):
            if key in update_functions_dict and key in options:
                update_functions_dict[key](options, key, value)
            else:
                warn("kwarg '" + str(key) + "' not recognized")
    return options


def _update_plot_function(options, key, value):
    """Validate the 'plot_function' kwarg (Helper function).

    Warnings
    --------
    This method is intended for internal use only.

    """
    allowed_plotting_functions = {"loglog", "semilogx", "semilogy", "plot"}
    if value not in allowed_plotting_functions:
        raise ValueError(str(key) + "must be one of the following: " +
                         str(allowed_plotting_functions))
    options[key] = value


def _update_labels(options, key, value):
    """Validate kwargs for axis labels (Helper function).

    Warnings
    --------
    This method is intended for internal use only.

    """
    # Check input types for option
    if not hasattr(value, '__iter__'):
        raise TypeError(str(key) + " must be an iterable")
    if isinstance(value, string_types):
        value = (value, options[key][1])
    if not isinstance(value[1], dict):
        raise TypeError("fontdict must be a dictionary")
    # Update options
    options[key] = value


def _update_limits(options, key, value):
    """Validate kwargs for axis limits (Helper function).

    Warnings
    --------
    This method is intended for internal use only.

    """
    if not hasattr(value, '__iter__') or len(value) != 2:
        raise TypeError(str(key) + " must be an iterable of length 2")
    elif not all(isinstance(v, (float, integer_types))for v in value):
        raise TypeError("Limits must be ints or floats")
    elif value[0] > value[1]:
        warn("{0}min is greater than {0}max, attempting to swap values"
             .format(key[0]))
        value = (value[1], value[0])
    else:
        value = tuple(value)
    # Update options
    options[key] = value


def _update_grid(options, key, value):
    """Validate the 'grid' kwarg (Helper function).

    Warnings
    --------
    This method is intended for internal use only.

    """
    which_args = ["major", "minor", "both"]
    axis_args = ["x", "y", "both"]

    if isinstance(value, string_types):
        if value in which_args:
            value = (value, "both")
        else:
            value = ("both", value)

    if isinstance(value, bool):
        pass
    elif len(value) == 2 and all(isinstance(v, string_types) for v in value):
        if value[0] not in which_args:
            raise ValueError("The first value must be one of the following: "
                             + str(which_args))
        if value[1] not in axis_args:
            raise ValueError("The second value must be one of the following: "
                             + str(axis_args))
    else:
        raise TypeError(str(key) + " must be string or a tuple of strings")

    options[key] = value


def _update_lines(options, key, value):
    """Validate kwargs for line properties (Helper function).

    Warnings
    --------
    This method is intended for internal use only.

    """
    is_color_like = mpl.colors.is_color_like
    if isinstance(value, string_types):
        value = [value]

    if isinstance(value, Cycler):
        value = mpl.rcsetup.cycler(value)
    elif (not isinstance(value, Iterable) or
          not all(isinstance(v, string_types) for v in value)):
        raise TypeError(key + " must be a string or an iterable of strings.")
    elif key is "linecolor" and not all(is_color_like(v) for v in value):
        raise ValueError("Invalid colors found: " + str(
                         [v for v in value if not is_color_like(v)]))
    else:
        pass
    options[key] = value


def _update_legend_properties(options, key, value):
    """Validate kwargs for legend properties (Helper function).

    Warnings
    --------
    This method is intended for internal use only.

    """
    allowed_dict = {"default_legend_loc": _LEGEND_LOCS,
                    "default_legend_fontsize": _FONTSIZES}
    if not isinstance(value, (integer_types, float)):
        if key is "default_legend_ncol":
            raise TypeError(key + " must be a numerical value.")
        if value not in allowed_dict[key]:
            raise ValueError(key + " must be one of the following: "
                             + str(allowed_dict[key]))
    elif value < 0:
        raise ValueError(key + " must be a non negative value")
    elif key is "default_legend_loc" and value > 10:
        raise ValueError(key + " must be a value between 0 and 10, inclusive.")
    else:
        pass

    options[key] = value


def _get_legend_values(legend, allowable, default, success):
    """Help set legend location, fontsize, and entries (Helper function).

    Warnings
    --------
    This method is intended for internal use only.

    """
    parameter = None
    if legend:
        if allowable:
            if isinstance(legend[-1], Hashable) and \
               (isinstance(legend[-1], (integer_types, float)) or
               legend[-1] in allowable):
                parameter = legend.pop(-1)
                success = True
        elif (not isinstance(legend[-1], Hashable) and
              len(legend[-1]) == len(default)):
                parameter = legend.pop(-1)
                success = True
        else:
            parameter = default
    else:
        parameter = default

    if parameter is None:
        parameter = default
    return parameter, success


def _get_legend_properties(options, loc, ncol, n_items):
    """Help set legend location and number of columns (Helper function).

    Warnings
    --------
    This method is intended for internal use only.

    """
    loc_anchor_dict = {"left outside": ("center right", (-0.15, 0.5)),
                       "right outside": ("center left", (1.05, 0.5)),
                       "lower outside": ("upper center", (0.5, -0.2)),
                       "upper outside": ("lower center", (0.5, 1.15))}
    if loc in loc_anchor_dict:
        loc, anchor = loc_anchor_dict[loc]
    else:
        anchor = None

    if ncol is None:
        ncol = int(math.ceil(math.sqrt(n_items)/3))

    return {"loc": loc, "bbox_to_anchor": anchor, "ncol": ncol}


def _update_tiles(options, key, value):
    """Validate kwargs for tile placement and color (Helper function).

    Warnings
    --------
    This method is intended for internal use only.

    """
    if key is "placement" and value not in {"both", "lower", "upper"}:
        raise ValueError(key + " must be one of the following: "
                         + str({'both', 'lower', 'upper'}))
    if "_color" in key and not mpl.colors.is_color_like(value):
        raise ValueError("Invalid color input: " + str(value))

    options[key] = value
