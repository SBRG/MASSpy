# -*- coding: utf-8 -*-
"""Contains internal functions common to various visualization functions.

Warnings
--------
The functions found in this module are NOT intended for direct use.

"""
import math
from collections.abc import Iterable
from warnings import warn

try:
    from cycler import cycler

    import matplotlib as mpl
    from matplotlib import rcsetup as rc
    import matplotlib.pyplot as plt
except ImportError:
    mpl = None
    rc = None
    cycler = None

import numpy as np

from six import integer_types, iteritems, iterkeys, itervalues, string_types

from mass.util.util import ensure_iterable

INSTALLED_VISUALIZATION_PACKAGES = {
    "matplotlib": bool(mpl),
    "cycler": bool(cycler),
}
"""dict: Package names and ``bool``\ s indicating if they are installed."""

OUTSIDE_LEGEND_LOCATION_AND_ANCHORS = {
    "left outside": ("center right", (-0.15, 0.5)),
    "right outside": ("center left", (1.05, 0.5)),
    "lower outside": ("upper center", (0.5, -0.2)),
    "upper outside": ("lower center", (0.5, 1.15)),
}
"""dict: Legend location and anchors for default outside legend locations."""


def _validate_visualization_packages(package):
    """Validate whether a visualization package has been installed.

    Warnings
    --------
    This method is intended for internal use only.

    """
    if not INSTALLED_VISUALIZATION_PACKAGES[package]:
        raise ValueError(
            "The " + package + " has not been installed. To utilize the "
            "mass visualization functions, all visualization packages must be"
            "installed.")


def _validate_axes_instance(ax):
    """Validate the given axes instance vector and return it.

    Warnings
    --------
    This method is intended for internal use only.

    """
    if ax is None:
        ax = plt.gca()
    elif not isinstance(ax, mpl.axes.Axes):
        raise TypeError("ax must be a matplotlib.axes.Axes instance.")

    return ax


def _validate_mass_solution(mass_solution):
    """Validate the mass solution input and return a copy of it.

    Warnings
    --------
    This method is intended for internal use only.

    """
    if mass_solution.__class__.__name__ != "MassSolution":
        raise TypeError("mass_solution must be a mass.MassSolution.")

    if not mass_solution:
        warn("The MassSolution is empty.")

    return mass_solution.copy()


def _validate_time_vector(time_points, default_time_points):
    """Validate the given vector of time points and return it.

    Warnings
    --------
    This method is intended for internal use only.

    """
    if time_points is None:
        return np.array(default_time_points)
    if isinstance(time_points, Iterable)\
       and not isinstance(time_points, string_types):
        return np.array(sorted(time_points))

    raise ValueError("Invalid input for 'time_points'.")


def _validate_plot_observables(mass_solution, time_points, observable):
    """Validate the given observables and return them and their solutions.

    Warnings
    --------
    This method is intended for internal use only.

    """
    # Return all items in the solution if no observables are provided.
    if observable is None:
        observable = list(iterkeys(mass_solution))
    else:
        # Otherwise ensure observable is iterable
        observable = ensure_iterable(observable)

    # Replace mass objects with their identifiers.
    observable = [getattr(x, "id", x) for x in observable]

    # Check to ensure specified observables are in the MassSolution
    if not set(observable).issubset(set(iterkeys(mass_solution))):
        raise ValueError("observable must keys from the mass_solution.")

    # Turn solutions into interpolating functions if the timepoints provided
    # are not identical to those used in the simulation.
    if not isinstance(time_points, np.ndarray):
        time_points = np.array(time_points)

    if not np.array_equal(mass_solution.time, time_points):
        mass_solution.interpolate = True

    observable = dict(
        (x, mass_solution[x]) if isinstance(mass_solution[x], np.ndarray)
        else (x, mass_solution[x](time_points)) for x in observable)

    return observable


def _validate_legend_input_fmt(legend, observable):
    """Validate whether the given ``legend`` can be interpreted.

    Warnings
    --------
    This method is intended for internal use only.

    """
    # Ensure legend is iterable and get default labels and legend location
    legend = ensure_iterable(legend)
    default_loc = rc.defaultParams['legend.loc'][0]

    # Get lengths of observables and legend arguments for comparisions.
    n_obs, n_leg = len(observable), len(legend)

    # Get the possible location value in the input.
    possible_loc = legend[-1]

    def _validate_legend_loc(legend, poss_loc, def_labels, def_loc):
        """Valdidate legend location in this function."""
        try:
            if poss_loc not in OUTSIDE_LEGEND_LOCATION_AND_ANCHORS:
                poss_loc = _validate_kwarg_input("legend_loc", poss_loc,
                                                 as_warning=False)
        except ValueError:
            items = legend, def_loc
        else:
            items = def_labels, poss_loc

        return items

    # Legend location is an integer
    if isinstance(possible_loc, (integer_types, float)):
        if legend[0] != possible_loc and\
           (isinstance(legend[0], list) or n_obs <= n_leg <= n_obs + 1):
            items = ensure_iterable(legend[0]), possible_loc
        elif legend[0] == possible_loc:
            items = list(iterkeys(observable)), possible_loc
        else:
            items = None, possible_loc

    # Legend format includes iterable of labels and a location as a str
    elif n_leg == 2 and (not isinstance(legend[0], string_types)
                         and isinstance(possible_loc, string_types)):
        items = legend[0], legend[-1]
    # Only one observable value or legend value
    elif 1 in [n_obs, n_leg]:
        # Legend has one label and one location
        if n_leg == n_obs + 1:
            items = ensure_iterable(legend[0]), possible_loc
        elif n_leg in [n_obs, 1]:
            items = _validate_legend_loc(
                legend, possible_loc, list(iterkeys(observable)), default_loc)
    # Either only labels provided or bad label input and location
    elif n_leg == n_obs:
        items = _validate_legend_loc(legend, possible_loc, None, default_loc)
    # Bad legend input, return None to set off warnings
    else:
        items = None, None

    items = ensure_iterable(items)
    # Ensure number of legend labels is equal to the number of observables
    valid_bool = bool(items[0] is not None and len(items[0]) == n_obs)
    items[0] = items[0] if valid_bool else None

    # Ensure legend location is valid
    try:
        items[1] = _validate_kwarg_input("legend_loc", items[1],
                                         as_warning=False)
    except (AttributeError, ValueError):
        valid_bool = bool(items[1] in OUTSIDE_LEGEND_LOCATION_AND_ANCHORS
                          or items[1] in list(range(0, 11)))
        items[1] = items[1] if valid_bool else None

    return items


def _get_legend_args(ax, legend, observable, **kwargs):
    """Validate whether the given ``legend`` can be interpreted.

    Warnings
    --------
    This method is intended for internal use only.

    """
    n_current = len([line for line in ax.get_lines()
                     if line.get_label()[:2] != "t="])

    # Validate legend input format
    legend_labels, legend_loc = _validate_legend_input_fmt(legend, observable)
    if legend_labels is None:
        _raise_kwarg_warning(
            "legend",
            msg="therefore utilizing keys from the MassSolution instead")
        legend_labels = list(iterkeys(observable))

    if legend_loc is None:
        _raise_kwarg_warning(
            "loc", msg="therefore utilizing default legend 'loc' instead",
            kwarg_prefix="legend")
        legend_loc = rc.defaultParams['legend.loc'][0]

    if legend_loc in OUTSIDE_LEGEND_LOCATION_AND_ANCHORS:
        legend_loc, anchor = OUTSIDE_LEGEND_LOCATION_AND_ANCHORS[legend_loc]
    else:
        anchor = None

    # Get the legend properties
    ncols = _validate_kwarg_input("legend_ncol", kwargs.get("legend_ncol"))
    if ncols is None:
        ncols = int(math.ceil(math.sqrt(len(observable) + n_current) / 3))

    legend_kwargs = {
        "loc": legend_loc, "bbox_to_anchor": anchor, "ncol": ncols}

    return legend_labels, legend_kwargs


def _map_labels_to_solutions(observable, legend_labels):
    """Map the legend labels to the observable solutions and return.

    Warnings
    --------
    This method is intended for internal use only.

    """
    # Map solutions to labels
    if legend_labels == list(iterkeys(observable)):
        legend_labels = iteritems(observable)
    else:
        legend_labels = zip(legend_labels, itervalues(observable))

    return {label.lstrip("_"): sol for label, sol in legend_labels}


def _validate_kwarg_input(arg_name, arg_value, prefix=None, as_warning=True):
    """Validate whether the given ``input`` can be interpreted.

    Warnings
    --------
    This method is intended for internal use only.

    """
    KWARG_VALIDATION_FUNCTIONS = {
        "cycler": rc.cycler,
        "color": rc.validate_color,
        "linestyle": rc._validate_linestyle,
        "linewidth": rc.validate_float,
        "marker": rc.validate_string,
        "markersize": rc.validate_float,
        "axis_locator": rc.validate_axis_locator,
        "grid_axis": rc.validate_grid_axis,
        "legend_loc": rc.validate_legend_loc,
        "legend_ncol": rc.validate_int_or_None,
    }
    # Allow None to be returned
    if arg_value is not None:
        # Check input validity based on the arg_name
        validator = KWARG_VALIDATION_FUNCTIONS[arg_name]
        try:
            arg_value = validator(arg_value)
        except ValueError as e:
            # Return as error if warning disabled
            if not as_warning:
                raise ValueError(e)
            # Otherwise raise warning and set arg_value to None
            _raise_kwarg_warning(arg_name, msg=e, kwarg_prefix=prefix)
            arg_value = None

    return arg_value


def _get_plotting_function(ax, plot_function_str, valid):
    """Get the plotting function to utilize.

    Warnings
    --------
    This method is intended for internal use only.

    """
    plotting_functions_dict = {
        "plot": ax.plot,
        "semilogx": ax.semilogx,
        "semilogy": ax.semilogy,
        "loglog": ax.loglog
    }

    if plot_function_str not in plotting_functions_dict\
       or plot_function_str not in valid:
        raise ValueError(
            "'{0}' not a valid plotting option. Must be one of the "
            "following options: '{1}'.".format(plot_function_str, valid))

    return plotting_functions_dict[plot_function_str]


def _set_line_properties(ax, n_new, **kwargs):
    """Set the line properties.

    Warnings
    --------
    This method is intended for internal use only.

    """
    # Get the prop_cycler
    prop_cycler = kwargs.get("prop_cycler")

    if prop_cycler is None:
        cycler_kwargs_values = {}
        # Get current number of lines on the plot, excluding time points.
        n_current = len([line for line in ax.get_lines()
                         if line.get_label()[:2] != "t="])
        default_cycler_values = {
            "color": _get_default_colors(n_current + n_new)[n_current:],
            "linestyle": [rc.defaultParams["lines.linestyle"][0]] * n_new,
            "linewidth": [rc.defaultParams["lines.linewidth"][0]] * n_new,
            "marker": [rc.defaultParams["lines.marker"][0]] * n_new,
            "markersize": [rc.defaultParams["lines.markersize"][0]] * n_new
        }

        for k in default_cycler_values:
            values = kwargs.get(k)
            # Don't add to cycler instance if it was not provided.
            if not values:
                continue

            values = ensure_iterable(values)
            msg = ""
            # Validate values using appropriate validation method.
            try:
                values = [_validate_kwarg_input(k, v, as_warning=False)
                          for v in values]
            except ValueError as e:
                # Catch error and append to warning message
                msg += str(e)

            # If only one value given, apply it to all lines.
            if len(values) == 1 and n_new != 1:
                values = values * n_new
            # Raise warning if the number of given values does not match the
            # number of new items. Do not need to include this message
            # if invalid kwarg values were found.
            if len(values) != n_new and not msg:
                msg += "Wrong number of values for '{0}' provided".format(k)

            if msg:
                # Raise warning if a warning message was created and utilize
                # default values for that kwarg.
                warn(msg + ". Therefore utilizing default values instead.")
                cycler_kwargs_values[k] = default_cycler_values[k]
            else:
                # Add to cycler initialization args.
                cycler_kwargs_values[k] = values

        if "color" not in cycler_kwargs_values:
            cycler_kwargs_values["color"] = default_cycler_values["color"]
        prop_cycler = cycler(**cycler_kwargs_values)
        # Validate cycler and ensure matplotlib format. Though the earlier
        # validation functions should have caught all issues, this extra
        # step catches potential edge cases that may have slipped through.
        prop_cycler = _validate_kwarg_input("cycler", prop_cycler)
    else:
        prop_cycler = _validate_kwarg_input("cycler", prop_cycler)

    # Set cycler if there is one.
    if prop_cycler:
        ax.set_prop_cycle(prop_cycler)


def _set_axes_labels(ax, **kwargs):
    """Set the axes labels and title.

    Warnings
    --------
    This method is intended for internal use only.

    """
    for label_type in ["xlabel", "ylabel", "title"]:
        label_values = kwargs.get(label_type)
        if label_values is None:
            # No label provided
            continue
        if isinstance(label_values, string_types):
            # Only label string provided
            label_str = label_values
            fontdict = {}
        elif len(label_values) == 2:
            # Label string and fontdict provided
            label_str, fontdict = label_values
        else:
            _raise_kwarg_warning(label_type)
            label_str = None
        # Try setting the label
        if label_str is not None:
            try:
                getattr(ax, "set_" + label_type)(label_str, fontdict=fontdict)
            except ValueError as e:
                warn("Could not set '{0}' due to the following: '{1}'.".format(
                    label_type, e))


def _set_axes_limits(ax, **kwargs):
    """Set the axes limits.

    Warnings
    --------
    This method is intended for internal use only.

    """
    for limit_type in ["xlim", "ylim"]:
        limit_values = kwargs.get(limit_type)
        if limit_values is None:
            # No limits provided
            continue
        if not hasattr(limit_values, "__iter__") or len(limit_values) != 2:
            # Invalid limit input provided
            _raise_kwarg_warning(limit_type, msg="'{0}'".format(limit_values))
            limit_values = None
        # Try setting the label
        if limit_values is not None:
            try:
                getattr(ax, "set_" + limit_type)(tuple(limit_values))
            except ValueError as e:
                warn("Could not set '{0}' due to the following: '{1}'".format(
                    limit_type, e))


def _set_axes_gridlines(ax, **kwargs):
    """Set the axes gridlines.

    Warnings
    --------
    This method is intended for internal use only.

    """
    grid = kwargs.get("grid")

    if grid is not None:
        # Validate additional grid argument inputs
        grid_options = {}
        for k in ["color", "linestyle", "linewidth"]:
            prefix = "grid_"
            v = kwargs.get(prefix + k)
            v = _validate_kwarg_input(k, v, prefix=prefix)
            if v is not None:
                grid_options[k] = v

        if isinstance(grid, bool):
            ax.grid(grid, **grid_options)

        elif isinstance(grid, Iterable) and len(grid) == 2:
            which, axis = grid
            # Validate which argument for gridlines
            msg = ""
            try:
                _validate_kwarg_input("axis_locator", which, as_warning=False)
            except ValueError as e:
                # Format warning message
                msg += str(e).split(":")[-1].strip()
                msg += " for grid 'which' argument"
                which = None

            # Validate axis argument for gridlines
            try:
                _validate_kwarg_input("grid_axis", axis, as_warning=False)
            except ValueError as e:
                # Format warning message
                if msg:
                    msg += " and "
                msg += str(e).split(":")[-1].strip()
                msg += " for grid 'which' argument"
                axis = None

            if which is None or axis is None:
                _raise_kwarg_warning("grid", msg=msg)
            else:
                ax.grid(True, which=which, axis=axis, **grid_options)

        else:
            _raise_kwarg_warning("grid")


def _get_default_cycler():
    """Return the default cycler instance.

    Warnings
    --------
    This method is intended for internal use only.

    """
    return rc.defaultParams['axes.prop_cycle'][0]


def _get_default_colors(n_items):
    """Return a list of the colors for a given number of items.

    Warnings
    --------
    This method is intended for internal use only.

    """
    default_color_cycler = _get_default_cycler()
    points = np.linspace(0, 1, 20)
    # Get additional colors if n_items exceeds number of colors in the cycler
    if len(default_color_cycler) < n_items <= 20:
        colors = mpl.cm.get_cmap("tab20")(points)
    elif len(default_color_cycler) < n_items <= 60:
        colors = np.vstack((
            mpl.cm.get_cmap("tab20")(points),
            mpl.cm.get_cmap("tab20b")(points),
            mpl.cm.get_cmap("tab20c")(points)
        ))
    elif len(default_color_cycler) < n_items:
        # A large number of items exists, cannot use distinct tab20 colors,
        # use a different colormap altogether
        colors = mpl.cm.get_cmap("nipy_spectral")(np.linspace(0, 1, n_items))
    else:
        # Utilize the default cycler instance if it contains enough colors
        colors = [c["color"] for c in default_color_cycler]

    return list(colors)


def _raise_kwarg_warning(kwarg_name, msg=None, kwarg_prefix=None):
    """Raise a UserWarning if the kwarg input is invalid.

    Warnings
    --------
    This method is intended for internal use only.

    """
    # Prefix kwarg_name if any prefix given
    if kwarg_prefix is not None:
        kwarg_name = "_".join((kwarg_prefix, kwarg_name))

    # Format warning message
    if msg is None:
        warning_msg = "Invalid '{0}' input{1}"
        msg = "."
    else:
        warning_msg = "Invalid '{0}' input: {1}."

    # Raise warnings
    warn(warning_msg.format(kwarg_name, str(msg)))


__all__ = ()
