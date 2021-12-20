# -*- coding: utf-8 -*-
"""Contains internal functions common to various visualization functions.

Warnings
--------
The functions found in this module are NOT intended for direct use.

"""
from collections.abc import Iterable
from warnings import warn


try:
    import matplotlib as mpl
    import matplotlib.pyplot as plt
    from cycler import cycler
    from matplotlib import rcsetup as rc
except ImportError:
    mpl = None
    rc = None
    cycler = None

import numpy as np
import pandas as pd
from scipy import interpolate, stats
from six import integer_types, iteritems, iterkeys, itervalues, string_types

from mass.core.mass_configuration import MassConfiguration
from mass.util.util import apply_decimal_precision, ensure_iterable


MASSCONFIGURATION = MassConfiguration()

INSTALLED_VISUALIZATION_PACKAGES = {
    "matplotlib": bool(mpl),
    "cycler": bool(cycler),
}
r"""dict: Package names and ``bool``\ s indicating if they are installed."""

L_PAD = 0.15
"""float: Padding between legend and figure for outside legend locations."""

UL_PAD = 0.2
"""float: Additional padding for outside upper and lower legend locations."""

OUTSIDE_LEGEND_LOCATION_AND_ANCHORS = {
    "upper right outside": ("center left", (1 + L_PAD, 1 + L_PAD)),
    "upper left outside": ("center right", (0 - L_PAD, 1 + L_PAD)),
    "lower left outside": ("center right", (0 - L_PAD, 0 - L_PAD)),
    "lower right outside": ("center left", (1 + L_PAD, 0 - L_PAD)),
    "right outside": ("center left", (1 + L_PAD, 0.5)),
    "left outside": ("center right", (0 - L_PAD, 0.5)),
    "lower outside": ("upper center", (0.5, 0 - L_PAD - UL_PAD)),
    "upper outside": ("lower center", (0.5, 1 + L_PAD + UL_PAD)),
    11: (6, (1 + L_PAD, 1 + L_PAD)),
    12: (7, (0 - L_PAD, 1 + L_PAD)),
    13: (7, (0 - L_PAD, 0 - L_PAD)),
    14: (6, (1 + L_PAD, 0 - L_PAD)),
    15: (6, (1 + L_PAD, 0.5)),
    16: (7, (0 - L_PAD, 0.5)),
    17: (9, (0.5, 0 - L_PAD - UL_PAD)),
    18: (8, (0.5, 1 + L_PAD + UL_PAD)),
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
            "installed."
        )


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
    """Validate the mass solution input and return it.

    Warnings
    --------
    This method is intended for internal use only.

    """
    if mass_solution.__class__.__name__ != "MassSolution":
        raise TypeError("mass_solution must be a mass.MassSolution.")

    if not mass_solution:
        warn(
            "MassSolution '" + mass_solution.id + "' does not contain " "any solutions."
        )

    return mass_solution


def _validate_time_vector(time_vector, default_time_vector):
    """Validate the given vector of time points and return it.

    Warnings
    --------
    This method is intended for internal use only.

    """
    # Ensure time_vector is valid
    if time_vector is None:
        return np.array(sorted(default_time_vector))
    if isinstance(time_vector, Iterable) and not isinstance(time_vector, string_types):
        return np.array(sorted(time_vector))

    raise ValueError("Invalid input for `time_vector`.")


def _validate_plot_observables(mass_solution, observable, **kwargs):
    """Validate the given observables and return them and their solutions.

    Warnings
    --------
    This method is intended for internal use only.

    """
    time_vector = _validate_time_vector(kwargs.get("time_vector"), mass_solution.time)
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
        raise ValueError("`observable` must keys from the mass_solution.")

    # Turn solutions into interpolating functions if the timepoints provided
    # are not identical to those used in the simulation.
    if not isinstance(time_vector, np.ndarray):
        time_vector = np.array(time_vector)

    # Turn observable into a copy of the MassSolution containing only
    # the observable values.
    observable = mass_solution.__class__(
        id_or_model=mass_solution.id,
        solution_type=mass_solution.solution_type,
        data_dict={x: mass_solution[x] for x in observable},
        time=mass_solution.time,
        interpolate=mass_solution.interpolate,
        initial_values={
            x: mass_solution.initial_values.get(x)
            for x in observable
            if mass_solution.initial_values.get(x, None) is not None
        },
    )

    # Change the time points and solutions if the time_vector has changed
    if not np.array_equal(observable.time, time_vector):
        observable.interpolate = True
        observable.time = time_vector

    # Turn interpolating functions into solutions.
    observable.interpolate = False

    if kwargs.get("deviation"):
        observable = _calculate_deviation_solutions(observable, **kwargs)

    return observable


def _validate_ensemble_plot_observables(mass_solution_list, observable, **kwargs):
    """Validate the given observables and observable solutions in DataFrames.

    Warnings
    --------
    This method is intended for internal use only.

    """

    # Get time vector with the most points
    time_vector_lens = [len(sol.time) for sol in mass_solution_list]
    idx = time_vector_lens.index(max(time_vector_lens))
    # Use time vector with most points as default time vector for all solutions
    time_vector = _validate_time_vector(
        kwargs.get("time_vector"), mass_solution_list[idx].time
    )

    # Replace mass objects with their identifiers.
    observable = [getattr(x, "id", x) for x in ensure_iterable(observable)]

    # Raise error if solutions are missing observables
    missing_observables = [
        mass_solution_list[idx]
        for idx, x in enumerate(
            [set(observable).issubset(set(iterkeys(sol))) for sol in mass_solution_list]
        )
        if not x
    ]
    if missing_observables:
        raise ValueError(
            "The following MassSolution objects are missing the specified "
            "observables:\n{0!r}".format([s.id for s in missing_observables])
        )

    # Turn observable into a copy of the MassSolution containing only
    # the observable values as interpolating functions
    mass_solution_list = [
        sol.__class__(
            id_or_model=sol.id,
            solution_type=sol.solution_type,
            data_dict={x: sol[x] for x in observable},
            time=sol.time,
            interpolate=True,
            initial_values={
                x: sol.initial_values.get(x)
                for x in observable
                if sol.initial_values.get(x) is not None
            },
        )
        for sol in mass_solution_list
    ]

    for sol in mass_solution_list:
        sol.time = time_vector
        sol.interpolate = False
        if kwargs.get("deviation"):
            sol = _calculate_deviation_solutions(sol, **kwargs)

    # Make DataFrames for each observable solution with timepoints as columns
    observable_dataframes = {
        sol_key: pd.DataFrame(
            data=[sol[sol_key] for sol in mass_solution_list if sol_key in sol],
            columns=time_vector,
        )
        for sol_key in observable
    }

    return observable_dataframes, time_vector


def _validate_interval_type(interval_type):
    """Validate the given `interval_type` value.

    Warnings
    --------
    This method is intended for internal use only.

    """
    if interval_type in [None, "range"]:
        return interval_type

    if interval_type.startswith("CI="):
        interval_type = float(interval_type[3:]) / 100
        if 0 <= interval_type <= 1:
            return interval_type

    raise ValueError("Invalid `interval_type`")


def _calculate_confidence_interval(ci_value, sol_df, CI_distribution):
    """Calculate the confidence interval bounds.

    Warnings
    --------
    This method is intended for internal use only.

    """
    n = len(sol_df)

    if CI_distribution not in ["t", "z"]:
        warn("Unrecognized `CI_distribution`, defaulting to t-distribution")
        CI_distribution = "t"

    options = {}
    if CI_distribution.startswith("z"):
        interval_func = stats.norm.interval
        ddof = 0
    else:
        interval_func = stats.t.interval
        ddof = 1
        options["df"] = n - ddof

    mean = sol_df.mean(axis=0)
    options["scale"] = sol_df.std(axis=0, ddof=ddof) / np.sqrt(n)
    options["loc"] = mean

    lower_bound, upper_bound = interval_func(ci_value, **options)

    return lower_bound, upper_bound, mean


def _calculate_deviation_solutions(observable, initial_values=None, **kwargs):
    """Calculate the new deviated solution using the deviation kwargs.

    Warnings
    --------
    This method is intended for internal use only.

    """
    if kwargs.get("deviation_normalization") in ["range", "initial value"]:
        deviation_normalization = kwargs.get("deviation_normalization")
    else:
        raise ValueError("Unrecognized `deviation_normalization` input")

    if initial_values is None:
        initial_values = observable.initial_values

    for label, sol in iteritems(observable.copy()):
        if kwargs.get("deviation_zero_centered"):
            # Shift solutions to be centered around 0.
            sol = sol - initial_values[label]

        if deviation_normalization == "initial value":
            # Avoid plotting solutions that produce a divide by 0 error
            if initial_values[label] == 0:
                warn(
                    "Initial value for '{0}' is 0, therefore cannot plot "
                    "deviation from initial value.".format(label)
                )
                del observable[label]
                continue
            sol = sol / initial_values[label]

        if deviation_normalization == "range":
            # Avoid plotting solutions that produce a divide by 0 error
            value_range = apply_decimal_precision(
                sol.max() - sol.min(), MASSCONFIGURATION.decimal_precision
            )
            if value_range == 0:
                sol.fill(0)
            else:
                sol = sol / value_range
        # Overwrite old solution with new one
        observable[label] = sol

    return observable


def _validate_legend_input_fmt(legend, observable):
    """Validate whether the given ``legend`` can be interpreted.

    Warnings
    --------
    This method is intended for internal use only.

    """
    # Ensure legend is iterable and get default labels and legend location
    legend = ensure_iterable(legend)

    # Get lengths of observables and legend arguments for comparisions.
    n_obs, n_leg = len(observable), len(legend)

    # Get the possible location value in the input.
    possible_loc = legend[-1]

    # Legend location is an integer
    if isinstance(possible_loc, (integer_types, float)):
        if legend[0] != possible_loc and (
            isinstance(legend[0], list) or n_obs <= n_leg <= n_obs + 1
        ):
            items = ensure_iterable(legend[0]), possible_loc
        elif legend[0] == possible_loc:
            items = list(observable), possible_loc
        else:
            items = None, possible_loc

    # Legend format includes iterable of labels and a location as a str
    elif n_leg == 2 and (
        not isinstance(legend[0], string_types)
        and isinstance(possible_loc, string_types)
    ):
        items = legend[0], legend[-1]
    # Only one observable value or legend value
    elif 1 in [n_obs, n_leg]:
        # Legend has one label and one location
        if n_leg == n_obs + 1:
            items = ensure_iterable(legend[0]), possible_loc
        elif n_leg in [n_obs, 1]:
            items = list(observable), possible_loc
    # Either only labels provided or bad label input and location
    elif n_leg == n_obs:
        items = legend, possible_loc
    # Bad legend input, return None to set off warnings
    else:
        items = None, None

    return items[0], items[1]


def _validate_annotate_time_points_input(*args):
    """Validate the inputs for kwargs of ``annotate_time_points`` and return.

    Warnings
    --------
    This method is intended for internal use only.

    """
    # Get args and determine minimum and maximum time points
    time_vector, time_points, colors, markers, marker_sizes = args
    t_min = min(time_vector)
    t_max = max(time_vector)

    # Get the start and finish points
    if isinstance(time_points, string_types) and time_points.lower() == "endpoints":
        time_points = [t_min, t_max]
        default_markers = ["o", "D"]
        default_colors = ["r", "b"]
    else:
        # Ensure the time points are iterable
        time_points = ensure_iterable(time_points)
        default_markers = ["o"] * len(time_points)
        default_colors = ["k"] * len(time_points)

    # Get the default marker sizes
    default_sizes = [rc.defaultParams["lines.markersize"][0]] * len(default_markers)

    # Ensure time points are in the time vector
    outside_t_range = []
    for t in time_points:
        if t_min <= t <= t_max:
            continue
        else:
            outside_t_range += [t]

    # If outside the time range, no need to check colors, markers or sizes
    if outside_t_range:
        outside_t_range += [[t_min, t_max]]
    else:
        # Valdiate markers
        if markers is not None and len(time_points) < len(markers):
            _raise_kwarg_warning(
                "marker",
                msg=", therefore utilizing default `markers` values instead",
                kwarg_prefix="annotate_time_points",
            )
            markers = None

        # Assign value dependning on whether validation was successful
        markers = default_markers if not markers else markers

        # Validate marker sizes
        marker_sizes = _validate_time_points_marker_properties(
            "markersize", marker_sizes, len(markers)
        )
        # Assign value dependning on whether validation was successful
        marker_sizes = default_sizes if not marker_sizes else marker_sizes

    # Validate time points colors
    colors = _validate_time_points_marker_properties("color", colors, len(time_points))
    if colors is None:
        colors = default_colors

    return ([time_points, colors, markers, marker_sizes], outside_t_range)


def _validate_time_points_marker_properties(marker_prop, values, num_expected):
    """Validate marker properties such as markersize and color.

    Warnings
    --------
    This method is intended for internal use only.

    """
    # Allow 'None' values
    if values is not None:
        values = ensure_iterable(values)
        try:
            # Check whether the values are valid
            # Ensure that the number of values is equal to the number expected
            if len(values) == 1:
                values = values * num_expected
            elif len(values) != num_expected:
                raise ValueError(
                    "wrong number of values for '{0}' provided".format(marker_prop)
                )
        except ValueError as e:
            # Raise warning and set marker defaults if invalid input.
            msg = str(e).lower()
            msg += ", therefore utilizing default {0} values instead".format(
                marker_prop
            )
            _raise_kwarg_warning(
                marker_prop, msg=msg, kwarg_prefix="annotate_time_points"
            )
            # Returning None causes defaults to be used.
            values = None

    return values


def _validate_tile_placement(tile_placement, prefix=None):
    """Validate the given ``tile_placement`` input.

    Warnings
    --------
    This method is intended for internal use only.

    """
    # Ensure tile placement value is valid
    valid = {"all", "lower", "upper"}
    if tile_placement.lower() not in valid:
        arg_name = "tile_placement"
        if prefix is not None:
            arg_name = "_".join(("plot", arg_name))
        raise ValueError(
            "Invalid `{0}` input '{1}'. Can only be one of the following: "
            "{2}.".format(arg_name, tile_placement, str(valid))
        )

    return tile_placement.lower()


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
        "loglog": ax.loglog,
    }

    # Ensure plotting function is valid
    if (
        plot_function_str not in plotting_functions_dict
        or plot_function_str not in valid
    ):
        raise ValueError(
            "'{0}' not a valid plotting option. Must be one of the "
            "following options: '{1}'.".format(plot_function_str, valid)
        )

    return plotting_functions_dict[plot_function_str]


def _get_legend_args(ax, legend, observable, **kwargs):
    """Validate whether the given ``legend`` can be interpreted.

    Warnings
    --------
    This method is intended for internal use only.

    """
    # Validate legend input format
    legend_labels, legend_loc = _validate_legend_input_fmt(legend, observable)
    if legend_labels is None:
        _raise_kwarg_warning(
            "legend", msg="therefore utilizing keys from the MassSolution instead"
        )
        legend_labels = list(observable)

    if legend_loc is None:
        _raise_kwarg_warning(
            "loc",
            msg="therefore utilizing default legend `loc` instead",
            kwarg_prefix="legend",
        )
        legend_loc = rc.defaultParams["legend.loc"][0]

    if legend_loc in OUTSIDE_LEGEND_LOCATION_AND_ANCHORS:
        legend_loc, anchor = OUTSIDE_LEGEND_LOCATION_AND_ANCHORS[legend_loc]
    else:
        anchor = None

    # Get the current number of lines on the axes
    n_current = len(_get_ax_current(ax))
    # Get the legend properties
    ncols = kwargs.get("legend_ncol")
    # Use a default number of columns based on the total number of items
    if ncols is None:
        ncols = int(np.ceil(np.sqrt(len(observable) + n_current) / 3))

    legend_kwargs = {"loc": legend_loc, "bbox_to_anchor": anchor, "ncol": ncols}

    return legend_labels, legend_kwargs


def _map_labels_to_solutions(obs, labels, initial_values=None):
    """Map the legend labels to the observable solutions and return.

    Warnings
    --------
    This method is intended for internal use only.

    """
    labels = [label.lstrip("_") for label in labels]
    if isinstance(obs, pd.DataFrame):
        obs.index = labels
    else:
        # Map solutions to labels
        for old, new in zip(list(obs), list(labels)):
            obs[new] = obs.pop(old)
            try:
                obs._initial_values[new] = obs._initial_values.pop(old)
            except (AttributeError, KeyError):
                pass
    return obs


def _get_line_property_cycler(n_current, n_new, kwarg_prefix=None, **kwargs):
    """Get the line propertie and return a cycler to be set.

    Warnings
    --------
    This method is intended for internal use only.

    """
    # Get the prop_cycler
    prop_cycler = kwargs.get("prop_cycler")
    n_total = n_current + n_new
    if prop_cycler is None:
        cycler_kwargs_values = {}
        # Get current number of lines on the plot, excluding time points.
        default_cycler_values = {
            "color": _get_default_colors(n_total)[n_current:n_total],
            "linestyle": [rc.defaultParams["lines.linestyle"][0]] * n_new,
            "linewidth": [rc.defaultParams["lines.linewidth"][0]] * n_new,
            "marker": [rc.defaultParams["lines.marker"][0]] * n_new,
            "markersize": [rc.defaultParams["lines.markersize"][0]] * n_new,
        }

        for k in default_cycler_values:
            values = ensure_iterable(kwargs.get(k))
            # Don't add to cycler instance if it was not provided.
            if not values:
                continue

            msg = ""

            # If only one value given, apply it to all lines.
            if len(values) == 1 and n_new != 1:
                values = values * n_new
            # Raise warning if the number of given values does not match the
            # number of new items. Do not need to include this message
            # if invalid kwarg values were found.
            if len(values) != n_new and not msg:
                msg += "Wrong number of values for `{0}` provided".format(k)

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

    return prop_cycler


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
                warn(
                    "Could not set `{0}` due to the following: '{1}'.".format(
                        label_type, e
                    )
                )


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
                warn(
                    "Could not set `{0}` due to the following: '{1}'".format(
                        limit_type, e
                    )
                )


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
            if v is not None:
                grid_options[k] = v

        if isinstance(grid, bool):
            if grid:
                ax.grid(grid, **grid_options)
            else:
                ax.grid(None)

        elif isinstance(grid, Iterable) and len(grid) == 2:
            which, axis = grid

            ax.grid(True, which=which, axis=axis, **grid_options)

        else:
            _raise_kwarg_warning("grid")


def _set_axes_margins(ax, x_default=None, y_default=None, tile=False, **kwargs):
    """Set the axes margins.

    Warnings
    --------
    This method is intended for internal use only.

    """
    if tile:
        prefix = "tile_"
    else:
        prefix = ""
    for arg, margin_default in zip(["x", "y"], [x_default, y_default]):
        margin_arg = arg + "margin"
        # Validate margin value
        margin_value = kwargs.get(prefix + margin_arg, None)
        # Use default value if None
        if margin_value is None:
            # Use default from rcsetup
            if margin_default is None:
                margin_default = rc.defaultParams["axes." + margin_arg][0]
            margin_value = margin_default
        # Set the margin value
        getattr(ax, "set_" + margin_arg)(margin_value)


def _set_annotated_time_points(
    ax, observable, type_of_plot=None, first_legend=None, time_range=None, **kwargs
):
    """Set the given ``time_points`` and kwargs.

    Warnings
    --------
    This method is intended for internal use only.

    """
    items = _validate_annotate_time_points_input(
        time_range,
        kwargs.get("annotate_time_points"),
        kwargs.get("annotate_time_points_color"),
        kwargs.get("annotate_time_points_marker"),
        kwargs.get("annotate_time_points_markersize"),
    )

    if items[-1]:
        # Invalid time points, end function early
        _raise_kwarg_warning(
            "annotate_time_points",
            msg="points '{1}' outside of time range {0}".format(
                items[-1].pop(-1), items[-1]
            ),
        )

        return ax

    time_points, colors, markers, marker_sizes = items[0]

    plot_function = _get_plotting_function(
        ax,
        plot_function_str=kwargs.get("plot_function"),
        valid={"plot", "semilogx", "semilogy", "loglog"},
    )

    if hasattr(observable, "time"):
        # Change the time points of the MassSolution to make the MassSolution
        # carry only solutions for time points to be plotted.
        observable.time = time_points
    else:
        # Interpolate mean solution at time points to annotate
        for label, df in iteritems(observable):
            func = interpolate.interp1d(
                time_range, df.mean(axis=0), kind="cubic", fill_value="extrapolate"
            )
            observable[label] = np.array(func(time_points))

    lines = _get_ax_current(ax)
    if type_of_plot == "phase_portrait":
        lines = lines[len(lines) - 1 :]
    else:
        lines = lines[len(lines) - len(observable) :]

    for line in lines:
        # Set up the prop_cycler arguments for the time points
        if colors is None:
            colors = [line.get_color()] * len(time_points)
        # Make and set the prop_cycler for the time points
        ax.set_prop_cycle(
            cycler(
                **{
                    "color": colors,
                    "linestyle": [" "] * len(time_points),
                    "marker": markers,
                    "markersize": marker_sizes,
                }
            ),
        )

        for i, t in enumerate(time_points):
            label = "t=" + str(t)
            # Check if phase portrait or time profile
            if type_of_plot == "phase_portrait":
                # Handle as a phase portrait, get solutions and plot
                items = _group_xy_items(observable, itervalues)
                x, y = items[0][i], items[1][i]
            else:
                # Handle as a time profile, get solutions and plot
                items = observable[line.get_label()]
                x, y = t, items[i]

            plot_function(
                x, y, label=label, zorder=kwargs.get("annotate_time_points_zorder")
            )
            if kwargs.get("annotate_time_points_labels"):
                ax.annotate(
                    label,
                    xy=(x, y),
                    xytext=(10, 10),
                    textcoords="offset pixels",
                    horizontalalignment="left",
                    verticalalignment="top",
                )

    # Add the legend to the axes.
    if kwargs.get("annotate_time_points_legend") is not None:
        items = _get_annotated_time_points_legend_args(
            ax,
            desired_loc=kwargs.get("annotate_time_points_legend"),
            taken_loc=first_legend[1],
        )

        ax = _set_additional_legend_box(ax, items, first_legend[0])

    return ax


def _check_second_legend_location(desired_loc, taken_loc):
    """Check the legend location and determine if it is being occupied already.

    Warnings
    --------
    This method is intended for internal use only.

    """
    # Get the taken location (check if it is an outside location too)
    if taken_loc:
        bbox_to_anchor = taken_loc["bbox_to_anchor"]
        taken_loc = taken_loc["loc"]
        for k, v in iteritems(OUTSIDE_LEGEND_LOCATION_AND_ANCHORS):
            if v == (taken_loc, bbox_to_anchor):
                taken_loc = k
                break

    # Make sure desired location is not already taken
    if desired_loc is not None and desired_loc == taken_loc:
        msg = " location already used, utilizing default value instead"
        _raise_kwarg_warning("legend_loc", msg=msg, kwarg_prefix="annotate_time_points")
        desired_loc = None

    return desired_loc, taken_loc


def _get_annotated_time_points_legend_args(ax, desired_loc=None, taken_loc=None):
    """Get the arguments needed for a legend box of the annotated time points.

    Warnings
    --------
    This method is intended for internal use only.

    """
    points, labels = _get_handles_and_labels(ax, time_points=True)
    labels_and_points = {}
    for point, label in zip(points, labels):
        # Only make a legend containing the latest set of time points if there
        # are duplicates.
        labels_and_points[label] = point
    # Get points and labels as list
    points = list(itervalues(labels_and_points))
    labels = list(iterkeys(labels_and_points))

    desired_loc, taken_loc = _check_second_legend_location(desired_loc, taken_loc)
    anchor = None
    # Get kwargs for legend location
    if desired_loc in OUTSIDE_LEGEND_LOCATION_AND_ANCHORS:
        desired_loc, anchor = OUTSIDE_LEGEND_LOCATION_AND_ANCHORS[desired_loc]

    if desired_loc is None and anchor is None:
        legend_kwargs = {}
    else:
        legend_kwargs = {"loc": desired_loc, "bbox_to_anchor": anchor}

    return points, labels, legend_kwargs


def _set_additional_legend_box(ax, legend_args, first_legend=None):
    """Set the legend box of the annotated time points.

    Warnings
    --------
    This method is intended for internal use only.

    """
    points, labels, legend_kwargs = legend_args
    # Remove old legends
    for old_legend in ax.get_children():
        if old_legend.__class__.__name__ == "Legend":
            old_legend.remove()

    # Set legend, adding the first legend back if needed
    ax.legend(handles=points, labels=labels, **legend_kwargs)
    if first_legend:
        ax.add_artist(first_legend)

    return ax


def _get_default_cycler():
    """Return the default cycler instance.

    Warnings
    --------
    This method is intended for internal use only.

    """
    return rc.defaultParams["axes.prop_cycle"][0]


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
        colors = np.vstack(
            (
                mpl.cm.get_cmap("tab20")(points),
                mpl.cm.get_cmap("tab20b")(points),
                mpl.cm.get_cmap("tab20c")(points),
            )
        )
    elif len(default_color_cycler) < n_items:
        # A large number of items exists, cannot use distinct tab20 colors,
        # use a different colormap altogether
        colors = mpl.cm.get_cmap("nipy_spectral")(np.linspace(0, 1, n_items))
    else:
        # Utilize the default cycler instance if it contains enough colors
        colors = [c["color"] for c in default_color_cycler]

    return list(colors)


def _get_ax_current(ax, time_points=False):
    """Return current lines or time points on the axis.

    Warnings
    --------
    This method is intended for internal use only.

    """
    lines_with_labels = [
        label for label in ax.get_lines() if "_line" not in label.get_label()
    ]
    if time_points:
        return [
            label for label in lines_with_labels if label.get_label().startswith("t=")
        ]

    lines_with_labels = [
        label
        for label in lines_with_labels
        if not (
            label.get_label().endswith("_lb")
            or label.get_label().endswith("_ub")
            or label.get_label().startswith("t=")
        )
    ]

    return lines_with_labels


def _get_handles_and_labels(ax, time_points=False, bounds=False):
    """Return legend handles and labels for the axis.

    Warnings
    --------
    This method is intended for internal use only.

    """
    lines = _get_ax_current(ax, time_points=time_points)
    labels = [line.get_label() for line in lines]
    return (lines, labels)


def _group_xy_items(xy, iterfunction):
    """Seperate items in xy, grouping either the keys or the values.

    Warnings
    --------
    This method is intended for internal use only.

    """
    return [i for i in list(iterfunction(xy))]


def _raise_kwarg_warning(kwarg_name, msg=None, kwarg_prefix=None):
    """Raise a UserWarning if the kwarg input is invalid.

    Warnings
    --------
    This method is intended for internal use only.

    """
    # Prefix kwarg_name if any prefix given
    if kwarg_prefix:
        kwarg_name = "_".join((kwarg_prefix, kwarg_name))

    # Format warning message
    if msg is None:
        warning_msg = "Invalid `{0}` input{1}"
        msg = "."
    else:
        warning_msg = "Invalid `{0}` input: {1}."

    # Raise warnings
    warn(warning_msg.format(kwarg_name, str(msg)))


def _get_values_as_series(obj, compare, name=None):
    """Get values from ``obj`` as a pandas.Series.

    Warnings
    --------
    This method is intended for internal use only.

    """
    cls_ = obj.__class__

    def _validate_compare_type(compare, valid):
        """Ensure ``compare`` is valid."""
        if compare is None or compare not in valid:
            raise ValueError(
                "Invalid `compare` '{0}' given, cannot access "
                "'{1}' values.".format(compare, cls_.__name__)
            )

    if cls_.__name__ == "Series":
        series = obj.copy()

    elif "MassModel" in [cls_.__name__, cls_.__base__.__name__]:
        _validate_compare_type(compare, ["concentrations", "Keqs", "fluxes", "kfs"])
        values = getattr(
            obj,
            {
                "concentrations": "initial_conditions",
                "fluxes": "steady_state_fluxes",
                "Keqs": "parameters",
                "kfs": "parameters",
            }[compare],
        )
        if compare in ["Keqs", "kfs"]:
            values = values[compare[:-1]]
        else:
            values = {x.id: v for x, v in iteritems(values)}
        series = pd.Series(values)

    elif cls_.__name__ == "ConcSolution":
        _validate_compare_type(compare, ["concentrations", "Keqs"])
        series = getattr(obj, compare).copy()

    elif cls_.__name__ == "Solution":
        _validate_compare_type(compare, ["fluxes"])
        series = getattr(obj, compare).copy()
    else:
        raise TypeError("Unrecognized object type.")

    if name:
        series.name = name

    return series


def _get_dataframe_of_observables(x, y, compare, observable):
    """Concat ``x`` and ``y`` into DataFrame containing ``observable`` values.

    Warnings
    --------
    This method is intended for internal use only.

    """
    # Concat x and y together
    xy = pd.concat([x, y], axis=1, sort=False)
    # Filter for plot observables
    if observable is not None:
        if compare in ["Keqs", "kfs"]:
            p_type = compare[:-1]
            observable = [
                getattr(x, p_type + "_str", "_".join((p_type, x)))
                for x in ensure_iterable(observable)
            ]
        else:
            observable = [getattr(x, "id", x) for x in ensure_iterable(observable)]
        if set(observable).difference(xy.index):
            raise ValueError(
                "Invalid `observable` values: '{0}'".format(
                    set(observable).difference(xy.index)
                )
            )
        else:
            xy = xy.loc[observable]

    # Check for NA values, raise warning
    diff = set(xy.index).symmetric_difference(set(xy.dropna().index))
    if diff:
        warn("Ignoring {0}, only in one set of given values".format(diff))
        xy.dropna(inplace=True)

    return xy


__all__ = ()
