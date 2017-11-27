# -*- coding: utf-8 -*-

# Compatibility with Python 2.7
from __future__ import absolute_import

# Import necesary packages
import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.patches as patches

from pandas.compat import range, lrange, zip
from pandas.io.formats.printing import pprint_thing
from pandas.plotting._style import _get_standard_colors
from pandas.plotting._tools import _subplots, _set_ticks_props

from scipy.interpolate import interp1d
from math import inf
from six import iterkeys, itervalues
from cycler import cycler

# from cobra
from cobra import DictList

# from mass
from mass.core.massmetabolite import MassMetabolite
from mass.core.massreaction import MassReaction
from mass.core.massmodel import MassModel


def plot_simulation(time, solution_profile, default_fontsize=15, **kwargs):
    """Generates a plot of the data in the solution_profile over time.

    ``kwargs`` are passed on to various matplotlib methods. 
    See get_default_options() for a full description.

    Parameters
    ----------
    time: list or tuple or np.array
        An array containing the time points over with the system was simulated
    solution_profile : np.array
        An array containing the simulated results for either concentration
        or flux
    default_fontsize: int
        The value of the default fontsize for much of the plot text

    Returns
    -------
    plt.gcf()
        A reference to the current figure instance. Shows plot when returned.
        Can be used to modify plot after initial generation.

    See Also:
    ---------
    get_default_options()
    """

    # Process time vector into numpy array
    time = _gen_np_time_vector(time)

    # Step 0: Get options if provided, else use defaults
    options = _get_options(**kwargs)
    options = _process_title_options(default_fontsize, **options)

    # Step 1: Index by time and get time range
    start, final = _get_time_range(**options)
    time = pd.DataFrame(time).loc[start:final].values

    # Step 2: Generate conc/flux arrays for plotting
    plt_args = _gen_plot_args(time, solution_profile)

    # Step 3: Make plot using options and vectors provided
    style, xgrid, ygrid = _set_style(**options)

    with matplotlib.style.context(style):
        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.plot(*plt_args)

        # Step 4: Add remaining plotting options
        plt.rc("axes", prop_cycle=(cycler("color", _get_colormap())))
        #FIXME
        #_add_custom_linecolors(fig, ax, legend_ids, **options)
        _add_plot_range(ax, **options)

        _plot_title_options(**options)
        _plot_figsize(fig, **options)
        #FIXME
        #_plot_legend(legend_ids, default_fontsize, **options)

        _set_log_scale(ax, default=True, **options)
        _plot_gridlines(ax, xgrid, ygrid)
        _option_savefig(**options)

        # Step 5: Return plot/show plot figure
        plt.show()
        return plt.gcf()




def get_default_options():
    """Returns the current default options as a dictionary. These options are
    used when ``kwargs`` are not provided for plot methods (and thus optional)
    Below is the list of default options for the plotting functions.
    These options correspond to the ``kwargs``

    ``kwargs``:
    -----------
    plot_function: str
        A string that controls the scaling of the plot.
        "LogLogPlot" sets the x and y axes to log scale
        "LogLinearPlot", sets only the x axis to log scale (y axis is linear)
        "LinearLogPlot", sets only the y axis to log scale (x axis is linear)
        "LinearPlot", sets the x and y axes to linear scale
    title: str or tuple
        The title of the figure. If tuple, the 2nd value determines fontsize
    xlabel: str or tuple
        The title of the x-axis. If tuple, the 2nd value determines fontsize
    ylabel: str or tuple
        The title of the y-axis. If tuple, the 2nd value determines fontsize
    plot_legend: bool
        Whether of not to display the legend
    observable: list
        A list containing mass.MassMetabolite or str elements
        This list is used to filter plot_simulation to only plot specific elements
    set_xlim: tuple
        A tuple containing 2 float elements
        This tuple is used to manually restrict the plotting range for the x-axis
    set_ylim: tuple
        A tuple containing 2 float elements
        This tuple is used to manually restrict the plotting range for the y-axis
    trange: tuple
        A tuple containing 2 float elements
        This tuple is used to restrict the range over which plotlines are drawn
    truncate: bool
        Whether or not to truncate floating time point labels on the plot
    linecolor: dict
        A dict containing mass.MassMetabolite or mass.MassReaction elements 
        as the key and matplotlib-compliant str elements as the values.
        This allows various metabolite/flux plotlines to be colored a specific
        value.

        See matplotlib's set_color() method for additional details
    figsize: tuple
        A tuple containing 2 float elements
        This tuple is used to manually change the figsize of the plot (in inches)
    style: str
        Used to determine the style of the plot. 
        Accepted str values are: 
            bmh, classic, dark_background, fivethirtyeight, ggplot, grayscale,
            seaborn-bright, seaborn-colorblind, seaborn-dark-palette, 
            seabon-dark, seabon-darkgrid, seabon-deep, seabon-muted, 
            seabon-notebook, seabon-paper, seabon-pastel, seabon-poster, 
            seabon-talk, seabon-seabon-ticks, seabon-white, seabon-whitegrid,
            seabon, _classic-test, default

        See matplotlib.style.use() for additional details
    grid: bool
        Used to turn on/off the gridlines in the plot.
        True turns gridlines on, False turns gridlines off
    savefig: dict
        A dict containing at least "fname" and "dpi" as keys
        Saves the file at fname with dpi of dpi

        See matplotlib's savefig() method for additional details

    See Also:
    ---------
    set_default_options(**custom)
    restore_default_options()
    """
    return default_options



def set_default_options(**custom):
    """Allows user to change the global default options (provided if 
    ``kwargs`` are not not specified)

    Parameters:
    -----------
    **custom: dict
        A dictionary of ``kwargs`` with options as the keys (str) and the
        desired option defaults as the values

    See Also:
    ---------
    get_default_options()
    restore_default_options()
    """
    default_options.update(custom)

def restore_default_options():
    """Restores default_options to their original values

    See Also:
    ---------
    get_default_options()
    set_default_options(**custom)
    restore_default_options()
    """
    default_options = _base_default_options

def get_tiled_options():
    """Returns the current default options as a dictionary. These options are
    used when ``kwds`` are not provided for plot methods (and thus optional)
    Below is the list of default options for the plotting functions.
    These options correspond to the ``kwds``

    ``kwds``:
    -----------
    trange: tuple
        A tuple containing 2 float elements
        This tuple is used to restrict the range over which plotlines are drawn
    truncate: bool
        Whether or not to truncate floating time point labels on the plot
    linecolor: dict
        A dict containing mass.MassMetabolite or mass.MassReaction elements 
        as the key and matplotlib-compliant str elements as the values.
        This allows various metabolite/flux plotlines to be colored a specific
        value.

        See matplotlib's set_color() method for additional details
    style: str
        Used to determine the style of the plot. 
        Accepted str values are: 
            bmh, classic, dark_background, fivethirtyeight, ggplot, grayscale,
            seaborn-bright, seaborn-colorblind, seaborn-dark-palette, 
            seabon-dark, seabon-darkgrid, seabon-deep, seabon-muted, 
            seabon-notebook, seabon-paper, seabon-pastel, seabon-poster, 
            seabon-talk, seabon-seabon-ticks, seabon-white, seabon-whitegrid,
            seabon, _classic-test, default

        See matplotlib.style.use() for additional details
    grid: bool
        Used to turn on/off the gridlines in the plot.
        True turns gridlines on, False turns gridlines off
    display_axes: bool
        Used to turn on/off the axes in the plot.
        True gives all the subplots their own axes with labels
        False removes them from the subplots
    savefig: dict
        A dict containing at least "fname" and "dpi" as keys
        Saves the file at fname with dpi of dpi

        See matplotlib's savefig() method for additional details

    See Also:
    ---------
    set_tiled_options(**custom)
    restore_tiled_options()
    """
    return tiled_default_options

def set_tiled_options(**custom):
    """Allows user to change the global default options (provided if 
    ``kwds`` are not not specified)

    Parameters:
    -----------
    **custom: dict
        A dictionary of ``kwargs`` with options as the keys (str) and the
        desired option defaults as the values

    See Also:
    ---------
    get_tiled_options()
    restore_tiled_options()
    """
    tiled_default_options.update(custom)

def restore_tiled_options():
    """Restores default_options to their original values

    See Also:
    ---------
    get_tiled_options()
    set_tiled_options(**custom)
    restore_tiled_options()
    """
    tiled_default_options = _base_tiled_default_options



# Internal Methods - Global
def _sort_df(sol_df):
    ds = DictList(sol_df.columns)
    ls = []

    for col in sol_df.columns:
        ls.append(col.id)
    ls = sorted(ls)

    ss = DictList()
    for i in range(len(ls)):
        ss.append(ds.get_by_id(ls[i]))

    return sol_df.reindex_axis(ss, axis=1)

def _gen_np_time_vector(time):
    if isinstance(time, list) or isinstance(time, tuple):
        return np.array(time)
    elif isinstance(time, np.ndarray):
        return time
    else:
        msg = "time must be of type list, tuple, or numpy.ndarray"
        raise TypeError(msg)

def _gen_plot_args(time, solution_profile):
    # Convert list of time points into list of lists of len(solution_profile)
    t_vec = [time] * len(solution_profile)

    # Convert dict of interpolating functions into list of lists
    v_list = list(solution_profile.values())
    v_list = [v_list[i](time) for i in range(len(v_list))]

    # Generate list of interleaving lists (for Axes.plot method)
    plt_args = [val for pair in zip(t_vec, v_list) for val in pair]

    return plt_args

def _get_options(tiled=False, **kwargs):
    options = {}
    default = tiled_default_options if tiled else default_options

    for key in default:
        options[key] = kwargs[key] if key in kwargs else default[key]
    
    return options

def _process_title_options(default_fontsize, **options):
    dict_of_titles = {
        "title": plt.title, 
        "xlabel": plt.xlabel, 
        "ylabel": plt.ylabel
    }
    for name in dict_of_titles:
        if isinstance(options[name], str):
            options.update({name: (options[name], default_fontsize)})

    return options

def _get_time_range(**options):
    start = options["trange"][0] if options["trange"][0] is not None else None
    final = options["trange"][1] if options["trange"][1] is not None else None

    return start, final

def _add_plot_range(axes, **options):
    if options["set_xlim"] is not None:
        axes.set_xlim(options["set_xlim"][0], options["set_xlim"][1])
        axes.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
    if options["set_ylim"] is not None:
        axes.set_ylim(options["set_ylim"][0], options["set_ylim"][1])
        axes.yaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())

def _plot_title_options(**options):
    dict_of_titles = {
        "title": plt.title, 
        "xlabel": plt.xlabel, 
        "ylabel": plt.ylabel
    }
    for name in dict_of_titles:
        if options[name][0] is not None:
            dict_of_titles[name](**{
                "s": options[name][0],
                "fontsize": options[name][1]
            })

def _plot_figsize(fig, **options):
    if options["figsize"] != default_options["figsize"]:
        fig.set_size_inches(options["figsize"][0], options["figsize"][1])

def _set_style(**options):    
    style = options["style"] if options["style"] is not None else "default"
    grid  = options["grid"]
    xgrid, ygrid = (None, None)

    if isinstance(grid, tuple):
        xgrid = grid[0]
        ygrid = grid[1]
    elif (grid is None) and (style is "default"):
        xgrid = True
        ygrid = True
    elif (grid is None) and (style is not "default"):
        xgrid = False
        ygrid = False
    elif isinstance(grid, bool):
        xgrid = grid
        ygrid = grid

    return style, xgrid, ygrid

def _plot_gridlines(ax, xgrid, ygrid):
    if xgrid:
        ax.xaxis.grid(True, linestyle="--")
    if ygrid:
        ax.yaxis.grid(True, linestyle="--")

def _option_savefig(**options):
    if options["savefig"] != default_options["savefig"]:
        fig.savefig(**options["savefig"])



# Internal Methods - plot_simulation
def _get_conc_flux_array(sol_df, start, final, **options):
    dictlist_metabs = DictList(sol_df.columns.tolist())
    legend_ids = None

    if options["observable"] is None or options["observable"] == []:
        df_conc_flux = sol_df
        legend_ids = [item.id for item in sol_df.columns.tolist()]
    else:
        if isinstance(options["observable"], MassMetabolite):
            options["observable"] = [options["observable"]]
        elif isinstance(options["observable"], str):
            options["observable"] = [dictlist_metabs.get_by_id(
                options["observable"])]
        elif isinstance(options["observable"], list):
            if isinstance(options["observable"][0], str):
                new_list = []
                for string in options["observable"]:
                    new_list.append(dictlist_metabs.get_by_id(string))
                options["observable"] = new_list
            elif isinstance(options["observable"][0], MassMetabolite):
                pass
            else:
                raise TypeError("Expected MassMetabolite, string, list of "\
                                "MassMetabolites, or list of strings")

        df_conc_flux = sol_df[options["observable"]]
        legend_ids = [item.id for item in options["observable"]]

    return df_conc_flux, legend_ids

def _get_colormap():
    cm = plt.cm.get_cmap("tab20")
    colors1 = [cm.colors[i] for i in range(len(cm.colors))]
    cm = plt.cm.get_cmap('tab20b')
    colors2 = [cm.colors[i] for i in range(len(cm.colors))]
    cm = plt.cm.get_cmap('tab20c')
    colors3 = [cm.colors[i] for i in range(len(cm.colors))]
    colors = colors1 + colors2 + colors3

    return colors

def _add_custom_linecolors(fig, axes, legend_ids, **options):
    if options["linecolor"] is not None:
        dict_of_legend_ids = dict(zip(legend_ids, np.arange(len(legend_ids))))
        for metab in options["linecolor"]:
            axes.get_lines()[dict_of_legend_ids[metab]].set_color(
                options["linecolor"][metab])

def _plot_legend(legend_ids, default_fontsize, **options):
    fontsize = default_fontsize-3
    if options["plot_legend"]:
        plt.legend(legend_ids, loc="center left", bbox_to_anchor=(1.1, 0.5), 
                   prop={"size":fontsize})

def _set_log_scale(ax, default=True, **options):
    log_xscale, log_yscale = _is_log_scale(default, **options)
    ax.set_xscale("log") if log_xscale else ax.set_xscale("linear")
    ax.set_yscale("log") if log_yscale else ax.set_yscale("linear")



# Internal Methods - plot_phase_portrait
def _get_conc_flux_vector(sol_df, mass_obj, start, final, **options):
    dictlist_metabs = DictList(sol_df.columns)

    if isinstance(mass_obj, MassMetabolite):
            pass
    elif isinstance(mass_obj, str):
            mass_obj = dictlist_metabs.get_by_id(mass_obj)
    else:
        raise TypeError("Expected MassMetabolite or string")

    df_xx = sol_df[[mass_obj]]
    df_xx = df_xx.loc[start:final]

    return df_xx, mass_obj

def _annotate_time_range(axes, sol_df, df_x, px, df_y, py, 
                         default_fontsize, **options):
    fontsize = default_fontsize - 3

    if options["truncate"] is False:
        t_value_i = sol_df.index.tolist()[0]
        t_value_f = sol_df.index.tolist()[-1]
    else:
        t_value_i = _truncate(sol_df.index.tolist()[0], 2)
        t_value_f = _truncate(sol_df.index.tolist()[-1], 2)

    annotate_start = "t="+str(t_value_i)
    annotate_final = "t="+str(t_value_f)

    axes.annotate(annotate_start, 
                  xy=(df_x.iloc[0][px], df_y.iloc[0][py]), 
                  xytext=(df_x.iloc[0][px], df_y.iloc[0][py]), 
                  size=fontsize)

    axes.annotate(annotate_final, 
                  xy=(df_x.iloc[-1][px], df_y.iloc[-1][py]), 
                  xytext=(df_x.iloc[-1][px], df_y.iloc[-1][py]), 
                  size=fontsize)

def _annotate_time_scales(axes, list_of_time_scales, np_time_vector, 
                          df_x, px, df_y, py, start, final):
    """Add rectangles for various x,y values corresponding to times where
    pooling occurs"""
    if list_of_time_scales is not None and list_of_time_scales != []:
        i = 0
        for time_scale in list_of_time_scales:
            i+=1
            if time_scale not in np_time_vector:
                x_val, y_val = _interpolate_points(df_x, df_y, 
                                                   time_scale, 
                                                   np_time_vector, 
                                                   start, final)
            else:
                x_val = df_x.loc[time_scale][px]
                y_val = df_y.loc[time_scale][py]

            cx, cy = _make_rectangle(axes, x_val, y_val, i, df_x, px, df_y, py)
            axes.annotate(i, (cx, cy), color="w", weight="bold", 
                          fontsize=6, ha="center", va="center")

def _label_poi(axes, points, ts, np_time_vector, df_x, px, 
               df_y, py, start, final, default_fontsize):
    """Labels points of interest (poi) on the phase portrait"""
    fontsize = default_fontsize-2
    width, height = _get_width_height(df_x, px, df_y, py)

    if points is not None and points != []:
        for point in points:
            label = "t="+str(point)
            if point not in np_time_vector:
                x_val, y_val = _interpolate_points(df_x, df_y, 
                                                   point, 
                                                   np_time_vector, 
                                                   start, final)
            else:
                x_val = df_x.loc[point][px]
                y_val = df_y.loc[point][py]

            if ts is None or ts == []:
                xytext_val = (x_val, y_val)
            else:
                xytext_val = (x_val+width, y_val+height)

            axes.annotate(label, 
                          xy=(x_val, y_val), 
                          xytext=xytext_val, 
                          size=fontsize)

def _validate_datapoints(df_x, px, df_y, py):
    threshold = 0.4e-6

    x_i, y_i = (df_x.iloc[0][px], df_y.iloc[0][py])
    x_f, y_f = (df_x.iloc[-1][px], df_y.iloc[-1][py])

    if abs(x_f-x_i) < threshold or abs(y_f-y_i) < threshold:
        msg =  "datapoints are too close together to make a plot."
        msg += " Try using points that are less than "+str(threshold)
        raise ValueError(msg)



# Internal Methods - plot_tiled_phase_portrait
def _validate_datapoints_tiled(df_x, px, df_y, py):
    try:
        _validate_datapoints(df_x, px, df_y, py)
    except ValueError as v:
        return False
    else:
        return True

def _get_boundaries(df, range_padding):
    boundaries_list = []
    for a in df.columns:
        values = df[a].values
        rmin_, rmax_ = np.min(values), np.max(values)
        rdelta_ext = (rmax_ - rmin_) * range_padding / 2.
        boundaries_list.append((rmin_ - rdelta_ext, rmax_ + rdelta_ext))

    return boundaries_list

def _axes_visibility(fig, ax, **options):
    if options["display_axes"] is False:
        ax.set_xticklabels([])
        ax.set_yticklabels([])
        fig.subplots_adjust(wspace=0.1, hspace=0.1)

def _label_poi_lists(poi, sol_df, xcol, ycol, ax, ts, time, df_x, px, df_y, py, 
                     start, final, default_fontsize):
    if poi is not None:
        if isinstance(poi, tuple):
            ax1, ax2, tiled_poi = _unpack_poi_tiled(sol_df, poi)
            _label_poi_tiled(sol_df, ax1, ax2, xcol, ycol, ax, 
                             tiled_poi, ts, time, df_x, px, df_y, py, 
                             start, final, default_fontsize)
        elif isinstance(poi, list):
            if isinstance(poi[0], tuple):
                for item in poi:
                    ax1, ax2, tiled_poi = _unpack_poi_tiled(sol_df, item)
                    _label_poi_tiled(sol_df, ax1, ax2, xcol, ycol, ax, 
                                     tiled_poi, ts, time, df_x, px, df_y, py, 
                                     start, final, default_fontsize)

def _unpack_poi_tiled(sol_df, poi):
    if poi is None:
        return None, None, None

    dictlist_metabs = DictList(sol_df.columns)
    x, y, list_of_values = poi

    xy = [x, y]
    for i in range(len(xy)):
        if isinstance(xy[i], MassMetabolite):
                pass
        elif isinstance(xy[i], str):
                xy[i] = dictlist_metabs.get_by_id(xy[i])
        else:
            raise TypeError("Expected MassMetabolite or string")

    x = xy[0]
    y = xy[1]

    return x, y, list_of_values

def _label_poi_tiled(sol_df, ax1, ax2, xcol, ycol, ax, tiled_poi, ts, time, 
                     df_x, px, df_y, py, start, final, default_fontsize):
    if (ax1 is not None) and (ax2 is not None) and (tiled_poi is not None):
        result = (ax1.id, ax2.id)
        condition1 = (xcol.id, ycol.id)
        condition2 = (ycol.id, xcol.id)

        if result == condition1 or result == condition2:
            _label_poi(ax, tiled_poi, ts, time, 
                       df_x, px, df_y, py, 
                       start, final, default_fontsize)
    # if ax1 is not None and ax2 is not None:
    #     if ax1.id is xcol.id and ax2.id is ycol.id:
    #         _label_poi(ax, tiled_poi, ts, time, 
    #                    df_x, px, df_y, py, start, final, 
    #                    default_fontsize)
    #     elif ax1.id is ycol.id and ax2.id is xcol.id:
    #         _label_poi(ax, tiled_poi, ts, time,
    #                    df_x, px, df_y, py, start, final,
    #                    default_fontsize)

def _annotate_tiles(ax, df_x, px, df_y, py, default_fontsize):
    fontsize = default_fontsize - 2

    startx, starty = (df_x.iloc[0][px], df_y.iloc[0][py])
    finalx, finaly = (df_x.iloc[-1][px], df_y.iloc[-1][py])

    slope = (finaly - starty) / (finalx - startx)

    str_slope = "slope of "+str(py)+" vs "+str(px)+": "+str(slope)

    ax.annotate(str_slope, 
                  xy=(0.1, 0.5), 
                  xytext=(0.1, 0.5), 
                  size=fontsize)




# Internal Methods - helper functions
def _is_log_scale(default=True, **options):
    plot_function = options["plot_function"].lower()
    log_xscale = default
    log_yscale = default
    if "loglog" in plot_function:
        log_xscale = True
        log_yscale = True
    if "loglinear" in plot_function:
        log_xscale = True
        log_yscale = False
    elif "linearlog" in plot_function:
        log_xscale = False
        log_yscale = True
    elif "linear" in plot_function:
        log_xscale = False
        log_yscale = False

    return log_xscale, log_yscale

def _truncate(f, n):
    """Truncates/pads a float f to n decimal places without rounding"""
    s = '{}'.format(f)
    if 'e' in s or 'E' in s:
        return '{0:.{1}f}'.format(f, n)

    i, p, d = s.partition('.')

    return '.'.join([i, (d+'0'*n)[:n]])

def _interpolate_points(df_x, df_y, time_scale, np_time_vector, start, final):
    df_time = pd.DataFrame(np_time_vector, index=np_time_vector)
    df_time = df_time.loc[start:final]
    interp_fn_x = interp1d(df_time.squeeze(), df_x.squeeze(), kind="cubic")
    interp_fn_y = interp1d(df_time.squeeze(), df_y.squeeze(), kind="cubic")

    x_val = interp_fn_x(time_scale)
    y_val = interp_fn_y(time_scale)

    return x_val, y_val

def _make_rectangle(axes, x, y, num, df_x, px, df_y, py):
    width, height = _get_width_height(df_x, px, df_y, py)

    rect = patches.Rectangle(xy=(x-0.5*width, y-0.5*height), 
                             width=width, height=height, 
                             facecolor="k", label=num)
    axes.add_patch(rect)
    rx, ry = rect.get_xy()
    cx = rx + rect.get_width()/2.0
    cy = ry + rect.get_height()/2.0

    return cx, cy

def _get_width_height(df_x, px, df_y, py):
    width  = 0.01*abs(df_x.iloc[0][px] - df_x.iloc[-1][px])
    height = 0.01*abs(df_y.iloc[0][py] - df_y.iloc[-1][py])

    if width > height:
        width = height
    else:
        height = width

    return width, height



default_options = {
    "plot_function": "LogLogPlot",
    "title": (None, 15),
    "xlabel": (None, 15),
    "ylabel": (None, 15),
    "plot_legend": True,
    "observable": None,
    "set_xlim": None,
    "set_ylim": None,
    "trange": (None, None),
    "truncate": True,
    "linecolor": None,
    "figsize": (None, None),
    "style": None,
    "grid": None,
    "savefig": {
        "fname": None,
        "dpi": None
    }
}

tiled_default_options = {
    "trange": (None, None),
    "truncate": True,
    "linecolor": None,
    "figsize": (None, None),
    "style": None,
    "grid": None,
    "display_axes": True,
    "savefig": {
        "fname": None,
        "dpi": None
    }
}

_base_default_options = {
    "plot_function": "LogLogPlot",
    "title": (None, 15),
    "xlabel": (None, 15),
    "ylabel": (None, 15),
    "plot_legend": True,
    "observable": None,
    "set_xlim": None,
    "set_ylim": None,
    "trange": (None, None),
    "truncate": True,
    "linecolor": None,
    "figsize": (None, None),
    "style": None,
    "grid": None,
    "savefig": {
        "fname": None,
        "dpi": None
    }
}

_base_tiled_default_options = {
    "trange": (None, None),
    "truncate": True,
    "linecolor": None,
    "figsize": (None, None),
    "style": None,
    "grid": None,
    "display_axes": True,
    "savefig": {
        "fname": None,
        "dpi": None
    }
}