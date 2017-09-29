# -*- coding: utf-8 -*-

# Compatibility with Python 2.7
from __future__ import absolute_import

# Import necesary packages
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.lines as mlines

from pandas.compat import range, lrange, zip
from pandas.io.formats.printing import pprint_thing
from pandas.plotting._style import _get_standard_colors
from pandas.plotting._tools import _subplots, _set_ticks_props

from math import inf
from six import iterkeys, itervalues
from cycler import cycler

# from cobra
import cobra
from cobra import DictList

# from mass
import mass
from mass import MassMetabolite, MassReaction, MassModel
from mass.core.simulation import *

# Class begins
def plot_simulation(solution_profile, **kwargs):
    """Generates a plot of the data in the solution_profile over time.

    ``kwargs`` are passed on to various matplotlib methods. 
    See options_plot_simulation for a full list.

    Parameters
    ----------
    solution_profile : np.ndarray
        An array containing the time vector and simulated results of solutions.

    Returns
    -------
    plt.gcf()
        A reference to the current figure instance. Shows plot when returned.
        Can be used to modify plot after initial generation.
    """
    sol_profile = dict(solution_profile)
    default_fontsize = 15

    # get options (kwargs)
    options = {}
    for k in options_plot_simulation:
        options[k] = kwargs[k] if k in kwargs else options_plot_simulation[k]

    # build x-axis (get time_vector from sol_profile)
    start = None
    end = None
    if options["tstart"] is not None:
        t = sol_profile["t"]
        start = np.where(t==options["tstart"])[0][0]
    else:
        start = 0
    if options["tfinal"] is not None:
        t = sol_profile["t"]
        end = np.where(t==options["tfinal"])[0][0] + 1
    else:
        end = len(sol_profile["t"])

    np_time_vector = np.array(sol_profile.pop("t"))[start:end]

    # build y-axis (process sol_profile)
    dictlist_metabs = DictList(sol_profile.keys())
    if options["observable"] is None or options["observable"] == []:
        np_sol_profile = np.array(
            [solution for solution in itervalues(sol_profile)]).T
        np_sol_profile = np_sol_profile[start:end]
        if options["plot_legend"]:
            legend_ids = [x.id for x in iterkeys(sol_profile)]
    else:
        observable_profile = {}
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

        for metab in options["observable"]:
            observable_profile[metab] = sol_profile[metab]
        np_sol_profile = np.array(
            [solution for solution in itervalues(observable_profile)]).T
        np_sol_profile = np_sol_profile[start:end]
        if options["plot_legend"]:
            legend_ids = [x.id for x in iterkeys(observable_profile)]

    # options: chart and axis titles
    dict_of_titles = {
        "title": plt.title, 
        "xlabel": plt.xlabel, 
        "ylabel": plt.ylabel
    }
    for name in dict_of_titles:
        if isinstance(options[name], str):
            options[name] = (options[name], default_fontsize)

    # options: log or linear scale
    plot_function = options["plot_function"].lower()
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

    # get colormap
    cm = plt.cm.get_cmap("tab20")
    colors1 = [cm.colors[i] for i in range(len(cm.colors))]
    cm = plt.cm.get_cmap('tab20b')
    colors2 = [cm.colors[i] for i in range(len(cm.colors))]
    cm = plt.cm.get_cmap('tab20c')
    colors3 = [cm.colors[i] for i in range(len(cm.colors))]
    colors = colors1 + colors2 + colors3

    # make plot
    plt.plot(np_time_vector, np_sol_profile)
    axes = plt.gca()
    fig = plt.gcf()

    # plot options: change metabolite color
    if options["linecolor"] is not None:
        dict_of_legend_ids = dict(zip(legend_ids, np.arange(len(legend_ids))))
        for metab in options["linecolor"]:
            axes.get_lines()[dict_of_legend_ids[metab]].set_color(
                options["linecolor"][metab])

    # plot options: fig size (in inches)
    if options["figsize"] != options_plot_simulation["figsize"]:
        fig.set_size_inches(options["figsize"][0], options["figsize"][1])


    # plot options: chart and axis titles
    for name in dict_of_titles:
        if options[name][0] is not None:
            dict_of_titles[name](**{
                "s": options[name][0],
                "fontsize": options[name][1]
            })

    #plot options: plot range
    if options["set_xlim"] is not None:
        axes.set_xlim(options["set_xlim"])
    if options["set_ylim"] is not None:
        axes.set_ylim(options["set_ylim"])

    # plot options: log or linear scale
    plt.xscale("log") if log_xscale else plt.xscale("linear")
    plt.yscale("log") if log_yscale else plt.yscale("linear")

    # plot options: plot legend
    if options["plot_legend"]:
        plt.legend(legend_ids, loc="center left", bbox_to_anchor=(1.1, 0.5), 
                   prop={"size":default_fontsize})

    # plot options: savefig
    if options["savefig"] != options_plot_simulation["savefig"]:
        fig.savefig(**options["savefig"])


    # show plot
    return plt.gcf()



def plot_phase_portrait(solution_profile, x, y, **kwargs):
    """Generates a phase portrait of x,y in the solution_profile over time.

    ``kwargs`` are passed on to various matplotlib methods. 
    See options_plot_simulation for a full list.

    Parameters
    ----------
    solution_profile : np.ndarray
        An array containing the time vector and simulated results of solutions.
    x: mass.MassMetabolite
        The mass metabolite to plot on the x-axis
    y: mass.MassMetabolite
        The mass metabolite to plot on the y-axis

    Returns
    -------
    plt.gcf()
        A reference to the current figure instance. Shows plot when returned.
        Can be used to modify plot after initial generation.
    """

    sol_profile = dict(solution_profile)
    default_fontsize = 15

    # get options (kwargs)
    options = {}
    for k in options_plot_simulation:
        options[k] = kwargs[k] if k in kwargs else options_plot_simulation[k]

    # pop time_vector from sol_profile
    start = None
    end = None
    if options["tstart"] is not None:
        t = sol_profile["t"]
        start = np.where(t==options["tstart"])[0][0]
    else:
        start = 0
    if options["tfinal"] is not None:
        t = sol_profile["t"]
        end = np.where(t==options["tfinal"])[0][0] + 1
    else:
        end = len(sol_profile["t"])

    np_time_vector = np.array(sol_profile.pop("t"))[start:end]

    # build x-axis
    dictlist_metabs = DictList(sol_profile.keys())
    if isinstance(x, MassMetabolite):
            pass
    elif isinstance(x, str):
            x = dictlist_metabs.get_by_id(x)
    else:
        raise TypeError("Expected MassMetabolite or string")

    np_x = np.array(sol_profile[x])[start:end]

    # build y-axis
    dictlist_metabs = DictList(sol_profile.keys())
    if isinstance(y, MassMetabolite):
            pass
    elif isinstance(y, str):
            y = dictlist_metabs.get_by_id(y)
    else:
        raise TypeError("Expected MassMetabolite or string")

    np_y = np.array(sol_profile[y])[start:end]

    # options: chart and axis titles
    dict_of_titles = {
        "title": plt.title, 
        "xlabel": plt.xlabel, 
        "ylabel": plt.ylabel
    }
    for name in dict_of_titles:
        if isinstance(options[name], str):
            options[name] = (options[name], default_fontsize)

    # options: log or linear scale
    plot_function = options["plot_function"].lower()
    log_xscale = False
    log_yscale = False
    if "loglinear" in plot_function:
        log_xscale = True
        log_yscale = False
    elif "linearlog" in plot_function:
        log_xscale = False
        log_yscale = True
    elif "loglog" in plot_function:
        log_xscale = True
        log_yscale = True

    # get colormap
    cm = plt.cm.get_cmap("tab20")
    colors1 = [cm.colors[i] for i in range(len(cm.colors))]
    cm = plt.cm.get_cmap('tab20b')
    colors2 = [cm.colors[i] for i in range(len(cm.colors))]
    cm = plt.cm.get_cmap('tab20c')
    colors3 = [cm.colors[i] for i in range(len(cm.colors))]
    colors = colors1 + colors2 + colors3

    # make plot
    plt.plot(np_x, np_y)
    axes = plt.gca()
    fig = plt.gcf()

    # annotate tstart, tfinal
    annotate_start = "t="+str(np_time_vector[0])
    axes.annotate(annotate_start, xy=(np_x[0], np_y[0]), 
                  xytext=(np_x[0], np_y[0]), size=default_fontsize)
    annotate_final = "t="+str(np_time_vector[-1])
    axes.annotate(annotate_final, xy=(np_x[-1], np_y[-1]), 
                  xytext=(np_x[-1], np_y[-1]), size=default_fontsize)

    # plot options: fig size (in inches)
    if options["figsize"] != options_plot_simulation["figsize"]:
        fig.set_size_inches(options["figsize"][0], options["figsize"][1])


    # plot options: chart and axis titles
    for name in dict_of_titles:
        if options[name][0] is not None:
            dict_of_titles[name](**{
                "s": options[name][0],
                "fontsize": options[name][1]
            })

    # plot options: plot range
    if options["set_xlim"] is not None:
        axes.set_xlim(options["set_xlim"])
    if options["set_ylim"] is not None:
        axes.set_ylim(options["set_ylim"])

    # plot options: log or linear scale
    plt.xscale("log") if log_xscale else plt.xscale("linear")
    plt.yscale("log") if log_yscale else plt.yscale("linear")

    # plot options: plot legend
    if options["plot_legend"]:
        plt.legend(legend_ids, loc="center left", bbox_to_anchor=(1.1, 0.5), 
                   prop={"size":default_fontsize})

    # plot options: savefig
    if options["savefig"] != options_plot_simulation["savefig"]:
        fig.savefig(**options["savefig"])


    # show plot
    return plt.gcf()



def plot_tiled_phase_portrait(frame, alpha=0.5, figsize=None, ax=None, 
                              grid=False, marker='.', range_padding=0.05, 
                              place_tiles="upper", **kwds):
    """Draw a matrix of scatter plots.

    Parameters
    ----------
    frame : DataFrame
    alpha : float, optional
        amount of transparency applied
    figsize : (float,float), optional
        a tuple (width, height) in inches
    ax : Matplotlib axis object, optional
    grid : bool, optional
        setting this to True will show the grid
    diagonal : {'hist', 'kde'}
        pick between 'kde' and 'hist' for
        either Kernel Density Estimation or Histogram
        plot in the diagonal
    marker : str, optional
        Matplotlib marker type, default '.'
    range_padding : float, optional
        relative extension of axis range in x and y
        with respect to (x_max - x_min) or (y_max - y_min),
        default 0.05
    kwds : other plotting keyword arguments
        To be passed to plot function

    Examples
    --------
    >>> df = DataFrame(np.random.randn(1000, 4), columns=['A','B','C','D'])
    >>> scatter_matrix(df, alpha=0.2)
    """
    frame.set_index(keys="t", inplace=True)
    df = frame._get_numeric_data()
    n = df.columns.size
    naxes = n * n
    fig, axes = _subplots(naxes=naxes, figsize=figsize, ax=ax,
                          squeeze=False)

    # no gaps between subplots
    fig.subplots_adjust(wspace=0, hspace=0)

    marker = _get_marker_compat(marker)

    boundaries_list = []
    for a in df.columns:
        values = df[a].values
        rmin_, rmax_ = np.min(values), np.max(values)
        rdelta_ext = (rmax_ - rmin_) * range_padding / 2.
        boundaries_list.append((rmin_ - rdelta_ext, rmax_ + rdelta_ext))

    place_tiles_dict = {"lower": float.__gt__, "upper": float.__lt__}
    for i, a in zip(lrange(n), df.columns):
        for j, b in zip(lrange(n), df.columns):
            ax = axes[i, j]

            if i == j:
                # do nothing along diagonals
                pass

            elif place_tiles_dict[place_tiles](float(i),float(j)):

                ax.plot(df[b], df[a],
                           marker=marker, alpha=alpha, **kwds)

                ax.set_xlim(boundaries_list[j])
                ax.set_ylim(boundaries_list[i])

            ax.set_xlabel(b)
            ax.set_ylabel(a)

            if j != 0:
                ax.yaxis.set_visible(False)
            if i != n - 1:
                ax.xaxis.set_visible(False)

    _set_ticks_props(axes, xlabelsize=8, xrot=90, ylabelsize=8, yrot=0)

    return axes



def _get_marker_compat(marker):
    if mpl.__version__ < '1.1.0' and marker == '.':
        return 'o'
    if marker not in mlines.lineMarkers:
        return 'o'
    return marker



options_plot_simulation = {
    "plot_function": "LogLogPlot",
    "plot_legend": False,
    "title":  (None, 15),
    "xlabel": (None, 15),
    "ylabel": (None, 15),
    "observable": None,
    "linecolor": None,
    "set_xlim": None,
    "set_ylim": None,
    "tstart": None,
    "tfinal": None,
    "figsize": (None, None),
    "savefig": {
        "fname": None,
        "dpi": None
    }
}





