# -*- coding: utf-8 -*-

# Compatibility with Python 2.7
from __future__ import absolute_import

# Import necesary packages
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.lines as mlines

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
    See default_options for a full list.

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
    # Generate seperate mutable copy of solution profile
    sol_profile = dict(solution_profile)
    sol_df = pd.DataFrame(solution_profile)
    default_fontsize = 15

    # Step 0: Get options if provided, else use defaults
    options = _get_options(**kwargs)
    options = _process_title_options(default_fontsize, **options)

    # Step 1: Get time vector for x-axis
    start, final = _get_time_range(t_sol_profile=sol_profile, **options)
    np_time_vector = np.array(sol_profile.pop("t"))[start:final]

    # Step 2: Get conc/flux array for y-axis
    np_sol_profile, legend_ids = _get_conc_flux_array(
        sol_profile, start, final, **options)

    # Step 3: Make plot using options and vectors provided
    plt.plot(np_time_vector, np_sol_profile)
    axes = plt.gca()
    fig = plt.gcf()

    # Step 4: Add remaining plotting options
    plt.rc("axes", prop_cycle=(cycler("color", _get_colormap())))
    _add_custom_linecolors(fig, axes, legend_ids, **options)
    _add_plot_range(axes, **options)
    _plot_title_options(**options)
    _plot_figsize(fig, **options)
    _plot_legend(legend_ids, default_fontsize, **options)
    _option_savefig(**options)

    log_xscale, log_yscale = _is_log_scale(**options)
    plt.xscale("log") if log_xscale else plt.xscale("linear")
    plt.yscale("log") if log_yscale else plt.yscale("linear")

    # Step 5: Return plot/show plot figure
    return plt.gcf()



def plot_phase_portrait(solution_profile, x, y, **kwargs):
    """Generates a phase portrait of x,y in the solution_profile over time.

    ``kwargs`` are passed on to various matplotlib methods. 
    See default_options for a full list.

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

    # Generate seperate mutable copy of solution profile
    sol_profile = dict(solution_profile)
    sol_df = pd.DataFrame(solution_profile)
    default_fontsize = 15

    # Step 0: Get options if provided, else use defaults
    options = _get_options(**kwargs)
    
    dict_of_titles = {
        "title": plt.title, 
        "xlabel": plt.xlabel, 
        "ylabel": plt.ylabel
    }

    # Step 1: Pop time vector from sol_profile
    start, final = _get_time_range(t_sol_profile=sol_profile, **options)
    np_time_vector = np.array(sol_profile.pop("t"))[start:final]

    # Step 2: Get conc/flux vectors for x-axis and y-axis metabolites
    np_x = _get_conc_flux_vector(sol_profile, x, start, final, **options)
    np_y = _get_conc_flux_vector(sol_profile, y, start, final, **options)

    # Step 3: Make plot using options and vectors provided
    plt.plot(np_x, np_y)
    axes = plt.gca()
    fig = plt.gcf()

    #Step 4: Add remaining plotting options

    # Step 5: Return plot/show plot figure
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

# Internal Methods
def _get_options(**kwargs):
    options = {}
    for key in default_options:
        options[key] = kwargs[key] if key in kwargs else default_options[key]
    
    return options

def _get_time_range(t_sol_profile, **options):
    start = 0
    final = len(t_sol_profile["t"])
    if options["tstart"] is not None:
        start = np.where(t_sol_profile["t"]==options["tstart"])[0][0]
    if options["tfinal"] is not None:
        final = np.where(t_sol_profile["t"]==options["tstart"])[0][0] + 1

    return start, final

def _get_conc_flux_array(sol_profile, start, final, **options):
    dictlist_metabs = DictList(sol_profile.keys())
    legend_ids = None

    if options["observable"] is None or options["observable"] == []:
        np_sol_profile = np.array(
            [solution for solution in itervalues(sol_profile)]).T
        np_sol_profile = np_sol_profile[start:final]
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
        np_sol_profile = np_sol_profile[start:final]
        if options["plot_legend"]:
            legend_ids = [x.id for x in iterkeys(observable_profile)]

    return np_sol_profile, legend_ids

def _get_conc_flux_vector(sol_profile, metab, start, final, **options):
    dictlist_metabs = DictList(sol_profile.keys())

    if isinstance(metab, MassMetabolite):
            pass
    elif isinstance(metab, str):
            metab = dictlist_metabs.get_by_id(metab)
    else:
        raise TypeError("Expected MassMetabolite or string")

    np_vector = np.array(sol_profile[metab])[start:final]
    return np_vector

def _is_log_scale(**options):
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

    return log_xscale, log_yscale

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

def _add_plot_range(axes, **options):
    if options["set_xlim"] is not None:
        axes.set_xlim(options["set_xlim"])
        axes.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
    if options["set_ylim"] is not None:
        axes.set_ylim(options["set_ylim"])
        axes.yaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())

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

def _plot_legend(legend_ids, default_fontsize, **options):
    if options["plot_legend"]:
        plt.legend(legend_ids, loc="center left", bbox_to_anchor=(1.1, 0.5), 
                   prop={"size":default_fontsize})

def _option_savefig(**options):
    if options["savefig"] != default_options["savefig"]:
        fig.savefig(**options["savefig"])



default_options = {
    "plot_function": "LogLogPlot",
    "title": (None, 15),
    "xlabel": (None, 15),
    "ylabel": (None, 15),
    "plot_legend": True,
    "observable": None,
    "set_xlim": None,
    "set_ylim": None,
    "tstart": None,
    "tfinal": None,
    "linecolor": None,
    "figsize": (None, None),
    "savefig": {
        "fname": None,
        "dpi": None
    }
}






