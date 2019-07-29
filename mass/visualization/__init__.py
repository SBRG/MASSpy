# -*- coding: utf-8 -*-
r"""Contains function for visualizing simulation results.

This module contains the various functions for visualization of solutions
returned in :class:`.MassSolution`\ s after simulation of models. Note that
to use the :mod:`mass.visualization` functions, additional packages
for visualziation included during the :mod:`mass` installation process as
follows::

    # Installs plotting visualization packages.
    pip install mass["plotting"]
    # Or to install all additional packages.
    pip install mass["all"]

In general, it is recommended to pass an :class:`matplotlib.axes.Axes` instance
into the visualization function in order to guaruntee more control over the
generated plot and ensure that the plot be placed onto a figure as expected.
If an :class:`~matplotlib.axes.Axes` is not passed into the function, then the
currrent axes instance will be used, accessed via::

    matplotlib.pyplot.gca()

All functions will return the axes instance utilzed in generating the plot.

The ``legend`` input format for the plotting function must be one of the
following:

    1. An iterable of legend labels as strings
       (e.g. ``legend=["A", "B", ...]``).

    2. A ``str`` representing the location of the legend, or an
       ``int`` between 0 and 14 (inclusive) corresponding to the
       legend location (e.g. ``legend="best"`` or ``legend=1``).

    3. An iterable of the format ``(labels, loc)`` to set both
       the legend labels and location, where the ``labels`` and ``loc``
       in the iterable follows the formats specified above in both
       **1** and **2**, respectively.
       (e.g. ``legend=(["A", "B", ...], "best")``)

Valid legend location strings and the integer corresponding to the legend
location codes are the following:

    =====================   ===================
    Location String (str)   Location Code (int)
    =====================   ===================
    'best'                  0
    'upper right'           1
    'upper left'            2
    'lower left'            3
    'lower right'           4
    'right'                 5
    'center left'           6
    'center right'          7
    'lower center'          8
    'upper center'          9
    'center'                10
    'right outside'         11
    'left outside'          12
    'upper outside'         13
    'lower outside'         14
    =====================   ===================

Some other important things to note about the legend include the following:

    * Note that the last four "outside" locations are NOT
      :mod:`matplotlib.legend` location values, but rather they utilize the
      legend kwarg ``bbox_to_anchor`` to place the legend outside of the plot.
    * If only legend labels are provided (i.e. format specified above in
      **1**), the default legend location specified in the
      :mod:`matplotlib.rcsetup` will be used.
    * If only the legend location is provided (i.e. format specified above
      in **2**), the corresponding keys in ``mass_solution`` will be used as
      default legend labels.
    * If invalid input is provided (e.g. too many legend labels provided,
      invalid legend location), a warning will be issued and the default
      values will be used.
    * See the :mod:`matplotlib.legend` documentation for additional control
      over the plot legend.

"""
from mass.visualization.time_profiles import plot_time_profile
from mass.visualization.phase_portraits import (
    plot_phase_portrait, plot_tiled_phase_portraits)


__all__ = (
  "plot_time_profile", "plot_phase_portrait", "plot_tiled_phase_portraits")
