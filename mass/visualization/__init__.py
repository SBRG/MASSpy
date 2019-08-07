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

The following are optional ``kwargs`` that can be passed to the functions
of the :mod:`~mass.visualiation` module.

  time_vector :
      ``iterable`` of values to treat as time points for the solutions. If
      provided, the original solutions in the :class:`.MassSolution` input
      will be converted into interpolating functions and the solutions are
      recalculated based on the provided ``time_vector``. If ``None`` then
      the current :attr:`.MassSolution.time` values will be used.

      Default is ``None``.
  plot_function :
      ``str`` representing the plotting function to use. Accepted values are
      the following:

          For functions in :mod:`~mass.visualization.time_profiles` and
          :mod:`~mass.visualization.phase_portraits` modules:

              * ``"plot"`` for a linear x-axis and a linear y-axis
                via :meth:`Axes.plot() <matplotlib.axes.Axes.plot>`
              * ``"semilogx``" for a logarithmic x-axis and a linear y-axis
                via :meth:`Axes.semilogx() <matplotlib.axes.Axes.semilogx>`
              * ``"semilogy"`` for a linear x-axis and a logarithmic y-axis
                via :meth:`Axes.semilogy() <matplotlib.axes.Axes.semilogy>`
              * ``"loglog"`` for a logarithmic x-axis and a logarithmic y-axis
                via :meth:`Axes.loglog() <matplotlib.axes.Axes.loglog>`

           Default is ``"plot"``.
  title :
      Either a ``str`` to set as the title or a ``tuple`` of length 2 where
      the first value is the title string and the second value is a ``dict``
      of font options. Arguments are passed to the corresponding setter method
      :meth:`Axes.set_title() <matplotlib.axes.Axes.set_title>`.

      Default is ``None``.
  xlabel :
      Either a ``str`` to set as the xlabel or a ``tuple`` of length 2 where
      the first value is the xlabel string and the second value is a ``dict``
      of font options. Arguments are passed to the corresponding setter method
      :meth:`Axes.set_xlabel() <matplotlib.axes.Axes.set_xlabel>`.

      Default is ``None``. Not valid for :func:`~.plot_tiled_phase_portraits`
  ylabel :
      Either a ``str`` to set as the ylabel or a ``tuple`` of length 2 where
      the first value is the ylabel string and the second value is a ``dict``
      of font options. Arguments are passed to the corresponding setter method
      :meth:`Axes.set_ylabel() <matplotlib.axes.Axes.set_ylabel>`.

      Default is ``None``. Not valid for :func:`~.plot_tiled_phase_portraits`
  xlim :
      ``tuple`` of form ``(xmin, xmax)`` containing numerical values
      specifying the limits of the x-axis.  Passed to the corresponding setter
      method :meth:`Axes.set_xlim() <matplotlib.axes.Axes.set_xlim>`.
      For :func:`~.plot_tiled_phase_portraits`, the limits will be applied to
      all tiles containing phase portrait plots.

      Default is ``None``.
  ylim :
      ``tuple`` of form ``(ymin, ymax)`` containing numerical values
      specifying the limits of the y-axis.  Passed to the corresponding setter
      method :meth:`Axes.set_ylim() <matplotlib.axes.Axes.set_ylim>`.
      For :func:`~.plot_tiled_phase_portraits`, the limits will be applied to
      all tiles containing phase portrait plots.

      Default is ``None``.
  xmargin :
      ``float`` value greater than -0.5 to set as the padding for the x-axis
      data limits prior to autoscaling. Passed to the corresponding setter
      method :meth:`Axes.set_xmargin() <matplotlib.axes.Axes.set_xmargin>`.
      For :func:`~.plot_tiled_phase_portraits`, the margins will be applied to
      all tiles containing phase portrait plots.

      Default is ``0.15`` for :func:`~.plot_tiled_phase_portraits`,
      otherwise ``None``.
  ymargin :
      ``float`` value greater than -0.5 to set as the padding for the y-axis
      data limits prior to autoscaling. Passed to the corresponding setter
      method :meth:`Axes.set_ymargin() <matplotlib.axes.Axes.set_ymargin>`.
      For :func:`~.plot_tiled_phase_portraits`, the margins will be applied to
      all tiles containing phase portrait plots.

      Default is ``0.15`` for :func:`~.plot_tiled_phase_portraits`,
      otherwise ``None``.
  color :
      Value or ``iterable`` of values representing valid
      :mod:`matplotlib.colors` values to use as line colors. If a single color
      is provided, that color will be applied to all solutions being plotted.
      If an ``iterable`` of color values is provided, it must be equal to the
      number of solutions to be plotted.
      For :func:`~.plot_tiled_phase_portraits`, the colors will be applied to
      all tiles containing phase portrait plots.

      Default is ``None``. Ignored if the kwarg ``prop_cycler`` is also
      provided.
  linestyle :
      Value or ``iterable`` of values representing valid :mod:`matplotlib`
      linestyles_. If a single linestyle is provided, that linestyle will be
      applied to all solutions being plotted. If an ``iterable`` is provided,
      it must be equal to the number of solutions to be plotted.
      For :func:`~.plot_tiled_phase_portraits`, the linestyles will be applied
      to all tiles containing phase portrait plots.

      Default is ``None`` to use default value in :mod:`matplotlib.rcsetup`.
      Ignored if the kwarg ``prop_cycler`` is also provided.
  linewidth :
      ``float`` value representing the linewidth (in points) to set.

      Default is ``None`` to use default value in :mod:`matplotlib.rcsetup`.
      Ignored if the kwarg ``prop_cycler`` is also provided.
  marker :
      Value or ``iterable`` of values representing valud :mod:`matplotlib`
      marker_ values to use as line markers. If a single marker is provided,
      that marker will be applied to all solutions being plotted.
      If an ``iterable`` is provided, it must be equal to the number of
      solutions to be plotted.
      For :func:`~.plot_tiled_phase_portraits`, the markers will be applied
      to all tiles containing phase portrait plots.

      Default is ``None`` to use default value in :mod:`matplotlib.rcsetup`.
      Ignored if the kwarg ``prop_cycler`` is also provided.
  markersize :
      ``float`` value representing the size of the marker (in points) to set.
      For :func:`~.plot_tiled_phase_portraits`, the markersizes will be applied
      to all tiles containing phase portrait plots. Ignored if the kwarg
      ``prop_cycler`` is also provided.

      Default is ``None`` to use default value in :mod:`matplotlib.rcsetup`.
      Ignored if the kwarg ``prop_cycler`` is also provided.
  grid :
      Either a ``bool`` or a ``tuple`` of form ``(which, axis)`` where the
      values for ``which`` and ``axis`` are one of the following:

          * ``which`` arguement must be ``"major"``, ``"minor"``, or ``"both"``
          * ``axis`` argument must be ``"x"``, ``"y"``, or ``"both"``

      If ``grid=False`` then grid lines will be removed.
      If ``grid=True``, then grid lines will be created based on the
      default values in :mod:`matplotlib.rcsetup` unless overridden by
      ``grid_color``, ``grid_linestyle``, or ``grid_linewidth`` kwargs. Passed
      to the setter method :meth:`Axes.grid() <matplotlib.axes.Axes.grid>`.
      For :func:`~.plot_tiled_phase_portraits`, the grid lines will be
      applied to all tiles containing phase portrait plots.

      Default is ``None``.
  grid_color :
      Value representing a valid :mod:`matplotlib.colors` value to use as
      the color of the grid lines. For :func:`~.plot_tiled_phase_portraits`,
      the color of the grid lines will be applied to all tiles containing
      phase portrait plots.

      Default is ``None`` to use default value in :mod:`matplotlib.rcsetup`.
      Ignored if the kwarg ``grid`` is ``None``.
  grid_linestyle :
      Value representing a valid :mod:`matplotlib` value to use as the
      linestyles_ of the grid lines. For :func:`~.plot_tiled_phase_portraits`,
      the linestyle of the grid lines will be applied to all tiles containing
      phase portrait plots.

      Default is ``None`` to use default value in :mod:`matplotlib.rcsetup`.
      Ignored if the kwarg ``grid`` is ``None``.
  grid_linewidth :
       ``float`` value representing the grid linewidth (in points) to set.

      Default is ``None`` to use default value in :mod:`matplotlib.rcsetup`.
      Ignored if the kwarg ``grid`` is ``None``.
  legend_ncol :
      ``int`` indicating the number of columns to use in the legend.
      If ``None`` the the following formula is applied::

          ncols = int(ceil(sqrt(N_total) / 3))

      where ``N_total`` is equal to the total number of solution lines.

      Default is ``None``. Not valid for :func:`~.plot_tiled_phase_portraits`
  annotate_time_points :
    Either the string ``"endpoints"`` or an ``iterable`` containing the
    numerical values for the time points of interest to be annotated by
    plotting points for the solutions at the given time points.
    If ``annotate_time_points="endpoints"`` then only the initial and final
    time points will be utilized. If ``None`` then no time points will be
    annotated. For :func:`~.plot_tiled_phase_portraits`, the time points
    will be applied to all tiles containing phase portrait plots.

    Default is ``endpoints`` for :func:`~.plot_tiled_phase_portraits`,
    otherwise ``None``.
  annotate_time_points_color :
      Value or ``iterable`` of values representing valid
      :mod:`matplotlib.colors` values to use as time point colors. If a single
      color is provided, that color will be applied to all time points being
      plotted. If an ``iterable`` of color values is provided, it must be
      equal to the number of time points to be plotted.
      For :func:`~.plot_tiled_phase_portraits`, the colors will be applied to
      all tiles containing phase portrait plots.

      Default is ``["red", "blue"]`` for :func:`~.plot_tiled_phase_portraits`,
      otherwise ``None``.
  annotate_time_points_marker :
      Value or ``iterable`` of values representing valud :mod:`matplotlib`
      marker_ values to use as time point markers. If a single marker is
      provided, that marker will be applied to all solutions being plotted.
      If an ``iterable`` is provided, it must be equal to the number of time
      points to be plotted.
      For :func:`~.plot_tiled_phase_portraits`, the markers will be applied
      to all tiles containing phase portrait plots.

      Default is ``["o", "D"]`` for :func:`~.plot_tiled_phase_portraits`,
      otherwise ``None``.
  annotate_time_points_markersize :
      ``float`` value representing the size of the marker (in points) to set.
      For :func:`~.plot_tiled_phase_portraits`, the markersizes will be applied
      to all tiles containing phase portrait plots.

      Default is ``None`` to use default value in :mod:`matplotlib.rcsetup`.
  annotate_time_points_legend_loc :
      A ``str`` representing the location of the legend, or an ``int``
      between 0 and 14 (inclusive) corresponding to the location to use for
      the legend of annotated time points. Cannot be the same location value
      as the plot's ``legend`` location.

      Default is ``None`` to use either the ``"right outside"`` location,
      or ``"left outside"`` location if ``"right outside"`` is occupied.
  prop_cycle :
      A valid :func:`matplotlib.rcsetup.cycler` instance to use in the plot.
      If provided, then the ``color``, ``linestyle``, ``linewidth``,
      ``marker``, and ``markersize`` kwargs are ignored.

      Default is ``None``.
  tile_ticks_on :
      ``bool`` indicating whether to leave tick marks on tiles containing
      phase portraits.
  
      Default is ``False``.
      Only valid for :func:`~.plot_tiled_phase_portraits`.
  tile_xlabel_fontdict :
      :class:`Font properties <matplotlib.font_manager.FontProperties>`
      to set using a ``dict``. Applied to all tile ylabels.

      Default is ``None``.
      Only valid for :func:`~.plot_tiled_phase_portraits`.
  tile_ylabel_fontdict :
      :class:`Font properties <matplotlib.font_manager.FontProperties>`
      to set using a ``dict``. Applied to all tile ylabels.

      Default is ``None``.
      Only valid for :func:`~.plot_tiled_phase_portraits`.
  data_tile_fontsize :
      Valid :class:`fontsize <matplotlib.font_manager.FontProperties>` value
      representing the fontsize to be used for the data values on tiles
      displaying additional data.

      Default is ``large``.
      Only valid for :func:`~.plot_tiled_phase_portraits`. Ignored if no
      additional data is provided (i.e. ``additional_data=None``).
  data_tile_color :
      Value representing a valid :mod:`matplotlib.colors` value to use as
      the facecolor of the tiles displaying additional data.

      Default is ``None`` to utilize the color ``"lightgray"``.
      Only valid for :func:`~.plot_tiled_phase_portraits`. Ignored if no
      additional data is provided (i.e. ``additional_data=None``).
  diag_tile_color :
      Value representing a valid :mod:`matplotlib.colors` value to use as
      the facecolor of the tiles on the diagonal.

      Default is ``None`` to utilize the color ``"black"``.
      Only valid for :func:`~.plot_tiled_phase_portraits`.
  empty_tile_color :
      Value representing a valid :mod:`matplotlib.colors` value to use as
      the facecolor of the empty tiles.

      Default is ``None`` to utilize the color ``"white"``.
      Only valid for :func:`~.plot_tiled_phase_portraits`.

.. _marker:  https://matplotlib.org/api/markers_api.html
.. _linestyles: https://matplotlib.org/gallery/lines_bars_and_markers/linestyles

"""  # noqa
from mass.visualization.time_profiles import plot_time_profile
from mass.visualization.phase_portraits import (
    plot_phase_portrait, plot_tiled_phase_portraits)


__all__ = (
    "plot_time_profile", "plot_phase_portrait", "plot_tiled_phase_portraits")