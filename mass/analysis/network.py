# -*- coding: utf-8 -*-

# Compatibility with Python 2.7
from __future__ import absolute_import

# Import necesary packages
import os
import escher
import ipywidgets as ipyw
from six import iteritems
from IPython.display import display, IFrame, HTML

# from mass
from mass.core import massmodel

# Class begins
## Hex codes for certain colors and default color schemes
l_gray = "#D3D3D3"


def visualize_current_state(model, map_name=None, map_json=None, **kwargs):
	"""Display the current concentrations and steady state fluxes of the model
	onto an Escher map."""
	# Check inputs
	if map_name is None and map_json is None:
		raise ValueError("Either map_name or map_json must be specified")
	if not isinstance(model, massmodel.MassModel):
		raise TypeError("model must be a mass.MassModel")

	concentrations = dict((m.id, ic)
					   for m, ic in iteritems(model.initial_conditions))
	fluxes = dict((r.id, r.ssflux)
					   for r in model.reactions if r.ssflux is not None)
	if concentrations == {}:
		concentrations = None
		ndc = None
	else:
		ndc = l_gray

	if fluxes == {}:
		fluxes = None
		ndf = None
	else:
		ndf = l_gray

	b = escher.Builder(map_name=map_name, map_json=map_json, model=model,
						metabolite_data=concentrations, reaction_data=fluxes,
						metabolite_no_data_color=ndc,
						reaction_no_data_color=ndf, **kwargs)
	return b.display_in_notebook()

def visualize_tsd(model, time_scales,
				map_name=None, map_json=None,
				metabolite_data=None, reaction_data=None,
				highlight_missing=False, **kwargs):
	"""Docstring"""
	# Check inputs
	if map_name is None and map_json is None:
		raise ValueError("Either map_name or map_json must be specified")
	if not isinstance(model, massmodel.MassModel):
		raise TypeError("model must be a mass.MassModel")

	n_ts = len(time_scales)
	if metabolite_data is None:
		metabolite_data = [None] * n_ts
	elif len(metabolite_data) != n_ts:
		raise ValueError("metabolite_data must have the same length as"
							" time_scales.")
	if reaction_data is None:
		reaction_data = [None] * n_ts
	elif len(reaction_data) != n_ts:
		raise ValueError("reaction_data must have the same length as"
							" time_scales.")

	directory = './Escher/%s' % model.id
	_ensure_dir(directory)
	for i in range(n_ts):
		b = escher.Builder(map_name=map_name, map_json=map_json, model=model,
						metabolite_data=metabolite_data[i],
						reaction_data=reaction_data[i],
						metabolite_no_data_color=l_gray,
						reaction_no_data_color=l_gray,
						highlight_missing=highlight_missing,
						**kwargs)

		b.save_html('%s/ts_%s' % (directory, i), overwrite=True, menu='none')

	def display_escher(i):
		display(HTML('<h3> Time Constant %s </h3>' % str(time_scales[i].real)))
		display(IFrame('%s/ts_%s.html'% (directory, i), width=700, height=700))
		return i

	return ipyw.interactive(display_escher,
								i=ipyw.IntSlider(min=0,max=n_ts-1,step=1))

def _ensure_dir(directory):
	if not os.path.exists(directory):
		os.makedirs(directory)
