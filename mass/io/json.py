# -*- coding: utf-8 -*-

# Compatibility with Python 2.7
from __future__ import absolute_import

# Import Necessary Packages
import io
from math import inf
from collections import OrderedDict
from six import iteritems, string_types
from numpy import bool_, float_, float64
from operator import attrgetter, itemgetter

try:
	import simplejson as json
except ImportError:
	import json

# from cobra
from cobra.core import Gene
# from mass
from mass.core import MassMetabolite, MassReaction, MassModel

# Class begins
## JSON_SPEC
JSON_SPEC = "1"
## Internal Attribute List for MassMetabolites
_REQUIRED_METABOLITE_ATTRIBUTES = ["id", "name", "initial_condition"]
_ORDERED_OPTIONAL_METABOLITE_KEYS = [
	"charge", "compartment", "formula", "notes", "annotation"]
_OPTIONAL_METABOLITE_ATTRIBUTES = {
	"charge": None,
	"compartment": None,
	"formula": None,
	"notes": {},
	"annotation": {}}

## Internal Attribute List for MassReactions
_REQUIRED_REACTION_ATTRIBUTES = [
	"id", "name", "metabolites", "_reversible", "_forward_rate_constant",
	"_equilibrium_constant", "_genes", "_gene_reaction_rule", "ssflux"]
_ORDERED_OPTIONAL_REACTION_KEYS = [
	"subsystem", "_reverse_rate_constant", "notes", "annotation"]
_OPTIONAL_REACTION_ATTRIBUTES = {
	"subsystem": "",
	"_reverse_rate_constant": 0,
	"notes": {},
	"annotation": {}}

## Internal Attribute List for Genes
_REQUIRED_GENE_ATTRIBUTES = ["id", "name"]
_ORDERED_OPTIONAL_GENE_KEYS = ["notes", "annotation"]
_OPTIONAL_GENE_ATTRIBUTES = {
	"notes": {},
	"annotation": {}}

## Internal Attribute List for MassModels
_ORDERED_OPTIONAL_MODEL_KEYS = [
	"name", "_rtype", "_custom_rates", "_custom_parameters",
	"fixed_concentrations", "compartments", "units", "_matrix_type", "_dtype",
	"notes", "annotation"]
_OPTIONAL_MODEL_ATTRIBUTES = {
	"name": None,
	"_rtype": 1,
	"_custom_rates": {},
	"_custom_parameters": {},
	"fixed_concentrations": {},
	"compartments": [],
	"units": {},
	"_matrix_type": "dense",
	"_dtype": float64,
	"notes": {},
	"annotation": {}}

## Public Methods
def model_to_json(model, **kwargs):
	"""Return the model as a JSON document.

	``kwargs`` are passed on to ``json.dumps``.

	Parameters
	----------
	model : mass.MassModel
		The mass model to represent.

	Returns
	-------
	String representation of the mass model as a JSON document.

	See Also
	--------
	write_json_model : Write directly to a file.
	json.dumps : Base function.
	"""
	obj = _model_to_dict(model)
	obj[u"version"] = JSON_SPEC
	return json.dumps(obj, allow_nan=False, **kwargs)

def parse_json_into_model(document):
	"""Load a mass model from a JSON document.

	Parameters
	----------
	document : str
		The JSON document representation of a mass model.

	Returns
	-------
	mass.MassModel
		The mass model as represented in the JSON document.

	See Also
	--------
	read_json_model : Load directly from a file.
	"""
	return _model_from_dict(json.loads(document))

def write_json_model(model, filename, pretty=False, **kwargs):
	"""Write the mass model to a file in JSON format.

	``kwargs`` are passed on to ``json.dump``.

	Parameters
	----------
	model : mass.MassModel
		The mass model to represent.
	filename : str or file-like
		File path or descriptor that the JSON representation should be
		written to.
	pretty : bool, optional
		Whether to format the JSON more compactly (default) or in a more
		verbose but easier to read fashion. Can be partially overwritten by the
		``kwargs``.

	See Also
	--------
	model_to_json : Return a string representation.
	json.dump : Base function.
	"""
	obj = _model_to_dict(model)
	obj[u"version"] = JSON_SPEC

	if pretty:
		dump_opts = {
			"indent": 4, "separators": (",", ": "), "sort_keys": True,
			"allow_nan": False}
	else:
		dump_opts = {
			"indent": 0, "separators": (",", ":"), "sort_keys": False,
			"allow_nan": False}
	dump_opts.update(**kwargs)

	if isinstance(filename, string_types):
		if ".json" not in filename:
			filename += ".json"
		with open(filename, "w") as file_handle:
			json.dump(obj, file_handle, **dump_opts)
	else:
		json.dump(obj, filename, **dump_opts)

def read_json_model(filename):
	"""Load a mass model from a file in JSON format.

	Parameters
	----------
	filename : str or file-like
		File path or descriptor that contains the JSON document describing the
		mass model.

	Returns
	-------
	mass.MassModel
		The mass model as represented in the JSON document.

	See Also
	--------
	parse_json_into_model : Load from a string.
	"""
	if isinstance(filename, string_types):
		with open(filename, "r") as file_handle:
			return _model_from_dict(json.load(file_handle))
	else:
		return _model_from_dict(json.load(filename))

# Internal Methods
def _fix_type(value):
	"""Convert possible types to str, float, and bool"""
	# Because numpy floats can not be pickled to json
	if isinstance(value, string_types):
		return str(value)
	if isinstance(value, float_):
		return float(value)
	if isinstance(value, bool_):
		return bool(value)
	if isinstance(value, set):
		return list(value)
	if value is None:
		return ""
	return value

def _update_optional(mass_object, new_dict, optional_attribute_dict,
					 ordered_keys):
	"""Update new_dict with optional attributes from mass_object"""
	for key in ordered_keys:
		default = optional_attribute_dict[key]
		value = getattr(mass_object, key)
		if value is None or value == default:
			continue
		new_dict[key] = _fix_type(value)

def _inf_to_string(dict_object):
	"""Converts math.inf (float) to "inf" (str) for json conversion"""
	for key in dict_object:
		if dict_object[key] == inf:
			dict_object[key] = "inf"
		elif dict_object[key] == -inf:
			dict_object[key] = "-inf"
	return dict_object

def _metabolite_to_dict(metabolite):
	"""Convert metabolite to a dict.

	Parameters
	----------
	metabolite : mass.MassMetabolite
		The metabolite to reformulate as a dict

	Returns
	-------
	OrderedDict
		A dictionary with elements, 'initial_condition', 'charge',
		'compartment', 'formula', 'notes', and 'annotation'

	See Also
	--------
	mass.io._metabolite_from_dict
	"""
	new_met = OrderedDict()
	for key in _REQUIRED_METABOLITE_ATTRIBUTES:
		new_met[key] = _fix_type(getattr(metabolite, key))
	_update_optional(metabolite, new_met, _OPTIONAL_METABOLITE_ATTRIBUTES,
					 _ORDERED_OPTIONAL_METABOLITE_KEYS)
	return new_met

def _metabolite_from_dict(metabolite_dict):
	"""Build a metabolite from a dict. Models stored in json are first
	formulated as a dict that can be read to mass model using this method
	to convert the metabolites.

	Parameters
	----------
	metabolite_dict : dict
		A dictionary with elements, 'initial_condition', 'charge',
		'compartment', 'formula', 'notes', and 'annotation'

	Returns
	-------
	mass.MassMetabolite
		The generated metabolite

	See Also
	--------
	mass.io._metabolite_to_dict
	"""
	new_metabolite = MassMetabolite()
	for k, v in iteritems(metabolite_dict):
		setattr(new_metabolite, k, v)
	return new_metabolite

def _reaction_to_dict(reaction):
	"""Convert reaction to a dict.

	Parameters
	----------
	reaction : mass.MassReaction
		The reaction to reformulate as a dict

	Returns
	-------
	OrderedDict
		A dictionary with elements, 'metabolites', '_reversible',
		'_forward_rate_constant', '_equilibrium_constant', '_genes',
		'_gene_reaction_rule', and 'ssflux'; where 'metabolites' and 'genes'
		are in turn lists with dictionaries holding all attributes to form the
		corresponding object.

	See Also
	--------
	mass.io._reaction_from_dict
	"""
	new_reaction = OrderedDict()
	for key in _REQUIRED_REACTION_ATTRIBUTES:
		if key != "metabolites":
			new_reaction[key] = _fix_type(getattr(reaction, key))
			continue
		mets = OrderedDict()
		for met in sorted(reaction.metabolites, key=attrgetter("id")):
			mets[str(met)] = reaction.metabolites[met]
		new_reaction["metabolites"] = mets
	_update_optional(reaction, new_reaction, _OPTIONAL_REACTION_ATTRIBUTES,
					 _ORDERED_OPTIONAL_REACTION_KEYS)
	return _inf_to_string(new_reaction)

def _reaction_from_dict(reaction, model):
	"""Build a reaction from a dict. Models stored in json are first formulated
	as a dict that can be read to mass model using this method to convert
	the reactions.

	Parameters
	----------
	reaction : dict
		A dictionary with elements, 'metabolites', '_reversible',
		'_forward_rate_constant', '_equilibrium_constant', '_genes',
		'_gene_reaction_rule', and 'ssflux'; where 'metabolites' and 'genes'
		are in turn lists with dictionaries holding all attributes to form the
		corresponding object.

	model : mass.MassModel
		The model of the reaction to reformulate as a dict

	Returns
	-------
	mass.MassReaction
		The generated reaction

	See Also
	--------
	mass.io._reaction_to_dict
	"""
	new_reaction = MassReaction(id=reaction["id"],
								name=reaction["name"],
								reversible=reaction["_reversible"])
	for k, v in iteritems(reaction):
		if k == 'metabolites':
			new_reaction.add_metabolites(OrderedDict(
				(model.metabolites.get_by_id(str(met)), coeff)
				for met, coeff in iteritems(v)))
		elif v == "inf":
			setattr(new_reaction, k, inf)
		elif v == "-inf":
			setattr(new_reaction, k, -inf)
		else:
			setattr(new_reaction, k, v)
	return new_reaction

def _gene_to_dict(gene):
	"""Converts Gene object to OrderedDict. Identical to cobra method."""
	new_gene = OrderedDict()
	for key in _REQUIRED_GENE_ATTRIBUTES:
		new_gene[key] = _fix_type(getattr(gene, key))
	_update_optional(gene, new_gene, _OPTIONAL_GENE_ATTRIBUTES,
					 _ORDERED_OPTIONAL_GENE_KEYS)
	return new_gene

def _gene_from_dict(gene):
	"""reformulates dict to Gene object. Identical to cobra method"""
	new_gene = Gene(gene["id"])
	for k, v in iteritems(gene):
		setattr(new_gene, k, v)
	return new_gene

def _model_to_dict(model):
	"""Convert model to a dict.

	Parameters
	----------
	model : mass.MassModel
		The model to reformulate as a dict

	Returns
	-------
	OrderedDict
		A dictionary with elements, 'genes', 'compartments', 'id',
		'metabolites', 'notes' and 'reactions'; where 'metabolites', 'genes'
		and 'metabolites' are in turn lists with dictionaries holding all
		attributes to form the corresponding object.

	See Also
	--------
	mass.io._model_from_dict
	"""
	obj = OrderedDict()
	obj["id"] = model.id
	obj["metabolites"] = sorted(
		(_metabolite_to_dict(metabolite) for metabolite in model.metabolites),
		key=itemgetter("id"))
	obj["reactions"] = sorted(
		(_reaction_to_dict(reaction) for reaction in model.reactions),
		key=itemgetter("id"))
	obj["genes"] = sorted(
		(_gene_to_dict(gene) for gene in model.genes), key=itemgetter("id"))
	ics = OrderedDict()
	for met in sorted(model.initial_conditions, key=attrgetter("id")):
		ics[str(met)] = _fix_type(model.initial_conditions[met])
	obj["initial_conditions"] = ics
	# Add custom_rates and custom_parameters here (change to optional)
	_update_optional(model, obj, _OPTIONAL_MODEL_ATTRIBUTES,
					 _ORDERED_OPTIONAL_MODEL_KEYS)
	return obj

def _model_from_dict(obj):
	"""Build a model from a dict. Models stored in json are first formulated as
	a dict that can be read to mass model using this function.

	Parameters
	----------
	obj : dict
		A dictionary with elements, 'genes', 'compartments', 'id',
		'metabolites', 'notes' and 'reactions'; where 'metabolites', 'genes'
		and 'metabolites' are in turn lists with dictionaries holding all
		attributes to form the corresponding object.

	Returns
	-------
	mass.MassModel
		The generated model

	See Also
	--------
	mass.io._model_to_dict
	"""
	if 'reactions' not in obj:
		raise ValueError('Object has no reactions attribute. Cannot load.')
	model = MassModel()
	model.add_metabolites([
		_metabolite_from_dict(metabolite) for metabolite in obj['metabolites']]
	)
	model.genes.extend([_gene_from_dict(gene) for gene in obj['genes']])
	model.add_reactions(
		[_reaction_from_dict(reaction, model) for reaction in obj['reactions']]
	)
	# Add update_initial_conditions
	ics = obj["initial_conditions"]
	ic_dict = {metab: ics[str(metab)] for metab in model.metabolites}
	model.update_initial_conditions(ic_dict)

	for k, v in iteritems(obj):
		if k not in {
		'metabolites', 'reactions', 'genes', 'initial_conditions'}:
			setattr(model, k, v)
	return model
