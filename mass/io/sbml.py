# -*- coding: utf-8 -*-

from __future__ import absolute_import

import re
from ast import And, BoolOp, Name, Or
from bz2 import BZ2File
from collections import defaultdict
from decimal import Decimal
from gzip import GzipFile
from tempfile import NamedTemporaryFile
from warnings import catch_warnings, simplefilter, warn

from six import iterkeys, iteritems, string_types

import sympy
from sympy import Basic, Symbol
from sympy.printing.mathml import mathml

from cobra.core import Gene
from cobra.core.gene import parse_gpr
from cobra.manipulation.modify import _renames
from cobra.manipulation.validate import check_metabolite_compartment_formula
# For Gene and GPR handling
from cobra.io.sbml3 import annotate_cobra_from_sbml, annotate_sbml_from_cobra

from mass.core import MassMetabolite, MassReaction, MassModel
from mass.core.expressions import strip_time
from mass.exceptions import MassSBMLError

try:
	from lxml.etree import (
		parse, Element, SubElement, ElementTree, register_namespace,
		ParseError, XPath, fromstring, tostring)

	_with_lxml = True
except ImportError:
	_with_lxml = False
	try:
		from xml.etree.cElementTree import (
			parse, Element, SubElement, ElementTree, register_namespace,
			ParseError, fromstring, tostring)
	except ImportError:
		XPath = None
		from xml.etree.ElementTree import (
			parse, Element, SubElement, ElementTree, register_namespace,
			ParseError, fromstring, tostring)

try:
	import libsbml
except ImportError:
	raise MassSBMLError("Need to install libsbml to proceed")

## Set a float infinity (Compatibility with Python 2.7)
inf = float('inf')
# deal with namespaces
namespaces = {"fbc": "http://www.sbml.org/sbml/level3/version1/fbc/version2",
			  "sbml": "http://www.sbml.org/sbml/level3/version1/core",
			  "rdf": "http://www.w3.org/1999/02/22-rdf-syntax-ns#",
			  "bqbiol": "http://biomodels.net/biology-qualifiers/",
			  "mml": "http://www.w3.org/1998/Math/MathML"}

for key in namespaces:
	register_namespace(key, namespaces[key])

def ns(query):
	"""replace prefixes with namespace"""
	for prefix, uri in iteritems(namespaces):
		query = query.replace(prefix + ":", "{" + uri + "}")
	return query

# XPATH query wrappers
fbc_prefix = "{" + namespaces["fbc"] + "}"
sbml_prefix = "{" + namespaces["sbml"] + "}"

SBML_DOT = "__SBML_DOT__"
# FBC TAGS
OR_TAG = ns("fbc:or")
AND_TAG = ns("fbc:and")
GENEREF_TAG = ns("fbc:geneProductRef")
GPR_TAG = ns("fbc:geneProductAssociation")
GENELIST_TAG = ns("fbc:listOfGeneProducts")
GENE_TAG = ns("fbc:geneProduct")
# XPATHS
COMPARTMENT_XPATH = ns("sbml:listOfCompartments/sbml:compartment")
GENES_XPATH = GENELIST_TAG + "/" + GENE_TAG
SPECIES_XPATH = ns("sbml:listOfSpecies/sbml:species[@boundaryCondition='%s']")
OBJECTIVES_XPATH = ns("fbc:objective[@fbc:id='%s']/"
					  "fbc:listOfFluxObjectives/"
					  "fbc:fluxObjective")

if _with_lxml:
	RDF_ANNOTATION_XPATH = ("sbml:annotation/rdf:RDF/"
							"rdf:Description[@rdf:about=$metaid]/"
							"*[self::bqbiol:isEncodedBy or self::bqbiol:is]/"
							"rdf:Bag/rdf:li/@rdf:resource")
	extract_rdf_annotation = XPath(RDF_ANNOTATION_XPATH, namespaces=namespaces,
								   smart_strings=False)
else:
	RDF_ANNOTATION_XPATH = ns("sbml:annotation/rdf:RDF/"
							  "rdf:Description[@rdf:about='%s']/"
							  "bqbiol:isEncodedBy/"
							  "rdf:Bag/rdf:li[@rdf:resource]")

	def extract_rdf_annotation(sbml_element, metaid):
		search_xpath = RDF_ANNOTATION_XPATH % metaid
		for i in sbml_element.iterfind(search_xpath):
			yield _get_attr(i, "rdf:resource")
		for i in sbml_element.iterfind(search_xpath.replace(
				"isEncodedBy", "is")):
			yield _get_attr(i, "rdf:resource")


# Public Methods
def parse_xml_into_model(xml, number=float):
	"""Load a mass model from an xml object.

	Parameters
	----------
	xml : xml object
		The string-like representation of a mass model in xml format.
	number : float or int
		The typecasting to use for certain number parsing instances

	Returns
	-------
	mass.MassModel
		The mass model as represented in the SBML/XML document.

	See Also
	--------
	read_sbml_model : Load directly from a file.
	"""
	# add model
	xml_model = xml.find(ns("sbml:model"))
	if _get_attr(xml_model, "fbc:required") == "false":
		warn('loading SBML model with fbc:required="false"')

	model_id = _get_attr(xml_model, "id")
	model = MassModel(model_id)
	model.name = xml_model.get("name")

	# add in units
	unit_dict = {}
	for unit_def in xml_model.findall(ns(
		"sbml:listOfUnitDefinitions/sbml:unitDefinition")):
		unit = _get_attr(unit_def, "id", require=True)
		if unit is None:
			continue
		if "mole" in unit.lower():
			unit_dict["N"] = unit
		if ("liter" in unit.lower()) or ("litre" in unit.lower()):
			unit_dict["Vol"] = unit
		if ("hour" in unit.lower()) or ("sec" in unit.lower()):
			unit_dict["Time"] = unit

	model.units = unit_dict

	# add compartments
	model.compartments = {c.get("id"): c.get("name") for c in
						  xml_model.findall(COMPARTMENT_XPATH)}
	# Detect fixed concentration metabolites (boundary metabolites)
	boundary_metabolites = {_clip(i.get("id"), "M_")
							for i in xml_model.findall(SPECIES_XPATH % 'true')}
	bound_dict = {}
	# add metabolites (species)
	for species in xml_model.findall(ns("sbml:listOfSpecies/sbml:species")):
		met = _get_attr(species, "id", require=True)
		met = MassMetabolite(_clip(met, "M_"))
		met.name = species.get("name")
		_annotate_mass_from_sbml(met, species)
		met.compartment = species.get("compartment")
		met.charge = _get_attr(species, "fbc:charge", float)
		met.formula = _get_attr(species, "fbc:chemicalFormula")
		model.add_metabolites([met])
		model.update_initial_conditions(
			{met: _get_attr(species, "initialConcentration", float)})
		if met.id in boundary_metabolites:
			bound_dict[met] = _get_attr(species, "initialConcentration",
										 float)
	# add fixed concentration metabolites (boundary metabolites)
	model.add_fixed_concentrations(bound_dict)

	# add genes
	for sbml_gene in xml_model.iterfind(GENES_XPATH):
		gene_id = _get_attr(sbml_gene, "fbc:id").replace(SBML_DOT, ".")
		gene = Gene(_clip(gene_id, "G_"))
		gene.name = _get_attr(sbml_gene, "fbc:name")
		if gene.name is None:
			gene.name = _get_attr(sbml_gene, "fbc:label")
		annotate_cobra_from_sbml(gene, sbml_gene)
		model.genes.append(gene)

	def process_gpr(sub_xml):
		"""recursively convert gpr xml to a gpr string"""
		if sub_xml.tag == OR_TAG:
			return "( " + ' or '.join(process_gpr(i) for i in sub_xml) + " )"
		elif sub_xml.tag == AND_TAG:
			return "( " + ' and '.join(process_gpr(i) for i in sub_xml) + " )"
		elif sub_xml.tag == GENEREF_TAG:
			gene_id = _get_attr(sub_xml, "fbc:geneProduct", require=True)
			return _clip(gene_id, "G_")
		else:
			raise Exception("unsupported tag " + sub_xml.tag)

	# add reactions
	custom_rate_dict = {}
	custom_param_dict = {}
	reactions = []
	for sbml_reaction in xml_model.iterfind(
			ns("sbml:listOfReactions/sbml:reaction")):
		reaction = _get_attr(sbml_reaction, "id", require=True)
		isReversible = _get_attr(sbml_reaction, "reversible")
		if isReversible == "false":
			isReversible = False
		else:
			isReversible = True
		reaction = MassReaction(_clip(reaction, "R_"), reversible=isReversible)
		reaction.name = sbml_reaction.get("name")
		_annotate_mass_from_sbml(reaction, sbml_reaction)
		reaction.subsystem = _get_attr(sbml_reaction, "compartment")

		# add mathml rate law extraction here
		result_from_mathml = ""
		for local_rate in sbml_reaction.findall(
			ns("sbml:kineticLaw/mml:math")):
			ast = libsbml.readMathMLFromString(tostring(local_rate,
										encoding="utf-8").decode("utf-8"))
			result_from_mathml = libsbml.formulaToL3String(ast)
			if result_from_mathml is None:
				warn("%s has imparsible MathML" % reaction.id)
				continue
			else:
				result_from_mathml = result_from_mathml.replace("^", "**")

		# add rate constants, ssflux, custom paramaters (local parameters)
		custom_param_ids = []
		custom_param_vals = []
		for local_param in sbml_reaction.findall(ns(
			"sbml:kineticLaw/sbml:listOfLocalParameters/sbml:localParameter")):
			pid = local_param.get("id")
			lpval = local_param.get("value")
			if pid == "ssflux_"+reaction.id and lpval is not None:
				reaction.ssflux = number(local_param.get("value"))
			elif pid == "kf_"+reaction.id and lpval is not None:
				reaction.kf = number(local_param.get("value"))
			elif pid == "Keq_"+reaction.id and reaction.reversible is True:
				if local_param.get("value") == "inf":
					reaction.Keq = inf
				elif local_param.get("value") == "-inf":
					reaction.Keq = -inf
				elif local_param.get("value") == None:
					pass
				else:
					reaction.Keq = number(local_param.get("value"))
			elif pid == "kr_"+reaction.id and reaction.reversible is True:
				reaction.kr = local_param.get("value")
			else:
				if lpval is None:
					pass
				else:
					custom_param_ids.append(pid)
					custom_param_vals.append(number(local_param.get("value")))
		custom_param_dict[reaction] = dict(zip(custom_param_ids,
												custom_param_vals))

		reactions.append(reaction)
		# add stoichiometry (reactants, products, modifiers)
		stoichiometry = defaultdict(lambda: 0)
		for species_reference in sbml_reaction.findall(
				ns("sbml:listOfReactants/sbml:speciesReference")):
			met_name = _clip(species_reference.get("species"), "M_")
			stoichiometry[met_name] -= \
				float(species_reference.get("stoichiometry"))
		for species_reference in sbml_reaction.findall(
				ns("sbml:listOfProducts/sbml:speciesReference")):
			met_name = _clip(species_reference.get("species"), "M_")
			stoichiometry[met_name] += \
				_get_attr(species_reference, "stoichiometry",
						   type=float, require=True)
		# needs to have keys of metabolite objects, not ids
		object_stoichiometry = {}
		for met_id in stoichiometry:
			try:
				metabolite = model.metabolites.get_by_id(met_id)
			except KeyError:
				warn("ignoring unknown metabolite '%s' in reaction %s" %
					 (met_id, reaction.id))
				continue
			object_stoichiometry[metabolite] = stoichiometry[met_id]
		reaction.add_metabolites(object_stoichiometry)
		# set gene reaction rule (gene product association)
		gpr_xml = sbml_reaction.find(GPR_TAG)
		if gpr_xml is not None and len(gpr_xml) != 1:
			warn("ignoring invalid geneAssociation for " + repr(reaction))
			gpr_xml = None
		gpr = process_gpr(gpr_xml[0]) if gpr_xml is not None else ''
		# remove outside parenthesis, if any
		if gpr.startswith("(") and gpr.endswith(")"):
			gpr = gpr[1:-1].strip()
		gpr = gpr.replace(SBML_DOT, ".")
		reaction.gene_reaction_rule = gpr

		# Test to see mathml rate matches standard massreaction rate
		# if no match, add mathml rate as custom rate law for reaction
		for rxn_rate_type in [1, 2, 3]:
			reaction._rtype = rxn_rate_type
			curr_string = strip_time(reaction.rate)[0]
			curr_string = _remove_trailing_zeroes(str(curr_string))
			if result_from_mathml is None:
				pass
			elif curr_string == _remove_trailing_zeroes(result_from_mathml):
				break

		curr_string = strip_time(reaction.rate)[0]
		curr_string = _remove_trailing_zeroes(str(curr_string))
		if result_from_mathml is None:
			pass
		elif curr_string != _remove_trailing_zeroes(result_from_mathml):
			custom_rate_dict[reaction] = result_from_mathml
	try:
		model.add_reactions(reactions)
	except ValueError as e:
		warn(str(e))

	for custom_rxn, custom_rate in iteritems(custom_rate_dict):
		model.add_custom_rate(
			custom_rxn, custom_rate, custom_param_dict[custom_rxn])

	# add external metabolites (parameters)
	ext_dict = {}
	for param in xml_model.findall(ns("sbml:listOfParameters/sbml:parameter")):
		ext_dict[_get_attr(param, "id", require=True)] = _get_attr(
														param, "value", float)

	# add fixed concentration metabolites (external metabolites)
	model.add_fixed_concentrations(ext_dict)

	return model

def model_to_xml(model):

	"""Return the model as an xml object.

	Parameters
	----------
	model : mass.MassModel
		The mass model to represent.

	Returns
	-------
	ElementTree.Element
		string-like representation of the mass model as an xml object.

	See Also
	--------
	write_sbml_model : Write directly to a file.
	"""

	# add in model
	xml = Element("sbml", xmlns=namespaces["sbml"], level="3", version="1",
				  sboTerm="SBO:0000624")
	_set_attrib(xml, "fbc:required", "true")
	xml_model = SubElement(xml, "model")
	_set_attrib(xml_model, "fbc:strict", "true")
	if model.id is not None:
		xml_model.set("id", model.id)
	if model.name is not None:
		xml_model.set("name", model.name)

	 # add in units
	if len(model.units) != 0:
		units_list = SubElement(xml_model, "listOfUnitDefinitions")
		for key, unit in iteritems(model.units):
			unit_def = SubElement(units_list, "unitDefinition", id=unit)
			scale = _get_scale(unit.lower())

			unit_str = ""
			special = ""
			for i in range(len(_SBML_base_units)):
				if _SBML_base_units[i] in unit.lower():
					unit_str = _SBML_base_units[i]
				elif "liter" in unit.lower():
					unit_str = "litre"
				elif "meter" in unit.lower():
					unit_str = "metre"
				elif "hour" in unit.lower():
					unit_str = "second"
					special = "hour"
				elif "minute" in unit.lower():
					unit_str = "minute"
					special = "minute"

			if (unit_str == "") or (unit_str is None):
				warn("Units are not SBML-compliant")
				unit_str = unit

			list_units = SubElement(unit_def, "listOfUnits")
			unit_elem = SubElement(list_units, "unit", kind=unit_str)
			if "minute" in special:
				_set_attrib(unit_elem, "multiplier", 60)
			elif "hour" in special:
				_set_attrib(unit_elem, "multiplier", 3600)
			else:
				_set_attrib(unit_elem, "multiplier", 1)

			_set_attrib(unit_elem, "exponent", 1)
			_set_attrib(unit_elem, "scale", scale)


	# add in compartments
	if len(model.compartments) != 0:
		compartments_list = SubElement(xml_model, "listOfCompartments")
		for compartment, name in iteritems(model.compartments):
			SubElement(compartments_list, "compartment", id=compartment,
						name=name, constant="true")

	# add in species (metabolites)
	for met in model.metabolites:
		species = SubElement(SubElement(xml_model, "listOfSpecies"), "species",
							 id="M_" + met.id,
							 # Required SBML parameter
							 hasOnlySubstanceUnits="false")
		_set_attrib(species, "name", met.name)
		_annotate_sbml_from_mass(species, met)
		_set_attrib(species, "compartment", met.compartment)
		_set_attrib(species, "fbc:charge", met.charge)
		_set_attrib(species, "fbc:chemicalFormula", met.formula)
		_set_attrib(species, "initialConcentration",
				   model.initial_conditions[met])
		# differentiate fixed concentration metabolites
		if met in model.fixed_concentrations:
			_set_attrib(species, "constant", "true")
			_set_attrib(species, "boundaryCondition", "true")
		else:
			_set_attrib(species, "constant", "false")
			_set_attrib(species, "boundaryCondition", "false")

	# add in parameters (external metabolites)

	for ext in model.get_external_metabolites:
		if ext in iterkeys(model.fixed_concentrations):
			param = SubElement(SubElement(xml_model, "listOfParameters"),
											"parameter", id=ext)
			_set_attrib(param, "value", model.fixed_concentrations[ext])
			_set_attrib(param, "constant", "true")

	# add in genes
	if len(model.genes) != 0:
		for gene in model.genes:
			gene_id = gene.id.replace(".", SBML_DOT)
			sbml_gene = SubElement(SubElement(xml_model, GENELIST_TAG),
						GENE_TAG)
			_set_attrib(sbml_gene, "fbc:id", "G_" + gene_id)
			if gene.name is None or len(gene.name) == 0:
				gene.name = gene.id
			_set_attrib(sbml_gene, "fbc:label", gene_id)
			_set_attrib(sbml_gene, "fbc:name", gene.name)
			annotate_sbml_from_cobra(sbml_gene, gene)

	# add in reactions
	rxns_list = SubElement(xml_model, "listOfReactions")
	for rxn, rate in iteritems(strip_time(model.rates)):
		rxn_id = "R_%s" % rxn.id
		sbml_rxn = SubElement(rxns_list, "reaction", id=rxn_id,
								   reversible=str(rxn.reversible).lower(),
								   # Required SBML parameter
								   fast="false")
		_set_attrib(sbml_rxn, "name", rxn.name)
		_set_attrib(sbml_rxn, "compartment", rxn.subsystem)
		_annotate_sbml_from_mass(sbml_rxn, rxn)
		# add in kinetic law (rate law + associated constants)
		sbml_kinetic_law = SubElement(sbml_rxn, "kineticLaw")
		# add in math (rate law in mathml)
		rate_str = str(rate)
		if "**" in rate_str:
			rate_str = rate_str.replace("**", "^")
		new_str = libsbml.writeMathMLToString(libsbml.parseL3Formula(rate_str))
		if new_str is None:
			raise MassSBMLError("Error: Can't parse rate law expression "
							"for %s with rate law:" % (rxn.name, str(rate)))
		else:
			new_str = new_str.replace(
				'<?xml version="1.0" encoding="UTF-8"?>\n', "")
			sbml_kinetic_law.append(fromstring(new_str))

		# add in local parameters (kf, Keq, kr, ssflux)
		sbml_param_list = SubElement(sbml_kinetic_law, "listOfLocalParameters")
		sbml_kf = SubElement(sbml_param_list, "localParameter",
						id="kf_%s" % rxn.id, name="kf_%s" % rxn.id)
		sbml_Keq = SubElement(sbml_param_list, "localParameter",
						id="Keq_%s" % rxn.id, name="Keq_%s" % rxn.id)
		sbml_kr = SubElement(sbml_param_list, "localParameter",
						id="kr_%s" % rxn.id, name="kr_%s" % rxn.id)
		sbml_ssflux = SubElement(sbml_param_list, "localParameter",
						id="ssflux_%s" % rxn.id, name="ssflux_%s" % rxn.id)
		_set_attrib(sbml_kf, "value", rxn.kf)
		_set_attrib(sbml_Keq, "value", rxn.Keq)
		_set_attrib(sbml_kr, "value", rxn.kr)
		_set_attrib(sbml_ssflux, "value", rxn.ssflux)

		# add in local parameters (custom parameters)
		c_param_list = list(iterkeys(model.custom_parameters))
		try:
			c_rate = model.custom_rates[rxn]
		except KeyError:
			pass
		else:
			for sym in c_rate.atoms(Symbol):
				if sym is Symbol("t"):
					continue
				elif str(sym) in c_param_list:
					v = model.custom_parameters[str(sym)]
					sbml_cparam = SubElement(sbml_param_list, "localParameter",
											 id=str(sym), name=str(sym))
					_set_attrib(sbml_cparam, "value", v)

		# add in reactants, products, modifiers (stoichiometry)
		reactants = {}
		products = {}
		for metabolite, stoichiomety in iteritems(rxn._metabolites):
			met_id = "M_" + metabolite.id
			if stoichiomety > 0:
				products[met_id] = _strnum(stoichiomety)
			else:
				reactants[met_id] = _strnum(-stoichiomety)
		if len(reactants) > 0:
			reactant_list = SubElement(sbml_rxn, "listOfReactants")
			for met_id, stoichiomety in sorted(iteritems(reactants)):
				SubElement(reactant_list, "speciesReference", species=met_id,
						   stoichiometry=stoichiomety, constant="true")
		if len(products) > 0:
			product_list = SubElement(sbml_rxn, "listOfProducts")
			for met_id, stoichiomety in sorted(iteritems(products)):
				SubElement(product_list, "speciesReference", species=met_id,
						   stoichiometry=stoichiomety, constant="true")

		# add in gene product association (gene reaction rule)
		gpr = rxn.gene_reaction_rule
		if gpr is not None and len(gpr) > 0:
			gpr = gpr.replace(".", SBML_DOT)
			gpr_xml = SubElement(sbml_rxn, GPR_TAG)
			try:
				parsed, _ = parse_gpr(gpr)
				_construct_gpr_xml(gpr_xml, parsed.body)
			except Exception as e:
				print("failed on '%s' in %s" %
					  (rxn.gene_reaction_rule, repr(rxn)))
				raise e

	return xml

def write_sbml_model(model, filename, use_fbc_package=True, **kwargs):
	"""Write the mass model to a file in SBML/XML format.

	``kwargs`` are passed on to ``ElementTree.write``.

	Parameters
	----------
	model : mass.MassModel
		The mass model to represent.
	filename : str or file-like
		File path or descriptor that the SBML/XML representation should be
		written to.
	use_fbc_package: bool
		Option to use fbc package for generating SBML/XML document

	See Also
	--------
	model_to_xml : Return an xml object.
	"""
	if not _with_lxml:
		warn("Install lxml for faster SBML I/O", ImportWarning)
	if not use_fbc_package:
		if libsbml is None:
			raise ImportError("libSBML required to write non-fbc models")
		return

	# create xml
	xml = model_to_xml(model, **kwargs)
	write_args = {"encoding": "UTF-8", "xml_declaration": True}
	if _with_lxml:
		write_args["pretty_print"] = True
		write_args["pretty_print"] = True
	else:
		_indent_xml(xml)

	# write xml to file
	xmlfile = filename
	should_close = True
	if hasattr(filename, "write"):
		should_close = False
	elif filename.endswith(".gz"):
		xmlfile = GzipFile(filename, "wb")
	elif filename.endswith(".bz2"):
		xmlfile = BZ2File(filename, "wb")
	elif isinstance(filename, str):
		if not filename.endswith(".xml") and not filename.endswith(".sbml"):
				filename += ".xml"
				xmlfile = filename
		xmlfile = open(filename, "wb")

	ElementTree(xml).write(xmlfile, **write_args)
	if should_close:
		xmlfile.close()

def read_sbml_model(filename):
	"""Load a mass model from a file in SBML/XML format.

	Parameters
	----------
	filename : str or file-like
		File path or descriptor that contains the SBML/XML document describing
		the mass model.

	Returns
	-------
	mass.MassModel
		The mass model as represented in the SBML/XML document.

	See Also
	--------
	parse_xml_into_model : Load from an xml object.
	"""
	return parse_xml_into_model(_parse_stream(filename))

# Internal Methods
def _parse_stream(filename):
	"""Parses filename or compressed stream to xml"""
	try:
		if hasattr(filename, "read"):
			return parse(filename)
		elif filename.endswith(".gz"):
			with GzipFile(filename) as infile:
				return parse(infile)
		elif filename.endswith(".bz2"):
			with BZ2File(filename) as infile:
				return parse(infile)
		else:
			return parse(filename)
	except ParseError as e:
		raise MassSBMLError("Malformed XML file: " + str(e))

def _get_attr(tag, attribute, type=lambda x: x, require=False):
	value = tag.get(ns(attribute))
	if require and value is None:
		msg = "required attribute '%s' not found in tag '%s'" % \
			  (attribute, tag.tag)
		if tag.get("id") is not None:
			msg += " with id '%s'" % tag.get("id")
		elif tag.get("name") is not None:
			msg += " with name '%s'" % tag.get("name")
		raise MassSBMLError(msg)
	return type(value) if value is not None else None

def _set_attrib(xml, attribute_name, value):
	if value is None or value == "":
		return
	xml.set(ns(attribute_name), str(value))

def _construct_gpr_xml(parent, expression):
	"""create gpr xml under parent node"""
	if isinstance(expression, BoolOp):
		op = expression.op
		if isinstance(op, And):
			new_parent = SubElement(parent, AND_TAG)
		elif isinstance(op, Or):
			new_parent = SubElement(parent, OR_TAG)
		else:
			raise Exception("unsupported operation " + op.__class__)
		for arg in expression.values:
			_construct_gpr_xml(new_parent, arg)
	elif isinstance(expression, Name):
		gene_elem = SubElement(parent, GENEREF_TAG)
		_set_attrib(gene_elem, "fbc:geneProduct", "G_" + expression.id)
	else:
		raise Exception("unsupported operation  " + repr(expression))

def _annotate_mass_from_sbml(mass_element, sbml_element):
	sbo_term = sbml_element.get("sboTerm")
	if sbo_term is not None:
		mass_element.annotation["SBO"] = sbo_term
	meta_id = _get_attr(sbml_element, "metaid")
	if meta_id is None:
		return
	annotation = mass_element.annotation
	for uri in extract_rdf_annotation(sbml_element, metaid="#" + meta_id):
		if not uri.startswith("http://identifiers.org/"):
			warn("%s does not start with http://identifiers.org/" % uri)
			continue
		try:
			provider, identifier = uri[23:].split("/", 1)
		except ValueError:
			warn("%s does not conform to http://identifiers.org/provider/id"
				 % uri)
			continue
		# handle multiple id's in the same database
		if provider in annotation:
			# make into a list if necessary
			if isinstance(annotation[provider], string_types):
				annotation[provider] = [annotation[provider]]
			annotation[provider].append(identifier)
		else:
			mass_element.annotation[provider] = identifier

def _annotate_sbml_from_mass(sbml_element, mass_element):
	if len(mass_element.annotation) == 0:
		return
	# get the id so we can set the metaid
	tag = sbml_element.tag
	if tag.startswith(sbml_prefix) or tag[0] != "{":
		prefix = ""
	elif tag.startswith(fbc_prefix):
		prefix = fbc_prefix
	else:
		raise ValueError("Can not annotate " + repr(sbml_element))
	_id = sbml_element.get(prefix + "id")
	if len(_id) == 0:
		raise ValueError("%s does not have id set" % repr(sbml_element))
	_set_attrib(sbml_element, "metaid", _id)
	annotation = SubElement(sbml_element, ns("sbml:annotation"))
	rdf_desc = SubElement(SubElement(annotation, ns("rdf:RDF")),
						  ns("rdf:Description"))
	_set_attrib(rdf_desc, "rdf:about", "#" + _id)
	bag = SubElement(SubElement(rdf_desc, ns("bqbiol:is")),
					 ns("rdf:Bag"))
	for provider, identifiers in sorted(iteritems(mass_element.annotation)):
		if provider == "SBO":
			_set_attrib(sbml_element, "sboTerm", identifiers)
			continue
		if isinstance(identifiers, string_types):
			identifiers = (identifiers,)
		for identifier in identifiers:
			li = SubElement(bag, ns("rdf:li"))
			_set_attrib(li, "rdf:resource", "http://identifiers.org/%s/%s" %
					   (provider, identifier))

# inspired by http://effbot.org/zone/element-lib.htm#prettyprint
def _indent_xml(elem, level=0):
	"""indent xml for pretty printing"""
	i = "\n" + level * "  "
	if len(elem):
		if not elem.text or not elem.text.strip():
			elem.text = i + "  "
		if not elem.tail or not elem.tail.strip():
			elem.tail = i
		for elem in elem:
			_indent_xml(elem, level + 1)
		if not elem.tail or not elem.tail.strip():
			elem.tail = i
	else:
		if level and (not elem.tail or not elem.tail.strip()):
			elem.tail = i

def _get_scale(string):

	scale = None

	for key in iterkeys(_SI_prefix_dict):
		if key in string:
			scale = _SI_prefix_dict[key]

	if scale is None:
		scale = 0

	return scale

# string utility methods
def _clip(string, prefix):
	"""_clips a prefix from the beginning of a string if it exists

	>>> _clip("R_pgi", "R_")
	"pgi"

	"""
	return string[len(prefix):] if string.startswith(prefix) else string

def _strnum(number):
	"""Utility function to convert a number to a string"""
	if isinstance(number, (Decimal, Basic, str)):
		return str(number)
	s = "%.15g" % number
	return s.rstrip(".")

def _remove_trailing_zeroes(string):
	# Remove whitespace
	string = string.replace(" ", "")
	original = string

	# Split by decimal
	str_list = string.split(".")

	#if len of 1, return original
	if len(str_list) == 1:
		return original

	# Remove leading zeroes
	str_list = [s.lstrip("0") for s in str_list]

	# If any string starts with a number, return original
	for s in str_list:
		if s[0].isdigit():
			return original

	# Else, join strings and return them
	result = ""
	for s in str_list:
		result += s
	#remove any whitespace
	result = result.replace(" ", "")
	return result

# Internal Variables
_SI_prefix_dict = {
	"atto": -18,
	"femto": -15,
	"pico": -12,
	"nano": -9,
	"micro": -6,
	"milli": -3,
	"centi": -2,
	"deci": -1,
	"deca": 1,
	"hecto": 2,
	"kilo": 3,
	"mega": 6,
	"giga": 9,
	"tera": 12,
	"peta": 15,
	"exa": 18
}

# litre <-> liter, metre <-> meter
_SBML_base_units = [
	"ampere",
	"avogadro",
	"becquerel",
	"candela",
	"coulomb",
	"dimensionless",
	"farad",
	"gram",
	"gray",
	"henry",
	"hertz",
	"item",
	"joule",
	"katal",
	"kelvin",
	"kilogram",
	"litre",
	"lumen",
	"lux",
	"metre",
	"mole",
	"newton",
	"ohm",
	"pascal",
	"radian",
	"second",
	"siemens",
	"sievert",
	"steradian",
	"tesla",
	"volt",
	"watt",
	"weber"

]
