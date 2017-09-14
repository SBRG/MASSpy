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
from math import inf

from six import iteritems, string_types

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

# deal with sbml2 here
# use sbml level 2 from sbml.py (which uses libsbml). Eventually, it would
# be nice to use the libSBML converters directly instead.
try:
    import libsbml
except ImportError:
    raise MassSBMLError("Need to install libsbml to proceed")
else:
    from cobra.io.sbml import create_cobra_model_from_sbml_file as read_sbml2
    from cobra.io.sbml import write_cobra_model_to_sbml_file as write_sbml2

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
BOUND_XPATH = ns("sbml:listOfParameters/sbml:parameter[@value]")
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
            yield get_attrib(i, "rdf:resource")
        for i in sbml_element.iterfind(search_xpath.replace(
                "isEncodedBy", "is")):
            yield get_attrib(i, "rdf:resource")

class MassSBMLError(Exception):
    pass



def get_attrib(tag, attribute, type=lambda x: x, require=False):
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



def set_attrib(xml, attribute_name, value):
    if value is None or value == "":
        return
    xml.set(ns(attribute_name), str(value))



def parse_stream(filename):
    """parses filename or compressed stream to xml"""
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



# string utility functions
def clip(string, prefix):
    """clips a prefix from the beginning of a string if it exists

    >>> clip("R_pgi", "R_")
    "pgi"

    """
    return string[len(prefix):] if string.startswith(prefix) else string


def strnum(number):
    """Utility function to convert a number to a string"""
    if isinstance(number, (Decimal, Basic, str)):
        return str(number)
    s = "%.15g" % number
    return s.rstrip(".")




def construct_gpr_xml(parent, expression):
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
            construct_gpr_xml(new_parent, arg)
    elif isinstance(expression, Name):
        gene_elem = SubElement(parent, GENEREF_TAG)
        set_attrib(gene_elem, "fbc:geneProduct", "G_" + expression.id)
    else:
        raise Exception("unsupported operation  " + repr(expression))



def annotate_mass_from_sbml(mass_element, sbml_element):
    sbo_term = sbml_element.get("sboTerm")
    if sbo_term is not None:
        mass_element.annotation["SBO"] = sbo_term
    meta_id = get_attrib(sbml_element, "metaid")
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



def annotate_sbml_from_mass(sbml_element, mass_element):
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
    id = sbml_element.get(prefix + "id")
    if len(id) == 0:
        raise ValueError("%s does not have id set" % repr(sbml_element))
    set_attrib(sbml_element, "metaid", id)
    annotation = SubElement(sbml_element, ns("sbml:annotation"))
    rdf_desc = SubElement(SubElement(annotation, ns("rdf:RDF")),
                          ns("rdf:Description"))
    set_attrib(rdf_desc, "rdf:about", "#" + id)
    bag = SubElement(SubElement(rdf_desc, ns("bqbiol:is")),
                     ns("rdf:Bag"))
    for provider, identifiers in sorted(iteritems(mass_element.annotation)):
        if provider == "SBO":
            set_attrib(sbml_element, "sboTerm", identifiers)
            continue
        if isinstance(identifiers, string_types):
            identifiers = (identifiers,)
        for identifier in identifiers:
            li = SubElement(bag, ns("rdf:li"))
            set_attrib(li, "rdf:resource", "http://identifiers.org/%s/%s" %
                       (provider, identifier))



def parse_xml_into_model(xml, number=float, importRatesAsCustom=False):
    # add model
    xml_model = xml.find(ns("sbml:model"))
    if get_attrib(xml_model, "fbc:required") == "false":
        warn('loading SBML model with fbc:required="false"')

    model_id = get_attrib(xml_model, "id")
    model = MassModel(model_id)
    model.name = xml_model.get("name")

    # add compartments
    model.compartments = {c.get("id"): c.get("name") for c in
                          xml_model.findall(COMPARTMENT_XPATH)}
    # Detect fixed concentration metabolites (boundary metabolites)
    boundary_metabolites = {clip(i.get("id"), "M_")
                            for i in xml_model.findall(SPECIES_XPATH % 'true')}
    bound_dict = {}
    # add metabolites (species)
    for species in xml_model.findall(ns("sbml:listOfSpecies/sbml:species")):
        met = get_attrib(species, "id", require=True)
        met = MassMetabolite(clip(met, "M_"))
        met.name = species.get("name")
        annotate_mass_from_sbml(met, species)
        met.compartment = species.get("compartment")
        met.charge = get_attrib(species, "fbc:charge", float)
        met.formula = get_attrib(species, "fbc:chemicalFormula")
        model.add_metabolites([met])
        model.update_initial_conditions(
            {met: get_attrib(species, "initialConcentration", float)})
        if met.id in boundary_metabolites:
            bound_dict[met] = get_attrib(species, "initialConcentration", 
                                         float)
    # add fixed concentration metabolites (boundary metabolites)
    model.add_fixed_concentrations(bound_dict)


    # add genes
    for sbml_gene in xml_model.iterfind(GENES_XPATH):
        gene_id = get_attrib(sbml_gene, "fbc:id").replace(SBML_DOT, ".")
        gene = Gene(clip(gene_id, "G_"))
        gene.name = get_attrib(sbml_gene, "fbc:name")
        if gene.name is None:
            gene.name = get_attrib(sbml_gene, "fbc:label")
        annotate_cobra_from_sbml(gene, sbml_gene)
        model.genes.append(gene)

    def process_gpr(sub_xml):
        """recursively convert gpr xml to a gpr string"""
        if sub_xml.tag == OR_TAG:
            return "( " + ' or '.join(process_gpr(i) for i in sub_xml) + " )"
        elif sub_xml.tag == AND_TAG:
            return "( " + ' and '.join(process_gpr(i) for i in sub_xml) + " )"
        elif sub_xml.tag == GENEREF_TAG:
            gene_id = get_attrib(sub_xml, "fbc:geneProduct", require=True)
            return clip(gene_id, "G_")
        else:
            raise Exception("unsupported tag " + sub_xml.tag)
    
    # add reactions
    custom_rate_dict = {}
    custom_param_dict = {}
    reactions = []
    for sbml_reaction in xml_model.iterfind(
            ns("sbml:listOfReactions/sbml:reaction")):
        reaction = get_attrib(sbml_reaction, "id", require=True)
        isReversible = get_attrib(sbml_reaction, "reversible")
        if isReversible == "false":
            isReversible = False
        else:
            isReversible = True
        reaction = MassReaction(clip(reaction, "R_"), reversible=isReversible)
        reaction.name = sbml_reaction.get("name")
        annotate_mass_from_sbml(reaction, sbml_reaction)
        reaction.subsystem = get_attrib(sbml_reaction, "compartment")

        # add mathml rate law extraction here
        result_from_mathml = ""
        for local_rate in sbml_reaction.findall(
            ns("sbml:kineticLaw/mml:math")):
            rate_xml_string = tostring(local_rate, encoding="utf-8").decode(
                "utf-8")
            ast = libsbml.readMathMLFromString(rate_xml_string)
            result_from_mathml = libsbml.formulaToL3String(ast)
            if result_from_mathml is None:
                msg = str(reaction.id) + " has imparsible MathML"
                warn(msg)
                continue
            else:
                result_from_mathml = result_from_mathml.replace("^", "**")

        # add rate constants, ssflux, custom paramaters (local parameters)
        custom_param_ids = []
        custom_param_values = []
        for local_parameter in sbml_reaction.findall(ns(
            "sbml:kineticLaw/sbml:listOfLocalParameters/sbml:localParameter")):
            pid = local_parameter.get("id")
            if pid == "ssflux_"+reaction.id:
                reaction.ssflux = number(local_parameter.get("value"))
            elif pid == "kf_"+reaction.id:
                reaction.kf = number(local_parameter.get("value"))
            elif pid == "Keq_"+reaction.id and reaction.reversible is True:
                if local_parameter.get("value") == "inf":
                    reaction.Keq = inf
                elif local_parameter.get("value") == "-inf":
                    reaction.Keq = -inf
                else:
                    reaction.Keq = number(local_parameter.get("value"))
            elif pid == "kr_"+reaction.id and reaction.reversible is True:
                reaction.kr = local_parameter.get("value")
            else:
                num_val = number(local_parameter.get("value"))
                custom_param_ids.append(pid)
                custom_param_values.append(num_val)
        custom_param_dict[reaction] = dict(zip(custom_param_ids, 
                                                custom_param_values))



        reactions.append(reaction)
        # add stoichiometry (reactants, products, modifiers)
        stoichiometry = defaultdict(lambda: 0)
        for species_reference in sbml_reaction.findall(
                ns("sbml:listOfReactants/sbml:speciesReference")):
            met_name = clip(species_reference.get("species"), "M_")
            stoichiometry[met_name] -= \
                float(species_reference.get("stoichiometry"))
        for species_reference in sbml_reaction.findall(
                ns("sbml:listOfProducts/sbml:speciesReference")):
            met_name = clip(species_reference.get("species"), "M_")
            stoichiometry[met_name] += \
                get_attrib(species_reference, "stoichiometry",
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

        test_dict = dict()
        test_dict[reaction] = reaction.rate_expression
        test_dict = strip_time(test_dict)

        if result_from_mathml is None:
            pass
        elif result_from_mathml.replace(" ", "") != str(
            test_dict[reaction]).replace(" ", ""):
            custom_rate_dict[reaction] = result_from_mathml



    try:
        model.add_reactions(reactions)
    except ValueError as e:
        warn(str(e))


    for custom_rxn, custom_rate in iteritems(custom_rate_dict):
        model.add_custom_rate(
            custom_rxn, custom_rate, custom_param_dict[custom_rxn])

    return model



def model_to_xml(mass_model):
    # add in model
    xml = Element("sbml", xmlns=namespaces["sbml"], level="3", version="1",
                  sboTerm="SBO:0000624")
    set_attrib(xml, "fbc:required", "false")
    xml_model = SubElement(xml, "model")
    set_attrib(xml_model, "fbc:strict", "true")
    if mass_model.id is not None:
        xml_model.set("id", mass_model.id)
    if mass_model.name is not None:
        xml_model.set("name", mass_model.name)

    # # add in parameters (fixed concentration metabolites)
    # parameter_list = SubElement(xml_model, "listOfParameters")
    # param_attr = {"constant": "true", "sboTerm": "SBO:0000285"}
    # for k, v in iteritems(mass_model.fixed_concentrations):
    #     param = SubElement(parameter_list, "parameter", id="M_" + k.id)
    #     set_attrib(param, "name", k.name)
    #     annotate_sbml_from_mass(param, k)
    #     set_attrib(param, "value", v)
    #     set_attrib(param, "constant", "true")
    #     set_attrib(param, "sboTerm", "SBO:0000285")

    # add in compartments
    compartments_list = SubElement(xml_model, "listOfCompartments")
    compartments = mass_model.compartments
    for compartment, name in iteritems(compartments):
        SubElement(compartments_list, "compartment", id=compartment, name=name,
                   constant="true")

    # add in species (metabolites)
    species_list = SubElement(xml_model, "listOfSpecies")
    for met in mass_model.metabolites:
        species = SubElement(species_list, "species", 
                             id="M_" + met.id,
                             # Useless required SBML parameter
                             hasOnlySubstanceUnits="false")
        set_attrib(species, "name", met.name)
        annotate_sbml_from_mass(species, met)
        set_attrib(species, "compartment", met.compartment)
        set_attrib(species, "fbc:charge", met.charge)
        set_attrib(species, "fbc:chemicalFormula", met.formula)
        set_attrib(species, "initialConcentration", 
                   mass_model.initial_conditions[met])
        # differentiate fixed concentration metabolites
        if met in mass_model.fixed_concentrations:
            set_attrib(species, "constant", "true")
            set_attrib(species, "boundaryCondition", "true")
        else:
            set_attrib(species, "constant", "false")
            set_attrib(species, "boundaryCondition", "false")

    # add in genes
    if len(mass_model.genes) > 0:
        genes_list = SubElement(xml_model, GENELIST_TAG)
        for gene in mass_model.genes:
            gene_id = gene.id.replace(".", SBML_DOT)
            sbml_gene = SubElement(genes_list, GENE_TAG)
            set_attrib(sbml_gene, "fbc:id", "G_" + gene_id)
            name = gene.name
            if name is None or len(name) == 0:
                name = gene.id
            set_attrib(sbml_gene, "fbc:label", gene_id)
            set_attrib(sbml_gene, "fbc:name", name)
            annotate_sbml_from_cobra(sbml_gene, gene)

    # add in reactions
    reactions_list = SubElement(xml_model, "listOfReactions")
    rate_dictionary = strip_time(mass_model.rate_expressions)
    for reaction, rate in iteritems(rate_dictionary):
        id = "R_" + reaction.id
        sbml_reaction = SubElement(reactions_list, "reaction", 
                                   id=id, 
                                   reversible=str(reaction.reversible).lower(),
                                   # Useless required SBML parameter
                                   fast="false")
        set_attrib(sbml_reaction, "name", reaction.name)
        set_attrib(sbml_reaction, "compartment", reaction.subsystem)
        annotate_sbml_from_mass(sbml_reaction, reaction)
        # add in kinetic law (rate law + associated constants)
        sbml_kinetic_law = SubElement(sbml_reaction, "kineticLaw")
        # add in math (rate law in mathml)
        rate_str = str(rate)
        if "**" in rate_str:
            rate_str = rate_str.replace("**", "^")
        new_mathml = libsbml.parseL3Formula(rate_str)
        new_str = libsbml.writeMathMLToString(new_mathml)
        if new_str is None:
            msg = "Error: Can't parse rate law expression for "
            raise MassSBMLError(msg+reaction.name+" with rate law: "+str(rate))
        else:
            new_str = new_str.replace(
                '<?xml version="1.0" encoding="UTF-8"?>\n', "")
            sbml_kinetic_law.append(fromstring(new_str))
        
        # add in local parameters (kf, Keq, kr, ssflux)
        sbml_param_list = SubElement(sbml_kinetic_law, "listOfLocalParameters")
        fwd_rate = "kf_"+reaction.id
        eq_const = "Keq_"+reaction.id
        rev_rate = "kr_"+reaction.id
        rxn_ssflux = "ssflux_"+reaction.id
        sbml_kf = SubElement(sbml_param_list, "localParameter", 
                             id=fwd_rate, name=fwd_rate)
        sbml_Keq = SubElement(sbml_param_list, "localParameter", 
                              id=eq_const, name=eq_const)
        sbml_kr = SubElement(sbml_param_list, "localParameter", 
                             id=rev_rate, name=rev_rate)
        sbml_ssflux = SubElement(sbml_param_list, "localParameter", 
                                 id=rxn_ssflux, name=rxn_ssflux)
        set_attrib(sbml_kf, "value", reaction.kf)
        set_attrib(sbml_Keq, "value", reaction.Keq)
        set_attrib(sbml_kr, "value", reaction.kr)
        set_attrib(sbml_ssflux, "value", reaction.ssflux)

        # add in local parameters (custom parameters)
        c_param_list = list(mass_model.custom_parameters.keys())
        try:
            c_rate = mass_model.custom_rates[reaction]
        except KeyError:
            pass
        else:
            for sym in c_rate.atoms(Symbol):
                if sym is Symbol("t"):
                    continue
                elif str(sym) in c_param_list:
                    v = mass_model.custom_parameters[str(sym)]
                    sbml_cparam = SubElement(sbml_param_list, "localParameter", 
                                             id=str(sym), name=str(sym))
                    set_attrib(sbml_cparam, "value", v)

        # add in reactants, products, modifiers (stoichiometry)
        reactants = {}
        products = {}
        for metabolite, stoichiomety in iteritems(reaction._metabolites):
            met_id = "M_" + metabolite.id
            if stoichiomety > 0:
                products[met_id] = strnum(stoichiomety)
            else:
                reactants[met_id] = strnum(-stoichiomety)
        if len(reactants) > 0:
            reactant_list = SubElement(sbml_reaction, "listOfReactants")
            for met_id, stoichiomety in sorted(iteritems(reactants)):
                SubElement(reactant_list, "speciesReference", species=met_id,
                           stoichiometry=stoichiomety, constant="true")
        if len(products) > 0:
            product_list = SubElement(sbml_reaction, "listOfProducts")
            for met_id, stoichiomety in sorted(iteritems(products)):
                SubElement(product_list, "speciesReference", species=met_id,
                           stoichiometry=stoichiomety, constant="true")

        # add in gene product association (gene reaction rule)
        gpr = reaction.gene_reaction_rule
        if gpr is not None and len(gpr) > 0:
            gpr = gpr.replace(".", SBML_DOT)
            gpr_xml = SubElement(sbml_reaction, GPR_TAG)
            try:
                parsed, _ = parse_gpr(gpr)
                construct_gpr_xml(gpr_xml, parsed.body)
            except Exception as e:
                print("failed on '%s' in %s" %
                      (reaction.gene_reaction_rule, repr(reaction)))
                raise e

    return xml



def write_sbml_model(mass_model, filename, use_fbc_package=True, **kwargs):
    if not _with_lxml:
        warn("Install lxml for faster SBML I/O", ImportWarning)
    if not use_fbc_package:
        if libsbml is None:
            raise ImportError("libSBML required to write non-fbc models")
        # ignore for now (deal with after sbml3 is finished)
        #write_sbml2(cobra_model, filename, use_fbc_package=False, **kwargs)
        return
    # create xml
    xml = model_to_xml(mass_model, **kwargs)
    write_args = {"encoding": "UTF-8", "xml_declaration": True}
    if _with_lxml:
        write_args["pretty_print"] = True
        write_args["pretty_print"] = True
    else:
        indent_xml(xml)
    # write xml to file
    should_close = True
    if hasattr(filename, "write"):
        xmlfile = filename
        should_close = False
    elif filename.endswith(".gz"):
        xmlfile = GzipFile(filename, "wb")
    elif filename.endswith(".bz2"):
        xmlfile = BZ2File(filename, "wb")
    else:
        xmlfile = open(filename, "wb")
    ElementTree(xml).write(xmlfile, **write_args)
    if should_close:
        xmlfile.close()



# inspired by http://effbot.org/zone/element-lib.htm#prettyprint
def indent_xml(elem, level=0):
    """indent xml for pretty printing"""
    i = "\n" + level * "  "
    if len(elem):
        if not elem.text or not elem.text.strip():
            elem.text = i + "  "
        if not elem.tail or not elem.tail.strip():
            elem.tail = i
        for elem in elem:
            indent_xml(elem, level + 1)
        if not elem.tail or not elem.tail.strip():
            elem.tail = i
    else:
        if level and (not elem.tail or not elem.tail.strip()):
            elem.tail = i





