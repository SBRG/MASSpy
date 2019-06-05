# -*- coding: utf-8 -*-
"""TODO Module Docstrings."""
from __future__ import absolute_import

import datetime
import logging
import operator
import re

import libsbml

from six import integer_types, iteritems, string_types

from sympy import Basic, Symbol, mathml

from cobra.io.sbml import (
    BOUND_MINUS_INF, BOUND_PLUS_INF, LOWER_BOUND_ID, SBO_DEFAULT_FLUX_BOUND,
    SBO_FLUX_BOUND, UPPER_BOUND_ID, ZERO_BOUND_ID, _create_bound,
    _sbase_annotations, _sbase_notes_dict)

from mass.core import (
    MassConfiguration, MassMetabolite, MassModel, MassReaction, Unit,
    UnitDefinition, _SBML_BASE_UNIT_KINDS_DICT)
from mass.enzyme_modules import (
    EnzymeModule, EnzymeModuleDict, EnzymeModuleForm, EnzymeModuleReaction,
    _ORDERED_ENZYMEMODULE_DICT_DEFAULTS, _make_bound_attr_str_repr)
from mass.exceptions import MassSBMLError
from mass.util import get_subclass_specific_attributes, strip_time

LOGGER = logging.getLogger(__name__)


# -----------------------------------------------------------------------------
# Defaults and constants for writing SBML
# -----------------------------------------------------------------------------
MASSCONFIGURATION = MassConfiguration()
# For SBML level and core version, and package extension versions
SBML_LEVEL_VERSION = (3, 2)
FBC_VERSION = 2
GROUPS_VERSION = 1

# Default compartment and boundary compartment definitions
DEFAULT_COMPARTMENT_DICT = MASSCONFIGURATION.default_compartment
BOUNDARY_COMPARTMENT_DICT = MASSCONFIGURATION.boundary_compartment

# Precompiled regex for mass and SBML moieties
MASS_MOIETY_RE = re.compile("\[\S+\]")
SBML_MOIETY_RE = re.compile("Moiety\S+$")

# For SBO terms
SBO_MODELING_FRAMEWORK = "SBO:0000062"

# COBRA Flux units
COBRA_UNIT_FLUX = UnitDefinition(
    id="mmol_per_gDW_per_hr", name="cobra_flux_unit", list_of_units=[
        Unit(kind="mole", exponent=1, scale=-3, multiplier=1),
        Unit(kind="gram", exponent=-1, scale=0, multiplier=1),
        Unit(kind="second", exponent=-1, scale=0, multiplier=3600),
    ]
)
# -----------------------------------------------------------------------------
# Functions for replacements (import/export)
# -----------------------------------------------------------------------------
SBML_DOT = "__SBML_DOT__"


def _clip(sid, prefix):
    """Clips a prefix from the beginning of a string if it exists."""
    return sid[len(prefix):] if sid.startswith(prefix) else sid


def _f_gene(sid, prefix="G_"):
    """Clip gene prefix from id."""
    sid = sid.replace(SBML_DOT, ".")
    return _clip(sid, prefix)


def _f_gene_rev(sid, prefix="G_"):
    """Add gene prefix to id."""
    return prefix + sid.replace(".", SBML_DOT)


def _f_specie(sid, prefix="M_"):
    """Clip specie/metabolite prefix from id."""
    return _clip(sid, prefix)


def _f_specie_rev(sid, prefix="M_"):
    """Add specie/metabolite prefix to id."""
    return prefix + sid


def _f_reaction(sid, prefix="R_"):
    """Clip reaction prefix from id."""
    return _clip(sid, prefix)


def _f_reaction_rev(sid, prefix="R_"):
    """Add reaction prefix to id."""
    return prefix + sid


F_GENE = "F_GENE"
F_GENE_REV = "F_GENE_REV"
F_SPECIE = "F_SPECIE"
F_SPECIE_REV = "F_SPECIE_REV"
F_REACTION = "F_REACTION"
F_REACTION_REV = "F_REACTION_REV"

F_REPLACE = {
    F_GENE: _f_gene,
    F_GENE_REV: _f_gene_rev,
    F_SPECIE: _f_specie,
    F_SPECIE_REV: _f_specie_rev,
    F_REACTION: _f_reaction,
    F_REACTION_REV: _f_reaction_rev,
}


def _f_moiety_formula(metabolite):
    """Fix moieties in formula from SBML compatible to be masspy compatible."""
    print("TODO _f_moiety_formula ", metabolite)
    return


def _f_moiety_formula_rev(metabolite):
    """Fix moieties in formula to be SBML compatible from masspy compatible."""
    sbml_formula = metabolite.formula
    moieties = list(filter(MASS_MOIETY_RE.match, metabolite.elements))
    if moieties:
        for moiety in moieties:
            replacement = "Moiety" + moiety.lstrip("[").rstrip("]").lower()
            # Add to end of formula prevents issues if moiety ends with number
            sbml_formula = sbml_formula.replace(moiety, "")
            sbml_formula += replacement
        # Replace any dashes
        sbml_formula = sbml_formula.replace("-", "")

    return sbml_formula


def _for_id(sid):
    """Return a string specifying the object id for logger messages."""
    return " for '{0}'".format(sid)


# -----------------------------------------------------------------------------
# MathML
# -----------------------------------------------------------------------------
# For MathML representation of kinetic laws and other sympy equations
MATHML_MML_TAG_RE = re.compile("\<mml:.*?\>|\<\/mml:.*?\>")
MATH_XML_FMT = '<math xmlns="http://www.w3.org/1998/Math/MathML">{0}</math>'


def _create_math_xml_str_from_sympy_expr(sympy_equation):
    """Create a MathML string from a sympy expression to be parsed by libsbml.

    This also requires removing 'mml' tags from XML representation of
        mathemtical equations.

    Sympy expressions can be converted into their equivalent MathML string
    through the `sympy.printing.mathml` module. However, the mathml module will
    interpret '_' as a subscript and '__' as a superscript. In order to prevent
    the module from misinterpreting underscores in identifiers, the underscores
    are converted into ampsersands. Additionally, all MathML presentation
    markup must be removed for similar reasons. After the XML string is made,
    the identifiers are returned to their original state.
    """
    underscore_replace = {str(arg): str(arg).replace("_", "&")
                          for arg in sympy_equation.atoms(Symbol)}
    math_xml_str = mathml(sympy_equation.subs(underscore_replace))
    math_xml_str = MATH_XML_FMT.format(math_xml_str.replace("&", "_"))
    math_xml_str = MATHML_MML_TAG_RE.sub("", math_xml_str)

    return math_xml_str

def _create_sympy_expr_from_math_xml_str(mathml_xml_str):
    """Create a sympy expression from a MathML string to be parsed by libsbml.
    """
# -----------------------------------------------------------------------------
# Read SBML
# -----------------------------------------------------------------------------
def read_sbml_model(filename):
    """TODO DOCSTRING."""
    print(filename)


def _parse_notes_dict(sbase):
    """ Create dictionary of MASS notes.

    Parameters
    ----------
    sbase : libsbml.SBase

    Returns
    -------
    dict of notes

    """
    mass_notes = {}
    notes = sbase.getNotesString()
    if notes:
        start_pos_list = (match.end() for match in re.finditer(r"<p>", notes))
        end_pos_list = (match.start() for match in re.finditer(r"</p>", notes))

        for start, end in zip(start_pos_list, end_pos_list):
            match = re.search(
                "\s*(\w+\s*\w*)\s*:\s*(...+)", notes[start:end], re.DOTALL)
            k, v = match.groups()
            mass_notes.update({k.strip(): v.strip()})

    return mass_notes


# -----------------------------------------------------------------------------
# Write SBML
# -----------------------------------------------------------------------------
def write_sbml_model(mass_model, filename, f_replace=F_REPLACE, **kwargs):
    """Write MassModel to filename in SBML format.

    The created model is SBML level 3 version 2 core (L3V2).

    If the given filename ends with the suffix ".gz" (for example,
    "myfile.xml.gz"), libSBML assumes the caller wants the file to be
    written compressed in gzip format. Similarly, if the given filename
    ends with ".zip" or ".bz2", libSBML assumes the caller wants the
    file to be compressed in zip or bzip2 format (respectively). Files
    whose names lack these suffixes will be written uncompressed. Special
    considerations for the zip format: If the given filename ends with
    ".zip", the file placed in the zip archive will have the suffix
    ".xml" or ".sbml".  For example, the file in the zip archive will
    be named "test.xml" if the given filename is "test.xml.zip" or
    "test.zip". Similarly, the filename in the archive will be
    "test.sbml" if the given filename is "test.sbml.zip".

    Parameters
    ----------
    mass_model : mass.core.MassModel
        MassModel instance which is written to SBML
    filename : string
        path to which the model is written
    f_replace: dict of replacement functions for id replacement
        A dict of replacement functions to apply on identifiers.
    units: bool
        Should the units be written into the SBMLDocument.
        Default is True.
    include_all: bool
        If True, all information in the model be written into the SBMLDocument,
        otherwise only the information required to simulate the model ODEs.
        Default is True.

    """
    doc = _model_to_sbml(mass_model, f_replace=f_replace, **kwargs)
    print()
    print(libsbml.writeSBMLToString(doc))
    if isinstance(filename, string_types):
        # Write to path
        libsbml.writeSBMLToFile(doc, filename)

    elif hasattr(filename, "write"):
        # Write to file handle
        sbml_str = libsbml.writeSBMLToString(doc)
        filename.write(sbml_str)


def _model_to_sbml(mass_model, f_replace=None, units=True, include_all=True):
    """Convert MassModel to SBMLDocument.

    Parameters
    ----------
    mass_model: mass.MassModel
        MassModel instance which is written to SBML
    f_replace: dict of replacement functions
        Replacement to apply on identifiers.
    units: bool
        Should the units be written into the SBMLDocument.
    include_all: bool
        Should all information in the model be written into the SBMLDocument,
        or only the information required to simulate the model ODEs.

    Returns
    -------
    doc: libsbml.SBMLDocument

    """
    if f_replace is None:
        f_replace = {}

    doc, models = _create_SBMLDocument_and_model_objects(mass_model,
                                                         include_all)
    model, model_fbc, model_group = models

    if include_all:
        # Write the model Meta Information (ModelHistory) into the SBMLDocument
        if hasattr(mass_model, "_sbml"):
            meta = mass_model._sbml
            _write_model_meta_info_to_sbml(model, meta)
            # Set meta annotation and notes
            if "annotation" in meta:
                _sbase_annotations(doc, meta["annotation"])
            if "notes" in meta:
                _sbase_notes_dict(doc, meta["notes"])

        # Write the units information into the SBMLDocument
        if units:
            _write_model_units_to_sbml(model, mass_model)

        # Write model notes and annotations into the SBMLDocument
        _sbase_notes_dict(model, mass_model.notes)
        _sbase_annotations(model, mass_model.annotation)

        # Write the flux upper/lower bound parameters into the SBMLDocument
        _write_model_flux_parameters_to_sbml(model, mass_model, units=units)

    # Write the rate law parameters into the SBMLDocument
    _write_model_parameters_to_sbml(
        model, mass_model, f_replace, units=units, include_all=include_all)

    # Write the compartment information into the SBMLDocument
    _write_model_compartments_to_sbml(model, mass_model)

    # Write the species information into the SBMLDocument
    # Includes MassMetabolites, EnzymeModuleForms, and their concentrations
    _write_model_species_to_sbml(
        model, mass_model, f_replace, include_all=include_all)
    # Write the boundary conditions into the SBMLDocument
    _write_model_boundary_conditions_to_sbml(model, mass_model, f_replace)
    # Write the gene information into the SBMLDocument
    if include_all:
        _write_model_genes_to_sbml(model_fbc, mass_model, f_replace)
    # Write the reaction information into the SBMLDocument
    # Includes MassReactions and EnzymeModuleReactions
    _write_model_reactions_to_sbml(
        model, mass_model, f_replace, units=units, include_all=include_all)

    if include_all:
        # Write the EnzymeModule specific information into the SBMLDocument
        if isinstance(mass_model, EnzymeModule):
            _write_enzyme_modules_to_sbml(model_group, mass_model, f_replace)
        # Write the EnzymeModuleDict specific information into the SBMLDocument
        if mass_model.enzyme_modules:
            for enzyme in mass_model.enzyme_modules:
                _write_enzyme_modules_to_sbml(model_group, enzyme, f_replace)

    return doc


def _create_SBMLDocument_and_model_objects(mass_model, include_all=None):
    """Create and return the SBMLDocument and model objects.

    Will also set package extensions, model ID, meta ID, and name

    Warnings
    --------
    This method is intended for internal use only.

    """
    # Try creating the SBMLDocument object
    sbml_ns = libsbml.SBMLNamespaces(*SBML_LEVEL_VERSION)  # SBML L3V2 Core
    try:
        doc = libsbml.SBMLDocument(sbml_ns)
    except ValueError as e:
        MassSBMLError("Could not create SBMLDocument due to:\n" + str(e))

    # Create model
    model = doc.createModel()
    _check(model, "create model")
    _check(doc.setSBOTerm(SBO_MODELING_FRAMEWORK),
           "set SBO term for modeling framework")

    sbml_plugin_dict = {"fbc": FBC_VERSION}
    if (isinstance(mass_model, EnzymeModule) or mass_model.enzyme_modules)\
       and include_all:
        sbml_plugin_dict.update({"groups": GROUPS_VERSION})

    plugin_models = []
    for plugin, version_info in iteritems(sbml_plugin_dict):
        # Enable package extension
        _check(
            doc.enablePackage(
                "http://www.sbml.org/sbml/level{0}/version1/{1}/version{2}"
                .format(SBML_LEVEL_VERSION[0], plugin, version_info),
                plugin, True),
            "enable package extension {0}-v{1}".format(plugin, version_info))
        # Set required to false
        _check(doc.setPackageRequired(plugin, False),
               "set {0} extension required to false".format(plugin))
        # Get plugin and return
        plugin_models += [model.getPlugin(plugin)]
        _check(plugin_models[-1], "get {0} plugin for model".format(plugin))

    if len(plugin_models) < 2:
        plugin_models += [None for i in range(2 - len(plugin_models))]
    model_fbc, model_group = plugin_models

    # Set strictness of fbc package
    _check(model_fbc.setStrict(False), "set fbc plugin strictness to false")

    # Set ID, meta ID, and name
    if mass_model.id is not None:
        _check(model.setId(mass_model.id), "set model id")
        _check(model.setMetaId("meta_" + mass_model.id), "set model meta id")
    else:
        _check(model.setMetaId("meta_model"), "set model meta id")
    if mass_model.name is not None:
        _check(model.setName(mass_model.name), "set model name")

    return doc, (model, model_fbc, model_group)


def _write_model_meta_info_to_sbml(model, meta):
    """Write MassModel Meta Information (ModelHistory) into SBMLDocument.

    Warnings
    --------
    This method is intended for internal use only.

    """
    # Make ModelHistory object and populate
    history = libsbml.ModelHistory()
    _check(history, "create ModelHistory")
    if meta.get("created", None):
        # Use previous creation date
        _check(history.setCreatedDate(meta["created"]), "set created date")
    else:
        # Get the current time and date and set for the model.
        time = datetime.datetime.now()
        timestr = time.strftime('%Y-%m-%dT%H:%M:%S')
        date = libsbml.Date(timestr)
        _check(date, "create the date")
        _check(history.setCreatedDate(date), "set created date")
        _check(history.setModifiedDate(date), "set modified date")

    # Make Creator objects and add to ModelHistory
    if meta.get("creators", None):
        for mass_creator in meta["creators"]:
            # Make creator object
            creator = libsbml.ModelCreator()
            _check(creator, "create model creator")
            # Add family name, given name, organisation, and email attributes
            for k in ["familyName", "givenName", "organisation", "email"]:
                if mass_creator.get(k, None):
                    # Make logger message
                    msg = k.replace("Name", " name") if "Name" in k else k
                    # Set creator attribute
                    _check(creator.__class__.__dict__[k].fset(
                        creator, mass_creator[k]), "set creator " + msg)
            # Add creator to the ModelHistory
            _check(history.addCreator(creator),
                   "adding creator to ModelHistory")
    # Set the ModelHistory
    _check(model.setModelHistory(history), 'set ModelHistory')


def _write_model_units_to_sbml(model, mass_model):
    """Write MassModel unit information into SBMLDocument.

    Currently, all units will be in terms of Millimole, Litres, and Hours.

    Warnings
    --------
    This method is intended for internal use only.

    """
    # Create nested function for making new SBML unit definitions
    def _write_unit_definition(id_str, name="", unit_def=None):
        """Create unit definition."""
        udef = model.createUnitDefinition()
        _check(udef, "create UnitDefinition" + _for_id(id_str))
        _check(udef.setId(id_str), "set UnitDefinition id" + _for_id(id_str))
        if name:
            _check(udef.setName(name),
                   "set UnitDefinition name" + _for_id(id_str))
        for u in unit_def:
            unit = udef.createUnit()
            _check(unit, "create Unit" + _for_id(id_str))
            _check(unit.setKind(_SBML_BASE_UNIT_KINDS_DICT[u.kind]),
                   "set Unit kind" + _for_id(id_str))
            _check(unit.setExponent(u.exponent),
                   "set Unit exponent" + _for_id(id_str))
            _check(unit.setScale(u.scale), "set Unit scale" + _for_id(id_str))
            _check(unit.setMultiplier(u.multiplier),
                   "set Unit multiplier" + _for_id(id_str))

    # Write MassModel units
    if mass_model.units:
        to_write = list(mass_model.units) + [COBRA_UNIT_FLUX]
    else:
        to_write = [COBRA_UNIT_FLUX]

    # Write flux bound units
    for unit_def in to_write:
        _write_unit_definition(
            id_str=unit_def.id, name=unit_def.name, unit_def=unit_def)


def _write_model_flux_parameters_to_sbml(model, mass_model, units=None):
    """Write MassModel flux bound information into SBMLDocument.

    This includes the default upper and lower bounds, positive and negative
    infinity values, and a flux bound of value 0.

    Warnings
    --------
    This method is intended for internal use only.

    """
    # Determine the minimum and maximum flux bound values
    if mass_model.reactions:
        min_value = min(mass_model.reactions.list_attr("lower_bound"))
        max_value = max(mass_model.reactions.list_attr("upper_bound"))
    else:
        # Use default values if no reactions in model.
        min_value = MASSCONFIGURATION.lower_bound
        max_value = MASSCONFIGURATION.upper_bound

    # Create the model flux parameters
    flux_udef = model.getUnitDefinition(COBRA_UNIT_FLUX.id)
    _create_parameter(
        model, pid=LOWER_BOUND_ID, value=min_value, sbo=SBO_DEFAULT_FLUX_BOUND,
        units=units, udef=flux_udef)
    _create_parameter(
        model, pid=UPPER_BOUND_ID, value=max_value, sbo=SBO_DEFAULT_FLUX_BOUND,
        units=units, udef=flux_udef)
    _create_parameter(
        model, pid=ZERO_BOUND_ID, value=0, sbo=SBO_DEFAULT_FLUX_BOUND,
        units=units, udef=flux_udef)
    _create_parameter(
        model, pid=BOUND_MINUS_INF, value=-float("Inf"), sbo=SBO_FLUX_BOUND,
        units=units, udef=flux_udef)
    _create_parameter(
        model, pid=BOUND_PLUS_INF, value=float("Inf"), sbo=SBO_FLUX_BOUND,
        units=units, udef=flux_udef)


def _write_model_parameters_to_sbml(model, mass_model, f_replace, units=None,
                                    include_all=None):
    """Write MassModel kinetic parameter information into SBMLDocument.

    The kinetic parameters include the reaction parameters kf, Keq, and kr, and
    any custom parameters in the model.

    Warnings
    --------
    This method is intended for internal use only.

    """
    for parameter_type, parameter_dict in iteritems(mass_model.parameters):
        # Boundary metabolites addressed later with metabolites
        if parameter_type == "Boundary":
            continue
        if parameter_type == "v":
            constant = False
        else:
            constant = True
        # Create kf, Keq, kr, and custom parameters, handling ID corrections
        # for recognized reaction IDs in the parameters
        for parameter_id, value in iteritems(parameter_dict):
            try:
                pid, rid = parameter_id.split("_", 1)
                rid = mass_model.reactions.get_by_id(rid).id
            except (ValueError, KeyError):
                pid = parameter_id
            else:
                if f_replace and F_REACTION_REV in f_replace:
                    rid = f_replace[F_REACTION_REV](rid)
                pid = "_".join((pid, rid))
            finally:
                _create_parameter(
                    model, pid=pid, value=value, sbo=None, constant=constant,
                    units=units, udef=None, include_all=include_all)


def _create_parameter(model, pid, value, sbo=None, constant=True, units=None,
                      udef=None, include_all=True):
    """Create a global model parameter to be written into the SBMLDocument.

    Warnings
    --------
    This method is intended for internal use only.

    """
    # Create parameter and set its ID, value, and whether it is constant
    parameter = model.createParameter()
    _check(parameter, "create model parameter" + _for_id(pid))
    _check(parameter.setId(pid), "set parameter id" + _for_id(pid))
    _check(parameter.setValue(value), "set parameter value" + _for_id(pid))
    _check(parameter.setConstant(constant),
           "set parameter constant" + _for_id(pid))
    if include_all:
        # Set SBO term and units if desired.
        if sbo:
            _check(parameter.setSBOTerm(sbo),
                   "set parameter sbo term" + _for_id(pid))
        if units and udef is not None:
            _check(parameter.setUnits(udef.getId()),
                   "set parameter units" + _for_id(pid))


def _write_model_compartments_to_sbml(model, mass_model):
    """Write MassModel compartment information into SBMLDocument.

    Warnings
    --------
    This method is intended for internal use only.

    """
    # TODO Revisit based on new cobra compartment implementation
    if not mass_model.compartments:
        LOGGER.warning(
            "No compartments found in model. Therefore creating compartment "
            "'%s' for entire model.", list(DEFAULT_COMPARTMENT_DICT)[0])
        compartment_dict = DEFAULT_COMPARTMENT_DICT
    else:
        compartment_dict = mass_model.compartments

    if mass_model.boundary_conditions:
        compartment_dict.update(BOUNDARY_COMPARTMENT_DICT)

    for cid, name in iteritems(compartment_dict):
        compartment = model.createCompartment()
        _check(compartment, "create model compartment" + _for_id(cid))
        _check(compartment.setId(cid), "set compartment id" + _for_id(cid))
        _check(compartment.setName(name),
               "set compartment name" + _for_id(cid))
        _check(compartment.setConstant(True),
               "set compartment constant" + _for_id(cid))


def _write_model_species_to_sbml(model, mass_model, f_replace,
                                 include_all=None):
    """Write MassModel species information into SBMLDocument.

    Warnings
    --------
    This method is intended for internal use only.

    """
    for mass_metabolite in mass_model.metabolites:
        # Create specie and set ID, name, and compartment
        mid = mass_metabolite.id
        if f_replace and F_SPECIE_REV in f_replace:
            mid = f_replace[F_SPECIE_REV](mid)
        specie = model.createSpecies()
        _check(specie, "create model specie" + _for_id(mid))
        _check(specie.setId(mid), "set specie id" + _for_id(mid))
        _check(specie.setName(mass_metabolite.name),
               "set specie name" + _for_id(mid))
        # Use a generic compartment if no species have set compartments
        if not mass_model.compartments:
            cid = list(DEFAULT_COMPARTMENT_DICT)[0]
        else:
            cid = mass_metabolite.compartment
        _check(specie.setCompartment(cid),
               "set specie compartment" + _for_id(mid))
        # Set metabolite initial condition value
        if mass_metabolite.initial_condition is not None:
            _check(
                specie.setInitialConcentration(
                    mass_metabolite.initial_condition),
                "set specie initial concentration" + _for_id(mid))

        # Set species constant, boundary condition, and unit restrictions
        _check(specie.setConstant(mass_metabolite.fixed),
               "set specie constant" + _for_id(mid))
        _check(specie.setBoundaryCondition(False),
               "set specie boundary condition" + _for_id(mid))
        _check(specie.setHasOnlySubstanceUnits(False),
               "set specie unit substance restrictions" + _for_id(mid))

        if include_all:
            # Set metabolite charge via fbc package
            specie_fbc = specie.getPlugin("fbc")
            _check(specie_fbc, "get fbc plugin for specie" + _for_id(mid))
            if mass_metabolite.charge is not None:
                _check(specie_fbc.setCharge(mass_metabolite.charge),
                       "set specie charge" + _for_id(mid))
            # Set metabolite formula via fbc package
            if mass_metabolite.formula is not None:
                sbml_formula = _f_moiety_formula_rev(mass_metabolite)
                _check(specie_fbc.setChemicalFormula(sbml_formula),
                       "set specie formula" + _for_id(mid))
            # Set metabolite annotation and notes
            if isinstance(mass_metabolite, EnzymeModuleForm):
                # Set enzyme information if EnzymeModuleForm
                _write_enzyme_attr_info_to_notes(
                    specie, mass_metabolite, f_replace=f_replace)
            else:
                # Otherwise just set the MassMetabolite notes regularly
                _sbase_notes_dict(specie, mass_metabolite.notes)
            # Set metabolite annotation
            _sbase_annotations(specie, mass_metabolite.annotation)


def _write_model_boundary_conditions_to_sbml(model, mass_model, f_replace):
    """Write MassModel boundary condition information into SBMLDocument.

    Warnings
    --------
    This method is intended for internal use only.

    """
    cid = str(list(BOUNDARY_COMPARTMENT_DICT)[0])
    for bmid, bc_value in iteritems(mass_model.boundary_conditions):
        # Create boundary specie and set ID, name, and compartment
        if f_replace and F_SPECIE_REV in f_replace:
            bmid = f_replace[F_SPECIE_REV](bmid)
        specie = model.createSpecies()
        _check(specie, "create model boundary specie" + _for_id(bmid))
        _check(specie.setId(bmid), "set boundary specie id" + _for_id(bmid))
        _check(specie.setCompartment(cid),
               "set boundary specie compartment" + _for_id(bmid))
        # Set boundary specie value and constant
        if isinstance(bc_value, (integer_types, float)):
            _check(specie.setInitialConcentration(bc_value),
                   "set boundary specie concentration" + _for_id(bmid))
            _check(specie.setConstant(True),
                   "set specie constant" + _for_id(bmid))
        else:
            # TODO handle functions of time for boundary conditions
            print("TODO handle functions of time for boundary conditions")
            _check(specie.setConstant(False),
                   "set specie constant" + _for_id(bmid))

        # Set species as boundary condition, and set unit restrictions
        _check(specie.setBoundaryCondition(True),
               "set specie boundary condition" + _for_id(bmid))
        _check(specie.setHasOnlySubstanceUnits(False),
               "set specie unit substance restrictions" + _for_id(bmid))


def _write_model_genes_to_sbml(model_fbc, mass_model, f_replace):
    """Write MassModel gene information into SBMLDocument.

    Warnings
    --------
    This method is intended for internal use only.

    """
    for mass_gene in mass_model.genes:
        # Create gene product and set ID
        gid = mass_gene.id
        if f_replace and F_GENE_REV in f_replace:
            gid = f_replace[F_GENE_REV](gid)
        gp = model_fbc.createGeneProduct()
        _check(gp, "create model gene product" + _for_id(gid))
        _check(gp.setId(gid), "set gene id " + _for_id(gid))
        # Set gene name and label
        gname = mass_gene.name if mass_gene.name else gid
        _check(gp.setName(gname), "set gene name" + _for_id(gid))
        _check(gp.setLabel(gname), "set gene label" + _for_id(gid))

        # Set gene annotation and notes
        _sbase_annotations(gp, mass_gene.annotation)
        _sbase_notes_dict(gp, mass_gene.notes)


def _write_model_reactions_to_sbml(model, mass_model, f_replace, units=None,
                                   include_all=None):
    """Write MassModel reaction information into SBMLDocument.

    Warnings
    --------
    This method is intended for internal use only.

    """
    for mass_reaction in mass_model.reactions:
        # Create reaction and set ID, name, and reversible
        rid = mass_reaction.id
        if f_replace and F_REACTION_REV in f_replace:
            rid = f_replace[F_REACTION_REV](rid)
        reaction = model.createReaction()
        _check(reaction, "create reaction" + _for_id(rid))
        _check(reaction.setId(rid), "set reaction id" + _for_id(rid))
        _check(reaction.setName(mass_reaction.name),
               "set reaction name" + _for_id(rid))
        _check(reaction.setReversible(mass_reaction.reversible),
               "set reaction reversible" + _for_id(rid))

        # Write specie references into reaction
        reaction = _write_reaction_specie_ref_to_sbml(
            reaction, mass_reaction=mass_reaction, f_replace=f_replace)

        # Set the reaction kinetic law
        _write_reaction_kinetic_law_to_sbml(
            reaction, mass_reaction=mass_reaction, f_replace=f_replace)

        if include_all:
            # Get plugin for reaction
            reaction_fbc = reaction.getPlugin("fbc")
            _check(reaction_fbc, "get fbc plugin for reaction" + _for_id(rid))

            # Set reaction flux bounds
            flux_udef = model.getUnitDefinition(COBRA_UNIT_FLUX.id)
            reaction_fbc.setLowerFluxBound(
                _create_bound(
                    model, reaction=mass_reaction, bound_type="lower_bound",
                    f_replace=f_replace, units=units, flux_udef=flux_udef))
            reaction_fbc.setUpperFluxBound(
                _create_bound(
                    model, reaction=mass_reaction, bound_type="upper_bound",
                    f_replace=f_replace, units=units, flux_udef=flux_udef))

            # Set the reaction GPR if it exists
            _write_reaction_gpr_to_sbml(reaction_fbc, mass_reaction, f_replace)

            if isinstance(mass_reaction, EnzymeModuleReaction):
                # Set enzyme information if EnzymeModuleReaction
                _write_enzyme_attr_info_to_notes(
                    reaction, mass_reaction, f_replace=f_replace)
            else:
                # Otherwise just set the MassReaction notes regularly
                _sbase_notes_dict(reaction, mass_reaction.notes)
            # Set reaction annotation
            _sbase_annotations(reaction, mass_reaction.annotation)


def _write_reaction_specie_ref_to_sbml(reaction, mass_reaction, f_replace):
    """Write MassReaction species reference information into SBMLDocument.

    Warnings
    --------
    This method is intended for internal use only.

    """
    # Helper function for creating specie references
    def _write_specie_ref_into_reaction(reaction, metabolite, stoich,
                                        comparison_op, constant, f_replace):
        """Write specie references into reactions."""
        # Get specie id
        sid = str(metabolite)
        if f_replace and F_SPECIE_REV in f_replace:
            sid = f_replace[F_SPECIE_REV](sid)
        # Comparison to determine whether metabolite is reactant or product
        if comparison_op(stoich, 0):
            sref = reaction.createReactant()
        else:
            sref = reaction.createProduct()
        _check(sref,
               "create specie reference " + sid + _for_id(reaction.getId()))
        # Set ID, stoichiometry, and whether species is constant
        _check(sref.setSpecies(sid), "set specie reference id" + _for_id(sid))
        _check(sref.setStoichiometry(abs(stoich)),
               "set specie reference stoichiometry" + _for_id(sid))
        _check(sref.setConstant(constant),
               "set specie reference constant" + _for_id(sid))

        return reaction

    # Set reaction species references and stoichiometry
    for metabolite, stoichiometry in iteritems(mass_reaction.metabolites):
        # Write metabolite into reaction
        reaction = _write_specie_ref_into_reaction(
            reaction=reaction, metabolite=metabolite, stoich=stoichiometry,
            comparison_op=operator.lt, constant=metabolite.fixed,
            f_replace=f_replace)

        # Set boundary specie reference if reaction on boundary
        bc_met = mass_reaction.boundary_metabolite
        if bc_met is not None \
           and bc_met in mass_reaction.model.boundary_conditions:
            # Boundary value is constant
            bc_value = mass_reaction.model.boundary_conditions.get(bc_met)
            if isinstance(bc_value, (integer_types, float)):
                constant = True
            # Boundary value is a function
            else:
                constant = False
            # Write boundary metabolite into reaction
            reaction = _write_specie_ref_into_reaction(
                reaction=reaction, metabolite=bc_met, stoich=stoichiometry,
                comparison_op=operator.gt, constant=constant,
                f_replace=f_replace)

    return reaction


def _write_reaction_kinetic_law_to_sbml(reaction, mass_reaction, f_replace):
    """Write MassReaction kinetic law information into SBMLDocument.

    Warnings
    --------
    This method is intended for internal use only.

    """
    rid = reaction.getId()
    rate_equation = strip_time(mass_reaction.rate)
    # If ID replacements were performed earlier then apply the ID
    # replacements for metabolite and parameter arguments in rate law also.
    met_id_subs = {}
    metabolite_ids_to_correct = []
    for arg in list(rate_equation.atoms(Symbol)):
        # Fix reaction ID in rate law arguments
        if mass_reaction.id in str(arg):
            met_id_subs.update({
                str(arg): str(arg).replace(mass_reaction.id, rid)})
        # Correct any boundary metabolite IDs in rate law
        elif str(arg) in mass_reaction.model.boundary_metabolites:
            metabolite_ids_to_correct.append(str(arg))
        # Correct any metabolite IDs in rate law
        else:
            # Ensure metabolite exists in model before replacing its ID
            try:
                met = mass_reaction.model.metabolites.get_by_id(str(arg))
            except KeyError:
                pass
            else:
                # Account for metabolites in the rate law that
                # are not considered reactants or products.
                if met not in mass_reaction.metabolites:
                    # Get specie id
                    msid = str(met)
                    if f_replace and F_SPECIE_REV in f_replace:
                        msid = f_replace[F_SPECIE_REV](msid)
                    # Create modifier specie reference and set the specie
                    msref = reaction.createModifier()
                    _check(msref,
                           "create modifier specie " + msid + _for_id(rid))
                    _check(msref.setSpecies(msid),
                           "set modifier specie species" + _for_id(msid))
                metabolite_ids_to_correct.append(str(met))

    # Match metabplite IDs in the rate equation to those in SBMLDocument
    met_id_subs.update(
        dict((sid, f_replace[F_SPECIE_REV](sid))
             if f_replace and F_SPECIE_REV in f_replace
             else (sid, sid) for sid in metabolite_ids_to_correct))

    # Make xml string of rate equation via sympy conversion to Math ML
    math_xml_str = _create_math_xml_str_from_sympy_expr(
        rate_equation.subs(met_id_subs))
    # Create kinetic law as AST math and write into the SBMLDocument
    kinetic_law = reaction.createKineticLaw()
    _check(kinetic_law, "create kinetic law" + _for_id(rid))
    _check(kinetic_law.setMath(libsbml.readMathMLFromString(math_xml_str)),
           "set math on kinetic law" + _for_id(rid))


def _write_reaction_gpr_to_sbml(reaction_fbc, mass_reaction, f_replace):
    """Write MassReaction gene reaction rule information into SBMLDocument.

    Warnings
    --------
    This method is intended for internal use only.

    """
    # Set the reaction GPR if it exists
    gpr = mass_reaction.gene_reaction_rule
    if gpr:
        rid = mass_reaction.id
        if f_replace and F_REACTION_REV in f_replace:
            rid = f_replace[F_REACTION_REV](rid)
        # Replace IDs in string
        if f_replace and F_GENE_REV in f_replace:
            gpr = gpr.replace("(", "( ")
            gpr = gpr.replace(")", " )")
            tokens = gpr.split(" ")
            for i, token in enumerate(tokens):
                if token not in [" ", "and", "or", "(", ")"]:
                    tokens[i] = f_replace[F_GENE_REV](token)
            gpr = " ".join(tokens)
        # Create assoication
        gpa = reaction_fbc.createGeneProductAssociation()
        _check(gpa, "create gene product association" + _for_id(rid))
        _check(gpa.setAssociation(gpr),
               "set gene product association" + _for_id(rid))


def _write_enzyme_modules_to_sbml(model_group, mass_obj, f_replace):
    """Write EnzymeModule and EnzymeModuleDict information into SBMLDocument.

    Utilizes the groups package extension.

    Warnings
    --------
    This method is intended for internal use only.

    """
    # Create group to represent EnzymeModule/EnzymeModuleDict information
    mass_notes = {}
    # Create groups for attributes, add them to enzyme group members list
    enzyme_group_members = []
    for attr_name, default in iteritems(_ORDERED_ENZYMEMODULE_DICT_DEFAULTS):
        attr_value = getattr(mass_obj, attr_name, None)
        # Ignore attributes set at default and the id, name, and S matrix.
        if attr_name in["id", "name", "S"] or attr_value in [default, None]:
            continue
        elif isinstance(attr_value, list):
            # Write attribute with list value as a group.
            enzyme_group_members += [
                _create_group(
                    model_group, gid=attr_name, members=attr_value,
                    f_replace=f_replace)]
        elif isinstance(attr_value, dict):
            # Write attribute with dict containing lists as a group of groups.
            members = []
            # Create enzyme group members (as groups) and add to list
            for i, (category, values) in enumerate(iteritems(attr_value)):
                members += [
                    _create_group(
                        model_group, gid="Category" + str(i), gname=category,
                        members=values, f_replace=f_replace)]
            enzyme_group_members += [
                _create_group(
                    model_group, gid=attr_name, members=members,
                    f_replace=f_replace)]
        elif isinstance(attr_value, Basic):
            # Write attribute with sympy equation to the group notes.
            mass_notes.update({
                attr_name: _create_math_xml_str_from_sympy_expr(attr_value.rhs)
            })
        else:
            # Add attribute to the group notes.
            mass_notes.update({attr_name: attr_value})

    # Create enzyme group
    enzyme_group = _create_group(
        model_group, gid=mass_obj.id, gname=mass_obj.name, kind="collection",
        members=enzyme_group_members, f_replace=f_replace)
    _sbase_notes_dict(enzyme_group, mass_notes)


def _write_enzyme_attr_info_to_notes(sbml_obj, mass_obj, f_replace=None):
    """Write the enzyme object attribute into SBMLDocument in the notes.

    Applies to the EnzymeModuleForms and EnzymeModuleReaction

    Warnings
    --------
    This method is intended for internal use only.

    """
    # Get original notes from mass object.
    mass_notes = mass_obj.notes.copy()
    # Get the subclass specific attributes to write into the notes.
    attributes = get_subclass_specific_attributes(mass_obj)
    default_values = [{} if "bound" in attr_name else ""
                      for attr_name in attributes]
    for attr_name, default_value in zip(attributes, default_values):
        attr_value = getattr(mass_obj, attr_name)
        # Skip attributes equal to their defaults
        if attr_value == default_value:
            continue
        elif "bound" in attr_name:
            # Create a string representation of the bound dict attr.
            bound_dict = {}
            for key, value in iteritems(attr_value):
                key = str(key)
                if f_replace and F_SPECIE_REV in f_replace:
                    key = f_replace[F_SPECIE_REV](key)
                bound_dict.update({key: value})
            # Add to notes to write
            mass_notes.update({
                attr_name: _make_bound_attr_str_repr(bound_dict)})
        else:
            # Add to notes to write
            mass_notes.update({attr_name: attr_value})
    # Write the notes to the SBML object.
    _sbase_notes_dict(sbml_obj, mass_notes)


def _create_group(model_group, gid, gname="", kind="collection", members=None,
                  f_replace=None):
    """Write an SBML group into SBMLDocument, and return the group.

    Additional group information is written into the notes attribute of the
    SBML object.

    Warnings
    --------
    This method is intended for internal use only.

    """
    group = model_group.createGroup()
    _check(group, "create group" + _for_id(gid))
    _check(group.setId(gid), "set group id" + _for_id(gid))
    _check(group.setName(gname), "set group name" + _for_id(gid))
    _check(group.setKind(kind), "set group kind" + _for_id(gid))

    if members:
        for mass_member in members:
            # Get member ID and name
            if not isinstance(mass_member, libsbml.SBase):
                mid = mass_member.id
                mname = mass_member.name
            else:
                mid = mass_member.getId()
                mname = ""
                if mass_member.isSetName():
                    mname = mass_member.getName()

            # ID replacements
            m_type = str(type(mass_member))
            if "Reaction" in m_type:
                if f_replace and F_REACTION_REV in f_replace:
                    mid = f_replace[F_REACTION_REV](mid)
            if "Metabolite" in m_type:
                if f_replace and F_SPECIE_REV in f_replace:
                    mid = f_replace[F_SPECIE_REV](mid)
            if "Gene" in m_type:
                if f_replace and F_GENE_REV in f_replace:
                    mid = f_replace[F_GENE_REV](mid)

            # Create member
            member = group.createMember()
            _check(member, "create group member " + mid + _for_id(gid))
            _check(member.setIdRef(mid),
                   "set group member id ref" + _for_id(mid))
            if mname:
                _check(member.setName(mname),
                       "set group member name" + _for_id(mid))

    return group


# -----------------------------------------------------------------------------
# Validation
# -----------------------------------------------------------------------------
def validate_sbml_model(filename):
    """TODO DOCSTRING."""
    print(filename)


def _check(value, message):
    """Check the libsbml return value and log error messages.

    If value is None, logs an error messaage constructed using 'message' and
    then exists with status code 1. If 'value' is an integer, it assumes it is
    an libSBML return status code. If the code value is
    LIBSBML_OPERATION_SUCCESS, returns without further action; if it is not,
    logs an error message constructed using 'message' along with text from
    libSBML explaining the meaning of the code, and exits with status code 1.
    """
    print(message)
    if value is None:
        LOGGER.error("Error: LibSBML returned a null value trying to "
                     "<%s>.", message)
    if isinstance(value, integer_types)\
       and value != libsbml.LIBSBML_OPERATION_SUCCESS:
        LOGGER.error("Error encountered trying to  <%s>.", message)
        LOGGER.error("LibSBML error code %s: %s", str(value),
                     libsbml.OperationReturnValue_toString(value).strip())
