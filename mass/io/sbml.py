# -*- coding: utf-8 -*-
"""TODO Module Docstrings.

CLASS IS EXPERIMENTAL. Features defined in class are necessary for
libRoadRunner compatibility, but not for full SBML. Support for full SBML will
be forthcoming.
"""
from __future__ import absolute_import

import datetime
import logging
import operator
import re

import libsbml

from six import integer_types, iteritems, string_types

from sympy import Symbol, mathml

from cobra.io.sbml import (_sbase_annotations, _sbase_notes_dict)

from mass.exceptions import MassSBMLError
from mass.util import strip_time

LOGGER = logging.getLogger(__name__)


# -----------------------------------------------------------------------------
# Defaults and constants for writing SBML
# -----------------------------------------------------------------------------
# For SBML level and core version, and package extension versions
SBML_LEVEL_VERSION = (3, 2)
FBC_VERSION = 2

DEFAULT_COMPARTMENT_STR = "default_compartment"  # Default compartment
BOUNDARY_COMPARTMENT_DICT = {"b": "boundary"}  # Boundary compartment

MASS_MOIETY_RE = re.compile("\[\S+\]")  # Precompiled regex for moieties
SBML_MOIETY_RE = re.compile("Moiety\S+$")

# For SBO terms
SBO_MODELING_FRAMEWORK = "SBO:0000062"

# For MathML representation of kinetic laws
MATHML_MML_TAG_RE = re.compile("\<mml:.*?\>|\<\/mml:.*?\>")
MATH_XML_FMT = '<math xmlns="http://www.w3.org/1998/Math/MathML">{0}</math>'
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


def _create_math_xml_from_expr(equation):
    """Remove 'mml' tags from XML representation of mathemtical equations.

    Sympy expressions can be converted into their equivalent MathML string
    through the `sympy.printing.mathml` module. However, the mathml module will
    interpret '_' as a subscript and '__' as a superscript. In order to prevent
    the module from misinterpreting underscores in identifiers, the underscores
    are converted into ampsersands. Additionally, all MathML presentation
    markup must be removed for similar reasons. After the XML string is made,
    the identifiers are returned to their original state.
    """
    underscore_replace = {str(arg): str(arg).replace("_", "&")
                          for arg in equation.atoms(Symbol)}
    math_xml_str = mathml(equation.subs(underscore_replace))
    math_xml_str = MATH_XML_FMT.format(math_xml_str.replace("&", "_"))
    math_xml_str = MATHML_MML_TAG_RE.sub("", math_xml_str)

    return math_xml_str


# -----------------------------------------------------------------------------
# Read SBML
# -----------------------------------------------------------------------------
def read_sbml_model(filename):
    """TODO DOCSTRING."""
    print(filename)


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


def _model_to_sbml(mass_model, f_replace=None, units=True):
    """Convert MassModel to SBMLDocument.

    Parameters
    ----------
    mass_model : mass.core.MassModel
        MassModel instance which is written to SBML
    f_replace : dict of replacement functions
        Replacement to apply on identifiers.
    units : boolean
        Should the units be written in the SBMLDocument. All units will be in
        terms of: {'Millimoles', 'Litre', 'Hours'}.

    Returns
    -------
    doc: libsbml.SBMLDocument

    """
    if f_replace is None:
        f_replace = {}

    doc, models = _create_SBMLDocument_and_model_objects(mass_model)
    model, model_fbc = models

    # Write model annotations into the SBMLDocument
    _sbase_annotations(model, mass_model.annotation)
    # Write the model Meta Information (ModelHistory) into the SBMLDocument
    if hasattr(mass_model, "_sbml"):
        meta = mass_model._sbml
        _write_model_meta_info_to_sbml(model, meta)
        # Set meta annotation and notes
        if "annotation" in meta:
            _sbase_annotations(doc, meta["annotation"])
        if "notes" in meta:
            _sbase_notes_dict(doc, meta["notes"])

    # TODO Write the units information into the SBMLDocument
    if units:
        _write_model_units_to_sbml(model, mass_model)

    # TODO Write the flux upper/lower bound parameters into the SBMLDocument
    # _write_flux_bounds_to_sbml()

    # TODO Write the rate law parameters into the SBMLDocument
    _write_model_parameters_to_sbml(model, mass_model, f_replace, units=units)

    # Write the compartment information into the SBMLDocument
    _write_model_compartments_to_sbml(model, mass_model)

    # Write the species information into the SBMLDocument
    # Includes MassMetabolites, EnzymeModuleForms, and their concentrations
    _write_model_species_to_sbml(model, mass_model, f_replace)
    # Write the boundary conditions into the SBMLDocument
    _write_model_boundary_conditions_to_sbml(model, mass_model, f_replace)
    # Write the gene information into the SBMLDocument
    _write_model_genes_to_sbml(model_fbc, mass_model, f_replace)
    # Write the reaction information into the SBMLDocument
    _write_model_reactions_to_sbml(model, mass_model, f_replace)

    # TODO Write the EnzymeModuleDict information to the SBMLDocument

    return doc


def _create_SBMLDocument_and_model_objects(mass_model):
    """Create and return the SBMLDocument and model objects.

    Will also set package extensions, model ID, meta ID, and name

    Warnings
    --------
    This method is intended for internal use only.

    """
    # Try creating the SBMLDocument object
    sbml_ns = libsbml.SBMLNamespaces(*SBML_LEVEL_VERSION)  # SBML L3V2 Core
    _check(sbml_ns.addPackageNamespace("fbc", FBC_VERSION),  # SBML L3 fbc-V2
           "adding fbc-V{0} package extension".format(FBC_VERSION))
    try:
        doc = libsbml.SBMLDocument(sbml_ns)
    except ValueError as e:
        MassSBMLError("Could not create SBMLDocument due to:\n" + str(e))

    # Set package extensions
    _check(doc.setPackageRequired("fbc", False),
           "set fbc package extension required to false")
    _check(doc.setSBOTerm(SBO_MODELING_FRAMEWORK),
           "set SBO term for modeling framework")

    # Create model
    model = doc.createModel()
    _check(model, "create model")

    # Set plugins
    model_fbc = model.getPlugin("fbc")
    _check(model_fbc, "get fbc plugin for model")
    _check(model_fbc.setStrict(False), "set fbc plugin strictness to true")

    # Set ID, meta ID, and name
    if mass_model.id is not None:
        _check(model.setId(mass_model.id), "set model id")
        _check(model.setMetaId("meta_" + mass_model.id), "set model meta id")
    else:
        _check(model.setMetaId("meta_model"), "set model meta id")
    if mass_model.name is not None:
        _check(model.setName(mass_model.name), "set model name")

    return doc, (model, model_fbc)


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


def _write_model_parameters_to_sbml(model, mass_model, f_replace, units=None):
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
                _create_parameter(model, pid=pid, value=value, sbo=None,
                                  constant=constant, units=units, udef=None)


def _create_parameter(model, pid, value, sbo=None, constant=True, units=None,
                      udef=None):
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
    # Set SBO term and units if desired.
    if sbo:
        _check(parameter.setSBOterm(sbo),
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
            "'%s' for entire model.", DEFAULT_COMPARTMENT_STR)
        compartment_dict = {DEFAULT_COMPARTMENT_STR: ""}
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


def _write_model_species_to_sbml(model, mass_model, f_replace):
    """Write MassModel species information into SBMLDocument.

    Warnings
    --------
    This method is intended for internal use only.

    """
    for metabolite in mass_model.metabolites:
        # Create specie and set ID, name, and compartment
        mid = metabolite.id
        if f_replace and F_SPECIE_REV in f_replace:
            mid = f_replace[F_SPECIE_REV](mid)
        specie = model.createSpecies()
        _check(specie, "create model specie" + _for_id(mid))
        _check(specie.setId(mid), "set specie id" + _for_id(mid))
        _check(specie.setName(metabolite.name),
               "set specie name" + _for_id(mid))
        # Use a generic compartment if no species have set compartments
        if not mass_model.compartments:
            cid = DEFAULT_COMPARTMENT_STR
        else:
            cid = metabolite.compartment
        _check(specie.setCompartment(cid),
               "set specie compartment" + _for_id(mid))
        # Set metabolite initial condition value
        if metabolite.initial_condition is not None:
            _check(
                specie.setInitialConcentration(metabolite.initial_condition),
                "set specie initial concentration" + _for_id(mid))

        # Set species constant, boundary condition, and unit restrictions
        _check(specie.setConstant(metabolite.fixed),
               "set specie constant" + _for_id(mid))
        _check(specie.setBoundaryCondition(False),
               "set specie boundary condition" + _for_id(mid))
        _check(specie.setHasOnlySubstanceUnits(False),
               "set specie unit substance restrictions" + _for_id(mid))

        # Set metabolite charge via fbc package
        specie_fbc = specie.getPlugin("fbc")
        _check(specie_fbc, "get fbc plugin for specie" + _for_id(mid))
        if metabolite.charge is not None:
            _check(specie_fbc.setCharge(metabolite.charge),
                   "set specie charge" + _for_id(mid))
        # Set metabolite formula via fbc package
        if metabolite.formula is not None:
            sbml_formula = _f_moiety_formula_rev(metabolite)
            _check(specie_fbc.setChemicalFormula(sbml_formula),
                   "set specie formula" + _for_id(mid))
        # Set metabolite annotation and notes
        _sbase_annotations(specie, metabolite.annotation)
        _sbase_notes_dict(specie, metabolite.notes)


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


def _write_model_reactions_to_sbml(model, mass_model, f_replace):
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

        # Set reaction annotation and notes
        _sbase_annotations(reaction, mass_reaction.annotation)
        _sbase_notes_dict(reaction, mass_reaction.notes)

        # Write specie references into reaction
        reaction = _write_reaction_specie_ref_to_sbml(
            reaction, mass_reaction=mass_reaction, f_replace=f_replace)
        # Get plugin for reaction
        reaction_fbc = reaction.getPlugin("fbc")
        _check(reaction_fbc, "get fbc plugin for reaction" + _for_id(rid))
        # TODO Flux bounds

        # Set the reaction GPR if it exists
        gpr = mass_reaction.gene_reaction_rule
        if gpr:
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

        # Set the reaction kinetic law
        _write_reaction_kinetic_law_to_sbml(
            reaction, mass_reaction=mass_reaction, f_replace=f_replace)


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
    math_xml_str = _create_math_xml_from_expr(rate_equation.subs(met_id_subs))
    # Create kinetic law as AST math and write into the SBMLDocument
    kinetic_law = reaction.createKineticLaw()
    _check(kinetic_law, "create kinetic law" + _for_id(rid))
    _check(kinetic_law.setMath(libsbml.readMathMLFromString(math_xml_str)),
           "set math on kinetic law" + _for_id(rid))


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


# -----------------------------------------------------------------------------
# Validation
# -----------------------------------------------------------------------------
def validate_sbml_model(filename):
    """TODO DOCSTRING."""
    print(filename)
