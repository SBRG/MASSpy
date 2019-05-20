# -*- coding: utf-8 -*-
"""TODO Module Docstrings.

CLASS IS EXPERIMENTAL. Features defined in class are necessary for
libRoadRunner compatibility, but not for full SBML. Support for full SBML will
be forthcoming.
"""
from __future__ import absolute_import

import datetime
import logging
import re

import libsbml

from six import integer_types, iteritems, string_types

from mass.exceptions import MassSBMLError

LOGGER = logging.getLogger(__name__)


# -----------------------------------------------------------------------------
# Defaults and constants for writing SBML
# -----------------------------------------------------------------------------
# For SBML level and core version, and package extension versions
SBML_LEVEL_VERSION = (3, 2)
FBC_VERSION = 2

DEFAULT_COMPARTMENT_STR = "default_compartment"  # Default compartment
MASS_MOIETY_RE = re.compile("\[\S+\]")  # Precompiled regex for moieties
SBML_MOIETY_RE = re.compile("Moiety\S+$")
# For SBO terms
SBO_MODELING_FRAMEWORK = "SBO:0000004"  # TODO Consider more specific term

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

    # TODO annotations for model
    # TODO notes for model

    # Write the model Meta Information (ModelHistory) to the SBMLDocument
    if hasattr(mass_model, "_sbml"):
        _write_model_meta_info_to_sbml(model, mass_model._sbml)
    # Write the model unit definitions to the SBMLDocument
    if units:
        _write_model_units_to_sbml()

    # TODO flux bounds for model

    # Write the compartment information
    _write_model_compartments_to_sbml(model, mass_model)

    # TODO enzyme module information

    # Write the species information
    _write_model_species_to_sbml(model, mass_model, f_replace)

    # Write the gene information
    _write_model_genes_to_sbml(model_fbc, mass_model, f_replace)

    # Write the reaction information
    _write_model_reactions_to_sbml(model, mass_model, f_replace)

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
    _check(model_fbc.setStrict(True), "set fbc plugin strictness to true")

    # Set  ID, meta ID, and name
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
    # TODO annotations for meta
    # TODO annotations for meta

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


def _write_model_units_to_sbml():
    """Write MassModel unit information into SBMLDocument.

    Currently, all units will be in terms of Millimole, Litres, and Hours.

    Warnings
    --------
    This method is intended for internal use only.

    """
    # TODO Enable other units and/or use of units from MassModel
    print("TODO: Units need to be implemented")


def _write_model_compartments_to_sbml(model, mass_model):
    """Write MassModel compartment information into SBMLDocument.

    Warnings
    --------
    This method is intended for internal use only.

    """
    # TODO Revisit based on new cobra compartment implementation
    if not mass_model.compartments:
        LOGGER.warning("No compartments found in MassModel. Therefore creating"
                       " compartment 'default_compartment' for entire model")
        compartment_dict = {DEFAULT_COMPARTMENT_STR: ""}
    else:
        compartment_dict = mass_model.compartments

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

        # TODO Set initial/fixed concentration, Clean up after IC/FC changes
        # Four possible cases
        # Metabolite varies, metabolite is on boundary
        # Metabolite varies, metabolite is not on boundary
        # Metabolite constant, metabolite is on boundary
        # Metabolite constant, metabolite is not on boundary
        boundary_condition = False
        constant = False

        # Set species constant, boundary condition, and unit restrictions
        _check(specie.setConstant(constant),
               "set specie constant" + _for_id(mid))
        _check(specie.setBoundaryCondition(boundary_condition),
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
        # TODO annotations for specie
        # TODO notes for specie


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

        # TODO annotations for gene
        # TODO notes for gene


def _write_model_reactions_to_sbml(model, mass_model, f_replace):
    """Write MassModel reaction information into SBMLDocument.

    Warnings
    --------
    This method is intended for internal use only.

    """
    # TODO Break into smaller functions
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
        # TODO annotations for reaction
        # TODO notes for reaction

        # Set reaction species references and stoichiometry
        for metabolite, stoichiometry in iteritems(mass_reaction.metabolites):
            # Get specie id
            sid = metabolite.id
            if f_replace and F_SPECIE_REV in f_replace:
                sid = f_replace[F_SPECIE_REV](sid)
            # Create specie reference as reactant or product
            # based on the reaction stoichiometry
            if stoichiometry < 0:
                sref = reaction.createReactant()
            else:
                sref = reaction.createProduct()
            _check(sref, "create specie reference " + sid + _for_id(rid))
            # Set ID, stoichiometry, and whether species is constant
            _check(sref.setSpecies(sid),
                   "set specie reference id" + _for_id(sid))
            _check(sref.setStoichiometry(abs(stoichiometry)),
                   "set specie reference stoichiometry" + _for_id(sid))
            _check(sref.setConstant(True),
                   "set specie reference constant" + _for_id(sid))

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
