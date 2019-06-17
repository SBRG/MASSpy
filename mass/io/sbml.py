# -*- coding: utf-8 -*-
"""TODO Module Docstrings."""
from __future__ import absolute_import

import datetime
import logging
import operator
import re
import traceback
import warnings
from collections import defaultdict
from io import StringIO

import libsbml

from six import integer_types, iteritems, itervalues, string_types

from sympy import Basic, Symbol, SympifyError, mathml, sympify

from cobra.core.dictlist import DictList
from cobra.core.gene import Gene
from cobra.io.sbml import (
    BOUND_MINUS_INF, BOUND_PLUS_INF, CobraSBMLError, LOWER_BOUND_ID,
    SBO_DEFAULT_FLUX_BOUND, SBO_FLUX_BOUND, UPPER_BOUND_ID, ZERO_BOUND_ID,
    _create_bound, _error_string,
    _get_doc_from_filename as _cobra_get_doc_from_filename, _parse_annotations,
    _sbase_annotations as _cobra_sbase_annotations, _sbase_notes_dict)

from mass.core.massconfiguration import MassConfiguration
from mass.core.massmetabolite import MassMetabolite
from mass.core.massmodel import MassModel
from mass.core.massreaction import MassReaction
from mass.core.units import (Unit, UnitDefinition, _SBML_BASE_UNIT_KINDS_DICT)
from mass.enzyme_modules.enzyme_module import EnzymeModule
from mass.enzyme_modules.enzyme_module_dict import (
    EnzymeModuleDict, _ORDERED_ENZYMEMODULE_DICT_DEFAULTS)
from mass.enzyme_modules.enzyme_module_form import (
    EnzymeModuleForm, _make_bound_attr_str_repr)
from mass.enzyme_modules.enzyme_module_reaction import EnzymeModuleReaction
from mass.exceptions import MassSBMLError
from mass.util.expressions import strip_time
from mass.util.util import get_subclass_specific_attributes

LOGGER = logging.getLogger(__name__)
# -----------------------------------------------------------------------------
# Defaults and constants for writing SBML
# -----------------------------------------------------------------------------
MASSCONFIGURATION = MassConfiguration()
# For SBML level and core version, and package extension versions
SBML_LEVEL_VERSION = (3, 1)
FBC_VERSION = 2
GROUPS_VERSION = 1

# Default compartment and boundary compartment definitions
DEFAULT_COMPARTMENT_DICT = MASSCONFIGURATION.default_compartment
BOUNDARY_COMPARTMENT_DICT = MASSCONFIGURATION.boundary_compartment

# Precompiled regex
CHAR_RE = re.compile(r"(?P<char>_Char\d*_)")           # For ASCII Char removal
MASS_MOIETY_RE = re.compile(r"\[\S+\]")                # For mass moieties
SBML_MOIETY_RE = re.compile(r"Moiety")                 # For SBML moieties
CATEGORY_GROUP_RE = re.compile(r"_Category\d$")        # For category groups
RATE_CONSTANT_RE = re.compile(r"k_(\S*)_(fwd|rev)\Z")  # For kf & kr
NOTES_RE = re.compile(r"\s*(\w+\s*\w*)\s*:\s*(...+)", re.DOTALL)  # For notes
# For parsing rate laws
KLAW_POW_RE = re.compile(r"pow\((?P<arg>...*?(?=,)), (?P<exp>\d*)\)")

# For SBO terms
SBO_MODELING_FRAMEWORK = "SBO:0000062"

# COBRA Flux units
COBRA_FLUX_UNIT = UnitDefinition(
    id="mmol_per_gDW_per_hr", name="cobra flux unit", list_of_units=[
        Unit(kind="mole", exponent=1, scale=-3, multiplier=1),
        Unit(kind="gram", exponent=-1, scale=0, multiplier=1),
        Unit(kind="second", exponent=-1, scale=0, multiplier=3600),
    ]
)
# -----------------------------------------------------------------------------
# Functions for replacements (import/export)
# -----------------------------------------------------------------------------
SBML_DOT = "__SBML_DOT__"
NUMBER_ID_PREFIX = "_"  # To prefix identifiers with starting with numbers


def _prefix_number_id(sid, obj_type=None):
    """Add a prefix to the beginning of the string if it starts with a number.

    Warnings
    --------
    This method is intended for internal use only.

    """
    needs_prefix = re.match(r"^\d", sid)
    if needs_prefix:
        if obj_type is not None:
            LOGGER.warning(
                "ID cannot begin with Unicode decimal digit, adding prefix "
                "'%s' to %s ID '%s'.", NUMBER_ID_PREFIX, obj_type, sid)
        sid = NUMBER_ID_PREFIX + sid

    return sid


def _clip(sid, prefix):
    """Clips a prefix from the beginning of a string if it exists.

    Warnings
    --------
    This method is intended for internal use only.

    """
    return sid[len(prefix):] if sid.startswith(prefix) else sid


def _f_gene(sid, prefix="G_"):
    """Clip gene prefix from id.

    Warnings
    --------
    This method is intended for internal use only.

    """
    sid = sid.replace(SBML_DOT, ".")
    return _clip(sid, prefix)


def _f_gene_rev(sid, prefix="G_"):
    """Add gene prefix to id.

    Warnings
    --------
    This method is intended for internal use only.

    """
    return prefix + sid.replace(".", SBML_DOT)


def _f_specie(sid, prefix="M_"):
    """Clip specie/metabolite prefix from id.

    Warnings
    --------
    This method is intended for internal use only.

    """
    return _clip(sid, prefix)


def _f_specie_rev(sid, prefix="M_"):
    """Add specie/metabolite prefix to id.

    Warnings
    --------
    This method is intended for internal use only.

    """
    return prefix + sid


def _f_reaction(sid, prefix="R_"):
    """Clip reaction prefix from id.

    Warnings
    --------
    This method is intended for internal use only.

    """
    return _clip(sid, prefix)


def _f_reaction_rev(sid, prefix="R_"):
    """Add reaction prefix to id.

    Warnings
    --------
    This method is intended for internal use only.

    """
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

ID_ASCII_REPLACE = {
    "_Char45_": "__",
    "_Char91_": "_",
    "_Char93_": "_",
}


def _f_moiety_formula(formula_str):
    """Fix moieties in formula from SBML compatible to be masspy compatible.

    Warnings
    --------
    This method is intended for internal use only.

    """
    mass_formula = None
    if formula_str is not None:
        # Split formula into moiety and regular formula
        formula_elements = SBML_MOIETY_RE.split(formula_str)
        mass_formula = formula_elements[0]
        if len(formula_elements) > 1:
            # Bracket moieties and add to beginning of formula
            formula_elements = [mass_formula] + [
                "[{0}]".format(e.upper()) for e in formula_elements[1:]]
            formula_elements.reverse()
            mass_formula = "-".join(formula_elements).rstrip("-")

    return mass_formula


def _f_moiety_formula_rev(metabolite):
    """Fix moieties in formula to be SBML compatible from masspy compatible.

    Warnings
    --------
    This method is intended for internal use only.

    """
    sbml_formula = metabolite.formula
    moieties = list(filter(MASS_MOIETY_RE.match, metabolite.elements))
    if moieties:
        for moiety in moieties:
            replacement = SBML_MOIETY_RE.pattern
            replacement += moiety.lstrip("[").rstrip("]").lower()
            if not re.match('^[A-Z]+[a-z]+$', replacement):
                LOGGER.warning(
                    "SBML user defined compounds must be in the form of a "
                    "single capital letter followed by zero or more lowercase "
                    "letters. Therefore removing all non-lowercase letters "
                    "from the Moiety.")
                replacement = re.sub("[0-9]", "", replacement)
            # Add to end of formula prevents issues if moiety ends with number
            sbml_formula = sbml_formula.replace(moiety, "")
            sbml_formula += replacement
        # Replace any dashes
        sbml_formula = sbml_formula.replace("-", "")

    return sbml_formula


def _fix_parameter_id(pid, rid=None, **kwargs):
    """Fix the rate constant ID, making it masspy compatible.

    Also removes ASCII characters if desired.

    """
    promote_local_parameters = kwargs["promote_local_parameters"]
    remove_char = kwargs["remove_char"]

    if promote_local_parameters and rid:
        if isinstance(rid, list):
            rid = rid[-1]
        pid = pid[len(rid) + 1:]

    rate_constant_fix = RATE_CONSTANT_RE.search(pid)
    if rate_constant_fix:
        rid, p_type = rate_constant_fix.groups()
        p_type = {"fwd": "kf_", "rev": "kr_"}[p_type]
        pid = pid[:rate_constant_fix.start()] + p_type + rid

    if remove_char:
        pid = _remove_char_from_id(pid)

    return pid


def _remove_char_from_id(sid):
    """Remove ASCII characters from an identifier."""
    for match in CHAR_RE.finditer(sid):
        char_key = match.group("char")
        if char_key in ID_ASCII_REPLACE:
            replacement_str = ID_ASCII_REPLACE[char_key]
        else:
            replacement_str = "_"
        sid = re.sub(char_key, replacement_str, sid)
    return sid


# -----------------------------------------------------------------------------
# MathML
# -----------------------------------------------------------------------------
# For MathML representation of kinetic laws and other sympy equations
MATHML_MML_TAG_RE = re.compile("\<mml:.*?\>|\<\/mml:.*?\>")
MATH_XML_FMT = '<math xmlns="http://www.w3.org/1998/Math/MathML">{0}</math>'
SBML_MATH_PACKAGE_NAME = "l{0}v{1}extendedmath".format(*SBML_LEVEL_VERSION)


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

    Warnings
    --------
    This method is intended for internal use only.

    """
    underscore_replace = {str(arg): str(arg).replace("_", "&")
                          for arg in sympy_equation.atoms(Symbol)}
    math_xml_str = mathml(sympy_equation.subs(underscore_replace))
    math_xml_str = MATH_XML_FMT.format(math_xml_str.replace("&", "_"))
    math_xml_str = MATHML_MML_TAG_RE.sub("", math_xml_str)

    return math_xml_str


# -----------------------------------------------------------------------------
# Read SBML
# -----------------------------------------------------------------------------
def read_sbml_model(filename, f_replace=None, **kwargs):
    """Read SBML model from the given filename into a MassModel.

    If the given filename ends with the suffix ''.gz'' (for example,
    ''myfile.xml.gz'),' the file is assumed to be compressed in gzip
    format and will be automatically decompressed upon reading. Similarly,
    if the given filename ends with ''.zip'' or ''.bz2',' the file is
    assumed to be compressed in zip or bzip2 format (respectively).  Files
    whose names lack these suffixes will be read uncompressed.  Note that
    if the file is in zip format but the archive contains more than one
    file, only the first file in the archive will be read and the rest ignored.
    To read a gzip/zip file, libSBML needs to be configured and linked
    with the zlib library at compile time.  It also needs to be linked
    with the bzip2 library to read files in bzip2 format. (Both of these
    are the default configurations for libSBML.)
    This function supports SBML with FBC-v1 and FBC-v2. FBC-v1 models
    are converted to FBC-v2 models before reading.

    Parameters
    ----------
    filename: path to SBML file, or SBML string, or SBML file handle
        SBML which is read into a MassModel.
    f_replace: dict of replacement functions for id replacement
        Dictionary of replacement functions for gene, specie, and reaction.
        By default the following id changes are performed on import:
        clip G_ from genes, clip M_ from species, clip R_ from reactions
        If no replacements should be performed, set f_replace={}.
    number: data type of stoichiometry: {float, int}
        In which data type should the stoichiometry be parsed.
        Default is float.
    set_missing_bounds: bool
        If True, missing bounds are set to default bounds from the
        mass.MassConfiguration.
        Default is False.
    remove_char: bool
        If True, remove ASCII characters from identifiers.
        Default is True.
    promote_local_parameters: bool
        If True, promote all local parameters in kinetic laws to global
        parameters before parsing model.
        Default is False.

    Returns
    -------
    mass.core.MassModel

    Notes
    -----
    Provided file handles cannot be opened in binary mode, i.e., use
        with open(path, "r" as f):
            read_sbml_model(f)
    File handles to compressed files are not supported yet.

    """
    # Try getting the SBMLDocument from the file 'filename' and if successful,
    # try creating a model from the SBMLDocument.
    try:
        doc = _get_doc_from_filename(filename)
        return _sbml_to_model(doc, f_replace=f_replace, **kwargs)

    except IOError as e:
        # Raise Import/export error if error not SBML parsing related.
        raise e

    except Exception:
        # Log traceback and raise a MassSBMLError for parsing related errors.
        LOGGER.error(traceback.print_exc())
        raise MassSBMLError(
            "Something went wrong reading the SBML model. Most likely the SBML"
            " model is not valid. Please check that your model is valid using "
            "the `cobra.io.sbml.validate_sbml_model` function or via the "
            "online validator at http://sbml.org/validator .\n"
            "\t`(model, errors) = validate_sbml_model(filename)`"
            "\nIf the model is valid and cannot be read please open an issue "
            "at https://github.com/SBRG/masspy/issues .")


def _get_doc_from_filename(filename):
    """Get SBMLDocument from given filename.

    Parameters
    ----------
    filename : path to SBML, or SBML string, or filehandle

    Returns
    -------
    libsbml.SBMLDocument

    """
    try:
        doc = _cobra_get_doc_from_filename(filename)
    # Reclassify cobra error as a mass error
    except CobraSBMLError as e:
        raise MassSBMLError(e)

    return doc


def _sbml_to_model(doc, f_replace=None, **kwargs):
    """Create MassModel from SBMLDocument.

    Parameters
    ----------
    doc: libsbml.SBMLDocument
        SBMLDocument which is read into a MassModel
    f_replace: dict of replacement functions for id replacement
        Dictionary of replacement functions for gene, specie, and reaction.
        By default the following id changes are performed on import:
        clip G_ from genes, clip M_ from species, clip R_ from reactions
        If no replacements should be performed, set f_replace={}.
    number: data type of stoichiometry: {float, int}
        In which data type should the stoichiometry be parsed.
        Default is float.
    set_missing_bounds: bool
        If True, missing bounds are set to default bounds from the
        mass.MassConfiguration.
        Default is False.
    remove_char: bool
        If True, remove ASCII characters from identifiers.
        Default is True.
    promote_local_parameters: bool
        If True, promote all local parameters in kinetic laws to global
        parameters before parsing model.
        Default is False.

    Returns
    -------
    mass.core.MassModel

    Warnings
    --------
    This method is intended for internal use only.

    """
    default_kwargs = {
        "number": float,
        "set_missing_bounds": False,
        "remove_char": False,
        "promote_local_parameters": False}
    # Use default kwargs or get missing kwargs with their default values
    if kwargs is None:
        kwargs = default_kwargs
    else:
        for key, value in iteritems(default_kwargs):
            if key not in kwargs:
                kwargs[key] = value
    if f_replace is None:
        f_replace = F_REPLACE

    doc, models = _get_sbml_models_from_doc(doc, **kwargs)
    model, model_fbc, model_groups = models

    # Get model ID
    model_id = _check_required(model, model.getIdAttribute(), "id")

    # Create MassModel
    mass_model = MassModel(model_id)
    mass_model.name = model.getName()

    # Read the MassModel meta information from the SBMLDocument
    mass_model._sbml = _read_model_meta_info_from_sbml(doc, model, model_id)

    # Read the MassModel notes and annotations from the SBMLDocument
    mass_model.notes = _parse_notes_dict(model)
    if "description" in mass_model.notes:
        # Set description and remove from the dictionary
        mass_model.description = mass_model.notes["description"]
        del mass_model.notes["description"]
    mass_model.annotation = _parse_annotations(model)

    # Read the MassModel unit information from the SBMLDocument
    mass_model.add_units(_read_model_units_from_sbml(model))

    # Read the MassModel compartment information from the SBMLDocument
    # TODO Revisit based on new cobra compartment implementation
    mass_model.compartments = _read_model_compartments_from_sbml(model)

    # Read the MassModel species information from the SBMLDocument
    # Includes MassMetabolites and EnzymeModuleForms, boundary conditions are
    # stored in order to be added after boundary reactions.
    metabolites, bc_values = _read_model_species_from_sbml(model, f_replace,
                                                           **kwargs)
    mass_model.add_metabolites(metabolites)
    if model_fbc:
        genes = _read_model_genes_from_sbml(model_fbc, f_replace=f_replace)
        mass_model.genes += genes
    else:
        # TODO read legacy genes here
        pass

    # Read the MassModel reaction information from the SBMLDocument
    # Includes MassReactions, EnzymeModuleReactions, kinetic laws,
    # and GPR information.
    reactions, custom_rates, local_params = _read_model_reactions_from_sbml(
        model, f_replace=f_replace, metabolites=DictList(metabolites),
        **kwargs)
    mass_model.add_reactions(reactions)

    # Read the MassModel parameter information from the SBMLDocument
    parameters, bounds = _read_model_parameters_from_sbml(model,
                                                          f_replace=f_replace,
                                                          **kwargs)
    # TODO assign reaction bounds

    # Add the parameters and boundary conditions to the model
    parameters.update(local_params)
    parameters.update(bc_values)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        mass_model.update_parameters(parameters)

    # Add custom rates
    mass_model.update_custom_rates(custom_rates)

    # TODO Revisit when cobra groups implemented, store in mass_groups
    if model_groups:
        enzyme_groups, mass_groups = _read_model_groups_from_sbml(model_groups)
        # Read the enzyme groups into EnzymeModuleDicts
        enzyme_module_dicts_list = _read_enzyme_modules_from_sbml(
            enzyme_groups, f_replace=f_replace,
            remove_char=kwargs["remove_char"])

        # Add enzyme modules to model
        for enzyme_module_dict in enzyme_module_dicts_list:
            enzyme_module_dict._update_object_pointers(mass_model)
            if enzyme_module_dict["model"] == "EnzymeModule":
                mass_model = EnzymeModule(mass_model)
                for attr, value in iteritems(enzyme_module_dict):
                    if attr in ["id", "name", "description", "S", "model"]:
                        continue
                    setattr(mass_model, attr, value)
            else:
                # The dict is an EnzymeModuleDict to add to the model.
                mass_model.enzyme_modules.append(enzyme_module_dict)
                enzyme_module_dict["model"] = mass_model
                enzyme_module_dict._make_enzyme_stoichiometric_matrix(True)

    return mass_model


def _get_sbml_models_from_doc(doc, **kwargs):
    """Return the SBML model objects, performing version conversions if needed.

    Warnings
    --------
    This method is intended for internal use only.

    """
    model = doc.getModel()
    if model is None:
        raise MassSBMLError("No SBML model detected in file.")

    # Convert to current SBML level
    doc_level = doc.getLevel()
    doc_version = doc.getVersion()
    # Determine whether conversion to a different SBML level/version is needed
    if (doc_level, doc_version) != SBML_LEVEL_VERSION:
        # Create conversion message
        conversion_msg = str(
            "from Level {0} Version {1} to Level {2} Version {3}".format(
                doc_level, doc_version, *SBML_LEVEL_VERSION))
        LOGGER.warning(
            "SBML model is Level %s Version %s, attempting to convert model "
            "%s before parsing.", doc_level, doc_version, conversion_msg)
        # Return exception if conversion fails
        result = doc.setLevelAndVersion(*SBML_LEVEL_VERSION)
        if result != libsbml.LIBSBML_OPERATION_SUCCESS:
            raise Exception("Conversion " + conversion_msg + " failed")

    def _convert_legacy_package_extensions(doc, plugin_str, desired_vers):
        """Convert packages from legacy to current versions."""
        # Get plugin and plugin version
        doc_plugin = doc.getPlugin(plugin_str)
        plugin_vers = doc_plugin.getPackageVersion()
        # Perform conversion if plugin version is less than the desired version
        if plugin_vers < desired_vers:
            LOGGER.warning(
                "Loading SBML with %s-v%s (models should be encoded using "
                "%s-v%s).", plugin_str, plugin_vers, plugin_str, desired_vers)
            # Get SBML conversion properties and add conversion option
            conversion_properties = libsbml.ConversionProperties()
            conversion_properties.addOption(
                "convert {0} v{1} to {0} v{2}".format(
                    plugin_str, plugin_vers, desired_vers), True,
                "Convert {0}-v{1} to {0}-v{2}".format(
                    plugin_str.upper(), plugin_vers, desired_vers))
            # Convert the doc and check if successful
            result = doc.convert(conversion_properties)
            if result != libsbml.LIBSBML_OPERATION_SUCCESS:
                raise Exception(
                    "Conversion of SBML {0} v{1} to {0} v{2} failed".format(
                        plugin_str, plugin_vers, desired_vers))
        return doc

    model_fbc = model.getPlugin("fbc")
    if not model_fbc:
        LOGGER.warning("Model does not contain SBML fbc package information.")
    else:
        # Convert fbc-v1 legacy models
        doc = _convert_legacy_package_extensions(doc, "fbc", FBC_VERSION)

    model_groups = model.getPlugin("groups")

    # Promote local parameters to global parameters
    if kwargs["promote_local_parameters"]:
        conversion_properties = libsbml.ConversionProperties()
        conversion_properties.addOption(
            "promoteLocalParameters", True,
            "Promotes all Local Parameters to Global ones")
        # Convert the doc and check if successful
        result = doc.convert(conversion_properties)
        if result != libsbml.LIBSBML_OPERATION_SUCCESS:
            raise Exception(
                "Conversion of local parameters to global parameters failed")

    return doc, (model, model_fbc, model_groups)


def _read_model_meta_info_from_sbml(doc, model, model_id=None):
    """Read SBML meta information to create a dict of meta information.

    Warnings
    --------
    This method is intended for internal use only.

    """
    # Make meta information as a dict
    meta_dict = {
        "model.id": model_id,
        "level": model.getLevel(),
        "version": model.getVersion()}
    # Read Model history, creation date, and creators
    created = None
    creators = []
    modified = []
    if model.isSetModelHistory():
        history = model.getModelHistory()
        # Get creation date from SBML, if one is set in the SBML.
        if history.isSetCreatedDate():
            date = history.getCreatedDate()
            created = date.getDateAsString()
        if history.isSetModifiedDate():
            for k in range(history.getNumModifiedDates()):
                date = history.getModifiedDate(k)
                modified += [date.getDateAsString()]
        # Get the list of creators, if any
        list_of_creators = history.getListCreators()
        # Make a dict of creator attributes and values
        for creator in list_of_creators:
            creator_dict = {}
            # Make a dict for the creator attributes initialized to None
            for key in ["familyName", "givenName", "organization", "email"]:
                attr = key[0].capitalize() + key[1:]
                if getattr(creator, "isSet" + attr)():
                    creator_dict[key] = getattr(creator, "get" + attr)()
                else:
                    creator_dict[key] = None
            creators.append(creator_dict)

    # Get version and package info
    info = "<{0}> SBML L{1}V{2}".format(
        model_id, model.getLevel(), model.getVersion())
    packages = {}
    for k in range(doc.getNumPlugins()):
        plugin = doc.getPlugin(k)
        key, value = (plugin.getPackageName(), plugin.getPackageVersion())
        packages[key] = value
        info += ", {0}-v{1}".format(key, value)
        if key not in [SBML_MATH_PACKAGE_NAME, "fbc", "groups"]:
            LOGGER.warning(
                "SBML package '%s' not supported by masspy. Therefore the "
                "package information is not parsed.", key)

    # Add meta information to the dict
    meta_dict["creators"] = creators
    meta_dict["created"] = created
    meta_dict["modified"] = modified
    meta_dict["notes"] = _parse_notes_dict(doc)
    meta_dict["annotation"] = _parse_annotations(doc)
    meta_dict["info"] = info
    meta_dict["packages"] = packages

    return meta_dict


def _read_model_units_from_sbml(model):
    """Read the SBML units, returning them as a list.

    UnitDefinitions are created and returned in a list.

    Warnings
    --------
    This method is intended for internal use only.

    """
    list_of_mass_unit_defs = []
    for sbml_unit_def in model.getListOfUnitDefinitions():
        # Get the ID and name attributes
        try:
            unit_def_id = sbml_unit_def.getIdAttribute()
            unit_def_id = _check_required(sbml_unit_def, unit_def_id, "id")
            unit_def_name = ""
            if sbml_unit_def.isSetName():
                unit_def_name = sbml_unit_def.getName()

            # Get list of units that make the UnitDefinition
            list_of_sbml_units = sbml_unit_def.getListOfUnits()

            # Make a list for holding mass UnitDefinitions
            list_of_mass_units = []
            for sbml_unit in list_of_sbml_units:
                unit_args = {
                    "kind": None,
                    "exponent": None,
                    "scale": None,
                    "multiplier": None}
                for key in unit_args:
                    value = getattr(sbml_unit, "get" + key.capitalize())()
                    if key == "kind" and value == 36:
                        raise MassSBMLError(
                            "Cannot set Unit 'kind' attribute as 'invalid'")
                    value = _check_required(sbml_unit, value, key)
                    unit_args[key] = value
                if None in list(itervalues(unit_args)):
                    raise MassSBMLError(
                        "Missing Unit requirements: " + ", ".join((
                            key for key, value in iteritems(unit_args)
                            if value is None)))
                list_of_mass_units += [Unit(**unit_args)]
            # Create the mass UnitDefinition and add it to a list
            mass_unit_def = UnitDefinition(unit_def_id, name=unit_def_name,
                                           list_of_units=list_of_mass_units)
            list_of_mass_unit_defs.append(mass_unit_def)
        except MassSBMLError as e:
            # Failed Units and UnitDefinitions should not prevent model
            # from loading, log as warnings instead
            LOGGER.warning(
                "Skipping UnitDefinition '%s'. %s", unit_def_id, str(e))
            continue

    return list_of_mass_unit_defs


def _read_model_compartments_from_sbml(model):
    """Read the SBML compartments, returning them as a dict.

    Compartments are created and returned in a dict.

    Warnings
    --------
    This method is intended for internal use only.

    """
    # TODO Revisit based on new cobra compartment implementation
    compartments = {}
    for compartment in model.getListOfCompartments():
        cid = _check_required(compartment, compartment.getIdAttribute(), "id")
        # Skip the default compartment
        if cid in MASSCONFIGURATION.default_compartment:
            continue
        # Add compartment to the compartments to be returned
        compartments[cid] = compartment.getName()

    return compartments


def _read_model_species_from_sbml(model, f_replace=None, **kwargs):
    """Read SBML species and boundary conditions and return them.

    MassMetabolites and EnzymeModuleForms are created and returned in a list.
    Boundary conditions are returned as a dict.

    Warnings
    --------
    This method is intended for internal use only.

    """
    metabolite_list = []
    if model.getNumSpecies() == 0:
        LOGGER.warning("No metabolites in model.")

    boundary_values_dict = {}
    for specie in model.getListOfSpecies():
        # Check the identifier
        sid = _check_required(specie, specie.getIdAttribute(), "id")
        if f_replace and F_SPECIE in f_replace:
            sid = f_replace[F_SPECIE](sid)
        sid = _prefix_number_id(sid, "Metabolite")
        if kwargs["remove_char"]:
            sid = _remove_char_from_id(sid)
        # Handle boundary conditions
        if specie.getBoundaryCondition():
            # Collect boundary metabolites and their values to be added after
            # boundary reactions have been read into the model.
            if specie.getConstant():
                boundary_values_dict[sid] = specie.getInitialConcentration()
            else:
                # TODO handle functions of time for boundary conditions
                boundary_values_dict[sid] = specie.getInitialConcentration()
            continue

        # Set MassMetabolite id, name, and fixed attributes
        metabolite = MassMetabolite(sid)
        metabolite.name = specie.getName()
        metabolite.fixed = specie.getConstant()
        # Set the MassMetabolite initial condition, if any
        if specie.isSetInitialConcentration():
            metabolite.initial_condition = specie.getInitialConcentration()
        # Set the MassMetabolite compartment if it is not the generic default
        cid = specie.getCompartment()
        if cid not in MASSCONFIGURATION.default_compartment:
            metabolite.compartment = cid
        # Get the fbc specie and set the charge and formula
        specie_fbc = specie.getPlugin("fbc")
        if specie_fbc:
            metabolite.charge = specie_fbc.getCharge()
            sbml_formula = specie_fbc.getChemicalFormula()
            metabolite.formula = _f_moiety_formula(sbml_formula)
        else:
            # TODO Read legacy SBML species here
            pass

        # Using notes information, determine whether the metabolite should
        # be added to the model as an EnzymeModuleForm.
        notes = _parse_notes_dict(specie)
        if "enzyme_module_id" in notes:
            metabolite = EnzymeModuleForm(metabolite)
            _read_enzyme_attr_info_from_notes(
                metabolite, notes, f_replace=f_replace)

        # Add the notes and annotations to the metabolite
        metabolite.notes = notes
        metabolite.annotation = _parse_annotations(specie)
        # Append metabolite to list of metabolites
        metabolite_list.append(metabolite)

    return metabolite_list, boundary_values_dict


def _read_model_genes_from_sbml(model_fbc, f_replace=None):
    """Read the SBML genes, returning them as a list.

    cobra Genes are created and returned in a list.

    Warnings
    --------
    This method is intended for internal use only.

    """
    gene_list = []
    for gp in model_fbc.getListOfGeneProducts():
        # Check the identifier
        gid = _check_required(gp, gp.getIdAttribute(), "id")
        if f_replace and F_GENE in f_replace:
            gid = f_replace[F_GENE](gid)
        gid = _prefix_number_id(gid, "Gene")

        # Set the gene ID and name
        gene = Gene(gid)
        gene.name = gp.getName() if gp.isSetName() else gid
        # Parse notes and annotations
        gene.notes = _parse_notes_dict(gp)
        gene.annotation = _parse_annotations(gp)
        # Add gene to list of genes
        gene_list.append(gene)

    return gene_list


def _read_model_reactions_from_sbml(model, f_replace=None, metabolites=None,
                                    **kwargs):
    """Read the SBML reactions, returning them as a list.

    MassReactions and EnzymeModuleReactions are created and returned in a list.

    Warnings
    --------
    This method is intended for internal use only.

    """
    reaction_list = []
    if model.getNumReactions() == 0:
        LOGGER.warning("No reactions in model.")

    custom_rates_dict = {}
    local_parameters_dict = {}
    for reaction in model.getListOfReactions():
        # Check the identifier
        rid = _check_required(reaction, reaction.getIdAttribute(), "id")
        if f_replace and F_REACTION in f_replace:
            rid = f_replace[F_REACTION](rid)
        rid = _prefix_number_id(rid, "Reaction")
        if kwargs["remove_char"]:
            rid = _remove_char_from_id(rid)
        # Set the reaction ID, name, and reversible
        mass_reaction = MassReaction(rid)
        mass_reaction.name = reaction.getName()
        mass_reaction.reversible = reaction.getReversible()

        # Determine reaction stoichiometry
        met_stoic = _read_reaction_stoichiometry_from_sbml(
            reaction, f_replace=f_replace, **kwargs)
        # Get metabolite objects and set reaction stoichiometry
        stoichiometry = {metabolites.get_by_id(mid): stoichiometry
                         for mid, stoichiometry in iteritems(met_stoic)
                         if mid in metabolites}
        mass_reaction.add_metabolites(stoichiometry)

        if reaction.isSetKineticLaw():
            rate_eq, local_parameters = _read_reaction_kinetic_law_from_sbml(
                reaction, mass_reaction, metabolites, f_replace=f_replace,
                **kwargs)
            # Get the rate equation, check if it is a mass action rate law.
            # If not a mass action rate law, assume it is a custom rate law
            mass_action_rates = [strip_time(mass_reaction.get_rate_law(x))
                                 for x in range(1, 4)]
            local_parameters_dict.update(local_parameters)
            if rate_eq in mass_action_rates:
                # Rate law is a mass action rate law identical to the one
                # automatically generated when reaction rate type is 1, 2, or 3
                for x, rate in enumerate(mass_action_rates):
                    if rate_eq == rate:
                        mass_reaction._rtype = int(x + 1)
            else:
                custom_rates_dict[mass_reaction] = rate_eq
        else:
            LOGGER.warning(
                "No kinetic law found for reaction '%s'. Therefore, assigning "
                "the reaction a rate law based on Mass Action Kinetics.",
                reaction.id)

        # TODO add GPR to the reaction
        # TODO add flux bounds to reaction, including missing
        # Using notes information, determine whether the reaction should
        # be added to the model as an EnzymeModuleReaction.
        notes = _parse_notes_dict(reaction)
        if "enzyme_module_id" in notes:
            mass_reaction = EnzymeModuleReaction(mass_reaction)
            _read_enzyme_attr_info_from_notes(
                mass_reaction, notes, f_replace=f_replace)

        # Add the notes and annotations to the reaction
        mass_reaction.notes = notes
        # mass_reaction.annotation = _parse_annotations(reaction)

        # Add reaction to list of reaction
        reaction_list.append(mass_reaction)

    return reaction_list, custom_rates_dict, local_parameters_dict


def _read_reaction_stoichiometry_from_sbml(reaction, f_replace=None, **kwargs):
    """Read the SBML reaction stoichiometry, returning it as a dict.

    Warnings
    --------
    This method is intended for internal use only.

    """
    number = kwargs["number"]
    remove_char = kwargs["remove_char"]
    stoichiometry = defaultdict(lambda: 0)
    for attr, sign in zip(["Reactants", "Products"], [-1, 1]):
        # Iterate through list
        for sref in getattr(reaction, "getListOf" + attr)():
            # Check specie ID
            sid = _check_required(sref, sref.getSpecies(), "species")
            if f_replace and F_SPECIE in f_replace:
                sid = f_replace[F_SPECIE](sid)
            sid = _prefix_number_id(sid)
            if remove_char:
                sid = _remove_char_from_id(sid)
            # Add specie IDs and stoichiometry to dict
            stoichiometry[sid] += sign * number(
                _check_required(sref, sref.getStoichiometry(),
                                "stoichiometry"))
    return stoichiometry


def _read_reaction_kinetic_law_from_sbml(reaction, mass_reaction, metabolites,
                                         f_replace=None, **kwargs):
    """Read the SBML reaction kinetic law and return it.

    Warnings
    --------
    This method is intended for internal use only.

    """
    # Get the reaction ID
    rid = reaction.getIdAttribute()
    # Get the kinetic law and the rate equation as a string.
    kinetic_law = reaction.getKineticLaw()
    rate_eq = _check_required(kinetic_law, kinetic_law.getFormula(),
                              "formula")
    # Perform substitution for certain mathematical operatiosn to sympify rate
    for match in KLAW_POW_RE.finditer(rate_eq):
        old = match.group(0)
        new = "(({0})**{1})".format(match.group("arg"), match.group("exp"))
        rate_eq = rate_eq.replace(old, new)
    # If ID replacements were performed then apply the ID replacements
    # for metabolite and parameter arguments in rate law also.
    id_subs = {}
    try:
        rate_eq = sympify(rate_eq)
    except SympifyError as e:
        raise MassSBMLError(e)
    # Get local parameters from the kinetic law
    local_parameters = {}
    if not kwargs["promote_local_parameters"]:
        for local_parameter in kinetic_law.getListOfLocalParameters():
            pid = _check_required(
                local_parameter, local_parameter.getIdAttribute(), "id")
            value = local_parameter.getValue()
            pid = _fix_parameter_id(pid, rid=rid, **kwargs)
            local_parameters[pid] = value

    for arg in list(rate_eq.atoms(Symbol)):
        new_arg = str(arg)
        # Fix reaction ID in rate law arguments
        if rid in new_arg:
            # Fix parameter identifier if needed.
            new_arg = _fix_parameter_id(new_arg, rid=rid, **kwargs)
            id_subs.update({arg: new_arg})
            continue

        # Correct any metabolite IDs in rate law
        if f_replace and F_SPECIE in f_replace:
            new_arg = f_replace[F_SPECIE](new_arg)
        new_arg = _prefix_number_id(new_arg)
        if kwargs["remove_char"]:
            new_arg = _remove_char_from_id(new_arg)
        # Ensure metabolite exists
        try:
            new_arg = metabolites.get_by_id(new_arg)
        except KeyError as e:
            if mass_reaction.boundary \
               and mass_reaction.boundary_metabolite == str(e).strip("'"):
                new_arg = mass_reaction.boundary_metabolite
            elif not mass_reaction.reactants or not mass_reaction.products:
                msg = "({0}.boundary: {1}, {0}.boundary_metabolite: {2}!={3})"
                msg = msg.format(
                    mass_reaction.id, mass_reaction.boundary, str(new_arg),
                    mass_reaction._make_boundary_metabolites())
                LOGGER.warning(
                    "Reaction '%s' appears to be on the boundary, but %s.",
                    mass_reaction, msg)
        else:
            new_arg = str(new_arg)
        finally:
            id_subs.update({arg: new_arg})
    # Make rate equation
    rate_eq = rate_eq.subs(id_subs)
    return rate_eq, local_parameters


def _read_model_parameters_from_sbml(model, **kwargs):
    """Read the SBML model parameters and return them in dicts.

    Warnings
    --------
    This method is intended for internal use only.

    """
    bound_pid_list = [
        LOWER_BOUND_ID, UPPER_BOUND_ID, ZERO_BOUND_ID, BOUND_MINUS_INF,
        BOUND_PLUS_INF, "_lower_bound", "_upper_bound"]

    parameters_dict = {}
    bounds_dict = {}
    reaction_ids = [
        reaction.getIdAttribute() for reaction in model.getListOfReactions()]
    for parameter in model.getListOfParameters():
        # Get parameter ID and value
        pid = _check_required(
            parameter, parameter.getIdAttribute(), "id")
        value = parameter.getValue()
        # Check if parameter is a flux bound
        if any([x in pid for x in bound_pid_list]):
            bounds_dict[pid] = value
            continue
        rid = list(filter(lambda x: pid.startswith(x), reaction_ids))
        pid = _fix_parameter_id(pid, rid=rid, **kwargs)
        parameters_dict.update({pid: value})

    return parameters_dict, bounds_dict


def _read_model_groups_from_sbml(model_groups):
    """Read the SBML groups, returning enzymes and mass groups seperately.

    Warnings
    --------
    This method is intended for internal use only.

    """
    # Read the MassModel enzyme modules information from the SBMLDocument
    # Also determine whether the MassModel should be an EnzymeModule
    other_groups = []
    enzyme_groups = defaultdict(list)
    for group in model_groups.getListOfGroups():
        # Get the group ID
        group_id = _check_required(group, group.getIdAttribute(), "id")
        # Determine whether the group represents an enzyme
        try:
            enzyme_module_id, attr = group_id.split("_", 1)
            if CATEGORY_GROUP_RE.search(attr):
                attr = attr[:CATEGORY_GROUP_RE.search(attr).start()]
            if attr not in list(_ORDERED_ENZYMEMODULE_DICT_DEFAULTS)\
               and "EnzymeModule" not in attr:
                raise ValueError
        except ValueError:
            # TODO revisit for other groups
            other_groups.append(group)
            LOGGER.warning("Group '%s' not parsed.", group_id)
        else:
            # Add enzyme group to the dictionary
            enzyme_groups[enzyme_module_id].append(group)

    return enzyme_groups, other_groups


def _read_enzyme_modules_from_sbml(enzyme_group_dict, f_replace=None,
                                   remove_char=None):
    """Read the SBML groups for enzymes, returning them as a list.

    EnzymeModuleDicts are created and returned in a list.

    Warnings
    --------
    This method is intended for internal use only.

    """
    def _get_group_members(group, gid=None, remove_char=None):
        """Get the group members of a group and performn ID corrections."""
        member_list = []
        for member in group.getListOfMembers():
            if member.isSetIdRef():
                obj_id = member.getIdRef()
            elif member.isSetMetaIdRef():
                obj_id = member.getMetaIdRef()
            # Correct object ID if necessary.
            replace_key = {
                "enzyme_module_ligands": F_SPECIE,
                "enzyme_module_forms": F_SPECIE,
                "enzyme_module_reactions": F_REACTION}.get(gid)
            if f_replace and replace_key in f_replace:
                obj_id = f_replace[replace_key](obj_id)
            obj_id = _prefix_number_id(obj_id)
            if remove_char:
                obj_id = _remove_char_from_id(obj_id)
            # Add object to member list
            member_list.append(obj_id)

        return member_list

    def _parse_enzyme_group(gid, gname, group, enzyme_module_dict,
                            remove_char=None):
        """Parse the various enzyme groups and add to the EnzymeModuleDict."""
        if "_categorized" in gid and CATEGORY_GROUP_RE.search(gid):
            attr = gid[:CATEGORY_GROUP_RE.search(gid).start()]
            gid = attr[:-len("_categorized")]
            member_list = _get_group_members(group, gid, remove_char)
            categorized_dicts[attr][gname].extend(member_list)
        elif "EnzymeModule" in gid:
            # Set attributes stored in the group notes
            group_notes = _parse_notes_dict(group)
            for attr, value in iteritems(group_notes):
                try:
                    value = float(value)
                except ValueError:
                    if attr == "enzyme_net_flux_equation":
                        value = sympify(value)
                finally:
                    enzyme_module_dict[attr] = value
            # Set the name for the EnzymeModuleDict
            enzyme_module_dict["name"] = gname
            # Temporary store object class
            enzyme_module_dict["model"] = gid
        else:
            # Get members for DictList attribute
            member_list = _get_group_members(group, gid, remove_char)
            # Set enzyme module ligands, forms,or reactions attribute list
            enzyme_module_dict[gid] = member_list

    enzyme_module_dicts = []
    for enzyme_module_id, group_list in iteritems(enzyme_group_dict):
        # Make EnzymeModuleDict
        enzyme_module_dict = EnzymeModuleDict(id_or_enzyme=enzyme_module_id)
        categorized_dicts = defaultdict(lambda: defaultdict(list))
        for group in group_list:
            # Get the group ID and name
            gid = _check_required(group, group.getIdAttribute(), "id")
            gid = gid.replace(enzyme_module_id + "_", "")
            gname = group.getName() if group.isSetName() else ""
            # Parse enzyme group
            _parse_enzyme_group(gid, gname, group, enzyme_module_dict,
                                remove_char=remove_char)

        for attr in ["ligands", "forms", "reactions"]:
            attr = "enzyme_module_" + attr + "_categorized"
            enzyme_module_dict[attr] = categorized_dicts[attr]

        enzyme_module_dicts.append(enzyme_module_dict)

    return enzyme_module_dicts


def _read_enzyme_attr_info_from_notes(mass_obj, notes, f_replace=None):
    """Read the enzyme object attributes from the notes into the mass object.

    Applies to the EnzymeModuleForms and EnzymeModuleReaction

    Warnings
    --------
    This method is intended for internal use only.

    """
    for enz_attr in get_subclass_specific_attributes(mass_obj):
        if enz_attr not in notes:
            continue
        # Handle bound dict attributes
        if "bound" in enz_attr:
            attr_value = {}
            for value_str in notes.get(enz_attr).split("; "):
                value, mid = value_str.split(" ")
                if f_replace and F_SPECIE in f_replace:
                    mid = f_replace[F_SPECIE](mid)
                mid = _prefix_number_id(mid)
                attr_value[mid] = value
            setattr(mass_obj, "_" + enz_attr, attr_value)
        else:
            setattr(mass_obj, enz_attr, notes.get(enz_attr))
        # Remove newly set attribute from the notes dict
        del notes[enz_attr]


# -----------------------------------------------------------------------------
# Write SBML
# -----------------------------------------------------------------------------
def write_sbml_model(mass_model, filename, f_replace=None, **kwargs):
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
        If no replacements should be performed, set f_replace={}.
    units: bool
        Should the units be written into the SBMLDocument.
        Default is True.
    local_parameters: bool
        If True, standard reaction parameters [kf, Keq, kr] be written
        as SBML local parameters of the kinetic law. Otherwise. all parameters
        are written as SBML global model parameters. Default is True.

    Notes
    -----
    Custom parameters are always written as SBML global parameters.

    """
    doc = _model_to_sbml(mass_model, f_replace=f_replace, **kwargs)
    if isinstance(filename, string_types):
        # Write to path
        libsbml.writeSBMLToFile(doc, filename)

    elif hasattr(filename, "write"):
        # Write to file handle
        sbml_str = libsbml.writeSBMLToString(doc)
        filename.write(sbml_str)


def _model_to_sbml(mass_model, f_replace=None, **kwargs):
    """Convert MassModel to SBMLDocument.

    Parameters
    ----------
    mass_model: mass.MassModel
        MassModel instance which is written to SBML
    f_replace: dict of replacement functions
        Replacement to apply on identifiers.
    units: bool
        Should the units be written into the SBMLDocument.

    Returns
    -------
    doc: libsbml.SBMLDocument

    Warnings
    --------
    This method is intended for internal use only.

    """
    # Create dict of default kwargs
    default_kwargs = {
        "units": True,
        "local_parameters": True}
    # Use default kwargs or get missing kwargs with their default values
    if kwargs is None:
        kwargs = default_kwargs
    else:
        for key, value in iteritems(default_kwargs):
            if key not in kwargs:
                kwargs[key] = value

    if f_replace is None:
        f_replace = F_REPLACE

    doc, models = _create_doc_and_model_objects(mass_model)
    model, model_fbc, model_groups = models

    # Write the model Meta Information (ModelHistory) into the SBMLDocument
    if hasattr(mass_model, "_sbml"):
        meta = mass_model._sbml
        _write_model_meta_info_to_sbml(model, meta)
        # Set meta annotation and notes
        if "annotation" in meta:
            _sbase_annotations(doc, meta["annotation"])
        if "notes" in meta:
            _sbase_notes_dict(doc, meta["notes"])
    else:
        # TODO write meta information for SBML using massconfiguration options
        pass

    # Write the units information into the SBMLDocument
    if kwargs["units"] and mass_model.units:
        _write_model_units_to_sbml(model, mass_model)

    # Write model notes and annotations into the SBMLDocument
    model_notes = mass_model.notes
    if mass_model.description:
        model_notes["description"] = mass_model.description
    _sbase_notes_dict(model, model_notes)
    _sbase_annotations(model, mass_model.annotation)

    # Write the flux bound parameters into the SBMLDocument
    _write_model_flux_parameters_to_sbml(
        model, mass_model, units=kwargs["units"])

    # Write the compartment information into the SBMLDocument
    _write_model_compartments_to_sbml(model, mass_model)

    # Write the species information into the SBMLDocument
    # Includes MassMetabolites, EnzymeModuleForms, and their concentrations
    _write_model_species_to_sbml(model, mass_model, f_replace=f_replace)
    # Write the boundary conditions into the SBMLDocument
    _write_model_boundary_conditions_to_sbml(
        model, mass_model, f_replace=f_replace)
    # Write the gene information into the SBMLDocument
    _write_model_genes_to_sbml(model_fbc, mass_model, f_replace=f_replace)
    # Write the reaction information into the SBMLDocument
    # Includes MassReactions, EnzymeModuleReactions, and kinetic parameters
    _write_model_reactions_to_sbml(
        model, mass_model, f_replace=f_replace, **kwargs)

    # Write the EnzymeModule specific information into the SBMLDocument
    if isinstance(mass_model, EnzymeModule):
        _write_enzyme_modules_to_sbml(
            model_groups, mass_model, f_replace=f_replace)
    # Write the EnzymeModuleDict specific information into the SBMLDocument
    if mass_model.enzyme_modules:
        for enzyme in mass_model.enzyme_modules:
            _write_enzyme_modules_to_sbml(
                model_groups, enzyme, f_replace=f_replace)

    return doc


def _create_doc_and_model_objects(mass_model):
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
        MassSBMLError(
            "Could not create SBMLDocument due to the following:\n" + str(e))

    # Create model
    model = doc.createModel()
    _check(model, "create model")
    _check(doc.setSBOTerm(SBO_MODELING_FRAMEWORK),
           "set SBO term for modeling framework")

    sbml_plugin_dict = {"fbc": FBC_VERSION}
    if (isinstance(mass_model, EnzymeModule) or mass_model.enzyme_modules):
        sbml_plugin_dict.update({"groups": GROUPS_VERSION})

    plugin_models = []
    for plugin, version_info in iteritems(sbml_plugin_dict):
        # Enable package extension
        _check(
            doc.enablePackage(
                "http://www.sbml.org/sbml/level{0}/version{1}/{2}/version{3}"
                .format(*SBML_LEVEL_VERSION, plugin, version_info),
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
    model_fbc, model_groups = plugin_models

    # Set strictness of fbc package
    _check(model_fbc.setStrict(False), "set fbc plugin strictness to False")

    # Set ID, meta ID, and name
    if mass_model.id is not None:
        _check(model.setIdAttribute(mass_model.id), "set model id")
        _check(model.setMetaId("meta_" + mass_model.id), "set model meta id")
    else:
        _check(model.setMetaId("meta_model"), "set model meta id")
    if mass_model.name is not None:
        _check(model.setName(mass_model.name), "set model name")

    return doc, (model, model_fbc, model_groups)


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
        timestr = meta["created"]
    else:
        # Get the current time and date and set for the model.
        time = datetime.datetime.now()
        timestr = time.strftime('%Y-%m-%dT%H:%M:%S')
    # Set the creation date
    date = libsbml.Date(timestr)
    _check(date, "create the created date")
    _check(history.setCreatedDate(date), "set created date")
    # Set the modified date
    if meta.get("modified", None) and meta["modified"]:
        for k, timestr in enumerate(meta["modified"]):
            date = libsbml.Date(timestr)
            _check(
                history.setModifiedDate(date), "set modified date " + str(k))
    else:
        _check(history.setModifiedDate(date), "set modified date")

    # Make Creator objects and add to ModelHistory
    if meta.get("creators", None):
        for mass_creator in meta["creators"]:
            # Make creator object
            creator = libsbml.ModelCreator()
            _check(creator, "create model creator")
            # Add family name, given name, organization, and email attributes
            for k in ["familyName", "givenName", "organization", "email"]:
                if mass_creator.get(k, None):
                    # Set creator attribute
                    setter = creator.__class__.__dict__.get(
                        "set" + k[0].capitalize() + k[1:])
                    _check(
                        setter(creator, mass_creator[k]), "set creator " + k)
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
        _check(udef.setIdAttribute(id_str),
               "set UnitDefinition id" + _for_id(id_str))
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

    # Get units in MassModel to write
    to_write = list(mass_model.units)
    # Add COBRA unit if it does not exist in model units
    if COBRA_FLUX_UNIT.id not in mass_model.units:
        to_write += [COBRA_FLUX_UNIT]

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

    # Get flux units if units are desired.
    if units:
        flux_udef = model.getUnitDefinition(COBRA_FLUX_UNIT.id)
    else:
        flux_udef = None

    # Create the model flux parameters
    _create_parameter(
        model, pid=LOWER_BOUND_ID, value=min_value, sbo=SBO_DEFAULT_FLUX_BOUND,
        udef=flux_udef, units=units)
    _create_parameter(
        model, pid=UPPER_BOUND_ID, value=max_value, sbo=SBO_DEFAULT_FLUX_BOUND,
        udef=flux_udef, units=units)
    _create_parameter(
        model, pid=ZERO_BOUND_ID, value=0, sbo=SBO_DEFAULT_FLUX_BOUND,
        udef=flux_udef, units=units)
    _create_parameter(
        model, pid=BOUND_MINUS_INF, value=-float("Inf"), sbo=SBO_FLUX_BOUND,
        udef=flux_udef, units=units)
    _create_parameter(
        model, pid=BOUND_PLUS_INF, value=float("Inf"), sbo=SBO_FLUX_BOUND,
        udef=flux_udef, units=units)


def _write_model_parameters_to_sbml(model, mass_model, **kwargs):
    """Write MassModel kinetic parameter information into SBMLDocument.

    The kinetic parameters include the reaction parameters kf, Keq, kr, v, and
    any custom parameters in the model.

    If local parameters is True, the assumption is made that
    all other parameters have already been written into the SBML Document as
    local reaction parameters and will therefore be ignored.

    Warnings
    --------
    This method is intended for internal use only.

    """
    units = kwargs["units"]
    local_parameters = kwargs["local_parameters"]
    for parameter_type, parameter_dict in iteritems(mass_model.parameters):
        # Skip over parameters already written into the model.
        if (local_parameters and parameter_type in ["kf", "Keq", "kr"]) \
           or parameter_type == "Boundary":
            continue
        if parameter_type == "v":
            constant = False
        else:
            constant = True
        # Create kf, Keq, kr, and custom parameters, handling ID corrections
        # for recognized reaction IDs in the parameters
        for pid, value in iteritems(parameter_dict):
            _create_parameter(
                model, pid=pid, value=value, sbo=None, constant=constant,
                udef=None, units=units)


def _create_parameter(sbml_obj, pid, value, sbo=None, constant=True, udef=None,
                      units=None, local_parameters=None):
    """Create a parameter to be written into the SBMLDocument.

    Warnings
    --------
    This method is intended for internal use only.

    """
    # Create parameter and set its ID, value, and whether it is constant
    if local_parameters:
        parameter = sbml_obj.createLocalParameter()
        p_type = " local "
    else:
        parameter = sbml_obj.createParameter()
        p_type = " model "

    _check(parameter, "create" + p_type + "parameter" + _for_id(pid))
    _check(parameter.setIdAttribute(pid),
           "set" + p_type + "parameter id" + _for_id(pid))
    _check(parameter.setValue(value),
           "set" + p_type + "parameter value" + _for_id(pid))

    # Local parameters are always assumed constant
    if not local_parameters:
        _check(parameter.setConstant(constant),
               "set" + p_type + "parameter constant" + _for_id(pid))
    else:
        pass

    # Set SBO term and units if desired.
    if sbo:
        _check(parameter.setSBOTerm(sbo),
               "set" + p_type + "parameter sbo term" + _for_id(pid))
    if units and udef is not None:
        _check(parameter.setUnits(udef.getIdAttribute()),
               "set" + p_type + "parameter units" + _for_id(pid))


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
        _check(compartment.setIdAttribute(cid),
               "set compartment id" + _for_id(cid))
        _check(compartment.setName(name),
               "set compartment name" + _for_id(cid))
        _check(compartment.setConstant(True),
               "set compartment constant" + _for_id(cid))


def _write_model_species_to_sbml(model, mass_model, f_replace=None):
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
        _check(specie.setIdAttribute(mid), "set specie id" + _for_id(mid))
        _check(specie.setName(metabolite.name),
               "set specie name" + _for_id(mid))
        # Use a generic compartment if no species have set compartments
        if not mass_model.compartments:
            cid = list(DEFAULT_COMPARTMENT_DICT)[0]
        else:
            cid = metabolite.compartment
        _check(specie.setCompartment(cid),
               "set specie compartment" + _for_id(mid))
        # Set metabolite initial condition value
        if metabolite.initial_condition is not None:
            _check(
                specie.setInitialConcentration(
                    metabolite.initial_condition),
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
        if isinstance(metabolite, EnzymeModuleForm):
            # Set enzyme information if EnzymeModuleForm
            _write_enzyme_attr_info_to_notes(
                specie, metabolite, f_replace=f_replace)
        else:
            # Otherwise just set the MassMetabolite notes regularly
            _sbase_notes_dict(specie, metabolite.notes)
        # Set metabolite annotation
        _sbase_annotations(specie, metabolite.annotation)


def _write_model_boundary_conditions_to_sbml(model, mass_model,
                                             f_replace=None):
    """Write MassModel boundary condition information into SBMLDocument.

    Warnings
    --------
    This method is intended for internal use only.

    """
    for bmid, bc_value in iteritems(mass_model.boundary_conditions):
        # Create boundary specie and set ID, name, and compartment
        if f_replace and F_SPECIE_REV in f_replace:
            bmid = f_replace[F_SPECIE_REV](bmid)
        specie = model.createSpecies()
        _check(specie, "create model boundary specie" + _for_id(bmid))
        _check(specie.setIdAttribute(bmid),
               "set boundary specie id" + _for_id(bmid))
        _check(specie.setCompartment(str(list(BOUNDARY_COMPARTMENT_DICT)[0])),
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


def _write_model_genes_to_sbml(model_fbc, mass_model, f_replace=None):
    """Write MassModel gene information into SBMLDocument.

    Warnings
    --------
    This method is intended for internal use only.

    """
    for gene in mass_model.genes:
        # Create gene product and set ID
        gid = gene.id
        if f_replace and F_GENE_REV in f_replace:
            gid = f_replace[F_GENE_REV](gid)
        gp = model_fbc.createGeneProduct()
        _check(gp, "create model gene product" + _for_id(gid))
        _check(gp.setIdAttribute(gid), "set gene id " + _for_id(gid))
        # Set gene name and label
        gname = gene.name if gene.name else gid
        _check(gp.setName(gname), "set gene name" + _for_id(gid))
        _check(gp.setLabel(gname), "set gene label" + _for_id(gid))

        # Set gene annotation and notes
        _sbase_annotations(gp, gene.annotation)
        _sbase_notes_dict(gp, gene.notes)


def _write_model_reactions_to_sbml(model, mass_model, f_replace=None,
                                   **kwargs):
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
        _check(reaction.setIdAttribute(rid), "set reaction id" + _for_id(rid))
        _check(reaction.setName(mass_reaction.name),
               "set reaction name" + _for_id(rid))
        _check(reaction.setReversible(mass_reaction.reversible),
               "set reaction reversible" + _for_id(rid))
        if SBML_LEVEL_VERSION[1] <= 1:
            _check(reaction.setFast(False), "set reaction fast" + _for_id(rid))

        # Write specie references into reaction
        reaction = _write_reaction_specie_ref_to_sbml(
            reaction, mass_reaction=mass_reaction, f_replace=f_replace)

        # Set the reaction kinetic law
        _write_reaction_kinetic_law_to_sbml(
            reaction, mass_reaction=mass_reaction, f_replace=f_replace,
            **kwargs)

        # Get plugin for reaction
        reaction_fbc = reaction.getPlugin("fbc")
        _check(reaction_fbc, "get fbc plugin for reaction" + _for_id(rid))

        # Set reaction flux bounds
        flux_udef = model.getUnitDefinition(COBRA_FLUX_UNIT.id)
        reaction_fbc.setLowerFluxBound(
            _create_bound(
                model, reaction=mass_reaction, bound_type="lower_bound",
                f_replace=f_replace, units=kwargs["units"],
                flux_udef=flux_udef))
        reaction_fbc.setUpperFluxBound(
            _create_bound(
                model, reaction=mass_reaction, bound_type="upper_bound",
                f_replace=f_replace, units=kwargs["units"],
                flux_udef=flux_udef))

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

    # Write the rate law parameters into the SBMLDocument
    _write_model_parameters_to_sbml(
        model, mass_model, f_replace=f_replace, **kwargs)


def _write_reaction_specie_ref_to_sbml(reaction, mass_reaction,
                                       f_replace=None):
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
               "create specie reference " + sid + _for_id(
                   reaction.getIdAttribute()))
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


def _write_reaction_kinetic_law_to_sbml(reaction, mass_reaction,
                                        f_replace=None, **kwargs):
    """Write MassReaction kinetic law information into SBMLDocument.

    Warnings
    --------
    This method is intended for internal use only.

    """
    units = kwargs["units"]
    local_parameters = kwargs["local_parameters"]

    rid = reaction.getIdAttribute()
    rate_eq = strip_time(mass_reaction.rate)
    # If ID replacements were performed earlier then apply the ID
    # replacements for metabolite and parameter arguments in rate law also.
    id_subs = {}
    for arg in list(rate_eq.atoms(Symbol)):
        new_arg = str(arg)
        if mass_reaction.id in new_arg:
            id_subs.update({arg: new_arg})
            continue
        # Ensure the metabolite argument exists in the model
        try:
            new_arg = mass_reaction.model.metabolites.get_by_id(new_arg)
            new_arg = str(new_arg)
        except KeyError:
            if new_arg in mass_reaction.model.custom_parameters:
                new_arg = None

        if new_arg is None:
            continue

        # Account for metabolites in the rate law that
        # are not considered reactants or products via modifiers
        rxn_mets = [str(met) for met in mass_reaction.metabolites]
        rxn_mets += [mass_reaction.boundary_metabolite]
        # Modify the ID if metabolite participates in the reaction
        if new_arg in rxn_mets and f_replace and F_SPECIE_REV in f_replace:
            new_arg = f_replace[F_SPECIE_REV](new_arg)
        elif new_arg not in rxn_mets:
            # Modify the ID and create a modifier species
            if f_replace and F_SPECIE_REV in f_replace:
                new_arg = f_replace[F_SPECIE_REV](new_arg)
            msref = reaction.createModifier()
            _check(msref, "create modifier specie " + new_arg + _for_id(rid))
            _check(msref.setSpecies(new_arg),
                   "set modifier specie species " + new_arg + _for_id(rid))

        id_subs.update({arg: new_arg})
    # Make xml string of rate equation via sympy conversion to Math ML
    math_xml_str = _create_math_xml_str_from_sympy_expr(rate_eq.subs(id_subs))
    # Create kinetic law as AST math and write into the SBMLDocument
    kinetic_law = reaction.createKineticLaw()
    _check(kinetic_law, "create kinetic law" + _for_id(rid))
    _check(kinetic_law.setMath(libsbml.readMathMLFromString(math_xml_str)),
           "set math on kinetic law" + _for_id(rid))

    # Set local reaction parameters if desired
    if local_parameters:
        for pid, value in iteritems(mass_reaction.parameters):
            _create_parameter(
                kinetic_law, pid=pid, value=value, sbo=None, udef=None,
                units=units, local_parameters=local_parameters)


def _write_reaction_gpr_to_sbml(reaction_fbc, mass_reaction, f_replace=None):
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


def _write_enzyme_modules_to_sbml(model_groups, mass_obj, f_replace=None):
    """Write EnzymeModule and EnzymeModuleDict information into SBMLDocument.

    Utilizes the groups package extension.

    Warnings
    --------
    This method is intended for internal use only.

    """
    # Create group to represent EnzymeModule/EnzymeModuleDict information
    notes = {}
    # Create groups for attributes, add them to enzyme group members list
    enzyme_group_members = []
    for attr_name, default in iteritems(_ORDERED_ENZYMEMODULE_DICT_DEFAULTS):
        attr_value = getattr(mass_obj, attr_name, None)
        # Ignore attributes set at default and the id, name, and S matrix.
        gid = "_".join((mass_obj.id, attr_name))
        if attr_name in["id", "name", "S", "model"] or attr_value == default:
            continue
        elif isinstance(attr_value, list):
            # Write attribute with list value as a group.
            enzyme_group_members += [
                _write_group_to_sbml(
                    model_groups, gid=gid, members=attr_value,
                    f_replace=f_replace)]
        elif isinstance(attr_value, dict):
            # Write attribute with dict containing lists as a group of groups.
            members = []
            # Create enzyme group members (as groups) and add to list
            for i, (category, values) in enumerate(iteritems(attr_value)):
                members += [
                    _write_group_to_sbml(
                        model_groups, gid=gid + "_Category" + str(i),
                        gname=category, members=values, f_replace=f_replace)]
            enzyme_group_members += [
                _write_group_to_sbml(
                    model_groups, gid=gid, members=members,
                    f_replace=f_replace)]
        elif isinstance(attr_value, Basic):
            # Write attribute with sympy equation to the group notes.
            notes.update({attr_name: str(attr_value.rhs)})
        else:
            # Add attribute to the group notes.
            notes.update({attr_name: str(attr_value)})
    # Create enzyme group
    gid = "_".join((mass_obj.id, mass_obj.__class__.__name__))
    enzyme_group = _write_group_to_sbml(
        model_groups, gid=gid, gname=mass_obj.name, kind="collection",
        members=enzyme_group_members, f_replace=f_replace)
    _sbase_notes_dict(enzyme_group, notes)


def _write_enzyme_attr_info_to_notes(sbml_obj, mass_obj, f_replace=None):
    """Write the enzyme object attribute into SBMLDocument in the notes.

    Applies to the EnzymeModuleForms and EnzymeModuleReaction

    Warnings
    --------
    This method is intended for internal use only.

    """
    # Get original notes from mass object.
    notes = mass_obj.notes.copy()
    # Get the subclass specific attributes to write into the notes.
    attributes = get_subclass_specific_attributes(mass_obj)
    default_values = [{} if "bound" in attr_name else ""
                      for attr_name in attributes]
    for attr_name, default_value in zip(attributes, default_values):
        attr_value = getattr(mass_obj, attr_name)
        # Skip attributes equal to their defaults
        if attr_value == default_value and attr_name != "enzyme_module_form":
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
            notes.update({
                attr_name: _make_bound_attr_str_repr(bound_dict)})
        else:
            # Add to notes to write
            notes.update({attr_name: attr_value})
    # Write the notes to the SBML object.
    _sbase_notes_dict(sbml_obj, notes)


def _write_group_to_sbml(model_groups, gid, gname="", kind="collection",
                         members=None, f_replace=None):
    """Write an SBML group into SBMLDocument, and return the group.

    Additional group information is written into the notes attribute of the
    SBML object.

    Warnings
    --------
    This method is intended for internal use only.

    """
    group = model_groups.createGroup()
    _check(group, "create group" + _for_id(gid))
    _check(group.setIdAttribute(gid), "set group id" + _for_id(gid))
    _check(group.setName(gname), "set group name" + _for_id(gid))
    _check(group.setKind(kind), "set group kind" + _for_id(gid))

    if members:
        for member in members:
            # Get member ID and name
            if not isinstance(member, libsbml.SBase):
                mid = member.id
                mname = member.name
            else:
                mid = member.getIdAttribute()
                mname = ""
                if member.isSetName():
                    mname = member.getName()

            # ID replacements
            m_type = str(member.__class__.__name__)
            if m_type in ["MassReaction", "EnzymeModuleReaction"]:
                if f_replace and F_REACTION_REV in f_replace:
                    mid = f_replace[F_REACTION_REV](mid)
            if m_type in ["MassMetabolite", "EnzymeModuleForm"]:
                if f_replace and F_SPECIE_REV in f_replace:
                    mid = f_replace[F_SPECIE_REV](mid)
            if m_type in ["Gene"]:
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
# Notes and Annotations
# -----------------------------------------------------------------------------
def _sbase_annotations(sbase, annotation):
    """Set SBase annotations based on mass annotations.

    Annotations are written for the SBML object using the cobra annotation
    writing function: `cobra.io._sbase_annotations`. Any raised CobraSBMLError
    will be reclassified and raised as a MassSBMLError.

    Parameters
    ----------
    sbase : libsbml.SBase
        SBML object to annotate
    annotation : mass annotation structure
        mass object with annotation information

    """
    try:
        _cobra_sbase_annotations(sbase, annotation)
    # Reclassify cobra error as a mass error
    except CobraSBMLError as e:
        raise MassSBMLError(e)


def _parse_notes_dict(sbase):
    """Read notes from SBML object and create dictionary of MASS notes.

    Expected notes format corresponds with the output format of the the
    cobra notes writing function: `cobra.io.sbml_sbase_notes_dict`

    Parameters
    ----------
    sbase : libsbml.SBase

    Returns
    -------
    dict of notes

    Warnings
    --------
    This method is intended for internal use only.

    """
    notes = {}
    notes_str = sbase.getNotesString()
    if notes_str:
        # Get starting and ending positions for notes
        start_list = (match.end() for match in re.finditer(r"<p>", notes_str))
        end_list = (match.start() for match in re.finditer(r"</p>", notes_str))
        # Iterate through notes string to get keys and values for notes
        for s, e in zip(start_list, end_list):
            match = NOTES_RE.search(notes_str[s:e])
            if match is not None:
                k, v = match.groups()
                k, v = (k.strip(), v.strip())
            else:
                LOGGER.warning(
                    "No key provided for note in %s '%s'. The note will be "
                    "added to the notes attribute with key 'NOTES'.",
                    sbase.__class__.__name__, sbase.getIdAttribute())
                k, v = ("NOTES", [])
                if k in notes:
                    v += notes["NOTES"]
                v += [notes_str[s:e]]
            notes.update({k: v})

    return notes


# -----------------------------------------------------------------------------
# Validation
# -----------------------------------------------------------------------------
def validate_sbml_model(filename, check_model=True, internal_consistency=True,
                        check_units_consistency=False,
                        check_modeling_practice=False, **kwargs):
    """Validate SBML model and returns the model along with a list of errors.

    Parameters
    ----------
    filename : str
        The filename (or SBML string) of the SBML model to be validated.
    check_model: boolean {True, False}
        Check some basic model properties. Default is True
    internal_consistency: bool
        Check internal consistency. Default is True.
    check_units_consistency: bool
        Check consistency of units. Default is False.
    check_modeling_practice: bool
        Check modeling practice. Default is False.

    Returns
    -------
    (model, errors)
    model : :class:`~mass.core.massmodel.MassModel` object
        The mass model if the file could be read successfully or None
        otherwise.
    errors : dict
        Warnings and errors grouped by their respective types.

    Raises
    ------
    MassSBMLError

    """
    # Errors and warnings are grouped based on their type. SBML_* types are
    # from the libsbml validator. MASS_* types are from the masspy SBML
    # parser.
    keys = (
        "SBML_FATAL",
        "SBML_ERROR",
        "SBML_SCHEMA_ERROR",
        "SBML_WARNING",

        "MASS_FATAL",
        "MASS_ERROR",
        "MASS_WARNING",
        "MASS_CHECK",
    )
    errors = {key: [] for key in keys}
    # [1] libsbml validation
    doc = _get_doc_from_filename(filename)

    # Set checking of units and modeling practice
    doc.setConsistencyChecks(libsbml.LIBSBML_CAT_UNITS_CONSISTENCY,
                             check_units_consistency)
    doc.setConsistencyChecks(libsbml.LIBSBML_CAT_MODELING_PRACTICE,
                             check_modeling_practice)
    # Check consistency
    if internal_consistency:
        doc.checkInternalConsistency()
    doc.checkConsistency()

    for k in range(doc.getNumErrors()):
        e = doc.getError(k)
        msg = _error_string(e, k=k)
        sev = e.getSeverity()
        if sev == libsbml.LIBSBML_SEV_FATAL:
            errors["SBML_FATAL"].append(msg)
        elif sev == libsbml.LIBSBML_SEV_ERROR:
            errors["SBML_ERROR"].append(msg)
        elif sev == libsbml.LIBSBML_SEV_SCHEMA_ERROR:
            errors["SBML_SCHEMA_ERROR"].append(msg)
        elif sev == libsbml.LIBSBML_SEV_WARNING:
            errors["SBML_WARNING"].append(msg)

    # [2] masspy validation (check that SBML can be read into model)
    # all warnings generated while loading will be logged as errors
    log_stream = StringIO()
    stream_handler = logging.StreamHandler(log_stream)
    formatter = logging.Formatter('%(levelname)s:%(message)s')
    stream_handler.setFormatter(formatter)
    stream_handler.setLevel(logging.INFO)
    LOGGER.addHandler(stream_handler)
    LOGGER.propagate = False

    try:
        # read model and allow additional parser arguments
        model = _sbml_to_model(doc, **kwargs)
    except MassSBMLError as e:
        errors["MASS_ERROR"].append(str(e))
        return None, errors
    except Exception as e:
        errors["MASS_FATAL"].append(str(e))
        return None, errors

    mass_errors = log_stream.getvalue().split("\n")
    for mass_error in mass_errors:
        tokens = mass_error.split(":")
        error_type = tokens[0]
        error_msg = ":".join(tokens[1:])

        if error_type == "WARNING":
            errors["MASS_WARNING"].append(error_msg)
        elif error_type == "ERROR":
            errors["MASS_ERROR"].append(error_msg)

    # remove stream handler
    LOGGER.removeHandler(stream_handler)
    LOGGER.propagate = True

    # [3] additional model tests
    if check_model:
        # TODO add checks
        pass

    for key in ["SBML_FATAL", "SBML_ERROR", "SBML_SCHEMA_ERROR"]:
        if errors[key]:
            LOGGER.error("SBML errors in validation, check error log "
                         "for details.")
            break
    for key in ["SBML_WARNING"]:
        if errors[key]:
            LOGGER.error("SBML warnings in validation, check error log "
                         "for details.")
            break
    for key in ["MASS_FATAL", "MASS_ERROR"]:
        if errors[key]:
            LOGGER.error("MASS errors in validation, check error log "
                         "for details.")
            break
    for key in ["MASS_WARNING", "MASS_CHECK"]:
        if errors[key]:
            LOGGER.error("MASS warnings in validation, check error log "
                         "for details.")
            break

    return model, errors


def _check_required(sbase, value, attribute):
    """Get required attribute from SBase.

    Parameters
    ----------
    sbase: libsbml.SBase
    value: existing value
    attribute: name of attribute

    Returns
    -------
    attribute value (or value if already set)

    """
    if (value is None) or (value == ""):
        msg = str("Required attribute '{0}' cannot be found or parsed in "
                  "'{1}'".format(attribute, sbase))
        if hasattr(sbase, "getId") and sbase.getId():
            msg += " with id '{0}'".format(sbase.getId())
        elif hasattr(sbase, "getName") and sbase.getName():
            msg += " with name '{0}'".format(sbase.getName())
        elif hasattr(sbase, "getMetaId") and sbase.getMetaId():
            msg += " with metaId '{0}'".format(sbase.getMetaId())
        raise MassSBMLError(msg)
    if attribute == "id":
        if not libsbml.SyntaxChecker.isValidSBMLSId(value):
            LOGGER.error("'%s' is not a valid SBML 'SId'.", value)

    return value


def _check(value, message):
    """Check the libsbml return value and log error messages.

    If value is None, logs an error messaage constructed using 'message' and
    then exists with status code 1. If 'value' is an integer, it assumes it is
    an libSBML return status code. If the code value is
    LIBSBML_OPERATION_SUCCESS, returns without further action; if it is not,
    logs an error message constructed using 'message' along with text from
    libSBML explaining the meaning of the code, and exits with status code 1.

    Warnings
    --------
    This method is intended for internal use only.

    """
    if value is None:
        LOGGER.error("Error: LibSBML returned a null value trying to "
                     "<%s>.", message)
    if isinstance(value, integer_types)\
       and value != libsbml.LIBSBML_OPERATION_SUCCESS:
        LOGGER.error("Error encountered trying to  <%s>.", message)
        LOGGER.error("LibSBML error code %s: %s", str(value),
                     libsbml.OperationReturnValue_toString(value).strip())


def _for_id(sid):
    """Return a string specifying the object id for logger messages."""
    return " for '{0}'".format(sid)
