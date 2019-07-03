# -*- coding: utf-8 -*-
"""TODO Module Docstrings."""
from __future__ import absolute_import

import datetime
import logging
import re
import traceback
from collections import defaultdict
from io import StringIO

import libsbml

from six import integer_types, iteritems, itervalues, string_types

from sympy import Basic, Symbol, SympifyError, mathml, sympify

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
from mass.core.units import Unit, UnitDefinition, _SBML_BASE_UNIT_KINDS_DICT
from mass.enzyme_modules.enzyme_module import EnzymeModule
from mass.enzyme_modules.enzyme_module_dict import (
    EnzymeModuleDict, _ORDERED_ENZYMEMODULE_DICT_DEFAULTS)
from mass.enzyme_modules.enzyme_module_form import (
    EnzymeModuleForm, _make_bound_attr_str_repr)
from mass.enzyme_modules.enzyme_module_reaction import EnzymeModuleReaction
from mass.exceptions import MassSBMLError
from mass.util.expressions import strip_time
from mass.util.util import _check_kwargs, get_subclass_specific_attributes

LOGGER = logging.getLogger(__name__)
# -----------------------------------------------------------------------------
# Defaults and constants for writing SBML
# -----------------------------------------------------------------------------
MASSCONFIGURATION = MassConfiguration()
# For SBML level and core version, and package extension versions
SBML_LEVEL_VERSION = (3, 1)
FBC_VERSION = 2
GROUPS_VERSION = 1

# Precompiled regex
CHAR_RE = re.compile(r"(?P<char>_Char\d*_)")           # For ASCII Char removal
MASS_MOIETY_RE = re.compile(r"\[\S+\]")                # For mass moieties
SBML_MOIETY_RE = re.compile(r"Moiety")                 # For SBML moieties
CATEGORY_GROUP_RE = re.compile(r"_Category\d$")        # For category groups
RATE_CONSTANT_RE = re.compile(r"k_(\S*)_(fwd|rev)\Z")  # For kf & kr
NOTES_RE = re.compile(r"\s*(\w+\s*\w*)\s*:\s*(.+)", re.DOTALL)  # For notes
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


def _remove_char_from_id(sid):
    """Remove ASCII characters from an identifier."""
    for match in CHAR_RE.finditer(sid):
        char_key = match.group("char")
        if ID_ASCII_REPLACE.get(char_key):
            replacement_str = ID_ASCII_REPLACE[char_key]
        else:
            replacement_str = "_"
        sid = re.sub(char_key, replacement_str, sid)
    return sid


def _get_corrected_id(item, f_items, obj_type, remove_char):
    """Return the id after performing any necessary corrections."""
    # Get id as string
    if isinstance(item, string_types):
        id_str = item
    else:
        id_str = _check_required(item, item.getIdAttribute(), "id")
    # Perform replacements
    if len(f_items) == 3:
        id_str = re.sub(*f_items)
    else:
        f_replace, f_key = f_items
        if f_replace and f_key in f_replace:
            id_str = f_replace[f_key](id_str)
    id_str = _prefix_number_id(id_str, obj_type)
    if remove_char:
        id_str = _remove_char_from_id(id_str)
    return id_str


# -----------------------------------------------------------------------------
# MathML
# -----------------------------------------------------------------------------
SBML_MATH_PACKAGE_NAME = "l{0}v{1}extendedmath".format(*SBML_LEVEL_VERSION)
MATHML_XML_NS = '<math xmlns="http://www.w3.org/1998/Math/MathML">{0}</math>'
REPLACE_TIME_RE = re.compile(r"\>\<ci\>t\<\/ci\>\<")
REMOVE_MML_TAGS_RE = re.compile(r"\<mml:.*?\>|\<\/mml:.*?\>")


def _create_math_xml_str_from_sympy_expr(sympy_equation):
    """Create a MathML string from a sympy expression to be parsed by libsbml.

    This also requires removing 'mml' tags from XML representation of
        mathematical equations.

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
    # Replace underscores before conversion
    underscore_replace = {str(arg): str(arg).replace("_", "&")
                          for arg in sympy_equation.atoms(Symbol)}
    sympy_equation = sympy_equation.subs(underscore_replace)
    # Convert equation into MathML and remove tags that libsbml cannot parse
    mathml_xml = mathml(sympy_equation)
    mathml_xml = MATHML_XML_NS.format(mathml_xml.replace("&", "_"))
    mathml_xml = REMOVE_MML_TAGS_RE.sub("", mathml_xml)

    # Fix the time symbol
    if REPLACE_TIME_RE.search(mathml_xml):
        time_symbol = '><csymbol encoding="text" definitionURL=' +\
                      '"http://www.sbml.org/sbml/symbols/time">' +\
                      't</csymbol><'
        mathml_xml = REPLACE_TIME_RE.sub(time_symbol, mathml_xml)
    return mathml_xml


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

    The parser tries to fall back to information in notes dictionaries
    if information is not available in the FBC packages, e.g.,
    CHARGE, FORMULA on species, or GENE_ASSOCIATION, SUBSYSTEM on reactions.

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
        Default is True.
    remove_char: bool
        If True, remove ASCII characters from identifiers.
        Default is True.
    stop_on_conversion_fail: bool
        Whether to stop trying to load the model if a conversion process fails.
        If False, then the loading of the model will be attempted anyways,
        despite a potential loss of information.
        Default is True.

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
    stop_on_conversion_fail: bool
        Whether to stop trying to load the model if a conversion process fails.
        If False, then the loading of the model will be attempted anyways,
        despite a potential loss of information.
        Default is True.

    Returns
    -------
    mass.core.MassModel

    Warnings
    --------
    This method is intended for internal use only.

    """
    default_kwargs = {
        "number": float,
        "set_missing_bounds": True,
        "remove_char": True,
        "stop_on_conversion_fail": True}
    # Check kwargs
    kwargs = _check_kwargs(default_kwargs, kwargs)

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
    mass_model._sbml = _read_model_meta_info_from_sbml(doc, model)

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
    metabolites, boundary_conditions = _read_model_species_from_sbml(model,
                                                                     f_replace,
                                                                     **kwargs)
    mass_model.add_metabolites(metabolites)

    # Read the MassModel gene information from the SBMLDocument.
    genes = _read_model_genes_from_sbml(model, model_fbc, f_replace, **kwargs)
    mass_model.genes += genes

    # Read the MassModel reaction information from the SBMLDocument
    # Includes MassReactions, EnzymeModuleReactions, kinetic laws,
    # and GPR information.
    reactions, local_params, custom_rates = _read_model_reactions_from_sbml(
        model, mass_model.metabolites, f_replace, **kwargs)
    mass_model.add_reactions(reactions)
    # Add the parameters and boundary conditions to the model
    # Read the MassModel parameter information from the SBMLDocument
    parameters = _read_global_parameters_from_sbml(model, mass_model.reactions,
                                                   f_replace, **kwargs)
    parameters.update(local_params)
    parameters.update(boundary_conditions)
    mass_model.update_parameters(parameters, verbose=False)

    # Add custom rates
    mass_model.update_custom_rates(custom_rates)

    # Add enzyme groups and other parsed groups to the model.
    # TODO Revisit when cobra groups implemented, store in mass_groups
    if model_groups:
        enzyme_groups, mass_groups = _read_model_groups_from_sbml(model_groups)
        # Read the enzyme groups into EnzymeModuleDicts
        enzyme_module_dicts_list = _read_enzyme_modules_from_sbml(
            enzyme_groups, f_replace, **kwargs)

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
    stop_on_conversion_fail = kwargs.get("stop_on_conversion_fail")
    model = doc.getModel()
    if model is None:
        raise MassSBMLError("No SBML model detected in file.")

    def check_conversion_success(result, msg, stop_on_conversion_fail):
        """Check whether a conversion was successful, and handle fails."""
        if result != libsbml.LIBSBML_OPERATION_SUCCESS:
            msg = "Conversion " + msg + " failed."
            if stop_on_conversion_fail:
                raise Exception(msg)
            else:
                LOGGER.warning(
                    "%s Loading of model will continue, but may yield "
                    "unexpected results and some information loss.", msg)

    def convert_package_extensions(plugin, new_version):
        """Convert the plugin package extension to a different version."""
        model_plugin = model.getPlugin(plugin)
        if model_plugin:
            cur_version = model_plugin.getPackageVersion()
            if cur_version != new_version:
                cur_version = "{0}-v{1}".format(plugin, cur_version)
                new_version = "{0}-v{1}".format(plugin, new_version)
                LOGGER.warning(
                    "Loading SBML with %s (models should be encoded using "
                    "%s).", cur_version, new_version)
                conversion_properties = libsbml.ConversionProperties()
                conversion_properties.addOption(
                    "convert {0} to {1}".format(
                        cur_version.replace("-", " "),
                        new_version.replace("-", " ")), True,
                    "Convert {0} to {1}".format(
                        cur_version.upper(), new_version.upper()))
                # Convert and check if successful
                result = doc.convert(conversion_properties)
                check_conversion_success(
                    result, "of {0} to {1} ".format(cur_version, new_version),
                    stop_on_conversion_fail)
        else:
            LOGGER.warning("Model does not contain SBML %s package", plugin)

        return model_plugin

    # Determine whether conversion to different SBML Levels & versions needed.
    current_level_version = (doc.getLevel(), doc.getVersion())
    if current_level_version != SBML_LEVEL_VERSION:
        # Create conversion message for logger.
        msg = "from Level {0} Version {1} to Level {2} Version {3}".format(
            *current_level_version, *SBML_LEVEL_VERSION)
        LOGGER.warning(
            "SBML model is written as Level %s Version %s, will attempt to "
            "convert the SBML model %s before parsing", *current_level_version,
            msg)
        # Attempt conversion
        result = doc.setLevelAndVersion(*SBML_LEVEL_VERSION)
        check_conversion_success(result, msg, stop_on_conversion_fail)

    model_fbc = convert_package_extensions("fbc", FBC_VERSION)
    model_groups = convert_package_extensions("groups", GROUPS_VERSION)

    return doc, (model, model_fbc, model_groups)


def _read_model_meta_info_from_sbml(doc, model):
    """Read SBML meta information to create a dict of meta information.

    Warnings
    --------
    This method is intended for internal use only.

    """
    # Make meta information as a dict
    meta = {
        "model.id": model.getIdAttribute(),
        "level": model.getLevel(),
        "version": model.getVersion()}

    # Read Model history, creation date, and creators
    created = ""
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
        # Make a dict of creator attributes and values
        for creator in history.getListCreators():
            creator_dict = {}
            # Make a dict for the creator attributes initialized to None
            for k in ["familyName", "givenName", "organization", "email"]:
                attr = k[0].capitalize() + k[1:]
                if getattr(creator, "isSet" + attr)():
                    creator_dict[k] = getattr(creator, "get" + attr)()
                else:
                    creator_dict[k] = None
            creators.append(creator_dict)

    # Add meta information to the dict
    meta["creators"] = creators
    meta["created"] = created
    meta["modified"] = modified
    meta["notes"] = _parse_notes_dict(doc)
    meta["annotation"] = _parse_annotations(doc)

    # Get version and package info
    meta["info"] = "<{0}> SBML L{1}V{2}".format(
        meta["model.id"], meta["level"], meta["version"])

    packages = {}
    for k in range(doc.getNumPlugins()):
        plugin = doc.getPlugin(k)
        k, v = (plugin.getPackageName(), plugin.getPackageVersion())
        packages[k] = v
        meta["info"] += ", {0}-v{1}".format(k, v)
        if k not in [SBML_MATH_PACKAGE_NAME, "fbc", "groups"]:
            LOGGER.warning(
                "SBML package '%s' not supported by masspy. Therefore the "
                "package information is not parsed.", k)

    meta["packages"] = packages

    return meta


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


def _read_global_parameters_from_sbml(model, reactions, f_replace, **kwargs):
    """Read the SBML model parameters and return them in dicts.

    Warnings
    --------
    This method is intended for internal use only.

    """
    bound_pid_list = [
        LOWER_BOUND_ID, UPPER_BOUND_ID, ZERO_BOUND_ID, BOUND_MINUS_INF,
        BOUND_PLUS_INF, "_lower_bound", "_upper_bound"]

    parameters_dict = {}
    for parameter in model.getListOfParameters():
        # Get parameter ID and value
        pid = _check_required(
            parameter, parameter.getIdAttribute(), "id")
        # Check if parameter is a flux bound
        if any([x in pid for x in bound_pid_list]):
            continue

        try:
            p_type, rid = pid.split("_", 1)
            if re.search(p_type, rid):
                p_type = re.sub(rid, "", p_type)
            rid = _get_corrected_id(rid, (f_replace, F_REACTION),
                                    None, kwargs.get("remove_char"))
            rid = getattr(reactions.get_by_id(rid), "id")
        except (ValueError, KeyError):
            pass
        else:
            pid = p_type + "_" + rid
        finally:
            value = parameter.getValue()

        parameters_dict[pid] = value

    return parameters_dict


def _read_model_species_from_sbml(model, f_replace, **kwargs):
    """Read SBML species and boundary conditions and return them.

    MassMetabolites and EnzymeModuleForms are created and returned in a list.
    Boundary conditions are returned as a dict.

    Warnings
    --------
    This method is intended for internal use only.

    """
    metabolite_list = []
    boundary_conditions = {}
    if model.getNumSpecies() == 0:
        LOGGER.warning("No metabolites in model.")

    for specie in model.getListOfSpecies():
        # Check the identifier and make corrections if necessary
        sid = _get_corrected_id(specie, (f_replace, F_SPECIE),
                                "Metabolite", kwargs.get("remove_char"))
        cid = specie.getCompartment()
        # Do not read species with boundary compartments as MassMetabolites
        # to be added but as boundary conditions to be added
        if cid == next(iter(MASSCONFIGURATION.boundary_compartment)):
            boundary_conditions[sid] = _read_specie_initial_value(model,
                                                                  specie,
                                                                  sid)
            continue
        # Set MassMetabolite id, name, and fixed attributes
        met = MassMetabolite(sid)
        met.name = specie.getName()
        met.fixed = specie.getConstant()
        # Set the MassMetabolite compartment if it is not the default
        if cid != next(iter(MASSCONFIGURATION.default_compartment)):
            met.compartment = cid

        met.initial_condition = _read_specie_initial_value(model, specie, sid)

        # Using notes information, determine whether the metabolite should
        # be added to the model as an EnzymeModuleForm.
        notes = _parse_notes_dict(specie)
        if "enzyme_module_id" in notes:
            met = EnzymeModuleForm(met)
            _read_enzyme_attr_info_from_notes(met, notes, f_replace)

        # Add the notes and annotations to the metabolite
        met.notes = notes
        met.annotation = _parse_annotations(specie)
        # Add the formula and charge to the metabolite
        met.formula, met.charge = _read_specie_formula_charge_from_sbml(specie,
                                                                        met)

        # Append metabolite to list of metabolites
        metabolite_list.append(met)

    return metabolite_list, boundary_conditions


def _read_specie_initial_value(model, specie, sid):
    """Read the specie as a boundary condition and return a dict.

    Warnings
    --------
    This method is intended for internal use only.

    """
    initial_condition = None
    # Try to set initial concentration or amount depending on rate law.
    if MASSCONFIGURATION.include_compartments_in_rates:
        attributes = ("InitialAmount", "InitialConcentration")
        msg = "include compartments"
    else:
        attributes = ("InitialConcentration", "InitialAmount")
        msg = "exclude compartments"
    for i, attr in enumerate(attributes):
        if getattr(specie, "isSet" + attr)():
            initial_condition = getattr(specie, "get" + attr)()
            if i == len(attributes) - 1:
                LOGGER.warning(
                    "Rate expressions are set to %s, but the initial "
                    "condition is an %s%s in the SBML model specie. ",
                    msg, attr, _for_id(sid))
            break

    if initial_condition is None:
        assignment_rule = model.getAssignmentRule(specie.getIdAttribute())
        if assignment_rule is not None:
            initial_condition = assignment_rule.getFormula()
            # Perform substitution for power law operations to sympify rate
            for match in KLAW_POW_RE.finditer(initial_condition):
                old = match.group(0)
                new = "(({0})**{1})".format(
                    match.group("arg"), match.group("exp"))
                initial_condition = initial_condition.replace(old, new)
            # Try to sympify the reaction rate
            try:
                initial_condition = sympify(initial_condition)
            except SympifyError as e:
                raise MassSBMLError(e)

    return initial_condition


def _read_specie_formula_charge_from_sbml(specie, metabolite):
    """Read the specie formula and charge values and return them.

     Warnings
    --------
    This method is intended for internal use only.

    """
    formula, charge = (None, None)
    # Get the fbc specie and set the charge and formula
    specie_fbc = specie.getPlugin("fbc")
    if specie_fbc:
        if specie_fbc.isSetChemicalFormula():
            formula = _f_moiety_formula(specie_fbc.getChemicalFormula())
        if specie_fbc.isSetCharge():
            charge = specie_fbc.getCharge()
    else:
        # Handle legacy species
        if specie.isSetCharge():
            LOGGER.warning(
                "Use of the species charge attribute is discouraged, use "
                "fbc:charge instead: %s", specie)
            charge = specie.getCharge()
        # Fallback to notes information for legacy charge
        elif "CHARGE" in metabolite.notes:
            LOGGER.warning(
                "Use of CHARGE in the notes element is discouraged, use "
                "fbc:charge instead: %s", specie)
            try:
                charge = int(metabolite.notes["CHARGE"])
            except ValueError:
                pass
            else:
                del metabolite.notes['CHARGE']
        # Fallback to notes information for legacy formula
        if 'FORMULA' in metabolite.notes:
            LOGGER.warning(
                "Use of FORMULA in the notes element is discouraged, use "
                "fbc:chemicalFormula instead: %s", specie)
            formula = metabolite.notes['FORMULA']
            del metabolite.notes['FORMULA']

    return formula, charge


def _read_model_genes_from_sbml(model, model_fbc, f_replace, **kwargs):
    """Read the SBML genes, returning them as a list.

    cobra Genes are created and returned in a list.

    Warnings
    --------
    This method is intended for internal use only.

    """
    gene_list = []
    # Read GPR via FBC package
    if model_fbc:
        for gp in model_fbc.getListOfGeneProducts():
            # Check the identifier and make corrections if necessary
            gid = _get_corrected_id(gp, (f_replace, F_GENE),
                                    None, kwargs.get("remove_char"))

            # Set the gene ID and name
            gene = Gene(gid)
            gene.name = gp.getName() if gp.isSetName() else gid
            # Parse notes and annotations
            gene.notes = _parse_notes_dict(gp)
            gene.annotation = _parse_annotations(gp)
            # Add gene to list of genes
            gene_list.append(gene)
    else:
        # Parse legacy GPRs
        for reaction in model.getListOfReactions():
            notes = _parse_notes_dict(reaction)
            if "GENE ASSOCIATION" in notes:
                gpr = notes['GENE ASSOCIATION']
            elif "GENE_ASSOCIATION" in notes:
                gpr = notes['GENE_ASSOCIATION']
            else:
                gpr = ''

            if gpr:
                gpr = gpr.replace("(", ";")
                gpr = gpr.replace(")", ";")
                gpr = gpr.replace("or", ";")
                gpr = gpr.replace("and", ";")
                gids = [t.strip() for t in gpr.split(';')]

                # Create missing genes
                for gid in gids:
                    gid = _get_corrected_id(gid, (f_replace, F_GENE),
                                            None, kwargs.get("remove_char"))

                    if gid not in gene_list:
                        gene = Gene(gid)
                        gene.name = gid
                        gene_list.append(gene)
    return gene_list


def _read_model_reactions_from_sbml(model, metabolites, f_replace, **kwargs):
    """Read the SBML reactions, returning them as a list.

    MassReactions and EnzymeModuleReactions are created and returned in a list.

    Warnings
    --------
    This method is intended for internal use only.

    """
    reaction_list = []
    local_parameters_dict = {}
    custom_rates_dict = {}
    if model.getNumReactions() == 0:
        LOGGER.warning("No reactions in model.")

    for reaction in model.getListOfReactions():
        # Check the identifier
        rid = _get_corrected_id(reaction, (f_replace, F_REACTION),
                                "Reaction", kwargs.get("remove_char"))
        # Set the reaction ID, name, and reversible
        mass_reaction = MassReaction(rid)
        mass_reaction.name = reaction.getName()
        mass_reaction.reversible = reaction.getReversible()

        # Read reaction species and determine stoichiometry
        metabolite_stoichiometry = _read_reaction_species_from_sbml(reaction,
                                                                    f_replace,
                                                                    **kwargs)
        metabolite_stoichiometry = {
            metabolites.get_by_id(mid): stoichiometry
            for mid, stoichiometry in iteritems(metabolite_stoichiometry)
            if mid in metabolites}
        mass_reaction.add_metabolites(metabolite_stoichiometry)

        # Read the reaction kinetic law
        rate_eq, local_parameters = _read_reaction_kinetic_law_from_sbml(
            reaction, mass_reaction, f_replace, **kwargs)
        local_parameters_dict.update(local_parameters)
        # Check whether rate law should be a custom rate law
        mass_action_rates = [
            strip_time(mass_reaction.get_mass_action_rate(x))
            for x in range(1, 4)]
        if rate_eq in mass_action_rates:
            # Rate law is a mass action rate law identical to the one
            # automatically generated when reaction rate type is 1, 2, or 3
            for x, rate in enumerate(mass_action_rates):
                if rate_eq == rate:
                    mass_reaction._rtype = int(x + 1)
        else:
            custom_rates_dict[mass_reaction] = str(rate_eq)

        # Parse the notes dict for fall back solutions
        mass_notes = _parse_notes_dict(reaction)
        mass_reaction.gpr = _read_reaction_gpr_from_sbml(reaction, mass_notes,
                                                         f_replace)
        # Add flux bounds to reaction, including missing bounds if desired.
        mass_reaction.bounds = _read_reaction_flux_bounds_from_sbml(
            model, reaction, kwargs.get("set_missing_bounds"))

        # Add subsystem and EnzymeModuleReaction attributes from note dict
        if "SUBSYSTEM" in mass_notes:
            mass_reaction.subsystem = mass_notes["SUBSYSTEM"]
            del mass_notes["SUBSYSTEM"]

        if "enzyme_module_id" in mass_notes:
            mass_reaction = EnzymeModuleReaction(mass_reaction)
            _read_enzyme_attr_info_from_notes(
                mass_reaction, mass_notes, f_replace=f_replace)

        # Add the notes and annotations to the reaction
        mass_reaction.notes = mass_notes
        mass_reaction.annotation = _parse_annotations(reaction)

        # Add reaction to list of reaction
        reaction_list.append(mass_reaction)

    return reaction_list, local_parameters_dict, custom_rates_dict


def _read_reaction_species_from_sbml(reaction, f_replace, **kwargs):
    """Read the SBML reaction species stoichiometry, returning it as a dict.

    Warnings
    --------
    This method is intended for internal use only.

    """
    # Get number type as float or int
    number = kwargs.get("number")

    metabolite_stoichiometry = defaultdict(lambda: 0)
    for specie_type, sign in zip(["Reactants", "Products"], [-1, 1]):
        # Iterate through list
        species_ref_list = getattr(reaction, "getListOf" + specie_type)()
        for sref in species_ref_list:
            # Check specie ID
            sid = _check_required(sref, sref.getSpecies(), "species")
            sid = _get_corrected_id(sid, (f_replace, F_SPECIE),
                                    None, kwargs.get("remove_char"))
            # Add specie IDs and stoichiometry to dict
            metabolite_stoichiometry[sid] += sign * number(
                _check_required(sref, sref.getStoichiometry(),
                                "stoichiometry"))

    return metabolite_stoichiometry


def _read_reaction_kinetic_law_from_sbml(reaction, mass_reaction, f_replace,
                                         **kwargs):
    """Read the SBML reaction kinetic law and return it.

    Warnings
    --------
    This method is intended for internal use only.

    """
    mass_rid = mass_reaction.id
    sbml_species = list(
        reaction.getListOfReactants()) + list(
            reaction.getListOfProducts()) + list(reaction.getListOfModifiers())
    sbml_species = [sref.getSpecies() for sref in sbml_species]

    local_parameters = {}
    if reaction.isSetKineticLaw():
        sbml_rid = reaction.getIdAttribute()
        kinetic_law = reaction.getKineticLaw()
        # Get the kinetic law and the rate equation as a string.
        kinetic_law = reaction.getKineticLaw()
        rate_eq = _check_required(kinetic_law, kinetic_law.getFormula(),
                                  "formula")

        # Perform substitution for power law operations to sympify rate
        for match in KLAW_POW_RE.finditer(rate_eq):
            old = match.group(0)
            new = "(({0})**{1})".format(match.group("arg"), match.group("exp"))
            rate_eq = rate_eq.replace(old, new)
        # Try to sympify the reaction rate
        try:
            rate_eq = sympify(rate_eq)
        except SympifyError as e:
            raise MassSBMLError(e)
        # If ID replacements were performed earlier then apply the ID
        # replacements for metabolite and parameter arguments in rate law also.
        id_subs = {}
        for arg in list(rate_eq.atoms(Symbol)):
            arg = str(arg)
            new_arg = arg
            # Check if reaction is in the name of the parameter
            if re.search(sbml_rid, arg) and sbml_rid != mass_rid:
                new_arg = _get_corrected_id(
                    new_arg, (sbml_rid, mass_rid, arg), "Parameter",
                    kwargs.get("remove_char"))
            elif arg in sbml_species:
                new_arg = _get_corrected_id(
                    new_arg, (f_replace, F_SPECIE), None,
                    kwargs.get("remove_char"))
            else:
                if kwargs.get("remove_char"):
                    new_arg = _remove_char_from_id(new_arg)
            id_subs[arg] = new_arg
        # Make rate equation
        rate_eq = rate_eq.subs(id_subs)
        for local_parameter in kinetic_law.getListOfLocalParameters():
            pid = _check_required(
                local_parameter, local_parameter.getIdAttribute(), "id")
            value = local_parameter.getValue()
            if re.search(sbml_rid, pid) and sbml_rid != mass_rid:
                pid = _get_corrected_id(
                    pid, (sbml_rid, mass_rid, pid), "Parameter",
                    kwargs.get("remove_char"))
            elif kwargs.get("remove_char"):
                pid = _remove_char_from_id(pid)
            local_parameters[pid] = value

    else:
        LOGGER.warning(
            "No kinetic law found for SBML reaction '%s'. Therefore, assigning"
            " the MassReaction '%s' a rate law based on Mass Action Kinetics.",
            reaction, mass_rid)
        rate_eq = mass_reaction.get_mass_action_rate(1)

    return rate_eq, local_parameters


def _read_reaction_gpr_from_sbml(reaction, mass_notes, f_replace):
    """Read the GPR information from SBMLDocument and return as a string.

    Warnings
    --------
    This method is intended for internal use only.

    """
    reaction_fbc = reaction.getPlugin("fbc")
    if reaction_fbc:
        # GPR rules
        def process_association(ass):
            """Recursively convert gpr association to a gpr string."""
            if ass.isFbcOr():
                return " ".join([
                    "(", ' or '.join(
                        process_association(c)
                        for c in ass.getListOfAssociations()), ")"])
            elif ass.isFbcAnd():
                return " ".join([
                    "(", ' and '.join(
                        process_association(c)
                        for c in ass.getListOfAssociations()), ")"])
            elif ass.isGeneProductRef():
                gid = ass.getGeneProduct()
                if f_replace and F_GENE in f_replace:
                    gid = f_replace[F_GENE](gid)
                return gid
            else:
                return ""

        gpr = ''
        gpa = reaction_fbc.getGeneProductAssociation()
        if gpa is not None:
            association = gpa.getAssociation()
            gpr = process_association(association)
    else:
        # Fallback to notes information
        if "GENE ASSOCIATION" in mass_notes:
            gpr = mass_notes['GENE ASSOCIATION']
        elif "GENE_ASSOCIATION" in mass_notes:
            gpr = mass_notes['GENE_ASSOCIATION']
        else:
            gpr = ''
        if gpr:
            LOGGER.warning(
                "Use of GENE ASSOCIATION or GENE_ASSOCIATION in the notes "
                "element is discouraged, use fbc:gpr instead: %s", reaction)
            if f_replace and F_GENE in f_replace:
                gpr = " ".join(f_replace[F_GENE](t) for t in gpr.split(' '))

    # Remove outside parenthesis, if any
    if gpr.startswith("(") and gpr.endswith(")"):
        gpr = gpr[1:-1].strip()

    return gpr


def _read_reaction_flux_bounds_from_sbml(model, reaction, set_missing_bounds):
    """Read the reaction lower and upper flux bounds, returned as a tuple.

    Warnings
    --------
    This method is intended for internal use only.

    """
    p_lb, p_ub = (None, None)
    reaction_fbc = reaction.getPlugin("fbc")
    if reaction_fbc:
        # bounds in fbc
        lb_id = reaction_fbc.getLowerFluxBound()
        if lb_id:
            p_lb = model.getParameter(lb_id)  # type: libsbml.Parameter
            if p_lb and p_lb.getConstant() and p_lb.getValue() is not None:
                lower_bound = p_lb.getValue()

        ub_id = reaction_fbc.getUpperFluxBound()
        if ub_id:
            p_ub = model.getParameter(ub_id)  # type: libsbml.Parameter
            if p_ub and p_ub.getConstant() and p_ub.getValue() is not None:
                upper_bound = p_ub.getValue()

    if p_lb is None:
        if set_missing_bounds:
            lower_bound = MASSCONFIGURATION.lower_bound
        else:
            lower_bound = -float("inf")

    if p_ub is None:
        if set_missing_bounds:
            upper_bound = MASSCONFIGURATION.upper_bound
        else:
            upper_bound = float("inf")

    return (lower_bound, upper_bound)


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


def _read_enzyme_modules_from_sbml(enzyme_group_dict, f_replace, **kwargs):
    """Read the SBML groups for enzymes, returning them as a list.

    EnzymeModuleDicts are created and returned in a list.

    Warnings
    --------
    This method is intended for internal use only.

    """
    def _get_group_members(group, gid=None, **kwargs):
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
            obj_id = _get_corrected_id(obj_id, (f_replace, replace_key),
                                       None, kwargs.get("remove_char"))
            # Add object to member list
            member_list.append(obj_id)

        return member_list

    def _parse_enzyme_group(gid, gname, group, enzyme_module_dict, **kwargs):
        """Parse the various enzyme groups and add to the EnzymeModuleDict."""
        if "_categorized" in gid and CATEGORY_GROUP_RE.search(gid):
            attr = gid[:CATEGORY_GROUP_RE.search(gid).start()]
            gid = attr[:-len("_categorized")]
            member_list = _get_group_members(group, gid, **kwargs)
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
            member_list = _get_group_members(group, gid, **kwargs)
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
                                **kwargs)

        for attr in ["ligands", "forms", "reactions"]:
            attr = "enzyme_module_" + attr + "_categorized"
            enzyme_module_dict[attr] = categorized_dicts[attr]

        enzyme_module_dicts.append(enzyme_module_dict)

    return enzyme_module_dicts


def _read_enzyme_attr_info_from_notes(mass_obj, notes, f_replace, **kwargs):
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
                mid = _get_corrected_id(mid, (f_replace, F_SPECIE), None,
                                        kwargs.get("remove_char"))
                attr_value[mid] = value
            setattr(mass_obj, "_" + enz_attr, attr_value)
        else:
            setattr(mass_obj, enz_attr, notes.get(enz_attr))
        # Remove newly set attribute from the notes dict
        del notes[enz_attr]


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


# -----------------------------------------------------------------------------
# Write SBML
# -----------------------------------------------------------------------------
def write_sbml_model(mass_model, filename, f_replace=None, **kwargs):
    """Write MassModel to filename in SBML format.

    The created model is SBML level 3 version 1 core (L3V1) using packages
    fbc-v2 and groups-v1 for optimal exporting. Not including these packages
    may result in some information loss when exporting the model.

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
    use_fbc_package: bool
        Should the fbc package be used. Defaukt is True.
    use_groups_package: bool
        Should the groups package be used. Defaukt is True.
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
    use_fbc_package: bool
        Should the fbc package be used. Defaukt is True.
    use_groups_package: bool
        Should the groups package be used. Defaukt is True.
    units: bool
        Should the units be written into the SBMLDocument.
        Default is True.
    local_parameters: bool
        If True, standard reaction parameters [kf, Keq, kr] be written
        as SBML local parameters of the kinetic law. Otherwise. all parameters
        are written as SBML global model parameters. Default is True.

    Returns
    -------
    doc: libsbml.SBMLDocument

    Warnings
    --------
    This method is intended for internal use only.

    """
    # Create dict of default kwargs
    default_kwargs = {
        "use_fbc_package": True,
        "use_groups_package": True,
        "units": True,
        "local_parameters": True}
    # Check kwargs
    kwargs = _check_kwargs(default_kwargs, kwargs)

    if f_replace is None:
        f_replace = F_REPLACE

    # Get the SBMLDocument and SBML model objects
    doc, (model, model_fbc, model_groups) = _create_doc_and_models(mass_model,
                                                                   **kwargs)
    # Write the model Meta Information (ModelHistory) into the SBMLDocument
    _write_model_meta_info_to_sbml(doc, mass_model, model)

    # Write model notes and annotations into the SBMLDocument
    model_notes = mass_model.notes
    if mass_model.description:
        model_notes["description"] = mass_model.description
    _sbase_notes_dict(model, model_notes)
    _sbase_annotations(model, mass_model.annotation)

    # Write the units information into the SBMLDocument
    if kwargs.get("units") and mass_model.units:
        _write_model_units_to_sbml(model, mass_model, **kwargs)

    # Write the compartment information into the SBMLDocument
    _write_model_compartments_to_sbml(model, mass_model)

    # Write the species information into the SBMLDocument
    # Includes MassMetabolites, EnzymeModuleForms, and their concentrations
    _write_model_species_to_sbml(model, mass_model, f_replace, **kwargs)
    if mass_model.boundary_conditions:
        _write_model_boundary_species_to_sbml(model, mass_model, f_replace)

    if kwargs.get("use_fbc_package"):
        # Write the flux bound parameters into the SBMLDocument
        _write_global_flux_bounds_to_sbml(model, mass_model, **kwargs)
        # Write the gene information into the SBMLDocument
        _write_model_genes_to_sbml(model_fbc, mass_model, f_replace)

    # Write the reaction information into the SBMLDocument
    # Includes MassReactions, EnzymeModuleReactions, kinetic laws,
    # GPR, reaction flux bounds, and kinetic parameters as local parameters
    _write_model_reactions_to_sbml(model, mass_model, f_replace, **kwargs)
    # Write kinetic parameters into the SBMLDocument as global parameters
    # if they are not already written as kinetic law local parameters
    _write_model_kinetic_parameters_to_sbml(model, mass_model, f_replace,
                                            **kwargs)
    # Write the EnzymeModule specific information into the SBMLDocument
    if kwargs.get("use_groups_package"):
        if isinstance(mass_model, EnzymeModule):
            _write_enzyme_modules_to_sbml(model_groups, mass_model, f_replace)
        # Write the EnzymeModuleDict specific information into the SBMLDocument
        if mass_model.enzyme_modules:
            for enzyme in mass_model.enzyme_modules:
                _write_enzyme_modules_to_sbml(model_groups, enzyme, f_replace)

    return doc


def _create_doc_and_models(mass_model, **kwargs):
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
    # Set model ID, meta ID, and name in SBMLDocument
    if mass_model.id is not None:
        _check(model.setIdAttribute(mass_model.id), "set model id")
        _check(model.setMetaId("meta_" + mass_model.id), "set model meta id")
    else:
        _check(model.setMetaId("meta_model"), "set model meta id")
    if mass_model.name is not None:
        _check(model.setName(mass_model.name), "set model name")

    # Function to enable package plugins.
    def _enable_package_plugin(plugin, version_info):
        """Enable the plugin for the package extension."""
        # Enable package extension
        uri_fmt = "http://www.sbml.org/sbml/level{0}/version{1}/{2}/version{3}"
        _check(
            doc.enablePackage(
                uri_fmt.format(*SBML_LEVEL_VERSION, plugin, version_info),
                plugin, True),
            "enable package extension {0}-v{1}".format(plugin, version_info))
        # Set package required to false
        # Set required to false
        _check(doc.setPackageRequired(plugin, False),
               "set {0} extension required to false".format(plugin))
        # Get plugin model and return
        model_plugin = model.getPlugin(plugin)
        return model_plugin

    # Get plugin models if package extensions are used.
    model_fbc, model_groups = (None, None)

    if kwargs.get("use_fbc_package"):
        # Add fbc package extension to SBMLDocument
        model_fbc = _enable_package_plugin("fbc", FBC_VERSION)
        # Set strictness of fbc package
        _check(model_fbc.setStrict(False),
               "set fbc plugin strictness to False")

    if kwargs.get("use_groups_package"):
        # Add groups package extension to SBMLDocument
        model_groups = _enable_package_plugin("groups", GROUPS_VERSION)

    return doc, (model, model_fbc, model_groups)


def _write_model_meta_info_to_sbml(doc, mass_model, model):
    """Write MassModel Meta Information (ModelHistory) into SBMLDocument.

    Warnings
    --------
    This method is intended for internal use only.

    """
    # Get the original meta information if it exists
    meta = mass_model._sbml if hasattr(mass_model, "_sbml") else {}

    # Check if there are creators
    creators = meta["creators"] if meta.get("creators") else []

    mass_creator = MASSCONFIGURATION.model_creator
    # Add creator in MassConfiguration if it doesn't already exist.
    if not any([mass_creator == c for c in creators]) and \
       all([v != "" for v in itervalues(mass_creator)]):
        creators += [mass_creator]

    # Make model history if creators are present (required for history)
    if creators:
        # Make ModelHistory object
        history = libsbml.ModelHistory()
        _check(history, "create ModelHistory")

        if meta.get("created", False):
            # Use the previous creation date
            timestr = meta["created"]
        else:
            timestr = datetime.datetime.now()
            timestr = timestr.strftime('%Y-%m-%dT%H:%M:%S')

        # Set the creation date
        date = libsbml.Date(timestr)
        _check(date, "create the created date")
        _check(history.setCreatedDate(date), "set created date")
        # Set the modified date
        if meta.get("modified", False):
            for k, timestr in enumerate(meta["modified"]):
                date = libsbml.Date(timestr)
                _check(history.setModifiedDate(date),
                       "set modified date " + str(k))
        else:
            date = libsbml.Date(timestr)
            _check(history.setModifiedDate(date), "set modified date")

        # Make Creator objects and add to ModelHistory
        for mass_creator in creators:
            # Make creator
            creator = libsbml.ModelCreator()
            _check(creator, "create model creator")
            # Add family name, given name, organization, and email attributes
            for k in ["familyName", "givenName", "organization", "email"]:
                if mass_creator.get(k, False):
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

    # Set meta annotation and notes in doc
    if "annotation" in meta:
        _sbase_annotations(doc, meta["annotation"])
    if "notes" in meta:
        _sbase_notes_dict(doc, meta["notes"])


def _write_model_units_to_sbml(model, mass_model, **kwargs):
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
                   "set Unit kind '" + u.kind + "'" + _for_id(id_str))
            _check(unit.setExponent(u.exponent),
                   "set Unit exponent" + _for_id(u.kind))
            _check(unit.setScale(u.scale), "set Unit scale" + _for_id(u.kind))
            _check(unit.setMultiplier(u.multiplier),
                   "set Unit multiplier" + _for_id(u.kind))

    # Get units in MassModel to write
    to_write = list(mass_model.units)
    # Add COBRA unit if it does not exist in model units
    if COBRA_FLUX_UNIT.id not in mass_model.units \
       and kwargs.get("use_fbc_package"):
        to_write += [COBRA_FLUX_UNIT]

    # Write flux bound units
    for unit_def in to_write:
        _write_unit_definition(
            id_str=unit_def.id, name=unit_def.name, unit_def=unit_def)


def _write_model_compartments_to_sbml(model, mass_model):
    """Write MassModel compartment information into SBMLDocument.

    Warnings
    --------
    This method is intended for internal use only.

    """
    # TODO Revisit based on new cobra compartment implementation
    if not mass_model.compartments:
        default_compartment_dict = MASSCONFIGURATION.default_compartment
        LOGGER.warning(
            "No compartments found in model. Therefore creating compartment "
            "'%s' for entire model.", list(default_compartment_dict)[0])
        compartment_dict = default_compartment_dict
    else:
        compartment_dict = mass_model.compartments

    if mass_model.boundary_conditions:
        compartment_dict.update(MASSCONFIGURATION.boundary_compartment)

    for cid, name in iteritems(compartment_dict):
        compartment = model.createCompartment()
        _check(compartment, "create model compartment" + _for_id(cid))
        _check(compartment.setIdAttribute(cid),
               "set compartment id" + _for_id(cid))
        _check(compartment.setName(name),
               "set compartment name" + _for_id(cid))
        _check(compartment.setSpatialDimensions(3),
               "set compartment spatial dimensions" + _for_id(cid))
        _check(compartment.setConstant(True),
               "set compartment constant" + _for_id(cid))
        _check(compartment.setSize(1), "set compartment size" + _for_id(cid))


def _write_global_flux_bounds_to_sbml(model, mass_model, **kwargs):
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
        min_value, max_value = MASSCONFIGURATION.bounds

    # Get flux units if units are desired.
    flux_udef = model.getUnitDefinition(COBRA_FLUX_UNIT.id)

    # Create the model flux parameters
    _create_parameter(
        model, pid=LOWER_BOUND_ID, value=min_value, sbo=SBO_DEFAULT_FLUX_BOUND,
        udef=flux_udef, **kwargs)
    _create_parameter(
        model, pid=UPPER_BOUND_ID, value=max_value, sbo=SBO_DEFAULT_FLUX_BOUND,
        udef=flux_udef, **kwargs)
    _create_parameter(
        model, pid=ZERO_BOUND_ID, value=0, sbo=SBO_DEFAULT_FLUX_BOUND,
        udef=flux_udef, **kwargs)
    _create_parameter(
        model, pid=BOUND_MINUS_INF, value=-float("Inf"), sbo=SBO_FLUX_BOUND,
        udef=flux_udef, **kwargs)
    _create_parameter(
        model, pid=BOUND_PLUS_INF, value=float("Inf"), sbo=SBO_FLUX_BOUND,
        udef=flux_udef, **kwargs)


def _write_model_kinetic_parameters_to_sbml(model, mass_model, f_replace,
                                            **kwargs):
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
    local_parameters = kwargs.get("local_parameters")
    for parameter_type, parameter_dict in iteritems(mass_model.parameters):
        # Skip over parameters already written into the model.
        if (local_parameters and parameter_type != "v") or \
           (not local_parameters and parameter_type == "Boundary"):
            continue

        if parameter_type == "v":
            constant = False
        else:
            constant = True
        # Create kf, Keq, kr, and custom parameters, handling ID corrections
        # for recognized reaction IDs in the parameters
        for pid, value in iteritems(parameter_dict):
            try:
                p_type, rid = pid.split("_", 1)
                rid = getattr(mass_model.reactions.get_by_id(rid), "id")
            except (ValueError, KeyError):
                pass
            else:
                if f_replace and F_REACTION_REV in f_replace:
                    rid = f_replace[F_REACTION_REV](rid)
                pid = p_type + "_" + rid
            _create_parameter(model, pid, value, constant=constant, sbo=None,
                              udef=None, **kwargs)


def _write_model_species_to_sbml(model, mass_model, f_replace, **kwargs):
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

        # Use the default compartment if no species have set compartments
        if metabolite.compartment is None:
            cid = next(iter(MASSCONFIGURATION.default_compartment))
            is_boundary_condition = False
        else:
            cid = metabolite.compartment
            is_boundary_condition = bool(
                cid == next(iter(MASSCONFIGURATION.boundary_compartment)))
        _check(specie.setCompartment(cid),
               "set specie compartment" + _for_id(mid))

        # Set species constant and whether it is a boundary condition
        _check(specie.setConstant(metabolite.fixed),
               "set specie constant" + _for_id(mid))
        _check(specie.setBoundaryCondition(is_boundary_condition),
               "set specie boundary condition" + _for_id(mid))

        if MASSCONFIGURATION.include_compartments_in_rates:
            has_only_substance_units = True
            setter = specie.setInitialAmount
        else:
            has_only_substance_units = False
            setter = specie.setInitialConcentration

        # Set metabolite initial condition value and unit restrictions
        _check(specie.setHasOnlySubstanceUnits(has_only_substance_units),
               "set specie unit substance restrictions" + _for_id(mid))
        if metabolite.initial_condition is not None:
            _check(setter(metabolite.initial_condition),
                   "set specie initial condition" + _for_id(mid))
        # Set specie formula and charge
        additional_notes = _write_specie_formula_charge_into_sbml(specie,
                                                                  metabolite,
                                                                  **kwargs)
        # Set metabolite annotation and notes
        if isinstance(metabolite, EnzymeModuleForm):
            # Set enzyme information if EnzymeModuleForm
            _write_enzyme_attr_info_to_notes(specie, metabolite, f_replace,
                                             additional_notes=additional_notes)
        else:
            # Otherwise just set the MassMetabolite notes normally
            additional_notes.update(metabolite.notes)
            _sbase_notes_dict(specie, additional_notes)
        # Set metabolite annotation
        _sbase_annotations(specie, metabolite.annotation)


def _write_specie_formula_charge_into_sbml(specie, metabolite, **kwargs):
    """Write the specie charge and formula into the SBMLDocument.

     Warnings
    --------
    This method is intended for internal use only.

    """
    notes = {}
    if kwargs.get("use_fbc_package"):
        # Set metabolite charge via fbc package
        specie_fbc = specie.getPlugin("fbc")
        mid = specie.getIdAttribute()
        _check(specie_fbc, "get fbc plugin for specie" + _for_id(mid))
        if metabolite.charge is not None:
            _check(specie_fbc.setCharge(metabolite.charge),
                   "set specie charge" + _for_id(mid))
        # Set metabolite formula via fbc package
        if metabolite.formula is not None:
            sbml_formula = _f_moiety_formula_rev(metabolite)
            _check(specie_fbc.setChemicalFormula(sbml_formula),
                   "set specie formula" + _for_id(mid))
    else:
        # Fallback to notes information
        if metabolite.charge is not None:
            notes["CHARGE"] = metabolite.charge
        if metabolite.formula is not None:
            notes["FORMULA"] = metabolite.formula

    return notes


def _write_model_boundary_species_to_sbml(model, mass_model, f_replace):
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
        _check(
            specie.setCompartment(
                next(iter(MASSCONFIGURATION.boundary_compartment))),
            "set boundary specie compartment" + _for_id(bmid))
        # Set boundary specie value
        if isinstance(bc_value, (integer_types, float)):
            constant = True
            if MASSCONFIGURATION.include_compartments_in_rates:
                has_only_substance_units = True
                setter = specie.setInitialAmount
            else:
                has_only_substance_units = False
                setter = specie.setInitialConcentration
            _check(setter(bc_value),
                   "set boundary specie concentration" + _for_id(bmid))
        else:
            constant = False
            assignment_rule = model.createAssignmentRule()
            _check(assignment_rule,
                   "create model assignment rule" + _for_id(bmid))
            _check(assignment_rule.setVariable(bmid),
                   "set assignment rule variable" + _for_id(bmid))

            # Create MathML string from the sympy expression
            math = _create_math_xml_str_from_sympy_expr(bc_value)
            # Set the math for the AssignmentRule
            _check(
                assignment_rule.setMath(libsbml.readMathMLFromString(math)),
                "set math on assignment rule" + _for_id(bmid))

        # Set specie constant, unit restrictions, and boundary condition
        _check(specie.setConstant(constant),
               "set specie constant" + _for_id(bmid))
        _check(specie.setHasOnlySubstanceUnits(has_only_substance_units),
               "set specie unit substance restrictions" + _for_id(bmid))
        _check(specie.setBoundaryCondition(True),
               "set specie boundary condition" + _for_id(bmid))


def _write_model_genes_to_sbml(model_fbc, mass_model, f_replace):
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


def _write_model_reactions_to_sbml(model, mass_model, f_replace, **kwargs):
    """Write MassModel reaction information into SBMLDocument.

    Warnings
    --------
    This method is intended for internal use only.

    """
    additional_notes = {}
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
        _check(reaction.setFast(False), "set reaction fast" + _for_id(rid))

        # Write specie references into SBML reaction
        _write_reaction_specie_ref_to_sbml(reaction, mass_reaction, f_replace)

        # Write reaction kinetic law into SBML reaction
        _write_reaction_kinetic_law_to_sbml(reaction, mass_reaction, f_replace,
                                            **kwargs)

        # Write bouns and GPR if FBC package enabled.
        if kwargs.get("use_fbc_package"):
            # Get plugin for reaction
            reaction_fbc = reaction.getPlugin("fbc")
            _check(reaction_fbc, "get fbc plugin for reaction" + _for_id(rid))

            # Set reaction flux bounds
            flux_udef = model.getUnitDefinition(COBRA_FLUX_UNIT.id)
            reaction_fbc.setLowerFluxBound(
                _create_bound(
                    model, reaction=mass_reaction, bound_type="lower_bound",
                    f_replace=f_replace, units=kwargs.get("units"),
                    flux_udef=flux_udef))
            reaction_fbc.setUpperFluxBound(
                _create_bound(
                    model, reaction=mass_reaction, bound_type="upper_bound",
                    f_replace=f_replace, units=kwargs.get("units"),
                    flux_udef=flux_udef))

            # Set the reaction GPR
            _write_reaction_gpr_to_sbml(reaction_fbc, mass_reaction, f_replace)

        # Write the subsystem attribute into the notes.
        if mass_reaction.subsystem:
            additional_notes["SUBSYSTEM"] = mass_reaction.subsystem
        if isinstance(mass_reaction, EnzymeModuleReaction):
            # Set enzyme information if EnzymeModuleReaction
            _write_enzyme_attr_info_to_notes(
                reaction, mass_reaction, f_replace,
                additional_notes=additional_notes)
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
    rid = reaction.getIdAttribute()
    metabolites = mass_reaction.metabolites
    boundary_conditions = mass_reaction.model.boundary_conditions
    # Write boundary metabolite in addition if boundary reaction
    if mass_reaction.boundary and not mass_reaction.products:
        metabolites.update({mass_reaction.boundary_metabolite: 1})
    if mass_reaction.boundary and not mass_reaction.reactants:
        metabolites.update({mass_reaction.boundary_metabolite: -1})

    for metabolite, stoichiometry in iteritems(metabolites):
        # Get specie id
        sid = str(metabolite)
        if f_replace and F_SPECIE_REV in f_replace:
            sid = f_replace[F_SPECIE_REV](sid)
        # Comparison to determine whether metabolite is reactant or product
        if stoichiometry < 0:
            sref = reaction.createReactant()
        else:
            sref = reaction.createProduct()
        _check(sref, "create specie reference " + sid + _for_id(rid))
        # Set ID and stoichiometry
        _check(sref.setSpecies(sid), "set specie reference id" + _for_id(sid))
        _check(sref.setStoichiometry(abs(stoichiometry)),
               "set specie reference stoichiometry" + _for_id(sid))

        if boundary_conditions.get(metabolite) and\
           not isinstance(boundary_conditions.get(metabolite), float):
            constant = False
        else:
            constant = True
        _check(sref.setConstant(constant),
               "set specie reference constant" + _for_id(sid))

        # Fixed metabolites are treated as modifier species in rate laws and
        # therefore cannot exist in the SBML reaction reactants or products.
        if isinstance(metabolite, MassMetabolite) and metabolite.fixed:
            msref = reaction.createModifier()
            _check(msref,
                   "create modifier specie " + sid + _for_id(rid))
            _check(msref.setSpecies(sid),
                   "set modifier specie species " + sid + _for_id(rid))


def _write_reaction_kinetic_law_to_sbml(reaction, mass_reaction, f_replace,
                                        **kwargs):
    """Write MassReaction kinetic law information into SBMLDocument.

    Warnings
    --------
    This method is intended for internal use only.

    """
    mass_rid = mass_reaction.id
    rid = reaction.getIdAttribute()
    sbml_species = list(
        reaction.getListOfReactants()) + list(
            reaction.getListOfProducts())
    sbml_species = [sref.getSpecies() for sref in sbml_species]

    rate_eq = strip_time(mass_reaction.rate)
    # Create kinetic law as AST math and write into the SBMLDocument
    kinetic_law = reaction.createKineticLaw()
    _check(kinetic_law, "create kinetic law" + _for_id(rid))
    if kwargs.get("local_parameters"):
        all_parameter_values = mass_reaction.model._get_all_parameters()

    def _create_local_parameter(key, pid, sbo, udef, **kwargs):
        """Create a local reaction paramter and write it into SBMLDocument."""
        try:
            value = all_parameter_values[key]
        except KeyError:
            pass
        else:
            _create_parameter(kinetic_law, pid, value,
                              sbo=sbo, udef=udef, **kwargs)

    # If ID replacements were performed earlier then apply the ID
    # replacements for metabolite and parameter arguments in rate law also.
    id_subs = {}
    for arg in list(rate_eq.atoms(Symbol)):
        arg = str(arg)
        new_arg = arg
        # Check if reaction is in the name of the parameter
        if re.search(mass_rid, arg):
            new_arg = re.sub(mass_rid, rid, new_arg)
            if kwargs.get("local_parameters"):
                _create_local_parameter(
                    arg, new_arg, sbo=None, udef=None, **kwargs)
            id_subs[arg] = new_arg
            continue

        try:
            if arg not in mass_reaction.model.boundary_conditions:
                new_arg = str(
                    mass_reaction.model.metabolites.get_by_id(arg))
        except KeyError:
            if arg in mass_reaction.model.custom_parameters:
                if kwargs.get("local_parameters"):
                    _create_local_parameter(
                        arg, new_arg, sbo=None, udef=None, **kwargs)
            else:
                # Account for metabolites in the rate law that
                # are not considered reactants or products via modifiers
                if f_replace and F_SPECIE_REV in f_replace:
                    new_arg = f_replace[F_SPECIE_REV](new_arg)
                msref = reaction.createModifier()
                _check(msref,
                       "create modifier specie " + new_arg + _for_id(rid))
                _check(msref.setSpecies(new_arg),
                       "set modifier specie species " + new_arg + _for_id(rid))
        else:
            if f_replace and F_SPECIE_REV in f_replace:
                new_arg = f_replace[F_SPECIE_REV](new_arg)
        id_subs[arg] = new_arg

    # Make xml string of rate equation via sympy conversion to MathML
    rate_eq = rate_eq.subs(id_subs)
    _check(
        kinetic_law.setMath(
            libsbml.readMathMLFromString(
                _create_math_xml_str_from_sympy_expr(rate_eq))),
        "set math on kinetic law" + _for_id(rid))


def _write_reaction_gpr_to_sbml(reaction_fbc, mass_reaction, f_replace):
    """Write MassReaction GPR information into SBMLDocument.

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


def _create_parameter(sbml_obj, pid, value, constant=True, sbo=None, udef=None,
                      **kwargs):
    """Create a parameter to be written into the SBMLDocument.

    Warnings
    --------
    This method is intended for internal use only.

    """
    local_parameters = bool(kwargs.get("local_parameters") and
                            sbml_obj.__class__.__name__ == "KineticLaw")
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

    # Set SBO term and units if desired.
    if sbo:
        _check(parameter.setSBOTerm(sbo),
               "set" + p_type + "parameter sbo term" + _for_id(pid))
    if kwargs.get("units") and udef is not None:
        _check(parameter.setUnits(udef.getIdAttribute()),
               "set" + p_type + "parameter units" + _for_id(pid))


def _write_enzyme_modules_to_sbml(model_groups, mass_obj, f_replace):
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


def _write_enzyme_attr_info_to_notes(sbml_obj, mass_obj, f_replace,
                                     additional_notes=None):
    """Write the enzyme object attribute into SBMLDocument in the notes.

    Applies to the EnzymeModuleForms and EnzymeModuleReaction

    Warnings
    --------
    This method is intended for internal use only.

    """
    # Get original notes from mass object.
    notes = mass_obj.notes.copy()
    if additional_notes:
        notes.update(additional_notes)
    # Get the subclass specific attributes to write into the notes.
    attributes = get_subclass_specific_attributes(mass_obj)
    default_values = [{} if "bound" in attr_name else ""
                      for attr_name in attributes]

    for attr_name, default_value in zip(attributes, default_values):
        attr_value = getattr(mass_obj, attr_name)
        # Skip attributes equal to their defaults
        if attr_value == default_value and attr_name != "enzyme_module_form":
            continue

        if "bound" in attr_name:
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
    model: :class:`~mass.core.massmodel.MassModel` object
        The mass model if the file could be read successfully or None
        otherwise.
    errors: dict
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
    # libsbml validation
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

    # masspy validation (check that SBML can be read into model)
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

    # Remove stream handler
    LOGGER.removeHandler(stream_handler)
    LOGGER.propagate = True

    # Additional model tests
    if check_model:
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


def validate_sbml_model_export(mass_model, filename, f_replace=None, **kwargs):
    """Validate export of model to SBML, returning a list of errors.

    If no SBML errors or MASS fatal errors occur, the model will be written to
    the 'filename'.

    Parameters:
    -----------
    mass_model : mass.core.MassModel
        MassModel instance which is written to SBML
    filename: path to SBML file, or SBML string, or SBML file handle
        SBML which is read into a MassModel.
    f_replace: dict of replacement functions for id replacement
        A dict of replacement functions to apply on identifiers.
        If no replacements should be performed, set f_replace={}.

    Returns
    -------
    (success, errors)
    success: bool
        A bool indicating whether the model was successfully exported to
        'filename'.
    errors: dict
        Warnings and errors grouped by their respective types.

    Raises
    ------
    MassSBMLError

    """
    all_kwargs = {
        "use_fbc_package": True,
        "use_groups_package": True,
        "units": True,
        "local_parameters": True,
        "number": float,
        "set_missing_bounds": True,
        "remove_char": True,
        "stop_on_conversion_fail": True,
        "check_model": True,
        "internal_consistency": True,
        "check_units_consistency": False,
        "check_modeling_practice": False}

    # Check kwargs
    kwargs = _check_kwargs(all_kwargs, kwargs)
    errors = {
        "SBML_FATAL": [],
        "SBML_ERROR": [],
        "SBML_SCHEMA_ERROR": [],
        "SBML_WARNING": [],
        "MASS_FATAL": [],
        "MASS_ERROR": [],
        "MASS_WARNING": [],
        "MASS_CHECK": [],
    }
    all_kwargs = {"export": {}, "validate": {}}
    for k, v in iteritems(kwargs):
        if k in ["use_fbc_package", "use_groups_package", "units",
                 "local_parameters"]:
            all_kwargs["export"][k] = v
        else:
            all_kwargs["validate"][k] = v
    try:
        doc = _model_to_sbml(mass_model, f_replace=f_replace,
                             **all_kwargs["export"])
        sbml_str = libsbml.writeSBMLToString(doc)
        model, errors = validate_sbml_model(sbml_str, **all_kwargs["validate"])

    except Exception as e:
        success = False
        model = None
        errors["MASS_FATAL"].append(str(e))
        return success, errors

    if model is not None:
        success = True
        for key in ["SBML_FATAL", "SBML_ERROR", "SBML_SCHEMA_ERROR",
                    "MASS_FATAL"]:
            if errors[key]:
                success = False
    else:
        success = False

    if success and isinstance(filename, string_types):
        # Write to path
        libsbml.writeSBMLToFile(doc, filename)
    elif success and hasattr(filename, "write"):
        # Write to file handle
        filename.write(sbml_str)

    return success, errors
