# -*- coding: utf-8 -*-
"""TODO Module Docstrings."""
from __future__ import absolute_import

import warnings
from collections import Counter
from math import ceil, floor

from mass import config as _config
from mass.util.util import ensure_iterable

from six import iteritems, iterkeys

import sympy as sym

from tabulate import tabulate

# Global
_T_SYM = sym.Symbol("t")
_ZERO_TOL = _config.ZERO_TOLERANCE


def qcqa_model(model, parameters=False, concentrations=False, fluxes=False,
               superfluous=False, elemental=False, thermodynamic=False,
               tol=None):
    """Check the model quality and print a summary of the results.

    Checking the model quality involves running a series of quality control and
    assessment tests to determine consistency (e.g. elemental, thermodynamic)
    in the model, missing values, and whether the model can be simulated.

    Parameters
    ----------
    model: mass.MassModel
        The MassModel to inspect.
    parameters: bool, optional
        If True, then check for undefined parameters in the model.
    concentrations: bool, optional
        If True, then check for undefined initial conditions and fixed
        concentrations in the model.
    fluxes: bool, optional
        If True, then check for undefined steady state fluxes in the model.
    superfluous: bool, optional
        If True, then check for superfluous parameters in the model and ensure
        existing parameters are consistent with one another if superfluous
        parameters are present.
    elemental: bool, optional
        If True, then check for elemental consistency in the model.
        Exchange reactions are ignored.
    thermodynamic: bool, optional
        If True, then check for thermodynamic consistency in the model.
        Exchange reactions are ignored.
    tol: float, optional
        The tolerance used in consistency checking. If None provided, the
        global zero tolerance is used.

    """
    tol = _set_tolerance(tol)
    # Set up empty lists for storing QC/QA report items.
    table_items = [[], [], []]
    # Get missing parameters
    if True in [parameters, fluxes]:
        results = _mk_parameter_content(model, tol, parameters, fluxes)
        for to_add, item_list in zip(results, table_items):
            item_list.extend(to_add)

    # Get missing initial and fixed concentrations
    if concentrations:
        results = _mk_concentration_content(model)
        for to_add, item_list in zip(results, table_items):
            item_list.extend(to_add)

    # Check for the desired consistencies in the values.
    if True in [superfluous, elemental, thermodynamic]:
        results = _mk_consistency_content(model, tol, superfluous, elemental,
                                          thermodynamic)
        for to_add, item_list in zip(results, table_items):
            item_list.extend(to_add)

    # Check if simulatable
    checks = is_simulatable(model)
    report = _format_table_for_print(table_items, checks, model.id)
    print(report)


def qcqa_simulation(simulation, model, parameters=False, concentrations=False,
                    superfluous=False, thermodynamic=False, tol=None):
    """Check the Simulation quality and print a summary of the results.

    Checking the Simulation quality involves running a series of quality
    control and assessment tests to determine consistency (e.g. numerical,
    thermodynamic) in the stored values, missing values, and whether the
    Simulation object contains what is needed in order to simulate the
    given model.

    Parameters
    ----------
    simulation: mass.Simulation
        The Simulation to inspect.
    model: mass.MassModel
        The MassModel to inspect.
    parameters: bool, optional
        If True, then check for undefined parameters in the model.
    concentrations: bool, optional
        If True, then check for undefined initial conditions and fixed
        concentrations in the model.
    superfluous: bool, optional
        If True, then check for superfluous parameters in the model and ensure
        existing parameters are consistent with one another if superfluous
        parameters are present.
    thermodynamic: bool, optional
        If True, then check for thermodynamic consistency in the model.
        Exchange reactions are ignored.
    tol: float, optional
        The tolerance used in consistency checking. If None provided, the
        global zero tolerance is used.

    """
    tol = _set_tolerance(tol)
    # Set up empty lists for storing QC/QA report items.
    table_items = [[], [], []]
    # Get missing parameters
    if parameters:
        results = _mk_parameter_content(model, tol, parameters, False,
                                        simulation=simulation)
        for to_add, item_list in zip(results, table_items):
            item_list.extend(to_add)

    # Get missing initial and fixed concentrations
    if concentrations:
        results = _mk_concentration_content(model, simulation=simulation)
        for to_add, item_list in zip(results, table_items):
            item_list.extend(to_add)

    # Check for the desired consistencies in the values.
    if True in [superfluous, thermodynamic]:
        results = _mk_consistency_content(model, tol, superfluous, False,
                                          thermodynamic, simulation=simulation)
        for to_add, item_list in zip(results, table_items):
            item_list.extend(to_add)
    # Check if simulatable
    checks = is_simulatable(model, simulation=simulation)
    report = _format_table_for_print(table_items, checks, model.id,
                                     sim_id=simulation.id)
    print(report)


def get_missing_reaction_parameters(model, simulation=None,
                                    reaction_list=None):
    """Identify the missing parameters for reactions in a model.

    Parameters
    ----------
    model: mass.MassModel
        The MassModel to inspect.
    reaction_list: list of mass.MassReaction, optional
        A list of mass.MassReaction objects in the model to be checked.
        If None provided, will use all reactions in the model.

    Returns
    -------
    missing: dict
        A dictionary with MassReactions objects as keys and a string
        identifying the missing parameters as values. Returns as
        empty if there are no missing values.

    Note
    ----
    Will include standard reaction parameters in custom rate laws. To get
        missing custom parameters for reactions with custom rate expressions,
        use qcqa.get_missing_custom_parameters instead.

    See Also
    --------
    qcqa_model

    """
    if reaction_list is None:
        reaction_list = model.reactions
    reaction_list = ensure_iterable(reaction_list)

    if simulation is not None:
        count = _mk_simulation_param_counter(model, simulation, reaction_list)
    else:
        count = []

    missing = {}
    for rxn in reaction_list:
        missing_params = []
        parameter_keys = [rxn.Keq_str, rxn.kf_str, rxn.kr_str]
        for key in parameter_keys:
            try:
                rxn.parameters[key]
            except KeyError:
                if not rxn.reversible and key is rxn.kr_str:
                    pass
                else:
                    missing_params.append(key)
        if not count:
            missing_params = "; ".join([k.split("_")[0]
                                        for k in missing_params]).rstrip("; ")
        else:
            missing_params = "; ".join([k.split("_")[0]
                                        for k in missing_params
                                        if count[rxn.id] < 2]).rstrip("; ")
        if missing_params:
            missing[rxn] = "{0}".format(missing_params.rstrip("; "))

    return missing


def get_missing_custom_parameters(model, simulation=None, reaction_list=None):
    """Identify the missing custom parameters in a model.

    Parameters
    ----------
    model: mass.MassModel
        The MassModel to inspect.
    reaction_list : list of mass.MassReaction, optional
        A list of mass.MassReaction objects in the model to be checked.
        If None provided, will use all reactions in the model.

    Returns
    -------
    missing: dict
        A dictionary with MassReactions objects as keys and a string
        identifying the missing custom parameters as values. Returns as
        empty if there are no missing values.

    Note
    ----
    Will not include standard reaction parameters. To get missing standard
        reaction parameters for reactions with custom rate laws, use
        qcqa.get_missing_reaction_parameters instead.

    See Also
    --------
    qcqa_model

    """
    if reaction_list is None:
        reaction_list = model.reactions
    reaction_list = ensure_iterable(reaction_list)

    if simulation is not None:
        existing_customs = simulation.view_parameter_values(model)[model.id]
    else:
        existing_customs = []
    missing = {}
    # Filter out reactions without custom rates
    reaction_list = [reaction for reaction in reaction_list
                     if reaction in model.custom_rates]

    for rxn in reaction_list:
        rate = model.custom_rates[rxn]
        symbols = [str(symbol) for symbol in list(rate.atoms(sym.Symbol))]
        customs = []
        for parameter in symbols:
            if parameter not in [rxn.Keq_str, rxn.kf_str, rxn.kr_str] \
               and parameter is not "t":
                try:
                    value = model.custom_parameters[parameter]
                    if value is None:
                        customs.append(parameter)
                except KeyError:
                    customs.append(parameter)
        if existing_customs:
            customs = [custom for custom in customs
                       if sym.Symbol(custom) not in existing_customs
                       or existing_customs[sym.Symbol(custom)] is None]
        missing[rxn] = "; ".join(customs)

    return missing


def get_missing_steady_state_fluxes(model, reaction_list=None):
    """Identify the missing steady state flux values for reactions in a model.

    Parameters
    ----------
    model: mass.MassModel
        The MassModel to inspect.
    reaction_list : list of mass.MassReaction, optional
        A list of mass.MassReaction objects in the model to be checked.
        If None provided, will use all reactions in the model.

    Returns
    -------
    missing: list
        A list of reactions with missing steady state fluxes values. Returns as
        empty if there are no missing values.

    See Also
    --------
    qcqa_model

    """
    if reaction_list is None:
        reaction_list = model.reactions
    reaction_list = ensure_iterable(reaction_list)

    missing = [rxn for rxn in reaction_list if rxn.steady_state_flux is None]

    return missing


def get_missing_initial_conditions(model, simulation=None,
                                   metabolite_list=None):
    """Identify the missing initial conditions for metabolites in a model.

    Parameters
    ----------
    model: mass.MassModel
        The MassModel to inspect.
    metabolite_list : list of mass.MassMetabolites, optional
        A list of mass.MassMetabolite objects in the model to be checked.
        If None provided, will use all metabolites in the model.

    Returns
    -------
    missing: list
        A list of metabolites with missing initial conditions. Returns as
        empty if there are no missing values.

    See Also
    --------
    qcqa_model

    """
    if metabolite_list is None:
        metabolite_list = model.metabolites
    metabolite_list = ensure_iterable(metabolite_list)

    # Filter out fixed concentration metabolites
    missing = [met for met in metabolite_list
               if met not in model.fixed_concentrations]

    missing = [met for met in missing if met not in model.initial_conditions
               or model.initial_conditions[met] is None]

    if simulation is not None:
        exist = simulation.view_initial_concentration_values(model)
        missing = [met for met in missing
                   if not sym.Function(str(met))(_T_SYM) in exist[model.id]
                   or exist[model.id][sym.Function(str(met))(_T_SYM)] is None]

    return missing


def get_missing_fixed_concentrations(model, simulation=None,
                                     metabolite_list=None):
    """Identify the missing fixed concentrations for metabolites in a model.

    Parameters
    ----------
    model: mass.MassModel
        The MassModel to inspect.
    metabolite_list : list of mass.MassMetabolites, optional
        A list of mass.MassMetabolite objects in the model to be checked.
        If None provided, will use all external metabolites in the model.

    Returns
    -------
    missing: list
        A list of metabolites with missing fixed concentrations. Returns as
        empty if there are no missing values.

    See Also
    --------
    qcqa_model

    """
    if metabolite_list is None:
        metabolite_list = model.external_metabolites
    metabolite_list = ensure_iterable(metabolite_list)

    # Filter out initial concentrations
    missing = [met for met in metabolite_list
               if met not in model.initial_conditions]

    missing = [met for met in missing if met not in model.fixed_concentrations
               or model.fixed_concentrations[met] is None]

    if simulation is not None:
        existing = simulation.view_parameter_values(model)
        missing = [met for met in missing
                   if not sym.Symbol(str(met)) in existing[model.id]
                   or existing[model.id][sym.Symbol(str(met))] is None]

    return missing


def check_superfluous_consistency(model, simulation=None, tol=None,
                                  reaction_list=None):
    """Check parameters of model reactions to ensure numerical consistentency.

    Parameter numerical consistency includes checking reaction rate constants
    and equilibrium constants to ensure they are mathematically consistent with
    one another. If there are no superfluous parameters, existing parameters
    are considered consistent.

    Parameters
    ----------
    model: mass.massmodel
        The MassModel to inspect
    tol: float, optional
        The tolerance for parameter consistency. Parameters are considered
        consistent if abs(rxn.kr - rxn.kf/rxn.Keq) <=tol. If None provided, the
        global zero tolerance is used.
    reaction_list: list of mass.MassReaction, optional
        A list of mass.MassReaction objects in the model to be checked.
        If None provided, will use all reactions in the model.

    Returns
    -------
    inconsistent: dict
        A dictionary with MassReactions objects as keys and a string
        identifying the inconsistencies as values. Returns as empty if there
        are no missing values.

    See Also
    --------
    qcqa_model

    """
    tol = _set_tolerance(tol)
    if simulation is not None:
        existing = simulation.view_parameter_values(model)[model.id]

    if reaction_list is None:
        reaction_list = model.reactions
    reaction_list = ensure_iterable(reaction_list)

    superfluous = {}
    for rxn in reaction_list:
        try:
            if simulation is not None:
                keys = [sym.Symbol("{0}_{1}".format(param_type, rxn.id))
                        for param_type in ["kf", "Keq", "kr"]]
                kf, Keq, kr = [existing[key] for key in keys]
            else:
                keys = [rxn.kf_str, rxn.Keq_str, rxn.kr_str]
                kf, Keq, kr = [rxn.parameters[key] for key in keys]
            superfluous[rxn] = _is_consistent(kf, Keq, kr, tol)
        except KeyError:
            pass

    return superfluous


def check_elemental_consistency(model, reaction_list=None):
    """Check the reactions in the model to ensure elemental consistentency.

    Elemental consistency includes checking reactions to ensure they are mass
    and charged balanced. Exchange reactions are ignored because they are
    typically unbalanced.

    Parameters
    ----------
    model: mass.MassModel
        The MassModel to inspect.
    reaction_list : list of mass.MassReaction, optional
        A list of mass.MassReaction objects in the model to be checked.
        If None provided, will use all reactions in the model.

    Returns
    -------
    inconsistent: dict
        A dictionary with MassReactions objects as keys and a string
        identifying the inconsistencies as values. Returns as empty if there
        are no missing values.

    See Also
    --------
    qcqa_model

    """
    if reaction_list is None:
        reaction_list = model.reactions
    reaction_list = ensure_iterable(reaction_list)

    inconsistent = {}
    for reaction in reaction_list:
        if not reaction.exchange and reaction.check_mass_balance():
            unbalanced = ""
            for elem, amount in iteritems(reaction.check_mass_balance()):
                unbalanced += "{0}: {1:.1f}; ".format(elem, amount)
            inconsistent[reaction] = unbalanced.rstrip("; ")

    return inconsistent


def check_thermodynamic_consistency(model, tol=None, reaction_list=None):
    """Check the model reactions for thermodynamic consistency.

    Parameters
    ----------
    model: mass.massmodel
        The MassModel to inspect
    tol: float, optional
        The tolerance for parameter consistency. Parameters are considered
        consistent if abs(rxn.kr - rxn.kf/rxn.Keq) <=tol. If None provided, the
        global zero tolerance is used.
    reaction_list: list of mass.MassReaction, optional
        A list of mass.MassReaction objects in the model to be checked.
        If None provided, will use all reactions in the model.

    Returns
    -------
    inconsistent: dict
        A dictionary with MassReactions objects as keys and a string
        identifying the inconsistencies as values. Returns as empty if there
        are no missing values.

    See Also
    --------
    qcqa_model

    """
    tol = _set_tolerance(tol)
    if reaction_list is None:
        reaction_list = model.reactions
    reaction_list = ensure_iterable(reaction_list)

    inconsistent = {}
    warnings.warn("Thermodynamic consistency checking will be implemented"
                  " in a future updated", FutureWarning)
    return inconsistent


def check_reaction_parameters(model, simulation=None, tol=None,
                              reaction_list=None):
    """Check the model reactions for missing and superfluous parameters.

    Parameters
    ----------
    model: mass.MassModel
        The MassModel to inspect.

    Returns
    -------
    missing: dict
        A dictionary with MassReactions objects as keys and missing parameters
        as values. Returns as empty if there are no missing values.
    superfluous: dict
        A dictionary with MassReactions objects as keys and superfluous
        parameters as values. Returns as empty if there are no superfluous
        parameters.
    reaction_list: list of mass.MassReaction, optional
        A list of mass.MassReaction objects in the model to be checked.
        If None provided, will use all reactions in the model.

    See Also
    --------
    qcqa_model

    """
    tol = _set_tolerance(tol)
    if reaction_list is None:
        reaction_list = model.reactions
    reaction_list = ensure_iterable(reaction_list)

    if simulation is not None:
        existing_parameters = simulation.view_parameter_values(model)[model.id]
        count = _mk_simulation_param_counter(model, simulation, reaction_list)
    else:
        count = []

    missing = []
    superfluous = []
    customs = {}
    for rxn in reaction_list:
        if rxn in model.custom_rates:
            missing_customs = _check_custom_for_standard(model, rxn)
            if simulation is not None:
                param_keys = [sym.Symbol("{0}_{1}".format(param, rxn.id))
                              for param in missing_customs[rxn].split("; ")]
                missing_customs = [str(param).split("_")[0]
                                   for param in param_keys
                                   if param not in existing_parameters
                                   or existing_parameters[param] is None]
            customs.update({rxn: "; ".join(missing_customs)})
        # Address reactions that are missing parameters
        elif (len(rxn.parameters) < 2 and not count) or \
             (isinstance(count, dict) and count[rxn.id] < 2):
            missing.append(rxn)
        # Address reactions that have superfluous parameters
        elif (len(rxn.parameters) > 2 and not count) or \
             (isinstance(count, dict) and count[rxn.id] > 2):
            superfluous.append(rxn)
        # Only two reaction parameters exist, no consistency check required
        else:
            pass

    if missing:
        missing = get_missing_reaction_parameters(model, simulation, missing)
    else:
        missing = {}

    if superfluous:
        superfluous = check_superfluous_consistency(model, simulation,
                                                    tol, superfluous)
    else:
        superfluous = {}

    missing.update(customs)

    return missing, superfluous


def is_simulatable(model, simulation=None):
    """Determine whether a model can be simulated.

    Parameters
    ----------
    model: mass.MassModel
        The MassModel to inspect.
    simulation: mass.Simulation, optional
        If provided, will check whether the model in the given Simulation
        object can be simulated.

    Returns
    -------
    simulate_check: bool
        True if the model can be simulated, False otherwise.
    consistency_check: bool
        True if the model has no numerical consistency issues, False otherwise.

    See Also
    --------
    qcqa_model
    qcqa_simulation

    """
    missing_params, superfluous = check_reaction_parameters(model, simulation)
    missing_concs = get_missing_initial_conditions(model, simulation)
    missing_concs += get_missing_fixed_concentrations(model, simulation)
    missing_params.update(get_missing_custom_parameters(model, simulation))
    consistency_check = True
    if superfluous:
        for rxn, consistency in iteritems(superfluous):
            if consistency is "Inconsistent":
                consistency_check = False

    if missing_params or missing_concs:
        simulate_check = False
    else:
        simulate_check = True

    return [simulate_check, consistency_check]


# Internal
def _mk_parameter_content(model, tol, parameters, fluxes, simulation=None):
    """Create the content for summarizing missing reaction parameters.

    Warnings
    --------
    This method is intended for internal use only.

    """
    missing = []
    headers = []
    # Check standard reaction parameters if desired.
    if parameters:
        headers.append("Reaction Parameters")
        missing_params = check_reaction_parameters(model, simulation, tol)[0]
        missing_params = ["{0}: {1}".format(rxn.id, params)
                          for rxn, params in iteritems(missing_params)]
        missing.append("\n".join(missing_params))
        # Check custom parameters
        headers.append("Custom Parameters")
        missing_params = get_missing_custom_parameters(model, simulation)
        missing_params = ["{0}: {1}".format(rxn.id, params)
                          for rxn, params in iteritems(missing_params)]
        missing.append("\n".join(missing_params))

    # Check steady state fluxes if desired.
    if fluxes:
        headers.append("S.S. Fluxes")
        missing_params = get_missing_steady_state_fluxes(model)
        missing.append("\n".join([r.id for r in missing_params]))

    section = "MISSING PARAMETERS"
    content_lists, columns, sections = _mk_content(missing, headers, section)
    return content_lists, columns, sections


def _mk_concentration_content(model, simulation=None):
    """Create the content for summarizing missing concentrations.

    Warnings
    --------
    This method is intended for internal use only.

    """
    missing = []
    for function in [get_missing_initial_conditions,
                     get_missing_fixed_concentrations]:
        missing_conc = [str(m) for m in function(model, simulation)]
        missing.append("\n".join(missing_conc))

    headers = ["Initial Conditions", "Fixed Concentrations"]
    section = "MISSING CONCENTRATIONS"
    content_lists, columns, sections = _mk_content(missing, headers, section)
    return content_lists, columns, sections


def _mk_consistency_content(model, tol, superfluous, elemental, thermodynamic,
                            simulation=None):
    """Create the content for summarizing missing reaction parameters.

    Warnings
    --------
    This method is intended for internal use only.

    """
    missing = []
    headers = []
    # Check superfluous parameters and their consistency if desired
    if superfluous:
        headers.append("Superfluous Parameters")
        inconsistent = check_reaction_parameters(model, simulation, tol)[1]
        inconsistent = ["{0}: {1}".format(rxn.id, consistency)
                        for rxn, consistency in iteritems(inconsistent)]
        missing.append("\n".join(inconsistent))
    # Check elemental consistency if desired
    if elemental:
        headers.append("Elemental")
        inconsistent = check_elemental_consistency(model)
        inconsistent = ["{0}: {{{1}}} unbalanced".format(reaction.id, unbal)
                        for reaction, unbal in iteritems(inconsistent)]
        missing.append("\n".join(inconsistent))

    # Check thermodynamic consistency if desired
    if thermodynamic:
        headers.append("Thermodynamic")
        inconsistent = check_thermodynamic_consistency(model, tol)
        missing.append("\n".join(inconsistent))

    section = "CONSISTENCY CHECKS"
    content_lists, columns, sections = _mk_content(missing, headers, section)
    return content_lists, columns, sections


def _mk_content(missing, headers, section):
    """Check if content exists and add to table setup lists if it does.

    Warnings
    --------
    This method is intended for internal use only.

    """
    content_lists = []
    columns = []
    sections = []
    for content, head in zip(missing, headers):
        if content:
            content_lists.append(content)
            columns.append(head)
    if content_lists and columns:
        content_lists = [content_lists]
        columns = [columns]
        sections.append(section)

    return content_lists, columns, sections


def _format_table_for_print(table_items, checks, model_id, sim_id=None):
    """Format qcqa report table such that it is ready to be printed.

    Warnings
    --------
    This method is intended for internal use only.

    """
    def make_formatted_table(content, header_list, table_format,
                             str_alignment):
        formatted_table = tabulate(content, headers=header_list,
                                   tablefmt=table_format,
                                   stralign=str_alignment)
        return formatted_table

    simulate_check, consistency_check = checks
    # Unpack table items
    content_lists, columns, sections = table_items

    # Create tables
    tables = [make_formatted_table([content], header, 'simple', u'left')
              for content, header in zip(content_lists, columns)]
    # Format based on longest string in the inner tables if content exists
    if tables:
        # Determine longest line in the table, minimum length of 42 characters
        max_len = max([len(table.split('\n')[1]) for table in tables] + [42])
        sections = [["{0}{1}{2}".format(" " * ceil((max_len - len(section))/2),
                    section, " " * floor((max_len - len(section))/2))]
                    for section in sections]

    # Format all indivual pieces of the report
    tables = [make_formatted_table([[table]], section, 'rst', u'left')
              for table, section in zip(tables, sections)]
    tables = [[table] for table in tables]
    report_head = ""
    if sim_id is not None:
        report_head += "SIMULATION OBJECT ID: {0}\nSPECIFIC ".format(sim_id)

    # Create and print report
    report_head += ("MODEL ID: {0}\nSIMULATABLE: {1};\nNUMERICAL CONSISTENCY: "
                    "{2}".format(model_id, simulate_check, consistency_check))
    report = make_formatted_table(tables, [report_head], 'fancy_grid', u'left')
    return report


def _check_custom_for_standard(model, reaction):
    """Check for missing standard reaction parameters in custom rate laws.

    Warnings
    --------
    This method is intended for internal use only.

    """
    customs = {}
    if reaction in model.custom_rates:
        symbols = list(model.custom_rates[reaction].atoms(sym.Symbol))
        symbols = sorted([str(s) for s in symbols
                         if str(s) in [reaction.Keq_str, reaction.kf_str,
                                       reaction.kr_str]])
        for param in symbols:
            try:
                reaction.parameters[param]
            except KeyError:
                if reaction not in customs:
                    customs[reaction] = "{0}; ".format(param.split("_")[0])
                else:
                    customs[reaction] += "{0}; ".format(param.split("_")[0])
        if reaction in customs:
            customs[reaction] = customs[reaction].rstrip("; ")
    return customs


def _mk_simulation_param_counter(model, simulation, reaction_list):
    """Make a reaction parameter counter for a model in a Simulation.

    Warnings
    --------
    This method is intended for internal use only.

    """
    count = []
    if simulation is not None:
        existing_parameters = simulation.view_parameter_values(model)[model.id]
        for param in iterkeys(existing_parameters):
            param_split = str(param).split("_", 1)
            if param_split[0] in ["kf", "Keq", "kr"] and \
               param_split[1] in [r.id for r in reaction_list]:
                count.append(param_split[1])
            elif str(param) not in model.external_metabolites:
                count.append(str(param))
        count = Counter(count)
    return count


def _is_consistent(kf, Keq, kr, tol):
    """Determine whether the reaction parameters are numerically consistency.

    Warnings
    --------
    This method is intended for internal use only.

    """
    return "Consistent" if (abs(kr - kf/Keq) <= tol) else "Inconsistent"


def _set_tolerance(tol):
    """Check value type and return the value for the zero tolerance.

    Warnings
    --------
    This method is intended for internal use only.

    """
    if tol is None:
        return _ZERO_TOL
    elif not isinstance(tol, float):
        raise TypeError("tol must be a float")
    else:
        return tol
