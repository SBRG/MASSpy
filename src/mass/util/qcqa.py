# -*- coding: utf-8 -*-
"""Module containing functions to assess the quality of a model."""
from math import ceil, floor

import sympy as sym
from cobra.util.util import format_long_string
from six import iteritems, itervalues, string_types
from tabulate import tabulate

from mass.core.mass_configuration import MassConfiguration
from mass.util.expressions import _mk_met_func
from mass.util.util import _check_kwargs, ensure_iterable


# Global
MASSCONFIGURTION = MassConfiguration()


def qcqa_model(model, **kwargs):
    """Check the model quality and print a summary of the results.

    Notes
    -----
    Checking the model quality involves running a series of quality control and
    assessment tests to determine consistency (e.g. elemental)
    in the model, missing values, and whether the model can be simulated.

    Parameters
    ----------
    model : MassModel
        The model to inspect.
    **kwargs
        parameters :
            ``bool`` indicating whether to check for undefined parameters in
            the model.

            Default is ``False.``
        concentrations :
            ``bool`` indicating whether to check for undefined initial and
            boundary conditions in the model.

            Default is ``False.``
        fluxes :
            ``bool`` indicating whether to check for undefined steady state
            fluxes in the model.

             Default is ``False.``
        superfluous :
            ``bool`` indicating whether to check for superfluous parameters
            in the model and ensure existing parameters are consistent with
            one another if superfluous parameters are present.

             Default is ``False.``
        elemental :
            ``bool`` indicating whether to check for elemental consistency in
            the model. Boundary reactions are ignored.

            Default is ``False.``
        simulation_only :
            Only check for undefined values necessary for simulating the model.

            Default is ``True``.

    """
    kwargs = _check_kwargs(
        {
            "parameters": False,
            "concentrations": False,
            "fluxes": False,
            "superfluous": False,
            "elemental": False,
            "simulation_only": True,
        },
        kwargs,
    )

    # Set up empty lists for storing QC/QA report items.
    table_items = [[], [], []]
    # Get missing parameters
    if any([kwargs.get(k) for k in ["parameters", "fluxes"]]):
        results = _mk_parameter_content(model, **kwargs)
        for to_add, item_list in zip(results, table_items):
            item_list.extend(to_add)

    # Get missing initial and fixed concentrations
    if kwargs.get("concentrations"):
        results = _mk_concentration_content(model, **kwargs)
        for to_add, item_list in zip(results, table_items):
            item_list.extend(to_add)

    # Check for the desired consistencies in the values.
    if any([kwargs.get(k) for k in ["superfluous", "elemental"]]):
        results = _mk_consistency_content(model, **kwargs)
        for to_add, item_list in zip(results, table_items):
            item_list.extend(to_add)

    # Check if simulatable
    checks = is_simulatable(model)
    report = _format_table_for_print(table_items, checks, model.id)
    print(report)


def get_missing_reaction_parameters(model, reaction_list=None, simulation_only=True):
    r"""Identify the missing parameters for reactions in a model.

    Notes
    -----
    Will include the default reaction parameters in custom rate laws. To get
    missing custom parameters for reactions with custom rate expressions,
    use :func:`get_missing_custom_parameters` instead.

    Parameters
    ----------
    model : MassModel
        The model to inspect.
    reaction_list : iterable
        An iterable of :class:`~.MassReaction`\ s in the model to be checked.
        If ``None`` then all reactions in the model will be utilized.
    simulation_only :
        Only check for undefined values necessary for simulating the model.

    Returns
    -------
    missing : dict
        A ``dict`` with :class:`~.MassReaction`\ s as keys and a string
        identifying the missing parameters as values. Will return as an
        empty ``dict`` if there are no missing values.

    See Also
    --------
    :attr:`.MassReaction.all_parameter_ids`
        List of default reaction parameters.

    """
    reaction_list = _get_objs_to_check(model, "reactions", reaction_list)
    missing = {}
    for rxn in reaction_list:
        missing_params = []
        parameter_keys = [rxn.Keq_str, rxn.kf_str, rxn.kr_str]
        for key in parameter_keys:
            try:
                rxn.parameters[key]
            except KeyError:
                if not rxn.reversible and key in [rxn.Keq_str, rxn.kr_str]:
                    pass
                else:
                    missing_params.append(key)
        missing_params = "; ".join([k.split("_")[0] for k in missing_params]).rstrip(
            "; "
        )

        # Remove missing equilibrium and reverse rate constants
        # for irreversible reactions
        for param in ["Keq", "kr"]:
            if not rxn.reversible and param in missing_params:
                missing_params = missing_params.replace(param, "")

        if missing_params:
            missing[rxn] = "{0}".format(missing_params.rstrip("; "))

    if simulation_only and missing:
        missing = _check_if_param_needed(model, missing)

    return missing


def get_missing_custom_parameters(model, reaction_list=None, simulation_only=True):
    r"""Identify the missing custom parameters in a model.

    Notes
    -----
    Will not include default reaction parameters. To get missing standard
    reaction parameters for reactions with custom rate laws, use
    :func:`get_missing_reaction_parameters` instead.

    Parameters
    ----------
    model : MassModel
        The model to inspect.
    reaction_list : iterable
        An iterable of :class:`~.MassReaction`\ s in the model to be checked.
        If ``None`` then all reactions in the model will be utilized.
    simulation_only :
        Only check for undefined values necessary for simulating the model.

    Returns
    -------
    missing : dict
        A ``dict`` with :class:`~.MassReaction`\ s as keys and a string
        identifying the missing custom parameters as values. Will return as an
        empty ``dict`` if there are no missing values.

    See Also
    --------
    :attr:`.MassReaction.all_parameter_ids`
        List of default reaction parameters.

    """
    reaction_list = _get_objs_to_check(model, "reactions", reaction_list)

    missing = {}
    # Filter out reactions without custom rates
    reaction_list = [
        reaction for reaction in reaction_list if reaction in model.custom_rates
    ]
    for rxn in reaction_list:
        rate = model.custom_rates[rxn]
        symbols = [
            str(symbol)
            for symbol in list(rate.atoms(sym.Symbol))
            if str(symbol) not in model.metabolites
        ]
        customs = []
        for parameter in symbols:
            if (
                parameter not in [rxn.Keq_str, rxn.kf_str, rxn.kr_str]
                and parameter != "t"
            ):
                try:
                    value = model.custom_parameters[parameter]
                    if value is None:
                        customs.append(parameter)
                except KeyError:
                    if parameter not in model.boundary_conditions:
                        customs.append(parameter)
        if customs:
            missing[rxn] = "; ".join(customs)

    if simulation_only and missing:
        missing = _check_if_param_needed(model, missing, customs=True)

    return missing


def get_missing_steady_state_fluxes(model, reaction_list=None):
    r"""Identify the missing steady state flux values for reactions in a model.

    Parameters
    ----------
    model : MassModel
        The model to inspect.
    reaction_list : iterable
        An iterable of :class:`~.MassReaction`\ s in the model to be checked.
        If ``None`` then all reactions in the model will be utilized.

    Returns
    -------
    missing : list
        List of :class:`~.MassReaction`\ s with missing steady state fluxes.
        Will return as an empty ``list`` if there are no missing values.

    """
    reaction_list = _get_objs_to_check(model, "reactions", reaction_list)

    missing = [rxn for rxn in reaction_list if rxn.steady_state_flux is None]

    return missing


def get_missing_initial_conditions(model, metabolite_list=None, simulation_only=True):
    r"""Identify the missing initial conditions for metabolites in a model.

    Notes
    -----
    Does not include boundary conditions.

    Parameters
    ----------
    model : MassModel
        The model to inspect.
    metabolite_list : iterable
        An iterable of :class:`~.MassMetabolite`\ s in the model to be checked.
        If ``None`` then all metabolites in the model will be utilized.
    simulation_only :
        Only check for undefined values necessary for simulating the model.

    Returns
    -------
    missing : list
        List of :class:`~.MassMetabolite`\ s with missing initial conditions.
        Will return as an empty ``list`` if there are no missing values.

    """
    metabolite_list = _get_objs_to_check(model, "metabolites", metabolite_list)

    # Filter out 'boundary metabolites'
    missing = [met for met in metabolite_list if met not in model.boundary_conditions]

    missing = [
        met
        for met in missing
        if met not in model.initial_conditions or model.initial_conditions[met] is None
    ]

    if simulation_only and missing:
        missing = _check_if_conc_needed(model, missing)

    return missing


def get_missing_boundary_conditions(model, metabolite_list=None, simulation_only=True):
    r"""Identify the missing boundary conditions for metabolites in a model.

    Parameters
    ----------
    model : MassModel
        The model to inspect.
    metabolite_list : iterable
        An iterable of 'boundary metabolites' or :class:`~.MassMetabolite`\ s
        in the model to be checked. If ``None`` then all 'boundary metabolites'
        in the model will be utilized.
    simulation_only :
        Only check for undefined values necessary for simulating the model.

    Returns
    -------
    missing : list
        List of metabolites with missing boundary conditions.
        Will return as an empty ``list`` if there are no missing values.

    See Also
    --------
    :attr:`.MassModel.boundary_metabolites`
        List of boundary metabolites found in the model.

    """
    if metabolite_list is None:
        metabolite_list = model.boundary_metabolites
    metabolite_list = ensure_iterable(metabolite_list)

    # Filter out initial concentrations
    missing = [met for met in metabolite_list if met not in model.initial_conditions]

    missing = [
        met
        for met in missing
        if met not in model.boundary_conditions
        or model.boundary_conditions[met] is None
    ]

    if simulation_only and missing:
        missing = _check_if_conc_needed(model, missing)

    return missing


def check_superfluous_consistency(model, reaction_list=None):
    r"""Check parameters of model reactions to ensure numerical consistentency.

    Parameter numerical consistency includes checking reaction rate and
    equilibrium constants to ensure they are mathematically consistent with
    one another. If there are no superfluous parameters, the existing
    parameters are considered consistent.

    Notes
    -----
    The `MassConfiguration.decimal_precision` is used to round the value of
    ``abs(rxn.kr - rxn.kf/rxn.Keq)`` before comparison.

    Parameters
    ----------
    model : MassModel
        The model to inspect.
    reaction_list : iterable
        An iterable of :class:`~.MassReaction`\ s in the model to be checked.
        If ``None`` then all reactions in the model will be utilized.

    Returns
    -------
    inconsistent : dict
        A ``dict`` with :class:`~.MassReaction`\ s as keys and a string
        identifying the incosistencies as values. Will return as an
        empty ``dict`` if there are no inconsistencies.

    """
    reaction_list = _get_objs_to_check(model, "reactions", reaction_list)

    superfluous = {}
    for rxn in reaction_list:
        try:
            args = [
                rxn.parameters[key] for key in [rxn.kf_str, rxn.Keq_str, rxn.kr_str]
            ]
            superfluous[rxn] = _is_consistent(*args)
        except KeyError:
            pass

    return superfluous


def check_elemental_consistency(model, reaction_list=None):
    r"""Check the reactions in the model to ensure elemental consistentency.

    Elemental consistency includes checking reactions to ensure they are mass
    and charged balanced. Boundary reactions are ignored because they are
    typically unbalanced.

    Parameters
    ----------
    model : MassModel
        The model to inspect.
    reaction_list : iterable
        An iterable of :class:`~.MassReaction`\ s in the model to be checked.
        If ``None`` then all reactions in the model will be utilized.

    Returns
    -------
    inconsistent : dict
        A ``dict`` with :class:`~.MassReaction`\ s as keys and a string
        identifying the incosistencies as values. Will return as an
        empty ``dict`` if there are no inconsistencies.

    """
    reaction_list = _get_objs_to_check(model, "reactions", reaction_list)

    inconsistent = {}
    for reaction in reaction_list:
        if not reaction.boundary and reaction.check_mass_balance():
            unbalanced = ""
            for elem, amount in iteritems(reaction.check_mass_balance()):
                unbalanced += "{0}: {1:.1f}; ".format(elem, amount)
            inconsistent[reaction] = unbalanced.rstrip("; ")

    return inconsistent


def check_reaction_parameters(model, reaction_list=None, simulation_only=True):
    r"""Check the model reactions for missing and superfluous parameters.

    Parameters
    ----------
    model : MassModel
        The model to inspect.
    reaction_list : iterable
        An iterable of :class:`~.MassReaction`\ s in the model to be checked.
        If ``None`` then all reactions in the model will be utilized.
    simulation_only :
        Only check for undefined values necessary for simulating the model.

    Returns
    -------
    tuple (missing, superfluous)
    missing : dict
        A ``dict`` with :class:`~.MassReaction`\ s as keys and a string
        identifying the missing parameters as values. Will return as an
        empty ``dict`` if there are no missing values.
    superfluous : dict
        A ``dict`` with :class:`~.MassReaction`\ s as keys and superfluous
        parameters as values. Will return as an empty ``dict`` if there are no
        superfluous values.

    """
    reaction_list = _get_objs_to_check(model, "reactions", reaction_list)

    missing = []
    superfluous = []
    customs = {}
    for rxn in reaction_list:
        if rxn in model.custom_rates:
            missing_customs = _check_custom_for_standard(model, rxn)
            if missing_customs:
                customs.update(
                    dict(
                        (rxn, "; ".join([missing]))
                        if isinstance(missing, string_types)
                        else (rxn, "; ".join(missing))
                        for rxn, missing in iteritems(missing_customs)
                    )
                )
        # Always check if forward rate constant defined
        elif rxn.forward_rate_constant is None:
            missing.append(rxn)
        # Reversible reaction without an equilibrium or reverse rate constant
        elif rxn.reversible and len(rxn.parameters) < 2:
            missing.append(rxn)
        elif rxn.reversible and len(rxn.parameters) > 2:
            superfluous.append(rxn)
        # Two reaction parameters exist for reversible reactions or
        # forward rate constant exists for an irreversible reaction
        else:
            pass

    if missing and simulation_only:
        missing = get_missing_reaction_parameters(model, missing, simulation_only)
    elif missing:
        missing = get_missing_reaction_parameters(model, None, simulation_only)
    else:
        missing = {}

    if superfluous:
        superfluous = check_superfluous_consistency(model, superfluous)
    else:
        superfluous = {}
    missing.update(customs)

    return missing, superfluous


def is_simulatable(model):
    """Determine whether a model can be simulated.

    Parameters
    ----------
    model : MassModel
        The model to inspect.

    Returns
    -------
    tuple (simulate_check, consistency_check)
    simulate_check : bool
        ``True`` if the model can be simulated, ``False`` otherwise.
    consistency_check : bool
        ``True`` if the model has no issues with numerical consistency,
        ``False`` otherwise.

    """
    missing_params, superfluous = check_reaction_parameters(model)
    missing_concs = get_missing_initial_conditions(model)
    missing_concs += get_missing_boundary_conditions(model)
    missing_params.update(get_missing_custom_parameters(model))
    consistency_check = True
    if superfluous:
        for consistency in itervalues(superfluous):
            if consistency == "Inconsistent":
                consistency_check = False

    if missing_params or missing_concs:
        simulate_check = False
    else:
        simulate_check = True

    return (simulate_check, consistency_check)


# Internal
def _mk_parameter_content(model, **kwargs):
    """Create the content for summarizing missing reaction parameters.

    Warnings
    --------
    This method is intended for internal use only.

    """
    parameters, fluxes = tuple(kwargs.get(k) for k in ["parameters", "fluxes"])

    missing = []
    headers = []
    # Check standard reaction parameters if desired.
    if parameters:
        headers.append("Reaction Parameters")
        missing_params = check_reaction_parameters(
            model, simulation_only=kwargs.get("simulation_only")
        )[0]
        missing_params = [
            "{0}: {1}".format(rxn.id, params)
            for rxn, params in iteritems(missing_params)
        ]
        missing.append("\n".join(missing_params))
        # Check custom parameters
        headers.append("Custom Parameters")
        missing_params = get_missing_custom_parameters(
            model, simulation_only=kwargs.get("simulation_only")
        )
        missing_params = [
            "{0}: {1}".format(rxn.id, params)
            for rxn, params in iteritems(missing_params)
        ]
        missing.append("\n".join(missing_params))

    # Check steady state fluxes if desired.
    if fluxes:
        headers.append("S.S. Fluxes")
        missing_params = get_missing_steady_state_fluxes(model)
        missing.append("\n".join([r.id for r in missing_params]))

    section = "MISSING PARAMETERS"
    content_lists, columns, sections = _mk_content(missing, headers, section)
    return content_lists, columns, sections


def _mk_concentration_content(model, **kwargs):
    """Create the content for summarizing missing concentrations.

    Warnings
    --------
    This method is intended for internal use only.

    """
    missing = []
    for i, function in enumerate(
        [get_missing_initial_conditions, get_missing_boundary_conditions]
    ):
        missing_conc = [
            m for m in function(model, simulation_only=kwargs.get("simulation_only"))
        ]
        for j, met in enumerate(missing_conc):
            if i == 0:
                # Identify reactions for missing initial conditions
                associated_rxns = sorted([r.id for r in met.reactions])
            else:
                # Identify reactions for missing boundary conditions
                associated_rxns = sorted(
                    [r.id for r in model.boundary if r.boundary_metabolite == met]
                )
            # Format string
            associated_rxn_str = ", ".join(associated_rxns)
            missing_conc[j] = "{0} (in {1})".format(
                str(met), format_long_string(associated_rxn_str, 30)
            )
        # Join all strings together.
        missing.append("\n".join(missing_conc))

    headers = ["Initial Conditions", "Boundary Conditions"]
    section = "MISSING CONCENTRATIONS"
    content_lists, columns, sections = _mk_content(missing, headers, section)
    return content_lists, columns, sections


def _mk_consistency_content(model, **kwargs):
    """Create the content for summarizing missing reaction parameters.

    Warnings
    --------
    This method is intended for internal use only.

    """
    superfluous, elemental = tuple(kwargs.get(k) for k in ["superfluous", "elemental"])

    missing = []
    headers = []
    # Check superfluous parameters and their consistency if desired
    if superfluous:
        headers.append("Superfluous Parameters")
        inconsistent = check_reaction_parameters(
            model, simulation_only=kwargs.get("simulation_only")
        )[1]
        inconsistent = [
            "{0}: {1}".format(rxn.id, consistency)
            for rxn, consistency in iteritems(inconsistent)
        ]
        missing.append("\n".join(inconsistent))
    # Check elemental consistency if desired
    if elemental:
        headers.append("Elemental")
        inconsistent = check_elemental_consistency(model)
        inconsistent = [
            "{0}: {{{1}}}".format(reaction.id, unbalanced)
            for reaction, unbalanced in iteritems(inconsistent)
        ]
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


def _format_table_for_print(table_items, checks, model_id):
    """Format qcqa report table such that it is ready to be printed.

    Warnings
    --------
    This method is intended for internal use only.

    """

    def make_formatted_table(content, header_list, table_format, str_alignment):
        formatted_table = tabulate(
            content, headers=header_list, tablefmt=table_format, stralign=str_alignment
        )
        return formatted_table

    simulate_check, consistency_check = checks
    # Unpack table items
    content_lists, columns, sections = table_items

    # Create tables
    tables = [
        make_formatted_table([content], header, "simple", "left")
        for content, header in zip(content_lists, columns)
    ]
    # Format based on longest string in the inner tables if content exists
    if tables:
        # Determine longest line in the table, minimum length of 42 characters
        max_l = max([len(table.split("\n")[1]) for table in tables] + [42])
        sections = [
            [
                "{0}{1}{2}".format(
                    " " * ceil((max_l - len(sect)) / 2),
                    sect,
                    " " * floor((max_l - len(sect)) / 2),
                )
            ]
            for sect in sections
        ]

    # Format all indivual pieces of the report
    tables = [
        make_formatted_table([[table]], section, "rst", "left")
        for table, section in zip(tables, sections)
    ]
    tables = [[table] for table in tables]
    report_head = ""

    # Create and print report
    report_head += (
        "MODEL ID: {0}\nSIMULATABLE: {1}\nPARAMETERS NUMERICALY CONSISTENT:"
        " {2}".format(model_id, simulate_check, consistency_check)
    )
    report = make_formatted_table(tables, [report_head], "fancy_grid", "left")
    return report


def _check_custom_for_standard(model, reaction):
    """Check for missing standard reaction parameters in custom rate laws.

    Warnings
    --------
    This method is intended for internal use only.

    """
    customs = {}
    if reaction in model.custom_rates and model.custom_rates[reaction] is not None:
        symbols = list(model.custom_rates[reaction].atoms(sym.Symbol))
        symbols = sorted(
            [
                str(s)
                for s in symbols
                if str(s) in [reaction.Keq_str, reaction.kf_str, reaction.kr_str]
            ]
        )
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


def _is_consistent(kf, Keq, kr):
    """Determine whether the reaction parameters are numerically consistency.

    Warnings
    --------
    This method is intended for internal use only.

    """
    if round(abs(kr - (kf / Keq)), MASSCONFIGURTION.decimal_precision) == 0:
        return "Consistent"

    return "Inconsistent"


def _check_if_conc_needed(model, missing):
    """Check whether the missing concentrations are needed for simulation.

    Warnings
    --------
    This method is intended for internal use only.

    """
    needed = set()
    for rate in itervalues(model.rates):
        needed.update(rate.atoms(sym.Function))
        needed.update(rate.atoms(sym.Symbol))

    missing = [
        met
        for met in missing
        if sym.Symbol(str(met)) in needed or _mk_met_func(str(met)) in needed
    ]

    return missing


def _check_if_param_needed(model, missing, customs=False):
    """Check whether the missing parameters are needed for simulation.

    Warnings
    --------
    This method is intended for internal use only.

    """
    needed = set()
    for rate in itervalues(model.rates):
        needed.update(rate.atoms(sym.Symbol))

    for reaction, missing_values_str in iteritems(missing.copy()):
        missing_params = missing_values_str.split("; ")
        if customs:
            missing_params = [
                param for param in missing_params if sym.Symbol(param) in needed
            ]
        else:
            missing_params = [
                param
                for param in missing_params
                if sym.Symbol("_".join((param, reaction.id))) in needed
            ]

        if missing_params:
            missing[reaction] = "; ".join(missing_params)
        else:
            del missing[reaction]

    return missing


def _get_objs_to_check(model, attribute, object_list):
    """Check whether the missing parameters are needed for simulation.

    Warnings
    --------
    This method is intended for internal use only.

    """
    attribute_dictlist = getattr(model, attribute)
    if object_list is not None:
        object_list = ensure_iterable(object_list)
        for i, obj in enumerate(object_list):
            try:
                obj = attribute_dictlist.get_by_id(getattr(obj, "_id", obj))
            except KeyError as e:
                raise ValueError("'{0}' not found in model.".format(str(e)))
            else:
                object_list[i] = obj
    else:
        object_list = attribute_dictlist

    return object_list


__all__ = (
    "qcqa_model",
    "get_missing_reaction_parameters",
    "get_missing_custom_parameters",
    "get_missing_steady_state_fluxes",
    "get_missing_initial_conditions",
    "get_missing_boundary_conditions",
    "check_superfluous_consistency",
    "check_elemental_consistency",
    "check_reaction_parameters",
    "is_simulatable",
)
