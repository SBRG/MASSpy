# -*- coding: utf-8 -*-
r"""Module to create and manage an ensemble of models.

The method of the :mod:`~.ensemble` submodule are designed to generate an
ensemble of models. It contains various methods to assist in generating
multiple models from existing :class:`~.MassModel`\ s, using flux data or
concentration data in :class:`pandas.DataFrame`\ s  (e.g. generated from
:mod:`~mass.thermo.conc_sampling`). There are also methods to help ensure that
models are thermodynamically feasible and can reach steady states with or
without perturbations applied.

In addition to containing various methods that can be combined into
an ensemble generation workflow, the  :mod:`~.ensemble` submodule contains the
:func:`generate_ensemble_of_models` function, which is optimized for
performance when generating a large number of models.

The :func:`generate_ensemble_of_models` function also ensures that the user
input is valid before generating models to reduce the likelihood of a user
error causing the model generation process to stop before completion. However,
there is time spent in function's setup, meaining that when generating a
smaller number of models, performance gains may not be seen.

"""
import logging
import warnings

import numpy as np
import pandas as pd
from six import iteritems, string_types

from mass.core.mass_model import MassModel
from mass.exceptions import MassEnsembleError
from mass.simulation.simulation import (
    STEADY_STATE_SOLVERS,
    Simulation,
    _get_sim_values_from_model,
)
from mass.util.util import _check_kwargs, _log_msg, _make_logger, ensure_iterable


# Set the logger
LOGGER = _make_logger(__name__)
"""logging.Logger: Logger for :mod:`~mass.thermo.ensemble` submodule."""


def create_models_from_flux_data(
    reference_model, data=None, raise_error=False, **kwargs
):
    """Generate ensemble of models for a given set of flux data.

    Parameters
    ----------
    reference_model : iterable, None
        A :class:`.MassModel` object to treat as the reference model.
    data : pandas.DataFrame
        A :class:`pandas.DataFrame` containing the flux data for generation
        of the models. Each row is a different set of flux values to
        generate a model for, and each column corresponds to the reaction
        identifier for the flux value.
    raise_error : bool
        Whether to raise an error upon failing to generate a model from a
        given reference. Default is ``False``.
    **kwargs
        verbose :
            ``bool`` indicating the verbosity of the function.

            Default is ``False``.
        suffix :
            ``str`` representing the suffix to append to generated models.

            Default is ``'_F'``.

    Returns
    -------
    new_models : list
        A ``list`` of successfully generated :class:`.MassModel` objects.

    Raises
    ------
    MassEnsembleError
        Raised if generation of a model fails and ``raise_error=True``.

    """
    kwargs = _check_kwargs(
        {
            "verbose": False,
            "suffix": "_F",
        },
        kwargs,
    )

    if not isinstance(reference_model, MassModel):
        raise TypeError("`reference_model` must be a MassModel")

    data, id_array = _validate_data_input(
        reference_model, data, "reactions", kwargs.get("verbose")
    )
    new_models = []
    for i, values in enumerate(data.values):
        # Create new model
        new_model = reference_model.copy()
        new_model.id += kwargs.get("suffix") + str(i)
        try:
            _log_msg(
                LOGGER,
                logging.INFO,
                kwargs.get("verbose"),
                "New model '%s' created",
                new_model.id,
            )

            # Update the model parameters
            new_model.update_parameters(dict(zip(id_array, values)), verbose=False)
            _log_msg(
                LOGGER,
                logging.INFO,
                kwargs.get("verbose"),
                "Updated flux values for '%s'",
                new_model.id,
            )

            # Add model to the ensemble
            new_models.append(new_model)

        except Exception as e:
            msg = str(
                "Could not create '{0}' for the ensemble due to the "
                "following error: {1!r}".format(new_model.id, str(e))
            )
            if raise_error:
                raise MassEnsembleError(msg)
            _log_msg(LOGGER, logging.ERROR, kwargs.get("verbose"), msg)

    return new_models


def create_models_from_concentration_data(
    reference_model, data=None, raise_error=False, **kwargs
):
    """Generate ensemble of models for a given set of concentration data.

    Parameters
    ----------
    reference_model : iterable, None
        A :class:`.MassModel` object to treat as the reference model.
    data : pandas.DataFrame
        A :class:`pandas.DataFrame` containing the concentration data for
        generation of the models. Each row is a different set of
        concentration values to generate a model for, and each column
        corresponds to the metabolite identifier for the concentraiton
        value.
    raise_error : bool
        Whether to raise an error upon failing to generate a model from a
        given reference. Default is ``False``.
    **kwargs
        verbose :
            ``bool`` indicating the verbosity of the function.

            Default is ``False``.
        suffix :
            ``str`` representing the suffix to append to generated models.

            Default is ``'_C'``.

    Returns
    -------
    new_models : list
        A ``list`` of successfully generated :class:`.MassModel` objects.

    Raises
    ------
    MassEnsembleError
        Raised if generation of a model fails and ``raise_error=True``.

    """
    kwargs = _check_kwargs(
        {
            "verbose": False,
            "suffix": "_C",
        },
        kwargs,
    )

    if not isinstance(reference_model, MassModel):
        raise TypeError("`reference_model` must be a MassModel")

    data, id_array = _validate_data_input(
        reference_model, data, "metabolites", kwargs.get("verbose")
    )
    new_models = []
    for i, values in enumerate(data.values):
        # Create new model
        new_model = reference_model.copy()
        new_model.id += kwargs.get("suffix") + str(i)
        try:
            _log_msg(
                LOGGER,
                logging.INFO,
                kwargs.get("verbose"),
                "New model '%s' created",
                new_model.id,
            )

            # Update the model parameters
            new_model.update_initial_conditions(
                dict(zip(id_array, values)), verbose=False
            )
            _log_msg(
                LOGGER,
                logging.INFO,
                kwargs.get("verbose"),
                "Updated initial conditions for '%s'",
                new_model.id,
            )

            # Add model to the ensemble
            new_models.append(new_model)

        except Exception as e:
            msg = str(
                "Could not create '{0}' for the ensemble due to the "
                "following error: {1!r}".format(new_model.id, str(e))
            )
            if raise_error:
                raise MassEnsembleError(msg)
            _log_msg(LOGGER, logging.ERROR, kwargs.get("verbose"), msg)

    return new_models


def ensure_positive_percs(
    models, reactions=None, raise_error=False, update_values=False, **kwargs
):
    """Seperate models based on whether all calculated PERCs are positive.

    Parameters
    ----------
    models : iterable
        An iterable of :class:`.MassModel` objects to use for PERC
        calculations.
    reactions : iterable
        An iterable of reaction identifiers to calculate the
        PERCs for. If ``None``, all reactions in the model will be used.
    raise_error : bool
        Whether to raise an error upon failing to generate a model from a
        given reference. Default is ``False``.
    update_values : bool
        Whether to update the PERC values for models that generate all
        positive PERCs. Default is ``False``.
    **kwargs
        verbose :
            ``bool`` indicating the verbosity of the function.

            Default is ``False``.
        at_equilibrium_default :
            ``float`` value to set the pseudo-order rate constant if the
            reaction is at equilibrium.

            Default is ``100,000``.

    Returns
    -------
    tuple (positive, negative)
    positive : list
        A ``list`` of :class:`.MassModel` objects whose calculated PERC
        values were postiive.
    negative : list
        A ``list`` of :class:`.MassModel` objects whose calculated PERC
        values were negative.

    Raises
    ------
    MassEnsembleError
        Raised if PERC calculation fails and ``raise_error=True``.

    """
    kwargs = _check_kwargs(
        {
            "verbose": False,
            "at_equilibrium_default": 100000,
        },
        kwargs,
    )
    positive = []
    negative = []

    models = ensure_iterable(models)
    if any([not isinstance(model, MassModel) for model in models]):
        raise TypeError("`models` must be an iterable of MassModels.")

    for model in models:
        model, is_positive = _ensure_positive_percs_for_model(
            model,
            reactions,
            kwargs.get("verbose"),
            raise_error,
            update_values,
            kwargs.get("at_equilibrium_default"),
        )
        if is_positive:
            positive.append(model)
        else:
            negative.append(model)

    _log_msg(
        LOGGER,
        logging.INFO,
        kwargs.get("verbose"),
        "Finished PERC calculations, returning seperated models.",
    )
    return positive, negative


def ensure_steady_state(
    models,
    strategy="simulate",
    perturbations=None,
    solver_options=None,
    update_values=False,
    **kwargs
):
    """Seperate models based on whether a steady state can be reached.

    All ``kwargs`` are passed to :meth:`~.Simulation.find_steady_state`.

    Parameters
    ----------
    models : MassModel, iterable
        A :class:`.MassModel` or an iterable of :class:`.MassModel` objects to
        find a steady state for.
    strategy : str
        The strategy for finding the steady state. Must be one of the
        following:

            * ``'simulate'``
            * ``'nleq1'``
            * ``'nleq2'``

    perturbations : dict
        A ``dict`` of perturbations to incorporate into the simulation.
        Models must reach a steady state with the given pertubration to be
        considered as feasible.
        See :mod:`~.simulation.simulation` documentation for more
        information on valid perturbations.
    solver_options : dict
        A `dict` of options to pass to the solver utilized in determining a
        steady state. Solver options should be for the
        :class:`roadrunner.Integrator` if ``strategy="simulate"``, otherwise
        options should correspond to the :class:`roadrunner.SteadyStateSolver`.
    update_values : bool
        Whether to update the model with the steady state results.
        Default is ``False``. Only updates models that reached steady state.
    **kwargs
        verbose :
            ``bool`` indicating the verbosity of the method.

            Default is ``False``.
        steps :
            ``int`` indicating number of steps at which the output is
            sampled where the samples are evenly spaced and
            ``steps = (number of time points) - 1.``
            Steps and number of time points may not both be specified.
            Only valid for ``strategy='simulate'``.

            Default is ``None``.
        tfinal :
            ``float`` indicating the final time point to use in when
            simulating to long times to find a steady state.
            Only valid for ``strategy='simulate'``.

            Default is ``1e8``.
        num_attempts :
            ``int`` indicating the number of attempts the steady state
            solver should make before determining that a steady state
            cannot be found. Only valid for ``strategy='nleq1'`` or
            ``strategy='nleq2'``.

            Default is ``2``.
        decimal_precision :
            ``bool`` indicating whether to apply the
            :attr:`~.MassBaseConfiguration.decimal_precision` attribute of
            the :class:`.MassConfiguration` to the solution values.

            Default is ``False``.

    Returns
    -------
    tuple (feasible, infeasible)
    feasible : list
        A ``list`` of :class:`.MassModel` objects that could successfully reach
        a steady state.
    infeasible : list
        A ``list`` of :class:`.MassModel` objects that could not successfully
        reach a steady state.

    """
    kwargs = _check_kwargs(
        {
            "verbose": False,
            "steps": None,
            "tfinal": 1e8,
            "num_attempts": 2,
            "decimal_precision": True,
        },
        kwargs,
    )

    models = ensure_iterable(models)
    if any([not isinstance(model, MassModel) for model in models]):
        raise TypeError("`models` must be an iterable of MassModels.")

    # Ensure strategy input is valid
    if strategy not in STEADY_STATE_SOLVERS and strategy != "simulate":
        raise ValueError("Invalid steady state strategy: '{0}'".format(strategy))

    simulation = _initialize_simulation(
        models[0], strategy, solver_options, kwargs.get("verbose")
    )

    if len(models) > 1:
        simulation.add_models(models[1:], verbose=kwargs.get("verbose"))

    conc_sol_list, flux_sol_list = simulation.find_steady_state(
        models, strategy, perturbations, update_values, **kwargs
    )

    feasible = []
    infeasible = []
    for i, model in enumerate(models):
        if len(models) == 1:
            conc_sol, flux_sol = conc_sol_list, flux_sol_list
        else:
            conc_sol, flux_sol = conc_sol_list[i], flux_sol_list[i]
        if conc_sol and flux_sol:
            ics, params = simulation.get_model_simulation_values(model)
            model.update_initial_conditions(ics)
            model.update_parameters(
                {
                    param: value
                    for param, value in params.items()
                    if param in model.reactions.list_attr("flux_symbol_str")
                }
            )
            feasible.append(model)
        else:
            infeasible.append(model)

    _log_msg(
        LOGGER,
        logging.INFO,
        kwargs.get("verbose"),
        "Finished finding steady states, returning seperated models.",
    )

    return feasible, infeasible


def generate_ensemble_of_models(
    reference_model,
    flux_data=None,
    conc_data=None,
    ensure_positive_percs=None,
    strategy=None,
    perturbations=None,
    **kwargs
):
    """Generate an ensemble of models for given data sets.

    This function is optimized for performance when generating a large
    ensemble of models when compared to the combination of various individual
    methods of the :mod:`ensemble` submodule used. However, this function may
    not provide as much control over the process when compared to utilizing a
    combination of other methods defined in the :mod:`ensemble` submodule.

    Notes
    -----
    * Only one data set is required to generate the ensemble, meaning that
      a flux data set can be given without a concentration data set, and vice
      versa.

    * If ``x`` flux data samples and ``y`` concentration data samples are
      provided, ``x * y`` total models will be generated.

    * If models deemed ``infeasible`` are to be returned, ensure the
      ``return_infeasible`` kwarg is set to ``True``.

    Parameters
    ----------
    reference_model : MassModel
        The reference model used in generating the ensemble.
    flux_data : pandas.DataFrame or None
        A :class:`pandas.DataFrame` containing the flux data for generation
        of the models. Each row is a different set of flux values to
        generate a model for, and each column corresponds to the reaction
        identifier for the flux value.
    conc_data : pandas.DataFrame or None
        A :class:`pandas.DataFrame` containing the concentration data for
        generation of the models. Each row is a different set of
        concentration values to generate a model for, and each column
        corresponds to the metabolite identifier for the concentraiton
        value.
    ensure_positive_percs :
        A ``list`` of reactions to calculate PERCs for, ensure they
        are postive, and update feasible models with the new PERC values.
        If ``None``, no PERCs will be checked.
    strategy : str, None
        The strategy for finding the steady state.
        Must be one of the following:

            * ``'simulate'``
            * ``'nleq1'``
            * ``'nleq2'``

        If a ``strategy`` is given, models must reach a steady state to be
        considered feasible. All feasible models are updated to steady state.
        If ``None``, no attempts will be made to determine whether a generated
        model can reach a steady state.
    perturbations : dict
        A ``dict`` of perturbations to incorporate into the simulation,
        or a list of perturbation dictionaries where each ``dict`` is applied
        to a simulation. Models must reach a steady state with all given
        pertubration dictionaries to be considered feasible.
        See :mod:`~.simulation.simulation` documentation for more
        information on valid perturbations.

        Ignored if ``strategy=None``.
    **kwargs
        solver_options :
            ``dict`` of options to pass to the solver utilized in determining
            a steady state. Solver options should be for the
            :class:`roadrunner.Integrator` if ``strategy="simulate"``,
            otherwise options should correspond to the
            :class:`roadrunner.SteadyStateSolver`.

            Default is ``None``.
        verbose :
            ``bool`` indicating the verbosity of the function.

            Default is ``False``.
        decimal_precision :
            ``bool`` indicating whether to apply the
            :attr:`~.MassBaseConfiguration.decimal_precision` attribute of
            the :class:`.MassConfiguration` to the solution values.

            Default is ``False``.
        flux_suffix :
            ``str`` representing the suffix to append to generated models
            indicating the flux data set used.

            Default is ``'_F'``.
        conc_suffix :
            ``str`` representing the suffix to append to generated models
            indicating the conc data set used.

            Default is ``'_C'``.
        at_equilibrium_default :
            ``float`` value to set the pseudo-order rate constant if the
            reaction is at equilibrium.

            Default is ``100,000``. Ignored if ``ensure_positive_percs=None``.
        return_infeasible :
            ``bool`` indicating whether to generate and return an
            :class:`Ensemble` containing the models deemed infeasible.

            Default is ``False``.

    Returns
    -------
    feasible : list
        A ``list`` containing the `MassModel` objects that are deemed
        `feasible` by sucessfully passing through all PERC and simulation
        checks in the ensemble building processes.
    infeasible : list
        A ``list`` containing the `MassModel` objects that are deemed
        `infeasible` by failing to passing through one of the PERC or
        simulation checks in the ensemble building processes.

    """
    # Check all inputs at beginning to ensure that ensemble generation is not
    # disrupted near the end due to invalid input format
    kwargs = _check_kwargs(
        {
            "verbose": False,
            "decimal_precision": False,
            "flux_suffix": "_F",
            "conc_suffix": "_C",
            "at_equilibrium_default": 100000,
            "solver_options": None,
            "return_infeasible": False,
        },
        kwargs,
    )
    verbose = kwargs.pop("verbose")
    _log_msg(LOGGER, logging.INFO, verbose, "Validating input")
    # Validate model input
    if not isinstance(reference_model, MassModel):
        raise TypeError("`reference_model` must be a MassModel.")

    # Validate DataFrame inputs, if any
    if flux_data is not None:
        # Validate flux data if provided
        flux_data, flux_ids = _validate_data_input(
            reference_model, flux_data, "reactions", verbose
        )
    else:
        # Set a value to allow for iteration
        flux_data = pd.DataFrame([0])
        flux_ids = np.array([])

    if conc_data is not None:
        # Validate conc data if provided
        conc_data, conc_ids = _validate_data_input(
            reference_model, conc_data, "metabolites", verbose
        )
    else:
        # Set a value to allow for iteration
        conc_data = pd.DataFrame([0])
        conc_ids = np.array([])

    if conc_ids.size == 0 and flux_ids.size == 0:
        raise ValueError(
            "No flux data or concentration data provided. "
            "At least one data set must be provided to generate "
            "the ensemble."
        )
    # Ensure strategy input and perturbation input is valid
    simulation = None
    if strategy is not None:
        simulation = _initialize_simulation(
            reference_model, strategy, kwargs.get("solver_options"), verbose
        )
        _log_msg(LOGGER, logging.INFO, verbose, "Creating Simulation")
        if strategy not in STEADY_STATE_SOLVERS and strategy != "simulate":
            raise ValueError("Invalid steady state strategy: '{0}'".format(strategy))

        if perturbations is not None:
            if isinstance(perturbations, dict):
                perturbations = [perturbations]

            for i, perturbation in enumerate(perturbations):
                # Parse the perturbations and format the input to be used
                perturbations[i] = simulation._format_perturbations_input(
                    perturbation, verbose
                )

    if perturbations is None:
        perturbations = []

    # Dictionary to track when models failed a check
    numbers = {}
    if ensure_positive_percs is not None:
        numbers.update({"Infeasible, negative PERCs": 0})

    feasible_models = []
    infeasible_models = []
    for i, flux_values in enumerate(flux_data.values):
        # Check if flux data provided
        if flux_ids.size == 0:
            # No flux data provided
            f_suffix = ""
            flux_values = {}
        else:
            f_suffix = kwargs.get("flux_suffix") + str(i)
            flux_values = dict(zip(flux_ids, flux_values))

        # Iterate through concentration data
        for j, conc_values in enumerate(conc_data.values):
            # Check if conc data provided
            if conc_ids.size == 0:
                # No conc data provided
                c_suffix = ""
                conc_values = {}
            else:
                c_suffix = kwargs.get("conc_suffix") + str(j)
                conc_values = dict(zip(conc_ids, conc_values))

            # Make a copy of the model
            model = reference_model.copy()
            model.id += f_suffix + c_suffix
            _log_msg(LOGGER, logging.INFO, verbose, "New model '%s' created", model.id)

            # Update the model steady state flux values
            if flux_values:
                model.update_parameters(flux_values, verbose=False)
                _log_msg(
                    LOGGER,
                    logging.INFO,
                    verbose,
                    "Updated flux values for '%s'",
                    model.id,
                )

            # Update model concentration values
            if conc_values:
                model.update_initial_conditions(conc_values, verbose=False)
                _log_msg(
                    LOGGER,
                    logging.INFO,
                    verbose,
                    "Updated initial conditions for '%s'",
                    model.id,
                )

            if ensure_positive_percs is None:
                feasible_models.append(model)
                continue

            # Ensure PERCs are positive, updating model if they are
            model, is_feasible = _ensure_positive_percs_for_model(
                model,
                ensure_positive_percs,
                verbose,
                False,
                True,
                kwargs.get("at_equilibrium_default"),
            )

            if is_feasible:
                # Add feasible model to ensemble
                feasible_models.append(model)
            else:
                # Add infeasible model to ensemble if specified
                infeasible_models.append(model)
                numbers["Infeasible, negative PERCs"] += 1

    # Ensure steady state exists if given a strategy
    if strategy is not None:
        # Add models to simulation, all models are copies of the original
        # and can be directly added for performance
        simulation._simulation_values += [
            _get_sim_values_from_model(model) for model in feasible_models
        ]

        _log_msg(
            LOGGER,
            logging.INFO,
            verbose,
            "Successfully loaded feasible models into Simulation.",
        )

        # Define helper function for determining steady state feasibility
        def ensure_steady_state_success(
            check_type, perturbation=None, update_values=False
        ):
            """Function to determine whether a steady state can be found."""
            models = [m for m in feasible_models if m.id in simulation.models]
            with warnings.catch_warnings():
                conc_sol_list, flux_sol_list = simulation.find_steady_state(
                    models,
                    strategy=strategy,
                    perturbations=perturbation,
                    update_values=update_values,
                    verbose=verbose,
                    decimal_precision=kwargs.get("decimal_precision"),
                )

            # Add failed count to final report
            numbers[check_type] = 0
            for i, model in enumerate(models):
                if len(models) == 1:
                    conc_sol, flux_sol = conc_sol_list, flux_sol_list
                else:
                    conc_sol, flux_sol = conc_sol_list[i], flux_sol_list[i]
                if conc_sol and flux_sol:
                    _log_msg(
                        LOGGER,
                        logging.INFO,
                        verbose,
                        "Successfully found steady state for '%s'",
                        model,
                    )
                    if update_values:
                        ics, params = simulation.get_model_simulation_values(model)
                        model.update_initial_conditions(ics)
                        model.update_parameters(
                            {
                                param: value
                                for param, value in params.items()
                                if param in model.reactions.list_attr("flux_symbol_str")
                            }
                        )
                        model.update_parameters(params)
                    continue

                # Remove from feasible models and add to infeasible models
                _log_msg(
                    LOGGER,
                    logging.INFO,
                    verbose,
                    "No steady state found for '%s', " "adding to infeasible models",
                    model,
                )
                simulation.remove_models([model], verbose=verbose)
                numbers[check_type] += 1

            return

        # Ensure models can reach a steady state
        ensure_steady_state_success(
            "Infeasible, no steady state found", update_values=True
        )
        # Ensure models can reach a steady state with the given perturbations
        for i, perturbation in enumerate(perturbations):
            perturbation_str = "pertubration " + str(i + 1)
            _log_msg(
                LOGGER,
                logging.INFO,
                verbose,
                "Attempting to find steady state with %s",
                perturbation_str,
            )
            ensure_steady_state_success(
                "Infeasible, no steady state with " + perturbation_str, perturbation
            )

    if simulation is not None:
        infeasible_models += [
            model for model in feasible_models if model.id not in simulation.models
        ]
        feasible_models = [
            model for model in feasible_models if model.id in simulation.models
        ]

    # Format numbers to display
    num_str = "\n" if verbose else ""
    num_str += "Total models generated: {0}".format(
        len(feasible_models) + len(infeasible_models)
    )
    if numbers:
        num_str += "\nFeasible: {0}\n".format(str(len(feasible_models)))
        num_str += "\n".join(["{0}: {1}".format(k, v) for k, v in iteritems(numbers)])

    _log_msg(LOGGER, logging.INFO, True, num_str)
    if kwargs.get("return_infeasible"):

        return feasible_models, infeasible_models

    return feasible_models


def _validate_data_input(reference_model, data, id_type, verbose):
    """Validate the DataFrame input.

    Warnings
    --------
    This method is intended for internal use only.

    """
    if isinstance(data, pd.Series):
        data = pd.DataFrame(data).T
    if not isinstance(data, pd.DataFrame):
        raise TypeError("`data` must be a pandas DataFrame or Series")

    values = data.columns.values
    # Ensure all DataFrame columns are identifier strings
    if not all([isinstance(v, string_types) for v in values]):
        values = [getattr(v, "_id", v) for v in values]

    obj_list = {
        "metabolites": reference_model.metabolites,
        "reactions": reference_model.reactions,
    }[id_type]

    if not all([obj_id in obj_list for obj_id in values]):
        msg = "Invalid {0}".format(id_type)
        if verbose:
            msg += " {0!r}".format(
                [getattr(obj, "_id", obj) for obj in values if obj not in obj_list]
            )
        msg += " in data columns."
        raise ValueError(msg)

    data.columns = values
    if id_type == "reactions":
        values = np.array(["v_" + rid for rid in values])

    return data, values


def _ensure_positive_percs_for_model(
    model, reactions, verbose, raise_error, update_values, at_equilibrium_default
):
    """Ensure calculated PERCs for a model are positive.

    Warnings
    --------
    This method is intended for internal use only.

    """
    if reactions is None:
        reactions = model.reactions

    reactions = [getattr(r, "_id", r) for r in reactions]
    try:
        _log_msg(LOGGER, logging.INFO, verbose, "Calculating PERCs for '%s'", model.id)
        # Calculate PERCs
        percs = model.calculate_PERCs(
            at_equilibrium_default,
            fluxes={
                r: v
                for r, v in iteritems(model.steady_state_fluxes)
                if r.id in reactions
            },
        )
        negative_percs = [kf for kf, v in iteritems(percs) if v < 0]
        if negative_percs:
            # Found negative percs
            _log_msg(
                LOGGER,
                logging.WARN,
                verbose,
                "Negative PERCs '%s' calculated for '%s'",
                str(negative_percs),
                model.id,
            )
            is_feasible = False
        else:
            # No negative PERCs
            _log_msg(
                LOGGER,
                logging.INFO,
                verbose,
                "All PERCs are positive for '%s'",
                model.id,
            )
            if update_values:
                _log_msg(
                    LOGGER, logging.INFO, verbose, "Updating PERCs for '%s'", model.id
                )
                model.update_parameters(percs, verbose=False)
            is_feasible = True

    except ValueError as e:
        msg = str(
            "Could not calculate PERCs for '{0}' due to the "
            "following error: {1!r}".format(model.id, str(e))
        )
        if raise_error:
            raise MassEnsembleError(msg)
        _log_msg(LOGGER, logging.ERROR, verbose, msg)

    return model, is_feasible


def _initialize_simulation(reference_model, strategy, solver_options, verbose):
    """Initialize the Simulation object used in ensemble creation.

    Warnings
    --------
    This method is intended for internal use only.

    """
    simulation = Simulation(reference_model=reference_model, verbose=verbose)

    if solver_options:
        if strategy == "simulate":
            solver = simulation.integrator
        else:
            solver = simulation.steady_state_solver
        for option, value in iteritems(solver_options):
            if option in solver.getSettings():
                solver.setSetting(option, value)
            else:
                warnings.warn(
                    "Unrecognized option '{0}' for '{1}' solver.".format(
                        option, solver.getName()
                    )
                )
    return simulation


__all__ = (
    "generate_ensemble_of_models",
    "create_models_from_flux_data",
    "create_models_from_concentration_data",
    "ensure_positive_percs",
    "ensure_steady_state",
)
