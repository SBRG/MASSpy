# -*- coding: utf-8 -*-
"""TODO DOCSTRING."""
import logging
from warnings import warn

from cobra.core.dictlist import DictList

import numpy as np

import pandas as pd

from six import iteritems, string_types

from mass.core.mass_model import MassModel
from mass.core.simulation import Simulation, _log_msg
from mass.exceptions import MassEnsembleError, MassSimulationError
from mass.util.util import _check_kwargs, _make_logger, ensure_iterable

# Set the logger
LOGGER = _make_logger(__name__)
"""logging.Logger: Logger for :mod:`~mass.thermo.ensemble` submodule."""


class Ensemble(Simulation):
    """TODO DOCSTRING."""

    def __init__(self, reference_model, id=None, **kwargs):
        """Initialize the Ensemble."""
        kwargs = _check_kwargs({
            "name": None,
            "verbose": False,
        }, kwargs)

        # Create ID if None
        if id is None:
            id = "{0}_Ensemble".format(str(reference_model))

        # Intiailize
        try:
            super(Ensemble, self).__init__(reference_model=reference_model,
                                           id=id, name=kwargs.get("name"),
                                           verbose=kwargs.get("verbose"))
        except MassSimulationError as e:
            raise MassEnsembleError(e)

    def create_models_from_flux_data(self, models=None, data=None,
                                     raise_error=False, **kwargs):
        """TODO DOCSTRING."""
        kwargs = _check_kwargs({
            "verbose": False,
            "suffix": "_F",
        }, kwargs)
        verbose = kwargs.pop("verbose")
        if models is None:
            models = [self.reference_model]
        else:
            models = ensure_iterable(models)

        if any([not isinstance(model, MassModel) for model in models]):
            raise TypeError("`models` must be an iterable of MassModels.")

        # Validate data input, ensuring all columns are identifiers
        data, id_array = _validate_data_input(self.reference_model, data,
                                              "reactions", verbose)
        ensemble = np.empty((len(models), data.shape[0])).astype(np.object_)
        # Loop through models
        for i, model in enumerate(models):
            try:
                # Loop through value sets.
                for j, values in enumerate(data.values):
                    # Create new model
                    new_model = model.copy()
                    new_model.id += kwargs.get("suffix") + str(j)
                    _log_msg(LOGGER, logging.INFO, verbose,
                             "New model '%s' created", new_model.id)

                    # Update the model parameters
                    new_model.update_parameters(dict(zip(id_array, values)),
                                                verbose=False)
                    _log_msg(LOGGER, logging.INFO, verbose,
                             "Updated flux values for '%s'", new_model.id)

                    # Add model to the ensemble
                    ensemble[i, j] = new_model

            except Exception as e:
                msg = str("Could not use '{0}' for ensemble creation due to "
                          "the following error: {1!r}".format(
                              model.id, str(e)))
                if raise_error:
                    raise MassEnsembleError(msg)
                _log_msg(LOGGER, logging.ERROR, verbose, msg)

        ensemble = [m for m in ensemble.flatten() if isinstance(m, MassModel)]

        return ensemble

    def create_models_from_concentration_data(self, models=None, data=None,
                                              raise_error=False, **kwargs):
        """TODO DOCSTRING."""
        kwargs = _check_kwargs({
            "verbose": False,
            "suffix": "_C",
        }, kwargs)
        verbose = kwargs.pop("verbose")
        if models is None:
            models = [self.reference_model]
        else:
            models = ensure_iterable(models)

        if any([not isinstance(model, MassModel) for model in models]):
            raise TypeError("`models` must be an iterable of MassModels.")

        # Validate data input, ensuring all columns are identifiers
        data, id_array = _validate_data_input(self.reference_model, data,
                                              "metabolites", verbose)

        ensemble = np.empty((len(models), data.shape[0])).astype(np.object_)
        # Loop through models
        for i, model in enumerate(models):
            try:
                # Loop through value sets.
                for j, values in enumerate(data.values):
                    # Create new model
                    new_model = model.copy()
                    new_model.id += kwargs.get("suffix") + str(j)
                    _log_msg(LOGGER, logging.INFO, verbose,
                             "New model '%s' created", new_model.id)

                    # Update the model parameters
                    new_model.update_initial_conditions(
                        dict(zip(id_array, values)), verbose=False)
                    _log_msg(LOGGER, logging.INFO, verbose,
                             "Updated initial conditions for '%s'", 
                             new_model.id)

                    # Add model to the ensemble
                    ensemble[i, j] = new_model

            except Exception as e:
                msg = str("Could not use '{0}' for ensemble creation due to "
                          "the following error: {1!r}".format(
                              model.id, str(e)))
                if raise_error:
                    raise MassEnsembleError(msg)
                _log_msg(LOGGER, logging.ERROR, verbose, msg)

        ensemble = [m for m in ensemble.flatten() if isinstance(m, MassModel)]

        return ensemble

    def ensure_positive_percs(self, models, raise_error=False,
                              update_values=False, **kwargs):
        """TODO DOCSTRING."""
        kwargs = _check_kwargs({
            "verbose": False,
            "at_equilibrium_default": 100000,
        }, kwargs)
        verbose = kwargs.pop("verbose")
        feasible = []
        infeasible = []

        models = ensure_iterable(models)
        if not models:
            warn("No models provided.")
            return feasible, infeasible

        if any([not isinstance(model, MassModel) for model in models]):
            raise TypeError("`models` must be an iterable of MassModels.")

        for model in models:
            model, is_feasible = _ensure_positive_percs_for_model(
                model, verbose, raise_error, update_values,
                kwargs.get("at_equilibrium_default"))
            if is_feasible:
                feasible.append(model)
            else:
                infeasible.append(model)

        _log_msg(LOGGER, logging.INFO, verbose,
                 "Finished PERC calculations, returning seperated models.")
        return feasible, infeasible

    def ensure_steady_state(self, models=None, strategy="nleq2",
                            perturbations=None, update_values=False,
                            **kwargs):
        """TODO DOCSTRING."""
        kwargs = _check_kwargs({
            "verbose": False,
            "selections": None,
            "boundary_metabolites": False,
            "steps": None,
            "tfinal": 1e8,
            "num_attempts": 2,
            "decimal_precision": True}, kwargs)
        # Get models, adding any models that may not already exist
        if models is None:
            models = self.models 
        elif any([str(m) not in self.models for m in models]):
            self.add_models([m for m in models if m not in self.models],
                            verbose=kwargs.get("verbose"))

        conc_sol_list, flux_sol_list = self.find_steady_state(
            models, strategy, perturbations, update_values, **kwargs)

        feasible = []
        infeasible = []
        for i, model in enumerate(models):
            conc_sol, flux_sol = conc_sol_list[i], flux_sol_list[i]
            if conc_sol and flux_sol:
                feasible.append(model)
            else:
                infeasible.append(model)

        _log_msg(LOGGER, logging.INFO, kwargs.get("verbose"),
                 "Finished finding steady states, returning seperated models.")

        return feasible, infeasible


def _validate_data_input(reference_model, data, id_type, verbose):
    """TODO DOCSTRING."""
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
        "reactions": reference_model.reactions
    }[id_type]

    if not all([obj_id in obj_list for obj_id in values]):
        msg = "Invalid {0}".format(id_type)
        if verbose:
            msg += " {0!r}".format([
                getattr(obj, "_id", obj) for obj in values
                if obj not in obj_list])
        msg += " in data columns."
        raise ValueError(msg)

    data.columns = values
    if id_type == "reactions":
        values = np.array(["v_" + rid for rid in values])

    return data, values


def _ensure_positive_percs_for_model(model, verbose, raise_error,
                                     update_values, at_equilibrium_default):
    """TODO DOCSTRING."""
    try:
        _log_msg(LOGGER, logging.INFO, verbose,
                 "Calculating PERCs for '%s'", model.id)
        # Calculate PERCs
        percs = model.calculate_PERCs(at_equilibrium_default)
        negative_percs = [kf for kf, v in iteritems(percs) if v < 0]
        if negative_percs:
            # Found negative percs
            _log_msg(LOGGER, logging.WARN, verbose,
                     "Negative PERCs '%s' calculated for '%s'",
                     str(negative_percs), model.id)
            is_feasible = False
        else:
            # No negative PERCs
            _log_msg(LOGGER, logging.INFO, verbose,
                     "All PERCs are positive for '%s'", model.id)
            if update_values:
                _log_msg(LOGGER, logging.INFO, verbose,
                         "Updating PERCs for '%s'", model.id)
                model.update_parameters(percs, verbose=False)
            is_feasible = True

    except ValueError as e:
        msg = str("Could not calculate PERCs for '{0}' due to the "
                  "following error: {1!r}".format(model.id, str(e)))
        if raise_error:
            raise MassEnsembleError(msg)
        _log_msg(LOGGER, logging.ERROR, verbose, msg)

    return model, is_feasible


def generate_ensemble(reference_model, flux_data=None, conc_data=None,
                      steady_state_strategy="simulate", perturbations=None,
                      **kwargs):
    """TODO DOCSTRING.

    Notes
    -----
    * This function is optimized for performance when generating a large
      ensemble of models when compared to the various methods of the
      :class:`Ensemble` class. However, this function may not provide as much
      control over the process as the methods defined in the :class:`Ensemble`
      class.

    """
    # Check all inputs at beginning to ensure that ensemble generation is not
    # disrupted near the end due to invalid input format
    kwargs = _check_kwargs({
        "verbose": False,
        "ensure_positive_percs": True,
        "flux_suffix": "_F",
        "conc_suffix": "_C",
        "at_equilibrium_default": 100000,
        "return_infeasible": False,
    }, kwargs)
    verbose = kwargs.pop("verbose")
    ensure_positive_percs = kwargs.pop("ensure_positive_percs")

    # Validate model input
    if not isinstance(reference_model, MassModel):
        raise TypeError("`reference_model` must be a MassModel.")

    # Validate DataFrame inputs, if any
    if flux_data is not None:
        # Validate flux data if provided
        flux_data, flux_ids = _validate_data_input(reference_model, flux_data,
                                                   "reactions", verbose)
    else:
        # Set a value to allow for iteration
        flux_data = pd.DataFrame([0])
        flux_ids = np.array([])

    if conc_data is not None:
        # Validate conc data if provided
        conc_data, conc_ids = _validate_data_input(reference_model, conc_data,
                                                   "metabolites", verbose)
    else:
        # Set a value to allow for iteration
        conc_data = pd.DataFrame([0])
        conc_ids = np.array([])

    if conc_ids.size == 0 and flux_ids.size == 0:
        raise ValueError("No flux data or concentration data provided. "
                         "At least one data set must be provided to generate "
                         "the ensemble.")

    # Generate ensemble for feasible models
    # Loading model into ensemble validates whether it can be simulated and
    # loaded into roadrunner. Allows for strategy and perturbation checking.
    feasible = Ensemble(reference_model)
    _log_msg(LOGGER, logging.INFO, verbose, "Creating feasible Ensemble")

    # Ensure strategy input and perturbation input is valid
    if steady_state_strategy is not None:
        valid = list(feasible.roadrunner.getRegisteredSteadyStateSolverNames())
        valid += ["simulate"]
        if steady_state_strategy not in valid:
            raise ValueError("Invalid steady state strategy")

        if perturbations is not None:
            if isinstance(perturbations, dict):
                perturbations = [perturbations]

            for i, perturbation in enumerate(perturbations):
                # Parse the perturbations and format the input to be used
                perturbations[i] = feasible._format_perturbations_input(
                    perturbation, verbose)
    else:
        perturbations = []

    # Dictionary to track when models were determined to be infeasible
    numbers = {"Infeasible, negative PERCs": 0}

    feasible_list = DictList([])
    infeasible_list = DictList([])
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
            _log_msg(LOGGER, logging.INFO, verbose,
                     "New model '%s' created", model.id)

            # Update the model steady state flux values
            if flux_values:
                model.update_parameters(flux_values, verbose=False)
                _log_msg(LOGGER, logging.INFO, verbose,
                         "Updated flux values for '%s'", model.id)

            # Update model concentration values
            if conc_values:
                model.update_initial_conditions(conc_values, verbose=False)
                _log_msg(LOGGER, logging.INFO, verbose,
                         "Updated initial conditions for '%s'", model.id)

            # Ensure PERCs are positive, updating model if they are
            if ensure_positive_percs:
                model, is_feasible = _ensure_positive_percs_for_model(
                    model, verbose, False, True,
                    kwargs.get("at_equilibrium_default"))

                if is_feasible:
                    # Add feasible model to ensemble
                    feasible_list.append(model)
                else:
                    # Add infeasible model to ensemble if specified
                    infeasible_list.append(model)
                    numbers["Infeasible, negative PERCs"] += 1

    # Add feasible models to ensemble
    feasible.add_models(feasible_list)

    # Ensure steady state exists if given a strategy
    if steady_state_strategy is not None:

        # Define helper function for determining steady state feasibility
        def ensure_steady_state_feasibility(num_key, perturbation=None):
            models = [m for m in feasible.models if m != reference_model.id]
            conc_sol_list, flux_sol_list = feasible.find_steady_state(
                models=models, strategy=steady_state_strategy,
                perturbations=perturbation,
                update_values=True, verbose=verbose)
            numbers[num_key] = 0
            for i, model in enumerate(models):
                conc_sol, flux_sol = conc_sol_list[i], flux_sol_list[i]
                if conc_sol and flux_sol:
                    _log_msg(LOGGER, logging.INFO, verbose,
                             "Successfully found steady state for '%s'",
                             model)
                    continue

                # Remove from feasible models and add to infeasible models
                _log_msg(LOGGER, logging.INFO, verbose,
                         "No steady state found for '%s', "
                         "adding to infeasible", model)
                feasible.remove_models([model], verbose=verbose)
                model = feasible_list.pop(feasible_list.index(model))
                infeasible_list.append(model)
                numbers[num_key] += 1

        # Ensure models can reach a steady state
        ensure_steady_state_feasibility("Infeasible, no steady state found")

        # Ensure models can reach a steady state with the given perturbations
        for i, perturbation in enumerate(perturbations):
            perturbation_str = "pertubration " + str(i + 1)
            _log_msg(LOGGER, logging.INFO, verbose,
                     "Attempting to find steady state with %s",
                     perturbation_str)
            ensure_steady_state_feasibility(
                "Infeasible, no steady state with " + perturbation_str,
                perturbation)

    num_str = "\nTotal Models generated: {0}\n".format(
        len(feasible.models) - 1 + len(infeasible_list))
    num_str += "Feasible: {0}\n".format(str(len(feasible.models) - 1))
    num_str += "\n".join(["{0}: {1}".format(k, v)
                         for k, v in iteritems(numbers)])

    if not kwargs.get("return_infeasible"):
        _log_msg(LOGGER, logging.INFO, verbose, num_str)

        return feasible

    # Alter the ID to make it clear which is the feasible Ensemble
    feasible.id += "_Feasible"

    # Create ensemble for infeasible models
    _log_msg(LOGGER, logging.INFO, verbose, "Creating infeasible Ensemble")
    infeasible = Ensemble(reference_model)
    infeasible.id += "_Infeasible"
    infeasible.add_models(infeasible_list)

    _log_msg(LOGGER, logging.INFO, verbose, num_str)

    return feasible, infeasible


__all__ = ("Ensemble", "generate_ensemble",)
