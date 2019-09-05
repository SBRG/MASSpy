# -*- coding: utf-8 -*-
"""TODO DOCSTRING."""
import logging
from warnings import warn

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
        data, id_list = self._validate_data_input(data, "reactions",
                                                  verbose)

        suffix = kwargs.get("suffix")

        ensemble = np.empty((len(models), data.shape[0])).astype(np.object_)
        # Loop through models
        for i, model in enumerate(models):
            try:
                # Loop through value sets.
                for j, values in enumerate(data.values):
                    # Create new model
                    new_model = model.copy()
                    new_model.id += suffix + str(j)
                    _log_msg(LOGGER, logging.INFO, verbose,
                             "New model '%s' created", new_model.id)

                    # Update the model parameters
                    new_model.update_parameters(dict(zip(id_list, values)),
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
        data, id_list = self._validate_data_input(data, "metabolites",
                                                  verbose)

        suffix = kwargs.get("suffix")

        ensemble = np.empty((len(models), data.shape[0])).astype(np.object_)
        # Loop through models
        for i, model in enumerate(models):
            try:
                # Loop through value sets.
                for j, values in enumerate(data.values):
                    # Create new model
                    new_model = model.copy()
                    new_model.id += suffix + str(j)
                    _log_msg(LOGGER, logging.INFO, verbose,
                             "New model '%s' created", new_model.id)

                    # Update the model parameters
                    new_model.update_initial_conditions(
                        dict(zip(id_list, values)), verbose=False)
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
            try:
                _log_msg(LOGGER, logging.INFO, verbose,
                         "Calculating PERCs for '%s'", model.id)
                percs = model.calculate_PERCs(**kwargs)
                negative_percs = [kf for kf, v in iteritems(percs) if v < 0]
                if not negative_percs:
                    _log_msg(LOGGER, logging.INFO, verbose,
                             "All PERCs are positive for '%s'", model.id)
                    if update_values:
                        _log_msg(LOGGER, logging.INFO, verbose,
                                 "Updating PERCs for '%s'", model.id)
                        model.update_parameters(percs, verbose=False)
                    feasible.append(model)
                else:
                    _log_msg(LOGGER, logging.WARN, verbose,
                             "Negative PERCs '%s' calculated for '%s'",
                             str(negative_percs), model.id)

                    infeasible.append(model)
            except ValueError as e:
                msg = str("Could not calculate PERCs for '{0}' due to the "
                          "following error: {1!r}".format(model.id, str(e)))
                if raise_error:
                    raise MassEnsembleError(msg)
                _log_msg(LOGGER, logging.ERROR, verbose, msg)

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

    def _validate_data_input(self, data, id_type, verbose):
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
            "metabolites": self.reference_model.metabolites,
            "reactions": self.reference_model.reactions
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
        return data, values


def generate_ensemble(reference_model, flux_data=None, concentration_data=None,
                      steady_state_strategy="nleq2", perturbation_list=None):
    """TODO DOCSTRING.

    Notes
    -----
    * This function is optimized for performance when generating an
      ensemble of models when compared to the various methods of the
      :class:`Ensemble` class. However, this function may not provide as much
      control over the process as the methods defined in the :class:`Ensemble`
      class.

    """


__all__ = ("Ensemble", "generate_ensemble",)
