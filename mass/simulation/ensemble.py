# -*- coding: utf-8 -*-
r"""Module to create and manage an ensemble of models.

The :mod:`~.ensemble` module contains the :class:`Ensemble` class and the
:func:`generate_ensemble` function.

The :class:`Ensemble` class is designed to manage an ensemble of models. It
contains various methods to assist in generating multiple models from existing
:class:`~.MassModel`\ s from either flux data or concentration data in
:class:`pandas.DataFrame`\ s  (e.g. generated from
:mod:`~mass.thermo.conc_sampling`).
There are also methods to help ensure that models are thermodynamically
feasible and can reach steady states with or without perturbations applied.

The benefits of using the the :class:`Ensemble` is that it provides more
control over the model creation and validation processes. However, if a large
number of models is to be generated and performance is a priority, the
:func:`generate_ensemble` function can be used.

The :func:`generate_ensemble` function is optimized for performance when
generating a large number of models. This function ensures valid user input is
before generating models to reduce the likelihood of a user error causing the
model generation process to stop before completion. However, there is time
spent in function's setup, so for generation of a smaller number of models,
performance gains may not be seen.

"""
import logging
import warnings

from cobra.core.dictlist import DictList

import numpy as np

import pandas as pd

from six import iteritems, string_types

from mass.core.mass_model import MassModel
from mass.exceptions import MassEnsembleError, MassSimulationError
from mass.simulation.simulation import Simulation, _log_msg
from mass.util.util import _check_kwargs, _make_logger, ensure_iterable

# Set the logger
LOGGER = _make_logger(__name__)
"""logging.Logger: Logger for :mod:`~mass.thermo.ensemble` submodule."""


class Ensemble(Simulation):
    """Class for managing a large ensemble of :mod:`mass` models.

    Parameters
    ----------
    reference_model : MassModel
        The model to load for simulation. The model will be set as the
        :attr:`Simulation.reference_model`.
    id : str or None
        An identifier to associate with the :class:`Simulation`. If ``None``
        then one is automatically created based on the model identifier.
    name : str
        A human readable name for the :class:`Simulation`.
    verbose : bool
        Whether to provide a QCQA report and more verbose messages when trying
        to load the model. Default is ``False``.
    *kwargs
        variable_step_size :
            ``bool`` indicating whether to initialize the integrator with a
            variable time step for simulations.

            Default is ``False``.

    """

    def __init__(self, reference_model, id=None, name=None, verbose=False,
                 **kwargs):
        """Initialize the Ensemble."""
        # Create ID if None
        if id is None:
            id = "{0}_Ensemble".format(str(reference_model))

        # Intiailize
        try:
            super(Ensemble, self).__init__(reference_model=reference_model,
                                           id=id, name=name, verbose=verbose,
                                           **kwargs)
        except MassSimulationError as e:
            raise MassEnsembleError(e)

    def create_models_from_flux_data(self, models=None, data=None,
                                     raise_error=False, **kwargs):
        """Generate ensemble of models for a given set of flux data.

        Notes
        -----
        * If ``x`` models are provided with ``y`` sets of flux data values,
          ``x * y`` total models will be generated.

        Parameters
        ----------
        models : iterable, None
            An iterable of :class:`.MassModel` objects to treat as references.
            If ``None`` provided, the :attr:`Ensemble.reference_model` is used.
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
        ensemble : list
            A ``list`` of successfully generated :class:`.MassModel` objects.

        Raises
        ------
        MassEnsembleError
            Raised if generation of a model fails and ``raise_error=True``.

        """
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
        """Generate ensemble of models for a given set of concentration data.

        Notes
        -----
        * If ``x`` models are provided with ``y`` sets of concentration data
          values, ``x * y`` total models will be generated.

        Parameters
        ----------
        models : iterable, None
            An iterable of :class:`.MassModel` objects to treat as references.
            If ``None`` provided, the :attr:`Ensemble.reference_model` is used.
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
        ensemble : list
            A ``list`` of successfully generated :class:`.MassModel` objects.

        Raises
        ------
        MassEnsembleError
            Raised if generation of a model fails and ``raise_error=True``.

        """
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

    def ensure_positive_percs(self, models, reactions=None, raise_error=False,
                              update_values=False, **kwargs):
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
        kwargs = _check_kwargs({
            "verbose": False,
            "at_equilibrium_default": 100000,
        }, kwargs)
        verbose = kwargs.pop("verbose")
        positive = []
        negative = []

        models = ensure_iterable(models)
        if not models:
            warnings.warn("No models provided.")
            return positive, negative

        if any([not isinstance(model, MassModel) for model in models]):
            raise TypeError("`models` must be an iterable of MassModels.")

        for model in models:
            model, is_positive = _ensure_positive_percs_for_model(
                model, reactions, verbose, raise_error, update_values,
                kwargs.get("at_equilibrium_default"))
            if is_positive:
                positive.append(model)
            else:
                negative.append(model)

        _log_msg(LOGGER, logging.INFO, verbose,
                 "Finished PERC calculations, returning seperated models.")
        return positive, negative

    def ensure_steady_state(self, models=None, strategy="nleq2",
                            perturbations=None, update_values=False,
                            **kwargs):
        """Seperate models based on whether a steady state can be reached.

        All ``kwargs`` are passed to :meth:`~.Simulation.find_steady_state`.

        Parameters
        ----------
        models : iterable, None
            An iterable of :class:`.MassModel` objects to find a steady state
            for. If ``None`` provided, all models in the :class:`Ensemble` are
            used.
        strategy : str
            The strategy for finding the steady state. Must be one of the
            following:

                * ``'simulate'``
                * ``'nleq1'``
                * ``'nleq2'``

        perturbations : dict
            A ``dict`` of perturbations to incorporate into the simulation.
            Models must reach a steady state with the given pertubration to be
            considered feasible.
            See :mod:`~.simulation.simulation` documentation for more
            information on valid perturbations.
        update_values : bool
            Whether to update the model with the steady state results.
            Default is ``False``.
        **kwargs
            verbose :
                ``bool`` indicating the verbosity of the method.

                Default is ``False``.
            selections :
                ``list`` of identifiers corresponding to the time course
                selections to return in the solutions. If pools or net fluxes
                are included, all variables for their formulas must be
                included.

                Default is ``None`` for all concentrations, fluxes, pools,
                and net fluxes.
            boundary_metabolites :
                ``bool`` indicating whether to include boundary metabolite
                concentrations in the output.

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
            A ``list`` of :class:`.MassModel` objects that could reach
            a steady state.
        infeasible : list
            A ``list`` of :class:`.MassModel` objects that could not reach
            a steady state.

        """
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


def _ensure_positive_percs_for_model(model, reactions, verbose, raise_error,
                                     update_values, at_equilibrium_default):
    """Ensure calculated PERCs for a model are positive.

    Warnings
    --------
    This method is intended for internal use only.

    """
    if reactions is None:
        reactions = model.reactions

    reactions = [getattr(r, "_id", r) for r in reactions]
    try:
        _log_msg(LOGGER, logging.INFO, verbose,
                 "Calculating PERCs for '%s'", model.id)
        # Calculate PERCs
        percs = model.calculate_PERCs(at_equilibrium_default, fluxes={
            r: v for r, v in iteritems(model.steady_state_fluxes)
            if r.id in reactions})
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
                      ensure_positive_percs=None, steady_state_strategy=None,
                      perturbations=None, **kwargs):
    """Generate an :class:`Ensemble` of feasible models for given data sets.

    This function is optimized for performance when generating a large
    ensemble of models when compared to the various methods of the
    :class:`Ensemble` class. However, this function may not provide as much
    control over the process as the methods defined in the :class:`Ensemble`
    class.

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
        The model will be set as the :attr:`Ensemble.reference_model`.
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
    steady_state_strategy : str, None
        The strategy for finding the steady state. Must be one of the
        following:

            * ``'simulate'``
            * ``'nleq1'``
            * ``'nleq2'``

        If ``None``, no attempts will be made to determine whether a generated
        model can reach a steady state.
    perturbations : dict
        A ``dict`` of perturbations to incorporate into the simulation,
        or a list of perturbation dictionaries where each ``dict`` is applied
        to a simulation. Models must reach a steady state with all given
        pertubration dictionaries to be considered feasible.
        See :mod:`~.simulation.simulation` documentation for more
        information on valid perturbations.

        Ignored if ``steady_state_strategy=None``.
    **kwargs
        verbose :
            ``bool`` indicating the verbosity of the function.

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
    feasible : Ensemble
        An :class:`Ensemble` containing the models deemed feasible.
    infeasible : Ensemble
        An :class:`Ensemble` containing the models deemed infeasible.
        Only returned if ``return_infeasible=True``.

    """
    # Check all inputs at beginning to ensure that ensemble generation is not
    # disrupted near the end due to invalid input format
    kwargs = _check_kwargs({
        "verbose": False,
        "flux_suffix": "_F",
        "conc_suffix": "_C",
        "at_equilibrium_default": 100000,
        "return_infeasible": False,
    }, kwargs)
    verbose = kwargs.pop("verbose")
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

    if perturbations is None:
        perturbations = []

    # Dictionary to track when models were determined to be infeasible
    numbers = {}

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

            if ensure_positive_percs is None:
                feasible_list.append(model)
            else:
                numbers.update({"Infeasible, negative PERCs": 0})
                if ensure_positive_percs is None:
                    ensure_positive_percs = []
                # Ensure PERCs are positive, updating model if they are
                model, is_feasible = _ensure_positive_percs_for_model(
                    model, ensure_positive_percs, verbose, False, True,
                    kwargs.get("at_equilibrium_default"))

                if is_feasible:
                    # Add feasible model to ensemble
                    feasible_list.append(model)
                else:
                    # Add infeasible model to ensemble if specified
                    infeasible_list.append(model)
                    numbers["Infeasible, negative PERCs"] += 1

    # Add feasible models to ensemble
    feasible.add_models(feasible_list, verbose=verbose)

    # Ensure steady state exists if given a strategy
    if steady_state_strategy is not None:
        # Define helper function for determining steady state feasibility
        def ensure_steady_state_feasibility(num_key, perturbation=None,
                                            update_values=False):
            models = [m for m in feasible.models if m != reference_model.id]
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                conc_sol_list, flux_sol_list = feasible.find_steady_state(
                    models=models, strategy=steady_state_strategy,
                    perturbations=perturbation,
                    update_values=update_values, verbose=verbose)
            numbers[num_key] = 0
            for i, model in enumerate(models):
                if len(models) == 1:
                    conc_sol, flux_sol = conc_sol_list, flux_sol_list
                else:
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
        ensure_steady_state_feasibility("Infeasible, no steady state found",
                                        update_values=True)

        # Ensure models can reach a steady state with the given perturbations
        for i, perturbation in enumerate(perturbations):
            perturbation_str = "pertubration " + str(i + 1)
            _log_msg(LOGGER, logging.INFO, verbose,
                     "Attempting to find steady state with %s",
                     perturbation_str)
            ensure_steady_state_feasibility(
                "Infeasible, no steady state with " + perturbation_str,
                perturbation)
    num_str = "\n" if verbose else ""
    num_str += "Total models generated: {0}".format(
        len(feasible.models) - 1 + len(infeasible_list))
    if numbers:
        num_str += "\nFeasible: {0}\n".format(str(len(feasible.models) - 1))
        num_str += "\n".join(["{0}: {1}".format(k, v)
                              for k, v in iteritems(numbers)])

    if not kwargs.get("return_infeasible"):
        _log_msg(LOGGER, logging.INFO, True, num_str)

        return feasible

    # Alter the ID to make it clear which is the feasible Ensemble
    feasible.id += "_Feasible"

    # Create ensemble for infeasible models
    _log_msg(LOGGER, logging.INFO, verbose, "Creating infeasible Ensemble")
    infeasible = Ensemble(reference_model)
    infeasible.id += "_Infeasible"
    infeasible.add_models(infeasible_list)

    _log_msg(LOGGER, logging.INFO, True, num_str)

    return feasible, infeasible


__all__ = ("Ensemble", "generate_ensemble",)
