# -*- coding: utf-8 -*-
r"""The Simulation module addresses the simulation of :class:`~.MassModel`\ s.

The :class:`~.Simulation` is designed to address all aspects related to the
simulation of one or more :class:`~.MassModel`\ s. These aspects include
initializing and compiling the model into a :class:`roadrunner.RoadRunner`
instance (the ODE integrator), storing the numerical values necessary for the
simulations, and handling the simulation results by creating and storing
:class:`~.MassSolution`\ s.

A :class:`~.Simulation` is initialized by providing a simulatable
:class:`~.MassModel` that can be converted into an SBML compliant model.

Multiple models can be added to the :class:`Simualtion` in order to simulate
an ensemble of models. To add an additional model to the ensemble, the model
must meet three criteria:

    1. The model must have equivalent ODEs to the
       :attr:`~Simulation.reference_model`
    2. The model must not have the same ID as the
       :attr:`~Simulation.reference_model`.
    3. The model must have all numerical values needed for simulation.

Perturbations can be implemented for a given simulation as long as they
follow the following guidelines:

    1. Perturbations are ``dict``\ s where the ``{key: value}`` pairs are the
       variables to be perturbed and the new numerical value or value change.
    2. To scale the current value of a variable, the value should be a ``str``
       representing the formula for altering perturbation variable, where the
       variable in the ``str`` is identical to the perturbation key.
    3. If providing a formula ``str`` as the perturbation value, it must be
       possible to 'sympify' the string using the
       :func:`~sympy.core.sympify.sympify` function. It must only have one
       variable, identical to the perturbation key.
    4. Only boundary conditions can be changed to have functions of time that
       represent the external concentration at that point in time.
       If a perturbation value is to be a string representing a function, it
       must be a function, it may only contain the time variable ``'t'`` and
       the boundary metabolite variable.

Some examples of perturbations following the guidelines for a model containing
the specie with ID ``'MID_c'``, boundary metabolite with ID ``'MID_b'``, and
reaction with ID ``'RID'``:

    * Altering :attr:`~.MassModel.initial_conditions` (ICs):

        * ``{'MID_c': 2}`` Change the IC value to 2.
        * ``{'MID_c': 'MID_c * 1.5'}`` Increase current IC value by 50%.

    * Altering :attr:`~.MassModel.parameters`:

        * ``{'kf_RID': 'kf_RID * 0.75'}`` Decrease kf parameter value by 25%.
        * ``{'Keq_RID': 100}`` Change Keq parameter value to 100.

    * Altering :attr:`~.MassModel.boundary_conditions` (BCs):

        * ``{'MID_b': 'sin(2 * pi * t)'}`` Change BC to a ``sin`` function.
        * ``{'MID_b': 'MID_b + cos(t)'}`` Add ``cos`` function to current BC
          value.

Note that perturbations using functions of time may take longer to implement
than other perturbations.

All simulation results are returned as :class:`.MassSolution`\ s.
Each simulated model has a corresponding :class:`.MassSolution`\ s stored
in the :class:`Simulation`. These solution objects are stored until being
replaced by a new :class:`~.MassSolution` upon resimulating the model. This
means that there can only be one concentration solution and one flux solution
per simulated model. A failed simulation of a model will return an empty
:class:`~.MassSolution`.

Though the :class:`Simulation` utilizes the :mod:`roadrunner` package, the
standard :mod:`logging` module will be used for :mod:`mass` logging purposes in
the :mod:`simulation` submodule. Therefore, the roadrunner logger is disabled
upon loading the :mod:`simulation` submodule. However, because the
:class:`Simulation` utilizes the :mod:`roadrunner` package for simulating
models, the :class:`roadrunner.Logger` can be accessed via the
:const:`RR_LOGGER` variable for those who wish to utilize it. See the
:mod:`roadrunner` documentation for more information on how to configure the
:class:`roadrunner.Logger`.
"""
import logging
from warnings import warn

from cobra.core.dictlist import DictList
from cobra.core.object import Object

from libsbml import writeSBMLToString

import numpy as np

import roadrunner

from six import iteritems

from sympy import Basic, Function, Symbol, sympify

from mass.core.mass_configuration import MassConfiguration
from mass.core.mass_model import MassModel
from mass.core.mass_solution import MassSolution, _CONC_STR, _FLUX_STR
from mass.exceptions import MassSimulationError
from mass.io.sbml import _model_to_sbml
from mass.util.dict_with_id import DictWithID
from mass.util.qcqa import is_simulatable, qcqa_model
from mass.util.util import (
    _check_kwargs, _make_logger, apply_decimal_precision, ensure_iterable)
# Set the logger
MASSCONFIGURATION = MassConfiguration()
# If working in Python application (e.g. iPython notebooks), enable logging

# Set the logger
LOGGER = _make_logger(__name__)
"""logging.Logger: Logger for :mod:`~mass.core.simulation` submodule."""

RR_LOGGER = roadrunner.Logger
"""roadrunner.Logger: The logger for the :mod:`roadrunner`."""
RR_LOGGER.disableLogging()

# SBML writing kwargs
_SBML_KWARGS = {"use_fbc_package": True, "use_groups_package": True,
                "units": False, "local_parameters": False,
                "write_objective": False}


class Simulation(Object):
    """Class for managing setup and result handling of simulations.

    The :class:`Simulation` class is designed to address all aspects related
    to the simulation of :class:`~.MassModel` objects, including setting the
    solver and solver options, perturbation of concentrations and parameters,
    simulation of a single model or an ensemble of models, and handling
    of simulation results.

    Parameters
    ----------
    reference_model : MassModel
        The model to load for simulation. The model will be set as the
        :attr:`Simulation.reference_model`.
    simulation_id : str or None
        An identifier to associate with the :class:`Simulation`. If ``None``
        then one is automatically created based on the model identifier.
    simulation_name : str
        A human readable name for the :class:`Simulation`.
    verbose : bool
        Whether to provide a QCQA report and more verbose messages when trying
        to load the model. Default is ``False``.

    """

    def __init__(self, reference_model, id=None, name=None, verbose=False):
        """Initialize the Simulation."""
        if not isinstance(reference_model, MassModel):
            raise TypeError(
                "'{0}' is not a valid MassModel instance".format(
                    str(reference_model)))

        if id is None:
            id = "{0}_Simulation".format(str(reference_model))

        try:
            # QCQA check model
            _assess_model_quality_for_simulation(
                reference_model, verbose, verbose)
            # Load model into RoadRunner
            rr = _load_model_into_roadrunner(reference_model, rr=None,
                                             verbose=verbose, **_SBML_KWARGS)
        except MassSimulationError as e:
            msg = "Could not load MassModel '{0}'".format(
                str(reference_model))
            if verbose:
                msg += ": " + str(e)
            raise MassSimulationError(msg)

        # Initialize Simulation
        super(Simulation, self).__init__(id, name)
        # Store the original model used to create the Simulation
        self._reference_model = reference_model
        # Set roadrunner
        self._roadrunner = rr

        # Storing model values for simulations
        self._simulation_values = DictList()
        self._simulation_values.add(
            _get_sim_values_from_model(reference_model))

        # Storing concentration and flux solutions in simulations
        self._concentration_solutions = DictList()
        self._flux_solutions = DictList()

    @property
    def reference_model(self):
        """Return the reference model of the :class:`Simulation`."""
        return getattr(self, "_reference_model")

    @property
    def models(self):
        """Return the IDs of models that exist in the :class:`Simulation`."""
        return [d.id.replace("_values", "") for d in self._simulation_values]

    @property
    def roadrunner(self):
        """Return the :class:`~roadrunner.RoadRunner` instance."""
        return self._roadrunner

    @property
    def concentration_solutions(self):
        r"""Get a copy of stored :class:`.MassSolution`\ s for concentrations.

        Returns
        -------
        ~cobra.core.dictlist.DictList
            Contains all :class:`.MassSolution` objects for concentrations.

        """
        return getattr(self, "_concentration_solutions").copy()

    @property
    def flux_solutions(self):
        r"""Get a copy of the stored :class:`.MassSolution`\ s for fluxes.

        Returns
        -------
        ~cobra.core.dictlist.DictList
            Contains all :class:`.MassSolution` objects for fluxes.

        """
        return getattr(self, "_flux_solutions").copy()

    def set_new_reference_model(self, model, verbose=False):
        """Set a new reference model for the :class:`Simulation`.

        To set a new reference model, the model must meet three criteria:

            1. The model must have equivalent ODEs to the
               :attr:`~Simulation.reference_model`.
            2. The model must not have the same ID as the
               :attr:`~Simulation.reference_model`
            3. The model must have all numerical values needed for simulation.

        If the criteria is not met, a warning is raised and the reference model
        will not change.

        After changing the reference model, the previous reference model will
        remain included in the :class:`Simulation`.

        Parameters
        ----------
        model : MassModel or str
            Either a new or existing :class:`~.MassModel`, or the string
            identifer of an existing model in the :class:`Simulation` to be set
            as the new reference model.
        verbose : bool
            Whether to output additional and more verbose messages.
            Default is ``False``.

        """
        # If a MassModel is provided, add the model to the Simulation,
        # replacing the model values if model already exists.
        if isinstance(model, MassModel):
            self.add_models(model, verbose=verbose)
            new_model = model
            if str(model) not in self.models:
                new_model = None
        else:
            # Use string identifier to get the model as an object.
            new_model = self.get_model_objects(model)
            try:
                new_model = new_model.pop()
            except IndexError:
                new_model = None

        # Only set the reference model if successful.
        if isinstance(new_model, MassModel):
            setattr(self, "_reference_model", new_model)
            self._roadrunner = _load_model_into_roadrunner(
                new_model, rr=self.roadrunner, verbose=verbose,
                **_SBML_KWARGS)

    def get_model_simulation_values(self, model):
        """Return two dictionaries containing initial and parameter values.

        Parameters
        ----------
        model : MassModel or str
            The model or its identifier whose values are to be returned.

        Returns
        -------
        tuple (init_conds, parameters)
        init_conds : :class:`~mass.util.dict_with_id.DictWithID`
            A :class:`~mass.util.dict_with_id.DictWithID` containing
            initial conditions.
        parameters : :class:`~mass.util.dict_with_id.DictWithID`
            A :class:`~mass.util.dict_with_id.DictWithID` containing
            model parameters.

        """
        try:
            values_dict = self._simulation_values.get_by_id(
                str(model) + "_values")
        except KeyError:
            raise ValueError(
                "MassModel '{0}' does not exist in the Simulation object."
                .format(str(model)))

        return (values_dict["init_conds"], values_dict["parameters"])

    def add_models(self, models, verbose=False):
        r"""Add the model values to the :class:`Simulation`.

        To add a model to the :class:`Simulation`, three criteria must be met:

            1. The model must have equivalent ODEs to the
               :attr:`~Simulation.reference_model`
            2. The model must not have the same ID as the
               :attr:`~Simulation.reference_model`.
            3. The model must have all numerical values needed for simulation.

        Notes
        -----
        * Only the model values are added to the :class:`Simulation`.
        * If a model already exists in the Simulation, it will be replaced.
        * To verify that the model has equivalent ODEs to the reference model,
          use :meth:`.MassModel.has_equivalent_odes`.

        Parameters
        ----------
        models : iterable of models
            An iterable containing the :class:`~.MassModel`\ s to add.
        verbose : bool
            Whether to print if loading of models succeeds or fails.
            Default is ``False``.

        """
        models = ensure_iterable(models)
        invalid_values = [
            model for model in models if not isinstance(model, MassModel)
            or str(model) == str(self.reference_model)]

        if invalid_values:
            raise ValueError(
                'Invalid models found: {0}. Ensure all models are MassModel '
                'objects with IDs different from the reference model before '
                'adding'.format(repr(invalid_values)))

        for model in models:
            if not self.reference_model.has_equivalent_odes(model, verbose):
                _log_msg(LOGGER, logging.WARN, verbose,
                         "Could not load MassModel '%s', ODEs are not "
                         "equivalent to the reference model.", str(model))
                continue
            try:
                _assess_model_quality_for_simulation(model, verbose, False)
            except MassSimulationError as e:
                _log_msg(LOGGER, logging.WARN, verbose,
                         "Could not load MassModel '%s': %s",
                         str(model), str(e))
                continue

            values = _get_sim_values_from_model(model)
            if values in self._simulation_values:
                _log_msg(LOGGER, logging.WARN, verbose,
                         "MassModel '%s' already exists, existing values will "
                         "be replaced", str(model))
                self._simulation_values._replace_on_id(values)
            else:
                self._simulation_values.add(values)

            _log_msg(LOGGER, logging.INFO, verbose,
                     "Successfully loaded MassModel '%s'.", str(model))

    def remove_models(self, models, verbose=False):
        r"""Remove the model values from the :class:`Simulation`.

        Notes
        -----
        The :attr:`~Simulation.reference_model` cannot be removed from the
        :class:`Simulation`. In order to remove the current reference model,
        the reference model must first be changed to a different model using
        the :meth:`set_new_reference_model` method.

        Parameters
        ----------
        models : iterable of models or their identifiers
            An iterable of :class:`.MassModel`\ s or their string identifiers
            to be removed.
        verbose : bool
            Whether to print if removal of models succeeds.
            Default is ``False``.

        """
        models = [str(m) for m in ensure_iterable(models)]
        for model in models:
            try:
                values_dict = self._simulation_values.get_by_id(
                    str(model) + "_values")
            except KeyError:
                warn("MassModel '{0}' does not exist.".format(model))
            else:
                self._simulation_values -= [values_dict]
                _log_msg(LOGGER, logging.INFO, verbose,
                         "Successfully removed MassModel '%s.", str(model))

    def get_model_objects(self, models=None):
        r"""Return the loaded models as :class:`~.MassModel`\ s.

        Notes
        -----
        With the exception of the :attr:`~Simulation.reference_model`, only
        the numerical values of a mdoel are stored in order to improve
        performance. Therefore, when using this method to retrieve the
        :class:`~.MassModel`\ s, all models are created anew, meaning that
        they will NOT be the same :class:`~.MassModel`\ s that were loaded
        into the Simulation.

        Parameters
        ----------
        models : iterable of model identifiers
            An iterable of strings containing the model identifiers of the
            desired :class:`~.MassModel`\ s to return. If ``None`` then all
            models in the :class:`Simulation` will be returned.

        Returns
        -------
        mass_models : ~cobra.core.dictlist.DictList
            A :class:`~cobra.core.dictlist.DictList` containing all of the
            :class:`~.MassModel`\ s.

        """
        if models is None:
            models = self.models
        models = [str(m) for m in ensure_iterable(models)]

        mass_models = DictList()
        for model in models:
            if model == self.reference_model:
                mass_models += [self.reference_model]
                continue

            # Copy reference model
            new_model = self.reference_model.copy()
            new_model.id = model

            try:
                # Try to update model values
                new_model = self._update_mass_model_with_values(new_model)
            except ValueError as e:
                warn(str(e))
            else:
                mass_models += [new_model]

        return mass_models

    def simulate(self, models=None, time=None, perturbations=None, **kwargs):
        r"""Simulate models and return results as :class:`~.MassSolution`\ s.

        A simulation is carried out by simultaneously integrating the ODEs
        of the ``models`` to compute their solutions over the time interval
        specified by ``time``, while temporarily incorporating events and
        changes specified in ``perturbations``.

        Parameters
        ----------
        models : iterable of models or their string identifiers, None
            The models to simulate. If ``None`` then all models loaded into
            the simulation object will be used. All models must already
            exist in the :class:`Simulation`.
        time : tuple
            Either a ``tuple`` containing the initial and final time points, or
            a ``tuple`` containing the initial time point, final time point,
            and the number of time points to use.
        perturbations : dict
            A ``dict`` of perturbations to incorporate into the simulation.
            See :mod:`~.simulation` documentation for more information on
            valid perturbations.
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

                Default is ``None``.
            interpolate :
                ``bool`` indicating whether simulation results should be
                returned to as interpolating functions.

                Default is ``False``.
            update_solutions :
                ``bool`` indicating whether to replace the stored solutions in
                the simulation with the new simulation results.

                Default is ``True``.
            decimal_precision :
                ``bool`` indicating whether to apply the
                :attr:`~.MassBaseConfiguration.decimal_precision` attribute of
                the :class:`.MassConfiguration` to the solution values.

                Default is ``False``.

        Returns
        -------
        tuple (conc_solutions, flux_solutions)
        conc_solutions : :class:`~cobra.core.dictlist.DictList`
            A :class:`~cobra.core.dictlist.DictList` of
            :class:`~.MassSolution`\ s containing the concentration solutions
            for successful simulations. If the simulation failed, the
            :class:`~.MassSolution` will be returned as empty.
        flux_solutions : :class:`~cobra.core.dictlist.DictList`
            A :class:`~cobra.core.dictlist.DictList` of
            :class:`~.MassSolution`\ s containing the flux solutions for
            successful simulations. If the simulation failed, the
            :class:`~.MassSolution` will be returned as empty.

        """
        # Check kwargs
        kwargs = _check_kwargs({
            "verbose": False,
            "selections": None,
            "boundary_metabolites": False,
            "steps": None,
            "interpolate": False,
            "update_solutions": True,
            "decimal_precision": False}, kwargs)
        # Set all models for simulation if None provided.
        if models is None:
            models = self.models
        models = [str(m) for m in ensure_iterable(models)]

        # Get roadrunner instance
        rr = self.roadrunner

        # Parse the time and format the input for the roadrunner
        time = _format_time_input(time, steps=kwargs.get("steps"),
                                  verbose=kwargs.get("verbose"))

        # Parse the perturbations and format the input to be used
        perturbations = self._format_perturbations_input(perturbations,
                                                         kwargs.get("verbose"))
        # Make the time course selection input and set the selections
        selections = self._make_rr_selections(kwargs.get("selections"),
                                              kwargs.get("verbose"))

        # Make DictLists for solution objects
        conc_sol_list = DictList()
        flux_sol_list = DictList()
        try:
            for model in models:
                try:
                    # Apply perturbations and set values in roadrunner
                    rr, reset = self._set_simulation_values(
                        model, perturbations, kwargs.get("verbose"))
                    # Simulate
                    rr.timeCourseSelections = selections
                    _log_msg(LOGGER, logging.INFO, kwargs.get("verbose"),
                             "Simulating '%s'", str(model))
                    results = rr.simulate(*time)

                    # Map results to their identifiers and return MassSolutions
                    _log_msg(LOGGER, logging.INFO, kwargs.get("verbose"),
                             "Simulation for '%s' successful", str(model))
                    solutions = self._make_mass_solutions(
                        model, default_selections=selections, results=results,
                        **kwargs)

                # Handle MassSimulationErrors
                except (RuntimeError, MassSimulationError) as e:
                    warn("One or more simulations failed. Check the log "
                         "for more details.")
                    if kwargs.get("verbose"):
                        print("Failed simulation for '{0}'".format(model))
                    LOGGER.error(
                        "Failed simulation for '%s' due the following error: "
                        "%s", model, str(e))
                    # Make empty MassSolutions
                    solutions = self._make_mass_solutions(
                        model, default_selections=selections, results=None,
                        **kwargs)
                    reset = False

                finally:
                    # Add solutions to overall simulation output
                    _log_msg(LOGGER, logging.INFO, kwargs.get("verbose"),
                             "Adding '%s' simulation solutions to output",
                             str(model))
                    conc_sol_list += [solutions[0]]
                    flux_sol_list += [solutions[1]]
                    # Reset the roadrunner state
                    self._reset_roadrunner(reset)

        # Handle unforseen errors as critical errors
        except Exception as e:
            # If an error occurs, raise it as a MassSimulationError
            raise MassSimulationError(
                "Critical simulation fail due to the following:\n"
                + str(e.__class__.__name__) + ": " + str(e))

        # Update the solutions stored in the Simulation
        if kwargs.get("update_solutions"):
            _log_msg(LOGGER, logging.INFO, kwargs.get("verbose"),
                     "Updating stored solutions")
            self._update_stored_solutions(_CONC_STR, conc_sol_list)
            self._update_stored_solutions(_FLUX_STR, flux_sol_list)

        return conc_sol_list, flux_sol_list

    def find_steady_state(self, models=None, strategy="nleq2",
                          perturbations=None, update_values=False, **kwargs):
        r"""Find steady states for models.

        The steady state is found by carrying out the provided strategy.

        * The ``'simulate'`` strategy will simulate the model for a long time
          (default ``1e8``), and ensure the absolute difference between
          solutions at the final two time points is less than the
          :attr:`~.MassBaseConfiguration.steady_state_threshold` in the
          :class:`.MassConfiguration`.
        * Other strategies involve using the
          :class:`roadrunner.roadrunner.SteadyStateSolver` class to determine
          the steady state through global Newtonian methods. The steady state
          is found when the sum of squares of the rates of change is less than
          the :attr:`~.MassBaseConfiguration.steady_state_threshold` in the
          :class:`.MassConfiguration`.

        Parameters
        ----------
        models : iterable of models or their string identifiers, None
            The models to simulate. If ``None`` then all models loaded into
            the simulation object will be used. All models must already
            exist in the :class:`Simulation`.
        strategy : str
            The strategy for finding the steady state. Must be one of the
            following:

                * ``'simulate'``
                * ``'nleq1'``
                * ``'nleq2'``

        perturbations : dict
            A ``dict`` of perturbations to incorporate into the simulation.
            See :mod:`~.simulation` documentation for more information on
            valid perturbations.
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
        tuple (conc_solutions, flux_solutions)
        conc_solutions : :class:`~cobra.core.dictlist.DictList`
            A :class:`~cobra.core.dictlist.DictList` of
            :class:`~.MassSolution`\ s containing the concentration solutions
            for successful simulations. If the simulation failed, the
            :class:`~.MassSolution` will be returned as empty.
        flux_solutions : :class:`~cobra.core.dictlist.DictList`
            A :class:`~cobra.core.dictlist.DictList` of
            :class:`~.MassSolution`\ s containing the flux solutions for
            successful simulations. If the simulation failed, the
            :class:`~.MassSolution` will be returned as empty.

        """
        # Check kwargs
        kwargs = _check_kwargs({
            "verbose": False,
            "selections": None,
            "boundary_metabolites": False,
            "steps": None,
            "tfinal": 1e8,
            "num_attempts": 2,
            "decimal_precision": False}, kwargs)
        # Set all models for simulation if None provided.
        if models is None:
            models = self.models
        models = [str(model) for model in ensure_iterable(models)]

        # Get roadrunner instance
        rr = self.roadrunner

        # Ensure strategy input is valid
        if strategy not in rr.getRegisteredSteadyStateSolverNames()\
           and strategy != "simulate":
            raise ValueError(
                "Invalid steady state strategy: '{0}'".format(strategy))
        if strategy == "simulate":
            steady_state_function = self._find_steady_state_simulate
        else:
            steady_state_function = self._find_steady_state_solver
            rr.setSteadyStateSolver(strategy)

        # Parse the perturbations and format the input to be used
        perturbations = self._format_perturbations_input(perturbations,
                                                         kwargs.get("verbose"))

        # Set species to use for steady state calculations
        selections = self._make_rr_selections(kwargs.get("selections"),
                                              kwargs.get("verbose"))
        # Remove time from selections
        selections.remove("time")

        # Make DictLists for solution objects
        conc_sol_list = DictList()
        flux_sol_list = DictList()
        try:
            for model in models:
                try:
                    # Apply perturbations and set values in roadrunner
                    rr, reset = self._set_simulation_values(
                        model, perturbations, kwargs.get("verbose"))
                    # Use simulate strategy
                    rr.steadyStateSelections = selections
                    results = steady_state_function(model, **kwargs)
                    # Map results to their identifiers and return MassSolutions
                    solutions = self._make_mass_solutions(
                        model, default_selections=selections, results=results,
                        update_values=update_values, **kwargs)

                # Handle MassSimulationErrors
                except (RuntimeError, MassSimulationError) as e:
                    warn("Unable to find a steady state for one or more "
                         "models. Check the log for more details.")
                    if kwargs.get("verbose"):
                        print("Failed to find steady state for {0}'".format(
                            model))
                    LOGGER.error(
                        "Unable to find a steady state for '%s' "
                        "using strategy '%s' due to the following: %s",
                        model, strategy, str(e))
                    # Make empty MassSolutions
                    solutions = self._make_mass_solutions(
                        model, default_selections=selections, results=None,
                        update_values=False, **kwargs)
                    reset = False

                finally:
                    # Add solutions to overall simulation output
                    _log_msg(LOGGER, logging.INFO, kwargs.get("verbose"),
                             "Adding '%s' simulation solutions to output",
                             str(model))
                    # Add solutions to output lists
                    conc_sol_list += [solutions[0]]
                    flux_sol_list += [solutions[1]]
                    # Reset the roadrunner state
                    self._reset_roadrunner(reset)

        # Handle unforseen errors as critical errors
        except Exception as e:
            # If an error occurs, raise it as a MassSimulationError
            raise MassSimulationError(
                "Critical simulation fail due to the following:\n"
                + str(e.__class__.__name__) + ": " + str(e))

        # Update reference model to have the new values
        if update_values:
            _log_msg(LOGGER, logging.INFO, kwargs.get("verbose"),
                     "Updating reference model values")
            model = self._update_mass_model_with_values(self.reference_model)
            setattr(self, "_reference_model", model)
            self._reset_roadrunner(True)

        return conc_sol_list, flux_sol_list

    def _make_rr_selections(self, selections=None, verbose=False):
        """Set the observable output of the simulation.

        Warnings
        --------
        This method is intended for internal use only.

        """
        # Get roadrunner executable model instance and reference model
        rr_model = self.roadrunner.model

        _log_msg(LOGGER, logging.INFO, verbose, "Setting output selections")
        rr_selections = ["time"]
        if selections is None:
            rr_selections += rr_model.getFloatingSpeciesConcentrationIds()
            rr_selections += rr_model.getReactionIds()
        else:
            # TODO Catch selection errors here
            rr_selections += selections

        return rr_selections

    def _format_perturbations_input(self, perturbations, verbose=False):
        """Check and format the perturbation input.

        Perturbations are checked before simulations are carried out to limit
        fails during the simulation due to bad syntax or values.

        Warnings
        --------
        This method is intended for internal use only.

        """
        _log_msg(LOGGER, logging.INFO, verbose, "Parsing perturbations")
        formatted_perturbations = {}
        # Check input type.
        if not isinstance(perturbations, dict) and perturbations is not None:
            raise TypeError("Perturbations must be a dict")

        # Get the reference model initial conditions and parameters
        if perturbations:
            sim_values = self._get_all_values_for_sim(self.reference_model)

            for key, value in iteritems(perturbations):
                # Ensure key exists in model values
                if key not in sim_values\
                   and _make_init_cond(key) not in sim_values:
                    raise ValueError(
                        "Invalid Perturbation: '{0}' not found in model "
                        "simulation values".format(key))
                try:
                    value = float(value)
                except (ValueError, TypeError):
                    if key in self.reference_model.boundary_metabolites:
                        value = sympify(value, locals={key: Symbol(key)})
                        if any([arg.atoms(Symbol) != {Symbol("t")}
                                for arg in list(value.atoms(Function))]):
                            raise ValueError(
                                "Function '{0}' for '{1}' has independent "
                                "variables other than 't' for time".format(
                                    value, key))
                    elif key in str(value):
                        value = sympify(value, locals={key: Symbol(key)})
                    else:
                        raise ValueError(
                            "Invalid Perturbation {{'{0}': '{1}'}}".format(
                                key, value))

                formatted_perturbations[key] = value
        return formatted_perturbations

    def _set_simulation_values(self, model, perturbations, verbose=False):
        """Set the simulation numerical values in the roadrunner instance.

        Warnings
        --------
        This method is intended for internal use only.

        """
        # Get roadrunner
        rr = self.roadrunner
        _log_msg(LOGGER, logging.INFO, verbose,
                 "Setting simulation values for '%s'", str(model))
        try:
            # Ensure values exist for the model
            sim_values_to_set = self._get_all_values_for_sim(model)
        except ValueError as e:
            raise MassSimulationError(e)

        # Apply perturbations to model values. Set flag for
        # reloading model into the roadrunner instance.
        reset = False
        try:
            for key, value in iteritems(perturbations):
                # Perturb value to a number if value is float
                if isinstance(value, float):
                    sim_values_to_set[key] = value
                    continue

                sim_values = self._get_all_values_for_sim(self.reference_model)
                # Determine type of perturbation
                if _make_init_cond(key) in sim_values:
                    accessor_key = _make_init_cond(key)
                else:
                    accessor_key = key

                value = value.subs({key: sim_values_to_set[accessor_key]})
                try:
                    value = float(value)
                except (ValueError, TypeError):
                    reset = True

                sim_values_to_set[accessor_key] = value
                continue

        except (ValueError, MassSimulationError) as e:
            raise MassSimulationError(e)

        # Set roadrunner to reflect given model variant
        rr = self._set_values_in_roadrunner(model, reset, sim_values_to_set)

        # Return the roadrunner instance and whether it will need a reload for
        return rr, reset

    def _set_values_in_roadrunner(self, model, reset, sim_values_to_set):
        """Set the roadrunner values to reflect the given model.

        Warnings
        --------
        This method is intended for internal use only.

        """
        rr = self.roadrunner
        if reset:
            new_model = self.reference_model.copy()
            new_model.id = model
            new_model = self._update_mass_model_with_values(
                new_model, sim_values_to_set)

            rr = _load_model_into_roadrunner(
                new_model, rr=self.roadrunner, verbose=False,
                **_SBML_KWARGS)
        else:
            for key, value in iteritems(sim_values_to_set):
                if isinstance(value, Basic):
                    continue
                rr.setValue(key, value)

        return rr

    def _make_mass_solutions(self, model, default_selections, results,
                             update_values=False, **kwargs):
        """Make the MassSolutions using the results of the Simulation.

        Warnings
        --------
        This method is intended for internal use only.

        """
        # Get roadrunner instance
        rr = self.roadrunner

        # Create dicts for containing concentrations and flux solutions.
        time = None
        conc_sol = {}
        flux_sol = {}
        if results is not None:
            # Get species and reaction lists
            species_list = rr.model.getFloatingSpeciesConcentrationIds()
            if kwargs.get("boundary_metabolites"):
                species_list += rr.model.getBoundarySpeciesConcentrationIds()
            reaction_list = rr.model.getReactionIds()

            selections = kwargs.get("selections")
            if not selections:
                selections = default_selections

            for sol_key in selections:
                # Get the time vector
                if sol_key == "time":
                    time = results[sol_key]
                elif sol_key in species_list:
                    # Add solution to concentration solutions if specie
                    result_array = _round_values(
                        results[sol_key], kwargs.get("decimal_precision"))
                    conc_sol[_strip_conc(sol_key)] = result_array
                elif sol_key in reaction_list:
                    # Add solution to flux solutions if reaction
                    result_array = _round_values(
                        results[sol_key], kwargs.get("decimal_precision"))
                    flux_sol[sol_key] = result_array
                else:
                    pass

        # Make a MassSolution object of the concentration solution dict.
        conc_sol = MassSolution(id_or_model=model, data_dict=conc_sol,
                                solution_type=_CONC_STR, time=time,
                                interpolate=kwargs.get("interpolate", False))

        # Make a MassSolution object of the concentration solution dict.
        flux_sol = MassSolution(id_or_model=model, data_dict=flux_sol,
                                solution_type=_FLUX_STR, time=time,
                                interpolate=kwargs.get("interpolate", False))

        if update_values:
            _log_msg(LOGGER, logging.INFO, kwargs.get("verbose"),
                     "Updating '%s' values", model)
            self._update_simulation_values(model, (conc_sol, flux_sol))

        return conc_sol, flux_sol

    def _update_stored_solutions(self, solution_type, solutions):
        """Update stored MassSolutions with new MassSolution objects.

        Warnings
        --------
        This method is intended for internal use only.

        """
        # Determine where solutions are stored
        stored_solution_dictlist = {
            "Conc": self._concentration_solutions,
            "Flux": self._flux_solutions}[solution_type]

        # Add solutions to their corresponding DictList
        for solution in solutions:
            if solution.id in stored_solution_dictlist:
                # Remove old solutions
                old_solution = stored_solution_dictlist.get_by_id(solution.id)
                old_solution._simulation = None
                # Replace with new solution
                stored_solution_dictlist._replace_on_id(solution)
            else:
                # Add the new solution
                stored_solution_dictlist.add(solution)
            # Set simulation reference in the MassSolution
            solution._simulation = self

    def _reset_roadrunner(self, reset):
        """Reset the RoadRunner to its the original state.

        Warnings
        --------
        This method is intended for internal use only.

        """
        if reset:
            _load_model_into_roadrunner(
                self.reference_model, rr=self.roadrunner, verbose=False,
                **_SBML_KWARGS)
        else:
            self.roadrunner.resetToOrigin()
            LOGGER.info("Reset roadrunner to origin.")

    def _find_steady_state_simulate(self, model, **kwargs):
        """Find the steady state of a model through simulation of the model.

        Warnings
        --------
        This method is intended for internal use only.

        """
        rr = self.roadrunner
        # Set simulation output to be the steady state selections
        rr.timeCourseSelections = rr.steadyStateSelections

        # Parse the time and format the input for the roadrunner
        steps = kwargs.get("steps")
        if steps is None:
            steps = int(max(kwargs.get("tfinal") / 1000, 10000))
        time = _format_time_input((0, kwargs.get("tfinal")), steps=steps,
                                  verbose=kwargs.get("verbose"))
        try:
            # Simulate for a long time
            _log_msg(LOGGER, logging.INFO, kwargs.get("verbose"),
                     "Simulating '%s'", str(model))
            simulation_results = rr.simulate(*time)
        except (RuntimeError, TypeError) as e:
            raise MassSimulationError(
                str(e.__class__.__name__) + ": " + str(e))

        final_points = np.array([simulation_results[key][-2:]
                                 for key in rr.timeCourseSelections])
        abs_diff = np.abs(final_points[:, 0] - final_points[:, -1])
        if kwargs.get("decimal_precision"):
            abs_diff = apply_decimal_precision(
                abs_diff, MASSCONFIGURATION.decimal_precision)

        abs_diff = abs_diff < MASSCONFIGURATION.steady_state_threshold
        if not all(abs_diff):
            raise MassSimulationError(
                'For MassModel "{0}", absolute difference for "{1!r}" is '
                'greater than the steady state threshold.'.format(model, [
                    key for i, key in enumerate(rr.timeCourseSelections)
                    if not abs_diff[i]]))

        _log_msg(LOGGER, logging.INFO, kwargs.get("verbose"),
                 "Found steady state for '%s'.", str(model))
        ss_results = {key: final_points[i, -1]
                      for i, key in enumerate(rr.timeCourseSelections)}

        return ss_results

    def _find_steady_state_solver(self, model, **kwargs):
        """Find the steady state of the model using a RoadRunner solver method.

        Warnings
        --------
        This method is intended for internal use only.

        """
        rr = self.roadrunner

        def is_steady_state(rr):
            """Solve for steady state via RoadRunner."""
            # Simulate for a long time
            is_ss_value = rr.steadyState()
            # Round value before comparision
            if kwargs.get("decimal_precision"):
                is_ss_value = apply_decimal_precision(
                    is_ss_value, MASSCONFIGURATION.decimal_precision)
            if is_ss_value <= MASSCONFIGURATION.steady_state_threshold:
                return True
            return False

        # See if a steady state can be reached
        i = 0
        success = False
        try:
            while not success and i < kwargs.get("num_attempts"):
                _log_msg(LOGGER, logging.INFO, kwargs.get("verbose"),
                         "Attempt %s to calculate steady state for '%s'",
                         str(i + 1), str(model))
                success = is_steady_state(rr)
                i += 1
        except RuntimeError as e:
            raise MassSimulationError(
                str(e.__class__.__name__) + ": " + str(e))

        if not success:
            raise MassSimulationError(
                "Could not find steady state for MassModel '{0}' after {1:d} "
                "attempts, steady state threshold always exceeded.".format(
                    model, i))

        _log_msg(LOGGER, logging.INFO, kwargs.get("verbose"),
                 "Found steady state for '%s'.", str(model))
        # Zip the solutions with their IDs into a dict and return
        results = zip(rr.steadyStateSelections, rr.getSteadyStateValues())
        return dict(results)

    def _update_simulation_values(self, model, values):
        """Update the stored values for the model with the new ones.

        Warnings
        --------
        This method is intended for internal use only.

        """
        sim_values = self.get_model_simulation_values(model)
        id_fix_dict = {model + "_init_conds": _make_init_cond,
                       model + "_parameters": _make_ss_flux}
        for current_value_dict, new_value_dict in zip(sim_values, values):
            id_fix_func = id_fix_dict[current_value_dict.id]
            for key, value in iteritems(new_value_dict):
                current_value_dict[id_fix_func(key)] = value

    def _update_mass_model_with_values(self, mass_model, value_dict=None):
        """Update the MassModel object with the stored model values.

        Warnings
        --------
        This method is intended for internal use only.

        """
        # Get the model values
        init_conds, parameters = self.get_model_simulation_values(mass_model)

        # Update with provided model values
        if value_dict is not None:
            init_conds = dict(
                (_strip_init_cond(key), value_dict[key]) if "init" in key
                else (key, value_dict[key]) for key in init_conds)
            parameters = {key: value_dict[key] for key in parameters}
        else:
            init_conds = dict(
                (_strip_init_cond(key), ic) if "init" in key
                else (key, ic) for key, ic in iteritems(init_conds))
        # Update the model initial conditions
        mass_model.update_initial_conditions(init_conds, verbose=False)
        # Update the model parameter values
        mass_model.update_parameters(parameters, verbose=False)

        return mass_model

    def _get_all_values_for_sim(self, mass_model):
        """Get all model values as a single dict.

        Warnings
        --------
        This method is intended for internal use only.

        """
        init_conds, parameters = self.get_model_simulation_values(mass_model)
        sim_values = {}
        sim_values.update(init_conds)
        sim_values.update(parameters)
        return sim_values


def _assess_model_quality_for_simulation(mass_model, verbose, report):
    """Assess the model quality to ensure it can be simulated.

    Warnings
    --------
    This method is intended for internal use only.

    """
    simulate_check, consistency_check = is_simulatable(mass_model)
    if not simulate_check:
        msg = ""
        if verbose:
            msg += "MassModel is missing numerical values that are " + \
                   "necessary for simulation."
        if report:
            qcqa_model(mass_model, **{
                "parameters": True,
                "concentrations": True,
                "fluxes": True,
                "superfluous": True,
                "elemental": True})

        raise MassSimulationError(msg)

    if not consistency_check:
        msg = "MassModel has numerical consistency issues, use with caution."
        if verbose:
            msg += " To help determine which values are not numerically " + \
                   "consistent, use the `qcqa_model` method in mass.util.qcqa."
        LOGGER.warning(msg)


def _load_model_into_roadrunner(mass_model, rr=None, verbose=False, **kwargs):
    """Create a RoadRunner instance for the given model.

    Warnings
    --------
    This method is intended for internal use only.

    """
    doc = _model_to_sbml(mass_model, f_replace={}, **kwargs)
    sbml_str = writeSBMLToString(doc)
    error = roadrunner.validateSBML(sbml_str)
    if error:
        msg = "Cannot load SBML Model '" + str(mass_model) + "' "
        if verbose:
            msg += str("into RoadRunner due to error:\n" + error + " Make "
                       "sure the model can be exported to a valid SBML model "
                       "using the 'validate_model_to_SBML' method in the SBML "
                       "submodule to ensure there are no SBML compliance "
                       "issues with the model.")
        raise MassSimulationError(msg)

    if rr is None:
        rr = roadrunner.RoadRunner(sbml_str)
    else:
        rr.clearModel()
        rr.load(sbml_str)

    msg = "Successfully loaded MassModel '{0}' into RoadRunner.".format(
        str(mass_model))
    LOGGER.info(msg)
    if verbose:
        print(msg)

    return rr


def _get_sim_values_from_model(mass_model):
    """Return the initial conditions and parameter values in a DictWithID.

    Warnings
    --------
    This method is intended for internal use only.

    """
    init_conds = mass_model.initial_conditions
    parameters = mass_model._get_all_parameters()

    values = {
        "init_conds": DictWithID(
            id=mass_model.id + "_init_conds", data_dict=dict(
                (_make_init_cond(met.id), ic)
                if not met.fixed else (met.id, ic)
                for met, ic in iteritems(init_conds))),
        "parameters": DictWithID(
            id=mass_model.id + "_parameters", data_dict={
                param: value
                for param, value in iteritems(parameters)}),
    }

    return DictWithID(id=mass_model.id + "_values", data_dict=values)


def _format_time_input(time, steps=None, verbose=False):
    """Format the time input and return it as a tuple for the RoadRunner.

    Warnings
    --------
    This method is intended for internal use only.

    """
    _log_msg(LOGGER, logging.INFO, verbose, "Getting time points")
    if len(time) == 2:
        t0, tf = time
        numpoints = None
    elif len(time) == 3:
        t0, tf, numpoints = time
        numpoints = int(numpoints)
    else:
        raise TypeError("The 'time' input mut be one of the following: "
                        "'(t0, tf)' or '(t0, tf, numpoints)'")
    time = (t0, tf, numpoints, None, steps)

    return time


def _round_values(values, decimal_precision=True):
    """Round a value to the decimal precision in the MassConfiguration.

    Warnings
    --------
    This method is intended for internal use only.

    """
    if decimal_precision:
        values = apply_decimal_precision(values,
                                         MASSCONFIGURATION.decimal_precision)

    return values


def _strip_conc(metabolite_str):
    """Strip the concentration brackets from a metabolite identifier.

    Warnings
    --------
    This method is intended for internal use only.

    """
    return metabolite_str.lstrip("[").rstrip("]")


def _make_conc(metabolite_str):
    """Add the concentration brackets from a metabolite identifier.

    Warnings
    --------
    This method is intended for internal use only.

    """
    return "[" + metabolite_str + "]"


def _strip_init_cond(metabolite_str):
    """Strip the initial condition format from a metabolite identifier.

    Warnings
    --------
    This method is intended for internal use only.

    """
    return metabolite_str[5:-1]


def _make_init_cond(metabolite_str):
    """Format metabolite identifier to match initial condition accesors.

    Warnings
    --------
    This method is intended for internal use only.

    """
    return "init(" + metabolite_str + ")"


def _make_ss_flux(reaction_str):
    """Format reaction identifier to match steady state flux parameter.

    Warnings
    --------
    This method is intended for internal use only.

    """
    return "v_" + reaction_str


def _log_msg(logger, level, verbose, msg, *args):
    """TODO DOCSTRING."""
    logger.log(level, msg, *args)
    if verbose:
        print(msg % args)


__all__ = ("Simulation", "LOGGER", "RR_LOGGER",)
