# -*- coding: utf-8 -*-
"""The Simulation module addresses the simulation of mass.MassModels.

The Simulation module is designed to address all aspects related to the
simulation of one or more MassModel objects. These aspects include initializing
and compiling the model into a RoadRunner instance (the ODE integrator),
storing the numerical values necessary for the simulations, and handling the
simulation results by creating and storing MassSolution objects.

A Simulation is initialized by providing a MassModel object

Perturbations can also be implemented for a given simulation as long as they
follow the following guidelines:
    1. Perturbations are dicts where key:value pairs are the variables to be
       perturbed and the corresponding value change.
    2. To scale the current value of a variable, the value should be a string
       representation of a formula altering the variable where the variable is
       identical to the perturbation key.
    3. If providing a formula string as the perturbation value, it must be
       possible to sympify the string via sympy.sympify. It must only have one
       variable that is identical to the perturbation key.
    4. Boundary metabolites can be changed to have functions of time that
       represent the external concentration at that point in time.
       If a perturbation value is to be a function, it must be a function of
       time and it must be possible to sympify the string via sympy.sympify.

Examples of perturbations following the guidelines for a model containing the
specie with ID 'MID_c', boundary metabolite with ID 'MID_b', and reaction
with ID 'RID':
    Altering initial conditions:
        {'init(MID_c)': float}, {'init(MID_c)': 'init(MID_c) * value'}
        {'init([MID_c])': float}, {'init([MID_c])': 'init([MID_c]) * value'}
    Setting a metabolite constant:
        {'fixed(MID_c)': float}, {'fixed(MID_c)': 'fixed(MID_c) * value'}
        {'fixed([MID_c])': float}, {'fixed([MID_c])': 'fixed([MID_c]) * value'}
        {'fixed(MID_b)': float}, {'fixed(MID_b)': 'fixed(MID_b) * value'}
        {'fixed([MID_b])': float}, {'fixed([MID_b])': 'fixed([MID_b]) * value'}
    Altering a reaction parameter:
        {'kf_RID': float}, {'kf_RID': 'kf_RID * value'}
        {'kr_RID': float}, {'kr_RID': 'kr_RID * value'}
        {'Keq_RID': float}, {'Keq_RID': 'Keq_RID * value'}
    Setting a function:
        {'func(MID_b)': 'sin(t)'}, {'func(MID_b)': 'sin(t) + func(MID_b)'}

All simulation results are returned as mass.MassSolution objects. MassSolutions
are specialized dictionaries with additional attributes and properties to help
with accessing, grouping, and plotting solutions. They can also behave as
booleans, with empty MassSolution objects returning as False, and those with
solutions returning as True.
"""
from __future__ import absolute_import

import re

from libsbml import writeSBMLToString

import numpy as np

import roadrunner

from six import iteritems

from sympy import sympify

from cobra.core.dictlist import DictList
from cobra.core.object import Object

from mass.core.massconfiguration import MassConfiguration
from mass.core.massmodel import MassModel
from mass.core.masssolution import MassSolution, _CONC_STR, _FLUX_STR
from mass.exceptions import MassSimulationError
from mass.io.sbml import _model_to_sbml
from mass.util.DictWithID import DictWithID
from mass.util.util import ensure_iterable
# Set the logger
MASSCONFIGURATION = MassConfiguration()
LOGGER = roadrunner.Logger
# Pre-compiled regex for perturbations
# re.compile(r"param\(([\s|\w]*)\)"),
PERTURBATIONS_RE_DICT = {
    "init": re.compile(r"init\([\s|\w]*\)"),
    "fixed": re.compile(r"fixed\((?P<metabolite>[\s|\w]*)\)"),
    "kf": re.compile(r"kf_[\s|\w]*"),
    "Keq": re.compile(r"Keq_[\s|\w]*"),
    "kr": re.compile(r"kr_[\s|\w]*"),
    "func": re.compile(r"func\((?P<metabolite>[\s|\w]*)\)"),
}


class Simulation(Object):
    """Class for managing setup and result handling of simulations.

    The mass.Simulation class is designed to address all aspects related to the
    simulation of MassModel objects, some of which includes setting solver
    and solver options, perturbation of concentrations and parameters,
    simulation of a single MassModel or an ensemble of models, and handling
    of simulation results.

    Parameters
    ----------
    mass_model: mass.MassModel
        The MassModel to load into the Simulation. The model will be set as the
        Simulation's reference model.
    simulation_id: str, optional
        An identifier to associate with the Simulation given as a string.
        If None provided, will default to "(MassModel.id)_Simulation."
    simulation_name: str, optional
        A human readable name for the Simulation.
    verbose: bool, optional
        Whether to provide a QCQA report and more verbose messages when trying
        to load the model.

    """

    def __init__(self, mass_model=None, simulation_id=None,
                 simulation_name=None, verbose=False):
        """Initialize the Simulation Object."""
        if not isinstance(mass_model, MassModel):
            raise TypeError(
                "'{0}' is not a valid MassModel instance".format(
                    str(mass_model)))

        if simulation_id is None:
            simulation_id = "{0}_Simulation".format(str(mass_model))
        try:
            # QCQA model
            _assess_model_quality_for_simulation(mass_model, verbose)
            # Load model into RoadRunner
            rr = _create_roadrunner_instance(
                mass_model, f_replace={}, verbose=verbose,
                **{"units": False, "local_parameters": False})
        except MassSimulationError as e:
            msg = "Could not load MassModel in Simulation object."
            if verbose:
                msg += ": " + str(e)
            raise MassSimulationError(msg)

        # Initialize Simulation
        super(Simulation, self).__init__(simulation_id, simulation_name)
        # Store the original model used to create the Simulation
        self._reference_massmodel = mass_model
        # Set roadrunner
        self._roadrunner = rr

        # Storing model values for simulations
        self._model_values = DictList()
        self._model_values.add(_get_sim_values_from_model(mass_model))

        # Storing concentration and flux solutions in simulations
        self._conc_solutions = DictList()
        self._flux_solutions = DictList()

    @property
    def reference_model(self):
        """Return the reference MassModel of the Simulation."""
        return getattr(self, "_reference_massmodel")

    @property
    def models(self):
        """Return the identifiers of models that exist in the Simulation."""
        return [d.id.replace("_values", "") for d in self._model_values]

    @property
    def roadrunner(self):
        """Return the roadrunner instance of the Simulation."""
        return self._roadrunner

    def get_model_values(self, model):
        """Return two DictWithIDs containing initial and parameter values.

        The first DictWithID contains the model initial condition values.
        The second DictWithID contains the model parameter values.

        """
        try:
            values_dict = self._model_values.get_by_id(str(model) + "_values")
        except KeyError:
            raise ValueError(
                "Model '{0}' does not exist in the Simulation object."
                .format(str(model)))

        return (values_dict["init_conds"], values_dict["parameters"])

    def simulate(self, models=None, time=None, perturbations=None, **kwargs):
        """Simulate models and return solutions as MassSolutions.

        A simulation is carried out by simultaneously integrating the ODEs
        of the 'models' to compute their solutions over the time interval
        specified by 'time', while temporarily incorporating events and changes
        specified in 'perturbations'. See the Simulation Module documentation
        for more information.

        Parameters
        ----------
        models: iterable of MassModels or iterable of MassModel IDs, None
            The models to simulate. If None provided, all models loaded into
            the simulation object will be used. All models must already
            exist in the Simulation.
        time: tuple of len 2 or len 3
            Either tuple of len 2 containing the start and end time points, or
            the a tuple of len 3 containing the starting time, ending time, and
            the number of time points to use.
        perturbations: dict, optional
            A dict of perturbations to incorporate into the simulation.
            See Simulation Module Documentation for more information on valid
            perturbations.
        **kwargs
            steps
            selections
            boundary_metabolites
            update_solutions

        Returns
        -------
        conc_solutions: DictList of mass.MassSolution, MassSolution
            A DictList of MassSolutions containing the concentration solutions
            for successful simulations. If the simulation failed, the
            MassSolution will be returned as empty. If there is only one
            model, a single MassSolution is returned instead of the DictList.
        flux_solutions: DictList of mass.MassSolution, MassSolution
            A DictList of MassSolutions containing the flux solutions for
            successful simulations. If the simulation failed, the MassSolution
            will be returned as empty. If there is only one model, a single
            MassSolution is returned instead of the DictList.

        """
        if kwargs:
            # TODO parse kwarg arguments
            pass

        # Set all models for simulation if None provided.
        if models is None:
            models = self.models
        models = [str(m) for m in ensure_iterable(models)]

        # Parse the time and format the input for the roadrunner
        time = _format_time_input(time, steps=kwargs.get("steps", None))
        # Parse the perturbations and format the input to be used
        perturbations = self._format_perturbations_input(perturbations)

        # Get roadrunner instance
        rr = self.roadrunner

        # Make the time course selection input and set the selections
        selections = self._make_rr_selections(
            selections=kwargs.get("selections", None),
            boundary_metabolites=kwargs.get("boundary_metabolites", False))

        # Make DictLists for solution objects
        conc_sol_list = DictList()
        flux_sol_list = DictList()
        try:
            for model in models:
                try:
                    # Apply perturbations and set model values in roadrunner
                    rr = self._set_simulation_values(model, perturbations)
                    # Simulate
                    rr.timeCourseSelections = selections
                    results = rr.simulate(*time)
                    # Map results to their identifiers and return MassSolutions
                    solutions = self._make_mass_solutions(
                        model, default_selections=selections, results=results,
                        **kwargs)

                # Handle MassSimulationErrors
                except MassSimulationError as e:
                    LOGGER.log(
                        LOGGER.LOG_ERROR, "Failed simulation for '{0}' due "
                        "the following MassSimulationError: {1}".format(
                            model, str(e)))
                    # Make empty MassSolutions
                    solutions = self._make_mass_solutions(
                        model, default_selections=selections, results=None,
                        **kwargs)
                finally:
                    # Add solutions to overall simulation output
                    conc_sol_list += [solutions[0]]
                    flux_sol_list += [solutions[1]]
                    # Reset the roadrunner state
                    rr = _reset_roadrunner_instance(rr)

        # Handle unforseen errors as critical errors
        except Exception as e:
            # If an error occurs, raise it as a MassSimulationError
            raise MassSimulationError(
                "Critical simulation fail due to the following:\n"
                + str(e.__class__.__name__) + ": " + str(e))

        # Update the solutions stored in the Simulation
        if kwargs.get("update_solutions", True):
            self._update_stored_solutions(_CONC_STR, conc_sol_list)
            self._update_stored_solutions(_FLUX_STR, flux_sol_list)

        # Return just a tuple of two MassSolutions if only one model simulated
        if len(models) == 1:
            return conc_sol_list[0], flux_sol_list[0]

        return conc_sol_list, flux_sol_list

    def find_steady_state(self, models=None, strategy="nleq2",
                          perturbations=None, update_values=False, **kwargs):
        """Find steady states for models and return solutions as MassSolutions.

        The steady state is found by carrying out the provided strategy.

        The 'simulate' strategy will simulate the model for a long time (1e8),
        and ensure the absolute difference between solutions at the final two
        time points is smaller than the steady state threshold specified in the
        MassConfiguration.

        All other strategies involve using the Roadrunner steady state solver
        to determine the steady state through root finding methods. The steady
        state is found when the sum of squares of the rates of change is less
        than the steady state threshold specified in the MassConfiguration.

        Parameters
        ----------
        models: iterable of MassModels or iterable of MassModel IDs, None
            The models to simulate. If None provided, all models loaded into
            the simulation object will be used. All models must already
            exist in the Simulation.
        strategy: {'simulate', 'nleq1', nleq2'}
            The strategy for finding the steady state.
        perturbations: dict, optional
            A dict of perturbations to incorporate into the simulation.
            See Simulation Module Documentation for more information on valid
            perturbations.
        update_values: bool, optional
            Whether to update the stored model values with the steady state
            results.
            Default is False.
        **kwargs
            selections
            boundary_metabolites
            tfinal
            steps
            num_attempts

        Returns
        -------
        conc_solutions: DictList of mass.MassSolution, MassSolution
            A DictList of MassSolutions containing the concentration solutions
            for successfully finding steady state. If no steady state is found,
            the MassSolution will be returned as empty. If there is only one
            model, a single MassSolution is returned instead of the
            DictList.
        flux_solutions: DictList of mass.MassSolution, MassSolution
            A DictList of MassSolutions containing the flux solutions for
            successfully finding steady state. If no steady state is found, the
            MassSolution will be returned as empty. If there is only one model,
            a single MassSolution is returned instead of the DictList.

        """
        if kwargs:
            # TODO parse kwarg arguments
            pass

        # Set all models for simulation if None provided.
        if models is None:
            models = self.models
        models = [str(model) for model in ensure_iterable(models)]

        # Get roadrunner and executable model instances
        rr = self.roadrunner
        # Ensure strategy input is valid
        if strategy not in rr.getRegisteredSteadyStateSolverNames()\
           and strategy != "simulate":
            raise ValueError(
                "Invalid steady state strategy: '{0}'".format(strategy))

        # Parse the perturbations and format the input to be used
        perturbations = self._format_perturbations_input(perturbations)
        if strategy == "simulate":
            steady_state_function = self._find_steady_state_simulate
        else:
            steady_state_function = self._find_steady_state_solver
            rr.setSteadyStateSolver(strategy)

        # Set species to use for steady state calculations
        selections = self._make_rr_selections(
            selections=kwargs.get("selections", None),
            boundary_metabolites=kwargs.get("boundary_metabolites", False))
        # Remove time from selections
        selections.remove("time")

        # Make DictLists for solution objects
        conc_sol_list = DictList()
        flux_sol_list = DictList()
        try:
            for model in models:
                try:
                    # Apply perturbations and set model values in roadrunner
                    rr = self._set_simulation_values(model, perturbations)
                    # Use simulate strategy
                    rr.steadyStateSelections = selections
                    results = steady_state_function(model, **kwargs)
                    # Map results to their identifiers and return MassSolutions
                    solutions = self._make_mass_solutions(
                        model, default_selections=selections, results=results,
                        **kwargs)
                    if update_values:
                        self._update_model_values(model, solutions)

                # Handle MassSimulationErrors
                except MassSimulationError as e:
                    LOGGER.log(
                        LOGGER.LOG_ERROR, "Unable to find a steady state for "
                        "Model '{0}' using strategy '{1}' due to the following"
                        " MassSimulationError: {2}".format(
                            model, strategy, str(e)))
                    # Make empty MassSolutions
                    solutions = self._make_mass_solutions(
                        model, default_selections=selections, results=results,
                        **kwargs)
                finally:
                    # Add solutions to output lists
                    conc_sol_list += [solutions[0]]
                    flux_sol_list += [solutions[1]]
                    # Reset the roadrunner state
                    rr = _reset_roadrunner_instance(rr)

        # Handle unforseen errors as critical errors
        except Exception as e:
            # If an error occurs, raise it as a MassSimulationError
            raise MassSimulationError(
                "Critical simulation fail due to the following:\n"
                + str(e.__class__.__name__) + ": " + str(e))

        # Return just a tuple of two MassSolutions if only one model simulated
        if len(models) == 1:
            return conc_sol_list[0], flux_sol_list[0]

        return conc_sol_list, flux_sol_list

    def _make_rr_selections(self, selections=None, boundary_metabolites=None):
        """Set the observable output of the simulation.

        Warnings
        --------
        This method is intended for internal use only.

        """
        # Get roadrunner executable model instance
        rr_model = self.roadrunner.model

        rr_selections = ["time"]
        if selections is None:
            rr_selections += rr_model.getFloatingSpeciesConcentrationIds()
            if boundary_metabolites:
                rr_selections += rr_model.getBoundarySpeciesConcentrationIds()
            rr_selections += rr_model.getReactionIds()
        else:
            # TODO Catch selection errors here
            rr_selections += selections
            pass

        return rr_selections

    def _format_perturbations_input(self, perturbations):
        """Check and format the perturbation input.

        Perturbations are checked before simulations are carried out to limit
        fails during the simulation due to bad syntax or values.

        Warnings
        --------
        This method is intended for internal use only.

        """
        formatted_perturbations = {}
        # Check input type.
        if not isinstance(perturbations, dict) and perturbations is not None:
            raise TypeError("Perturbations must be a dict")

        if perturbations:
            def validate_perturbation(key, value, pattern_re, item_dict):
                """Validate the perturbation for the given pattern."""
                perturbation = {}
                # Check if the perturbation matches the pattern.
                match = pattern_re.match(key)
                if match:
                    item = match.group()
                    # Try to convert the perturbation to a value, and ensure
                    # perturbation formula is valid.
                    if item in item_dict:
                        value = _convert_perturbation_str_to_number(key, value,
                                                                    pattern_re)
                    else:
                        # Raise ValueError if perturbed variable not found.
                        raise ValueError("Invalid Perturbation: '{0}' not "
                                         "found in model values".format(item))
                    perturbation[key] = value

                return perturbation

            # Get the reference model initial conditions and parameters
            inits, params = self.get_model_values(self.reference_model)
            # Iterate through perturbations to ensure they are valid.
            for key, value in iteritems(perturbations):
                for pert_type, pattern_re in iteritems(PERTURBATIONS_RE_DICT):
                    if pert_type in ["kf", "Keq", "kr"]:
                        args = (key, value, pattern_re, params)
                    elif pert_type in ["fixed", "func"]:
                        # TODO fixed concentrations and functions
                        pass

                    else:
                        args = (key, value, pattern_re, inits)
                    perturbation = validate_perturbation(*args)
                    if not perturbation:
                        continue

                    formatted_perturbations.update(perturbation)
                    break

        return formatted_perturbations

    def _set_simulation_values(self, model, perturbations):
        """Set the simulation numerical values in the roadrunner instance.

        Warnings
        --------
        This method is intended for internal use only.

        """
        # Get roadrunner
        rr = self.roadrunner
        try:
            # Ensure values exist for the model
            init_conds, parameters = self.get_model_values(model)
        except ValueError as e:
            raise MassSimulationError(e)

        # Gather model values in one dict
        model_values_to_set = {}
        model_values_to_set.update(init_conds)
        model_values_to_set.update(parameters)

        # Apply perturbations to model values.
        try:
            for key, value in iteritems(perturbations):
                # Perturb value to a number if value is float
                if isinstance(value, float):
                    model_values_to_set[key] = value
                    continue

                # Perturb value to scale current variable value
                if not isinstance(value, float) and key in value:
                    # Determine perturbation type and regex pattern
                    for k, v in iteritems(PERTURBATIONS_RE_DICT):
                        if v.search(value):
                            p_type, pattern_re = (k, v)
                            break
                    if p_type in ["kf", "Keq", "kr"]:
                        value = pattern_re.sub(str(parameters[key]), value)
                    elif p_type in ["fixed", "func"]:
                        # TODO fixed concentrations and functions
                        pass

                    else:
                        value = pattern_re.sub(str(init_conds[key]), value)
                    # Convert perturbation to a numerical value for simulation
                    value = sympify(value)
                    value = _convert_perturbation_str_to_number(key, value,
                                                                pattern_re)
                    model_values_to_set[key] = value
                    continue
                # Raise error if perturbation not applied.
                raise MassSimulationError(
                    "Could not apply perturbation {{'{0}': '{1}'}}".format(
                        key, value))
        except (ValueError, MassSimulationError) as e:
            raise MassSimulationError(e)

        # Set roadrunner to reflect given model variant
        for key, value in iteritems(model_values_to_set):
            rr.setValue(key, value)

        # Return the roadrunner instance
        return rr

    def _make_mass_solutions(self, model, default_selections, results,
                             **kwargs):
        """Make the MassSolutions using the results of the Simulation.

        Warnings
        --------
        This method is intended for internal use only.

        """
        selections = kwargs.get("selections", default_selections)
        interpolate = kwargs.get("interpolate", False)
        # Get roadrunner and roadrunner executable model instances
        rr = self.roadrunner
        rr_model = rr.model

        # Create dicts for containing concentrations and flux solutions.
        time = None
        conc_sol = {}
        flux_sol = {}
        if results is not None:
            # Get species and reaction lists
            species_list = rr_model.getFloatingSpeciesConcentrationIds()
            if kwargs.get("boundary_metabolites", False):
                species_list += rr_model.getBoundarySpeciesConcentrationIds()
            reaction_list = rr_model.getReactionIds()

            for sol_key in selections:
                # Get the time vector
                if sol_key == "time":
                    time = results[sol_key]
                elif sol_key in species_list:
                    # Add solution to concentration solutions if specie
                    result_array = _round_values(results[sol_key])
                    conc_sol[_strip_conc(sol_key)] = result_array
                elif sol_key in reaction_list:
                    # Add solution to flux solutions if reaction
                    result_array = _round_values(results[sol_key])
                    flux_sol[sol_key] = result_array
                else:
                    pass

        # Make a MassSolution object of the concentration solution dict.
        conc_sol = MassSolution(id_or_model=model, data_dict=conc_sol,
                                solution_type=_CONC_STR, time=time,
                                interpolate=interpolate)

        # Make a MassSolution object of the concentration solution dict.
        flux_sol = MassSolution(id_or_model=model, data_dict=flux_sol,
                                solution_type=_FLUX_STR, time=time,
                                interpolate=interpolate)

        return conc_sol, flux_sol

    def _update_stored_solutions(self, solution_type, solutions):
        """Update stored MassSolutions with new MassSolution objects.

        Warnings
        --------
        This method is intended for internal use only.

        """
        # Determine where solutions are stored
        stored_solution_dictlist = {
            "Conc": self._conc_solutions,
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
        steps = kwargs.get("steps", None)
        if steps is None:
            time = _format_time_input((0, kwargs.get("tfinal", 1e8), 1001))
        else:
            time = _format_time_input((0, kwargs.get("tfinal", 1e8)),
                                      steps=steps)
        # Simulate for a long time
        simulation_results = rr.simulate(*time)

        def is_steady_state(abs_diff):
            """Compare the absolute diff. to the steady state threshold."""
            # Round value before comparision
            if MASSCONFIGURATION.decimal_precision is not None:
                abs_diff = round(abs_diff, MASSCONFIGURATION.decimal_precision)
            if abs_diff <= MASSCONFIGURATION.steady_state_threshold:
                return True
            return False

        ss_results = {}
        for sol_key in rr.timeCourseSelections:
            second_to_last, last = simulation_results[sol_key][-2:]
            success = is_steady_state(abs(second_to_last - last))
            if success:
                ss_results[sol_key] = last
            else:
                # If no steady state found, raise MassSimulationError
                raise MassSimulationError(
                    "Absolute difference for '{0}' in Model '{1}' is greater "
                    "than the steady state threshold.".format(sol_key, model))

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
            is_ss_value = rr.steadyState()
            # Round value before comparision
            if MASSCONFIGURATION.decimal_precision is not None:
                is_ss_value = round(is_ss_value,
                                    MASSCONFIGURATION.decimal_precision)
            if is_ss_value <= MASSCONFIGURATION.steady_state_threshold:
                return True
            return False

        # See if a steady state can be reached
        i = 0
        success = False
        num_attempts = kwargs.get("num_attempts", 2)
        while not success and i < num_attempts:
            success = is_steady_state(rr)
            i += 1
        
        if not success:
            raise MassSimulationError(
                "Could not find steady state for Model '{0}' after {1:d} "
                "attempts, steady state threshold always exceeded.".format(
                    model, i))

        # Zip the solutions with their IDs into a dict and return
        results = zip(rr.steadyStateSelections, rr.getSteadyStateValues())
        return dict(results)

    def _update_model_values(self, model, values):
        """Update the stored values for the model with the new ones.

        Warnings
        --------
        This method is intended for internal use only.

        """
        model_values = self.get_model_values(model)
        id_fix_dict = {model + "_init_conds": _make_init_cond,
                       model + "_parameters": _make_ss_flux}
        for current_value_dict, new_value_dict in zip(model_values, values):
            id_fix_func = id_fix_dict[current_value_dict.id]
            for key, value in iteritems(new_value_dict):
                current_value_dict[id_fix_func(key)] = value


def _assess_model_quality_for_simulation(mass_model, verbose):
    """TODO DOCSTRING."""
    # TODO utilize via QCQA class

    pass


def _create_roadrunner_instance(mass_model, f_replace=None, verbose=False,
                                **kwargs):
    """Create a RoadRunner instance for the given model.

    Warnings
    --------
    This method is intended for internal use only.

    """
    doc = _model_to_sbml(mass_model, f_replace=f_replace, **kwargs)
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

    if verbose:
        print("Successfully loaded Model '" + str(mass_model) + "'")

    return roadrunner.RoadRunner(sbml_str)


def _reset_roadrunner_instance(roadrunner):
    """Reset a Simulation object's roadrunner instance.

    Warnings
    --------
    This method is intended for internal use only.

    """
    roadrunner.resetToOrigin()
    return roadrunner


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
            id=mass_model.id + "_init_conds", data_dict={
                _make_init_cond(met.id): ic
                for met, ic in iteritems(init_conds)}),
        "parameters": DictWithID(
            id=mass_model.id + "_parameters", data_dict={
                param: value
                for param, value in iteritems(parameters)}),
    }

    return DictWithID(id=mass_model.id + "_values", data_dict=values)


def _format_time_input(time, steps=None):
    """Format the time input and return it as a tuple for the RoadRunner.

    Warnings
    --------
    This method is intended for internal use only.

    """
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


def _convert_perturbation_str_to_number(key, value, pattern_re):
    """Attempt to convert perturbation value string to a numerical value.

    Warnings
    --------
    This method is intended for internal use only.

    """
    try:
        value = float(value)
    except (ValueError, TypeError):
        match = pattern_re.search(str(value))
        if not match or match.group().strip() != key:
            raise ValueError(
                "Invalid Perturbation {{'{0}': '{1}'}}".format(
                    key, value))
    return value


def _round_values(values):
    """Round a value to the decimal precision in the MassConfiguration.

    Warnings
    --------
    This method is intended for internal use only.

    """
    return_as_array = True

    if not hasattr(values, "__iter__"):
        values = ensure_iterable(values)
        return_as_array = False

    if MASSCONFIGURATION.decimal_precision is not None:
        values = np.array([
            round(v, MASSCONFIGURATION.decimal_precision) for v in values])

    if not return_as_array:
        values = values[0]

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

