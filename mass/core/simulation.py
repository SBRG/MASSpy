# -*- coding: utf-8 -*-
"""The Simulation module addresses the simulation of mass.MassModels.

The Simulation module is designed to address all aspects related to the
simulation of one or more MassModel objects. These aspects include storing the
model(s) to be simulated, the numerical values of associated parameters and
initial concentrations for quick access, and handling simulation results.
Perturbations can also be implemented in simulations as long as they follow the
following guidelines:

Valid keys for perturbation dictionary can include the following:

    "RxnID.forward_rate_constant" or "RxnID.kf"
    "RxnID.equilibrium_constant" or "RxnID.Keq"
    "RxnID.reverse_rate_constant" or "RxnID.kr"
    "MetID.initial_condition" or "MetID.ic"
    "MetID.fixed" or "MetID.fix"
    "MetID.function" (Only valid for exchange "boundary" metabolites)

where RxnID represents the reaction identifier and MetID represents the
metabolite identifier. To perform a perturbation where the current value is
increased or decreased by a scalar, placing "[key]" into the corresponding
value. Some examples of valid perturbations include:

    perturbation_dictionary = {
        "RxnID.kf": "[RxnID.kf]*1.5",
        "RxnID.equilibrium_constant": 5.51,
        "MetID.initial_condition": "[MetID.initial_condition]" + 11,
        "MetID.fixed": 0.0045,
        "MetID.fixed": "[MetID.ic]" / 10,
        "MetID.function" "(70 + 30*sin(120*pi*t))*2.8684*1e-4"
        }

Note that functions defined for any "MetID.function" perturbations must be
able to be interpreted by sympy.sympify. Additionally, a parameter or an
initial condition cannot be perturbed twice in the same simulation.

All simulation results are returned as mass.Solution objects. Solution objects
can behave as booleans, with Solution objects containining an empty solution
dictionary returning as False, and those with solutions returning as True.
"""
from __future__ import absolute_import

import re
from collections import OrderedDict
from copy import deepcopy
from math import log10
from warnings import filterwarnings, warn

import numpy as np

from scipy.integrate import solve_ivp
from scipy.optimize import root

from six import iteritems, iterkeys, itervalues, string_types

import sympy as sym

from cobra.core.dictlist import DictList
from cobra.core.object import Object

from mass.core.massmodel import MassModel
from mass.core.masssolution import (
    MassSolution, _CONC_STR, _FLUX_STR, _NETFLUX_STR, _POOL_STR)
from mass.exceptions import MassSimulationError
from mass.util.DictWithID import DictWithID
from mass.util.expressions import _mk_met_func, strip_time
from mass.util.qcqa import is_simulatable, qcqa_model, qcqa_simulation
from mass.util.util import ensure_iterable


_ZERO_TOL = 1e-8
# Pre-compiled regular expressions for perturbations
_KF_RE = re.compile("forward_rate_constant|kf")
_KEQ_RE = re.compile("equilibrium_constant|Keq")
_KR_RE = re.compile("reverse_rate_constant|kr")
_IC_RE = re.compile("initial_condition|ic")
_FIX_RE = re.compile("fix|fixed")
_FUNC_RE = re.compile("func|function")
_CUSTOM_RE = re.compile("custom")
# Global symbol for time
_T_SYM = sym.Symbol("t")
_LAMBDIFY_MODULE = ["numpy"]
# Define default option dicts for solvers

_SCIPY_DEFAULT_OPTIONS = DictWithID(
    id="scipy", dictionary={
        "method": "LSODA",
        "dense_output": True,
        "t_eval": None,
        "events": None,
        "max_step": np.inf,
        "rtol": 1e-6,
        "atol": 1e-9,
        "jac_sparsity": None,
        "lband": None,
        "uband": None,
        "min_step": 0.,
        "first_step": None,
    })
_ALL_DEFAULT_OPTIONS = DictList([_SCIPY_DEFAULT_OPTIONS])
try:
    from scikits.odes import ode as _sundials_wrapper
    _SUNDIALS_DEFAULT_OPTIONS = DictWithID(
        id="SUNDIALS", dictionary={
            
        })
    _ALL_DEFAULT_OPTIONS.add(_SUNDIALS_DEFAULT_OPTIONS)
except ImportError:
    pass

filterwarnings(action="ignore", module="scipy", message="^internal gelsd")


class Simulation(Object):
    """Class for managing setup and results of simulations for MassModels.

    The mass.Simulation class is designed to address all aspects related to the
    simulation of MassModel objects, some of which includes setting solver
    and solver options, perturbation of concentrations and parameters,
    simulation of a single MassModel or an ensemble of models, and handling
    of simulation results.

    Parameters
    ----------
    models: mass.MassModel, iterable of mass.MassModels
        A mass.MassModel or an iterable of mass.MassModels to load into
        the simulation object. 
    reference_model: mass.MassModel, optional
        The MassModel object to be considered as the main reference model for
        the simulation. If None provided, the first model of the provided
        iterable of models will be used.
    simulation_id: str, optional
        An identifier to associate with the Simulation given as a string. If
        None provided, will default to "(MassModel.id)_Simulation."
    name: str, optional
        A human readable name for the Simulation.

    Attributes
    ----------
    solver: str
        Return the solver that is currently set for the Simulation. Note that
        "solver" refers to the solver package to use (e.g. SUNDIALS, scipy), 
        and the "method" refers to the integration algorithm used. 
        (e.g. CVODE, RK45) A full list of possible solvers can be viewed using
        the Simulation,all_solvers property.
    models: cobra.DictList
        A cobra.DictList where the keys are the model identifiers and the
        values are the associated mass.MassModel objects.
    reference_model: mass.MassModel
        The MassModel to be the reference. Must exist in the Simulation models.
    units: dict
        A dictionary to store the units used in the simulation for referencing.
        Example: {'N': 'Millimoles', 'Vol': 'Liters', 'Time': 'Hours'}

    Warnings
    --------
    The Simulation object will not automatically track or convert units.
        Therefore, it is up to the user to ensure unit consistency. The
        Simulation.units attribute is provided as a way to inform the user or
        others what units are used in the simulation results.

    """

    def __init__(self, models, reference_model=None, simulation_id=None, 
                 name=None, solver="scipy"):
        """Initialize the Simulation Object."""
        models = ensure_iterable(models)
        for model in models:
            if not isinstance(model, MassModel):
                raise TypeError("'{0}' not a valid MassModel or subclass of "
                                "MassModel.".format(str(model)))
        if reference_model is None:
            reference_model = models[0]
        if not isinstance(reference_model, MassModel):
            raise TypeError("reference_model must be a mass.MassModel")

        if simulation_id is None:
            simulation_id = "{0}_Simulation".format(reference_model.id)
        Object.__init__(self, simulation_id, name)
        self._reference_model = reference_model
        self._models = DictList(models)

        # Get initial condition and parameter values, and store.
        all_param_vals = list(map(self._get_parameters_from_model, models))
        all_ic_vals = list(map(self._get_ics_from_model, models))

        self._values = DictList() + all_param_vals + all_ic_vals

        self._solutions = DictList()
        self._all_solver_option_dicts = deepcopy(_ALL_DEFAULT_OPTIONS)
        self._solver = solver

    # Public
    @property
    def reference_model(self):
        """Return the reference Massmodel for the Simulation."""
        return getattr(self, "_reference_model", None)

    @reference_model.setter
    def reference_model(self, value):
        """Set the reference MassModel for the Simulation.

        Notes
        -----
        If attempting to set a reference model that is not currently in
            self.models, the model will also added to self.models.

        """
        if not isinstance(value, MassModel):
            raise TypeError("reference_model must be a mass.MassModel")

        if value not in self.models:
            self.add_models(value)

        self._reference_model = value

    @property
    def solver(self):
        """Return the solver that is currently set for the Simulation."""
        return getattr(self, "_solver", None)

    @solver.setter
    def solver(self, value):
        """Set the solver in Simulation to use in simulations."""
        if value not in self.all_solvers:
            raise ValueError("solver '{0}' not recognized. Acceptable solvers "
                             "include the following: {1!r}"
                             .format(value, self.all_solvers))
        else:
            setattr(self, "_solver", value)

    @property
    def all_solvers(self):
        """Return a list of possible solvers that can be used."""
        return sorted(set([sol_opts.id for sol_opts in _ALL_DEFAULT_OPTIONS]))

    @property
    def models(self):
        """Return a DictList of models associated with this Simulation."""
        return getattr(self, "_models", None)

    def add_models(self, models):
        """Add a list of models and associated values to the Simulation.

        Parameters
        ----------
        models: mass.MassModel, list of mass.MassModels,
            The model or models to add to the Simulation.

        Notes
        -----
        If a model already exists in the Simulation, its associated numerical 
            values will be updated instead.

        """
        models = ensure_iterable(models)
        # Check whether a model is a MassModel object, then check if
        # the model already exists in the Simulation, ignoring those that do.
        for model in models:
            if not isinstance(model, MassModel):
                warn("Skipping {0}, not a MassModel object".format(model.id))
                models.remove(model)

            if model in self.models:
                warn("Updating {0} as it already exists in the Simulation "
                     "object.".format(model.id))
                self.update_values(model)
                models.remove(model)

        for model in models:
            self._values.add(self._get_parameters_from_model(model))
            self._values.add(self._get_ics_from_model(model))
            self.models.add(model)

    def remove_models(self, models):
        """Remove a list of models and associated values from the Simulation.

        Parameters
        ----------
        models: mass.MassModel, list of mass.MassModels
            The model or models to remove from the Simulation. The reference
            model cannot be removed.

        Notes
        -----
        The reference model will not be removed. To remove the reference model,
            first change the reference model.

        """
        models = ensure_iterable(models)
        # Check whether a model already exists in the Simulation,  
        # ignoring those that do not.
        model_list = []
        for model in ensure_iterable(models):
            try:
                model_list.append(self.models.get_by_id(str(model)))
            except KeyError:
                pass

        for model in model_list:
            if model != self.reference_model:
                for value_type in ["_parameters", "_ics"]:
                    self._values.remove(model.id + value_type)
                self.models.remove(model)

    def get_concentration_solutions(self, models=None):
        """Return the Conc. Solutions for a list of models.

        Parameters
        ----------
        models: mass.MassModel, list of mass.MassModels, None
            The model or models to lookup. If None provided, all models in
            the Simulation object will be used. Otherwise, any provided
            models must already be present in the Simulation object.

        Returns
        -------
        solutions: mass.Solution, dict
            Either a Solution object if there is only one model, or a dict
            with model identifiers as keys and the corresponding Solution
            objects as values for multiple models.

        """
        return self._lookup_solutions(models, _CONC_STR)

    def get_flux_solutions(self, models=None):
        """Return the Flux Solutions for a list of models.

        Models must already exist in the Simulation object.

        Parameters
        ----------
        models: mass.MassModel, list of mass.MassModels, None
            The model or models to lookup. If None provided, all models in
            the Simulation object will be used. Otherwise, any provided
            models must already be present in the Simulation object.

        Returns
        -------
        solutions: mass.Solution, dict
            Either a Solution object if there is only one model, or a dict
            with model identifiers as keys and the corresponding Solution
            objects as values for multiple models.

        """
        return self._lookup_solutions(models, _FLUX_STR)

    def get_pool_solutions(self, models=None):
        """Return the Pool Solutions for a list of models.

        Parameters
        ----------
        models: mass.MassModel, list of mass.MassModels, None
            The model or models to lookup. If None provided, all models in
            the Simulation object will be used. Otherwise, any provided
            models must already be present in the Simulation object.

        Returns
        -------
        solutions: mass.Solution, dict
            Either a Solution object if there is only one model, or a dict
            with model identifiers as keys and the corresponding Solution
            objects as values for multiple models.

        """
        return self._lookup_solutions(models, _POOL_STR)

    def get_net_flux_solutions(self, models=None):
        """Return the NetFlux Solutions for a list of models.

        Parameters
        ----------
        models: mass.MassModel, list of mass.MassModels, None
            The model or models to lookup. If None provided, all models in
            the Simulation object will be used. Otherwise, any provided
            models must already be present in the Simulation object.

        Returns
        -------
        solutions: mass.Solution, dict
            Either a Solution object if there is only one model, or a dict
            with model identifiers as keys and the corresponding Solution
            objects as values for multiple models.

        """
        return self._lookup_solutions(models, _NETFLUX_STR)

    def view_model_values(self, model):
        """Return copies of stored numerical values associated with a model.

        Parameters
        ----------
        model: mass.MassModel
            The model to lookup in the Simulation object.

        Return
        ------
        parameters: dict
            A copy of the stored parameters for the given model.
        initial_conditions: dict
            A copy of the stored initial comnditions for the given model.

        """
        # Check for the presence of the model in the Simulation
        self._check_for_one_model(model)
        parameters, initial_conditions = self._values.get_by_any(
            [model.id + "_parameters", model.id + "_ics"])
        return parameters.copy(), initial_conditions.copy()

    def view_parameter_values(self, models=None):
        """Return a copy of stored parameters for a list of models.

        Parameters
        ----------
        models: mass.MassModel, list of mass.MassModels, None
            The model or models to lookup. If None provided, all models in
            the Simulation object will be used. Otherwise, any provided
            models must already be present in the Simulation object.

        Returns
        -------
        value_dict: dict
            A copy of the stored dictionary with model identifiers as keys
            and the corresponding parameters as values.

        """
        # Check for models if a list provided, remove models from the list
        # that are not currently in the Simulation
        models = self._get_models_for_method(models, False)
        parameters = self._values.get_by_any([model.id + "_parameters"
                                              for model in models])
        value_dict = {model.id: parameter
                      for model, parameter, in zip(models, parameters)}
        return value_dict.copy()

    def view_initial_concentration_values(self, models=None):
        """Return a copy of stored initial concentrations for a list of models.

        Parameters
        ----------
        models: mass.MassModel, list of mass.MassModels, None
            The model or models to lookup. If None provided, all models in
            the Simulation object will be used. Otherwise, any provided
            models must already be present in the Simulation object.

        Returns
        -------
        value_dict: dict
            A copy of the stored dictionary with model identifiers as keys
            and the corresponding initial concentrations as values.

        """
        # Check for models if a list provided, remove models from the list
        # that are not currently in the Simulation
        models = self._get_models_for_method(models, False)
        ics = self._values.get_by_any([model.id + "_ics" for model in models])
        value_dict = {model.id: ic for model, ic, in zip(models, ics)}
        return value_dict.copy()

    def get_solver_options(self, solver=None):
        """Return a copy of the current solver options for a given solver.

        Parameters
        ----------
        solver: str, optional:
            If provided, will return the options for the given solver.
            Otherwise return the options for the current solver.

        Return
        ------
        options: dict
            A dictionary with the solver options.

        """
        if not solver:
            solver = self.solver
        if solver not in self.all_solvers:
            raise ValueError("solver '{0}' not recognized. Acceptable solvers "
                             "include the following: {1!r}"
                             .format(solver, self.all_solvers))

        return self._all_solver_option_dicts.get_by_id(solver).copy()

    def set_solver_options(self, solver=None, **options):
        """Set the current default solver options for a given solver.

        Parameters
        ----------
        solver: str, optional:
           The solver that will have its default use options updated. If None
           provided, will update the current solver.
        **options: 
            The options that correspond to the solver. To view possible solver
            options, use the Simulation.get_solver_options method.

        """
        if solver is None:
            solver = self.solver
        if solver not in self.all_solvers:
            raise ValueError("solver '{0}' not recognized. Acceptable solvers "
                             "include the following: {1!r}"
                             .format(solver, self.all_solvers))
        current_options = self._all_solver_option_dicts.get_by_id(solver)
        for key, value in iteritems(options):
            if key in current_options:
                current_options[key] = value
            else:
                warn("Unrecognized solver option: '{0}'.".format(str(key)))

    def reset_solver_options(self, solver=None):
        """Reet the current solver options for a given solver.

        Parameters
        ----------
        solvers: str:
           The solver that will have its default options reset. If None 
           provided, the current solver will have its default options reset.

        """
        if solver is None:
            solver = self.solver
        if solver not in self.all_solvers:
            raise ValueError("solver '{0}' not recognized. Acceptable solvers "
                             "include the following: {1!r}"
                             .format(solver, self.all_solvers))

        # Remove old solver dict and replace with a copy of the new one. 
        self._all_solver_option_dicts.remove(solver)
        new_options = deepcopy(_ALL_DEFAULT_OPTIONS.get_by_id(solver))
        self._all_solver_option_dicts.add(new_options)

        print("Default options for solver '{0}' reset".format(solver))

    def simulate_model(self, model, time=None, perturbations=None,
                       interpolate=True, verbose=True, return_obj=False,
                       update_solutions=True, **options):
        """Simulate a single MassModel and return the solution profiles.

        The "model" is simulated by integrating the ODEs using the solver
        set in self.solver to compute the solution at time points between the
        initial and final time points given in "time" while incorporating
        events specified in the "perturbations." See Simulation Module
        documentation for more information.

        Parameters
        ----------
        model: mass.MassModel or str
            The MassModel or the string identifier of the MassModel to
            simulate.
        time: tuple of floats
            A tuple containing the start and end time points of the simulation.
        perturbations: dict, optional
            A dict of events to incorporate into the simulation, where keys are
            the object and type of event to incorporate and values are the
            change to be made.
        interpolate: bool, optional
            If True, then solutions in both Solution objects are scipy
            interpolating functions, otherwise solutions are arrays of solution
            values. Default is True
        verbose: bool, optional
            If True, print a detailed report of why the simulation failed.
            Default is True.
        return_obj: bool, optional
            If True, then the original solution object is also returned in
            addition to the mass.Solution objects.
        update_solutions: bool, optional
            If True, then the mass.Solution objects stored within
            self.concentration_solutions and self.flux_solutions are updated.
        **options
            Options for the current solver to be used. Options must correspond
            to the set solver, otherwise they will be ignored.

        Returns
        -------
        conc_sol: mass.Solution
            A mass.Solution object containing the dict of concentrations
            solutions for successful simulations, and an empty dict for failed
            simulations.
        flux_sol: mass.Solution
            A mass.Solution object containing the dict of flux solutions for
            successful simulations, and an empty dict for failed simulations.
        sol_obj: object
            The original solution object returned by the solver directly after
            integration. Only returned if ``return_obj=True``.

        """
        if time is None or len(time) != 2 or isinstance(time, string_types):
            raise ValueError("time must be a tuple of 2 floats.")
        args = zip(["interpolate", "verbose", "return_obj", update_solutions],
                   [interpolate, verbose, return_obj, update_solutions])
        for arg, val in args:
            if not isinstance(val, bool):
                raise TypeError("'{0}' must be a bool. Found '{1}' instead."
                                .format(arg, val))

        # Helper function to update, store, and return solutions
        def _update_and_return(csol, fsol, sol_obj, update, return_obj):
            # Store the solutions if desired.
            if update:
                self._update_solution_storage([csol, fsol])
            # Return solutions and return solution object if desired.
            if return_obj:
                return csol, fsol, sol_obj
            else:
                return csol, fsol

        try:
            # Check for the presence of the model in the Simulation
            self._check_for_one_model(model, verbose)
            # Check whether a model can be simulated.
            self._assess_quality(model, verbose)
            # Set appropriate rate laws for each reaction
            self._set_rate_type(model)
            # Implement any perturbations and get the value substituion dicts.
            parameters, ics, functions = self._apply_perturbations(
                model, perturbations)
            # Setup ODEs for integration
            [odes, jacb, ics] = self._make_odes_functions(model, parameters,
                                                          ics, functions)
            # Integrate ODEs to obtain concentration solutions
            t, conc_sol, sol_obj = self._integrate_odes(time, odes, jacb,
                                                        ics, options)
            # Create Solution object for concentrations
            conc_sol = MassSolution(model, _CONC_STR, conc_sol, t)

            # Calculate flux solutions using the concentration solutions
            flux_sol = self._calculate_flux_solutions(model, parameters,
                                                      conc_sol, t)

            # Turn solutions into interpolating functions if desired
            if interpolate:
                conc_sol.interpolate = True
                flux_sol.interpolate = True

        # Return empty Solution objects if simulation cannot proceed.
        except MassSimulationError as e:
            warn(str(e))
            # Create empty solution objects
            conc_sol = MassSolution(model, _CONC_STR)
            flux_sol = MassSolution(model, _FLUX_STR)
            sol_obj = None

        return _update_and_return(conc_sol, flux_sol, sol_obj,
                                  update_solutions, return_obj)

    def simulate(self, models=None, time=None, perturbations=None,
                 interpolate=True, verbose=True, return_obj=False,
                 update_solutions=True, **options):
        """Simulate MassModel(s) and return the solution profiles.

        Each model provided in "models" is simulated by integrating the ODEs
        using the solver set in self.solver to compute the solution at time
        points between the initial and final time points given in "time" while
        incorporating events specified in the "perturbations." See Simulation
        Module documentation for more information.

        Parameters
        ----------
        models: mass.MassModel, list of mass.MassModels, None
            The model or models to simulate. If None provided, all models in
            the Simulation object will be simulated. Otherwise, any provided
            models must already be present in the Simulation object.
        time: tuple of floats
            A tuple containing the start and end time points of the simulation.
        perturbations: dict, optional
            A dict of events to incorporate into the simulation, where keys are
            the object and type of event to incorporate and values are the
            change to be made.
        interpolate: bool, optional
            If True, then solutions in both Solution objects are scipy
            interpolating functions, otherwise solutions are arrays of solution
            values. Default is True
        verbose: bool, optional
            If True, print a detailed report of why the simulation failed.
            Default is True.
        return_obj: bool, optional
            If True, then the original solution object is also returned in
            addition to the mass.Solution objects.
        update_solutions: bool, optional
            If True, then the mass.Solution objects stored within
            self.concentration_solutions and self.flux_solutions are updated.
        **options
            Options for the current solver to be used. Options must correspond
            to the set solver, otherwise they will be ignored.

        Returns
        -------
        conc_solutions: DictList of mass.Solution objects
            A DictList of mass.Solution objects, each containing the dict of
            concentrations solutions for successful simulations, and an
            empty dict for failed simulations.
        flux_solutions: DictList of mass.Solution objects
            A DictList of mass.Solution objects, each containing the dict of
            flux solutions for successful simulations, and an empty dict for
            failed simulations.
        sol_objects: dict of objects
            A dict of the original solution objects returned by the solver,
            where keys are model identifiers and values are solution objects.
            Only returned if ``return_obj=True``.

        Notes
        -----
        A mass.Solution object will contain a dict of solutions for successful
            simulations, and an empty dict for failed simulations.
            Empty solution objects also return as a bool False.

        """
        # Check for models if a list provided, remove models from the list
        # that are not currently in the Simulation
        models = self._get_models_for_method(models, verbose)
        conc_solutions = DictList()
        flux_solutions = DictList()
        sol_objects = {}

        # Simulate models that are present.
        for model in models:
            if verbose:
                print("Simulating " + str(model) + "...")
            sols = self.simulate_model(model=model, time=time,
                                       perturbations=perturbations,
                                       interpolate=interpolate,
                                       verbose=verbose, return_obj=return_obj,
                                       update_solutions=False,
                                       **options)
            sols = list(sols)
            if return_obj:
                sol_obj = sols.pop(2)
                sol_objects[model.id] = sol_obj

            if update_solutions:
                self._update_solution_storage(sols)
            for sol, storage in zip(sols, [conc_solutions, flux_solutions]):
                storage += [sol]

        if return_obj:
            return conc_solutions, flux_solutions, sol_objects
        else:
            return conc_solutions, flux_solutions

    def find_steady_state_model(self, model, strategy="simulate",
                                perturbations=None, verbose=True,
                                update_initial_conditions=False,
                                update_reactions=False, **options):
        """Find the steady state for a single model using the given strategy.

        Parameters
        ----------
        model: mass.MassModel or str
            The MassModel or the string identifier of the MassModel to
            simulate.
        strategy : 'simulate' or 'roots', tuple of len 2
            A string representing the strategy to use to solve for the 
            steady state. The strategy string can be 'simulate' to simulate the 
            model to steady state or it can be 'roots' to calculate the roots. 
            Default method for roots is lm.
        perturbations: dict, optional
            A dict of events to incorporate into the simulation, where keys are
            the object and type of event to incorporate and values are the
            change to be made.
        verbose: bool, optional
            If True, print a detailed report of why the method failed.
            Default is True.
        update_reactions : bool, optional
            If True, update the steady state flux attribute in each reaction.
        update_initial_conditions : bool, optional
            If True, update the model initial conditions for each metabolite,
            and update the Simulation.
        **options
            Options for the solver. Options must correspond to the set solver
            for the given stratgey, otherwise they will be ignored.

        Returns
        -------
        conc_sol: mass.Solution
            A mass.Solution object containing the dict of concentrations
            solutions for successful simulations, and an empty dict for failed
            simulations.
        flux_sol: mass.Solution
            A mass.Solution object containing the dict of flux solutions for
            successful simulations, and an empty dict for failed simulations.

        """
        # Check for the presence of the model in the Simulation
        self._check_for_one_model(model, verbose)
        strategy_dict = {"simulate": self._find_steady_state_simulate,
                         "roots": self._find_steady_state_roots}
        if strategy not in strategy_dict:
            raise ValueError("Unrecognized strategy. strategy must be "
                             "'simulate' or 'roots'.")

        chop = abs(int(log10(_ZERO_TOL)))
        update = (update_initial_conditions, update_reactions)
        try:
            # Check Simulation to determine whether simulation can proceed.
            self._assess_quality(model, verbose)
            conc_sol, flux_sol = strategy_dict[strategy](
                model, perturbations, verbose, chop, update, **options)
            conc_sol = MassSolution(
                model, _CONC_STR, solution_dictionary=conc_sol)
            flux_sol = MassSolution(
                model, _FLUX_STR, solution_dictionary=flux_sol)
        except MassSimulationError as e:
            warn(str(e))
            conc_sol = MassSolution(model, _CONC_STR)
            flux_sol = MassSolution(model, _FLUX_STR)

        if update_initial_conditions:
            self.update_values(model)

        return conc_sol, flux_sol

    def find_steady_state(self, models=None, strategy="simulate",
                          perturbations=None, verbose=True,
                          update_initial_conditions=False,
                          update_reactions=False, **options):
        """Find the steady state for MassModel(s) using the given strategy.

        Parameters
        ----------
        models: mass.MassModel, list of mass.MassModels, None
            The model or models to simulate. If None provided, all models in
            the Simulation object will be simulated. Otherwise, any provided
            models must already be present in the Simulation object.
        strategy : 'simulate' or 'roots'
            A string representing the strategy to use to solve for the steady
            state. Can be 'simulate' to simulate the model to steady state, or
            it can be 'roots' to calculate the roots.
        perturbations: dict, optional
            A dict of events to incorporate into the simulation, where keys are
            the object and type of event to incorporate and values are the
            change to be made.
        verbose: bool, optional
            If True, print a detailed report of why the method failed.
            Default is True.
        update_reactions : bool, optional
            If True, update the steady state flux attribute in each reaction.
        update_initial_conditions : bool, optional
            If True, update the model initial conditions for each metabolite,
            and update the Simulation.
        **options
            Options for the solver. Options must correspond to the set solver
            for the given stratgey, otherwise they will be ignored.

        Returns
        -------
        conc_solutions: DictList of mass.Solution objects
            A DictList of mass.Solution objects, each containing the dict of
            concentrations solutions for successful simulations, and an
            empty dict for failed simulations.
        flux_solutions: DictList of mass.Solution objects
            A DictList of mass.Solution objects, each containing the dict of
            flux solutions for successful simulations, and an empty dict for
            failed simulations.

        """
        # Check for models if a list provided, remove models from the list
        # that are not currently in the Simulation
        models = self._get_models_for_method(models, verbose)
        conc_solutions = DictList()
        flux_solutions = DictList()

        # Find the steady state for models that are present.
        for model in models:
            if verbose:
                print("Finding steady state for MassModel: " + str(model))
            sols = self.find_steady_state_model(
                model, strategy=strategy, perturbations=perturbations,
                verbose=verbose,
                update_initial_conditions=update_initial_conditions,
                update_reactions=update_reactions,
                **options)

            sols = list(sols)
            for sol, storage in zip(sols, [conc_solutions, flux_solutions]):
                storage += [sol]

        if update_initial_conditions:
            self.update_values(models=models)

        return conc_solutions, flux_solutions

    def update_values(self, models=None):
        """Update the Simulation with the models current numerical values.

        Values include initial conditions and all parameters.

        Parameters
        ----------
        models: mass.MassModel, list of mass.MassModels, None
            The model or models to update. If None provided, all models in
            the Simulation object will be updated. Otherwise, any provided
            models must already be present in the Simulation object.

        """
        models = self._get_models_for_method(models, False)
        new_values = DictList()
        for model in models:
            keys = [model.id + "_parameters", model.id + "_ics"]
            funcs = [self._get_parameters_from_model, self._get_ics_from_model]
            for key, func in zip(keys, funcs):
                self._values.remove(key)
                new_values.add(func(model))

        self._values += new_values

    def make_pools(self, pools, parameters=None, verbose=True):
        """Create Pool Solutions for a given list of pools.

        Example: For the reaction v1: x1 <=> x2 with Keq = 2,  a conservation
        pool and a disequilibrium pool can be defined by providing the
        following input for pools and parameters:

            pools = ['x1 + x2', 'x1 - x2/Keq_v1']
            parameters = {'Keq_v1' : 2}

        Parameters
        ----------
        pools : str, list of str, dict
            Either a string or a list of strings defining the pooling formula,
            or a dict where keys are pool identifiers and values are the
            corresponding pools. All metabolites to be pooled must exist in
            the Solution objects found in self.get_concentration_solutions().
        parameters : dict, optional
            A dictionary of aditional parameters to be used in the pools,
            where the key:value pairs are the parameter identifiers and its
            numerical value.
        verbose: bool, optional
            If True, provide warnings when pools cannot be created using a
            particular Solution object. Default is True.

        Returns
        -------
        pool_solutions: mass.Solution, DictList of mass.Solution objects
            A DictList of mass.Solution object each containing the dict of
            pool solutions, or a single mass.Solution object if there is only
            one model in the Simulation object.

        Notes
        -----
        If a pool cannot be created for a model due to including metabolites in 
            the pool formula that are not part of the model, the creation of
            that particular pool will be skipped. To enable a warning for when
            pool creation is skipped, set 'verbose' equal to True.

        """
        group_ids = None

        if isinstance(pools, string_types):
            pools = [pools]
        elif isinstance(pools, dict):
            group_ids = list(iterkeys(pools))
            pools = list(itervalues(pools))
        else:
            pools = ensure_iterable(pools)

        sols = self.get_concentration_solutions()
        if isinstance(sols, Solution):
            sols = {sols.id.replace("_ConcSol", ""): sols}
        if group_ids is None:
            group_ids = ["p{0}".format(str(i + 1)) for i in range(len(pools))]

        pool_solutions = self._create_group(
            sols=sols, to_create=pools, parameters=parameters, 
            new_ids=group_ids, sol_type=_POOL_STR, verbose=verbose)

        return pool_solutions

    def make_net_fluxes(self, net_fluxes, parameters=None, verbose=True):
        """Create NetFlux Solutions for a given list of flux summations.

        Example: To sum the fluxes of v1 and v2 scaled,
            net_fluxes = ['v1 + scalar * v2']
            parameters = {'scalar': 10}

        Parameters
        ----------
        net_fluxes : str, list of str, dict
            Either a string or a list of strings defining the net flux formula,
            or a dict where keys are net flux identifiers and values are the
            corresponding net fluxes. All fluxes to be combined must exist in
            the Solution objects found in self.get_flux_solutions().
        parameters : dictionary, optional
            A dictionary of aditional parameters to be used in the net fluxes,
            where the key:value pairs are the parameter identifiers and its
            numerical value.
        verbose: bool, optional
            If True, provide warnings when net fluxes cannot be created using a 
            particular Solution object. Default is True.

        Returns
        -------
        net_flux_solutions: mass.Solution, DictList of mass.Solution objects
            A DictList of mass.Solution object each containing the dict of
            net flux solutions, or a single mass.Solution object if there is
            only one model in the Simulation object.

        Notes
        -----
        If a net flux cannot be created for a model due to including reactions
            in the net flux formula that are not part of the model, the 
            creation of that particular net flux will be skipped. To enable a 
            warning for when net flux creation is skipped, set 'verbose' equal
            to True.

        """
        group_ids = None
        if isinstance(net_fluxes, string_types):
            net_fluxes = [net_fluxes]
        elif isinstance(net_fluxes, dict):
            group_ids = list(iterkeys(net_fluxes))
            net_fluxes = list(itervalues(net_fluxes))
        else:
            net_fluxes = ensure_iterable(net_fluxes)

        sols = self.get_flux_solutions()
        if isinstance(sols, Solution):
            sols = {sols.id.replace("_FluxSol", ""): sols}

        if group_ids is None:
            group_ids = ["net{0}".format(str(i + 1))
                         for i in range(len(net_fluxes))]

        net_fluxes_solutions = self._create_group(
            sols=sols, to_create=net_fluxes, parameters=parameters, 
            new_ids=group_ids, sol_type=_NETFLUX_STR, verbose=verbose)

        return net_fluxes_solutions

    # Internal
    def _set_rate_type(self, model):
        for rxn in model.reactions:
            # Ignore reactions with custom rates
            if rxn in model.custom_rates:
                continue
            # Reversible reactions depend on at least two parameters.
            elif rxn.reversible:
                parameters = rxn.parameters
                # Forward rate and equilibrium constants for type 1
                if rxn.kf_str in parameters and rxn.Keq_str in parameters:
                    rxn._rtype = 1
                # Forward and reverse rate constants for type 2
                elif rxn.kf_str in parameters and rxn.kr_str in parameters:
                    rxn._rtype = 2
                # Equilibrium and reverse rate constants for type 3
                elif rxn.Keq_str in parameters and rxn.kr_str in parameters:
                    rxn._rtype = 3
                # Default to rtype=1
                else:
                    rxn._rtype = 1
            # Irreversible reactions only depend on a kf
            else:
                rxn._rtype = 1

    def _create_group(self, sols, to_create, parameters, new_ids, sol_type, 
                      verbose):
        """Create a group and compute its solution.

        Warnings
        --------
        This method is intended for internal use only.

        """
        group_solutions = DictList()
        for sol in itervalues(sols):
            groups_sol_dict = {}
            items = [key for key in iterkeys(sol)]
            interpolate = sol.interpolate
            if interpolate:
                sol.interpolate = False
            groups_id_dict = dict(zip(new_ids, to_create))
            for g_id, group in iteritems(groups_id_dict):
                if isinstance(group, sym.Basic):
                    group = str(strip_time(group))
                args = sorted([arg for arg in items if arg in group])
                local_syms = {str(arg): sym.Symbol(arg) for arg in args}
                if parameters is not None:
                    local_syms.update({str(param): sym.Symbol(param)
                                       for param in iterkeys(parameters)})
                else:
                    parameters = {}
                try:
                    expr = sym.sympify(group, locals=local_syms)
                    expr = expr.subs(parameters)
                    extras = [exp_sym for exp_sym in expr.atoms(sym.Symbol)
                              if exp_sym not in list(itervalues(local_syms))]
                    if extras:
                        raise MassSimulationError

                except MassSimulationError:
                    if verbose:
                        warn("For Solution '{0}', extra arguments {1} found. "
                             "Therefore skipping creation of '{2}': {3}."
                             .format(sol.id, extras, g_id, group))
                    continue

                except sym.SympifyError:
                    raise ValueError("Could not convert expression '{0}: {1}' "
                                     "into a group.".format(g_id, group))

                func = sym.lambdify(args=[sym.Symbol(arg) for arg in args],
                                    expr=expr, modules=_LAMBDIFY_MODULE)
                values = np.array([sol[arg] for arg in args])
                groups_sol_dict.update({g_id: func(*values)})
            group_sol = MassSolution(
                sol.model, sol_type, groups_sol_dict, sol.t)
            group_sol._groups = groups_id_dict
            if interpolate:
                sol.interpolate = interpolate
                group_sol.interpolate = interpolate

            group_solutions.add(group_sol)

        self._update_solution_storage(group_solutions)

        if len(group_solutions) == 1:
            group_solutions = group_solutions[0]

        return group_solutions

    def _lookup_solutions(self, models, sol_type):
        """Return a dictionary solutions for the given models and 'sol_type'.

        Warnings
        --------
        This method is intended for internal use only.

        """
        # Check for models if a list provided, remove models from the list
        # that are not currently in the Simulation
        models = self._get_models_for_method(models, False)
        # Set up loop
        retry = True
        while retry:
            # Try to access solutions
            try:
                if not models:
                    sols = models
                    break
                sols = self._solutions.get_by_any([
                    "{0}_{1}Sol".format(str(model), sol_type)
                    for model in models])
                retry = False
            # If a solution does not exist for a model, remove the model from
            # the list and try again.
            except KeyError as e:
                model = str(e).replace("_{0}Sol".format(sol_type), "")
                model = self.models.get_by_id(model.strip("'"))
                warn("No '{0}Sol' found for '{1}'".format(sol_type, model.id))
                models.remove(model)

        solutions = {str(model): sol for model, sol in zip(models, sols)}
        if len(solutions) == 1:
            return solutions[str(models[0])]
        else:
            return solutions

    def _check_for_one_model(self, model, verbose=False):
        """Determine whether model is present in the Simulation object.

        Warnings
        --------
        This method is intended for internal use only.

        """
        try:
            model = self.models.get_by_id(str(model))
        except KeyError:
            msg = "Could not find '{0}' in self.models. ".format(str(model))
            if verbose:
                msg += ("Check the method input and ensure that model '{0}' "
                        "has been added to the Simulation through "
                        "self.add_models.".format(str(model)))
            raise MassSimulationError(msg)

    def _check_for_multiple_models(self, models, verbose=False):
        """Determine whether models are present in the Simulation object.

        Warnings
        --------
        This method is intended for internal use only.

        """
        models = ensure_iterable(models)
        for model in models:
            try:
                self._check_for_one_model(model, verbose)
            except MassSimulationError as e:
                warn(e)
                models.remove(model)
        return models

    def _get_parameters_from_model(self, model):
        """Get a dict of parameters as sympy symbols and their values.

        Warnings
        --------
        This method is intended for internal use only.

        """
        values = {sym.Symbol(str(k)): v 
                  for k, v in iteritems(model._get_all_parameters())}

        return DictWithID(model.id + "_parameters", dictionary=values)

    def _get_ics_from_model(self, model):
        """Get a dict of initial conditions as sympy symbols and their values.

        Warnings
        --------
        This method is intended for internal use only.

        """
        values = {_mk_met_func(k): v
                  for k, v in iteritems(model.initial_conditions)}
        return DictWithID(model.id + "_ics", dictionary=values)

    def _assess_quality(self, model, verbose, obj="Simulation"):
        """Assess the model quality to determine whether it is simulatable.

        The "obj" determines whether to check the model or the stored model
        values in the Simulation object.

        Warnings
        --------
        This method is intended for internal use only.

        """
        if obj == "Model":
            [sim_check, consistency_check] = is_simulatable(model)
            if verbose and (not sim_check or not consistency_check):
                # TODO Add thermodynamic consistency once implemented
                qcqa_model(model, parameters=not sim_check,
                           concentrations=not sim_check,
                           superfluous=not sim_check,
                           thermodynamic=False)

        elif obj == "Simulation":
            [sim_check, consistency_check] = is_simulatable(model, self)
            if verbose and (not sim_check or not consistency_check):
                # TODO Add thermodynamic consistency once implemented
                qcqa_simulation(self, model, parameters=not sim_check,
                                concentrations=not sim_check,
                                superfluous=not consistency_check,
                                thermodynamic=False)
        else:
            raise ValueError("obj must be one of the following {0!r}"
                             .format({"Simulation", "Model"}))

        if not consistency_check:
            warn("MassModel '{0}' has numerical consistency issues. Simulation"
                 " results may not be accurate.".format(model.id))

        if not sim_check:
            raise MassSimulationError("Unable to simulate MassModel '{0}'"
                                      .format(model.id))

    def _apply_perturbations(self, model, perturbations):
        """Apply the given perturbations to the value substituion dictionaries.

        Warnings
        --------
        This method is intended for internal use only.

        """
        # Get value substituion dictionaries
        functions = {}
        parameters, ics = self.view_model_values(model)
        if not perturbations:
            return parameters, ics, functions

        # Validate the perturbations to ensure they can be interpreted
        perturbations = self._validate_perturbations(perturbations)
        # Iterate through perturbations
        for item_str, value in iteritems(perturbations):
            # Perturb ICs and fixed concentrations
            if ".ic" in item_str[-3:] or ".fix" in item_str[-4:]:
                # Define correct value dict
                if sym.Symbol(item_str.split(".")[0]) in parameters:
                    value_dict = parameters
                else:
                    value_dict = ics
                is_ic = _IC_RE.search(item_str)
                item_sym = self._perturbation_string_replace(item_str, value,
                                                             "\.ic|\.fix",
                                                             value_dict, True)
                # Switch dictionaries for ICs changed into fixed concentrations
                if is_ic:
                    ics[item_sym] = value_dict[item_sym]
                else:
                    parameters[item_sym] = value_dict[item_sym]
                if item_sym in ics and not is_ic:
                    del ics[item_sym]

            # Perturb a fixed concentration to be a function
            elif ".func" in item_str[-5:]:
                item_sym = sym.Symbol(item_str[:-5])
                ics[item_sym] = parameters[item_sym]
                item_sym = self._perturbation_string_replace(
                    item_str, value, "\.func$", parameters, False)

                functions.update({item_sym: parameters[item_sym]})
                parameters[item_sym] = _mk_met_func(item_str[:-5])
            # Perturb a custom parameter
            elif ".custom" in item_str[-7:]:
                self._perturbation_string_replace(
                    item_str, value, "\.custom$", parameters, True)
            # Otherwise perturbation is on a rate parameter
            else:
                self._perturbation_string_replace(
                    item_str, value, "", parameters, True)

        return parameters, ics, functions

    def _validate_perturbations(self, perturbations):
        """Validate the given perturbations and raise errors if not valid.

        Warnings
        --------
        This method is intended for internal use only.

        """
        if not isinstance(perturbations, dict):
            raise TypeError("perturbations must be a dict.")
        validated = {}
        _re_list = [_KF_RE, _KEQ_RE, _KR_RE, _IC_RE,
                    _FIX_RE, _FUNC_RE, _CUSTOM_RE]
        _key_fixes = ["kf", "Keq", "kr", ".ic", ".fix", ".func", ".custom"]

        for old_key, value in iteritems(perturbations):
            item, pert_type = old_key.split(".")
            try:
                for _re, key_fix in zip(_re_list, _key_fixes):
                    if _re.match(pert_type):
                        new_key = ("{0}_{1}".format(key_fix, item)
                                   if key_fix in _key_fixes[:3]
                                   else item + key_fix)
                        if old_key in str(value) \
                           and "[{0}]".format(old_key) in str(value):
                            value = value.replace(old_key, new_key)
                        else:
                            # If value cannot be converted, ensure perturbation
                            # is allowed and can be interpreted later.
                            allow = [_FIX_RE.search(new_key)
                                     and _IC_RE.search(str(value)),
                                     _CUSTOM_RE.search(new_key),
                                     _FUNC_RE.search(new_key)]
                            if [True for b in allow if b is not None]:
                                pass
                            else:
                                # Otherwise try to convert value to a float
                                value = float(value)
                        validated[new_key.strip()] = value
                    elif [True for other in _re_list if other != _re
                          and other.match(pert_type)]:
                        pass
                    else:
                        raise ValueError
            except ValueError:
                raise MassSimulationError("Perturbation '{0!r} : {1!r}' not "
                                          "recognized.".format(old_key, value))
        return validated

    def _perturbation_string_replace(self, item_str, value, elim_pattern,
                                     value_dict, as_num):
        """Modify the perturbation strings to help format perturbations.

        Warnings
        --------
        This method is intended for internal use only.

        """
        item_str = re.sub(elim_pattern, "", item_str)
        if _mk_met_func(item_str) in value_dict:
            item_sym = _mk_met_func(item_str)
        else:
            item_sym = sym.Symbol(item_str)
        if item_sym not in value_dict:
            raise MassSimulationError("Could not find the '{0}' in defined "
                                      "Simulation values".format(item_str))
        elif item_str in str(value):
            if re.search(elim_pattern, str(value)):
                value = re.sub(elim_pattern, "", str(value))
            value = str(value).replace("[{0}]".format(item_str),
                                       str(value_dict[item_sym]))
        if as_num:
            value_dict[item_sym] = float(sym.sympify(value))
        else:
            value_dict[item_sym] = sym.sympify(value)
        return item_sym

    def _make_odes_functions(self, model, parameters, ics, functions):
        """Create lambda functions of the ODEs and the jacobian matrix.

        Will also ensure the order of the initial conditions matches the
        order of the ODEs.

        Warnings
        --------
        This method is intended for internal use only.

        """
        ordered_ics = OrderedDict()
        equations = OrderedDict()

        # Set up matrix of ODEs
        for met, equation in iteritems(model.odes):
            met_func = _mk_met_func(met)
            if met_func in ics:
                equations[met_func.diff(_T_SYM)] = equation.subs(parameters)
                ordered_ics[met_func] = ics[met_func]

        # Account for functions
        if functions:
            for met, func in iteritems(functions):
                met_func = _mk_met_func(met)
                equations[met_func.diff(_T_SYM)] = func.diff(_T_SYM)
                ordered_ics[met_func] = ics[met]

        metabolite_matrix = sym.Matrix(list(iterkeys(ordered_ics)))
        equations = metabolite_matrix.diff(_T_SYM).subs(equations)
        # Turn ODEs into lambda function
        lambda_odes = sym.lambdify(args=(_T_SYM, metabolite_matrix),
                                   expr=equations, modules=_LAMBDIFY_MODULE)

        # Determine Jacobian and turn into lambda function if possible
        try:
            jacb = equations.jacobian(metabolite_matrix)
            lambda_jacb = sym.lambdify(args=(_T_SYM, metabolite_matrix),
                                       expr=jacb, modules=_LAMBDIFY_MODULE)
        except TypeError:
            lambda_jacb = None

        return [lambda_odes, lambda_jacb, ordered_ics]

    def _integrate_odes(self, time, odes, jacb, ics, new_options):
        """Integrate the ODEs using the set solver and its options.

        Warnings
        --------
        This method is intended for internal use only.

        """
        # TODO Add additional solvers here, fix for universal input and output
        # once additional solvers implemented
        integrator = {
            "scipy": self._integrate_with_scipy}

        t, concs, sol_obj = integrator[self.solver](time, odes, jacb, ics,
                                                    new_options)
        # Map identifiers to their solutions, and store in a Solution object.
        id_list = strip_time(list(iterkeys(ics)))
        concs = {str(_id): sol for _id, sol in zip(id_list, concs)}

        return t, concs, sol_obj

    def _calculate_flux_solutions(self, model, parameters, conc_sol, t):
        """Calculate fluxes using rate equations and concentrations solutions.

        Warnings
        --------
        This method is intended for internal use only.

        """
        fluxes = {}
        for reaction, rate in iteritems(model.rates):
            rate = strip_time(rate.subs(parameters))
            args = tuple(rate.atoms(sym.Symbol))
            # Fluxes dependent on concentration solutions need to be calculated
            if args:
                concs = np.array([conc_sol[str(arg)] for arg in args])
                args = tuple([_T_SYM]) + args
                rate_function = sym.lambdify(args=args, expr=rate,
                                             modules=_LAMBDIFY_MODULE)
                flux = np.array([rate_function(t[i], *concs[:, i])
                                 for i in range(len(t))])
            # Constant flux, therefore flux is identical at each t
            else:
                flux = np.array([float(rate)] * len(t))
            fluxes[reaction.id] = flux

        flux_sol = MassSolution(model, _FLUX_STR, fluxes, t)

        return flux_sol

    def _update_solution_storage(self, solutions):
        """Update the solution storage DictList with the new solutions.

        Warnings
        --------
        This method is intended for internal use only.

        """
        for sol in solutions:
            if sol.id in self._solutions:
                self._solutions.remove(sol.id)
            if sol:
                self._solutions.add(sol)

    def _integrate_with_scipy(self, time, odes, jacb, ics, new_options):
        """Integrate the ODEs using the scipy.

        Warnings
        --------
        This method is intended for internal use only.

        """
        options = self.get_solver_options("scipy")
        for key, value in iteritems(new_options):
            if key in options:
                options[key] = value
            else:
                warn("Unknown kwarg '{0}' provided.".format(key))

        # Set jacobian
        options["jac"] = jacb
        # Remove options not relevant for solver method
        options_to_remove = {
            "RK45": ["lband", "uband", "min_step", 
                     "first_step", "jac_sparsity", "jac"],
            "RK23": ["lband", "uband", "min_step", 
                     "first_step", "jac_sparsity", "jac"],
            "Radau": ["lband", "uband", "min_step", "first_step"],
            "BDF": ["lband", "uband", "min_step", "first_step"],
            "LSODA": ["jac_sparsity"]}

        for option in options_to_remove[options["method"]]:
            del options[option]

        options["vectorized"] = True

        ic_vals = list(itervalues(ics))
        sol_obj = solve_ivp(odes, t_span=time, y0=ic_vals, **options)

        time = sol_obj.t
        unmapped_concs = sol_obj.y
        return time, unmapped_concs, sol_obj

    def _find_steady_state_simulate(self, model, perts, verbose, chop, 
                                    update, **options):
        """Find the steady state solution of a model using 'simulate' strategy.

        Warnings
        --------
        This method is intended for internal use only.

        """
        # Try simulating using a final time of 10^3 up to 10^6.
        power = 3
        fail = 7
        while power <= fail:
            retry = False
            solutions = self.simulate_model(model, time=(0, 10**power),
                                            perturbations=perts,
                                            interpolate=False, verbose=verbose,
                                            return_obj=False,
                                            update_solutions=False, 
                                            **options)
            # Check to see if concentrations reached a steady state.
            for sol in itervalues(solutions[0]):
                if not round(abs(sol[-1] - sol[-2]), chop) <= _ZERO_TOL:
                    break
            if not retry:
                break
            power += 1

        if power > fail:
            if verbose:
                warn("Unable to find a steady state for '{0}' using strategy "
                     "'simulate'.".format(model.id))
            return {}, {}
        else:
            return self._chop_store_sols(model, solutions, chop, update)

    def _find_steady_state_roots(self, model, perts, verbose, chop, 
                                 update, **options):
        """Find the steady state solution of a model using 'roots' strategy.

        Warnings
        --------
        This method is intended for internal use only.

        """
        # Implement perturbations and get the value substitution dicts.
        if perts:
            parameters, ics = self._apply_perturbations(model, perts)
        # Just get value substituion dicts if no perturbations provided.
        else:
            parameters, ics = self.view_model_values(model)

        ordered_ics = OrderedDict()
        equations = OrderedDict()
        for met, equation in iteritems(model.odes):
            met_func = _mk_met_func(met)
            if met_func in ics:
                equations[met_func] = equation.subs(parameters)
                ordered_ics[met_func] = ics[met_func]

        args = strip_time(list(iterkeys(equations)))
        lambda_eqs = sym.lambdify(args, strip_time(itervalues(equations)),
                                  modules=_LAMBDIFY_MODULE)

        # Make lambda function callable
        def root_func(ics):
            return lambda_eqs(*ics)

        input_dict = {"args": (), "method": "lm", 
                      "jac": None, "tol": None, "callback": None}
        for key in iterkeys(input_dict):
            if key in options:
                input_dict[key] = options.pop(key)

        ics = list(itervalues(ordered_ics))
        sol = root(fun=root_func, x0=ics, args=input_dict["args"],
                   method=input_dict["method"], jac=input_dict["jac"], 
                   tol=input_dict["callback"], callback=input_dict["callback"],
                   options=options)
        # Warn if no steady state was reached and return empty sols.
        if not sol.success:
            if verbose:
                warn("Unable to find a steady state for '{0}' using strategy "
                     "'roots'.".format(model.id))
            return {}, {}

        # Create solutions and put in Solution objects to return
        conc_sol = {met: sol for met, sol in zip(args, sol.x)}
        flux_sol = {rxn.id: [rate.subs(parameters).subs(conc_sol)]
                    for rxn, rate in iteritems(strip_time(model.rates))}
        conc_sol = {str(met): [sol] for met, sol in iteritems(conc_sol)}

        return self._chop_store_sols(model, [conc_sol, flux_sol], chop, update)

    def _chop_store_sols(self, model, solutions, chop, update):
        """Chop steady state solutions and store in dictionaries.

        Warnings
        --------
        This method is intended for internal use only.

        """
        
        conc_sol = {}
        flux_sol = {}
        update_initial_conditions, update_reactions = update
        for met, sol in iteritems(solutions[0]):
            conc_sol[met] = round(sol[-1], chop)
            if update_initial_conditions:
                met = model.metabolites.get_by_id(met)
                met.initial_condition = round(sol[-1], chop)

        for rxn, sol in iteritems(solutions[1]):
            flux_sol[rxn] = round(sol[-1], chop)
            if update_reactions:
                rxn = model.reactions.get_by_id(rxn)
                rxn.steady_state_flux = round(sol[-1], chop)

        return conc_sol, flux_sol

    def _get_models_for_method(self, models, verbose):
        """Return a copy of the DictList containing valid models.

        Warnings
        --------
        This method is intended for internal use only.

        """
        if models is None:
            models = self.models
        else:
            models = self._check_for_multiple_models(models, verbose=verbose)
        return models.copy()
