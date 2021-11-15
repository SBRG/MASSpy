# -*- coding: utf-8 -*-
"""Module handling :class:`optlang.interface.Model` for concentration problems.

The purpose of the :class:`ConcSolver` is to provide an interface to assist
with setting up various problems involving optimization-like problems
involving metabolite concentrations. Note that all internal solver variables
exist in logarithmic space and therefore all associated numerical values are
transformed from linear space into log space before being added to the solver.
Unless specified otherwise, all numerical solutions will be transformed back
into linear space from logspace before being returned.

Upon initialization, a generic problem is created, represented by an
:class:`optlang.interface.Model`. The generic problem includes variables for
the metabolite concentration and reaction equilibrium constants, along with
thermodynamic constraints for each reaction with the direction of the
constraint (i.e. greater/less than) dependent on the sign of the steady
state flux.

In addition to creating a generic problem, the :class:`ConcSolver` has the
following functions available to create a specific type of problem:

    * :meth:`ConcSolver.setup_feasible_qp_problem`
    * :meth:`ConcSolver.setup_sampling_problem`

If custom objectives, constraints and/or variables are to be used with the
solver, it is recommended to run one of the above functions first to setup the
problem, and then tailor the solver with customizations for that problem.

Notes
-----
* An :class:`optlang.interface.Model` represents an optimization problem and
  contains the variables, constraints, and objectives that make up the problem.
  See the :mod:`optlang` documentation for more information.
* All numerical values to be utilized by the solver must be defined. This
  includes:

        - Metabolite concentrations for all included metabolites, accessed via
          :attr:`.MassMetabolite.initial_condition`.
        - Reaction equilibrium constants for all included reactions, accessed
          via :attr:`.MassReaction.equilibrium_constant`.
        - Reaction steady state fluxes for all included reactions, accessed
          :attr:`.MassReaction.steady_state_flux`.

"""
import types
from collections import namedtuple
from functools import partial
from warnings import warn

import numpy as np
import optlang
import pandas as pd
from cobra.core.dictlist import DictList
from cobra.exceptions import (
    OPTLANG_TO_EXCEPTIONS_DICT,
    OptimizationError,
    SolverNotFound,
)
from cobra.util.context import get_context, resettable
from cobra.util.solver import get_solver_name, interface_to_str, qp_solvers, solvers
from optlang.symbolics import Zero
from scipy.sparse import dok_matrix, lil_matrix
from six import iteritems, string_types
from sympy import Basic, Matrix, eye

from mass.core.mass_configuration import MassConfiguration
from mass.thermo.conc_solution import (
    get_concentration_solution,
    update_model_with_concentration_solution,
)
from mass.util.qcqa import (
    get_missing_initial_conditions,
    get_missing_reaction_parameters,
    get_missing_steady_state_fluxes,
)
from mass.util.util import (
    _check_kwargs,
    _make_logger,
    apply_decimal_precision,
    ensure_iterable,
    ensure_non_negative_value,
    get_public_attributes_and_methods,
)


LOGGER = _make_logger(__name__)
"""logging.Logger: Logger for :mod:`~.conc_hr_sampler` submodule."""

MASSCONFIGURATION = MassConfiguration()


class ConcSolver:
    """Class providing an interface for concentration mathematical problems.

    Upon initialization, a generic problem is created, represented by an
    :class:`optlang.Model`. The generic problem includes variables for the
    metabolite concentration and reaction equilibrium constants, along with
    thermodynamic constraints for each reaction with the direction of the
    constraint (i.e. greater/less than) dependent on the sign of the steady
    state flux.

    Notes
    -----
    * All internal solver variables exist in logarithmic space and therefore
      all associated numerical values are transformed from linear space into
      log space before being added to the solver.
    * Unless specified otherwise, all numerical solutions will be transformed
      back into linear space from logspace before being returned.
    * Boundary reactions (a.k.a. reactions with only one metabolite involved)
      are excluded automatically. To work with reactions that cross compartment
      boundaries, :class:`.MassMetabolite` objects need to be defined for
      the metabolites in both compartments.

    Parameters
    ----------
    model : MassModel
        The :mod:`mass` model to associated with the :class:`ConcSolver`
        instance. The model is used to populate the solver with typical
        variables and constraints upon initialization.
    excluded_metabolites : iterable or None
        An iterable of model :class:`.MassMetabolite` objects or their
        identifiers in populating the solver with metabolite concentration
        variables and reaction constraints. If ``None``, no metabolites are
        excluded.
    excluded_reactions : iterable or None
        An iterable of model :class:`.MassReaction` objects or their
        identifiers in populating the solver with reaction equilibrium constant
        variables and reaction constraints. If ``None``, no reactions are
        excluded.
    equilibrium_reactions : iterable or None
        An iterable of model :class:`.MassReaction` objects or their
        identifiers that are intended to be at equilibrium. Reactions with
        steady state flux values equal to 0. are typically ignored unless
        they are specified in the :class:`ConcSolver.equilibrium_reactions`
    constraint_buffer : float or None
        A ``float`` value to use as a constraint buffer for all constraints.

        Default is `0.`.
    **kwargs
        exclude_infinite_Keqs :
            ``bool`` indicating whether to exclude reactions with equilibrium
            constant values of infinity from the concentration solver.

            Default is ``True``.
        fixed_conc_bounds :
            An iterable containing metabolites whose concentrations are
            to be set as fixed variables, meaning that their lower and
            upper bounds are equal to the base value.
        fixed_Keq_bounds :
            An iterable containing reactions whose equilibrium constants
            are to be set as fixed variables, meaning that their lower and
            upper bounds are equal to the base value.
        decimal_precision :
            ``bool`` indicating whether to apply the
            :attr:`~.MassBaseConfiguration.decimal_precision` attribute of
            the :class:`.MassConfiguration` to the bound values.

            Default is ``False``.
        zero_value_log_substitute :
            ``float`` value to substitute for 0 when trying to take the
            logarithm of 0 to avoid a domain error.

            Default is ``1e-10``.

    Attributes
    ----------
    problem_type : str
        The type of mathematical problem that the concentration solver has
        been setup to solve.
    excluded_metabolites : list
        A ``list`` of metabolite identifiers for model metabolites to exclude
        in populating the solver with metabolite concentration variables and
        reaction constraints.
    excluded_reactions : list
        A ``list`` of reaction identifiers for model reactions to exclude
        in populating the solver with reaction equilibrium constant variables
        and reaction constraints.
    equilibrium_reactions : list
        A ``list`` of reaction identifiers for model reactions that are
        intended to be at equilibrium.
    constraint_buffer : float
        A value to utilize when setting a constraint buffer.

    """

    def __init__(
        self,
        model,
        excluded_metabolites=None,
        excluded_reactions=None,
        equilibrium_reactions=None,
        constraint_buffer=0,
        **kwargs
    ):
        """Initialize the ConcSolver."""
        kwargs = _check_kwargs(
            {
                "exclude_infinite_Keqs": True,
                "fixed_conc_bounds": [],
                "fixed_Keq_bounds": [],
                "decimal_precision": False,
                "zero_value_log_substitute": 1e-10,
            },
            kwargs,
        )
        self._model = model
        model._conc_solver = self

        interface = MASSCONFIGURATION.solver
        self._solver = interface.Model()
        self.solver.objective = interface.Objective(Zero)
        self.problem_type = ""

        self._tolerance = None
        self.tolerance = MASSCONFIGURATION.tolerance

        self.excluded_metabolites = []
        self.excluded_reactions = []
        self.equilibrium_reactions = []

        self.constraint_buffer = ensure_non_negative_value(constraint_buffer)

        self._zero_value_log_substitute = ensure_non_negative_value(
            kwargs.pop("zero_value_log_substitute"), exclude_zero=True
        )

        # Try setting excluded and equilibrium attributes specified in kwargs
        if excluded_reactions is not None:
            excluded_reactions = [getattr(r, "_id", str(r)) for r in excluded_reactions]
        else:
            excluded_reactions = []
        excluded_reactions += [
            r.id for r in model.boundary if r.id not in excluded_reactions
        ]
        if kwargs.pop("exclude_infinite_Keqs"):
            excluded_reactions += [
                r.id
                for r in model.reactions
                if r.Keq == float("inf") and r.id not in excluded_reactions
            ]

        exclude_and_equilibrium = {
            "excluded_metabolites": excluded_metabolites,
            "excluded_reactions": excluded_reactions,
            "equilibrium_reactions": equilibrium_reactions,
        }
        for key, value in iteritems(exclude_and_equilibrium):
            try:
                if value is not None:
                    getattr(self, "add_" + key)(value)
            except ValueError as e:
                warn(
                    "Could not set `{0}` due to the following: {1}".format(key, str(e))
                )

        # Ensure model has Keq parameters and initial conditions defined
        missing = self._check_for_missing_values(model)
        for key, missing_values in iteritems(missing):
            if missing_values:
                LOGGER.warning("Missing %s for the following: %r", key, missing_values)

        self._initialize_solver(**kwargs)

    @property
    def model(self):
        """Return the model associated with the :class:`ConcSolver`."""
        return getattr(self, "_model")

    @property
    def solver(self):
        """Get or set the attached solver of the :class:`ConcSolver`.

        When using a :class:`~cobra.util.context.HistoryManager` context,
        this attribute can be set temporarily, reversed when the exiting
        the context.

        Parameters
        ---------
        value : str
            The optimization solver for problems concerning metabolite
            concentrations. The solver choices are the ones provided by
            :mod:`optlang` and solvers installed in your environment. Valid
            solvers typically include: ``"glpk"``, ``"cplex"``, ``"gurobi"``

        Notes
        -----
        * Like the :attr:`Model.solver <cobra.core.model.Model.solver>`
          attribute, the concentration solver instance is the associated
          solver object, which manages the interaction with the associated
          solver, e.g. glpk.

        * This property is useful for accessing the concentration optimization
          problem directly and for defining additional constraints.

        """
        return getattr(self, "_solver")

    @solver.setter
    @resettable
    def solver(self, value):
        """Set the attached solver of the :class:`ConcSolver`."""
        not_valid_interface = SolverNotFound(
            "{0} is not a valid solver interface. Pick from {1}.".format(
                value, str(list(solvers))
            )
        )

        if isinstance(value, string_types):
            try:
                interface = solvers[interface_to_str(value)]
            except KeyError:
                raise not_valid_interface
        elif isinstance(value, types.ModuleType) and hasattr(value, "Model"):
            interface = value
        elif isinstance(value, optlang.interface.Model):
            interface = value.interface
        else:
            raise not_valid_interface

        if self.problem != interface:
            setattr(self, "_solver", interface.Model.clone(self.solver))

    @property
    def tolerance(self):
        """Get or set the tolerance for the solver of the :class:`ConcSolver`.

        When using a :class:`~cobra.util.context.HistoryManager` context,
        this attribute can be set temporarily, reversed when the exiting
        the context.

        Parameters
        ----------
        value : float
            The tolerance of the concentration solver.

        """
        return getattr(self, "_tolerance")

    @tolerance.setter
    @resettable
    def tolerance(self, value):
        """Set the tolerance solver of the :class:`ConcSolver`."""
        solver_tolerances = self.solver.configuration.tolerances

        try:
            solver_tolerances.feasibility = value
        except AttributeError:
            LOGGER.info(
                "The current solver doesn't allow setting" "feasibility tolerance."
            )

        try:
            solver_tolerances.optimality = value
        except AttributeError:
            LOGGER.info(
                "The current solver doesn't allow setting" "optimality tolerance."
            )

        try:
            solver_tolerances.integrality = value
        except AttributeError:
            LOGGER.info(
                "The current solver doesn't allow setting" "integrality tolerance."
            )

        setattr(self, "_tolerance", value)

    @property
    def objective(self):
        """Get or set the solver objective.

        When using a :class:`~cobra.util.context.HistoryManager` context,
        this attribute can be set temporarily, reversed when the exiting
        the context.

        Parameters
        ----------
        value : dict, str, int, MassMetabolite, ~optlang.interface.Objective, or ~sympy.core.basic.Basic
            The following are allowable values to set as the objective.

                - ``dict`` where metabolites are keys, linear coefficients as
                  values.
                - ``str`` identifier of a :class:`.MassMetabolite` or the
                  metabolite object itself.
                - ``int`` metabolite index in :attr:`.MassModel.metabolites`
                - An :class:`~optlang.interface.Objective` or a :mod:`sympy`
                  expression to be directly interpreted as objectives.

        """  # noqa: E501
        return self.solver.objective

    @objective.setter
    def objective(self, value):
        """Set the solver objective."""
        if isinstance(value, Basic):
            value = self.problem.Objective(value, sloppy=False)
        if not isinstance(value, (dict, optlang.interface.Objective)):
            try:
                metabolites = self.model.metabolites.get_by_any(value)
            except KeyError:
                raise ValueError("Invalid objective.")
            value = {m: 1 for m in metabolites}
        self.set_objective(value, additive=False)

    @property
    def objective_direction(self):
        """Get or set the objective direction.

        When using a :class:`~cobra.util.context.HistoryManager` context,
        this attribute can be set temporarily, reversed when the exiting
        the context.

        Parameters
        ----------
        value : str
            The objective direction. Can be either ``"max"`` for the maximum,
            or ``"min"`` for the minimum.

        """
        return self.solver.objective.direction

    @objective_direction.setter
    @resettable
    def objective_direction(self, value):
        """Set the objective direction."""
        value = value.lower()
        if value.startswith("max"):
            self.solver.objective.direction = "max"
        elif value.startswith("min"):
            self.solver.objective.direction = "min"
        else:
            raise ValueError("Unknown objective direction '{}'.".format(value))

    @property
    def problem(self):
        """Return the interface to the underlying mathematical problem.

        Solutions to the :class:`ConcSolver` are obtained by formulating a
        mathematical problem and solving it. The :mod:`optlang` package is
        used to accomplish that and with this property, the problem interface
        can be accessed directly.

        Returns
        -------
        optlang.interface
            The problem interface that defines methods for interacting with
            the problem and associated solver directly.

        """
        return self.solver.interface

    @property
    def variables(self):
        """Return the mathematical variables in the :class:`ConcSolver`.

        In a :class:`ConcSolver`, most variables are metabolites and reaction
        equilibrium constants. However, for specific use cases, it may also be
        useful to have other types of variables. This property defines all
        variables currently associated with the underlying problem of the
        :class:`ConcSolver`.

        Notes
        -----
        * All variables exist in logspace.

        Returns
        -------
        optlang.container.Container
            A container with all associated variables.

        """
        return self.solver.variables

    @property
    def constraints(self):
        """Return the constraints in the :class:`ConcSolver`.

        In a :class:`ConcSolver`, most constraints are thermodynamic
        constraints relating the reaction equilibrium constant to the reaction
        metabolite concentrations.However, for specific use cases, it may also
        be useful to have other types of constraints. This property defines
        all constraints currently associated with the underlying problem of
        the :class:`ConcSolver`.

        Notes
        -----
        * All constraints exist in logspace.

        Returns
        -------
        optlang.container.Container
            A container with all associated constraints.

        """
        return self.solver.constraints

    @property
    def included_metabolites(self):
        """Return a ``list`` of metabolite identifiers included in the solver.

        These are the metabolites not in the
        :attr:`ConcSolver.excluded_metabolites` attribute.

        """
        return [
            m.id
            for m in self.model.metabolites
            if m.id not in self.excluded_metabolites
        ]

    @property
    def included_reactions(self):
        """Return a ``list`` of reaction identifiers included in the solver.

        These are the reactions not in the
        :attr:`ConcSolver.excluded_reactions` attribute.

        """
        return [
            r.id for r in self.model.reactions if r.id not in self.excluded_reactions
        ]

    @property
    def zero_value_log_substitute(self):
        """Get or set the a value to substitute for 0 when taking the log of 0.

        A value of ``1e-10`` means that instead of attempting ``log(0)`` which
        causes a :class:`ValueError`, it will be instead calculated as
        ``log(1e-10)``.

        Parameters
        ----------
        value : float
            A positive value to use instead of 0 when taking the logarithm.

        """
        return getattr(self, "_zero_value_log_substitute")

    @zero_value_log_substitute.setter
    def zero_value_log_substitute(self, value):
        """Set the value to substitute for 0. when taking the log of 0."""
        value = ensure_non_negative_value(value, exclude_zero=True)
        setattr(self, "_zero_value_log_substitute", value)

    def setup_sampling_problem(
        self,
        metabolites=None,
        reactions=None,
        conc_percent_deviation=0.2,
        Keq_percent_deviation=0.2,
        **kwargs
    ):
        """Set up the solver's mathematical problem for concentraiton sampling.

        Notes
        -----
        * This involves changing solver variable bounds based on the percent
          deviation of the base value, removing the objective, and setting the
          :attr:`~ConcSolver.problem_type` to ``"sampling"``.
        * If a percent deviation value is large enough to create a negative
          lower bound, it is set as the
          :attr:`~.ConcSolver.zero_value_log_substitute` value
          to ensure that returned values are not negative.

        Parameters
        ----------
        metabolites : iterable or None
            An iterable of metabolites whose concentration variable  bounds
            are to be changed. If ``None``, all metabolites except those that
            are in the :attr:`ConcSolver.excluded_metabolites` list are used.
        reactions : iterable or None
            An iterable of reactions whose equilibrium constant variable bounds
            are to be changed. If ``None``, all reactions except those that
            are in the :attr:`ConcSolver.excluded_reactions` list are used.
        conc_percent_deviation : float
            A non-negative number indicating the percent to deviate from the
            initial concentration to set as the lower and upper bounds for
            sampling.

            If a value of 0. is given, all given reaction equilibrium
            constants are set as fixed variables. Default is ``0.2`` for
            a 20% deviation from the base value.
        Keq_percent_deviation : float
            A non-negative number indicating the percent to deviate from the
            base reaction equilibrium constant to set as the lower and
            upper bounds for sampling.

            If a value of 0. is given, all given reaction equilibrium
            constants are set as fixed variables. Default is ``0.2`` for
            a 20% deviation from the base value.
        **kwargs
            fixed_conc_bounds :
                An iterable containing metabolites whose concentrations are
                to be set as fixed variables, meaning that their lower and
                upper bounds are equal to the base value.
            fixed_Keq_bounds :
                An iterable containing reactions whose equilibrium constants
                are to be set as fixed variables, meaning that their lower and
                upper bounds are equal to the base value.
            decimal_precision :
                ``bool`` indicating whether to apply the
                :attr:`~.MassBaseConfiguration.decimal_precision` attribute of
                the :class:`.MassConfiguration` to the bound values.

                Default is ``False``.

        """
        kwargs = _check_kwargs(
            {
                "fixed_conc_bounds": [],
                "fixed_Keq_bounds": [],
                "decimal_precision": False,
            },
            kwargs,
        )

        # Get fixed variables
        fixed_conc_bounds = kwargs.pop("fixed_conc_bounds")
        fixed_Keq_bounds = kwargs.pop("fixed_Keq_bounds")

        # Get list of included metabolites, reactions, and Keq_strs
        metabolites = self._get_included_metabolites(metabolites)
        reactions = self._get_included_reactions(reactions)
        Keq_strs = reactions.list_attr("Keq_str")

        for var in self.variables:
            if var.name in metabolites:
                # Set bounds if metabolite
                met = metabolites.get_by_id(var.name)
                if met.initial_condition is None:
                    bounds = (0, np.inf)
                    LOGGER.info(
                        "No initial conditions defined for '%s', no "
                        "variable bounds set",
                        met.id,
                    )
                elif any([m in fixed_conc_bounds for m in [met, met.id]]):
                    bounds = (met.initial_condition, met.initial_condition)
                else:
                    bounds = (conc_percent_deviation, conc_percent_deviation)
                    bounds = _get_deviation_values(met.initial_condition, *bounds)
                upper_def = np.inf
            elif var.name in Keq_strs:
                # Set bounds if reaction Keq variable
                rxn = reactions[Keq_strs.index(var.name)]
                if rxn.Keq is None:
                    bounds = (0, np.inf)
                    LOGGER.info(
                        "No Keq defined for '%s', no variable bounds " "set", rxn.id
                    )
                elif any([r in fixed_Keq_bounds for r in [rxn, rxn.id, rxn.Keq_str]]):
                    bounds = (rxn.Keq, rxn.Keq)
                else:
                    bounds = (Keq_percent_deviation, Keq_percent_deviation)
                    bounds = _get_deviation_values(rxn.Keq, *bounds)
                upper_def = MASSCONFIGURATION.irreversible_Keq
            else:
                continue
            # Get log values of bounds and set
            bounds = _get_log_bounds(
                bounds,
                upper_def,
                kwargs.get("decimal_precision"),
                self.zero_value_log_substitute,
            )
            var.lb, var.ub = bounds

        # Ensure objective is zero for sampling
        self.objective = Zero
        self.problem_type = "sampling"

    def setup_feasible_qp_problem(self, metabolites=None, reactions=None, **kwargs):
        """Set up the solver's mathematical problem for feasible conc. QP.

        Notes
        -----
        * This involves changing solver variable bounds to ``[0, inf]``,
          setting the objective as a QP problem, and setting the
          :attr:`~ConcSolver.problem_type` to ``"feasible_qp"``.
        * The :attr:`ConcSolver.solver` must have QP capabilities.

        Parameters
        ----------
        metabolites : iterable or None
            An iterable of metabolites whose concentration variable bounds
            are to be changed. If ``None``, all metabolites except those that
            are in the :attr:`ConcSolver.excluded_metabolites` list are used.
        reactions : iterable or None
            An iterable of reactions whose equilibrium constant variable bounds
            are to be changed. If ``None``, all reactions except those that
            are in the :attr:`ConcSolver.excluded_reactions` list are used.
        **kwargs
            fixed_conc_bounds :
                An iterable containing metabolites whose concentrations are
                to be set as fixed variables, meaning that their lower and
                upper bounds are equal to the base value.
            fixed_Keq_bounds :
                An iterable containing reactions whose equilibrium constants
                are to be set as fixed variables, meaning that their lower and
                upper bounds are equal to the base value.
            decimal_precision :
                ``bool`` indicating whether to apply the
                :attr:`~.MassBaseConfiguration.decimal_precision` attribute of
                the :class:`.MassConfiguration` to the bound values.

                Default is ``False``.

        Raises
        ------
        TypeError
            Raised when the current solver does not have QP capabilities

        See Also
        --------
        choose_solver
            Method to choose a solver with QP capabilities

        """
        # Ensure solver is QP capable
        if interface_to_str(self.solver.interface) not in qp_solvers:
            raise TypeError(
                "Solver does not have QP capabilities. Utilize the "
                "`ConcSolver.choose_solver` method to set a QP "
                "capable solver."
            )

        kwargs = _check_kwargs(
            {
                "fixed_conc_bounds": [],
                "fixed_Keq_bounds": [],
                "decimal_precision": False,
            },
            kwargs,
        )

        # Get fixed variables
        fixed_conc_bounds = kwargs.pop("fixed_conc_bounds")
        fixed_Keq_bounds = kwargs.pop("fixed_Keq_bounds")

        # Get list of included metabolites, reactions, and Keq_strs
        metabolites = self._get_included_metabolites(metabolites)
        reactions = self._get_included_reactions(reactions)
        Keq_strs = reactions.list_attr("Keq_str")

        # Initialize lists for x_vars and x_data for QP objective
        x_vars = []
        x_data = []
        for var in self.variables:
            if var.name in metabolites:
                # Set bounds if metabolite
                met = metabolites.get_by_id(var.name)
                if met.initial_condition is None:
                    LOGGER.info(
                        "No initial conditions defined for '%s', will "
                        "not be included in QP objective",
                        met.id,
                    )
                    continue
                if any([m in fixed_conc_bounds for m in [met, met.id]]):
                    bounds = (met.initial_condition, met.initial_condition)
                else:
                    bounds = (0, np.inf)
                x_data.append(met.initial_condition)
                x_vars.append(var)
            elif var.name in Keq_strs:
                # Set bounds if reaction Keq variable
                rxn = reactions[Keq_strs.index(var.name)]
                if any([r in fixed_Keq_bounds for r in [rxn.id, rxn.Keq_str]]):
                    bounds = (rxn.Keq, rxn.Keq)
                else:
                    bounds = (0, np.inf)
                x_data.append(rxn.equilibrium_constant)
                x_vars.append(var)
            else:
                continue
            # Get log values of bounds and set
            bounds = _get_log_bounds(
                bounds,
                np.inf,
                kwargs.get("decimal_precision"),
                self.zero_value_log_substitute,
            )
            var.lb, var.ub = bounds

        # Set objective for feasible qp problem
        x_vars = Matrix(x_vars)
        x_data = Matrix(np.log(x_data))
        F = 2 * eye(len(x_vars))

        # QP Objective must be in form of 0.5 * x.T * F * x - c.T * x
        obj = 0.5 * x_vars.T * F * x_vars - (2 * x_data).T * x_vars

        self.objective = obj[0]
        self.objective_direction = "min"
        self.problem_type = "feasible_qp"

    def choose_solver(self, solver=None, qp=False):
        """Choose a solver given a solver name.

        This will choose a solver compatible with the
        :class:`.ConcSolver` and required capabilities.

        Also respects :attr:`ConcSolver.solver` where it can.

        Parameters
        ----------
        solver : str
            The name of the solver to be used.
        qp : boolean
            Whether the solver needs Quadratic Programming capabilities.
            Default is ``False``.

        Returns
        -------
        solver : optlang.interface.Model
            Returns a valid solver for the problem.

        Raises
        ------
        SolverNotFound
            If no suitable solver could be found.

        """
        if solver is None:
            solver = self.problem
        else:
            self.solver = solver

        # Check for QP, raise error if no QP solver found
        if qp and interface_to_str(solver) not in qp_solvers:
            solver = solvers[get_solver_name(qp=True)]

        return solver

    def add_cons_vars(self, what, **kwargs):
        """Add constraints and variables to the solver's problem.

        Useful for variables and constraints that cannot be expressed through
        metabolite concentrations and reaction equilibrium constants.

        Additions are reversed upon exit if the solver itself is used as
        context.

        Parameters
        ----------
        what : list, tuple
           Either a ``list`` or a ``tuple`` of variables or constraints to
           add to the solver. Must be of :class:`optlang.interface.Variable`
           or :class:`optlang.interface.Constraint`.
        **kwargs : keyword arguments
           Passed to solver.add().

        """
        context = get_context(self)

        self.solver.add(what, **kwargs)
        if context:
            context(partial(self.solver.remove, what))

    def remove_cons_vars(self, what):
        """Remove constraints and variables from the solver's problem.

        Remove variables and constraints that were added directly to the
        solver's underlying mathematical problem. Removals are reversed
        upon exit if the model itself is used as context.

        Parameters
        ----------
        what : list, tuple
           Either a ``list`` or a ``tuple`` of variables or constraints to
           remove from the solver. Must be of
           :class:`optlang.interface.Variable` or
           :class:`optlang.interface.Constraint`.

        """
        context = get_context(self)

        self.solver.remove(what)
        if context:
            context(partial(self.solver.add, what))

    def reset_constraints(self):
        """Reset the constraints.

        Ensures constraints are updated if solver's problem changes in any way.

        """
        self.remove_cons_vars(list(self.constraints))
        # Get reaction variables for the solver,
        # filtering out those to be excluded.
        reactions = self._get_included_reactions()
        for rxn in reactions:
            if rxn.Keq_str in self.variables:
                # Add concentration Keq log constraint
                self.add_concentration_Keq_cons_to_problem(rxn)

    def add_excluded_metabolites(self, metabolites, reset_problem=False):
        """Add metabolites to the exclusion list for problem creation.

        Note that this will not remove metabolites from the current problem.
        The problem must first be reset in order for changes to take effect.

        Parameters
        ----------
        metabolites : iterable
            An iterable of :class:`.MassMetabolite` objects or their
            identifiers to be added to the
            :attr:`~ConcSolver.excluded_metabolites`.
        reset_problem : bool
            Whether to reset the underlying mathematical problem of the solver
            to a generic one after adding additional metabolites to exclude.
            If ``False`` then it is incumbent upon the user to remove the
            newly excluded metabolites from the solver.

        """
        metabolites = ensure_iterable(metabolites)
        # First check whether the metabolites exist in the model
        bad_ids = [m for m in metabolites if str(m) not in self.model.metabolites]
        if bad_ids:
            raise ValueError("Invalid metabolite identifiers in {0!r}".format(bad_ids))

        metabolites = self.model.metabolites.get_by_any(metabolites)
        # Add new metabolites to excluded list, ignoring those that exist
        for met in metabolites:
            if met.id in self.excluded_metabolites:
                warn("Metabolite '{0}' already excluded.".format(met.id))
            else:
                self.excluded_metabolites.append(met.id)

        if reset_problem:
            self._initialize_solver(reset_problem=reset_problem)

    def remove_excluded_metabolites(self, metabolites, reset_problem=False):
        """Remove metabolites from the exclusion list for problem creation.

        Note that this will not add metabolites to the current problem.
        The problem must first be reset in order for changes to take effect.

        Parameters
        ----------
        metabolites : iterable
            An iterable of :class:`.MassMetabolite` objects or their
            identifiers to be removed from the
            :attr:`~ConcSolver.excluded_metabolites`.
        reset_problem : bool
            Whether to reset the underlying mathematical problem of the solver
            to a generic one after removing additional metabolites to exclude.
            If ``False`` then it is incumbent upon the user to add the
            newly included metabolites to the solver.

        """
        metabolites = ensure_iterable(metabolites)
        # First check whether the metabolites exist in the model
        bad_ids = [m for m in metabolites if str(m) not in self.excluded_metabolites]
        if bad_ids:
            raise ValueError("Invalid metabolite identifiers in {0!r}".format(bad_ids))

        metabolites = self.model.metabolites.get_by_any(metabolites)
        # Remove metabolites from excluded list
        for met in metabolites:
            self.excluded_metabolites.remove(met.id)

        if reset_problem:
            self._initialize_solver(reset_problem=reset_problem)

    def add_excluded_reactions(self, reactions, reset_problem=False):
        """Add reactions to the exclusion list for problem creation.

        Note that this will not remove reaction equilibrium constants or
        constraints from the current problem. The problem must first be reset
        in order for changes to take effect.

        Parameters
        ----------
        reactions : iterable
            An iterable of :class:`.MassReaction` objects or their
            identifiers to be added to the
            :attr:`~ConcSolver.excluded_reactions`.
        reset_problem : bool
            Whether to reset the underlying mathematical problem of the solver
            to a generic one after adding additional reactions to exclude.
            If ``False`` then it is incumbent upon the user to remove the
            newly excluded reactions from the solver.

        """
        reactions = ensure_iterable(reactions)
        # First check whether the reactions exist in the model
        bad_ids = [
            r for r in reactions if getattr(r, "id", str(r)) not in self.model.reactions
        ]
        if bad_ids:
            raise ValueError("Invalid reaction identifiers in {0!r}".format(bad_ids))

        reactions = self.model.reactions.get_by_any(reactions)

        # Add new reactions to excluded list, ignoring those that exist
        for rxn in reactions:
            if rxn.id in self.excluded_reactions:
                warn("Reaction '{0}' already excluded.".format(rxn.id))
            else:
                self.excluded_reactions.append(rxn.id)

        if reset_problem:
            self._initialize_solver(reset_problem=reset_problem)

    def remove_excluded_reactions(self, reactions, reset_problem=False):
        """Remove reactions from the exclusion list for problem creation.

        Note that this will not add reaction equilibrium constants or
        constraints to the current problem. The problem must first be reset
        in order for changes to take effect.

        Parameters
        ----------
        reactions : iterable
            An iterable of :class:`.MassReaction` objects or their
            identifiers to be removed from the
            :attr:`~ConcSolver.excluded_reactions`.
        reset_problem : bool
            Whether to reset the underlying mathematical problem of the solver
            to a generic one after removing additional reactions to exclude.
            If ``False`` then it is incumbent upon the user to add the
            newly included reaction equilibrium constants and constraints
            to the solver.

        """
        reactions = ensure_iterable(reactions)
        # First check whether the reactions exist in the model
        bad_ids = [
            r
            for r in reactions
            if getattr(r, "id", str(r)) not in self.excluded_reactions
        ]

        if bad_ids:
            raise ValueError("Invalid reaction identifiers in {0!r}".format(bad_ids))

        reactions = self.model.metabolites.get_by_any(reactions)

        # Remove reactions from excluded list
        for rxn in reactions:
            self.excluded_reactions.remove(rxn.id)

        if reset_problem:
            self._initialize_solver(reset_problem=reset_problem)

    def add_equilibrium_reactions(self, reactions, reset_problem=False):
        """Add additional reaction to the equilibrium reaction list.

        The problem must first be reset in order for changes to take effect
        as a result of adding additional equilibrium reactions.

        Parameters
        ----------
        reactions : iterable
            An iterable of :class:`.MassReaction` objects or their
            identifiers to be added to the
            :attr:`~ConcSolver.equilibrium_reactions`.
        reset_problem : bool
            Whether to reset the underlying mathematical problem of the solver
            to a generic one after adding additional equilibrium reactions.
            If ``False`` then it is incumbent upon the user to make the
            changes necessary for the mathematical problem of the solver.

        """
        reactions = ensure_iterable(reactions)
        # First check whether the reactions exist in the model
        bad_ids = [
            r for r in reactions if getattr(r, "id", str(r)) not in self.model.reactions
        ]
        if bad_ids:
            raise ValueError("Invalid reaction identifiers in {0!r}".format(bad_ids))

        reactions = self.model.reactions.get_by_any(reactions)
        # Add new reactions to excluded list, ignoring those that exist
        for rxn in reactions:
            if rxn.id in self.equilibrium_reactions:
                warn(
                    "Reaction '{0}' already exists as an equilibrium "
                    "reaction.".format(rxn.id)
                )
            else:
                self.equilibrium_reactions.append(rxn.id)

        if reset_problem:
            self._initialize_solver(reset_problem=reset_problem)

    def remove_equilibrium_reactions(self, reactions, reset_problem=False):
        """Remove reactions from the equilibrium reaction list.

        The problem must first be reset in order for changes to take effect
        as a result of removing equilibrium reactions.

        Parameters
        ----------
        reactions : iterable
            An iterable of :class:`.MassReaction` objects or their
            identifiers to be removed from the
            :attr:`~ConcSolver.equilibrium_reactions`.
        reset_problem : bool
            Whether to reset the underlying mathematical problem of the solver
            to a generic one after removing equilibrium reactions.
            If ``False`` then it is incumbent upon the user to make the
            changes necessary for the mathematical problem of the solver.

        """
        reactions = ensure_iterable(reactions)
        # First check whether the reactions exist in the model
        bad_ids = [
            r
            for r in reactions
            if getattr(r, "id", str(r)) not in self.equilibrium_reactions
        ]

        if bad_ids:
            raise ValueError("Invalid reaction identifiers in {0!r}".format(bad_ids))

        reactions = self.model.metabolites.get_by_any(reactions)

        # Remove reactions from equilibrium list
        for rxn in reactions:
            self.equilibrium_reactions.remove(rxn.id)

        if reset_problem:
            self._initialize_solver(reset_problem=reset_problem)

    def optimize(self, objective_sense=None, raise_error=False, **kwargs):
        """Optimize the :class:`ConcSolver`.

        Notes
        -----
        Only the most commonly used parameters are presented here.  Additional
        parameters for solvers may be available and specified with the
        appropriate keyword argument.

        Parameters
        ----------
        objective_sense : str or None
            Either ``"maximize"`` or ``"minimize"`` indicating whether
            variables should be maximized or minimized. In case of ``None``,
            the previous direction is used.
        raise_error : bool
            If ``True``, raise an :class:`~cobra.exceptions.OptimizationError`
            if solver status is not optimal.
        **kwargs
            decimal_precision :
                ``bool`` indicating whether to apply the
                :attr:`~.MassBaseConfiguration.decimal_precision` attribute of
                the :class:`.MassConfiguration` to the solution values.

                Default is ``False``.

        """
        kwargs = _check_kwargs(
            {
                "decimal_precision": False,
            },
            kwargs,
        )
        original_direction = self.objective.direction
        self.objective.direction = {"maximize": "max", "minimize": "min"}.get(
            objective_sense, original_direction
        )
        # Optimize
        self.slim_optimize(**kwargs)
        # Format the solution as a ConcSolution object
        solution = get_concentration_solution(self, raise_error=raise_error, **kwargs)
        # Reset the objective direction
        self.objective.direction = original_direction
        return solution

    def slim_optimize(self, error_value=float("nan"), message=None, **kwargs):
        """Optimize model without creating a :class:`.ConcSolution` object.

        Creating a full solution object implies fetching shadow prices and
        solution values for all variables and constraints in the solver.

        This necessarily takes some time and in cases where only one
        or two values are of interest, it is recommended to instead use this
        function which does not create a solution object, returning only the
        value of the objective.

        Note however that the :meth:`~ConcSolver.slim_optimize` method uses
        efficient means to fetch values so if you need solution values or
        shadow prices for more than say 4 metabolites/reactions, then the
        total speed increase of :meth:`~ConcSolver.slim_optimize` versus
        :meth:`~ConcSolver.optimize` is expected to be small or even negative
        depending on other methods of fetching values after optimization.

        Parameters
        ----------
        error_value : float or None
           The value to return if optimization failed due to e.g.
           infeasibility. If ``None``, raise
           :class:`~cobra.exceptions.OptimizationError` if the
           optimization fails.
        message : string
           Error message to use if the optimization did not succeed.
        **kwargs
            decimal_precision :
                ``bool`` indicating whether to apply the
                :attr:`~.MassBaseConfiguration.decimal_precision` attribute of
                the :class:`.MassConfiguration` to the solution values.

                Default is ``False``.

        Returns
        -------
        float
            The objective value.

        """
        kwargs = _check_kwargs(
            {
                "decimal_precision": False,
            },
            kwargs,
        )
        self.solver.optimize()
        status = self.solver.status
        if status == optlang.interface.OPTIMAL:
            value = np.exp(self.solver.objective.value)
            if kwargs.get("decimal_precision"):
                value = apply_decimal_precision(
                    value, MASSCONFIGURATION.decimal_precision
                )
            return value

        if error_value is None:
            exception_cls = OPTLANG_TO_EXCEPTIONS_DICT.get(status, OptimizationError)
            raise exception_cls("{0} ({1})".format(message, status))

        return error_value

    def set_objective(self, value, additive=False):
        """Set the model objective.

        Parameters
        ----------
        value : optlang.interface.Objective, ~sympy.core.basic.Basic, or dict
            If the objective is linear, the value can be a new
            :class:`optlang.interface.Objective` object or a ``dict`` with
            linear coefficients where each key is a metabolite and the element
            the new coefficient (``float``).
        additive : bool
            If ``True``, add the terms to the current objective, otherwise
            start with an empty objective.

        """
        interface = self.problem
        reverse_value = self.solver.objective.expression
        reverse_value = interface.Objective(
            reverse_value, direction=self.solver.objective.direction, sloppy=True
        )

        if isinstance(value, dict):
            if not self.objective.is_Linear:
                raise ValueError(
                    "Can only update non-linear objectives additively using "
                    "object of class self.problem.Objective, not {0}".format(
                        type(value)
                    )
                )

            if not additive:
                self.solver.objective = interface.Objective(
                    Zero, direction=self.solver.objective.direction
                )

            self.solver.objective.set_linear_coefficients(value)

        elif isinstance(value, (Basic, optlang.interface.Objective)):
            if isinstance(value, Basic):
                value = interface.Objective(
                    value, direction=self.solver.objective.direction, sloppy=False
                )
            # Check whether expression only uses variables from current
            # clone the objective if not, faster than cloning without checking
            atoms = value.expression.atoms(optlang.interface.Variable)
            if all(a.problem is self.solver for a in atoms):
                value = interface.Objective.clone(value, model=self.solver)

            if not additive:
                self.solver.objective = value
            else:
                self.solver.objective += value.expression
        else:
            raise TypeError(
                "{0!r} is not a valid objective for {1!r}.".format(value, self.solver)
            )

        context = get_context(self)
        if context:

            def reset():
                self.model.solver.objective = reverse_value
                self.solver.objective.direction = reverse_value.direction

            context(reset)

    def add_metabolite_var_to_problem(
        self, metabolite, lower_bound=None, upper_bound=None, **kwargs
    ):
        """Add a metabolite concentration variable to the problem.

        The variable in linear space is represented as::

            log(x_lb) <= log(x) <= log(x_ub)

        where

            * ``x`` is the metabolite concentration variable.
            * ``x_lb`` is the lower bound for the concentration.
            * ``x_ub`` is the upper bound for the concentration.

        Parameters
        ----------
        metabolite : MassMetabolite
            The metabolite whose concentration should be added as a variable
            to the solver.
        lower_bound : float
            A non-negative number for the lower bound of the variable. If
            ``bound_type=deviation`` then value is treated as a percentage.
            Otherwise value is treated as the lower bound in linear space.
        upper_bound : float
            A non-negative number for the upper bound of the variable. If
            ``bound_type=deviation`` then value is treated as a percentage.
            Otherwise value is treated as the lower bound in linear space.
        **kwargs
            bound_type :
                Either ``"deviation"`` to indicate that bound values are
                percentages representing deviations of the base concentration
                value, or ``"absolute"`` to indicate that the bound values
                should be not be treated as percentages but as the bound values
                themselves (in linear space).

                Default is ``deviation``.
            concentration :
                A non-negative number to treat as the base concentration value
                in setting percent deviation bounds. Ignored if
                ``bound_type=absolute``

                Default is the current metabolite concentration accessed via
                :attr:`.MassMetabolite.initial_condition`.
            decimal_precision :
                ``bool`` indicating whether to apply the
                :attr:`~.MassBaseConfiguration.decimal_precision` attribute of
                the :class:`.MassConfiguration` to the bound values.

                Default is ``False``.

        """
        kwargs = _check_kwargs(
            {
                "concentration": metabolite.initial_condition,
                "bound_type": "deviation",
                "decimal_precision": False,
            },
            kwargs,
        )
        # Create and return log variable
        variable = self._create_variable(metabolite, lower_bound, upper_bound, **kwargs)
        self.add_cons_vars([variable])

        return variable

    def add_reaction_Keq_var_to_problem(
        self, reaction, lower_bound=None, upper_bound=None, **kwargs
    ):
        """Add a reaction equilibrium constant variable to the problem.

        The variable in linear space is represented as::

            log(Keq_lb) <= log(Keq) <= log(Keq_ub)

        where

            * ``Keq`` is the equilibrium constant variable for the reaction.
            * ``Keq_lb`` is the lower bound for the equilibrium constant.
            * ``Keq_ub`` is the upper bound for the equilibrium constant.

        Parameters
        ----------
        reaction : MassReaction
            The reaction whose equilibrium constant should be added as a
            variable to the solver.
        lower_bound : float
            A non-negative number for the lower bound of the variable. If
            ``bound_type=deviation`` then value is treated as a percentage.
            Otherwise value is treated as the lower bound in linear space.
        upper_bound : float
            A non-negative number for the upper bound of the variable. If
            ``bound_type=deviation`` then value is treated as a percentage.
            Otherwise value is treated as the lower bound in linear space.
        **kwargs
            bound_type :
                Either ``"deviation"`` to indicate that bound values are
                percentages representing deviations of the base equilibrium
                constant value, or ``"absolute"`` to indicate that the
                bound values should be not be treated as percentages but as
                the bound values themselves (in linear space).

                Default is ``deviation``.
            Keq :
                A non-negative number to treat as the base equilibrium
                constant value in setting percent deviation bounds. Ignored if
                ``bound_type=absolute``

                Default is the current reaction equilibrium constant accessed
                via :attr:`.MassReaction.equilibrium_constant`.
            steady_state_flux :
                The steady state flux value of the reaction. If set as 0.,
                the creation of the equilibrium constant variable will depend
                on whether the reaction is defined as an equilibrium reaction.

                Default is the current reaction steady state flux accessed
                via :attr:`.MassReaction.steady_state_flux`.
            decimal_precision :
                ``bool`` indicating whether to apply the
                :attr:`~.MassBaseConfiguration.decimal_precision` attribute of
                the :class:`.MassConfiguration` to the bound values.

                Default is ``False``.

        """
        kwargs = _check_kwargs(
            {
                "Keq": reaction.Keq,
                "steady_state_flux": reaction.steady_state_flux,
                "decimal_precision": False,
                "bound_type": "deviation",
            },
            kwargs,
        )

        steady_state_flux = kwargs.pop("steady_state_flux")
        if kwargs.get("decimal_precision"):
            steady_state_flux = round(
                steady_state_flux, MASSCONFIGURATION.decimal_precision
            )
        if steady_state_flux == 0 and reaction.id not in self.equilibrium_reactions:
            raise ValueError(
                "Steady state flux is at equilibrium for reaction '{0}' and "
                "not in `equilibrium_reactions` attribute".format(reaction.id)
            )

        # Create and return log variable
        variable = self._create_variable(reaction, lower_bound, upper_bound, **kwargs)
        self.add_cons_vars([variable])

        return variable

    def add_concentration_Keq_cons_to_problem(self, reaction, epsilon=None, **kwargs):
        """Add constraint using the reaction metabolite stoichiometry and Keq.

        The constraint in linear space is represented as::

            S^T * log(x) <= log(Keq) - epsilon if v > 0
            S^T * log(x) >= log(Keq) + epsilon if v < 0

        where

            * ``S^T`` is the transposed stoichiometry for the reaction.
            * ``x`` is the vector of concentration variables.
            * ``Keq`` is the equilibrium constant variable for the reaction.
            * ``v`` is the steady state flux value for the reaction.
            * ``epsilon`` is a buffer for the constraint .

        Parameters
        ----------
        reaction : MassReaction
            The reaction whose metabolite stoichiometry and equilibrium
            constant is used to create the constraint to be added.
        epsilon : float
            The buffer for the constraint.

            Default is ``None`` to use :attr:`ConcSolver.constraint_buffer`.
        **kwargs
            bound_type :
                Either ``"deviation"`` to indicate that bound values are
                percentages representing deviations of the base equilibrium
                constant value, or ``"absolute"`` to indicate that the
                bound values should be not be treated as percentages but as
                the bound values themselves (in linear space). Only used
                if the variables necessary for the constraint do not already
                exist.

                Default is ``deviation``.
            steady_state_flux :
                The steady state flux value of the reaction. Determines whether
                the contraint is set up as less than or as greater than.
                If set as 0.,  the creation of the reaction constraint will
                depend on whether the reaction is defined as an equilibrium
                reaction.

                Default is the current reaction steady state flux accessed
                via :attr:`.MassReaction.steady_state_flux`.
            decimal_precision :
                ``bool`` indicating whether to apply the
                :attr:`~.MassBaseConfiguration.decimal_precision` attribute of
                the :class:`.MassConfiguration` to the bound values.

                Default is ``False``.

        """
        kwargs = _check_kwargs(
            {
                "bound_type": "deviation",
                "steady_state_flux": reaction.steady_state_flux,
                "decimal_precision": False,
            },
            kwargs,
        )

        steady_state_flux = kwargs.get("steady_state_flux")
        if steady_state_flux is None:
            warn(
                "No flux defined for '{0}', no constraint "
                "can be made.".format(reaction.id)
            )
            return None

        if kwargs.get("decimal_precision"):
            steady_state_flux = round(
                steady_state_flux, MASSCONFIGURATION.decimal_precision
            )

        if epsilon is None:
            epsilon = self.constraint_buffer

        if steady_state_flux > 0:
            lb, ub = (None, -1 * epsilon)
        elif steady_state_flux < 0:
            lb, ub = (1 * epsilon, None)
        elif reaction.id in self.equilibrium_reactions:
            lb, ub = (0, 0)
        else:
            raise ValueError(
                "Steady state flux is at equilibrium for reaction '{0}' and "
                "not in `equilibrium_reactions` attribute".format(reaction.id)
            )

        var_and_coeffs = {}
        for metabolite, coefficient in iteritems(reaction.metabolites):
            if metabolite.id in self.excluded_metabolites:
                continue
            try:
                var = self.variables[metabolite.id]
            except KeyError:
                var = self.add_metabolite_var_to_problem(
                    metabolite, bound_type=kwargs.get("bound_type")
                )
            var_and_coeffs[var] = coefficient

        try:
            var = self.variables[reaction.Keq_str]
        except KeyError:
            var = self.add_reaction_Keq_var_to_problem(
                reaction, bound_type=kwargs.get("bound_type")
            )
        var_and_coeffs[var] = -1

        constraint = self.problem.Constraint(Zero, name=reaction.id, lb=lb, ub=ub)
        self.add_cons_vars([constraint], sloppy=True)
        constraint = self.constraints[reaction.id]
        constraint.set_linear_coefficients(var_and_coeffs)

        return constraint

    def update_model_with_solution(self, concentration_solution, **kwargs):
        """Update :attr:`ConcSolver.model` using a :class:`.ConcSolution`.

        Parameters
        ----------
        concentration_solution : ConcSolution
            The :class:`.ConcSolution` containing the solution values to use
            in updating the model.
        **kwargs
            concentrations :
                ``bool`` indicating whether to update the metabolite
                concentrations of the model (the
                :attr:`.MassMetabolite.initial_condition` values).

                Default is ``True``.
            Keqs :
                ``bool`` indicating whether to update the reaction equilibrium
                constants of the model (the
                :attr:`.MassReaction.equilibrium_constant` values).

                Default is ``True``.
            inplace :
                ``bool`` indicating whether to modify the current
                :attr:`ConcSolver.model` or to copy the model, then modify and
                set the model copy as the new :attr:`ConcSolver.model`. If
                ``False``, the association with the old model is removed and
                an association with the new model is created.

                Default is ``True``.

        """
        kwargs = _check_kwargs(
            {"concentrations": True, "Keqs": True, "inplace": True}, kwargs
        )
        if not kwargs.get("inplace"):
            # Remove association to prevent excessive replication
            setattr(self.model, "_conc_solver", None)
            self._model = update_model_with_concentration_solution(
                self.model, concentration_solution, **kwargs
            )
            setattr(self.model, "_conc_solver", self)
        else:
            update_model_with_concentration_solution(
                self.model, concentration_solution, **kwargs
            )

    def _create_variable(self, mass_obj, lower_bound, upper_bound, **kwargs):
        """Create an optlang variable for the solver.

        Warnings
        --------
        This method is intended for internal use only.

        """
        # Determine what type of variable to make
        if "Reaction" in mass_obj.__class__.__name__:
            var_id = mass_obj.Keq_str
            var_type = "Keq"
            upper_default = MASSCONFIGURATION.irreversible_Keq
        else:
            var_id = mass_obj.id
            var_type = "concentration"
            upper_default = np.inf

        # Validate bound type
        bound_type = _validate_bound_type(kwargs.get("bound_type"))
        if bound_type in {"deviation", "dev"}:
            value = kwargs.get(var_type)
            if value < 0 or value is None:
                raise ValueError(
                    "'{0}' {1} must be a non-negative number".format(var_id, var_type)
                )
        # Get deviation values
        bounds = lower_bound, upper_bound
        if bound_type in {"deviation", "dev"}:
            bounds = _get_deviation_values(value, *bounds)

        bounds = _get_log_bounds(
            bounds,
            upper_default,
            kwargs.get("decimal_precision"),
            self.zero_value_log_substitute,
        )

        # Create variable and return
        variable = self.problem.Variable(var_id, lb=bounds[0], ub=bounds[1])
        return variable

    def _initialize_solver(self, reset_problem=False, **kwargs):
        """Initialize the solver as a generic problem involving concentraitons.

        Warnings
        --------
        This method is intended for internal use only.

        """
        if reset_problem:
            interface = MASSCONFIGURATION.solver
            self._solver = interface.Model()
            self.solver.objective = interface.Objective(Zero)
            self.problem_type = ""

        # Get fixed variables
        fixed_conc_bounds = kwargs.pop("fixed_conc_bounds")
        fixed_Keq_bounds = kwargs.pop("fixed_Keq_bounds")

        # Get concentration variables for the solver,
        # filtering out those to be excluded.
        metabolites = self._get_included_metabolites()
        for met in metabolites:
            if any([m in fixed_conc_bounds for m in [met, met.id]]):
                lb, ub = (0, 0)
                bound_type = "deviation"
            else:
                lb, ub = (0, np.inf)
                bound_type = "absolute"
            kwargs["bound_type"] = bound_type
            self.add_metabolite_var_to_problem(met, lb, ub, **kwargs)

        # Get reaction variables for the solver,
        # filtering out those to be excluded.
        ignored_rxns = []
        reactions = self._get_included_reactions()
        for rxn in reactions:
            if rxn.id in fixed_Keq_bounds or rxn.Keq_str in fixed_Keq_bounds:
                lb, ub = (0, 0)
                bound_type = "deviation"
            else:
                lb, ub = (0, np.inf)
                bound_type = "absolute"
            kwargs["bound_type"] = bound_type
            try:
                # Add Keq variable to problem
                self.add_reaction_Keq_var_to_problem(rxn, lb, ub, **kwargs)
                # Add concentration Keq log constraint
                self.add_concentration_Keq_cons_to_problem(rxn)
            except ValueError:
                ignored_rxns += [rxn.id]

        if ignored_rxns:
            LOGGER.warning(
                "No Keq variables or reaction constraints created for the "
                "following reactions not defined as equilibrium reactions but "
                "with steady state fluxes of 0.:\n{0!r}".format(ignored_rxns)
            )

        self.problem_type = "generic"

    def _get_included_metabolites(self, metabolites=None):
        """Return a list of MassMetabolite objects for included metabolites.

        Warnings
        --------
        This method is intended for internal use only.

        """
        # If None, set up for all included model metabolites
        if metabolites is None:
            metabolites = self.included_metabolites
        else:
            # Get concentration variables for the solver,
            # filtering out those to be excluded.
            metabolites = [
                m
                for m in ensure_iterable(metabolites)
                if m.id not in self.excluded_metabolites
            ]
        metabolites = self.model.metabolites.get_by_any(self.included_metabolites)

        return DictList(metabolites)

    def _get_included_reactions(self, reactions=None):
        """Return a list of MassReaction objects for included reactions.

        Warnings
        --------
        This method is intended for internal use only.

        """
        # If None, set up for all included model reactions
        if reactions is None:
            reactions = self.included_reactions
        else:
            # Get Keq variables for the solver,
            # filtering out those to be excluded.
            reactions = [
                r
                for r in ensure_iterable(reactions)
                if r.id not in self.excluded_reactions
            ]
        reactions = self.model.reactions.get_by_any(self.included_reactions)

        return DictList(reactions)

    def _check_for_missing_values(
        self,
        model,
        concentrations=True,
        equilibrium_constants=True,
        steady_state_fluxes=True,
    ):
        """Determine missing values that prevent setup of a problem.

        Warnings
        --------
        This method is intended for internal use only.

        """
        missing = {}
        if concentrations:
            missing["concentrations"] = [
                m.id
                for m in get_missing_initial_conditions(
                    model,
                    metabolite_list=[
                        met
                        for met in model.metabolites
                        if met.id not in self.excluded_metabolites
                    ],
                )
            ]
        if equilibrium_constants:
            missing["equilibrium constants"] = [
                r.id
                for r, v in iteritems(
                    get_missing_reaction_parameters(
                        model,
                        reaction_list=[
                            rxn
                            for rxn in model.reactions
                            if rxn.id not in self.excluded_reactions
                        ],
                    )
                )
                if "Keq" in v
            ]

        if steady_state_fluxes:
            missing["steady state fluxes"] = [
                r.id
                for r in get_missing_steady_state_fluxes(
                    model,
                    reaction_list=[
                        rxn
                        for rxn in model.reactions
                        if rxn.id not in self.excluded_reactions
                    ],
                )
            ]

        return missing

    def __dir__(self):
        """Override default dir() implementation to list only public items.

        Warnings
        --------
        This method is intended for internal use only.

        """
        return get_public_attributes_and_methods(self)

    def __repr__(self):
        """Override default repr.

        Warnings
        --------
        This method is intended for internal use only.

        """
        return "<%s %s at 0x%x>" % (self.__class__.__name__, self.model.id, id(self))

    def __enter__(self):
        """Record all future changes, undoing them when __exit__ is called.

        Warnings
        --------
        This method is intended for internal use only.

        """
        self.model.__enter__()

    def __exit__(self, type, value, traceback):
        """Pop the top context manager and trigger the undo functions.

        Warnings
        --------
        This method is intended for internal use only.

        """
        self.model.__exit__(type, value, traceback)


def concentration_constraint_matricies(
    concentration_solver, array_type="dense", zero_tol=1e-6
):
    """Create a matrix representation of the problem.

    This is used for alternative solution approaches that do not use optlang.
    The function will construct the equality matrix, inequality matrix and
    bounds for the complete problem.

    Parameters
    ---------
    concentration_solver : ConcSolver
        The :class:`.ConcSolver` containing the mathematical problem.
    array_type : str
        A string identifiying the desired format for the returned matrix.
        Valid matrix types include ``'dense'``, ``'dok'``, ``'lil'``,
        ``'DataFrame'``, and ``'symbolic'``
        Default is the current ``dense``. See the :mod:`~.matrix`
        module documentation for more information on the ``array_type``.
    zero_tol : float
        The zero tolerance used to judge whether two bounds are the same.

    Returns
    -------
    collections.namedtuple
        A named tuple consisting of 6 matrices and 2 vectors:

        - ``"equalities"`` is a matrix ``S`` such that ``S*vars = b``.
          It includes a row for each equality constraint and a column for
          each variable.
        - ``"b"`` the right side of the equality equation such that
          ``S*vars = b.``
        - ``"inequalities"`` is a matrix ``M`` such that
          ``lb <= M*vars <= ub``. It contains a row for each inequality and
          as many columns as variables.
        - ``"bounds"`` is a compound matrix ``[lb ub]`` containing the lower
          and upper bounds for the inequality constraints in ``M``.
        - ``"variable_fixed"`` is a boolean vector indicating whether the
          variable at that index is fixed (``lower bound == upper_bound``) and
          is thus bounded by an equality constraint.
        - ``"variable_bounds"`` is a compound matrix ``[lb ub]`` containing
          the lower and upper bounds for all variables.

    """
    try:
        matrix_constructor = {
            "dense": np.array,
            "dok": dok_matrix,
            "lil": lil_matrix,
            "DataFrame": pd.DataFrame,
            "symbolic": Matrix,
        }[array_type]
    except KeyError as e:
        raise ValueError("Unrecognized `array_type` '{0}'.".format(str(e)))

    problem = namedtuple(
        "Problem",
        [
            "equalities",
            "b",
            "inequalities",
            "bounds",
            "variable_fixed",
            "variable_bounds",
        ],
    )
    equality_rows = []
    b = []
    inequality_rows = []
    inequality_bounds = []

    for constraint in concentration_solver.constraints:
        bounds = [
            -np.inf if constraint.lb is None else constraint.lb,
            np.inf if constraint.ub is None else constraint.ub,
        ]
        coefs = constraint.get_linear_coefficients(concentration_solver.variables)
        coefs = [coefs[v] for v in concentration_solver.variables]
        if (bounds[-1] - bounds[0]) < zero_tol:
            b.append(bounds[0] if abs(bounds[0]) > zero_tol else 0.0)
            equality_rows.append(coefs)
        else:
            inequality_rows.append(coefs)
            inequality_bounds.append(bounds)

    var_bounds = np.array([[v.lb, v.ub] for v in concentration_solver.variables])
    fixed = var_bounds[:, 1] - var_bounds[:, 0] < zero_tol

    return problem(
        equalities=matrix_constructor(equality_rows),
        b=np.array(b),
        inequalities=matrix_constructor(inequality_rows),
        bounds=matrix_constructor(inequality_bounds),
        variable_fixed=np.array(fixed),
        variable_bounds=matrix_constructor(var_bounds),
    )


def _validate_bound_type(bound_type):
    """Valdiate the bound type input.

    Warnings
    --------
    This method is intended for internal use only.

    """
    # Validate bound type
    if bound_type not in {"absolute", "abs", "deviation", "dev"}:
        raise ValueError("`bound_type` must be 'absolute' or 'deviation'")

    return bound_type


def _validate_bounds(lower_bound, upper_bound, prefix="", suffix="", exclude_zero=True):
    """Valdiate the bound value input.

    Warnings
    --------
    This method is intended for internal use only.

    """
    for value in [lower_bound, upper_bound]:
        # Verify value is valid
        try:
            value = ensure_non_negative_value(value, exclude_zero=exclude_zero)
        except (ValueError, TypeError) as e:
            raise e.__class__(
                "The {1}bounds{2} {0}.".format(str(e).lower(), prefix, suffix)
            )

    return lower_bound, upper_bound


def _get_deviation_values(value, lower_bound_percent, upper_bound_percent):
    """Calculate the deviation values from the given `value`.

    Warnings
    --------
    This method is intended for internal use only.

    """
    if lower_bound_percent == 0 and upper_bound_percent == 0:
        return value, value

    if lower_bound_percent is None:
        lower_bound_percent = MASSCONFIGURATION.percent_deviation
    if upper_bound_percent is None:
        upper_bound_percent = MASSCONFIGURATION.percent_deviation

    lower_bound, upper_bound = _validate_bounds(
        lower_bound_percent, upper_bound_percent, suffix=" percentage"
    )
    if value == np.inf:
        lower_bound = 0
        upper_bound = np.inf
    else:
        lower_bound = value - (value * lower_bound)
        if lower_bound < 0:
            lower_bound = 0
        upper_bound = value + (value * upper_bound)

    return lower_bound, upper_bound


def _get_log_bounds(
    bounds, upper_default, decimal_precision, zero_value_log_substitute
):
    """Calculate the bound values in logspace and return them.

    Warnings
    --------
    This method is intended for internal use only.

    """
    # Set defaults for None values
    bounds = [
        b if b is not None else default
        for b, default in zip(bounds, [0, upper_default])
    ]
    # Substitute for zero values to prevent math domain error in logspace
    bounds = [b if b != 0 else zero_value_log_substitute for b in bounds]
    # Validate bounds before taking log
    bounds = _validate_bounds(*bounds, exclude_zero=True)
    # Take log of bound values
    bounds = np.log(bounds)
    # Apply decimal precision if desired.
    if decimal_precision:
        bounds = [
            apply_decimal_precision(b, MASSCONFIGURATION.decimal_precision)
            for b in bounds
        ]

    return bounds


__all__ = ("ConcSolver", "concentration_constraint_matricies")
