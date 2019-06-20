# -*- coding: utf-8 -*-
"""TODO Module Docstrings."""
from __future__ import absolute_import

import logging

from six import integer_types, iteritems, with_metaclass

from cobra.core.configuration import Configuration
from cobra.core.singleton import Singleton
from cobra.util.solver import interface_to_str

__all__ = ("MassConfiguration",)

LOGGER = logging.getLogger(__name__)
COBRA_CONFIGURATION = Configuration()


class MassBaseConfiguration(object):
    """Define global configuration values that to be honored by mass functions.

    Attributes for model construction:
        ['boundary_compartment', 'default_compartment'. 'irreversible_Keq',
         'irreversible_kr', 'exclude_from_rates']

    Attributes for model simulation:
        ['decimal_precision',]

    Attributes for flux balance analysis (FBA):
        ['optimization_solver', 'optimization_tolerance', 'lower_bound',
         'upper_bound', 'bounds', 'processes']

    Attributes
    ----------
    boundary_compartment: dict
        A dictionary containing the identifier of the boundary compartment
        mapped to the name of the boundary compartment.
        Default value is {"b": "boundary"}.
    default_compartment: dict
        A dictionary containing the identifier of the default compartment
        mapped to the name of the desired name of default compartment. Used for
        writing models to SBML when there are no set compartments in the model.
        Default value is {"default": "default_compartment"}.
    irreversible_Keq: float
        The default value to assign to equilibrium constants (Keq) for
        irreversible reactions. Must be a non-negative value.
        Default value is the infinity (=float("inf")).
    irreversible_kr: float
        The default value to assign to equilibrium constants (Keq) for
        irreversible reactions. Must be a non-negative value.
        Default value is the 0.
    exclude_metabolites_from_rates: dict
        A dict where keys should correspond to a metabolite attrubute to
        utilize for filtering, and values are lists that contain the items to
        exclude that would be returned by the metabolite attribute. Does not
        apply to boundary reactions (MassReaction.boundary==True).
        Default is dict("elements", [{"H": 2, "O": 1}, {"H": 1}]) to remove
        the hydrogen and water metabolites using the 'elements' attribute
        to filter out the hydrogen and water in all rates except the hydrogen
        and water exchange reactions on the boundary.
    include_compartments_in_rates: bool
        Whether to include the compartment volumes in rate expressions.
        The boundary compartment will always excluded.
        Default is False.
    decimal_precision: int, None
        An integer indicating the decimal precision to use for rounding
        numerical values. Positive numbers indicated digits to the right of the
        decimal place, negative numbers indicate digits to the left of the
        decimal place. If None provided, no solutions will be rounded.
        Default is None.
    optimization_solver: {"glpk", "cplex", "gurobi"}
        The default optimization solver. The solver choices are the ones
        provided by `optlang` and solvers installed in your environment.
        Identical to the inherited `solver` attribute.
    optimization_tolerance: float
        The default tolerance for the optimization solver being used.
        Default value is 1e-7.
        Identical to the inherited `tolerance` attribute.
    lower_bound: float
        The standard lower bound for reversible reactions.
        Default value is -1000.
    upper_bound: float
        The standard upper bound for all reactions.
        Default value is 1000.
    bounds: tuple of floats
        The default reaction bounds for newly created reactions. The bounds
        are in the form of lower_bound, upper_bound.
        Default values are -1000.0, 1000.0.
    processes: int
        A default number of processes to use where multiprocessing is
        possible. The default number corresponds to the number of available
        cores (hyperthreads).

    Notes
    -----
    This object extends the BaseConfiguration class from cobra. However, in
    addition to the optimization solvers from cobrapy package, the masspy
    package utilizes ODE solvers. This may lead to confusion when trying to
    change solver options such as tolerances, since an optimization solver
    may need to utilize a different tolerance than the ODE solver. Therefore,
    the `solver` and `tolerance` attributes of the inherited BaseConfiguration
    class are renamed to `optimization_solver` and `optimization_tolerance` to
    help prevent confusion.

    """

    def __init__(self):
        """Initialize MassBaseConfiguration object."""
        # Model construction configuration options
        self._boundary_compartment = {"b": "boundary"}
        self._default_compartment = {"default": "default_compartment"}
        self._irreversible_Keq = float("inf")
        self._irreversible_kr = 0
        self.exclude_metabolites_from_rates = {
            "elements": [{"H": 2, "O": 1}, {"H": 1}]}
        self.include_compartments_in_rates = False

        # Model simulation options
        self._decimal_precision = None

        # For cobra configuration synchronization
        self._shared_state = COBRA_CONFIGURATION.__dict__

    @property
    def boundary_compartment(self):
        """Return the default value for the boundary compartment."""
        return getattr(self, "_boundary_compartment")

    @boundary_compartment.setter
    def boundary_compartment(self, value):
        """Set the default value for the boundary compartment.

        Parameters
        ----------
        value: dict
            A dictionary containing the identifier of the boundary compartment
            mapped to the name of the boundary compartment.

        """
        setattr(self, "_boundary_compartment", value)

    @property
    def default_compartment(self):
        """Return the default value for the default compartment."""
        return getattr(self, "_default_compartment")

    @default_compartment.setter
    def default_compartment(self, value):
        """Set the default value for the default compartment.

        Parameters
        ----------
        value: dict
            A dictionary containing the identifier of the default compartment
            mapped to the name of the default compartment.

        """
        setattr(self, "_default_compartment", value)

    @property
    def irreversible_Keq(self):
        """Return the default value for Keq of an irreversible reaction."""
        return getattr(self, "_irreversible_Keq")

    @irreversible_Keq.setter
    def irreversible_Keq(self, value):
        """Set the default value for Keq of an irreversible reaction.

        Parameters
        ----------
        value: float
            A non-negative number for the equilibrium constant (Keq)
            of the reaction.

        Warnings
        --------
        Equilibrium constants cannot be negative.

        """
        if value is None:
            pass
        elif not isinstance(value, (integer_types, float)):
            raise TypeError("Must be an int or float")
        elif value < 0.:
            raise ValueError("Must be a non-negative number")
        setattr(self, "_irreversible_Keq", value)

    @property
    def irreversible_kr(self):
        """Return the default value for kr of an irreversible reaction."""
        return getattr(self, "_irreversible_kr")

    @irreversible_kr.setter
    def irreversible_kr(self, value):
        """Set the default value for kr of an irreversible reaction.

        Parameters
        ----------
        value: float
            A non-negative number for the reverse rate constant (kr)
            of the reaction.

        Warnings
        --------
        Reverse rate constants cannot be negative.

        """
        if value is None:
            pass
        elif not isinstance(value, (integer_types, float)):
            raise TypeError("Must be an int or float")
        elif value < 0.:
            raise ValueError("Must be a non-negative number")
        setattr(self, "_irreversible_kr", value)

    @property
    def decimal_precision(self):
        """Return the default decimal precision when rounding."""
        return getattr(self, "_decimal_precision")

    @decimal_precision.setter
    def decimal_precision(self, value):
        """Set the default decimal precision when rounding.

        Parameters
        ----------
        value: int, None
            An integer indicating how many digits from the decimal should
            rounding occur. If None, no rounding will occur.

        """
        if value is not None and not isinstance(value, integer_types):
            raise TypeError("value must be an int.")

        setattr(self, "_decimal_precision", value)

    @property
    def optimization_solver(self):
        """Return the solver utilized for optimization."""
        return COBRA_CONFIGURATION.solver

    @optimization_solver.setter
    def optimization_solver(self, value):
        """Set the solver utilized for optimization.

        Parameters
        ----------
        value: {"glpk", "cplex", "gurobi"}
            The solver to utilize in optimizations.

        """
        COBRA_CONFIGURATION.solver = value

    @property
    def optimization_tolerance(self):
        """Return the tolerance value utilized by the optimization solver."""
        return COBRA_CONFIGURATION.tolerance

    @optimization_tolerance.setter
    def optimization_tolerance(self, value):
        """Set the tolerance value utilized by the optimization solver.

        Parameters
        ----------
        value: float
            The tolerance value to set.

        """
        COBRA_CONFIGURATION.tolerance = value

    @property
    def lower_bound(self):
        """Return the default value of the lower bound for reactions."""
        return COBRA_CONFIGURATION.lower_bound

    @lower_bound.setter
    def lower_bound(self, value):
        """Set the default value of the lower bound for reactions.

        Parameters
        ----------
        value: float
            The default bound value to set.

        """
        COBRA_CONFIGURATION.lower_bound = value

    @property
    def upper_bound(self):
        """Return the default value of the lower bound for reactions."""
        return COBRA_CONFIGURATION.upper_bound

    @upper_bound.setter
    def upper_bound(self, value):
        """Set the default value of the lower bound for reactions.

        Parameters
        ----------
        value: float
            The default bound value to set.

        """
        COBRA_CONFIGURATION.upper_bound = value

    @property
    def bounds(self):
        """Return the default lower and upper bounds for reactions."""
        return COBRA_CONFIGURATION.bounds

    @bounds.setter
    def bounds(self, bounds):
        """Set the default lower and upper bounds for reactions."""
        COBRA_CONFIGURATION.bounds = bounds

    @property
    def processes(self):
        """Return the default number of processes to use when possible."""
        return COBRA_CONFIGURATION.processes

    @property
    def shared_state(self):
        """Return a read-only dict for shared configuration attributes."""
        shared_state = {}
        for k, v in iteritems(self._shared_state):
            if k in ["_solver", "tolerance"]:
                k = "optimization_" + k.strip("_")
            shared_state[k] = v

        return shared_state

    def __repr__(self):
        """Return the representation of the MassConfiguration."""
        return """MassConfiguration:
        boundary compartment: {boundary_compartment}
        default compartment: {default_compartment}
        irreversible reaction Keq: {irreversible_Keq}
        irreversible reaction kr: {irreversible_kr}
        optimization solver: {optimization_solver}
        optimization solver tolerance: {optimization_tolerance}
        lower_bound: {lower_bound}
        upper_bound: {upper_bound}
        processes: {processes}""".format(
            boundary_compartment=[
                "{0} ({1})".format(v, k) if v else k for k, v in iteritems(
                    self.boundary_compartment)][0],
            default_compartment=[
                "{0} ({1})".format(v, k) if v else k for k, v in iteritems(
                    self.default_compartment)][0],
            irreversible_Keq=self.irreversible_Keq,
            irreversible_kr=self.irreversible_kr,
            optimization_solver=interface_to_str(self.optimization_solver),
            optimization_tolerance=self.optimization_tolerance,
            lower_bound=self.lower_bound,
            upper_bound=self.upper_bound,
            processes=self.processes)

    def _repr_html_(self):
        return """
        <table>
            <tr><tr>
                <td><strong>Boundary Compartment</strong></td>
                <td>{boundary_compartment}</td>
            </tr><tr>
                <td><strong>Default Compartment</strong></td>
                <td>{default_compartment}</td>
            </tr><tr>
                <td><strong>Irreversible Reaction Keq</strong></td>
                <td>{irreversible_Keq}</td>
            </tr><tr>
                <td><strong>Irreversible Reaction kr</strong></td>
                <td>{irreversible_kr}</td>
            </tr>
                <td><strong>Optimization solver</strong></td>
                <td>{optimization_solver}</td>
            </tr><tr>
                <td><strong>Optimization solver tolerance</strong></td>
                <td>{optimization_tolerance}</td>
            </tr><tr>
                <td><strong>Lower bound</strong></td>
                <td>{lower_bound}</td>
            </tr><tr>
                <td><strong>Upper bound</strong></td>
                <td>{upper_bound}</td>
            </tr><tr>
                <td><strong>Processes</strong></td>
                <td>{processes}</td>
            </tr>
        </table>""".format(
            boundary_compartment=[
                "{0}: {1}".format(k, v) if v else k for k, v in iteritems(
                    self.boundary_compartment)][0],
            default_compartment=[
                "{0}: {1}".format(k, v) if v else k for k, v in iteritems(
                    self.default_compartment)][0],
            irreversible_Keq=self.irreversible_Keq,
            irreversible_kr=self.irreversible_kr,
            optimization_solver=interface_to_str(self.optimization_solver),
            optimization_tolerance=self.optimization_tolerance,
            lower_bound=self.lower_bound,
            upper_bound=self.upper_bound,
            processes=self.processes)


class MassConfiguration(with_metaclass(Singleton, MassBaseConfiguration)):
    """Define the configuration to be singleton based."""
