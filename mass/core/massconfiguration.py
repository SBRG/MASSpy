# -*- coding: utf-8 -*-
"""TODO Module Docstrings."""
from __future__ import absolute_import

import logging

from six import with_metaclass

from cobra.core.configuration import BaseConfiguration
from cobra.core.singleton import Singleton
from cobra.util.solver import interface_to_str

from mass.util.util import ensure_non_negative_value

__all__ = ("MassConfiguration",)

LOGGER = logging.getLogger(__name__)


class MassBaseConfiguration(BaseConfiguration):
    """Define global configuration values that to be honored by mass functions.

    Attributes
    ----------
    irreversible_Keq: float
        The default value to assign to equilibrium constants (Keq) for
        irreversible reactions. Must be a non-negative value.
        Default value is the infinity (=float("inf")).
    irreversible_kr: float
        The default value to assign to equilibrium constants (Keq) for
        irreversible reactions. Must be a non-negative value.
        Default value is the 0.
    boundary_compartment: dict
        A dictionary containing the identifier of the boundary compartment
        mapped to the name of the boundary compartment.
        Default value is {"b": "boundary"}.
    default_compartment: dict
        A dictionary containing the identifier of the default compartment
        mapped to the name of the desired name of default compartment. Used for
        writing models to SBML when there are no set compartments in the model.
        Default value is {"default": "default_compartment"}.
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
        super(MassBaseConfiguration, self).__init__()
        self._boundary_compartment = {"b": "boundary"}
        self._default_compartment = {"default": "default_compartment"}
        self._irreversible_Keq = float("inf")
        self._irreversible_kr = 0

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
        ensure_non_negative_value(value)
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
        ensure_non_negative_value(value)
        setattr(self, "irreversible_kr", value)

    @property
    def optimization_solver(self):
        """Return the solver utilized for optimization."""
        return self.solver

    @optimization_solver.setter
    def optimization_solver(self, value):
        """Set the solver utilized for optimization."""
        self.solver = value

    @property
    def optimization_tolerance(self):
        """Return the tolerance value utilized by the optimization solver."""
        return self.tolerance

    @optimization_tolerance.setter
    def optimization_tolerance(self, value):
        """Set the tolerance value utilized by the optimization solver."""
        self.tolerance = value

    def __repr__(self):
        """Return the representation of the MassConfiguration."""
        return """
        optimization solver: {optimization_solver}
        optimization solver tolerance: {optimization_tolerance}
        lower_bound: {lower_bound}
        upper_bound: {upper_bound}
        processes: {processes}""".format(
            optimization_solver=interface_to_str(self.optimization_solver),
            optimization_tolerance=self.optimization_tolerance,
            lower_bound=self.lower_bound,
            upper_bound=self.upper_bound,
            processes=self.processes)

    def _repr_html_(self):
        return """
        <table>
            <tr>
                <td><strong>Optimization solver</strong></td>
                <td>{optimization_solver}</td>
            </tr>
            <tr>
                <td><strong>Optimization solver tolerance</strong></td>
                <td>{optimization_tolerance}</td>
            </tr>
            <tr>
                <td><strong>Lower bound</strong></td>
                <td>{lower_bound}</td>
            </tr>
            <tr>
                <td><strong>Upper bound</strong></td>
                <td>{upper_bound}</td>
            </tr>
            <tr>
                <td><strong>Processes</strong></td>
                <td>{processes}</td>
            </tr>
        </table>""".format(
            optimization_solver=interface_to_str(self.optimization_solver),
            optimization_tolerance=self.optimization_tolerance,
            lower_bound=self.lower_bound,
            upper_bound=self.upper_bound,
            processes=self.processes)


class MassConfiguration(with_metaclass(Singleton, MassBaseConfiguration)):
    """Define the configuration to be singleton based."""

    pass
