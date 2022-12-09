# -*- coding: utf-8 -*-

"""Provide unified interfaces for optimization solutions for concentrations.

Based on solution implementations in :mod:`cobra.core.solution`

"""
import pandas as pd
from cobra.util.solver import check_solver_status
from numpy import array, exp, nan
from optlang.interface import OPTIMAL

from mass.core.mass_configuration import MassConfiguration
from mass.util.util import (
    _check_kwargs,
    apply_decimal_precision,
    get_public_attributes_and_methods,
)


MASSCONFIGURATION = MassConfiguration()


class ConcSolution:
    """A unified interface to a :class:`.ConcSolver` optimization solution.

    Notes
    -----
    The :class:`.ConcSolution` is meant to be constructed by
    :func:`get_concentration_solution` please look at that function to fully
    understand the :class:`ConcSolution` class.

    Attributes
    ----------
    objective_value : float
        The (optimal) value for the objective function.
    status : str
        The solver status related to the solution.
    concentrations : pandas.pd.Series
        Contains the metabolite concentrations which are the primal values
        of metabolite variables.
    concentration_reduced_costs : pandas.pd.Series
        Contains metabolite reduced costs, which are the dual values of
        metabolites variables.
    Keqs : pandas.pd.Series
        Contains the reaction equilibrium constant values, which are primal
        values of Keq variables.
    Keq_reduced_costs : pandas.pd.Series
        Contains reaction equilibrium constant reduced costs, which are the
        dual values of Keq variables.
    shadow_prices : pandas.pd.Series
        Contains reaction shadow prices (dual values of constraints).

    """

    def __init__(
        self,
        objective_value,
        status,
        concentrations,
        Keqs,
        concentration_reduced_costs=None,
        Keq_reduced_costs=None,
        shadow_prices=None,
    ):
        """Initialize the ConcSolution."""
        super(ConcSolution, self).__init__()
        # For solver objective value and status
        self.objective_value = objective_value
        self.status = status

        # For solver variables
        self.concentrations = concentrations
        self.Keqs = Keqs
        # For variable reduced costs
        self.concentration_reduced_costs = concentration_reduced_costs
        self.Keq_reduced_costs = Keq_reduced_costs
        # For constraint shadow prices
        self.shadow_prices = shadow_prices

    def concentrations_to_frame(self):
        """Get a :class:`pandas.pd.DataFrame` of concs. and reduced costs."""
        return pd.DataFrame(
            {
                "concentrations": self.concentrations,
                "reduced_costs": self.concentration_reduced_costs,
            }
        )

    def Keqs_to_frame(self):
        """Get a :class:`pandas.pd.DataFrame` of Keqs and reduced costs."""
        return pd.DataFrame(
            {"Keqs": self.Keqs, "reduced_costs": self.Keq_reduced_costs}
        )

    def to_frame(self):
        """Get a :class:`pandas.pd.DataFrame` of variables and reduced costs."""
        return pd.DataFrame(
            {
                "variables": pd.concat((self.concentrations, self.Keqs)),
                "reduced_costs": pd.concat(
                    (self.concentration_reduced_costs, self.Keq_reduced_costs)
                ),
            }
        )

    def _repr_html_(self):
        """HTML representation of the overview for the ConcSolution.

        Warnings
        --------
        This method is intended for internal use only.

        """
        if self.status == OPTIMAL:
            with pd.option_context("display.max_rows", 10):
                html = (
                    "<strong><em>Optimal</em> solution with objective "
                    "value {:.3f}</strong><br>{}".format(
                        self.objective_value, self.to_frame()._repr_html_()
                    )
                )
        else:
            html = "<strong><em>{}</em> solution</strong>".format(self.status)
        return html

    def __repr__(self):
        """Set string representation of the solution instance.

        Warnings
        --------
        This method is intended for internal use only.

        """
        if self.status != OPTIMAL:
            return "<Solution {0:s} at 0x{1:x}>".format(self.status, id(self))
        return "<Solution {0:.3f} at 0x{1:x}>".format(self.objective_value, id(self))

    def __getitem__(self, variable):
        """Return the value of a metabolite concentration or reaction Keq.

        Parameters
        ----------
        variable : str
            A variable ID for a variable in the solution.

        Warnings
        --------
        This method is intended for internal use only.

        """
        try:
            return self.concentrations[str(variable)]
        except KeyError:
            pass

        try:
            return self.Keqs[str(variable)]
        except KeyError as e:
            raise ValueError(
                "{0!r} is not a str ID of a ConcSolution variable.".format(str(e))
            )

    def __dir__(self):
        """Override default dir() implementation to list only public items.

        Warnings
        --------
        This method is intended for internal use only.

        """
        return get_public_attributes_and_methods(self)

    get_primal_by_id = __getitem__


def get_concentration_solution(
    concentration_solver, metabolites=None, reactions=None, raise_error=False, **kwargs
):
    """Generate a solution representation of a :class:`.ConcSolver` state.

    Parameters
    ---------
    concentration_solver : ConcSolver
        The :class:`.ConcSolver` containing the mathematical problem solved.
    metabolites : list
        An iterable of :class:`.MassMetabolite` objects.
        Uses :attr:`.ConcSolver.included_metabolites` by default.
    reactions : list
        An iterable of :class:`.MassReaction` objects.
        Uses :attr:`.ConcSolver.included_reactions` by default.
    raise_error : bool
        Whether to raise an OptimizationError if solver status is not optimal.
    **kwargs
        decimal_precision :
            ``bool`` indicating whether to apply the
            :attr:`~.MassBaseConfiguration.decimal_precision` attribute of
            the :class:`.MassConfiguration` to the solution values.

            Default is ``False``.

    Returns
    -------
    ConcSolution
        The solution of the optimization as a :class:`ConcSolution` object.

    """
    kwargs = _check_kwargs(
        {
            "decimal_precision": False,
        },
        kwargs,
    )

    check_solver_status(concentration_solver.solver.status, raise_error=raise_error)

    # Get included metabolites and reactions
    metabolites = concentration_solver._get_included_metabolites(metabolites)
    reactions = concentration_solver._get_included_reactions(reactions)

    # Get variable IDs, metabolites, and reactions for Keqs and constraints
    metabolites = [m.id for m in metabolites if m.id in concentration_solver.variables]
    Keq_ids = [
        r.Keq_str for r in reactions if r.Keq_str in concentration_solver.variables
    ]
    reactions = [r.id for r in reactions if r.id in concentration_solver.constraints]

    # Get metabolite and Keq primal values
    concs = array([concentration_solver.solver.primal_values[m] for m in metabolites])
    Keqs = array([concentration_solver.solver.primal_values[Keq] for Keq in Keq_ids])

    if concentration_solver.solver.is_integer:
        # Fill irrelevant arrays with nan
        reduced_concs = array([nan] * len(metabolites))
        reduced_Keqs = array([nan] * len(Keq_ids))
        shadow = array([nan] * len(reactions))
    else:
        # Get reduced cost values and shadow prices
        reduced_concs = array(
            [concentration_solver.solver.reduced_costs[m] for m in metabolites]
        )
        reduced_Keqs = array(
            [concentration_solver.solver.reduced_costs[Keq] for Keq in Keq_ids]
        )
        shadow = array(
            [concentration_solver.solver.shadow_prices[r] for r in reactions]
        )

    def transform_values(arr, **kwargs):
        """Transform array from logs to linear space and round if desired."""
        if kwargs.get("decimal_precision"):
            arr = apply_decimal_precision(arr, MASSCONFIGURATION.decimal_precision)
        return arr

    objective_value = transform_values(
        exp(concentration_solver.solver.objective.value), **kwargs
    )
    concs = transform_values(exp(concs), **kwargs)
    Keqs = transform_values(exp(Keqs), **kwargs)
    reduced_concs = transform_values(reduced_concs, **kwargs)
    reduced_Keqs = transform_values(reduced_Keqs, **kwargs)
    shadow = transform_values(shadow, **kwargs)

    return ConcSolution(
        objective_value,
        concentration_solver.solver.status,
        pd.Series(concs, metabolites, name="concentrations"),
        pd.Series(Keqs, Keq_ids, name="Keqs"),
        pd.Series(reduced_concs, metabolites, name="concentration_reduced_costs"),
        pd.Series(reduced_Keqs, Keq_ids, name="Keq_reduced_costs"),
        pd.Series(shadow, reactions, name="shadow_prices"),
    )


def update_model_with_concentration_solution(
    model, concentration_solution, concentrations=True, Keqs=True, inplace=True
):
    """Update a :mod:`mass` model with values from a :class:`ConcSolution`.

    Parameters
    ----------
    model : MassModel
        A :mod:`mass` model to update with the new solution values.
    concentration_solution : ConcSolution
        The :class:`ConcSolution` containing the solution values to use.
    concentrations : bool
        Whether to update the metabolite concentrations of the model
        (the :attr:`.MassMetabolite.initial_condition` values).
    Keqs : bool
        Whether to update the reaction equilibrium constants of the model
        (the :attr:`.MassReaction.equilibrium_constant` values).
    inplace : bool
        Whether to modify the given model or to modify a copy of the model.

    Returns
    -------
    MassModel
        Either the given model if ``inplace=True``, or a new copy of the model
        ``inplace=False``.

    """
    if not isinstance(concentration_solution, ConcSolution):
        raise TypeError("Must be a ConcSolution object.")

    if not inplace:
        model = model.copy()

    if concentrations:
        model.update_initial_conditions(
            concentration_solution.concentrations.to_dict(), verbose=False
        )
    if Keqs:
        model.update_parameters(concentration_solution.Keqs.to_dict(), verbose=False)

    return model


__all__ = (
    "ConcSolution",
    "get_concentration_solution",
    "update_model_with_concentration_solution",
)
