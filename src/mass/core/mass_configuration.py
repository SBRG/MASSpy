# -*- coding: utf-8 -*-
"""
Define the global configuration values through the :class:`MassConfiguration`.

Involved in model construction:
    * :meth:`~MassBaseConfiguration.boundary_compartment`
    * :meth:`~MassBaseConfiguration.default_compartment`
    * :meth:`~MassBaseConfiguration.irreversible_Keq`
    * :meth:`~MassBaseConfiguration.irreversible_kr`
    * :meth:`~MassBaseConfiguration.exclude_metabolites_from_rates`
    * :meth:`~MassBaseConfiguration.exclude_compartment_volumes_in_rates`
    * :meth:`~MassBaseConfiguration.model_creator`
    * :meth:`~MassBaseConfiguration.decimal_precision`

Involved in model simulation:
    * :meth:`~MassBaseConfiguration.decimal_precision`
    * :meth:`~MassBaseConfiguration.steady_state_threshold`

Involved in flux balance analysis (FBA):
    * :meth:`~MassBaseConfiguration.solver`
    * :meth:`~MassBaseConfiguration.tolerance`
    * :meth:`~MassBaseConfiguration.processes`
    * :meth:`~MassBaseConfiguration.lower_bound`
    * :meth:`~MassBaseConfiguration.upper_bound`
    * :meth:`~MassBaseConfiguration.bounds`

Involved in thermodynamics:
    * :meth:`~MassBaseConfiguration.decimal_precision`
    * :meth:`~MassBaseConfiguration.solver`
    * :meth:`~MassBaseConfiguration.tolerance`
    * :meth:`~MassBaseConfiguration.processes`

Notes
-----
* Some attributes are used in multiple :mod:`mass` submodules, such as
  :attr:`~MassBaseConfiguration.decimal_precision`.

* The :class:`MassConfiguration` shares the configuration attributes of the
  :class:`~.Configuration` class from :mod:`cobra`.

"""
from cobra.core.configuration import Configuration
from cobra.core.singleton import Singleton
from cobra.util.solver import interface_to_str
from six import integer_types, iteritems, itervalues, string_types


COBRA_CONFIGURATION = Configuration()


class MassBaseConfiguration:
    """Define global configuration values honored by :mod:`mass` functions.

    Notes
    -----
    The :class:`MassConfiguration` should be always be used over the
    :class:`MassBaseConfiguration` in order for global configuration to work
    as intended.

    """

    def __init__(self):
        """Initialize MassBaseConfiguration."""
        # Model construction configuration options
        self._boundary_compartment = {"b": "boundary"}
        self._default_compartment = {"compartment": "default_compartment"}
        self._irreversible_Keq = float("inf")
        self._irreversible_kr = None
        self._exclude_metabolites_from_rates = {
            "elements": [{"H": 2, "O": 1}, {"H": 1}]
        }
        self._exclude_compartment_volumes_in_rates = True
        self._model_creator = {
            "familyName": "",
            "givenName": "",
            "organization": "",
            "email": "",
        }

        # Model simulation options
        self._decimal_precision = None
        self._steady_state_threshold = 1e-6

        # For cobra configuration synchronization
        self._shared_state = COBRA_CONFIGURATION.__dict__

    @property
    def boundary_compartment(self):
        """Get or set the default value for the boundary compartment.

        Parameters
        ----------
        compartment_dict : dict
            A ``dict`` containing the identifier of the boundary compartment
            mapped to the name of the boundary compartment.

        """
        return getattr(self, "_boundary_compartment")

    @boundary_compartment.setter
    def boundary_compartment(self, compartment_dict):
        """Set the default value for the boundary compartment."""
        setattr(self, "_boundary_compartment", compartment_dict)

    @property
    def default_compartment(self):
        """Get or set the default value for the default compartment.

        Parameters
        ----------
        compartment_dict : dict
            A ``dict`` containing the identifier of the default compartment
            mapped to the name of the default compartment.

        """
        return getattr(self, "_default_compartment")

    @default_compartment.setter
    def default_compartment(self, compartment_dict):
        """Set the default value for the default compartment."""
        setattr(self, "_default_compartment", compartment_dict)

    @property
    def irreversible_Keq(self):
        r"""Get or set the default ``Keq`` value of an irreversible reaction.

        Notes
        -----
        :attr:`.MassReaction.equilibrium_constant`\ s cannot be negative.

        Parameters
        ----------
        value : float
            A non-negative number for the equilibrium constant (Keq)
            of the reaction.

        Raises
        ------
        ValueError
            Occurs when trying to set a negative value.

        """
        return getattr(self, "_irreversible_Keq")

    @irreversible_Keq.setter
    def irreversible_Keq(self, value):
        """Set the default value for Keq of an irreversible reaction."""
        if value is not None:
            if not isinstance(value, (integer_types, float)):
                raise TypeError("Must be an int or float")
            if value < 0.0:
                raise ValueError("Must be a non-negative number")
        setattr(self, "_irreversible_Keq", value)

    @property
    def irreversible_kr(self):
        r"""Get or set the default ``kr`` value of an irreversible reaction.

        Notes
        -----
        :attr:`.MassReaction.reverse_rate_constant`\ s cannot be negative.

        Parameters
        ----------
        value : float
            A non-negative number for the reverse rate constant (kr)
            of the reaction.

        Raises
        ------
        ValueError
            Occurs when trying to set a negative value.

        """
        return getattr(self, "_irreversible_kr")

    @irreversible_kr.setter
    def irreversible_kr(self, value):
        """Set the default value for kr of an irreversible reaction."""
        if value is not None:
            if not isinstance(value, (integer_types, float)):
                raise TypeError("Must be an int or float")
            if value < 0.0:
                raise ValueError("Must be a non-negative number")
        setattr(self, "_irreversible_kr", value)

    @property
    def model_creator(self):
        """Get or set values for the ``dict`` representing the model creator.

        Notes
        -----
        * A read-only copy of the ``dict`` is returned.
        * To successfully export a model creator, all keys must have non-empty
          string values.

        Parameters
        ----------
        creator_dict : dict
            A ``dict`` containing the model creator information where keys
            are SBML creator fields and values are strings. Keys can only
            be the following:

            * 'familyName'
            * 'givenName'
            * 'organization'
            * 'email'

            Values must be strings or ``None``.

        """
        return self._model_creator.copy()

    @model_creator.setter
    def model_creator(self, creator_dict):
        """Set the information in the dict representing the model creator."""
        valid = {"familyName", "givenName", "organization", "email"}
        for k, v in iteritems(creator_dict):
            if k not in valid:
                raise ValueError(
                    "Invalid key '{0}'. Keys can only be the"
                    " following: {1:r}".format(k, str(valid))
                )
            if v is not None and not isinstance(v, string_types):
                raise TypeError(
                    "'{0}' not a string. Values must be strings or"
                    " None.".format(str(v))
                )

        self._model_creator.update(creator_dict)

    @property
    def exclude_metabolites_from_rates(self):
        """Get or set the metabolites that should be excluded from rates.

        Default is ``dict("elements", [{"H": 2, "O": 1}, {"H": 1}])`` to
        remove the hydrogen and water metabolites using the
        :attr:`~.MassMetabolite.elements` attribute to filter out the hydrogen
        and water in all rates except the hydrogen and water exchange
        reactions.

        Parameters
        ----------
        to_exclude : dict
            A dict where keys should correspond to a metabolite attribute to
            utilize for filtering, and values are lists that contain the items
            to exclude that would be returned by the metabolite attribute.
            Does not apply to boundary reactions.

        """
        return getattr(self, "_exclude_metabolites_from_rates")

    @exclude_metabolites_from_rates.setter
    def exclude_metabolites_from_rates(self, to_exclude):
        """Set the metabolites that should be excluded from rates."""
        if not isinstance(to_exclude, dict):
            raise TypeError("Must be a ``dict``.")
        setattr(self, "_exclude_metabolites_from_rates", to_exclude)

    @property
    def exclude_compartment_volumes_in_rates(self):
        """Get or set whether to exclude the compartment volumes in rates.

        The boundary compartment will always excluded.

        Parameters
        ----------
        value : bool
            Whether to exclude the compartment volumes in rate expressions.

        """
        return getattr(self, "_exclude_compartment_volumes_in_rates")

    @exclude_compartment_volumes_in_rates.setter
    def exclude_compartment_volumes_in_rates(self, value):
        """Set whether to exclude the compartment volumes in rates."""
        if not isinstance(value, bool):
            raise TypeError("Must be a ``bool``.")
        setattr(self, "_exclude_compartment_volumes_in_rates", value)

    @property
    def decimal_precision(self):
        """Get or set the default decimal precision when rounding.

        Positive numbers indicated digits to the right of the
        decimal place, negative numbers indicate digits to the left of the
        decimal place.

        Notes
        -----
        The :attr:`decimal_precison` is applied as follows::

            new_value = round(value, decimal_precison)

        Parameters
        ----------
        precision : int or None
            An integer indicating how many digits from the decimal should
            rounding occur. If ``None``, no rounding will occur.

        """
        return getattr(self, "_decimal_precision")

    @decimal_precision.setter
    def decimal_precision(self, precision):
        """Set the default decimal precision when rounding."""
        if precision is not None and not isinstance(precision, integer_types):
            raise TypeError("precision must be an int.")

        setattr(self, "_decimal_precision", precision)

    @property
    def steady_state_threshold(self):
        """Get or set the steady state threshold when using roadrunner solvers.

        A threshold for determining whether the RoadRunner steady state solver
        is at steady state. The steady state solver returns a value indicating
        how close the solution is to the steady state, where smaller values
        are better. Values less than the threshold indicate steady state.

        Notes
        -----
        * With simulations. the absolute difference between the last two points
          must be less than the steady state threshold.
        * With steady state solvers, the sum of squares of the steady state
          solution must be less than the steady state threshold.
        * Steady state threshold values cannot be negative.

        Parameters
        ----------
        threshold : float
            The threshold for determining whether a steady state occurred.

        Raises
        ------
        ValueError
            Occurs when trying to set a negative value.

        """
        return getattr(self, "_steady_state_threshold")

    @steady_state_threshold.setter
    def steady_state_threshold(self, threshold):
        """Set the default decimal precision when rounding."""
        if not isinstance(threshold, (integer_types, float)):
            raise TypeError("Must be an int or float")
        if threshold < 0.0:
            raise ValueError("Must be a non-negative number")
        setattr(self, "_steady_state_threshold", threshold)

    @property
    def solver(self):
        """Get or set the solver utilized for optimization.

        The solver choices are the ones
        provided by `optlang` and solvers installed in the environment.

        Parameters
        ----------
        solver : str
            The solver to utilize in optimizations. Valid solvers typically
            include:

                * ``"glpk"``
                * ``"cplex"``
                * ``"gurobi"``

        Raises
        ------
        :class:`cobra.exceptions.SolverNotFound`
            Occurs for invalid solver values.

        """
        return COBRA_CONFIGURATION.solver

    @solver.setter
    def solver(self, solver):
        """Set the solver utilized for optimization."""
        COBRA_CONFIGURATION.solver = solver

    @property
    def tolerance(self):
        """Get or set the tolerance value utilized by the optimization solver.

        Parameters
        ----------
        tol : float
            The tolerance value to set.

        """
        return COBRA_CONFIGURATION.tolerance

    @tolerance.setter
    def tolerance(self, tol):
        """Set the tolerance value utilized by the optimization solver."""
        COBRA_CONFIGURATION.tolerance = tol

    @property
    def lower_bound(self):
        """Get or set the default value of the lower bound for reactions.

        Parameters
        ----------
        bound : float
            The default bound value to set.

        """
        return COBRA_CONFIGURATION.lower_bound

    @lower_bound.setter
    def lower_bound(self, bound):
        """Set the default value of the lower bound for reactions."""
        COBRA_CONFIGURATION.lower_bound = bound

    @property
    def upper_bound(self):
        """Get or set the default value of the lower bound for reactions.

        Parameters
        ----------
        bound : float
            The default bound value to set.

        """
        return COBRA_CONFIGURATION.upper_bound

    @upper_bound.setter
    def upper_bound(self, bound):
        """Set the default value of the lower bound for reactions."""
        COBRA_CONFIGURATION.upper_bound = bound

    @property
    def bounds(self):
        """Get or set the default lower and upper bounds for reactions.

        Parameters
        ----------
        bounds : tuple of floats
            A tuple of floats to set as the new default bounds in the form
            of ``(lower_bound, upper_bound)``.

        Raises
        ------
        AssertionError
            Occurs when lower bound is greater than the upper bound.

        """
        return COBRA_CONFIGURATION.bounds

    @bounds.setter
    def bounds(self, bounds):
        """Set the default lower and upper bounds for reactions."""
        COBRA_CONFIGURATION.bounds = bounds

    @property
    def processes(self):
        """Return the default number of processes to use when possible.

        Number of processes to use where multiprocessing is
        possible. The default number corresponds to the number of available
        cores (hyperthreads).
        """
        return COBRA_CONFIGURATION.processes

    @property
    def shared_state(self):
        """Return a read-only ``dict`` for shared configuration attributes."""
        return {k.lstrip("_"): v for k, v in iteritems(getattr(self, "_shared_state"))}

    def _repr_html_(self):
        """Return the HTML representation of the MassConfiguration.

        Warnings
        --------
        This method is intended for internal use only.

        """
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
            </tr><tr>
                <td><strong>Metabolites excluded in rates</strong></td>
                <td>{excluded_metabolites_in_rates}</td>
            </tr><tr>
                <td><strong>Compartments in rates</strong></td>
                <td>{exclude_comp_vols}</td>
            </tr><tr>
                <td><strong>Model creator set</strong></td>
                <td>{model_creator}</td>
            </tr><tr>
                <td><strong>Decimal precision</strong></td>
                <td>{decimal_precision}</td>
            </tr><tr>
                <td><strong>Steady state threshold</strong></td>
                <td>{steady_state_threshold}</td>
            </tr>
                <td><strong>Solver</strong></td>
                <td>{solver}</td>
            </tr><tr>
                <td><strong>Solver tolerance</strong></td>
                <td>{tolerance}</td>
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
                "{0} ({1})".format(v, k) if v else k
                for k, v in iteritems(self.boundary_compartment)
            ][0],
            default_compartment=[
                "{0} ({1})".format(v, k) if v else k
                for k, v in iteritems(self.default_compartment)
            ][0],
            irreversible_Keq=self.irreversible_Keq,
            irreversible_kr=self.irreversible_kr,
            excluded_metabolites_in_rates=bool(self.exclude_metabolites_from_rates),
            exclude_comp_vols=self.exclude_compartment_volumes_in_rates,
            model_creator=bool(any(itervalues(self.model_creator))),
            decimal_precision=self.decimal_precision,
            steady_state_threshold=self.steady_state_threshold,
            solver=interface_to_str(self.solver),
            tolerance=self.tolerance,
            lower_bound=self.lower_bound,
            upper_bound=self.upper_bound,
            processes=self.processes,
        )

    def __repr__(self):
        """Override default :func:`repr` for the MassConfiguration.

        Warnings
        --------
        This method is intended for internal use only.

        """
        return """MassConfiguration:
        boundary compartment: {boundary_compartment}
        default compartment: {default_compartment}
        irreversible reaction Keq: {irreversible_Keq}
        irreversible reaction kr: {irreversible_kr}
        metabolites excluded in rates: {excluded_metabolites_in_rates}
        include compartments in rates: {exclude_comp_vols}
        model creator set: {model_creator}
        decimal_precision: {decimal_precision}
        steady_state_threshold: {steady_state_threshold}
        solver: {solver}
        solver tolerance: {tolerance}
        lower_bound: {lower_bound}
        upper_bound: {upper_bound}
        processes: {processes}""".format(
            boundary_compartment=[
                "{0} ({1})".format(v, k) if v else k
                for k, v in iteritems(self.boundary_compartment)
            ][0],
            default_compartment=[
                "{0} ({1})".format(v, k) if v else k
                for k, v in iteritems(self.default_compartment)
            ][0],
            irreversible_Keq=self.irreversible_Keq,
            irreversible_kr=self.irreversible_kr,
            excluded_metabolites_in_rates=bool(self.exclude_metabolites_from_rates),
            exclude_comp_vols=self.exclude_compartment_volumes_in_rates,
            model_creator=bool(any(itervalues(self.model_creator))),
            decimal_precision=self.decimal_precision,
            steady_state_threshold=self.steady_state_threshold,
            solver=interface_to_str(self.solver),
            tolerance=self.tolerance,
            lower_bound=self.lower_bound,
            upper_bound=self.upper_bound,
            processes=self.processes,
        )


class MassConfiguration(MassBaseConfiguration, metaclass=Singleton):
    """Define the configuration to be :class:`.Singleton` based."""


__all__ = (
    "MassConfiguration",
    "MassBaseConfiguration",
)
