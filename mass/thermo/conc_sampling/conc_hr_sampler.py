# -*- coding: utf-8 -*-
"""Provide base class for Hit-and-Run concentration samplers.

New samplers should derive from the abstract :class:`ConcHRSampler` class
where possible to provide a uniform interface.""

Based on sampling implementations in :mod:`cobra.sampling.hr_sampler`.

"""
from copy import deepcopy
from time import time

from cobra.exceptions import OptimizationError
from cobra.sampling.hr_sampler import Problem, shared_np_array

import numpy as np

from optlang.interface import OPTIMAL

from six import iterkeys

from mass.core.mass_configuration import MassConfiguration
from mass.thermo.conc_solver import (
    ConcSolver, concentration_constraint_matricies)
from mass.util.matrix import nullspace
from mass.util.util import _make_logger, ensure_non_negative_value

LOGGER = _make_logger(__name__)
"""logging.Logger: Logger for :mod:`~.conc_hr_sampler` submodule."""

MAX_TRIES = 100
"""int: Maximum number of retries for sampling."""

MASSCONFIGURATION = MassConfiguration()


class ConcHRSampler:
    """The abstract base class for hit and run concentration samplers.

    Parameters
    ----------
    concentration_solver : ConcSolver
        The :class:`.ConcSolver` to use in generating samples.
    thinning : int
        The thinning factor for the generated sampling chain as a positive
        ``int`` > 0. A thinning factor of 10 means samples are returned every
        10 steps.
    nproj : int or None
        A positive ``int`` > 0 indicating how often to reporject the sampling
        point into the feasibility space. Avoids numerical issues at the cost
        of lower samplimg. If ``None`` then the value is determined via the
        following::

            nproj = int(min(len(self.concentration_solver.variables)**3, 1e6))

        Default is ``None``
    seed : int or None
        A positive ``int`` > 0 indiciating random number seed that should be
        used. If ``None`` provided, the current time stamp is used.

        Default is ``None``.

    Attributes
    ----------
    concentration_solver : ConcSolver
        The :class:`.ConcSolver` used to generate samples.
    feasibility_tol : float
        The tolerance used for checking equalities feasibility.
    bounds_tol : float
        The tolerance used for checking bounds feasibility.
    thinning : int
        The currently used thinning factor.
    n_samples : int
        The total number of samples that have been generated by this
        sampler instance.
    retries : int
        The overall of sampling retries the sampler has observed. Larger
        values indicate numerical instabilities.
    problem : collections.namedtuple
        A :class:`~collections.namedtuple` whose attributes define the entire
        sampling problem in matrix form. See docstring of
        :class:`~cobra.sampling.hr_sampler.Problem` for more information.
    warmup : numpy.matrix
        A matrix of with as many columns as reactions in the model and more
        than 3 rows containing a warmup sample in each row. None if no warmup
        points have been generated yet.
    nproj : int
        How often to reproject the sampling point into the feasibility space.
    met_var_idx : numpy.ndarray
        Has one entry for each metabolite in the model containing the index of
        the respective metabolite variable.

        Does not included :attr:`.ConcSolver.excluded_metabolites`.
    Keq_var_idx : numpy.ndarray
        Has one entry for each reaction in the model containing the index of
        the respective reaction Keq variable.

        Does not included :attr:`.ConcSolver.excluded_reactions`.

    """

    # pylint: disable=too-many-instance-attributes
    def __init__(self, concentration_solver, thinning, nproj=None, seed=None):
        """Initialize a new sampler."""
        if not isinstance(concentration_solver, ConcSolver):
            raise TypeError(
                "`concentration_solver` must be a ConcSolver instance")
        if concentration_solver.solver_problem_type != "sampling":
            raise ValueError(
                "The given `concentration_solver` is not set up for sampling. "
                "To set up the concentration solver for sampling, utilize the "
                "ConcSolver.setup_sampling_problem.")
        # Set a copy of the solver to prevent changes to the original
        self.concentration_solver = deepcopy(concentration_solver)

        # Set tolerances
        self.feasibility_tol = concentration_solver.tolerance
        self.bounds_tol = concentration_solver.tolerance

        # Set thinning factor
        self.thinning = thinning

        # Set nproj
        self._nproj = None
        self._nproj = nproj

        # Set n_samples, retries, and solver problem
        self.n_samples = 0
        self.retries = 0
        self.problem = self.__build_problem()

        # Set variable indicies
        var_idx = {
            v: idx for idx, v in enumerate(
                iterkeys(concentration_solver.variables))}
        self.met_var_idx = np.array([
            var_idx[m.id]
            for m in concentration_solver._get_included_metabolites()])

        self.Keq_var_idx = np.array([
            var_idx[r.Keq_str]
            for r in concentration_solver._get_included_reactions()])

        # Set warmup variables
        self.warmup = None
        self.n_warmup = 0

        # Set seed
        self._seed = None
        self.seed = seed

    @property
    def nproj(self):
        """Get or set :attr:`~ConcHRSampler.nproj` value.

        Parameters
        ----------
        value : int or None
            A positive ``int`` > 0 indicating how often to reporject the
            sampling point into the feasibility space. Avoids numerical issues
            at the cost of lower sampling. If ``None`` then the value is
            determined via the following::

                nproj = int(min(len(self.concentration_solver.variables)**3, 1e6))

        """  # noqa: E501
        return getattr(self, "_nproj")

    @nproj.setter
    def nproj(self, value):
        """Set :attr:`~ConcHRSampler.nproj` value."""
        # Verify value is valid
        try:
            value = ensure_non_negative_value(value, exclude_zero=True)
        except (ValueError, TypeError) as e:
            raise e.__class__("'nproj' {0}.".format(str(e).lower()))

        # Set value
        if value is None:
            value = int(min(len(self.concentration_solver.variables)**3, 1e6))

        setattr(self, "_nproj", value)

    @property
    def seed(self):
        """Get or set :attr:`~ConcHRSampler.nproj` value.

        Parameters
        ----------
        value : int or None
            A positive ``int`` > 0 indiciating random number seed that should
            be used. If ``None`` provided, the current time stamp is used.

        """
        return getattr(self, "_seed")

    @seed.setter
    def seed(self, value):
        """Set :attr:`~ConcHRSampler.seed` value."""
        # Verify value is valid
        try:
            value = ensure_non_negative_value(value, exclude_zero=True)
        except (ValueError, TypeError) as e:
            raise e.__class__("'seed' {0}.".format(str(e).lower()))

        # Set value
        if value is None:
            value = int(time())

        # Avoid overflow
        value = value % np.iinfo(np.int32).max
        setattr(self, "_seed", value)

    def generate_cva_warmup(self):
        """Generate the warmup points for the sampler.

        Generates warmup points by setting each concentration as the sole
        objective and minimizing/maximizing it. Also caches the projection of
        the warmup points into the nullspace for non-homogenous problems.

        """
        self.n_warmup = 0
        c_solver = self.concentration_solver
        metabolites = c_solver._get_included_metabolites()

        self.warmup = np.zeros((2 * len(metabolites), len(c_solver.variables)))
        for sense in ("min", "max"):
            c_solver.objective_direction = sense

            for i, met in enumerate(metabolites):
                variable = c_solver.variables[self.met_var_idx[i]]

                # Omit fixed metabolites if they are non-homogeneous
                if variable.ub - variable.lb < self.bounds_tol:
                    LOGGER.info("Skipping fixed metabolite %s", met.id)
                    continue
                # Set the objective
                c_solver.objective.set_linear_coefficients({variable: 1})

                c_solver.slim_optimize()

                if not c_solver.solver.status == OPTIMAL:
                    LOGGER.info("Can not %simize metabolite %s, skipping it",
                                sense, met.id)
                    c_solver.objective.set_linear_coefficients({variable: 0})
                    continue

                primals = c_solver.solver.primal_values
                sol = [primals[v.name] for v in c_solver.variables]
                self.warmup[self.n_warmup, ] = sol
                self.n_warmup += 1

                # Reset objective
                c_solver.objective.set_linear_coefficients({variable: 0})

        # Shrink to measure
        self.warmup = self.warmup[0:self.n_warmup, :]
        if self.warmup.size == 0:
            raise OptimizationError(
                "CVA warmup found no feasible solutions. Ensure the systems "
                "has the appropriate variables and constraints by excluding "
                "certain metabolites (e.g. hydrogen) and reactions "
                "(e.g. boundary reactions), and by indicating the "
                "equilibrium reactions.")
        # Remove redundant search directions
        keep = np.logical_not(self._is_redundant(self.warmup))
        self.warmup = self.warmup[keep, :]
        self.n_warmup = self.warmup.shape[0]

        # Catch some special cases
        if len(self.warmup.shape) == 1 or self.warmup.shape[0] == 1:
            raise ValueError(
                "The concentration cone consists only of a single point.")
        if self.n_warmup == 2:
            if not self.problem.homogeneous:
                raise ValueError("Can not sample from an inhomogenous problem"
                                 " with only 2 search directions.")
            LOGGER.info("All search directions on a line, adding another one.")
            newdir = self.warmup.T.dot([0.25, 0.25])
            self.warmup = np.vstack([self.warmup, newdir])
            self.n_warmup += 1

        # Shrink warmup points to measure
        self.warmup = shared_np_array(
            (self.n_warmup, len(c_solver.variables)), self.warmup)

    def sample(self, n, concs=True):
        """Abstract sampling function.

        Should be overwritten by child classes.

        """

    def batch(self, batch_size, batch_num, concs=True):
        """Create a batch generator.

        This is useful to generate ``n`` batches of ``m`` samples each.

        Parameters
        ----------
        batch_size : int
            The number of samples contained in each batch (``m``).
        batch_num : int
            The number of batches in the generator (``n``).
        concs : boolean
            Whether to return concentrations or the internal solver variables.
            If ``False`` will return a variable for each metabolite and
            reaction Keq as well as all additional variables that may have
            been defined in the model.

        Yields
        ------
        pandas.core.frame.DataFrame
            A :class:`pandas.DataFrame <pandas.core.frame.DataFrame>` with
            dimensions ``(batch_size x n_m)`` containing a valid concentration
            sample for a total of n_m metabolites (or variables if
            ``concs=False``) in each row.

        """
        for _ in range(batch_num):
            yield self.sample(batch_size, concs=concs)

    def _reproject(self, p):
        """Reproject a point into the feasibility region.

        This function is guarunteed to return a new feasible point. However,
        no guaruntees in terms of proximity to the original point can be made.

        Parameters
        ----------
        p : numpy.ndarray
            The current sample point.

        Returns
        -------
        numpy.ndarray
            A new feasible point. If `p` was feasible it wil return p.

        Warnings
        --------
        This method is intended for internal use only.

        """
        nulls = self.problem.nullspace
        equalities = self.problem.equalities

        # Don't reproject if the point is feasible
        if np.allclose(equalities.dot(p), self.problem.b,
                       rtole=0, atol=self.feasibility_tol):
            new = p
        else:
            LOGGER.info(
                "feasibility violated in sample %d, trying to reproject",
                self.n_samples)
            new = nulls.dot(nulls.T.dot(p))

        # If projection may violate bounds, set to random point in space
        if any(new != p):
            LOGGER.info(
                "Reprojection failed in sample %d, using random point in "
                "space", self.n_samples)
            new = self._random_point()

        return new

    def _random_point(self):
        """Find an approximately random point in the concentration cone.

        Warnings
        --------
        This method is intended for internal use only.

        """
        idx = np.random.randint(self.n_warmup,
                                size=min(2, np.ceil(np.sqrt(self.n_warmup))))
        return self.warmup[idx, :].mean(axis=0)

    def _is_redundant(self, matrix, cutoff=None):
        """Identify redundant rows in a matrix that can be removed.

        Warnings
        --------
        This method is intended for internal use only.

        """
        cutoff = 1.0 - self.feasibility_tol

        # Avoid zero variances
        extra_col = matrix[:, 0] + 1

        # Avoid zero rows being correlated with constant rows
        extra_col[matrix.sum(axis=1) == 0] = 2
        corr = np.corrcoef(np.c_[matrix, extra_col])
        corr = np.tril(corr, -1)

        return (np.abs(corr) > cutoff).any(axis=1)

    def _bounds_dist(self, p):
        """Get the lower and upper bound distances. Negative is bad.

        Warnings
        --------
        This method is intended for internal use only.

        """
        problem = self.problem
        lb_dist = (p - problem.variable_bounds[0, ]).min()
        ub_dist = (problem.variable_bounds[1, ] - p).min()

        if problem.bounds.shape[0] > 0:
            const = problem.inequalities.dot(p)
            const_lb_dist = (const - problem.bounds[0, ]).min()
            const_ub_dist = (problem.bounds[1, ] - const).min()

            lb_dist = min(lb_dist, const_lb_dist)
            ub_dist = min(ub_dist, const_ub_dist)

        return np.array([lb_dist, ub_dist])

    def __build_problem(self):
        """Build the matrix representation of the sampling problem.

        Warnings
        --------
        This method is intended for internal use only.

        """
        problem = concentration_constraint_matricies(
            self.concentration_solver, zero_tol=self.feasibility_tol)

        equalities = problem.equalities
        b = problem.b
        inequalities = problem.inequalities
        bounds = np.atleast_2d(problem.bounds).T
        variable_fixed = problem.variable_fixed
        variable_bounds = np.atleast_2d(problem.variable_bounds).T

        fixed_non_zero = np.abs(problem.variable_bounds[:, 1]) > \
            self.feasibility_tol
        fixed_non_zero &= variable_fixed

        # Add any non-zero fixed variables as equalities to the matrix.
        if any(fixed_non_zero):
            n_fixed = fixed_non_zero.sum()
            rows = np.zeros((n_fixed, equalities.shape[1]))
            rows[range(n_fixed), np.where(fixed_non_zero)] = 1.0
            # Add to equalities
            equalities = np.vstack([equalities, rows])
            b = np.hstack([b, problem.variable_bounds[:, 1][fixed_non_zero]])

        # Set up a projection that can cast point into the nullspace
        nulls = nullspace(self.concentration_solver.model.S.T)
        # Create the problem
        return Problem(
            equalities=shared_np_array(equalities.shape, equalities),
            b=shared_np_array(b.shape, b),
            inequalities=shared_np_array(inequalities.shape, inequalities),
            bounds=shared_np_array(bounds.shape, bounds),
            variable_fixed=shared_np_array(variable_fixed.shape,
                                           variable_fixed, integer=True),
            variable_bounds=shared_np_array(variable_bounds.shape,
                                            variable_bounds),
            nullspace=shared_np_array(nulls.shape, nulls),
            homogeneous=False  # Conc sampling is never homogeneous
        )


# Required by ACHRSampler and OptGPSampler
# Has to be declared outside of class to be used for multiprocessing
def step(sampler, x, delta, fraction=None, tries=0):
    """Sample new feasible point from point ``x`` in the direction ``delta``.

    Has to be declared outside of class to be used for multiprocessing

    Parameters
    ----------
    sampler : ConcHRSampler
        The sampler instance.
    x : float
        The starting point from which to sample.
    delta : float
        The direction to travel from the point at ``x``.
    fraction : float or None
        The fraction of the alpha range to use in determining alpha.
        If ``None`` then the :func:`np.random.uniform` function to get alpha.
    tries : int
        Number of tries. If the number of tries is greater than the
        :const:`MAX_TRIES`, a :class:`RuntimeError` will be raised.

    Returns
    -------
    float
        The new feasible point.

    Raises
    ------
    RunTimeError
        Raised when ``tries > MAX_TRIES``

    """
    problem = sampler.problem
    valid = ((np.abs(delta) > sampler.feasibility_tol) &
             np.logical_not(problem.variable_fixed))

    # Permisible alphas for staying in variable bounds
    valphas = (
        (1.0 - sampler.bounds_tol) * problem.variable_bounds - x)[:, valid]
    valphas = (valphas / delta[valid]).flatten()

    if problem.bounds.shape[0] > 0:
        # Permissible alphas for staying in constraint bounds
        ineqs = problem.inequalities.dot(delta)
        valid = np.abs(ineqs) > sampler.feasibility_tol
        balphas = ((1.0 - sampler.bounds_tol) * problem.bounds -
                   problem.inequalities.dot(x))[:, valid]
        balphas = (balphas / ineqs[valid]).flatten()

        # Combined alphas
        alphas = np.hstack([valphas, balphas])
    else:
        alphas = valphas
    # Split positive and negative alphas and get alpha range
    pos_alphas = alphas[alphas > 0.0]
    neg_alphas = alphas[alphas <= 0.0]

    alpha_range = np.array([neg_alphas.max() if neg_alphas.size > 0 else 0,
                            pos_alphas.min() if pos_alphas.size > 0 else 0])

    if fraction:
        alpha = alpha_range[0] + fraction * (alpha_range[1] - alpha_range[0])
    else:
        alpha = np.random.uniform(alpha_range[0], alpha_range[1])

    p = x + alpha * delta

    # Numerical instabilities may cause bounds invalidation. If so,
    # reset sampler and sample from one of the original warmup directions
    # if that occurs. Also reset if it got stuck.
    if np.any(sampler._bounds_dist(p) < -sampler.bounds_tol) or\
       np.abs(np.abs(alpha_range).max() * delta).max() < sampler.bounds_tol:
        if tries > MAX_TRIES:
            raise RuntimeError("Can not escape sampling region, model seems "
                               "numerically unstable. Reporting the model to "
                               "https://github.com/SBRG/masspy/issues "
                               "will help us to fix this.")
        LOGGER.info(
            "found bounds infeasibility in sample, resetting to center")
        newdir = sampler.warmup[np.random.randint(sampler.n_warmup)]
        sampler.retries += 1

        return step(sampler, sampler.center, newdir - sampler.center, None,
                    tries + 1)
    return p


__all__ = ("ConcHRSampler", "LOGGER", "MAX_TRIES", "step")
