# -*- coding: utf-8 -*-
"""Module implementing concentration sampling for :mod:`mass` models.

Based on sampling implementations in :mod:`cobra.sampling.sampling`

"""
import pandas as pd

from mass.thermo.conc_sampling.conc_achr import ConcACHRSampler
from mass.thermo.conc_sampling.conc_optgp import ConcOptGPSampler


def sample_concentrations(
    concentration_solver, n, method="optgp", thinning=100, processes=1, seed=None
):
    """Sample valid concentration distributions from a :mod:`mass` model.

    This function samples valid concentration distributions from a
    :mod:`mass` model using a :class:`.ConcSolver`.

    Currently supports two methods.

    1. ``'optgp'`` which uses the :class:`~.ConcOptGPSampler` that supports
       parallel sampling :cite:`MHM14`. Requires large numbers of samples to be
       performant (n > 1000). For smaller samples, ``'achr'`` might be better
       suited.
    2. ``'achr'`` which uses artificial centering hit-and-run via the
       :class:`.ConcACHRSampler`. This is a single process method with
       good convergence :cite:`KS98`.

    Parameters
    ----------
    concentration_solver : ConcSolver
        The :class:`.ConcSolver` to use in generating samples.
    n : int
        The number of samples to obtain. When using ``'method=optgp'``, this
        must be a multiple of ``processes``, otherwise a larger number of
        samples will be returned.
    method : str
        The sampling algorithm to use. Default is ``'optgp'``.
    thinning : int
        The thinning factor for the generated sampling chain as a positive
        ``int`` > 0. A thinning factor of 10 means samples are returned every
        10 steps. If set to one, all iterates are returned.

        Default is ``100``.
    processes : int or None
        The number of processes used to generate samples. If ``None`` the
        number of processes specified in the :class:`.MassConfiguration` is
        utilized. Only valid for ``method='optgp'``.

        Default is ``1``.
    seed : int or None
        A positive ``int`` > 0 indiciating random number seed that should be
        used. If ``None`` provided, the current time stamp is used.

        Default is ``None``.

    Returns
    -------
    pandas.DataFrame
        The generated concentration samples. Each row corresponds to a sample
        of the concentrations and the columns are the metabolites.

    See Also
    --------
    :meth:`.ConcSolver.setup_sampling_problem`
        For setup of the sampling problem in the given :class:`.ConcSolver`.

    """
    if method == "optgp":
        sampler = ConcOptGPSampler(
            concentration_solver, processes=processes, thinning=thinning, seed=seed
        )
    elif method == "achr":
        sampler = ConcACHRSampler(concentration_solver, thinning=thinning, seed=seed)
    else:
        raise ValueError("method must be either 'optgp' or 'achr'.")

    return pd.DataFrame(
        columns=concentration_solver.included_metabolites, data=sampler.sample(n)
    )


__all__ = ("sample_concentrations",)
