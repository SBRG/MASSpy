# -*- coding: utf-8 -*-
"""This module contains the Hit-and-Run concentration samplers.

Key components in the :mod:`.thermo.conc_sampling` module are the following:

    1. The :func:`~.conc_sampling.sample_concentrations` function, a function
       to call one of the concentration sampler methods and valid
       concentration distributions from the :associated :mod:`mass` model
       of the :class:`.ConcSolver`.
       Currently provides two sampling methods:

        1.  ``'optgp'`` to utilize the :class:`.ConcOptGPSampler`
        2.  ``'achr'`` to utilize the :class:`.ConcACHRSampler`

    2. The :class:`.ConcOptGPSampler` is parallel optimized sampler with
       fast convergence and parallel execution based on :cite:`MHM14`, with its
       implementation similar to the Python :mod:`cobra` package.
       See the :mod:`.conc_optgp` documentation for more information.

    3. The :class:`.ConcACHRSampler` is a sampler that utilizes an
       Artifial Centering Hit-and-Run (ACHR) sampler for a low memory
       footprint and good convergence based on :cite:`KS98`, with its
       implementation similar to the Python :mod:`cobra` package.
       See the :mod:`.conc_achr` documentation for more information.

    4. The :class:`.ConcHRSampler` is the base class for the
       samplers. All current samplers are derived from this class and all new
       samplers should be derived from this class to provide a unified
       interface for concentration sampling.
       See the :mod:`.conc_hr_sampler` documentation for more information.

To properly use the concentration samplers and associated functions, note
the following:

    * It is required that a model has been loaded into a :class:`.ConcSolver`
      instance and that the solver has been setup for concentraiton sampling
      via the :meth:`.ConcSolver.setup_sampling_problem` method.

    * All numerical values to be utilized by the solver must be defined. This
      includes:

        - Metabolite concentrations for all included metabolites, accessed via
          :attr:`.MassMetabolite.initial_condition`.
        - Reaction equilibrium constants for all included reactions, accessed
          via :attr:`.MassReaction.equilibrium_constant`.
        - Reaction steady state fluxes for all included reactions, accessed
          :attr:`.MassReaction.steady_state_flux`.

    * In order to perform the sampling, all numerical values are transformed
      from linear space into logarithmic space before the solver is populated.
      Therefore, all variables and constraints in the solver will exist in
      logspace.

      However, all numerical solution values are transformed from logspace
      back into a linear space before being returned, unless otherwise
      specified.

"""
from mass.thermo.conc_sampling.conc_achr import ConcACHRSampler
from mass.thermo.conc_sampling.conc_hr_sampler import ConcHRSampler
from mass.thermo.conc_sampling.conc_optgp import ConcOptGPSampler
from mass.thermo.conc_sampling.conc_sampling import sample_concentrations


__all__ = ()
