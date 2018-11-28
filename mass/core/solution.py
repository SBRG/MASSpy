# -*- coding: utf-8 -*-
"""TODO Module Docstrings."""
from __future__ import absolute_import

from copy import copy, deepcopy
import re

from mass.core.massmodel import MassModel

import pandas as pd

from scipy.interpolate import interp1d

from six import iteritems, string_types


# Strings of valid solution types
_CONC_STR = "Conc"
_FLUX_STR = "Flux"
_POOL_STR = "Pool"
_NETFLUX_STR = "NetFlux"

# Pre-compiled regular expressions for valid solution types
_conc_re = re.compile(_CONC_STR)
_flux_re = re.compile(_FLUX_STR)
_pool_re = re.compile(_POOL_STR)
_netflux_re = re.compile(_NETFLUX_STR)


class _DictWithID(dict):
    """A standard python dictionary with an ID attribute.

    Parameters
    ----------
    id: str, None
        The identifier to associate with the object.
    dictionary: dict, optional
        If provided, the new _DictWithID will contain the key:value pairs from
        "dictionary". Otherwise the new _DictWithID is initialized as empty.

    Warnings
    --------
    This object is primarily intended for internal use only.

    """

    def __init__(self, id, dictionary=None):
        super(_DictWithID, self).__init__(self)
        self._id = id
        if isinstance(dictionary, dict):
            self.update(dictionary)
        elif dictionary is None:
            pass
        else:
            raise TypeError("dictionary must be a dict")

    @property
    def id(self):
        return getattr(self, "_id", None)

    @id.setter
    def id(self, value):
        if value == self.id:
            pass
        elif not isinstance(value, string_types):
            raise TypeError("ID must be a string")
        elif getattr(self, "_model", None) is not None:
            self._set_id_with_model(value)
        else:
            self._id = value

    def _set_id_with_model(self, value):
        self._id = value


class Solution(_DictWithID):
    """Container to store the solutions for the simulation of a MassModel.

    Solution containers are given an ID of the following form:
        Solution.id = Solution_{id_or_model}_{solution_type}

    Parameters
    ----------
    id_or_model: str, MassModel
        A string identifier or a MassModel to associate with the solution. If
        a MassModel, the MassModel identifier will be used, and a reference to
        the model is stored.
    solution_type: str
        The type of solution being stored. Can be one of the following:
        {"concentration", "flux"}


    Notes
    -----
    Solution objects can behave as booleans, with empty Solution objects
        returning as False, and those with solutions returning as True.

    """

    def __init__(self, id_or_model, solution_type, solution_dictionary=None,
                 time=None, interpolate=False):
        """Container with an identifier for "solution_type" solutions."""
        valid_check = [True if _re.match(solution_type) else False
                       for _re in [_conc_re, _flux_re, _pool_re, _netflux_re]]
        if True not in valid_check:
            raise ValueError("'{0}' is not a valid solution type."
                             .format(solution_type))
        if isinstance(id_or_model, MassModel):
            self._model = id_or_model
            id_or_model = id_or_model.id
        else:
            self._model = None

        id_or_model = "{0}_{1}Sol".format(id_or_model, solution_type)

        _DictWithID.__init__(self, id_or_model, dictionary=solution_dictionary)
        self._solution_type = solution_type
        self._time = time
        self._interpolate = interpolate
        self._groups = None

    @property
    def t(self):
        """Return time points used in computing the solution."""
        return getattr(self, "_time", None)

    @property
    def df(self):
        """Return the Solution object as a pandas.DataFrame."""
        if self.interpolate:
            sols = {k: v(self.t) for k, v in iteritems(self.solutions)}
        else:
            sols = self.solutions
        df = pd.DataFrame.from_dict(sols)
        df.index = pd.Series(self.t, name="Time")
        return df

    @property
    def solutions(self):
        """Return a normal dict containing solutions of the Solution object."""
        return dict(self)

    @property
    def solution_type(self):
        """Return the solution type of the Solution object."""
        return getattr(self, "_solution_type", None)

    @property
    def t0(self):
        """Return the initial time point used in computing the solution."""
        if self.t is not None and (len(self.t) >= 2):
            return self.t[0]

    @property
    def tf(self):
        """Return the final time point used in computing the solution."""
        if self.t is not None and (len(self.t) >= 2):
            return self.t[-1]

    @property
    def numpoints(self):
        """Return the number of time points used to compute the solutions."""
        return int(len(self.t))

    @property
    def groups(self):
        """Return the groups used to create grouped solutions if they exist."""
        return getattr(self, "_groups", None)

    @property
    def interpolate(self):
        """Return whether solutions are given as interpolating functions."""
        return getattr(self, "_interpolate", None)

    @interpolate.setter
    def interpolate(self, value):
        """Set whether solutions are stored as interpolating functions.

        Parameters
        ----------
        value: bool
            If True, solutions are stored in the Solution object as
            scipy interpolating functions. Otherwise store solutions as arrays.

        """
        if not isinstance(value, bool):
            raise ValueError("value must be a bool")
        if value != self.interpolate:
            for key, sol in iteritems(self):
                if value:
                    if not isinstance(sol, interp1d):
                        self[key] = interp1d(self.t, sol, kind='cubic',
                                             fill_value='extrapolate')
                else:
                    self[key] = sol(self.t)

        self._interpolate = value

    @property
    def model(self):
        """Return the model used in creating the Solution if one was used."""
        return getattr(self, "_model", None)

    def copy(self):
        """Copy a Solution object."""
        return copy(self)

    def __copy__(self):
        """Create a copy of the Solution object."""
        return copy(super(Solution, self))

    def __repr__(self):
        """Representation of the Solution object."""
        return "<%s %s at 0x%x>" % (self.__class__.__name__, self.id, id(self))
