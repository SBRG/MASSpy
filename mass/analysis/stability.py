# -*- coding: utf-8 -*-

# Compatibility with Python 2.7
from __future__ import absolute_import

# Import necesary packages
import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns

from sympy.core import Symbol
from scipy.interpolate import interp1d
from math import inf, log10
from six import iterkeys, itervalues

# from cobra
from cobra import DictList

# from mass
from mass.core.massmetabolite import MassMetabolite
from mass.core.massreaction import MassReaction
from mass.core.massmodel import MassModel
from mass.core.expressions import _sort_symbols
from mass.core.simulation import *
from mass.analysis.linear import *



# Methods to generate perturbed eigenvalues
def get_Jx_sym(massmodel, jacobian_type="metabolite"):
    """Get symbolic Jacobian

    Parameters
    ----------
    massmodel: mass.MassModel
        The system model for which to find the Jacobian of
    jacobian_type : "metabolite" or "reaction"
        Whether to obtain the jacobian with respect to the metabolites (Jx)
        or to obbtain the jacobian with respect to the reactions (Jv)

    Returns
    -------
    sympy.MutableDenseMatrix
        A sympy matrix representation of the Jacobian (in symbolic form)
    """
    return jacobian(massmodel, jacobian_type=jacobian_type, 
                    sub_parameters=False, sub_concentrations=False,
                    matrix_type="symbolic")

def fill_ssconc(Jx, massmodel, strategy="simulate"):
    """Fill Jx with steady-state conc data

    Parameters
    ----------
    Jx: sympy.MutableDenseMatrix
        A sympy matrix representation of the Jacobian (in symbolic form)
    massmodel: mass.MassModel
        The system model for which to find the Jacobian
    strategy : "simulate" or "find_roots"
        The strategy to use to solve for the steady state

    Returns
    -------
    sympy.MutableDenseMatrix
        A sympy matrix representation of the Jacobian (in symbolic form)
        with numeric values for everything except parameters
    """
    cf, ff = find_steady_state(massmodel, strategy)

    #fill in ssconc values
    new_keys = [Symbol(item.id) for item in cf.keys()]
    new_dict = dict(zip(new_keys, cf.values()))
    Jx = Jx.subs(new_dict)

    #fill in fixed ssconc values
    Jx = Jx.subs(massmodel.fixed_concentrations)
    new_keys = [Symbol(item) for item in massmodel.fixed_concentrations.keys()]
    new_dict = dict(zip(new_keys, massmodel.fixed_concentrations.values()))
    Jx = Jx.subs(new_dict)

    return Jx

def get_params(massmodel):
    """Get parameters to be tweaked

    Parameters
    ----------
    massmodel: mass.MassModel
        The system model which contains the relevent parameter values

    Returns
    -------
    list(syms[1]): list of sympy.Symbol
        A list of parameters in the symbolic Jacobian which can be perturbed
    react_dict: dict
        A dictionary containing parameters as keys (sympy.Symbol)
        and their model values as values (int, float, np.int64, np.float64)
    param_dict: dict
        A dictionary with kf parameters as keys (sympy.Symbol)
        and Keq parameters as values (sympy.Symbol)
    """
    odes, rates, syms = _sort_symbols(massmodel)

    react_dict = {}
    param_dict = {}
    a = None
    b = None

    massmodel.generate_rate_laws(rate_type=1)
    for react in massmodel.reactions:
        for item in list(syms[1]):
            if react.id in str(item):

                if "kf" in str(item):
                    a = item
                    react_dict[item] = react.kf

                elif "Keq" in str(item):
                    b = item
                    if str(a)[2:] in str(b)[3:]:
                        param_dict[a] = b
                    react_dict[item] = react.Keq

                elif "kr" in str(item):
                    react_dict[item] = react.kr

    return list(syms[1]), react_dict, param_dict

def gen_eigvalues(
    Jx, param, react_dict, scale="logspace", *custom_values, **values):
    """Generate dataframe of eigenvalues given a parameter to perturb
       and (optional) a range over which to perturb it

    Parameters
    ----------
    Jx: sympy.MutableDenseMatrix
        A sympy matrix representation of the Jacobian (in symbolic form)
        with numeric values for everything except parameters
    param: sympy.Symbol or str
        The specific parameter to perturb over a particular range in [0, inf]
        default range is [1e-9, 1e6] for linspace and [1e-9, 1e9] for logspace
    react_dict: dict
        A dictionary containing parameters as keys (sympy.Symbol)
        and their model values as values (int, float, np.int64, np.float64)
    scale: "linspace" or "logspace" or "custom"
        Whether to use np.linspace() or np.logspace() to generate parameter
        range for perturbing the value of param. If custom, a custom vector
        from *custom_values is used instead
    *custom_values: tuple
        A tuple of values to be used as a custom vector for perturbation
        of param. Can be a list or np.ndarray as long as "*" is present and
        attached to vector
        Examples:
            *list(1,2,3)
                Range is 1 to 3
            *np.array([1,2,3])
                Range is 1 to 3
    **values: dict
        keyword arguments for changing default range of linspace/logspace
        methods. These arguments are passed in as parameters into these methods
        directly. See documentation on np.linspace and/or np.logspace for
        details. 
        Examples:
            linspace: **{"start": 10, "stop": 100, "num": 100}
                Range is np.linspace(start=10, stop=100, num=100)
            logspace: **{"start": 1, "stop": 10, "num": 100}
                Range is np.logspace(start=1, stop=10, num=100)

    Returns
    -------
    eigs_df: pd.DataFrame
        A pandas dataframe, indexed by the values of a pertubed parameter and
        containing the eigenvalues for each value of that perturbed parameter
    orig: float
        The original steady-state value of the perturbed parameter
    """
    Jx, orig = _pop_param(Jx, param, react_dict)

    prange = _set_scale(scale, *custom_values, **values)
    eigs_df = pd.DataFrame()
    Jcopy = Jx.copy()

    for i in prange:
        eigs = np.real(eigenvalues(jj))
        eigs_df[i] = eigs

    eigs_df = eigs_df.T
    eigs_df.index.name = str(param)

    return eigs_df, orig



# Methods to process generated eigenvalues
def normalize_eigenvalues(Jx, param, react_dict, eigs_df):
    """Divide each column of eigs_df by each eigenvalue's original 
       steady-state value

    Parameters
    ----------
    Jx: sympy.MutableDenseMatrix
        A sympy matrix representation of the Jacobian (in symbolic form)
        with numeric values for everything except parameters
    param: sympy.Symbol or str
        The specific parameter to perturb over a particular range in [0, inf]
        default range is [1e-9, 1e6] for linspace and [1e-9, 1e9] for logspace
    react_dict: dict
        A dictionary containing parameters as keys (sympy.Symbol)
        and their model values as values (int, float, np.int64, np.float64)
    eigs_df: pd.DataFrame
        A pandas dataframe, indexed by the values of a pertubed parameter and
        containing the eigenvalues for each value of that perturbed parameter

    Returns
    -------
    norm_eigs_df: pd.DataFrame
        A pandas dataframe, indexed by the values of a pertubed parameter and
        containing the normalized eigenvalues for each value of that perturbed
        parameter
    """
    Jx, orig = _pop_param(Jx, param, react_dict)
    Jcopy = Jx.copy()

    jj = np.array(Jcopy.subs({param: orig}).tolist()).astype(np.float64)
    eigs = np.real(eigenvalues(jj))

    norm_eigs_df = eigs_df.T.apply(lambda x: x/eigs)
    norm_eigs_df = norm_eigs_df.T

    return norm_eigs_df

def applymap_log_wsp(eigs_df):
    """Applies log_wsp function to each element in eigs_df

    Parameters
    ----------
    norm_eigs_df: pd.DataFrame
        A pandas dataframe, indexed by the values of a pertubed parameter and
        containing the log_wsp of the normalized eigenvalues for each value 
        of that perturbed parameter
    Returns
    -------
    pd.DataFrame
        A pandas dataframe, indexed by the values of a pertubed parameter and
        containing the log_wsp of the normalized eigenvalues for each value 
        of that perturbed parameter
    """
    return eigs_df.applymap(lambda x: _log_wsp(x))



# Methods to plot dataframes of eigenvalues
def plot_eigvalues(eigs_df, size=(15,15), linewidth=None, fname=None):
    """Generate heatmap of eigs_df"""
    fig, ax = plt.subplots()
    fig.set_size_inches(size[0],size[1])

    if linewidth is not None:
        sns.heatmap(eigs_df, ax=ax, linewidth=linewidth)
        plt.show()
    else:
        sns.heatmap(eigs_df, ax=ax)
        plt.show()

    if fname is not None:
        fig.savefig(fname, format="png")

def gen_perturbed_parameter_eigenvalue_plots(massmodel, 
                                             jacobian_type="metabolite", 
                                             strategy="simulate", 
                                             scale="logspace", 
                                             size=(15,15), 
                                             linewidth=0.1,
                                             dname=None, 
                                             custom_values=None, **values):
    """FIXME"""
    #generate symbolic Jacobian
    Jx = get_Jx_sym(massmodel, jacobian_type)
    Jx = fill_ssconc(Jx, massmodel, strategy)

    #find list of parameters that can be perturbed
    list_syms, react_dict, param_dict = get_params(massmodel)
    infs = [item for item in react_dict.keys() if react_dict[item] == inf]

    #generate, process, and plot eigenvalue data
    for param in list_syms:
        if param not in infs:
            eigs_df, orig = gen_eigvalues(
                Jx, param, react_dict, scale, *custom_values, **values)

            norm_eigs_df = normalize_eigenvalues(
                Jx, param, react_dict, eigs_df)

            log_wsp_df = applymap_log_wsp(norm_eigs_df)

            if dname is not None and isinstance(dname, str):
                if dname[-1] not in "/":
                    name = dname+str(param)
                elif dname[-1] in "/":
                    name = dname+"/"+str(param)
            else:
                name=None
            plot_eigvalues(
                log_wsp_df, size, linewidth, name)



# Internal Methods
def _pop_param(Jx, param, react_dict):
    if isinstance(param, str):
        param = Symbol(param)
    elif isinstance(param, Symbol):
        pass
    else:
        msg = "param must be of type str or sympy.Symbol"
        raise TypeError(msg)

    subs_dict = dict(react_dict)
    orig = subs_dict.pop(param)
    Jx = Jx.subs(subs_dict)

    return Jx, orig

def _set_scale(scale, *custom_values, **values):
    """Generates a linear or logarithmic scaled vector of values for
       perturbation. A custom vector can be added as well using *custom_values.
       Default options can be tweaked using **values.

    Parameters
    ----------
    scale: "linspace" or "logspace" or "custom"
        Whether to use np.linspace() or np.logspace() to generate parameter
        range for perturbing the value of param. If custom, a custom vector
        from *custom_values is used instead
    *custom_values: tuple
        A tuple of values to be used as a custom vector for perturbation
        of param
    **values: dict
        keyword arguments for changing default range of linspace/logspace
        methods. These arguments are passed in as parameters into these methods
        directly. See documentation on np.linspace and/or np.logspace for
        details. 
        Examples:
            linspace: **{"start": 10, "stop": 100, "num": 100}
                Range is np.linspace(start=10, stop=100, num=100)
            logspace: **{"start": 1, "stop": 10, "num": 100}
                Range is np.logspace(start=1, stop=10, num=100)

    Returns
    -------
    prange: np.ndarray
        A numpy array of values to be used for perturbing a parameter

    """

    if scale.lower() in "logspace":
        if values == {}:
            values = {"start": -9, "stop": 9, "num": 100}
        prange = np.logspace(**values)

    elif scale.lower() in "linspace":
        if values == {}:
            a = np.linspace(1e-9, 1, 34, endpoint=False)
            b = np.linspace(1, 1e3, 33, endpoint=False)
            c = np.linspace(1e3, 1e6, 33)
            prange = np.append(a, np.append(b, c))
        else:
            prange = np.linspace(**values)

    elif scale.lower() in "custom":
        if custom_values == () or custom_values is None:
            msg = "'custom' requires a non-empty '*custom_values' parameter"
            raise ValueError(msg)
        else:
            prange = np.array(custom_values)

    else:
        msg = "scale must be either linspace or logspace"
        raise ValueError(msg)

    return prange

def _log_wsp(x):
    if x == 0:
        return 0
    else:
        result = log10(x) if x > 0 else -1*log10(-x)
        return result

















