# -*- coding: utf-8 -*-
"""TODO Module Docstrings."""
import colorsys
import warnings

from mass.analysis import linear
from mass.core.visualization import plot_tiled_phase_portrait

import matplotlib as mpl

import numpy as np

from scipy import linalg

from six import integer_types, iteritems, itervalues

import sympy as sym

_T_SYM = sym.Symbol("t")


def calculate_modal_matrix(model, percentages=False, atol=1e-10,
                           time_invariants=False):
    """TODO DOCTSTRING."""
    if not isinstance(percentages, bool):
        raise TypeError("percentages must be a bool")
    # Calculate the jacobian
    J = linear.jacobian(model, jacobian_type="Jx", use_parameter_values=True,
                        use_concentration_values=True, matrix_type="dense")
    # Get eigenvectors and eigenvalues
    w, vr = linear.eig(J, right=True)
    # Get the modal matrix M
    vr_inv = linalg.inv(vr)

    if time_invariants:
        end = len(J)
    else:
        end = linear.matrix_rank(J, atol=atol)

    ts = np.real_if_close(-1/np.sort(w)[:end])
    M = np.array([vr_inv[i] for i in np.argsort(w)[:end]])

    # Normalize rows by largest weight
    M = np.array([np.real_if_close(r/max(abs(r))) for r in M])
    M = np.array([[0 if abs(v) < atol else v for v in r] for r in M])
    if percentages:
        M = np.array([r/sum(abs(r)) for r in M])

    return ts, M


def dynamic_pairwise_angles(M):
    """TODO DOCSTRINGS."""
    # if (not isinstance(cutoff, (float, integer_types)) or
    #    (cutoff > 1 or cutoff < 0)):
    #     raise ValueError("cutoff must be a float value between 0 and 1")
    def calc_angle(u, v):
        with warnings.catch_warnings():
            warnings.simplefilter(action="ignore",
                                  category=RuntimeWarning)
            return np.dot(u, v) / (linalg.norm(u) * linalg.norm(v))

    correlation_matricies = []
    M_i = M.copy()
    n_ts, n_mets = M_i.shape
    for k in range(n_ts):
        if k > 0:
            M_i[k-1] = np.zeros(n_mets)
        angle_matrix = np.zeros((n_mets, n_mets))
        for i in range(n_mets):
            for j in range(n_mets):
                angle_matrix[i, j] = calc_angle(M_i[:, i], M_i[:, j])

        correlation_matricies.append(angle_matrix)

    return np.array(correlation_matricies)


def tiled_time_correlated_portraits(solution, correlations, cutoff=.95,
                                    ax=None, cmap="RdYlGn", **kwargs):
    """TODO DOCSTRING."""
    if (not isinstance(cutoff, (float, integer_types)) or
       (cutoff > 1 or cutoff < 0)):
        raise ValueError("cutoff must be a float value between 0 and 1")

    if (len(solution), len(solution)) != correlations[0].shape:
        raise ValueError("Each correlation matrix must be an N x N matrix "
                         "where N is the number of solutions in the given "
                         "Solution object. ")

    cmap = mpl.cm.get_cmap(cmap)
    tiled_default_args = {"time_poi": None,
                          "poi_labels": True,
                          "empty_tiles": "upper"}

    if kwargs is not None:
        for kwarg, value in iteritems(kwargs.copy()):
            if kwarg in ["time_poi", "poi_labels", "empty_tiles"]:
                tiled_default_args[kwarg] = value
                del kwargs[kwarg]

    time_poi, poi_labels, empty_tiles = tuple(itervalues(tiled_default_args))

    n_sols = len(solution)
    data = np.zeros((n_sols, n_sols), dtype=np.int64)
    for k, angle_matrix in enumerate(correlations):
        for i, row in enumerate(angle_matrix):
            for j, value in enumerate(row):
                if abs(value) >= cutoff and i != j and data[i, j] == 0:
                    data[i, j] = k

    ax = plot_tiled_phase_portrait(solution, ax=ax, time_poi=time_poi,
                                   poi_labels=poi_labels,
                                   display_data=data, empty_tiles=empty_tiles,
                                   **kwargs)

    # Set face colors based on time scales
    _set_correlation_face_color(ax, cmap, k)

    return ax


# Internal
def _set_correlation_face_color(ax, cmap, k):
    """Set the face colors for data tiles based on the number of timescales.

    Warnings
    --------
    This method is intended for internal use only.

    """
    colors = cmap(np.linspace(0, 1, k+1))
    for subax in ax.child_axes:
        if subax.texts:
            value = float(subax.texts[0].get_text())
            c = colors[int(value)]
            c = colorsys.rgb_to_hls(*mpl.colors.to_rgb(c))
            c = colorsys.hls_to_rgb(c[0], 1 - .8 * (1 - c[1]), c[2])
            subax.set_facecolor(c)


# def get_modes(model, percentages=False, atol=1e-10, cutoff=None,
#               sympy_expr=True):
#     """TODO DOCSTRINGS."""
#     if cutoff is not None and (cutoff > 1 or cutoff < 0):
#         raise ValueError("cutoff must be a value between 0 and 1")
#     # Get modal matrix and metabolite function symbols
#     ts, M = calculate_modal_matrix(model, percentages=percentages, atol=atol)
#     metabolites = [sym.Function(m.id)(_T_SYM) for m in model.metabolites]
#     modes = []
#     # Iterate through the modal matrix
#     for r in M:
#         if cutoff is not None:
#
#             abs_r = abs(r).tolist()
#             total = sum(abs_r)
#             current = 0
#             indicies = []
#             # Determine the next largest weight for the row, sum with previous
#             # largest weights and compare to cutoff.
#             while current/total < cutoff:
#                 current += max(abs_r)
#                 ind = np.argmax(abs_r)
#                 indicies.append(ind)
#                 abs_r[ind] = 0
#                 # Break out if all entries are 0.
#                 if current == 0:
#                     indicies = list(range(len(r)))
#                     break
#             # Get the values and metabolites that sum to the cutoff
#             r = sym.Matrix([v for i, v in enumerate(r) if i in indicies])
#             mets = sym.Matrix([m for i, m in enumerate(metabolites)
#                                if i in indicies])
#         else:
#             # Use all values and metabolites
#             r = sym.Matrix([v for i, v in enumerate(r)])
#             mets = sym.Matrix([m for i, m in enumerate(metabolites)])
#         # Create the mode ODE
#         modes.extend(list(r.T * mets))
#     # Turn into strings if desired, otherwise leave them as sympy expressions
#     modes = np.array([mode if sympy_expr else str(strip_time(mode))
#                       for mode in modes])
#     return modes
