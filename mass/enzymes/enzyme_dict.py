# -*- coding: utf-8 -*-
"""TODO Module Docstrings."""
from __future__ import absolute_import

from collections import OrderedDict
from copy import copy, deepcopy

import numpy as np

import pandas as pd

from six import iteritems, iterkeys, itervalues

from cobra.core.dictlist import DictList

from mass.util.DictWithID import OrderedDictWithID
from mass.util.util import (
    _get_matrix_constructor, _mk_new_dictlist, convert_matrix)


class EnzymeDict(OrderedDictWithID):
    """Container to store the attributes of an EnzymeModel.

    This object is intended to represent the EnzymeModel once merged into
    a MassModel in order to retain EnzymeModel specific attributes of the 
    EnzymeModel without the need of storing the EnzymeModel object itself.

    The EnzymeDict class is essentially a subclass of an OrderedDict with an id
    and attribute accessors. One can see which attributes are accessible using
    the keys() method.

    Parameters
    ----------
    enzyme: mass.EnzymeModel, mass.EnzymeDict
        The EnzymeModel to be converted into the EnzymeDict, or an existing 
        EnzymeDict object. If an existing EnzymeDict object is
        provided, a new EnzymeDict object is instantiated with the same
        properties as the original EnzymeDict.

    Notes
    -----
    EnzymeDict objects can behave as booleans, with empty EnzymeDict objects
        returning as False.

    """

    def __init__(self, id_or_enzyme=None):
        """Initialize EnzymeDict."""
        # Instiantiate a new EnzymeDict object if a EnzymeDict is given.
        if isinstance(id_or_enzyme, (EnzymeDict, dict)):
            super(EnzymeDict, self).__init__(
                id_or_enzyme["id"], dictionary=dict(id_or_enzyme))
        # Initialize an EnzymeDict using an EnzymeModel
        elif hasattr(id_or_enzyme, "_convert_self_into_enzyme_dict"):
            super(EnzymeDict, self).__init__(id_or_enzyme.id)
            for key, value in iteritems(id_or_enzyme.__dict__):
                new_key = key.lstrip("_")
                if new_key not in iterkeys(_ORDERED_ENZYMEDICT_DEFAULTS):
                    continue
                elif new_key is "S":
                    self[new_key] = id_or_enzyme._mk_stoich_matrix(
                        matrix_type="DataFrame", update_model=False)
                else:
                    self[new_key] = value
        # Initialize EnzymeDict with an id and defaults for everything else.
        else:
            self["id"] = id_or_enzyme
        # Ensure all required attributes make it into the EnzymeDict
        self._set_missing_to_defaults()
        # Move entries to keep a consisent order for ordered EnzymeDicts.
        self._fix_order()

    def copy(self):
        """Copy an EnzymeDict object."""
        # No references to the MassModel when copying the EnzymeDict
        model = self.model
        setattr(self, "model", None)
        for attr in ["ligands", "enzyme_forms", "enzyme_reactions"]:
            for i in getattr(self, attr):
                setattr(i, "_model", None)
            for dictlist in itervalues(self["categorized_" + attr]):
                for i in dictlist:
                    setattr(i, "_model", None)

        # The EnzymeDict can now be copied
        enzyme_dict_copy = deepcopy(self)
        # Restore references for the original EnzymeDict
        setattr(self, "model", model)
        for attr in ["ligands", "enzyme_forms", "enzyme_reactions"]:
            for i in getattr(self, attr):
                setattr(i, "_model", model)
            for dictlist in itervalues(self["categorized_" + attr]):
                for i in dictlist:
                    setattr(i, "_model", model)

        return enzyme_dict_copy

    # Internal
    def _set_missing_to_defaults(self):
        """Set all of the missing attributes to their default values.

        Warnings
        --------
        This method is intended for internal use only. 

        """
        for key, value in iteritems(_ORDERED_ENZYMEDICT_DEFAULTS):
            if key not in self:
                self[key] = value

    def _update_object_pointers(self, model=None):
        """Update objects in the EnzymeDict to come from the associated model.

        Warnings
        --------
        This method is intended for internal use only. 

        """
        if model is None:
            model = self.model
        for attr in ["ligands", "enzyme_forms", "enzyme_reactions"]:
            model_dictlist = {
                "ligands": model.metabolites, 
                "enzyme_forms": model.metabolites,
                "enzyme_reactions": model.reactions}.get(attr)
            self[attr] = _mk_new_dictlist(model_dictlist, getattr(self, attr))
            attr = "categorized_" + attr
            self[attr] = {
                key: _mk_new_dictlist(model_dictlist, old_dictlist)
                for key, old_dictlist in iteritems(getattr(self, attr))}

    def _make_enzyme_stoichiometric_matrix(self, update=False):
        """Return the S matrix based on enzyme forms, reactions, and ligands.

        Warnings
        --------
        This method is intended for internal use only. 

        """
        # Set up for matrix construction.
        (matrix_constructor, matrix_type, dtype) = _get_matrix_constructor(
            matrix_type="DataFrame", dtype=np.float_)

        metabolites = DictList([
            met for attr in ["ligands", "enzyme_forms"] for met in self[attr]])

        stoich_mat = matrix_constructor((len(metabolites),
                                         len(self.enzyme_reactions)))
        # Get the indicies for the species and reactions
        m_ind = metabolites.index
        r_ind = self.enzyme_reactions.index

        # Build the matrix
        for rxn in self.enzyme_reactions:
            for met, stoich in iteritems(rxn.metabolites):
                stoich_mat[m_ind(met), r_ind(rxn)] = stoich

        # Convert the matrix to the desired type
        stoich_mat = convert_matrix(
            stoich_mat, matrix_type=matrix_type, dtype=dtype, 
            row_ids=[m.id for m in metabolites],
            col_ids=[r.id for r in self.enzyme_reactions])

        if update:
            self["S"] = stoich_mat

        return stoich_mat

    def _fix_order(self):
        """Fix the order of the items in the EnzymeDict.

        Warnings
        --------
        This method is intended for internal use only. 

        """
        for key in iterkeys(_ORDERED_ENZYMEDICT_DEFAULTS):
            if key in self:
                self.move_to_end(key)

    def _repr_html_(self):
        """HTML representation of the overview for the EnzymeDict."""
        try:
            dim_S = "{0}x{1}".format(self.S.shape[0], self.S.shape[1])
            rank = np.linalg.matrix_rank(self.S)
        except np.linalg.LinAlgError:
            dim_S = "0x0"
            rank = 0

        return """
            <table>
                <tr>
                    <td><strong>Name</strong></td><td>{name}</td>
                </tr><tr>
                    <td><strong>Memory address</strong></td><td>{address}</td>
                </tr><tr>
                    <td><strong>Stoichiometric Matrix</strong></td>
                    <td>{dim_stoich_mat}</td>
                </tr><tr>
                    <td><strong>Matrix Rank</strong></td>
                    <td>{mat_rank}</td>
                </tr><tr>
                    <td><strong>Subsystem</strong></td>
                    <td>{subsystem}</td>
                </tr><tr>
                    <td><strong>Number of Ligands</strong></td>
                    <td>{num_ligands}</td>
                </tr><tr>
                    <td><strong>Number of EnzymeForms</strong></td>
                    <td>{num_enz_forms}</td>
                </tr><tr>
                    <td><strong>Number of Enzyme Reactions</strong></td>
                    <td>{num_enz_reactions}</td>
                </tr><tr>
                    <td><strong>Total Enzyme Concentration</strong></td>
                    <td>{enz_conc}</td>
                </tr><tr>
                    <td><strong>Enzyme Net Flux</strong></td>
                    <td>{enz_flux}</td>
                </tr>
            </table>
        """.format(name=self.id, address='0x0%x' % id(self),
                   dim_stoich_mat=dim_S, mat_rank=rank,
                   subsystem=self.subsystem,
                   num_ligands=len(self.ligands),
                   num_enz_forms=len(self.enzyme_forms),
                   num_enz_reactions=len(self.enzyme_reactions),
                   enz_conc=self.enzyme_concentration_total, 
                   enz_flux=self.enzyme_net_flux)

    # Dunders
    def __getattr__(self, name):
        """Override of default getattr() implementation."""
        if name == "_id" or name == "_model":
            name = name.lstrip("_")
        if name in self:
            return self[name]
        else:
            raise AttributeError("Attribute '" + name + "' does not exist.")

    def __setattr__(self, name, value):
        """Override of default setattr() implementation."""
        if name == "_id" or name == "_model":
            name = name.lstrip("_")
        self[name] = value

    def __delattr__(self, name):
        """Override of default delattr() implementation."""
        if name in self:
            del self[name]
        else:
            raise AttributeError("Attribute '" + name + "' does not exist.")

    def __dir__(self):
        """Override of default dir() implementation to include the keys."""
        return sorted(set(
            list(iterkeys(self)) + super(EnzymeDict, self).__dir__()))

    def __copy__(self):
        """Create a copy of the EnzymeDict."""
        return copy(super(EnzymeDict, self))

    def __deepcopy__(self, memo):
        """Create a deepcopy of the EnzymeDict."""
        return deepcopy(super(EnzymeDict, self), memo)

    def __repr__(self):
        """Override of default repr() implementation."""
        return "<%s %s at 0x%x>" % (
            self.__class__.__name__[:-4], self.id, id(self))


_ORDERED_ENZYMEDICT_DEFAULTS = OrderedDict({
    "id": None,
    "name": None,
    "subsystem": "",
    "ligands": DictList(),
    "enzyme_forms": DictList(),
    "enzyme_reactions": DictList(),
    "categorized_ligands": {"Undefined": DictList()},
    "categorized_enzyme_forms": {"Undefined": DictList()},
    "categorized_enzyme_reactions": {"Undefined": DictList()},
    "enzyme_concentration_total": None,
    "enzyme_net_flux": None,
    "enzyme_net_flux_equation": None,
    "description": "",
    "S": pd.DataFrame(),
    "model": None,
})
