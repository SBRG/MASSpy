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


class EnzymeModuleDict(OrderedDictWithID):
    """Container to store the attributes of an EnzymeModule.

    This object is intended to represent the EnzymeModule once merged into
    a MassModel in order to retain EnzymeModule specific attributes of the
    EnzymeModule without the need of storing the EnzymeModule object itself.

    The EnzymeModuleDict class is essentially a subclass of an OrderedDict with
    an id and attribute accessors. One can see which attributes are accessible
    using the keys() method.

    Parameters
    ----------
    enzyme: EnzymeModule, EnzymeModuleDict
        The EnzymeModule to be converted into the EnzymeModuleDict, or an
        existing EnzymeModuleDict object. If an existing EnzymeModuleDict
        object is provided, a new EnzymeModuleDict object is instantiated with
        the same properties as the original EnzymeModuleDict.

    Notes
    -----
    EnzymeModuleDict objects can behave as booleans, with empty
        EnzymeModuleDict objects returning as False.

    """

    def __init__(self, id_or_enzyme=None):
        """Initialize EnzymeModuleDict."""
        # Instiantiate new EnzymeModuleDict object if given an EnzymeModuleDict
        if isinstance(id_or_enzyme, (EnzymeModuleDict, dict)):
            super(EnzymeModuleDict, self).__init__(
                id_or_enzyme["id"], dictionary=dict(id_or_enzyme))
        # Initialize an EnzymeModuleDict using an EnzymeModule
        elif hasattr(id_or_enzyme, "_convert_self_into_enzyme_dict"):
            super(EnzymeModuleDict, self).__init__(id_or_enzyme.id)
            for key, value in iteritems(id_or_enzyme.__dict__):
                nkey = key.lstrip("_")
                if nkey not in iterkeys(_ORDERED_ENZYMEMODULE_DICT_DEFAULTS):
                    continue
                elif nkey == "S":
                    self[nkey] = id_or_enzyme._mk_stoich_matrix(
                        matrix_type="DataFrame", update_model=False)
                elif nkey == "enzyme_net_flux_equation":
                    self[nkey] = id_or_enzyme.enzyme_net_flux_equation
                else:
                    self[nkey] = value
        # Initialize EnzymeModuleDict with an id and defaults for everything.
        else:
            self["id"] = id_or_enzyme
        # Ensure all required attributes make it into the EnzymeModuleDict
        self._set_missing_to_defaults()
        # Move entries to keep a consisent order for ordered EnzymeModuleDicts.
        self._fix_order()

    def copy(self):
        """Copy an EnzymeModuleDict object."""
        # No references to the MassModel when copying the EnzymeModuleDict
        model = self.model
        setattr(self, "model", None)
        for attr in ["enzyme_module_ligands", "enzyme_module_forms",
                     "enzyme_module_reactions"]:
            for i in getattr(self, attr):
                setattr(i, "_model", None)
            for dictlist in itervalues(self[attr + "_categorized"]):
                for i in dictlist:
                    setattr(i, "_model", None)

        # The EnzymeModuleDict can now be copied
        enzyme_dict_copy = deepcopy(self)
        # Restore references for the original EnzymeModuleDict
        setattr(self, "model", model)
        for attr in ["enzyme_module_ligands", "enzyme_module_forms",
                     "enzyme_module_reactions"]:
            for i in getattr(self, attr):
                setattr(i, "_model", model)
            for dictlist in itervalues(self[attr + "_categorized"]):
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
        for key, value in iteritems(_ORDERED_ENZYMEMODULE_DICT_DEFAULTS):
            if key not in self:
                self[key] = value

    def _update_object_pointers(self, model=None):
        """Update objects in the EnzymeModuleDict to point to the given model.

        Warnings
        --------
        This method is intended for internal use only.

        """
        if model is None:
            model = self.model
        for attr in ["enzyme_module_ligands", "enzyme_module_forms",
                     "enzyme_module_reactions"]:
            model_dictlist = {
                "enzyme_module_ligands": model.metabolites,
                "enzyme_module_forms": model.metabolites,
                "enzyme_module_reactions": model.reactions}.get(attr)
            self[attr] = _mk_new_dictlist(model_dictlist, getattr(self, attr))
            attr += "_categorized"
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
            met for attr in ["enzyme_module_ligands", "enzyme_module_forms"]
            for met in self[attr]])

        stoich_mat = matrix_constructor((len(metabolites),
                                         len(self.enzyme_module_reactions)))
        # Get the indicies for the species and reactions
        m_ind = metabolites.index
        r_ind = self.enzyme_module_reactions.index

        # Build the matrix
        for rxn in self.enzyme_module_reactions:
            for met, stoich in iteritems(rxn.metabolites):
                stoich_mat[m_ind(met), r_ind(rxn)] = stoich

        # Convert the matrix to the desired type
        stoich_mat = convert_matrix(
            stoich_mat, matrix_type=matrix_type, dtype=dtype,
            row_ids=[m.id for m in metabolites],
            col_ids=[r.id for r in self.enzyme_module_reactions])

        if update:
            self["S"] = stoich_mat

        return stoich_mat

    def _fix_order(self):
        """Fix the order of the items in the EnzymeModuleDict.

        Warnings
        --------
        This method is intended for internal use only.

        """
        for key in iterkeys(_ORDERED_ENZYMEMODULE_DICT_DEFAULTS):
            if key in self:
                self.move_to_end(key)

    def _repr_html_(self):
        """HTML representation of the overview for the EnzymeModuleDict."""
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
                    <td>{num_enzyme_module_ligands}</td>
                </tr><tr>
                    <td><strong>Number of EnzymeForms</strong></td>
                    <td>{num_enz_forms}</td>
                </tr><tr>
                    <td><strong>Number of EnzymeModuleReactions</strong></td>
                    <td>{num_enz_reactions}</td>
                </tr><tr>
                    <td><strong>Enzyme Concentration Total</strong></td>
                    <td>{enz_conc}</td>
                </tr><tr>
                    <td><strong>Enzyme Net Flux</strong></td>
                    <td>{enz_flux}</td>
                </tr>
            </table>
        """.format(name=self.id, address='0x0%x' % id(self),
                   dim_stoich_mat=dim_S, mat_rank=rank,
                   subsystem=self.subsystem,
                   num_enzyme_module_ligands=len(self.enzyme_module_ligands),
                   num_enz_forms=len(self.enzyme_module_forms),
                   num_enz_reactions=len(self.enzyme_module_reactions),
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
            list(iterkeys(self)) + super(EnzymeModuleDict, self).__dir__()))

    def __copy__(self):
        """Create a copy of the EnzymeModuleDict."""
        return copy(super(EnzymeModuleDict, self))

    def __deepcopy__(self, memo):
        """Create a deepcopy of the EnzymeModuleDict."""
        return deepcopy(super(EnzymeModuleDict, self), memo)

    def __repr__(self):
        """Override of default repr() implementation."""
        return "<%s %s at 0x%x>" % (
            self.__class__.__name__[:-4], self.id, id(self))


_ORDERED_ENZYMEMODULE_DICT_DEFAULTS = OrderedDict({
    "id": None,
    "name": None,
    "subsystem": "",
    "enzyme_module_ligands": DictList(),
    "enzyme_module_forms": DictList(),
    "enzyme_module_reactions": DictList(),
    "enzyme_module_ligands_categorized": {"Undefined": DictList()},
    "enzyme_module_forms_categorized": {"Undefined": DictList()},
    "enzyme_module_reactions_categorized": {"Undefined": DictList()},
    "enzyme_concentration_total": None,
    "enzyme_net_flux": None,
    "enzyme_net_flux_equation": None,
    "description": "",
    "S": pd.DataFrame(),
    "model": None,
})
