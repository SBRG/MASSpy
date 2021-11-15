# -*- coding: utf-8 -*-
r"""EnzymeModuleDict is a class representing the EnzymeModule after merging.

This object is intended to represent the :class:`.EnzymeModule` once merged
into a :class:`.MassModel` in order to retain :class:`.EnzymeModule` specific
attributes of the :class:`.EnzymeModule` without the need of storing the
:class:`.EnzymeModule` object itself.

When merging an :class:`.EnzymeModule` into another model, the
:class:`.EnzymeModule` is converted into an :class:`EnzymeModuleDict`, allowing
for most of the enzyme specific attribute information of the
:class:`EnzymeModule` to be preserved during the merging process, and accessed
after the merging. All keys of the :class:`EnzymeModuleDict` can be used as
attribute accessors.  Additionally :class:`EnzymeModuleDict` is a subclass of
an :class:`~.OrderedDictWithID` which in turn is a subclass of an
:class:`~.collections.OrderedDict`, thereby inheriting its methods and behavior.

The :class:`.EnzymeModule` attributes preserved in the
:class:`EnzymeModuleDict` are the following:

    * :attr:`.EnzymeModule.id`
    * :attr:`.EnzymeModule.name`
    * :attr:`.EnzymeModule.subsystem`
    * :attr:`.EnzymeModule.enzyme_module_ligands`
    * :attr:`.EnzymeModule.enzyme_module_forms`
    * :attr:`.EnzymeModule.enzyme_module_reactions`
    * :attr:`.EnzymeModule.enzyme_module_ligands_categorized`
    * :attr:`.EnzymeModule.enzyme_module_forms_categorized`
    * :attr:`.EnzymeModule.enzyme_module_reactions_categorized`
    * :attr:`.EnzymeModule.enzyme_concentration_total`
    * :attr:`.EnzymeModule.enzyme_rate`
    * :attr:`.EnzymeModule.enzyme_concentration_total_equation`
    * :attr:`.EnzymeModule.enzyme_rate_equation`
    * :attr:`.EnzymeModule.S`

If one of the above attributes has not been set, it will be added to the
:class:`EnzymeModuleDict` as its default value. This means that the above
attributes can *always* be found in an :class:`EnzymeModuleDict`.

Note that this class is not intended to be used for construction of an
:class:`.EnzymeModule`, but rather a representation of one after construction.
See the :mod:`.enzyme_module` documentation for more information on
constructing :class:`.EnzymeModule`\ s.
"""
from collections import OrderedDict
from copy import copy, deepcopy

import numpy as np
import pandas as pd
from cobra.core.dictlist import DictList
from cobra.core.group import Group
from six import iteritems, iterkeys, itervalues

from mass.util.dict_with_id import OrderedDictWithID
from mass.util.matrix import _get_matrix_constructor, convert_matrix, matrix_rank
from mass.util.util import _mk_new_dictlist


class EnzymeModuleDict(OrderedDictWithID):
    """Container to store :class:`.EnzymeModule` information after merging.

    Parameters
    ----------
    enzyme : EnzymeModule or EnzymeModuleDict
        The :class:`.EnzymeModule` to be converted into an
        :class:`EnzymeModuleDict`, or an existing :class:`EnzymeModuleDict`.
        If an existing :class:`EnzymeModuleDict` is provided, a new
        :class:`EnzymeModuleDict` is instantiated with the same information
        as the original.

    """

    def __init__(self, id_or_enzyme=None):
        """Initialize EnzymeModuleDict."""
        # Instiantiate new EnzymeModuleDict object if given an EnzymeModuleDict
        if isinstance(id_or_enzyme, (EnzymeModuleDict, dict)):
            super(EnzymeModuleDict, self).__init__(
                id_or_enzyme["id"], data_dict=dict(id_or_enzyme)
            )
        # Initialize an EnzymeModuleDict using an EnzymeModule
        elif id_or_enzyme.__class__.__name__ == "EnzymeModule":
            super(EnzymeModuleDict, self).__init__(id_or_enzyme.id)
            items = {"enzyme_concentration_total_equation": None}
            items.update(id_or_enzyme.__dict__)
            for key, value in iteritems(items):
                nkey = key.lstrip("_")
                if nkey not in iterkeys(_ORDERED_ENZYMEMODULE_DICT_DEFAULTS):
                    continue
                elif nkey == "S":
                    self[nkey] = id_or_enzyme._mk_stoich_matrix(
                        array_type="DataFrame", update_model=False
                    )
                elif "_equation" in nkey:
                    self[nkey] = getattr(id_or_enzyme, nkey, None)
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
        """Copy an :class:`EnzymeModuleDict`."""
        # No references to the MassModel when copying the EnzymeModuleDict
        model = self.model
        setattr(self, "model", None)
        for attr in [
            "enzyme_module_ligands",
            "enzyme_module_forms",
            "enzyme_module_reactions",
        ]:
            for i in getattr(self, attr):
                setattr(i, "_model", None)
            for dictlist in itervalues(self[attr + "_categorized"]):
                for i in dictlist:
                    setattr(i, "_model", None)

        # The EnzymeModuleDict can now be copied
        enzyme_dict_copy = deepcopy(self)
        # Restore references for the original EnzymeModuleDict
        setattr(self, "model", model)
        for attr in [
            "enzyme_module_ligands",
            "enzyme_module_forms",
            "enzyme_module_reactions",
        ]:
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
        for attr in [
            "enzyme_module_ligands",
            "enzyme_module_forms",
            "enzyme_module_reactions",
        ]:
            model_dictlist = {
                "enzyme_module_ligands": model.metabolites,
                "enzyme_module_forms": model.metabolites,
                "enzyme_module_reactions": model.reactions,
            }.get(attr)
            attr_value = getattr(self, attr)
            self[attr] = _mk_new_dictlist(model_dictlist, attr_value)
            attr += "_categorized"
            attr_value = getattr(self, attr)
            if isinstance(attr_value, dict):
                attr_value = [
                    Group(category, members=model_dictlist.get_by_any(obj_ids))
                    for category, obj_ids in iteritems(attr_value)
                ]
                model.add_groups(attr_value)
            self[attr] = _mk_new_dictlist(model.groups, attr_value)

        for enzyme_module_form in self.enzyme_module_forms:
            enzyme_module_form._repair_bound_obj_pointers()

    def _make_enzyme_stoichiometric_matrix(self, update=False):
        """Return the S matrix based on enzyme forms, reactions, and ligands.

        Warnings
        --------
        This method is intended for internal use only.

        """
        # Set up for matrix construction.
        (matrix_constructor, array_type, dtype) = _get_matrix_constructor(
            array_type="DataFrame", dtype=np.float_
        )

        metabolites = DictList(
            [
                met
                for attr in ["enzyme_module_ligands", "enzyme_module_forms"]
                for met in self[attr]
            ]
        )

        stoich_mat = matrix_constructor(
            (len(metabolites), len(self.enzyme_module_reactions))
        )
        # Get the indicies for the forms and reactions
        m_ind = metabolites.index
        r_ind = self.enzyme_module_reactions.index

        # Build the matrix
        for rxn in self.enzyme_module_reactions:
            for met, stoich in iteritems(rxn.metabolites):
                stoich_mat[m_ind(met), r_ind(rxn)] = stoich

        # Convert the matrix to the desired type
        stoich_mat = convert_matrix(
            stoich_mat,
            array_type=array_type,
            dtype=dtype,
            row_ids=[m.id for m in metabolites],
            col_ids=[r.id for r in self.enzyme_module_reactions],
        )

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
        """HTML representation of the overview for the EnzymeModuleDict.

        Warnings
        --------
        This method is intended for internal use only.

        """
        try:
            dim_S = "{0}x{1}".format(self.S.shape[0], self.S.shape[1])
            rank = matrix_rank(self.S)
        except (np.linalg.LinAlgError, ValueError):
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
        """.format(
            name=self.id,
            address="0x0%x" % id(self),
            dim_stoich_mat=dim_S,
            mat_rank=rank,
            subsystem=self.subsystem,
            num_enzyme_module_ligands=len(self.enzyme_module_ligands),
            num_enz_forms=len(self.enzyme_module_forms),
            num_enz_reactions=len(self.enzyme_module_reactions),
            enz_conc=self.enzyme_concentration_total,
            enz_flux=self.enzyme_rate,
        )

    # Dunders
    def __getattr__(self, name):
        """Override of default ``getattr`` implementation.

        Warnings
        --------
        This method is intended for internal use only.

        """
        if name in ["_id", "_model"]:
            name = name.lstrip("_")
        if name not in self:
            raise AttributeError("Attribute '" + name + "' does not exist.")

        return self[name]

    def __setattr__(self, name, value):
        """Override of default ``setattr`` implementation.

        Warnings
        --------
        This method is intended for internal use only.

        """
        if name in ["_id", "_model"]:
            name = name.lstrip("_")
        self[name] = value

    def __delattr__(self, name):
        """Override of default ``delattr`` implementation.

        Warnings
        --------
        This method is intended for internal use only.

        """
        if name in self:
            del self[name]
        else:
            raise AttributeError("Attribute '" + name + "' does not exist.")

    def __dir__(self):
        """Override of default ``dir`` implementation to include the keys.

        Warnings
        --------
        This method is intended for internal use only.

        """
        return sorted(
            set(list(iterkeys(self)) + super(EnzymeModuleDict, self).__dir__())
        )

    def __copy__(self):
        """Create a copy of the EnzymeModuleDict.

        Warnings
        --------
        This method is intended for internal use only.

        """
        return copy(super(EnzymeModuleDict, self))

    def __deepcopy__(self, memo):
        """Create a deepcopy of the EnzymeModuleDict.

        Warnings
        --------
        This method is intended for internal use only.

        """
        return deepcopy(super(EnzymeModuleDict, self), memo)


_ORDERED_ENZYMEMODULE_DICT_DEFAULTS = OrderedDict(
    {
        "id": None,
        "name": "",
        "subsystem": "",
        "enzyme_module_ligands": DictList(),
        "enzyme_module_forms": DictList(),
        "enzyme_module_reactions": DictList(),
        "enzyme_module_ligands_categorized": DictList(),
        "enzyme_module_forms_categorized": DictList(),
        "enzyme_module_reactions_categorized": DictList(),
        "enzyme_concentration_total": None,
        "enzyme_rate": None,
        "enzyme_concentration_total_equation": None,
        "enzyme_rate_equation": None,
        "S": pd.DataFrame(),
        "model": None,
    }
)
