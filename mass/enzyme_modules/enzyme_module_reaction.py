# -*- coding: utf-8 -*-
"""TODO Module Docstrings."""
from __future__ import absolute_import

from collections import defaultdict

from mass.core.massreaction import MassReaction
from mass.enzyme_modules.enzyme_module_form import EnzymeModuleForm


class EnzymeModuleReaction(MassReaction):
    """Class for representing a reaction containing an EnzymeModuleForm.

    Parameters
    ----------
    id_or_reaction: str, MassReaction, EnzymeModuleReaction
        Either a string identifier to associate with the EnzymeModuleReaction,
        or an existing EnzymeModuleReaction object. If an existing
        EnzymeModuleReaction or MassReaction object is provided, a new
        EnzymeModuleReaction object is instantiated with the same properties as
        the original object.
    name: str, optional
        A human readable name for the reaction.
    subsystem: str, optional
        The subsystem where the reaction is meant to occur.
    reversible: bool, optional
        The kinetic reversibility of the reaction. Irreversible reactions have
        an equilibrium constant of infinity and a reverse rate constant of 0.
        If not provided, the reaction is assumed to be reversible.

    Attributes
    ----------
    enzyme_module_id: str, optional
        The identifier of the EnzymeModule represented by the EnzymeModuleForm.
    steady_state_flux: float, optional
        The stored (typically steady state) flux for the reaction. Stored flux
        values can be accessed for operations such as PERC calculations.

    """

    def __init__(self, id_or_reaction=None, name="", subsystem="",
                 reversible=True, steady_state_flux=None, enzyme_module_id=""):
        """Initialize the MassReaction Object."""
        super(EnzymeModuleReaction, self).__init__(
            id_or_reaction=id_or_reaction, name=name, subsystem=subsystem,
            reversible=reversible, steady_state_flux=steady_state_flux)
        if isinstance(id_or_reaction, EnzymeModuleReaction):
            self.__dict__.update(id_or_reaction.__dict__)
        else:
            # Set the id of the enzyme represented by the EnzymeModuleReaction
            self.enzyme_module_id = enzyme_module_id

    def generate_enzyme_module_reaction_name(self, update_enzyme=False):
        """Generate a name for the EnzymeModuleReaction based on bound ligands.

        The name is generated based on, bound_catalytic and bound_effector
        attributes of the EnzymeModuleForms involved in the enzyme module
        reactions.

        Parameters
        ----------
        update_enzyme: bool, optional
            If True, update the name attribute of the enzyme form in
            addition to returning the automatically generated name of the
            enzyme form as a str. Default is False.

        Returns
        -------
        name: A str representing the name of the enzyme reaction.

        """
        name = ""
        items = defaultdict(list)
        for met in self.metabolites:
            key = "Enz" if isinstance(met, EnzymeModuleForm) else "Lig"
            key += " React" if met in self.reactants else " Prod"
            items[key].append(met)

        for attr in ["bound_catalytic", "bound_effectors"]:
            for enz_r, enz_p in zip(items["Enz React"], items["Enz Prod"]):
                r_dict, p_dict = (getattr(enz_r, attr), getattr(enz_p, attr))
                diff = {}
                for k in list(set(p_dict).union(set(r_dict))):
                    if k in p_dict and k in r_dict:
                        coeff = abs(p_dict[k] - r_dict[k])
                    elif k in p_dict or k in r_dict:
                        coeff = [d[k] for d in [r_dict, p_dict]
                                 if k in d].pop()
                    if coeff != 0:
                        diff[k] = coeff

                if diff:
                    if list(diff) == list(items["Lig React"]):
                        mets = "-".join([m._remove_compartment_from_id_str()
                                         for m in [enz_r] + list(diff)])
                        action = " binding"
                    elif list(diff) == list(items["Lig Prod"]):
                        mets = "-".join([m._remove_compartment_from_id_str()
                                         for m in [enz_r] + list(diff)])
                        action = " release"
                    else:
                        mets = enz_r._remove_compartment_from_id_str()
                        action = " catalyzation"
                    name = mets + action

                if not name:
                    name = "-".join([enz_form._remove_compartment_from_id_str()
                                     for enz_form in [enz_r, enz_p]])
                    name += " transition"

        # Assign the new name to the name attribute
        if update_enzyme:
            self.name = name

        return name
