# -*- coding: utf-8 -*-
r"""
EnzymeModuleReaction is a class for holding information regarding enzymatic reactions.

The :class:`EnzymeModuleReaction` class inherits and extends the
:class:`~.MassReaction` class. It is designed to represent the reactions
and transitions involving :class:`~.EnzymeModuleSpecies`\ s represented in the
:class:`~.EnzymeModule` class.

The enzyme specific attributes on the :class:`EnzymeModuleReaction` are the
following:

    * :attr:`~EnzymeModuleReaction.enzyme_module_id`

Some other important points about the :class:`EnzymeModuleReaction` include:

    * If the :attr:`name` attribute is not set upon initializing, it is
      automatically generated using the enzyme specific attributes of the
      associated :class:`~.EnzymeModuleSpecies`\ s.
    * Even though :class:`~.MassReaction`\ s are catalyzed by enzymes, an
      enzymatic reaction in the context of this module will refer to reactions
      that involve :class:`~.EnzymeModuleSpecies`\ (s) and are associated with
      an :class:`~.EnzymeModule`.

""" # noqa
from collections import defaultdict

from mass.core.mass_reaction import MassReaction
from mass.enzyme_modules.enzyme_module_species import EnzymeModuleSpecies


class EnzymeModuleReaction(MassReaction):
    """Class representing enzymatic reaction in an :class:`~.EnzymeModule`.

    Parameters
    ----------
    id_or_reaction: str, MassReaction, EnzymeModuleReaction
        A string identifier to associate with the enzymatic reaction, or an
        existing reaction object. If an existing reaction object is
        provided, a new :class:`EnzymeModuleReaction` is instantiated with the
        same properties as the original reaction.
    name : str
        A human readable name for the enzymatic reaction.
    subsystem : str
        The subsystem where the enzymatic reaction is meant to occur.
    reversible : bool
        The kinetic reversibility of the reaction. Irreversible reactions have
        an equilibrium constant and a reverse rate constant as set in the
        :attr:`~.MassBaseConfiguration.irreversible_Keq` and
        :attr:`~.MassBaseConfiguration.irreversible_kr` attributes of the
        :class:`~.MassConfiguration`. Default is ``True``.

    Attributes
    ----------
    enzyme_module_id : str
        The identifier of the associated :class:`~.EnzymeModule`.

    """

    def __init__(self, id_or_reaction=None, name="", subsystem="",
                 reversible=True, steady_state_flux=None, enzyme_module_id=""):
        """Initialize the EnzymeModuleReaction."""
        super(EnzymeModuleReaction, self).__init__(
            id_or_reaction=id_or_reaction, name=name, subsystem=subsystem,
            reversible=reversible, steady_state_flux=steady_state_flux)
        if isinstance(id_or_reaction, EnzymeModuleReaction):
            self.__dict__.update(id_or_reaction.__dict__)
        else:
            # Set the id of the enzyme represented by the EnzymeModuleReaction
            self.enzyme_module_id = enzyme_module_id

    def generate_enzyme_module_reaction_name(self, update_enzyme=False):
        r"""Generate a name for the enzymatic reaction based on bound ligands.

        Notes
        -----
        * The :attr:`~.EnzymeModuleSpecies.bound_catalytic` and
          :attr:`~.EnzymeModuleSpecies.bound_effectors` attributes of the
          associated :class:`~.EnzymeModuleSpecies` are used in generating
          the name.

        Parameters
        ----------
        update_enzyme : bool
            If ``True``, update the :attr:`name` attribute of the enzymatic
            reaction in addition to returning the generated name.
            Default is ``False``.

        Returns
        -------
        str
            String representing the name of the :class:`EnzymeModuleReaction`.

        """
        name = ""
        items = defaultdict(list)
        for met in self.metabolites:
            key = "Enz" if isinstance(met, EnzymeModuleSpecies) else "Lig"
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


__all__ = ("EnzymeModuleReaction",)
