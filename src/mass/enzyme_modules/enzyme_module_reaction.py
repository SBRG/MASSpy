# -*- coding: utf-8 -*-
r"""
EnzymeModuleReaction is a class for holding information regarding enzyme module reactions.

The :class:`EnzymeModuleReaction` class inherits and extends the
:class:`~.MassReaction` class. It is designed to represent the reactions
and transitions involving :class:`~.EnzymeModuleForm`\ s represented in the
:class:`~.EnzymeModule` class.

The enzyme specific attributes on the :class:`EnzymeModuleReaction` are the
following:

    * :attr:`~EnzymeModuleReaction.enzyme_module_id`

Some other important points about the :class:`EnzymeModuleReaction` include:

    * If the :attr:`name` attribute is not set upon initializing, it is
      automatically generated using the enzyme specific attributes of the
      associated :class:`~.EnzymeModuleForm`\ s.

    * Even though :class:`~.MassReaction`\ s are also catalyzed by enzymes, an
      enzyme module reaction in the context of this module will refer to
      reactions that involve :class:`~.EnzymeModuleForm`\ (s) and are
      associated with an :class:`~.EnzymeModule`.

"""
from collections import defaultdict

from mass.core.mass_reaction import MassReaction
from mass.enzyme_modules.enzyme_module_form import EnzymeModuleForm


class EnzymeModuleReaction(MassReaction):
    """Class representing an enzyme reaction in an :class:`~.EnzymeModule`.

    Accepted ``kwargs`` are passed to the initialization method of the base
    class, :class:`.MassReaction`.

    Parameters
    ----------
    id_or_reaction : str, MassReaction, EnzymeModuleReaction
        A string identifier to associate with the enzyme module reaction,
        or an existing reaction object. If an existing reaction object is
        provided, a new :class:`EnzymeModuleReaction` is instantiated with the
        same properties as the original reaction.
    enzyme_module_id : str
        The identifier of the associated :class:`~.EnzymeModule`.
    **kwargs
        name :
            ``str`` representing a human readable name for the enzyme
            module reaction.
        subsystem :
            ``str`` representing the subsystem where the enzyme
            module reaction is meant to occur.
        reversible :
            ``bool`` indicating the the kinetic reversibility of the reaction.
            Irreversible reactions have an equilibrium constant and a
            reverse rate constant as set in the
            :attr:`~.MassBaseConfiguration.irreversible_Keq` and
            :attr:`~.MassBaseConfiguration.irreversible_kr` attributes of the
            :class:`~.MassConfiguration`. Default is ``True``.
        steady_state_flux :
            ``float`` representing the stored (typically steady state) flux
            for the reaction.

    """

    def __init__(self, id_or_reaction=None, enzyme_module_id="", **kwargs):
        """Initialize the EnzymeModuleReaction."""
        super(EnzymeModuleReaction, self).__init__(
            id_or_reaction=id_or_reaction,
            name=kwargs.get("name", ""),
            subsystem=kwargs.get("subsystem", ""),
            reversible=kwargs.get("reversible", True),
            steady_state_flux=kwargs.get("steady_state_flux", None),
        )
        if isinstance(id_or_reaction, EnzymeModuleReaction):
            self.__dict__.update(id_or_reaction.__dict__)
        else:
            # Set the id of the enzyme represented by the EnzymeModuleReaction
            self.enzyme_module_id = enzyme_module_id

    def generate_enzyme_module_reaction_name(self, update_enzyme=False):
        r"""Generate name for an enzyme module reaction based on bound ligands.

        Notes
        -----
        * The :attr:`~.EnzymeModuleForm.bound_metabolites` attribute of the
          associated :class:`~.EnzymeModuleForm` is used in generating
          the name.

        Parameters
        ----------
        update_enzyme : bool
            If ``True``, update the :attr:`name` attribute of the enzyme
            module reaction in addition to returning the generated name.
            Default is ``False``.

        Returns
        -------
        str
            String representing the name of the :class:`EnzymeModuleReaction`.

        """
        name = ""
        items = defaultdict(list)
        for met in self.metabolites:
            key = "Enz" if isinstance(met, EnzymeModuleForm) else "Lig"
            key += " React" if met in self.reactants else " Prod"
            items[key].append(met)

        for enz_r, enz_p in zip(items["Enz React"], items["Enz Prod"]):
            r_dict, p_dict = (
                getattr(enz_r, "bound_metabolites"),
                getattr(enz_p, "bound_metabolites"),
            )
            diff = {}
            for key in list(set(p_dict).union(set(r_dict))):
                if key in p_dict and key in r_dict:
                    coeff = abs(p_dict[key] - r_dict[key])
                elif key in p_dict or key in r_dict:
                    coeff = [d[key] for d in [r_dict, p_dict] if key in d].pop()
                if coeff != 0:
                    diff[key] = coeff

            if diff:
                if list(diff) != list(items["Lig React"]) and list(diff) != list(
                    items["Lig Prod"]
                ):
                    name_str = enz_r._remove_compartment_from_id_str()
                    name_str += " catalyzation"
                else:
                    name_str = "-".join(
                        [
                            m._remove_compartment_from_id_str()
                            for m in [enz_r] + list(diff)
                        ]
                    )
                    name_str += str(
                        " binding"
                        if list(diff) == list(items["Lig React"])
                        else " release"
                    )
                name = name_str

            if not name:
                name = "-".join(
                    [
                        enz_form._remove_compartment_from_id_str()
                        for enz_form in [enz_r, enz_p]
                    ]
                )
                name += " transition"

        # Assign the new name to the name attribute
        if update_enzyme:
            self.name = name

        return name


__all__ = ("EnzymeModuleReaction",)
