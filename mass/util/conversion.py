# -*- coding: utf-8 -*-

# Compatibility with Python 2.7
from __future__ import absolute_import

# Import necesary packages

# from cobra
from cobra.core.metabolite import Metabolite
from cobra.core.reaction import Reaction

# from mass
from mass.core.massmetabolite import MassMetabolite
from mass.core.massreaction import MassReaction

def to_cobra_metabolite(MassMetab=None, cobra_id=None):
	"""To create a cobra Metabolite object from a MassMetabolite Object

	Parameters
	----------
	cobra_id : string or None
		id for the cobra Metabolite object. If no id is specified,
		one will automatically be generated with the
		Metabolite object's current ID + _mass at the end

	Warnings
	--------
	All other fields will initialize to default values.
	"""
	if MassMetab == None:
		warn("No MassMetabolite Object was given")
		return None

	if not isinstance(MassMetab, Metabolite):
		raise TypeError("Metabolite must be a  MassMetabolite Object")

	if cobra_id == None:
		cobra_id = MassMetab.id + "_cobra"

	new_cobra_metab = Metabolite(id=cobra_id, name=MassMetab.name,
							formula=MassMetab._formula,
                            charge=MassMetab._charge,
							compartment=MassMetab._compartment)

	return new_cobra_metab

def from_cobra_metabolite(CobraMetab=None, mass_id=None):
	"""To create a MassMetabolite object from a cobra Metabolite Object

	Parameters
	----------
	CobraMetabolite : Metabolite Object from cobra.core.metabolite
		The cobra Metabolite Object for creating the MassMetabolite object
	mass_id : string or None
		id for the MassMetabolite object. If no id is specified,
		one will automatically be generated with the
		MassMetabolite object's current ID + _mass at the end

	Warnings
	--------
	All other fields will initialize to default values.

	A Metabolite Object from cobra.core.metabolite must be imported into
		the enviroment in order for this method to work properly.
	"""

	if CobraMetab == None:
		warn("No cobra Metabolite Object was given")
		return None

	if not isinstance(CobraMetab, Metabolite):
		raise TypeError("Metabolite must be a cobra Metabolite Object")

	if mass_id == None:
		mass_id = CobraMetab.id + "_mass"

	new_mass_metab = MassMetabolite(id=mass_id, name=CobraMetab.name,
								formula=CobraMetab.formula,
								charge=CobraMetab.charge,
								compartment=CobraMetab.compartment)
	return new_mass_metab

def to_cobra_reaction(MassRxn=None, cobra_id=None,
                        lower_bound=None,
                        upper_bound=None):
    """To convert a MassReaction object into a cobra Reaction object

    If the lower and/or upper bounds are not specified,the reversibility
    will be used to determine the bounds for initializing the reaction.

    For reversible MassReaction objects:
        lower_bound =-1000 and/or upper_bound=1000
    For irreversible MassReaction objects:
        lower_bound =0 and upper_bound=1000

    Warnings
    --------
    All other fields in a cobra Reaction will initialize to their defaults
    """
    if MassRxn == None:
        warn("No MassReaction Object was given")
        return None

    if not isinstance(MassRxn, MassReaction):
        raise TypeError("Reaction must be a MassReaction Object")

    if cobra_id == None:
        cobra_id = MassRxn._id + "_cobra"

    if lower_bound == None:
        if MassRxn._reversibility == True:
            lb = -1000
        else:
            lb = 0.

    if upper_bound == None:
        ub = 1000

    cobra_rxn = Reaction(id=cobra_id, name=MassRxn.name,
                        subsystem=MassRxn.subsystem, lower_bound=lb,
                        upper_bound=ub, objective_coefficient=0.)

    cobra_metab_dict = {to_cobra_metabolite(metab) : coefficient
                    for metab, coefficient in iteritems(MassRxn._metabolites)}

    cobra_rxn.add_metabolites(cobra_metab_dict)
    cobra_rxn._genes = MassRxn._genes
    cobra_rxn._gene_reaction_rule = MassRxn._gene_reaction_rule
    cobra_rxn._update_awareness()

    return cobra_rxn

def from_cobra_reaction(CobraRxn=None, mass_id=None,
                        kinetic_reversibility=None):
    """To convert a cobra Reaction object into a MassReaction object

    If kinetic_reversibility is not specifiied, will try to infer
    reversibility from the upper and lower bounds of the cobra object.

    Warnings
    --------
    All other fields in a MassReaction will initialize to their defaults

    A Reaction Object from cobra.core.reaction must be imported into
        the enviroment in order for this method to work properly.
    """
    if CobraRxn == None:
        warn("No cobra Reaction Object was given")
        return None

    if not isinstance(CobraRxn, Reaction):
        raise TypeError("Reaction must be a cobra Reaction Object")

    if mass_id == None:
        mass_id = CobraRxn._id + "_mass"

    if kinetic_reversibility == None:
        kinetic_reversibility = CobraRxn.reversibility

    mass_rxn = MassReaction(id=mass_id, name=CobraRxn.name,
                        subsystem=CobraRxn.subsystem,
                        reversibility=kinetic_reversibility)

    mass_metab_dict = {from_cobra_metabolite(cobra_metab): coeff
        for cobra_metab, coeff in iteritems(CobraRxn._metabolites)
        }

    mass_rxn.add_metabolites(mass_metab_dict)
    mass_rxn._genes = CobraRxn._genes
    mass_rxn._gene_reaction_rule = CobraRxn._gene_reaction_rule
    mass_rxn._rate_law = mass_rxn.generate_rate_law()
    mass_rxn._rate_law_expr = mass_rxn.generate_rate_law_expr()
    mass_rxn._update_awareness()

    return mass_rxn
