# -*- coding: utf-8 -*-

# Compatibility with Python 2.7
from __future__ import absolute_import

# Import necesary packages
from six import iteritems

# cobra
from cobra.core.metabolite import Metabolite
from cobra.core.model import Model
from cobra.core.reaction import Reaction

# mass
from mass.core import massmetabolite
from mass.core import massmodel
from mass.core import massreaction

# Class begins
## Public
def to_cobra_metabolite(mass_metabolite, cobra_id=None):
	"""To create a cobra Metabolite from a mass MassMetabolite.

	Parameters
	----------
	mass_metabolite : mass.MassMetabolite
		The MassMetabolite for creating the cobra.Metabolite object
	cobra_id : string, optional
		id for the new cobra Metabolite. If no id is specified,
		one will automatically be generated with the
		MassMetabolite object's current id.

	Returns
	-------
	cobra.Metabolite
		The new cobra Metabolite

	Warnings
	--------
	All similar fields will initialize to be identical to the input object.
	All other fields will initialize to default values.
	"""
	# Check the input
	if not isinstance(mass_metabolite, massmetabolite.MassMetabolite):
		raise TypeError("Must be a mass MassMetabolite.")

	# Generate a new ID if none is specified
	if cobra_id is None:
		cobra_id = mass_metabolite.id

	# Generate the cobra Metabolite
	cobra_metab = Metabolite(id=cobra_id, name=mass_metabolite.name,
							formula=mass_metabolite.formula,
							charge=mass_metabolite.charge,
							compartment=mass_metabolite.compartment)

	cobra_metab._constraint_sense = mass_metabolite._constraint_sense
	cobra_metab._bound = mass_metabolite._bound

	return cobra_metab

def to_massmetabolite(cobra_metabolite, mass_id=None):
	"""To create a mass MassMetabolite from a cobra Metabolite.

	Parameters
	----------
	cobra_metabolite : cobra.Metabolite
		The Metabolite for creating the mass.MassMetabolite object
	mass_id : string, optional
		id for the new mass MassMetabolite. If no id is specified,
		one will automatically be generated with the
		Metabolite object's current id.

	Returns
	-------
	mass.MassMetabolite
		The new mass MassMetsbolite object

	Warnings
	--------
	All similar fields will initialize to be identical to the input object.
	All other fields will initialize to default values.
	"""
	# Check the input
	if not isinstance(cobra_metabolite, Metabolite):
		raise TypeError("Must be a cobra Metabolite.")

	# Generate a new ID if none is specified
	if mass_id is None:
		mass_id = cobra_metabolite.id

	# Generate the mass MassMetabolite
	mass_metab = massmetabolite.MassMetabolite(id=mass_id,
				name=cobra_metabolite.name, formula=cobra_metabolite.formula,
				charge=cobra_metabolite.charge,
				compartment=cobra_metabolite.compartment)

	mass_metab._constraint_sense = cobra_metabolite._constraint_sense
	mass_metab._bound = cobra_metabolite._bound

	return mass_metab

def to_cobra_reaction(mass_reaction, cobra_id=None,
					upper_bound=None, lower_bound=None):
	"""To create a cobra Reaction from a mass MassReaction.

	If the lower and/or upper bounds are not specified, the reversibility
	will be used to determine the bounds for initializing the reaction.

	For reversible MassReaction objects:
		upper_bound=1000, lower_bound=-1000,
	For irreversible MassReaction objects:
		lower_bound=0, upper_bound=1000

	Parameters
	----------
	mass_reaction : mass.MassReacion
		The MassReaction for creating the cobra.Reaction object
	cobra_id : string, optional
		id for the new cobra Reaction. If no id is specified,
		one will automatically be generated with the
		MassReaction object's current id.
	lower_bound : float, optional
		The initialized lower bound of the cobra Reaction.
	upper_bound : float, optional
		The initialized upper bound of the cobra Reaction

	Returns
	-------
	cobra.Reaction
		The new cobra Reaction object

	Warnings
	--------
	All similar fields will initialize to be identical to the input object.
	All other fields will initialize to default values.
	"""
	# Check the input
	if not isinstance(mass_reaction, massreaction.MassReaction):
		raise TypeError("Must be a mass MassReaction.")

	# Generate a new ID if none is specified
	if cobra_id is None:
		cobra_id = mass_reaction.id

	# Generate the bounds
	if upper_bound is None:
		ub = mass_reaction.upper_bound
	if lower_bound is None:
		lb = mass_reaction.lower_bound

	# Generate the cobra Reaction
	cobra_rxn = Reaction(id=cobra_id, name=mass_reaction.name,
						subsystem=mass_reaction.subsystem, lower_bound=lb,
						upper_bound=ub, objective_coefficient=0.)
	cobra_rxn.variable_kind = mass_reaction.variable_kind
	# Generate and add cobra Metabolites
	cobra_metabs = {to_cobra_metabolite(metab) : coefficient
					for metab, coefficient in \
					iteritems(mass_reaction._metabolites)}
	cobra_rxn.add_metabolites(cobra_metabs)

	# Generate and add new genes
	for gene in mass_reaction._genes:
		cobra_rxn._genes.add(gene.copy())
	# Add the gene reaction rule
	cobra_rxn._gene_reaction_rule = mass_reaction._gene_reaction_rule
	# Make new metaboltites and genes aware they are involved in this reaction
	cobra_rxn._update_awareness()
	return cobra_rxn

def to_massreaction(cobra_reaction, mass_id=None,
					kinetic_reversibility=None):
	"""To create a MassReaction from a cobra Reaction.

	If kinetic_reversibility is not specifiied, will try to infer
	reversibility from the upper and lower bounds of the cobra object.

	Reversible if lower_bound < 0 < upper_bound. Otherwise irreversible

	Parameters
	----------
	cobra_reaction : cobra.Reaction
		The cobra Reaction for creating the mass.MassReaction object
	mass_id : string, optional
		id for the new mass MassReaction. If no id is specified,
		one will automatically be generated with the
		Reaction object's current id.
	kinetic_reversibility : bool, optional
		The reversibility of the mass MassReaction.

	Returns
	-------
	mass.MassReaction
		The new mass MassReaction object

	Warnings
	--------
	All similar fields will initialize to be identical to the input object.
	All other fields will initialize to default values.
	"""
	# Check the input
	if not isinstance(cobra_reaction, Reaction):
		raise TypeError("Must be a cobra Reaction.")

	# Generate a new ID if none is specified
	if mass_id is None:
		mass_id = cobra_reaction.id

	# Infer kinetic reversibility from bounds if none is specified
	if kinetic_reversibility is None:
		kinetic_reversibility = cobra_reaction.reversibility

	# Generate the mass MassReaction
	mass_rxn = massreaction.MassReaction(id=mass_id, name=cobra_reaction.name,
										subsystem=cobra_reaction.subsystem,
										reversible=kinetic_reversibility)
	# Store old cobra information
	mass_rxn.variable_kind = cobra_reaction.variable_kind
	mass_rxn.objective_coefficient = 0.
	mass_rxn.lower_bound = cobra_reaction.lower_bound
	mass_rxn.upper_bound = cobra_reaction.upper_bound

	# Generate and add mass MassMetsabolites
	mass_metabs = {to_massmetabolite(cobra_metab): coeff
			for cobra_metab, coeff in iteritems(cobra_reaction._metabolites)}
	mass_rxn.add_metabolites(mass_metabs)

	# Generate and add new genes
	for gene in cobra_reaction._genes:
		mass_rxn._genes.add(gene.copy())
	# Add the gene reaction rule
	mass_rxn._gene_reaction_rule = cobra_reaction._gene_reaction_rule
	# Make new metaboltites and genes aware they are involved in this reaction
	mass_rxn._update_awareness()
	return mass_rxn

def to_cobra_model(mass_model, cobra_id=None):
	"""To create a cobra Model from a mass MassModel.

	Parameters
	----------
	mass_model : mass.MassModel
		The MassModel for creating the cobra.Model object
	cobra_id : string, optional
		id for the new cobra Model. If no id is specified,
		one will automatically be generated with the
		MassModel object's current id + _cobra.

	Returns
	-------
	cobra.Model
		The new cobra Model object

	Warnings
	--------
	All similar fields will initialize to be identical to the input object.
	All other fields will initialize to default values.
	"""
	# Check the input
	if not isinstance(mass_model, massmodel.MassModel):
		raise TypeError("Must be a mass MassModel.")

	# Generate a new ID if none is specified
	if cobra_id is None:
		cobra_id = mass_model.id + "_mass"

	# Generate the cobra Model
	cobra_model = Model(id_or_model=cobra_id, name=mass_model.name)

	# Add reactions and metabolites
	cobra_rxns = [to_cobra_reaction(rxn) for rxn in mass_model.reactions]
	cobra_model.add_reactions(cobra_rxns)
	# Add compartments
	cobra_model.compartments = mass_model.compartments.copy()
	# Repair the new model
	cobra_model.repair(rebuild_index=True, rebuild_relationships=True)
	return cobra_model

def to_massmodel(cobra_model, mass_id=None):
	"""To create a mass MassModel from a cobra Model.

	Parameters
	----------
	cobra_model : cobra.Model
		The Model for creating the mass.MassModel object
	mass_id : string, optional
		id for the new mass MassModel. If no id is specified,
		one will automatically be generated with the
		Model object's current id + _mass.

	Returns
	-------
	mass.MassModel
		The new mass MassModel object

	Warnings
	--------
	All similar fields will initialize to be identical to the input object.
	All other fields will initialize to default values.
	"""
	# Check the input
	if not isinstance(cobra_model, Model):
		raise TypeError("Must be a cobra Model.")

	# Generate a new ID if none is specified
	if mass_id is None:
		mass_id = cobra_model.id + "_mass"

	# Generate the mass MassModel
	mass_model = massmodel.MassModel(id_or_massmodel=mass_id,
									name=cobra_model.name)

	# Add reactions and metabolites
	mass_rxns = [to_massreaction(rxn) for rxn in cobra_model.reactions]
	mass_model.add_reactions(mass_rxns)

	# Add compartments
	mass_model.compartments = cobra_model.compartments.copy()
	# Repair the new model
	mass_model.repair(rebuild_index=True, rebuild_relationships=True)
	return mass_model
