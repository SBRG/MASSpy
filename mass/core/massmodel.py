# -*- coding: utf-8 -*-

# Compatibility with Python 2.7
from __future__ import absolute_import

# Import necesary packages
import logging
import re
import sympy as sp
import pandas as pd
from six import string_types, iteritems
from copy import copy, deepcopy
from functools import partial
from warnings import warn

# from cobra
from cobra.core.object import Object
from cobra.core.dictlist import DictList
from cobra.util.context import HistoryManager, resettable, get_context

# from mass
from mass.core.massmetabolite import MassMetabolite
from mass.core.massreaction import MassReaction
from mass.util.array import *

# Class begins
## Set the logger
LOGGER = logging.getLogger(__name__)

# Class definition
class MassModel(Object):
	"""MassModel is a class for storing MassReactions, MassMetabolites, and
	all other information necessary to create a mass model for simulation.

	Parameters
	----------
	id_or_massmodel : MassModel, string
		Either an existing MassModel object in which case a new MassModel
		object is instantiated with the same properties as the original
		massmodel, or an identifier to associate with the massmodel as a string

	Attributes
	----------
	reactions : Dictlist
		A DictList where the keys are the reaction identifiers and the values
		are the associated MassReaction objects
	metabolites : DictList
		A DictList where the keys are the metabolite identifiers and the values
		are the associated MassMetabolite objects
	genes : DictList
		A DictList where the keys are the gene identifiers and the values
		are the associated Gene objects
	initial_conditions : dict
		A dictionary to store the initial conditions for the metabolites,
		where keys are metabolite objects and values are initial conditions.
		Can have different initial conditions from the metabolites if desired.
	custom_rates : dict
		A dictionary to store custom rates for specific reactions, where
		keys are reaction objects and values are the custom rate expressions.
		Custom rates will always have preference over rate laws in reactions.
	compartments : dict
		A dictionary to store the compartment shorthands and their full names.
		Keys are the shorthands and values are the full names.
		Example: {'c': 'cytosol'}
	"""
	def __init__(self, id_or_massmodel=None, name=None,
				matrix_type=None, dtype=None):
		"""Initialize the MassModel Object"""
		if isinstance(id_or_massmodel, MassModel):
			Object.__init__(self, id_or_massmodel, name=name)
			self.__setstate__(id_or_massmodel.__dict__)
			if not hasattr(self, "name"):
				self.name = None
		else:
			Object.__init__(self, id_or_massmodel, name=name)
			# A DictList of MassReactions
			self.reactions = DictList()
			# A DictList of MassMetabolites
			self.metabolites = DictList()
			# A DictList of cobra Genes
			self.genes = DictList()
			# A dictionary of initial conditions for MassMetabolites
			self.initial_conditions = dict()
			#For storing of custom rate laws
			self.custom_rates= dict()
			# A dictionary of the compartments in the model
			self.compartments = dict()

			# Internal storage of S matrix and data types for updating S
			self._matrix_type = matrix_type
			self._dtype = dtype

			# For storing the stoichiometric matrix.
			self._S = create_stoichiometric_matrix(self,
                                    matrix_type=self._matrix_type,
                                    dtype=self._dtype,
            						update_model=True)

			# For storing the HistoryManager contexts
			self._contexts = []

	# Properties
	@property
	def S(self):
		"""Get the Stoichiometric Matrix of the MassModel"""
		return update_S(self, reaction_list=None,
					matrix_type=self._matrix_type, dtype=self._dtype,
					update_model=False)

	@property
	def rates(self):
		""""Return the rate laws for the reactions as human readable strings
		in a dictionary where keys are the reaction objects and values are the
		rate laws
		"""
		rate_dict =  {rxn: rxn.rate_law for rxn in self.reactions}
		if self.custom_rates != {}:
			rate_dict.update(self.custom_rates)
		return rate_dict

	@property
	def rate_expressions(self):
		"""Get the rate laws for the reactions as sympy expressions in a
		dictionary where keys are the reaction objects and values are the
		rate law expressions
		"""
		rate_dict =  {rxn: rxn.rate_law_expr for rxn in self.reactions}
		if self.custom_rates != {}:
			rate_dict.update(self.custom_rates)
		return rate_dict

	@property
	def odes(self):
		"""Get the ODEs for the metabolites as sympy expressions where
		keys are the metabolite objects and values are the ODE expressions
		"""
		return {metab: metab.ode for metab in self.metabolites}

	@property
	def exchanges(self):
		"""Get the exchange reactions in the MassModel"""
		return [rxn for rxn in self.reactions if rxn.exchange]

	@property
	def get_metabolite_compartments(self):
		"""Return all metabolites' compartments

		Identical to the method in cobra.core.model
		"""
		return {metab.compartment for metab in self.metabolites
				if metab.compartment is not None}

	@property
	def get_irreversible_reactions(self):
		"""Return a list of all irreversible reactions in the model."""
<<<<<<< HEAD
		return [rxn for rxn in self.reactions if not rxn.reversible]
=======
		return [rxn for rxn in self.reactions if not rxn.reversibility]
>>>>>>> 0dd510179a2c170b1a5a99ba7b85c4fe71b0e12d

	# Methods
	def add_metabolites(self, metabolite_list, add_initial_conditons=False):
		"""Will add a list of metabolites to the MassModel object and add
		the MassMetabolite initial conditions accordingly.

		The change is revereted upon exit when using the MassModel as a context.

		Parameters
		----------
		metabolite_list : list
			A list of MassMetabolite objects to add to the MassModel
		add_initial_conditons : bool
			If True, will also add the initial conditions associated with each
			metabolite to the model. Otherwise just add metabolites without
			their intial conditions
		"""
		# If the metabolite list is not a list
		if not hasattr(metabolite_list, '__iter__'):
			metabolite_list = [metabolite_list]
		metabolite_list = DictList(metabolite_list)
		if len(metabolite_list) == 0:
			return None

		# Check whether a metabolite is a MassMetabolite object. Then check if
		# metabolites already exist in the massmodel, and ignore those that do
		for metab in metabolite_list:
			if not isinstance(metab, MassMetabolite):
				warn("Skipping %s, not a MassMetabolite" % metab)
				metabolite_list.remove(metab)
		existing_metabs = [metab for metab in metabolite_list
							if metab.id in self.metabolites]
		for metab in existing_metabs:
			LOGGER.info("Skipped MassMetabolite %s as it already exists"
						" in MassModel" % metab.id)
		metabolite_list = [metab for metab in metabolite_list
							if metab.id not in self.metabolites]

		# Have metabolites point to model object
		for metab in metabolite_list:
			metab._model = self
		# Add metabolites, and add initial conditions if True
		self.metabolites += metabolite_list
		if add_initial_conditons:
<<<<<<< HEAD
			self.set_initial_conditions(metabolite_list)
=======
			self.set_initial_conditons(metabolite_list)
>>>>>>> 0dd510179a2c170b1a5a99ba7b85c4fe71b0e12d

		context = get_context(self)
		if context:
			context(partial(self.metabolites.__isub__, metabolite_list))
			for metab in metabolite_list:
				context(partial(setattr, metab, '_model', None))

	def remove_metabolites(self, metabolite_list, destructive=False):
		"""Remove a list of metabolites from the MassModel object, and
		its associated initial condition

		The change is reverted upon exit when using the MassModel as a context.

		Parameters
		----------
		metabolite_list : list
			A list of MassMetabolite objects to remove from the MassModel.
		destructive : bool
			If False, then the MassMetabolite and its initial condition are
			removed from the associated MassReactions. If True, also remove
			the associated MassReactions and their rate laws from the MassModel
		"""
		# If the metabolite list is not a list
		if not hasattr(metabolite_list, '__iter__'):
			metabolite_list = [metabolite_list]
		metabolite_list = DictList(metabolite_list)
		if len(metabolite_list) == 0:
			return None

		# Check whether a metabolite already exists in the massmodel, and
		# ignore those that do not
		metabolite_list = [metab for metab in metabolite_list
							if metab.id in self.metabolites]

		# Remove assoication to model
		for metab in metabolite_list:
			metab._model = None
			if not destructive:
				# Remove metabolites from reactions
				for rxn in list(metab._reaction):
					coefficient = rxn._metabolites[metab]
					rxn.subtract_metabolites({metab : coefficient})
			else:
				# Remove associated reactions if destructive
				for rxn in list(metab._reaction):
					rxn.remove_from_model()

		# Remove initial conditions and then metabolites
		self.remove_initial_conditions(metabolite_list)
		self.metabolites -= metabolite_list

		context = get_context(self)
		if context:
			context(partial(self.metabolites.__iadd__, metabolite_list))
			for metab in metabolite_list:
				context(partial(setattr, metab, '_model', self))

<<<<<<< HEAD
	def set_initial_conditions(self, metabolite_list):
=======
	def set_initial_conditons(self, metabolite_list):
>>>>>>> 0dd510179a2c170b1a5a99ba7b85c4fe71b0e12d
		"""Set the initial conditions for a list of metabolites in the model.

		The metabolite must already exist in the model in order to set its
		initial condition. Initial conditions stored in the metabolite are
		accessed and added to the model. Any existing initial condition for
		a metabolite in the model is replaced.

		The change is reverted upon exit when using the MassModel as a context.

		Parameters
		----------
		metabolite_list : list
			A list of MassMetabolite objects
		"""
		# If the metabolite list is not a list
		if not hasattr(metabolite_list, '__iter__'):
			metabolite_list = [metabolite_list]
		metabolite_list = DictList(metabolite_list)
		if len(metabolite_list) == 0:
			return None

		# Check whether a metabolite already exists in the massmodel, and
		# ignore those that do not
		metabolite_list = [metab for metab in metabolite_list
							if metab.id in self.metabolites]
		# Add the initial conditions
		self.update_initial_conditions({metab: metab.ic
										for metab in metabolite_list})

	def remove_initial_conditions(self, metabolite_list):
		"""Remove initial conditions for a list of metabolites in the model.

		The change is reverted upon exit when using the MassModel as a context.

		Parameters
		----------
		metabolite_list : list
			A list of MassMetabolite objects
		"""
		# If the metabolite list is not a list
		if not hasattr(metabolite_list, '__iter__'):
			metabolite_list = [metabolite_list]
		metabolite_list = DictList(metabolite_list)
		if len(metabolite_list) == 0:
			return None

		# Check whether a metabolite already exists in the massmodel, and
		# ignore those that do not
		metabolite_list = [metab for metab in metabolite_list
							if metab.id in self.metabolites]

		# Keep track of existing initial conditions for HistoryManager
		context = get_context(self)
		if context:
			existing_ics = {metab : self.initial_conditions[metab]
							for metab in metabolite_list
							if metab in self.initial_conditions}
		# Remove the initial conditions
		for metab in metabolite_list:
			del self.initial_conditions[metab]

		if context:
			context(partial(self.initial_conditions.update, existing_ics))

	def update_initial_conditions(self, ic_dict, update_metabolites=False):
		"""Update the initial conditions in the model using a dictionary where
		MassMetabolites are keys and the initial conditions are the values.

		The metabolite must already exist in the model in order to set its
		initial condition. If a metabolites initial conditions already exists
		in the model, it is replaced by the new initial condition.

		The change is reverted upon exit when using the MassModel as a context.

		Parameters
		----------
		ic_dict : dict
			A dictionary where MassMetabolites are keys and the
			initial conditions are the values.
		update_metabolites : bool
			If True, will update the initial conditions in the MassMetabolite
			objects as well. Otherwise, only update the model initial conditons
		"""
		if not isinstance(ic_dict, dict):
			raise TypeError("Input must be a dictionary where keys are the "
						"MassMetabolites, and values are initial conditions")
		if len(ic_dict) == 0:
			return None

		# Check whether a metabolite already exists in the massmodel, and
		# ignore those that do not
		ic_dict = {metab : ic for metab, ic in iteritems(ic_dict)
					if metab in self.metabolites and ic is not None}

		# Keep track of existing initial conditions for HistoryManager
		context = get_context(self)
		if context:
			existing_ics = {metab : self.initial_conditions[metab]
					for metab in ic_dict if metab in self.initial_conditions}

		# Update initial conditions
		self.initial_conditions.update(ic_dict)
		if update_metabolites:
			# Update the initial condition stored in the MassMetabolite.
			for metab, ic_value in iteritems(ic_dict):
				metab._initial_condition = ic_value

		if context:
			for key in ic_dict.keys():
				if key not in existing_ics.keys():
					context(partial(self.initial_conditions.pop, key))
				context(partial(self.initial_conditions.update, existing_ics))
			if update_metabolites:
				for metab, ic_value in iteritems(ic_dict):
					context(partial(setattr, metab, '_initial_condition',
									existing_ics[metab]))

	def add_reactions(self, reaction_list, update_stoichiometry=False):
		"""Add MassReactions and their rates to the MassModel.

		MassReactions with identifiers identical to a reaction already in the
		MassModel are ignored.

		The change is reverted upon exit when using the MassModel as a context.

		Similar to the method in cobra.core.model
		Parameters
		----------
		reaction_list : list
			A list of MassReaction objects to add to the MassModel
		"""
		# If the reaction list is not a list
		if not hasattr(reaction_list, '__iter__'):
			reaction_list = [reaction_list]
		reaction_list = DictList(reaction_list)
		if len(reaction_list) == 0:
			return None

		# Check whether a reaction is a MassReaction object. Then check if
		# reactions already exist in the massmodel, and ignore those that do
		for rxn in reaction_list:
			if not isinstance(rxn, MassReaction):
				warn("Skipping %s, not a MassReaction" % rxn)
				reaction_list.remove(rxn)
		existing_rxns = [rxn for rxn in reaction_list
						if rxn.id in self.reactions]
		for rxn in existing_rxns:
			LOGGER.info("Skipped MassReaction %s as it already exists"
						" in MassModel" % rxn.id)
		reaction_list = [rxn for rxn in reaction_list
							if rxn.id not in self.reactions]

		context = get_context(self)

		# Add reactions, and have reactions point to the model
		for rxn in reaction_list:
			rxn._model = self
			# Loop through metabolites in a reaction
			for metab in list(rxn._metabolites.keys()):
				# If metabolite doesn't exist in the model, add it
				# with its associated initial condition
				if metab not in self.metabolites:
					self.add_metabolites(metab, add_initial_conditons=True)
					metab._reaction.add(rxn)
					if context:
						context(partial(metab._reaction.remove, rxn))
				# Otherwise have the reaction point to the metabolite
				# in the model.
				else:
					coefficient = rxn._metabolites.pop(metab)
					model_metab = self.metabolites.get_by_id(metab.id)
					rxn._metabolites[model_metab] = coefficient
					model_metab._reaction.add(rxn)
					if context:
						context(partial(model_metab._reaction.remove, rxn))

			# Loop through genes associated with a reaction
			for gene in list(rxn._genes):
				# If gene is not in the model, add and have it point to model
				if not self.genes.has_id(gene.id):
					self.genes += [gene]
					gene._model = self

					if context:
						context(partial(self.genes.__isub__, [gene]))
						context(partial(setattr, gene, '_model', None))
				# Otherwise, make the gene point to the one in the model
				else:
					model_gene = self.genes.get_by_id(gene.id)
					if model_gene is not gene:
						rxn._dissociate_gene(gene)
						rxn._associate_gene(model_gene)

		# Add reactions to the model
		self.reactions += reaction_list
		if update_stoichiometry:
			update_S(massmodel, reaction_list=reaction_list, update_model=True)

		if context:
			context(partial(self.reactions.__isub__, reaction_list))
			for rxn in reaction_list:
				context(partial(setattr, rxn, '_model', None))
			if update_stoichiometry:
				context(partial(update_S, massmodel, None, None, None, True))


	def remove_reactions(self, reaction_list, remove_orphans=False,
						update_stoichiometry=False):
		"""Remove MassReactions from the MassModel

		The change is reverted upon exit when the MassModel as a context.

		Parameters
		----------
		reaction_list : list
			A list of MassReaction objects to remove from the MassModel.
		remove_orphans : bool
			Remove orphaned genes and MassMetabolites from the
			MassModel as well.
		"""
		# If the reaction list is not a list
		if not hasattr(reaction_list, '__iter__'):
			reaction_list = [reaction_list]
		reaction_list = DictList(reaction_list)
		if len(reaction_list) == 0:
			return None

		# Check whether a reaction already exists in the massmodel, and
		# ignore those that do not
		reaction_list = [rxn for rxn in reaction_list
							if rxn.id in self.reactions]

		context = get_context(self)

		# Remove model association
		for rxn in reaction_list:
			rxn._model = None
			# Remove metabolite association
			for metab in rxn._metabolites:
				if rxn in metab._reaction:
					metab._reaction.remove(rxn)
					# Remove orphaned metabolites and their initial conditions
					if remove_orphans and len(metab._reaction) == 0:
						self.remove_metabolites(metab)
					if context:
						context(partial(metab._reaction.add, rxn))

			# Remove gene association
			for gene in rxn._genes:
				if rxn in gene._reaction:
					gene._reaction.remove(rxn)
					# Remove orphaned genes
					if remove_orphans and len(gene._reaction) == 0:
						self.genes.remove(gene)
						if context:
							context(partial(self.genes.add, gene))
					if context:
						context(partial(gene._reaction.add, rxn))
			if context:
				context(partial(setattr, rxn, '_model', self))

		# Remove reactions from the model
		self.reactions -= reaction_list
		if update_stoichiometry:
			update_S(massmodel, reaction_list=None, update_model=True)

		if context:
			context(partial(self.reactions.__iadd__, reaction_list))
			for rxn in reaction_list:
				context(partial(setattr, rxn, '_model', self))
			if update_stoichiometry:
				context(partial(update_S, massmodel, reaction_list,
							None, None, True))

	def add_exchange(self, metabolite, exchange_type="source",
					reversible=True, update_stoichiometry=False):
		"""Add an exchange reaction for a given metabolite using the
		pre-defined exchange types "source" for into the compartment
		and "sink" for exiting the compartment.

		The change is reverted upon exit when the MassModel as a context.

		Parameters
		----------
		metabolite : MassMetabolite
			Any given metabolite to create an exchange for.
		exchange_type : string, {"sink", "source"}
			The type of exchange reaction to create are not case sensitive.
		reversible : bool
			If True, exchange is reversible. When using a user-defined type,
			must specify the reversiblity
		"""
		# Check whether metabolite is a MassMetabolite object
		if not isinstance(metabolite, MassMetabolite):
			raise TypeError("metabolite must be a MassMetabolite object")

		type_dict = {"source":1,"sink":-1}

		# Set the type of exchange
		if not isinstance(exchange_type, string_types):
			raise TypeError("rxn_type must be a string")
		else:
			exchange_type = exchange_type.lower()

		if exchange_type in type_dict:
			rxn_id = "{}_{}".format("EX", metabolite.id)
			if rxn_id in self.reactions:
				warn("Reaction %s already exists in model" % rxn_id)
				return None
			c = type_dict[exchange_type]
		else:
			raise TypeError("Exchange type must be either 'source' or 'sink'")
		rxn_name = "{} {}".format(metabolite.name, exchange_type)

		rxn = MassReaction(id=rxn_id, name=rxn_name,
<<<<<<< HEAD
					subsystem="Transport/Exchange",reversible=reversible)
=======
					subsystem="Transport/Exchange",reversibility=reversible)
>>>>>>> 0dd510179a2c170b1a5a99ba7b85c4fe71b0e12d
		rxn.add_metabolites({metabolite: c})
		self.add_reactions([rxn], update_stoichiometry)


	def repair(self, rebuild_index=True, rebuild_relationships=True):
		"""Update all indexes and pointers in a MassModel

		Identical to the method in cobra.core.model

		Parameters
		----------
		rebuild_index : bool
			If True, rebuild the indecies kept in reactions,
			metabolites and genes.
		rebuild_relationships : bool
			If True, reset all associations between the reactions, metabolites,
			genes, and the MassModel and re-add them
		"""
		if not isinstance(rebuild_index, bool) or \
			not isinstance(rebuild_relationships, bool):
			raise TypeError("rebuild_index and rebuild_relationships "
							"must be True or False")
		# Rebuild DictList indicies
		if rebuild_index:
			self.reactions._generate_index()
			self.metabolites._generate_index()
			self.genes._generate_index()
		# Rebuild relationships between reactions and their associated
		# genes and metabolites
		if rebuild_relationships:
			for metab in self.metabolites:
				metab._reaction.clear()
			for gene in self.genes:
				gene._reaction.clear()
			for rxn in self.reactions:
				for metab in rxn.metabolites:
					metab._reaction.add(rxn)
				for gene in rxn.genes:
					gene._reaction.add(rxn)
		# Make all objects point to model
		for dictlist in (self.reactions, self.genes, self.metabolites):
			for item in dictlist:
				item._model = self

	def copy(self):
		"""Provides a partial 'deepcopy' of the MassModel. All of
		the MassMetabolite, MassReaction and Gene objects, the
		initial conditions, custom rates, and the S matrix are created anew
		but in a faster fashion than deepcopy
		"""
		# Define a new model
		new_model = self.__class__()
		# Define items to not copy by their references
		do_not_copy_by_ref = {"metabolites", "reactions", "genes",
							"initial_conditions","custom_rates", "_S"
							"notes", "annotations"}
		for attr in self.__dict__:
			if attr not in do_not_copy_by_ref:
				new_model.__dict__[attr] = self.__dict__[attr]
		new_model.notes = deepcopy(self.notes)
		new_model.annotation = deepcopy(self.annotation)

		# Copy metabolites
		new_model.metabolites = DictList()
		do_not_copy_by_ref = {"_reaction", "_model"}
		for metab in self.metabolites:
			new_metab = metab.__class__()
			for attr, value in iteritems(metab.__dict__):
				if attr not in do_not_copy_by_ref:
					new_metab.__dict__[attr] = copy(
							value) if attr == "formula" else value
			new_metab._model = new_model
			new_model.metabolites.append(new_metab)
			# Copy the initial condition
			if metab in self.initial_conditions.keys():
				ic = self.initial_conditions[metab]
				new_model.initial_conditions[new_metab] = ic

		# Copy the genes
		new_model.genes = DictList()
		for gene in self.genes:
			new_gene = gene.__class__(None)
			for attr, value in iteritems(gene.__dict__):
				if attr not in do_not_copy_by_ref:
					new_gene.__dict__[attr] = copy(
							value) if attr == "formula" else value
			new_gene._model = new_model
			new_model.genes.append(new_gene)

		# Copy the reactions
		new_model.reactions = DictList()
		do_not_copy_by_ref = {"_model", "_metabolites", "_genes"}
		for rxn in self.reactions:
			new_rxn = rxn.__class__()
			for attr, value in iteritems(rxn.__dict__):
				if attr not in do_not_copy_by_ref:
					new_rxn.__dict__[attr] = copy(value)
			new_rxn._model = new_model
			new_model.reactions.append(new_rxn)
			# Copy the custom rates
			if rxn in self.custom_rates.keys():
				custom_rate_expr = self.custom_rates[rxn]
				new_model.custom_rates[new_rxn] = custom_rate_expr
			# Update awareness
			for metab, stoic in iteritems(rxn._metabolites):
				new_metab = new_model.metabolites.get_by_id(metab.id)
				new_rxn._metabolites[new_metab] = stoic
				new_metab._reaction.add(new_rxn)
			for gene in rxn._genes:
				new_gene = new_model.genes.get_by_id(gene.id)
				new_rxn._genes.add(new_gene)
				new_gene._reaction.add(new_rxn)


		# Create the new stoichiometric matrix for the model
		massmodel._S = create_stoichiometric_matrix(massmodel,
						matrix_type=self._matrix_type,
						dtype=self._dtype, update_model=True)


		# Refresh contexts for the new model copy
		new_model._contexts = []
		return new_model

	def merge_models(self, second_model, prefix_existing=None, inplace=False,
					new_model_id=None):
		"""Merge two massmodels to create one MassModel object with
		the reactions and metabolites from both massmodels.

		Initial conditions and custom rate laws will also be added
		from the second model into the first model. However, initial conditions
		and custom rate laws are assumed to be the same if they have the same
		identifier and therefore will not be added.

		Parameters
		----------
		second_model : MassModel
			The other MassModel to add reactions and metabolites from
		prefix_existing : string_types
			Use the string to prefix the reaction identifier of a reaction
			in the second_model if that reaction already exists within
			the first model.
		inplace : bool
			If True, add the contents directly into the first model.
			If False, a new MassModel object is created and the first model is
			left untouched. When done within the model as context, changes to
			the models are reverted upon exit.
		new_id : String or None
			Will create a new model ID for the merged model if a string
			is given. Otherwise will just use the model IDs of the first model
			if inplace is True or create a combined ID if inplace is false.
		"""
		# Check inputs to ensure they are correct types
		if not isinstance(second_model, MassModel):
			raise TypeError("The second model to merge must be a MassModel")
		if not isinstance(prefix_existing, string_types) and \
			prefix_existing is not None:
			raise TypeError("prefix_existing must be a string or none")
		if not isinstance(inplace, bool):
			raise TypeError("inplace must be a bool")
		if not isinstance(new_model_id, string_types) and \
			new_model_id is not None:
			raise TypeError("new_model_id must be a string or none")

		if inplace:
			merged_model = self
		else:
			merged_model = self.copy()

		# Set the model ID
		if new_model_id is None:
			if inplace:
				new_model_id = self.id
			else:
				new_model_id = "{} & {}".format(self.id, second_model.id)
		merged_model.id = new_model_id

		# Add the reactions of the second model to the first model
		new_reactions = deepcopy(second_model.reactions)
		if prefix_existing is not None:
			existing_reactions = new_reactions.query(
				lambda rxn: rxn.id in self.reactions)
			for rxn in existing_reactions:
				rxn.id = "{}_{}".format(prefix_existing, rxn.id)
		merged_model.add_reactions(new_reactions, True)

		return merged_model

	# HTML representation
	def _repr_html_(self):
		return """
			<table>
				<tr>
					<td><strong>Name</strong></td><td>{name}</td>
				</tr><tr>
					<td><strong>Memory address</strong></td><td>{address}</td>
				</tr><tr>
					<td><strong>Stoichiometric Matrix</strong></td>
					<td>{dim_S_matrix}</td>
				</tr><tr>
					<td><strong>Matrix Type</strong></td>
					<td>{S_type}</td>
				</tr><tr>
					<td><strong>Number of Metabolites</strong></td>
					<td>{num_metabolites}</td>
				</tr><tr>
					<td><strong>Number of Reactions</strong></td>
					<td>{num_reactions}</td>
				</tr><tr>
					<td><strong>Number of Genes</strong></td>
					<td>{num_genes}</td>
				</tr><tr>
					<td><strong>Number of Parameters</strong></td>
					<td>{num_param}</td>
				</tr><tr>
					<td><strong>Number of Initial Conditions</strong></td>
					<td>{num_ic}</td>
				</tr><tr>
					<td><strong>Number of Exchanges</strong></td>
					<td>{num_exchanges}</td>
				</tr><tr>
					<td><strong>Number of Irreversible Reactions</strong></td>
					<td>{num_irreversible}</td>
				</tr><tr>
					<td><strong>Matrix Rank</strong></td>
					<td>{mat_rank}</td>
				</tr><tr>
					<td><strong>Dimension of Null Space</strong></td>
					<td>{dim_null}</td>
				</tr><tr>
					<td><strong>Dimension of Left Null Space</strong></td>
					<td>{dim_left_null}</td>
				</tr><tr>
					<td><strong>Number of Custom Rates</strong></td>
					<td>{num_custom_rates}</td>
				</tr><tr>
					<td><strong>Compartments</strong></td>
					<td>{compartments}</td>
				</tr>
			</table>
		""".format(name=self.id, address='0x0%x' % id(self),
					dim_S_matrix="{}x{}".format(self.S.shape[0],
												self.S.shape[1]),
					S_type="{}, {}".format(self._matrix_type,
									 self._dtype.__name__),
					num_metabolites=len(self.metabolites),
					num_reactions=len(self.reactions),
					num_genes=len(self.genes),
					num_param='FIXME: Get from a parameter summary',
					num_ic= len(self.initial_conditions),
					num_exchanges=len(self.exchanges),
					num_irreversible=len(self.get_irreversible_reactions),
					mat_rank=matrix_rank(self.S),
					dim_null=nullspace(self.S).shape[1],
					dim_left_null=left_nullspace(self.S).shape[1],
					num_custom_rates=len(self.custom_rates),
					compartments=", ".join(v if v else k for \
										k,v in iteritems(self.compartments))
					)
	# Module Dunders
	def __enter__(self):
		"""Record all future changes to the MassModel, undoing them when a
		call to __exit__ is received

		Identical to the method in cobra.core.model
		"""
		# Create a new context and add it to the stack
		try:
			self._contexts.append(HistoryManager())
		except AttributeError:
			self._contexts = [HistoryManager()]

		return self


	def __exit__(self, type, value, traceback):
		"""Pop the top context manager and trigger the undo functions

		Identical to the method in cobra.core.model
		"""
		context = self._contexts.pop()
		context.reset()

	def __setstate__(self, state):
		"""Make sure all Objects in the MassModel point to the model

		Similar to the method in cobra.core.model
		"""
		self.__dict__.update(state)
		for attr in ['reactions', 'metabolites', 'genes']:
			for x in getattr(self, attr):
				x._model = self
		if not hasattr(self, "name"):
			self.name = ""

	def __getstate__(self):
		"""Get the state for serialization.

		Ensures that the context stack is cleared prior to serialization,
		since partial functions cannot be pickled reliably

		Identical to the method in cobra.core.model
		"""
		odict = self.__dict__.copy()
		odict['_contexts'] = []
		return odict
