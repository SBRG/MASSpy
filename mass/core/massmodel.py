# -*- coding: utf-8 -*-

# Compatibility with Python 2.7
from __future__ import absolute_import

# Import necesary packages
import logging
import re
import pandas as pd
import numpy as np
import sympy as sp
from six import string_types, integer_types, iteritems, iterkeys, itervalues
from copy import copy, deepcopy
from functools import partial
from warnings import warn

# from cobra
from cobra.core.object import Object
from cobra.core.dictlist import DictList
from cobra.util.context import HistoryManager, resettable, get_context

# from mass
from mass.util import array
from mass.util import qcqa
from mass.core import expressions
from mass.core.massmetabolite import MassMetabolite
from mass.core.massreaction import MassReaction

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
	custom_parameters : dict
		A dictionary to store custom parameters for custom rates,
		keys are parameters and values are the parameter value.
		Custom rates will always have preference over rate laws in reactions.
	fixed_concentrations: dict
		A dictionary to store fixed concentrations for metabolites, where
		keys are metabolite objects or "external metabolites" of exchange
		reactions as strings, and values are the fixed concentrations.
		Fixed concentrations will always have preference over the metabolite
		ODEs representing concentrations.
	compartments : dict
		A dictionary to store the compartment shorthands and their full names.
		Keys are the shorthands and values are the full names.
		Example: {'c': 'cytosol'}
	units : dict
		A dictionary to store the units used in the model for referencing.

		WARNING: Note that the model will not track the units,
		Therefore all unit conversions must be manually in order to ensure
		numerical consistency in the model. It is highly recommended to stick
		with the following units:

		{'N': 'Millimoles', 'Vol': 'Liters', 'Time': 'Hours'}
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
			#For storing of custom rate laws and fixed concentrations
			self._rtype = 1
			self._custom_rates= dict()
			self._custom_parameters = dict()
			self.fixed_concentrations = dict()
			# A dictionary of the compartments in the model
			self.compartments = dict()
			# A dictionary to store the units utilized in the model.
			self.units = dict()
			# Internal storage of S matrix and data types for updating S
			self._matrix_type = matrix_type
			self._dtype = dtype
			# For storing the stoichiometric matrix.
			self._S = array.create_stoichiometric_matrix(self,
                                    matrix_type=self._matrix_type,
                                    dtype=self._dtype,
            						update_model=True)
			# For storing the HistoryManager contexts
			self._contexts = []

	# Properties
	@property
	def S(self):
		"""Get the Stoichiometric Matrix of the MassModel"""
		return self.update_S(matrix_type=self._matrix_type, dtype=self._dtype,
							update_model=False)

	@property
	def rates(self):
		""""Return the rate laws for the reactions as human readable strings
		in a dictionary where keys are the reaction objects and values are the
		rate laws
		"""
		rate_dict =  {rxn: rxn.generate_rate_law(rate_type=self._rtype,
								sympy_expr=False, update_reaction=True)
								for rxn in self.reactions}
		if self.custom_rates != {}:
			for rxn, custom_expression in iteritems(self.custom_rates):
				rate_dict.update({rxn : str(custom_expression)})
		return rate_dict

	@property
	def rate_expressions(self):
		"""Get the rate laws for the reactions as sympy expressions in a
		dictionary where keys are the reaction objects and values are the
		sympy rate law expressions
		"""
		rate_dict =  {rxn: rxn.generate_rate_law(rate_type=self._rtype,
								sympy_expr=True, update_reaction=True)
								for rxn in self.reactions}
		if self.custom_rates != {}:
			rate_dict.update(self.custom_rates)
		return rate_dict

	@property
	def odes(self):
		"""Get the ODEs for the metabolites as sympy expressions where
		keys are the metabolite objects and values are the ODE expressions
		"""
		return {metab: metab.ode for metab in self.metabolites
				if metab not in self.fixed_concentrations}

	@property
	def exchanges(self):
		"""Get the exchange reactions in the MassModel"""
		return [rxn for rxn in self.reactions if rxn.exchange]

	@property
	def get_external_metabolites(self):
		"""Get all 'external' metabolites in the reaction. Primarily used for
		setting fixed concentrations"""
		external_set = {rxn.get_external_metabolite for rxn in self.reactions
				if rxn.exchange}
		return list(sorted(external_set))

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
		return [rxn for rxn in self.reactions if not rxn.reversible]

	@property
	def steady_state_fluxes(self):
		"""Return all steady state fluxes stored in the model reactions"""
		return {rxn: rxn.ssflux for rxn in self.reactions
				if rxn.ssflux is not None}

	@property
	def custom_rates(self):
		"""Get the sympy custom rate expressions in the MassModel"""
		return self._custom_rates

	@property
	def custom_parameters(self):
		"""Get the custom rate parameters in the MassModel"""
		return self._custom_parameters

	# Methods
	## Public
	def update_S(self, reaction_list=None, matrix_type=None, dtype=None,
				update_model=True):
		"""For internal use only. Update the S matrix of the model.

		NOTE: reaction_list is assumed to be at the end of self.reactions.

		Parameters
		----------
		model : mass.MassModel
			The MassModel object to construct the matrix for
		reaction_list : list of MassReactions or None
			The list of MassReactions to add to the current stoichiometric matrix.
			Reactions must already exist in the model in order to update.
			If None, the entire stoichiometric matrix is reconstructed
		matrix_type: string {'dense', 'dok', 'lil', 'DataFrame'}, or None
			If None, will utilize the matrix type initialized with the original
			model. Otherwise reconstruct the S matrix with the specified type.
			Types can include 'dense' for a standard  numpy.array, 'dok' or
			'lil' to obtain the scipy sparse matrix of the corresponding type, and
			DataFrame for a pandas 'Dataframe' where species (excluding genes)
			are row indicies and reactions are column indicices
		dtype : data-type
			The desired data-type for the array. If None, defaults to float

		Returns
		-------
		matrix of class 'dtype'
			The stoichiometric matrix for the given MassModel
		"""
		return array._update_S(self, reaction_list=reaction_list,
						matrix_type=matrix_type, dtype=dtype,
						update_model=update_model)

	def add_metabolites(self, metabolite_list, add_initial_conditons=False):
		"""Will add a list of metabolites to the MassModel object and add
		the MassMetabolite initial conditions accordingly.

		The change is reverted upon exit when using the MassModel as a context.

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
			self.set_initial_conditions(metabolite_list)

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

	def set_initial_conditions(self, metabolite_list=None):
		"""Set the initial conditions for a list of metabolites in the model.

		The metabolite must already exist in the model in order to set its
		initial condition. Initial conditions stored in the metabolite are
		accessed and added to the model. Any existing initial condition for
		a metabolite in the model is replaced.

		The change is reverted upon exit when using the MassModel as a context.

		Parameters
		----------
		metabolite_list : list or None
			A list of MassMetabolite objects. If None, will use all metabolites
			in the model
		"""
		if metabolite_list is None:
			metabolite_list = self.metabolites
		# If the metabolite list is not a list
		elif not hasattr(metabolite_list, '__iter__'):
			metabolite_list = DictList([metabolite_list])
		else:
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
			for key in iterkeys(ic_dict):
				if key not in iterkeys(existing_ics):
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
			for metab in list(iterkeys(rxn._metabolites)):
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
			self.update_S(reaction_list=reaction_list, update_model=True)

		if context:
			context(partial(self.reactions.__isub__, reaction_list))
			for rxn in reaction_list:
				context(partial(setattr, rxn, '_model', None))
			if update_stoichiometry:
				context(partial(self.update_S, None, None, None, True))


	def remove_reactions(self, reaction_list, remove_orphans=False,
						update_stoichiometry=False):
		"""Remove MassReactions from the MassModel

		The change is reverted upon exit when using the MassModel as a context.

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
			self.update_S(reaction_list=None, update_model=True)

		if context:
			context(partial(self.reactions.__iadd__, reaction_list))
			for rxn in reaction_list:
				context(partial(setattr, rxn, '_model', self))
			if update_stoichiometry:
				context(partial(self.update_S, reaction_list,
							None, None, True))

	def add_exchange(self, metabolite, exchange_type="exchange",
						external_concentration=0.,
						update_stoichiometry=False):
		"""Add an exchange reaction for a given metabolite using the
		pre-defined exchange types "exchange" for reversibly into or exiting
		the compartment, "source" for irreversibly into the compartment,
		and "demand" for irreversibly exiting the compartment.

		The change is reverted upon exit when using the MassModel as a context.

		Parameters
		----------
		metabolite : MassMetabolite
			Any given metabolite to create an exchange for.
		exchange_type : string, {"demand", "source", "exchange"}
			The type of exchange reaction to create are not case sensitive.
		reversible : bool
			If True, exchange is reversible. When using a user-defined type,
			must specify the reversiblity
		"""
		# Check whether metabolite is a MassMetabolite object
		if not isinstance(metabolite, MassMetabolite):
			raise TypeError("metabolite must be a MassMetabolite object")

		type_dict = {"source":["S", 1, False],
					"demand":["DM", -1, False],
					"exchange":["EX", -1, True]}

		# Set the type of exchange
		if not isinstance(exchange_type, string_types):
			raise TypeError("exchange_type must be a string")
		else:
			exchange_type = exchange_type.lower()

		if exchange_type in type_dict:
			values = type_dict[exchange_type]
			rxn_id = "{}_{}".format(values[0], metabolite.id)
			if rxn_id in self.reactions:
				warn("Reaction %s already exists in model" % rxn_id)
				return None
			c = values[1]
			reversible = values[2]
		else:
			raise TypeError("Exchange type must be either "
							"'exchange''source' or 'sink'")
		rxn_name = "{} {}".format(metabolite.name, exchange_type)

		rxn = MassReaction(id=rxn_id, name=rxn_name,
					subsystem="Transport/Exchange",reversible=reversible)
		rxn.add_metabolites({metabolite: c})
		self.add_reactions([rxn], update_stoichiometry)
		self.add_fixed_concentrations(fixed_conc_dict={
			rxn.get_external_metabolite : external_concentration})

	def generate_rate_laws(self, reaction_list=None, rate_type=1,
							sympy_expr=False, update_reactions=False):
		"""Get the rate laws for a list of reactions in a MassModel and return
		them as human readable strings or as sympy expressions for simulations.

		The type determines which rate law format to return.
		For example: A <=> B

		type=1: kf*(A - B/Keq)
		type=2: kf*A - kr*B
		type=3: kr*(Keq*A - B)

		Parameters
		----------
		reaction_list = list or None
			The list of reactions to obtain the rates for. If none specified,
			will return the rates for all reactions in the MassModel
		rate_type : int {1, 2, 3}
			The type of rate law to display. Must be 1, 2, of 3.
			type 1 will utilize kf and Keq,
			type 2 will utilize kf and kr,
			type 3 will utilize kr and Keq.
		sympy_expr : bool
			If True, will output a sympy expression, otherwise
			will output a human readable string.

		Returns
		-------
		dict of reaction rates where keys are reaction identifiers and
			values are the strings or sympy expressions
		"""
		# Check inputs
		if not isinstance(rate_type, (integer_types, float)):
			raise TypeError("rate_type must be an int or float")
		elif not isinstance(sympy_expr, bool):
			raise TypeError("sympy_expr must be a bool")
		else:
			rate_type = int(rate_type)
			if rate_type not in {1, 2, 3}:
				raise ValueError("rate_type must be 1, 2, or 3")

		# Use massmodel reactions if no reaction list is given
		if reaction_list is None:
			reaction_list = self.reactions
		# If the reaction list is not a list
		elif not hasattr(reaction_list, '__iter__'):
			reaction_list = DictList([reaction_list])
		else:
			reaction_list = DictList(reaction_list)

		if len(reaction_list) == 0:
			return None
		if update_reactions:
			self._rtype = rate_type
		# Get the rates
		rates = {rxn :
				rxn.generate_rate_law(rate_type, sympy_expr, update_reactions)
				for rxn in reaction_list}
		if self.custom_rates != {}:
			rates.update(self.custom_rates)
		return rates

	def get_mass_action_ratios(self, reaction_list=None,sympy_expr=False):
		"""Get the mass action ratios for a list of reactions in a MassModel
		and return them as human readable strings or as sympy expressions
		for simulations

		Parameters
		----------
		reaction_list = list of MassReactions or None
			The list of MassReactions to obtain the disequilibrium ratios for.
			If None, will return the rates for all reactions in the MassModel
		sympy_expr : bool
			If True, will output sympy expressions, otherwise
			will output a human readable strings.
		Returns
		-------
		dict of disequilibrium ratios where keys are reaction identifiers and
			values are mass action ratios as strings or sympy expressions
		"""
		# Use massmodel reactions if no reaction list is given
		if reaction_list is None:
			reaction_list = self.reactions
		# If the reaction list is not a list
		elif not hasattr(reaction_list, '__iter__'):
			reaction_list = DictList([reaction_list])
		else:
			reaction_list = DictList(reaction_list)

		if len(reaction_list) == 0:
			return None
		# Get the mass action ratios
		return {rxn : rxn.get_mass_action_ratio(sympy_expr)
		 		for rxn in reaction_list}

	def get_disequilibrium_ratios(self, reaction_list=None, sympy_expr=False):
		"""Get the disequilibrium ratios for a list of reactions in a MassModel
		and return them as human readable strings or as sympy expressions
		for simulations

		Parameters
		----------
		reaction_list = list of MassReactions or None
			The list of MassReactions to obtain the disequilibrium ratios for.
			If None, will return the rates for all reactions in the MassModel
		sympy_expr : bool
			If True, will output sympy expressions, otherwise
			will output a human readable strings.
		Returns
		-------
		dict of disequilibrium ratios where keys are reaction identifiers and
			values are disequilibrium ratios as strings or sympy expressions
		"""
		# Use massmodel reactions if no reaction list is given
		if reaction_list is None:
			reaction_list = self.reactions
		# If the reaction list is not a list
		elif not hasattr(reaction_list, '__iter__'):
			reaction_list = DictList([reaction_list])
		else:
			reaction_list = DictList(reaction_list)

		if len(reaction_list) == 0:
			return None
		# Get the disequilibrium ratios
		return {rxn : rxn.get_disequilibrium_ratio(sympy_expr)
				for rxn in reaction_list}

	def add_custom_rate(self, reaction, custom_rate,
						custom_parameters=None):
		"""Add a custom rate to the MassModel for a reaction.

		Note: Metabolites must already be in the MassReaction

		The change is reverted upon exit when using the MassModel as a context.

		Parameters
		----------
		reaction : mass.MassReaction
			The MassReaction which the custom rate associated with
		custom_rate_law :  string
			The custom rate law as a string. The string representation of the
			custom rate lawwill be used to create a sympy expression that
			represents the custom rate law
		custom_parameters :  dictionary of strings
			A dictionary of where keys are custom parameters in the custom rate
			as strings, and values are their numerical value. The string
			representation of the custom parameters will be used to create the
			symbols in the sympy expressions of the custom rate law.
			If None, parameters are assumed to be already in the MassModel,
			or one of the MassReaction rate or equilibrium constants.
		"""
		# Get the custom parameters
		if custom_parameters is not None:
			custom_param_list = list(iterkeys(custom_parameters))
		else:
			custom_parameters = {}
			custom_param_list = []
		# Use existing ones if they are in the rate law
		existing_customs = self.custom_parameters
		if len(existing_customs) != 0:
			for custom_param in iterkeys(existing_customs):
				if re.search(custom_param, custom_rate) and \
					custom_param not in custom_param_list:
					custom_param_list.append(custom_param)

		custom_rate_expression = expressions.create_custom_rate(reaction,
												custom_rate, custom_param_list)

		self._custom_rates.update({reaction : custom_rate_expression})
		self._custom_parameters.update(custom_parameters)

		context = get_context(self)
		if context:
			context(partial(self._custom_rates.pop, reaction))
			for key in custom_param_list:
				if key in iterkeys(self._custom_parameters):
					context(partial((self._custom_parameters.pop, key)))
			context(partial(self._custom_parameters.update, existing_customs))


	def remove_custom_rate(self, reaction):
		"""Remove a custom rate to the MassModel for a reaction If no other custom
		rates rely on those custom parameters, remove those custom parameters from
		the model as well.

		The change is reverted upon exit when using the MassModel as a context.

		Parameters
		----------
		reaction : mass.MassReaction
			The MassReaction which the custom rate associated with
		"""
		# Remove the custom rate law
		custom_rate_to_remove = self.custom_rates[reaction]
		del self.custom_rates[reaction]

		# Remove custom parameters if they are not associated with any other
		# custom rate expression
		symbols = custom_rate_to_remove.atoms(sp.Symbol)
		if len(self.custom_rates) != 0:
			other_syms = set()
			for custom_rate in itervalues(self.custom_rates):
				for sym in list(custom_rate.atoms(sp.Symbol)):
					other_syms.add(sym)
			for sym in other_syms:
				if sym in symbols.copy():
					symbols.remove(sym)

		context = get_context(self)
		if context:
			existing = dict((str(sym), self._custom_parameters[str(sym)])
							for sym in symbols)

		for sym in symbols:
			del self._custom_parameters[str(sym)]

		if context:
			context(partial(self._custom_rates.update, {reaction:
											custom_rate_to_remove}))
			context(partial(self._custom_parameters.update, existing))

	def reset_custom_rates(self):
		"""Reset all custom rate laws and parameters in a model.

		Warnings
		--------
		Will remove all custom rates and custom rate parameters in the
		MassModel. To remove a specific rate(s) without affecting the others,
		use the remove_custom_rate method instead.
		"""
		self._custom_rates = {}
		self._custom_parameters = {}
		print("Reset all custom rate laws")

	def get_elemental_matrix(self, matrix_type=None, dtype=None):
		"""Get the elemental matrix of a model

		Parameters
		----------
		matrix_type: string {'dense', 'dok', 'lil', 'DataFrame'}, or None
			If None, will utilize the matrix type initialized with the original
			model. Otherwise reconstruct the S matrix with the specified type.
			Types can include 'dense' for a standard  numpy.array, 'dok' or
			'lil' to obtain the scipy sparse matrix of the corresponding type, and
			DataFrame for a pandas 'Dataframe' where species (excluding genes)
			are row indicies and reactions are column indicices
		dtype : data-type
			The desired data-type for the array. If None, defaults to float

		Returns
		-------
		matrix of class 'dtype'
			The elemental matrix for the given MassModel
		"""
		# Set defaults for the elemental matrix
		if matrix_type is None:
			matrix_type = 'dataframe'
		if dtype is None:
			dtype = np.int64

		# No need to construct a matrix if there are no metabolites
		if len(self.metabolites) == 0:
			return None

		CHOPNSq = ['C', 'H', 'O', 'P', 'N', 'S', 'q' ]

		(matrix_constructor, dtype) = array._setup_matrix_constructor(
												self, matrix_type, dtype)

		e_matrix = matrix_constructor((len(CHOPNSq),len(self.metabolites)),
									dtype=dtype)
		# Get index for elements and metabolites
		e_ind = CHOPNSq.index
		m_ind = self.metabolites.index

		# Build the matrix
		for metab in self.metabolites:
			for element in CHOPNSq:
				if element in iterkeys(metab.elements):
					amount = metab.elements[element]
				elif element == 'q' and metab.charge is not None:
					amount = metab.charge
				else:
					amount = 0
				e_matrix[e_ind(element), m_ind(metab)] = amount
		# Convert matrix to dataframe if matrix type is a dataframe
		if matrix_type == 'dataframe':
			metab_ids = [metab.id for metab in self.metabolites]
			e_matrix = pd.DataFrame(e_matrix, index=CHOPNSq, columns=metab_ids)

		return e_matrix

	def add_fixed_concentrations(self, fixed_conc_dict=None):
		"""Add fixed concentrations for metabolites, setting their ODEs
		to a constant value during simulation of the MassModel.

		The metabolites must already exist in the model or be an
		"external" metabolite for an exchange reaction"

		Parameters
		----------
		fixed_conc_dict : dictionary
			A dictionary of fixed concentrations where
			metabolites are keys and fixed concentrations are values
		"""
		# Check inputs
		if fixed_conc_dict is None:
			return None
		if not isinstance(fixed_conc_dict, dict):
			raise TypeError("fixed_conc_dict must be a dictionary")
		for metab, fixed_conc in iteritems(fixed_conc_dict):
			if metab not in self.get_external_metabolites and \
				metab not in self.metabolites:
				raise ValueError("Did not find %s in model metabolites"
								" or exchanges" % metab)
			if not isinstance(fixed_conc, (integer_types, float)):
				raise TypeError("Fixed concentration must be an int or float")
			elif fixed_conc < 0.:
				raise ValueError("External concentration must be non-negative")
			else:
				fixed_conc = float(fixed_conc)


		# Keep track of existing initial conditions for HistoryManager
		context = get_context(self)
		if context:
			existing_ics = {metab : self.fixed_concentrations[metab]
							for metab in fixed_conc_dict
							if metab in self.fixed_concentrations}

		self.fixed_concentrations.update(fixed_conc_dict)

		if context:
			for key in iterkeys(fixed_conc_dict):
				if key not in iterkeys(existing_ics):
					context(partial(self.fixed_concentrations.pop, key))
				context(partial(self.fixed_concentrations.update,
								existing_ics))

	def remove_fixed_concentrations(self, metabolites=None):
		"""Remove a fixed concentration for a specific metabolite

		Parameters
		----------
		metabolites : list of metabolites identifiers
			A list containing MassMetabolites or their identifiers. Can also
			be an "external" metabolite for an exchange reaction
		"""
		# Check inputs
		if metabolites is None:
			return None
		if not hasattr(metabolites, '__iter__'):
			metabolites = [metabolites]

		for metab in metabolites:
			if metab not in self.get_external_metabolites and \
				metab not in self.metabolites:
				raise ValueError("Did not find %s in model metabolites"
								" or exchanges" % metab)

		# Keep track of existing initial conditions for HistoryManager
		context = get_context(self)
		if context:
			existing_ics = {metab : self.fixed_concentrations[metab]
							for metab in metabolites
							if metab in self.fixed_concentrations}

		for metab in metabolites:
			del self.fixed_concentrations[metab]

		if context:
			context(partial(self.fixed_concentrations.update, existing_ics))

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
			if metab in iterkeys(self.initial_conditions) :
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
			if rxn in iterkeys(self.custom_rates):
				new_model.custom_rates[new_rxn] = self.custom_rates[rxn]
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
		new_model._S = array.create_stoichiometric_matrix(self,
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
		if not isinstance(prefix_existing, (string_types, type(None))):
			raise TypeError("prefix_existing must be a string or none")
		if not isinstance(inplace, bool):
			raise TypeError("inplace must be a bool")
		if not isinstance(new_model_id, (string_types, type(None))):
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

	def calc_PERCS(self, steady_state_concentrations=None,
					steady_state_fluxes=None, at_equilibrium_default=100000.,
					update_reactions=False):
		"""Calculate the pseudo rate constants (rxn.forward_rate_constant)
		for reactions in the MassModel using steady state concentrations and
		steady state fluxes.

		Parameters
		----------
		steady_state_concentrations : dict or None
			A dictionary of steady state concentrations where MassMetabolites are
			keys and the concentrations are the values. If None, will utilize
			the initial conditions in the MassModel.
		steady_state_fluxes : dict or None
			A dictionary of steady state fluxes where MassReactions are the keys
			and the fluxes are the values. If None, will utilize the
			steady state reactions stored in each reaction in the model.
		at_equilibrium_default : float or None
			The value to set the pseudo order rate constant if the reaction is
			at equilibrium. If None, will default to 100,000
		update_parameters : bool
			Whether to update the forward rate constants in the MassReactions.
			If True, will update the forward rate constants inside the
			MassReactions with the calculated pseudo order rate constants
		"""
		# Check inputs
		if steady_state_concentrations is None:
			steady_state_concentrations = self.initial_conditions
		if not isinstance(steady_state_concentrations, dict):
			raise TypeError("Steady state concentrations must be a dictionary"
					"where keys are MassMetabolites and values are concentrations")

		if not isinstance(steady_state_fluxes, (dict, type(None))):
			raise TypeError("Steady state fluxes must be a dictionary where"
							" keys are MassReactions and values are fluxes")

		if not isinstance(at_equilibrium_default, (integer_types, float)):
			raise TypeError("at_equilibrium_default must be an int or float")

		if not isinstance(update_reactions, bool):
			raise TypeError("update_reactions must be a bool")

		# Use model steady state fluxes and check for missing parameters
		if steady_state_fluxes is None:
			steady_state_fluxes = self.steady_state_fluxes
			missing_params = qcqa.get_missing_parameters(self, Keq=True,
										ssflux=True, custom_parameters=True)
		# Use the given steady state fluxes and check for missing parameters
		else:
			missing_params = qcqa.get_missing_parameters(self, Keq=True,
										ssflux=False, custom_parameters=True)
			for rxn in self.reactions:
				if rxn not in iterkeys(steady_state_fluxes):
					if rxn in iterkeys(missing_params):
						missing = missing_params[rxn]
						missing.append( "ssflux")
						missing_params[rxn] = missing
					else:
						missing_params[rxn] = ["ssflux"]

		# Use model initial conditions for the steady state concentratiosn
		# and check for missing initial conditions
		if steady_state_concentrations is None:
			missing_concs = qcqa.get_missing_initial_conditions(self)
		# Use the given steady state concentrations and
		# check for missing concentrations
		else:
			missing_concs = [m for m in self.metabolites
							if m not in iterkeys(steady_state_concentrations)]

		# If parameters or concentrations are missing, print a warning,
		# inform what values are missing, and return none
		if len(missing_params) != 0 or len(missing_concs) != 0:
			warn("\nCannot calculate PERCs due to missing values")
			reports = qcqa._qcqa_summary([missing_concs, missing_params])
			for report in reports:
				print("%s\n" % report)
			return None

		#  Group symbols in rate_laws
		if self._rtype != 1:
			self._rtype = 1
		odes, rates, symbols = expressions._sort_symbols(self)
		metabolites = symbols[0]
		rate_params = symbols[1]
		fixed_concs = symbols[2]
		custom_params = symbols[3]

		# Strip the time dependency
		rates = expressions.strip_time(rates)
		percs_dict = {}
		# Get values to subsitute into equation.
		for rxn, rate in iteritems(rates):
			values = {}
			symbols = rate.atoms(sp.Symbol)
			for sym in symbols:
				if sym in rate_params:

					if re.search("Keq", str(sym)):
						values.update({sym : rxn.Keq})
					else:
						perc = sym
				elif sym in custom_params:
					values.update({sym : self.custom_parameters[str(sym)]})
				elif sym in fixed_concs:
					values.update({sym : self.fixed_concentrations[str(sym)]})
				else:
					metab = self.metabolites.get_by_id(str(sym))
					values.update({sym : self.initial_conditions[metab]})

			flux = steady_state_fluxes[rxn]

			# Set equilbrium default if no flux
			if flux == 0:
				sol = {float(at_equilibrium_default)}
			# Otherwise calculate the PERC
			else:
				equation = sp.Eq(steady_state_fluxes[rxn], rate.subs(values))
				sol = set(sp.solveset(equation, perc, domain=sp.S.Reals))
			percs_dict.update({str(perc): float(sol.pop())})

			if update_reactions:
				rxn.kf = percs_dict[str(perc)]

		return percs_dict

	## Internal
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
				</tr><tr>
					<td><strong>Units</strong></td>
					<td>{units}</td>
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
					num_param=sum([len(rxn.parameters)
									for rxn in self.reactions] + \
									[len(self.fixed_concentrations)]),
					num_ic= len(self.initial_conditions),
					num_exchanges=len(self.exchanges),
					num_irreversible=len(self.get_irreversible_reactions),
					mat_rank=array.matrix_rank(self.S),
					dim_null=array.nullspace(self.S,'row').shape[0],
					dim_left_null=array.left_nullspace(self.S, 'row').shape[0],
					num_custom_rates=len(self.custom_rates),
					compartments=", ".join(v if v else k for \
										k,v in iteritems(self.compartments)),
					units=", ".join(v if v else k for \
										k,v in iteritems(self.units))
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
