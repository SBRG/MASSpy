# -*- coding: utf-8 -*-

# Import necessary packages

### Compatibility with Python 2.7
from __future__ import absolute_import

### Other necesary packages
import re

### cobra packages
from cobra.core.species import Species

# Class begins

### Regular expression for element parsing
element_re = re.compile("([A-Z][a-z]?)([0-9.]+[0-9.]?|(?=[A-Z])?)")

### Class definition
class Metabolite(Species):
	"""Metabolite is a class for holding information regarding a metabolite
	in mass.KineticReaction Object

	Parameters
	----------
	id: string
		The identifier associated with the metabolite
	name: string
		A human readable name
	formula: string or None
		Chemical formula associated with the metabolite
	charge: float or None
		The charge number associated with the metabolite
	compartment: string or None
		The compartment where the metabolite is located"""

	def __init__(self, id=None, name = "", formula= None, charge=None, compartment=None):
		'''Initialize the Metabolite Object'''
		Species.__init__(self, id, name)

		self._formula = formula
		self._charge = charge
		self._compartment = compartment
		self._initial_condition = None
		self._gibbs_formation= None ### Gibbs energy of formation associated with this metabolite



	### Properties
	@property
	def formula(self):
		'''Get the metabolite formula'''
		return self._formula

	@property
	def charge(self):
		'''Get the metabolite charge'''
		return self._charge

	@property
	def compartment(self):
		'''Get the metabolite compartment'''
		return self._compartment

	@property
	def elements(self):
		""" Dictionary of elements as keys and their count in the metabolite
		as integers. Identical with cobra.core.metabolite"""
		tmp_formula = self.formula
		if tmp_formula is None:
			return {}
		# necessary for some old pickles which use the deprecated
		# Formula class
		tmp_formula = str(self.formula)
		# commonly occuring characters in incorrectly constructed formulas
		if "*" in tmp_formula:
			warn("invalid character '*' found in formula '%s'" % self.formula)
			tmp_formula = tmp_formula.replace("*", "")
		if "(" in tmp_formula or ")" in tmp_formula:
			warn("invalid formula (has parenthesis) in '%s'" % self.formula)
			return None
		composition = {}
		parsed = element_re.findall(tmp_formula)
		for (element, count) in parsed:
			if count == '':
				count = 1
			else:
				try:
					count = float(count)
					int_count = int(count)
					if count == int_count:
						count = int_count
					else:
						warn("%s is not an integer (in formula %s)" %
							 (count, self.formula))
				except ValueError:
					warn("failed to parse %s (in formula %s)" %
						 (count, self.formula))
					return None
			if element in composition:
				composition[element] += count
			else:
				composition[element] = count
		return composition

	@property
	def initial_condition(self):
		'''Get initial condition of the metabolite

		Warnings
		--------
		In most cases, accessing the initial condition through the massmodel
		is more accurate than accessing it through the metabolite.
		'''
		return self._initial_condition

	@initial_condition.setter
	def initial_condition(self, value):
		'''Set the initial condition of the metabolite'''
		self._initial_condition = value

	@property
	def gibbs_formation(self):
		return self._gibbs_formation

	@gibbs_formation.setter
	def gibbs_formation(self, value):
		self._gibbs_formation = value

	@property
	def formula_weight(self):
		"""Calculate the formula weight. Identical with cobra.core.metabolite"""
		try:
			return sum([count * elements_and_molecular_weights[element]
						for element, count in self.elements.items()])
		except KeyError as e:
			warn("The element %s does not appear in the peridic table" % e)



	### Methods
	def _set_id_with_model(self, value):
		'''Set the internal id of the Metabolite object to
		the associated model. Mirrors the method in cobra.core.metabolite'''
		if value in self.massmodel.metabolites:
			raise ValueError("You stupid")
		self._id = value
		self.massmodel.metabolites._generate_index()


	def remove_from_model(self, destructive=False):
		"""Removes the metabolite association from self.massmodel

		The change is reverted upon exit when using the model as a context.
		Mirrors the method in cobra.core.metabolite

		Parameters
		----------
		destructive : bool
			If False then the metabolite is removed from all
			associated reactions.  If True then all associated
			reactions are removed from the Model."""
		return self._model.remove_metabolites(self, destructive)

	### Shorthands
	@property
	def ic(self):
		'''Shorthand getter for initial condition'''
		return initial_conditions(self)

	@ic.setter
	def ic(self, value):
		'''Shorthand setter for initial condition'''
		initial_conidition(self, value)

	@property
	def gf(self):
		'''Shorthand getter for Gibb's energy of formation'''
		return gibbs_formation(self)

	@gf.setter
	def gf(self, value):
		'''Shorthand setter for Gibb's energy of formation'''
		gibbs_formation(self, value)

	### Dictionary of elements and their corresponding atomic mass
	elements_and_molecular_weights = {
		'H':   1.007940,
		'He':  4.002602,
		'Li':  6.941000,
		'Be':  9.012182,
		'B':   10.811000,
		'C':   12.010700,
		'N':   14.006700,
		'O':   15.999400,
		'F':   18.998403,
		'Ne':  20.179700,
		'Na':  22.989770,
		'Mg':  24.305000,
		'Al':  26.981538,
		'Si':  28.085500,
		'P':   30.973761,
		'S':   32.065000,
		'Cl':  35.453000,
		'Ar':  39.948000,
		'K':   39.098300,
		'Ca':  40.078000,
		'Sc':  44.955910,
		'Ti':  47.867000,
		'V':   50.941500,
		'Cr':  51.996100,
		'Mn':  54.938049,
		'Fe':  55.845000,
		'Co':  58.933200,
		'Ni':  58.693400,
		'Cu':  63.546000,
		'Zn':  65.409000,
		'Ga':  69.723000,
		'Ge':  72.640000,
		'As':  74.921600,
		'Se':  78.960000,
		'Br':  79.904000,
		'Kr':  83.798000,
		'Rb':  85.467800,
		'Sr':  87.620000,
		'Y':   88.905850,
		'Zr':  91.224000,
		'Nb':  92.906380,
		'Mo':  95.940000,
		'Tc':  98.000000,
		'Ru':  101.070000,
		'Rh':  102.905500,
		'Pd':  106.420000,
		'Ag':  107.868200,
		'Cd':  112.411000,
		'In':  114.818000,
		'Sn':  118.710000,
		'Sb':  121.760000,
		'Te':  127.600000,
		'I':   126.904470,
		'Xe':  131.293000,
		'Cs':  132.905450,
		'Ba':  137.327000,
		'La':  138.905500,
		'Ce':  140.116000,
		'Pr':  140.907650,
		'Nd':  144.240000,
		'Pm':  145.000000,
		'Sm':  150.360000,
		'Eu':  151.964000,
		'Gd':  157.250000,
		'Tb':  158.925340,
		'Dy':  162.500000,
		'Ho':  164.930320,
		'Er':  167.259000,
		'Tm':  168.934210,
		'Yb':  173.040000,
		'Lu':  174.967000,
		'Hf':  178.490000,
		'Ta':  180.947900,
		'W':   183.840000,
		'Re':  186.207000,
		'Os':  190.230000,
		'Ir':  192.217000,
		'Pt':  195.078000,
		'Au':  196.966550,
		'Hg':  200.590000,
		'Tl':  204.383300,
		'Pb':  207.200000,
		'Bi':  208.980380,
		'Po':  209.000000,
		'At':  210.000000,
		'Rn':  222.000000,
		'Fr':  223.000000,
		'Ra':  226.000000,
		'Ac':  227.000000,
		'Th':  232.038100,
		'Pa':  231.035880,
		'U':   238.028910,
		'Np':  237.000000,
		'Pu':  244.000000,
		'Am':  243.000000,
		'Cm':  247.000000,
		'Bk':  247.000000,
		'Cf':  251.000000,
		'Es':  252.000000,
		'Fm':  257.000000,
		'Md':  258.000000,
		'No':  259.000000,
		'Lr':  262.000000,
		'Rf':  261.000000,
		'Db':  262.000000,
		'Sg':  266.000000,
		'Bh':  264.000000,
		'Hs':  277.000000,
		'Mt':  268.000000,
		'Ds':  281.000000,
		'Rg':  272.000000,
		'Cn':  285.000000,
		'Nh':  286.000000,
		'Fl':  289.000000,
		'Mc':  290.000000,
		'Lv':  293.000000,
		'Ts':  294.000000,
		'Og':  294.000000,
	}
