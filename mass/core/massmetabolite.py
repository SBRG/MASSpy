# -*- coding: utf-8 -*-

# Compatibility with Python 2.7
from __future__ import absolute_import

# Import necesary packages
import re
from six import string_types, integer_types
from warnings import warn

# cobra packages
from cobra.core.species import Species

# Class begins

# Regular expression for element parsing
element_re = re.compile("([A-Z][a-z]?)([0-9.]+[0-9.]?|(?=[A-Z])?)")

# Class definition
class MassMetabolite(Species):
	"""MassMetabolite is a class for holding information regarding a metabolite
	in mass.MassReaction Object

	Parameters
	----------
	id : string
		The identifier associated with the metabolite
	name : string
		A human readable name
	formula : string or None
		Chemical formula associated with the metabolite
	charge : float or None
		The charge number associated with the metabolite
	compartment : string or None
		The compartment where the metabolite is located
	"""

	def __init__(self, id=None, name ="", formula=None, charge=None,
					compartment=None):
		# Check inputs to ensure they are the correct types
		"""Initialize the MassMetabolite Object"""
		if not isinstance(name, string_types):
			raise TypeError("name must be a string type")
		if not isinstance(formula, string_types) and formula != None:
			raise TypeError("formula must be a string type")
		if not isinstance(charge, integer_types)  \
			and not isinstance(charge, float) and charge != None:
				raise TypeError("charge must be a number")
		if not isinstance(compartment, string_types) and compartment != None:
			raise TypeError("compartment must be a string or none")

		Species.__init__(self, id, name)
		# Chemical formula of the metabolite
		self._formula = formula
		# Charge associated with the metabolite
		self._charge = charge
		# Compartment in which metabolite is located
		self._compartment = compartment
		# Initial condition associated with this metabolite=
		self._initial_condition = None
		# Gibbs energy of formation associated with this metabolite
		self._gibbs_formation = None
		# Ordinary differential equation for the metabolite concentration
		self._ode = None

	# Properties
	@property
	def formula(self):
		"""Returns the metabolite formula"""
		return self._formula

	@property
	def charge(self):
		"""Returns the metabolite charge"""
		return self._charge

	@property
	def compartment(self):
		"""Returns the metabolite compartment"""
		return self._compartment

	@property
	def elements(self):
		"""Dictionary of elements as keys and their count in the
		metabolite as integers.

		Identical to the method in cobra.core.metabolite
		"""
		tmp_formula = self.formula
		if tmp_formula is None:
			return {}
		if not isinstance(tmp_formula, string_types):
			raise TypeError("Formula must be a string")
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
		"""Returns the initial condition of the metabolite

		Warnings
		--------
		In most cases, accessing initial conditions through the massmodel
		is more accurate than accessing it through the metabolite.
		"""
		return self._initial_condition

	@initial_condition.setter
	def initial_condition(self, value):
		"""Set the initial condition of the metabolite.

		Warnings
		--------
		Initial concentrations of metabolites cannot be negative.
		"""
		if not isinstance(value, integer_types) and \
			not isinstance(value, float):
			raise TypeError("Initial condition must be an integer or float")
		if value < 0.:
			raise ValueError("Initial condition must be a non-negative number")
		self._initial_condition = value

	@property
	def gibbs_formation(self):
		"""Returns the Gibbs formation energy of the metabolite"""
		return self._gibbs_formation

	@gibbs_formation.setter
	def gibbs_formation(self, value):
		"""Set the Gibbs formation energy of the metabolite"""
		if not isinstance(value, integer_types) and \
			not isinstance(value, float):
			raise TypeError("Initial condition must be an integer or float")

		self._gibbs_formation = value

	@property
	def formula_weight(self):
		"""Calculate the formula weight.

		Identical with cobra.core.metabolite
		"""
		try:
			return sum([count * elements_and_molecular_weights[element]
						for element, count in self.elements.items()])
		except KeyError as e:
			warn("The element %s does not appear in the peridic table" % e)

	@property
	def ode(self):
		"""Return the ordinary differential equation associated for the
		concentration of the metabolite

		Will return None if metabolite is not associated with a MassReaction
		"""
		return self._ode

	# Methods
	def _set_id_with_model(self, value):
		"""Set the id of the MassMetabolite object to the associated massmodel.

		Similar to the method in cobra.core.metabolite.
		"""
		if value in self.massmodel.metabolites:
			raise ValueError("The massmodel already contains a metabolite with \
								the id:")
		self._id = value
		self.massmodel.metabolites._generate_index()

	def remove_from_model(self, destructive=False):
		"""Removes the metabolite association from self.massmodel

		The change is reverted upon exit when using the massmodel as a context.

		Similar to the method in cobra.core.metabolite

		Parameters
		----------
		destructive : bool
			If False then the metabolite is removed from all associated
			reactions.  If True then all associated reactions are removed from
			the MassModel.
		"""
		return self._model.remove_metabolites(self, destructive)

	# Shorthands
	@property
	def ic(self):
		"""Shorthand getter for initial condition"""
		return self._initial_condition

	@ic.setter
	def ic(self, value):
		"""Shorthand setter for initial condition"""
		self.initial_condition = value

	@property
	def gf(self):
		"""Shorthand getter for Gibb's energy of formation"""
		return self._gibbs_formation

	@gf.setter
	def gf(self, value):
		"""Shorthand setter for Gibb's energy of formation"""
		self.gibbs_formation = value

	# Compatibility functions
	def to_cobra_metabolite(self, cobra_id=None):
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
		try:
			from cobra.core.metabolite import Metabolite
		except:
			raise ImportError("Failed to import the Metabolite Object from \
				cobra.core.metabolite. Ensure cobra is installed properly")

		if cobra_id == None:
			cobra_id = self.id + "_cobra"

		new_cobra_metab = Metabolite(id=cobra_id, name=self.name,
								formula=self._formula, charge=self._charge,
								compartment=self._compartment)

		return new_cobra_metab

	def from_cobra_metabolite(CobraMetabolite=None, mass_id=None):
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

		try:
			from cobra.core.metabolite import Metabolite
		except:
			raise ImportError("Failed to import the Metabolite Object from \
				cobra.core.metabolite. Ensure cobra is installed properly")

		if CobraMetabolite == None:
			warn("No cobra Metabolite Object was given")
			return None

		if not isinstance(CobraMetabolite, Metabolite):
			raise TypeError("Metabolite must be a cobra Metabolite Object")

		if mass_id == None:
			mass_id = CobraMetabolite.id + "_mass"

		new_mass_metab = MassMetabolite(id=mass_id, name=CobraMetabolite.name,
								formula=CobraMetabolite.formula,
								charge=CobraMetabolite.charge,
								compartment=CobraMetabolite.compartment)
		return new_mass_metab

	# HTML representation
	def _repr_html_(self):
		return """
		<table>
			<tr>
				<td><strong>MassMetabolite identifier</strong></td>
				<td>{id}</td>
			</tr><tr>
				<td><strong>Name</strong></td>
				<td>{name}</td>
			</tr><tr>
				<td><strong>Memory address</strong></td>
				<td>{address}</td>
			</tr><tr>
				<td><strong>Formula</strong></td>
				<td>{formula}</td>
			</tr><tr>
				<td><strong>Compartment</strong></td>
				<td>{compartment}</td>
			</tr><tr>
				<td><strong>Initial Condition</strong></td>
				<td>{ic}</td>
			</tr><tr>
				<td><strong>Gibbs formation energy</strong></td>
				<td>{gibbs}</td>
			</tr><tr>
				<td><strong>In {n_reactions} reaction(s)</strong></td>
				<td>{reactions}</td>
			</tr>
		<table>""".format(id=self.id, name=self.name, formula=self._formula,
							address='0x0%x' % id(self),
							compartment=self._compartment,
							ic=self._initial_condition,
							gibbs= self._gibbs_formation,
							n_reactions=len(self.reactions),
							reactions='. '.join(r.id for r in self.reactions))

# Dictionary of elements and their corresponding atomic mass
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
