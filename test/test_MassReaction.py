### Ask about type for _forward_rate_constant/_reverse_rate_constant
### They seem to be strings but also are inititlaized as numbers in __init__

### Ask about gene_reaction_rule, see what's up with exceptions

### Ask about shorthand S. What is it a shorthand for?

### Ask regarding 2 assertion errors below (about kr)
"""Unit test for MassReaction Class

Purpose: 
	To Test MassReaction attributes, properties, and methods 
	to ensure they work as intended and raise intended
	exceptions and warnings where appropriate


Intended Processes for Attributes:
	sym_kf:
		Should be a string or None, raise TypeError if not
	sym_kr:
		Should be a string or None, raise TypeError if not
	sym_Keq:
		Should be a string or None, raise TypeError if not
	_reversiblity(getter/setter):
		Should be a bool, raise TypeError if not
	subsystem:
		Should be a string or None, raise TypeError if not
	_forward_rate_constant(getter/setter):
		
	_reverse_rate_constant(getter/setter):
		
	_equilibrium_constant(getter/setter):
		
	_gibbs_formation:
		Should be a float or None, raise TypeError if not
	_rate_law:
		
	_rate_law_expression:
		
	_metabolites:
		Should be a dict, raise TypeError if not
	metabolites(getter):
		Should be a read-only shallow copy, raise TypeError if not
	_compartments(getter):
		Should be a list, raise TypeError if not
	_model:
		
	_genes:
		Should be a set, raise TypeError if not
	genes(getter)
		Should be a frozenset, raise TypeError if not
	_gene_reaction_rule(getter/setter):
		Should be a string, raise TypeError if not
		Context Manager should be implemented, raise UserWarning if not
		
	rate_constants(getter):
		Should be a list [kf, kr], raise TypeError if not
	reactants(getter):
		Should be a list, raise TypeError if not
	products(getter):
		Should be a list, raise TypeError if not
	stoichiometry(getter):
		Should be a list, raise TypeError if not
	forward_rate(getter):
		Should be a string, raise TypeError if not
	forward_rate_expr(getter):
		Should be a sympy expression, raise TypeError if not
	reverse_rate(getter):
		Should be a string, raise TypeError if not
	reverse_rate_expr(getter):
		Should be a sympy expression, raise TypeError if not
	rate_law(getter):
		Should be a string, raise TypeError if not
	rate_law_expr(getter):
		Shoulbe be a sympy expression, raise TypeError if not
	model(getter):
		Should be a MassModel object, raise TypeError if not
	reaction(getter/setter):
		Should be a string, raise TypeError if not
	boundary(getter):
		Should be a bool, raise TypeError if not
	gene_name_reaction_rule(getter):
		Should be a string, raise TypeError if not
	functional(getter):
		Should be a bool, raise TypeError if not


Intended Process for Methods:
	remove_from_model:
		Function the same ways as in cobrapy
	get_coefficient:
		Similar to cobrapy
	get_coefficients:
		Identical to cobrapy
	add_metabolites:
		Similar to cobrapy
	subtract_metabolites:
		Similar to cobrapy
	_set_id_with_model:
		Similar to cobrapy
	_update_awareness:
		
	build_reaction_string:
		Similar to cobrapy
	check_mass_balance:
		Identical to cobrapy
	get_compartments:
		Identical to cobrapy
	build_reaction_from_string:
		Identical to cobrapy
	_associate_gene:
		Identical to cobrapy
	_disassociate_gene:
		Identical to cobrapy
	knock_out:
		Similar to cobrapy
	_repr_html_:
		returns nice table format
	to_cobra_reaction:
		Function to generate cobra Reaction from MassReaction Object
		Should load cobra.core.metabolite, raise ImportError if not
	from_cobra_metabolite:
		Function to generate MassReaction from cobra Reaction Object
		Should take in cobra Reaction Object, raise TypeError if not
		Should take in non-None cobra Reaction Object, raise UserWarning if not


Intended Processes for Shorthands:
	kf:
		Getter should return _forward_rate_constant
		Setter should implement forward_rate_constant.setter
	kr:
		Getter should return _reverse_rate_constant
		Setter should implement reverse_rate_constant.setter
	Keq:
		Getter should return _equilibrium_constant
		Setter should implement equilibrium_constant.setter
	S:
		
Python Magic Methods also tested as appropriate
"""

### Import necessary packages
import unittest, doctest, pytest
### FIXME: replace ipy magic method with proper import from masspy
%run ./MassReaction_Class.ipynb
%run ./MassMetabolite_Class.ipynb

### Testing class begins here
class test_MassReaction(unittest.TestCase, MassReaction):
	
	### setUp/tearDown
	def setUp(self):
		
		### reversibility and other relevents
		self.react1 = MassReaction("TestIDr1") ## basic react for testing
		self.react2 = MassReaction("TestIDr2", reversibility=False) ## react with irreversibility
		
		### compartment and other relevents
		self.react3 = MassReaction("TestIDr3", subsystem = "Techoic Acid Biosynthesis") ## react w/subsystem
		self.react4 = MassReaction("TestIDr4", subsystem = "Transport/Exchange") ## react w/different subsystem
		
		### cobra Reaction (for compatibility functions)
		from cobra.core.reaction import Reaction
		self.cr = Reaction("TestIDr_cobra")
	
	def tearDown(self):
		del self.react1
		del self.react2
		
		del self.react3
		del self.react4
		
		del self.cr
	
	### Test IsInstance
	def test_IsInstance(self):
		self.assertIsInstance(self.react1, MassReaction)
		self.assertIsInstance(self.react2, MassReaction)
		self.assertIsInstance(self.react3, MassReaction)
		self.assertIsInstance(self.react4, MassReaction)
	
	### Test Equal/NotEqual
	def test_reversibility(self):
		expected_value = True
		self.assertEqual(self.react1.reversibility, expected_value)
		self.assertNotEqual(self.react2.reversibility, expected_value)
		self.assertEqual(self.react3.reversibility, expected_value)
		self.assertEqual(self.react4.reversibility, expected_value)
		
		expected_value = False
		self.assertNotEqual(self.react1.reversibility, expected_value)
		self.assertEqual(self.react2.reversibility, expected_value)
		self.assertNotEqual(self.react3.reversibility, expected_value)
		self.assertNotEqual(self.react4.reversibility, expected_value)
	
	def test_forward_rate_constant(self):
		expected_value = "kf_TestIDr1"
		self.assertEqual(self.react1.forward_rate_constant, expected_value)
		self.assertNotEqual(self.react2.forward_rate_constant, expected_value)
		self.assertNotEqual(self.react3.forward_rate_constant, expected_value)
		self.assertNotEqual(self.react4.forward_rate_constant, expected_value)
		
		expected_value = "kf_TestIDr2"
		self.assertNotEqual(self.react1.forward_rate_constant, expected_value)
		self.assertEqual(self.react2.forward_rate_constant, expected_value)
		self.assertNotEqual(self.react3.forward_rate_constant, expected_value)
		self.assertNotEqual(self.react4.forward_rate_constant, expected_value)
		
		expected_value = "kf_TestIDr3"
		self.assertNotEqual(self.react1.forward_rate_constant, expected_value)
		self.assertNotEqual(self.react2.forward_rate_constant, expected_value)
		self.assertEqual(self.react3.forward_rate_constant, expected_value)
		self.assertNotEqual(self.react4.forward_rate_constant, expected_value)
		
		expected_value = "kf_TestIDr4"
		self.assertNotEqual(self.react1.forward_rate_constant, expected_value)
		self.assertNotEqual(self.react2.forward_rate_constant, expected_value)
		self.assertNotEqual(self.react3.forward_rate_constant, expected_value)
		self.assertEqual(self.react4.forward_rate_constant, expected_value)
	
	def test_kf(self):
		expected_value = "kf_TestIDr1"
		self.assertEqual(self.react1.kf, expected_value)
		self.assertNotEqual(self.react2.kf, expected_value)
		self.assertNotEqual(self.react3.kf, expected_value)
		self.assertNotEqual(self.react4.kf, expected_value)
		
		expected_value = "kf_TestIDr2"
		self.assertNotEqual(self.react1.kf, expected_value)
		self.assertEqual(self.react2.kf, expected_value)
		self.assertNotEqual(self.react3.kf, expected_value)
		self.assertNotEqual(self.react4.kf, expected_value)
		
		expected_value = "kf_TestIDr3"
		self.assertNotEqual(self.react1.kf, expected_value)
		self.assertNotEqual(self.react2.kf, expected_value)
		self.assertEqual(self.react3.kf, expected_value)
		self.assertNotEqual(self.react4.kf, expected_value)
		
		expected_value = "kf_TestIDr4"
		self.assertNotEqual(self.react1.kf, expected_value)
		self.assertNotEqual(self.react2.kf, expected_value)
		self.assertNotEqual(self.react3.kf, expected_value)
		self.assertEqual(self.react4.kf, expected_value)
	
	def test_reverse_rate_constant(self):
		expected_value = "kr_TestIDr1"
		self.assertEqual(self.react1.reverse_rate_constant, expected_value)
		self.assertNotEqual(self.react2.reverse_rate_constant, expected_value)
		self.assertNotEqual(self.react3.reverse_rate_constant, expected_value)
		self.assertNotEqual(self.react4.reverse_rate_constant, expected_value)
		
		expected_value = 0
		self.assertNotEqual(self.react1.reverse_rate_constant, expected_value)
		self.assertEqual(self.react2.reverse_rate_constant, expected_value)
		self.assertNotEqual(self.react3.reverse_rate_constant, expected_value)
		self.assertNotEqual(self.react4.reverse_rate_constant, expected_value)
		
		expected_value = "kr_TestIDr3"
		self.assertNotEqual(self.react1.reverse_rate_constant, expected_value)
		self.assertNotEqual(self.react2.reverse_rate_constant, expected_value)
		self.assertEqual(self.react3.reverse_rate_constant, expected_value)
		self.assertNotEqual(self.react4.reverse_rate_constant, expected_value)
		
		expected_value = "kr_TestIDr4"
		self.assertNotEqual(self.react1.reverse_rate_constant, expected_value)
		self.assertNotEqual(self.react2.reverse_rate_constant, expected_value)
		self.assertNotEqual(self.react3.reverse_rate_constant, expected_value)
		self.assertEqual(self.react4.reverse_rate_constant, expected_value)
	
	def test_kr(self):
		expected_value = "kr_TestIDr1"
		self.assertEqual(self.react1.kr, expected_value)
		self.assertNotEqual(self.react2.kr, expected_value)
		self.assertNotEqual(self.react3.kr, expected_value)
		self.assertNotEqual(self.react4.kr, expected_value)
		
		expected_value = 0
		self.assertNotEqual(self.react1.kr, expected_value)
		self.assertEqual(self.react2.kr, expected_value)
		self.assertNotEqual(self.react3.kr, expected_value)
		self.assertNotEqual(self.react4.kr, expected_value)
		
		expected_value = "kr_TestIDr3"
		self.assertNotEqual(self.react1.kr, expected_value)
		self.assertNotEqual(self.react2.kr, expected_value)
		self.assertEqual(self.react3.kr, expected_value)
		self.assertNotEqual(self.react4.kr, expected_value)
		
		expected_value = "kr_TestIDr4"
		self.assertNotEqual(self.react1.kr, expected_value)
		self.assertNotEqual(self.react2.kr, expected_value)
		self.assertNotEqual(self.react3.kr, expected_value)
		self.assertEqual(self.react4.kr, expected_value)
	
	def test_equilibrium_constant(self):
		expected_value = "Keq_TestIDr1"
		self.assertEqual(self.react1.equilibrium_constant, expected_value)
		self.assertNotEqual(self.react2.equilibrium_constant, expected_value)
		self.assertNotEqual(self.react3.equilibrium_constant, expected_value)
		self.assertNotEqual(self.react4.equilibrium_constant, expected_value)
		
		expected_value = "Keq_TestIDr2"
		self.assertNotEqual(self.react1.equilibrium_constant, expected_value)
		self.assertEqual(self.react2.equilibrium_constant, expected_value)
		self.assertNotEqual(self.react3.equilibrium_constant, expected_value)
		self.assertNotEqual(self.react4.equilibrium_constant, expected_value)
		
		expected_value = "Keq_TestIDr3"
		self.assertNotEqual(self.react1.equilibrium_constant, expected_value)
		self.assertNotEqual(self.react2.equilibrium_constant, expected_value)
		self.assertEqual(self.react3.equilibrium_constant, expected_value)
		self.assertNotEqual(self.react4.equilibrium_constant, expected_value)
		
		expected_value = "Keq_TestIDr4"
		self.assertNotEqual(self.react1.equilibrium_constant, expected_value)
		self.assertNotEqual(self.react2.equilibrium_constant, expected_value)
		self.assertNotEqual(self.react3.equilibrium_constant, expected_value)
		self.assertEqual(self.react4.equilibrium_constant, expected_value)
	
	def test_Keq(self):
		expected_value = "Keq_TestIDr1"
		self.assertEqual(self.react1.Keq, expected_value)
		self.assertNotEqual(self.react2.Keq, expected_value)
		self.assertNotEqual(self.react3.Keq, expected_value)
		self.assertNotEqual(self.react4.Keq, expected_value)
		
		expected_value = "Keq_TestIDr2"
		self.assertNotEqual(self.react1.Keq, expected_value)
		self.assertEqual(self.react2.Keq, expected_value)
		self.assertNotEqual(self.react3.Keq, expected_value)
		self.assertNotEqual(self.react4.Keq, expected_value)
		
		expected_value = "Keq_TestIDr3"
		self.assertNotEqual(self.react1.Keq, expected_value)
		self.assertNotEqual(self.react2.Keq, expected_value)
		self.assertEqual(self.react3.Keq, expected_value)
		self.assertNotEqual(self.react4.Keq, expected_value)
		
		expected_value = "Keq_TestIDr4"
		self.assertNotEqual(self.react1.Keq, expected_value)
		self.assertNotEqual(self.react2.Keq, expected_value)
		self.assertNotEqual(self.react3.Keq, expected_value)
		self.assertEqual(self.react4.Keq, expected_value)
	
	def test_rate_constants(self):
		expected_value = ["kf_TestIDr1", "kr_TestIDr1"]
		self.assertEqual(self.react1.rate_constants, expected_value)
		self.assertNotEqual(self.react2.rate_constants, expected_value)
		self.assertNotEqual(self.react3.rate_constants, expected_value)
		self.assertNotEqual(self.react4.rate_constants, expected_value)
		
		expected_value = ["kf_TestIDr2", 0]
		self.assertNotEqual(self.react1.rate_constants, expected_value)
		self.assertEqual(self.react2.rate_constants, expected_value)
		self.assertNotEqual(self.react3.rate_constants, expected_value)
		self.assertNotEqual(self.react4.rate_constants, expected_value)
		
		expected_value = ["kf_TestIDr3", "kr_TestIDr3"]
		self.assertNotEqual(self.react1.rate_constants, expected_value)
		self.assertNotEqual(self.react2.rate_constants, expected_value)
		self.assertEqual(self.react3.rate_constants, expected_value)
		self.assertNotEqual(self.react4.rate_constants, expected_value)
		
		expected_value = ["kf_TestIDr4", "kr_TestIDr4"]
		self.assertNotEqual(self.react1.rate_constants, expected_value)
		self.assertNotEqual(self.react2.rate_constants, expected_value)
		self.assertNotEqual(self.react3.rate_constants, expected_value)
		self.assertEqual(self.react4.rate_constants, expected_value)
	
	### ### Requires add_metabolites and subract_metabolites
	
	def test_metabolites(self):
		
		### Testing with add_metabolites
		metab1 = MassMetabolite("TestIDm1", formula="CH4")
		react1.add_metabolites({metab1: 5})
		
		expected_value = {metab1: 5}
		self.assertEqual(self.react1.metabolites, expected_value)
		self.assertNotEqual(self.react2.metabolites, expected_value)
		self.assertNotEqual(self.react3.metabolites, expected_value)
		self.assertNotEqual(self.react4.metabolites, expected_value)
		
		### ### Testing combine=True
		react1.add_metabolites({metab1: 10})
		expected_value = {metab1: 15}
		self.assertEqual(self.react1.metabolites, expected_value)
		self.assertNotEqual(self.react2.metabolites, expected_value)
		self.assertNotEqual(self.react3.metabolites, expected_value)
		self.assertNotEqual(self.react4.metabolites, expected_value)
		
		### ### Testing combine=False
		
		### ### Testing reversibly=False
	
	#def test_reactants(self):
	
	#def test_products(self):
	
	#def test_stoichiometry(self):
	
	
	### Test Exceptions and Warnings
	
	### ### Parameter inputs (from __init__)
	def test_name_not_string_has_typeerror_int(self):
		with self.assertRaisesRegex(TypeError, "name must be a string type") as t:
			self.react_temp = MassReaction(id="TestIDr_temp", name=4)
			del self.react_temp
	
	def test_name_not_string_has_typeerror_float(self):
		with self.assertRaisesRegex(TypeError, "name must be a string type") as t:
			self.react_temp = MassReaction(id="TestIDr_temp", name=4.5)
			del self.react_temp
	
	def test_name_not_string_has_typeerror_bool(self):
		with self.assertRaisesRegex(TypeError, "name must be a string type") as t:
			self.react_temp = MassReaction(id="TestIDr_temp", name=True)
			del self.react_temp
	
	def test_subsystem_not_string_has_typeerror_int(self):
		with self.assertRaisesRegex(TypeError, "subsystem must be a string type") as t:
			self.react_temp = MassReaction(id="TestIDr_temp", subsystem=4)
			del self.react_temp
	
	def test_subsystem_not_string_has_typeerror_float(self):
		with self.assertRaisesRegex(TypeError, "subsystem must be a string type") as t:
			self.react_temp = MassReaction(id="TestIDr_temp", subsystem=4.5)
			del self.react_temp
	
	def test_subsystem_not_string_has_typeerror_bool(self):
		with self.assertRaisesRegex(TypeError, "subsystem must be a string type") as t:
			self.react_temp = MassReaction(id="TestIDr_temp", subsystem=True)
			del self.react_temp
	
	def test_reversibility_not_bool_has_typeerror_int(self):
		with self.assertRaisesRegex(TypeError, "reversibility must be a boolean") as t:
			self.react_temp = MassReaction(id="TestID_temp", reversibility=4)
			del self.react_temp
	
	def test_reversibility_not_bool_has_typeerror_float(self):
		with self.assertRaisesRegex(TypeError, "reversibility must be a boolean") as t:
			self.react_temp = MassReaction(id="TestIDr_temp", reversibility=4.5)
			del self.react_temp
	
	def test_reversibility_not_bool_has_typeerror_string(self):
		with self.assertRaisesRegex(TypeError, "reversibility must be a boolean") as t:
			self.react_temp = MassReaction(id="TestIDr_temp", reversibility="True")
			del self.react_temp
	
	### ### Properties and Attributes
	def test_reverse_rate_constant_with_irreversible_reaction_has_userwarning(self):
		with self.assertRaisesRegex(UserWarning, "No reverse rate constant for irreversible reactions") as w:
			self.react2.reverse_rate_constant
	
	def test_kr_with_irreversible_reaction_has_userwarning(self):
		with self.assertRaisesRegex(UserWarning, "No reverse rate constant for irreversible reactions") as w:
			self.react2.kr
	
if __name__ == '__main__':
	unittest.main(argv=['first-arg-is-ignored'], exit=False)