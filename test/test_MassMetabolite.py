"""Unit test for MassMetabolite Class

Purpose:
    To Test MassMetabolite attributes, properties, and methods
    to ensure they work as intended and raise intended
    exceptions and warnings where appropriate


Intended Processes for Attributes/Properties:
    _formula:
        Should be a string or None, raise TypeError if not
    _charge:
        Should be a float or None, raise TypeError if not
    _compartment:
        Should be a string or None, raise TypeError if not
    elements:
        Should be a dictionary, raise TypeError if not
    _initial_condition:
        Should be a float or None, raise TypeError if not
        Should be >= 0, raise ValueError if not
    _gibbs_formation:
        Should be a float or None, raise TypeError if not
    formula_weight:
        Should be a float, raise TypeError if not
        Should be > 0, raise ValueError if not


Intended Process for Methods:
    _set_id_with_model:
        Function the same way as in cobrapy
    remove_from_model:
        Function the same ways as in cobrapy


Intended Processes for Shorthands:
    ic:
        Getter should return _initial_condition
        Setter should implement initial_condition.setter
    gf:
        Getter should return _gibbs_formation
        Setter should implement gibb_formation.setter"""

### Import necessary packages
from __future__ import absolute_import

import unittest, doctest, pytest, warnings

from mass import MassMetabolite, MassReaction, MassModel ## won't work yet

### Testing class begins here
class test_MassMetabolite(unittest.TestCase, MassMetabolite):

    ### setUp/tearDown
    def setUp(self):

        ### formula, elements, and formula_weight
        self.metab1 = MassMetabolite('TestID1', formula='CH4') ## basic metab for testing
        self.metab2 = MassMetabolite('TestID2', formula='CH3') ## metab w/different formula
        self.metab3 = MassMetabolite('TestID3', formula='C6H12.5O6') ## metab with decimal formula

        ### charge
        self.metab4 = MassMetabolite('TestID4', charge= 2) ## basic metab for testing
        self.metab5 = MassMetabolite('TestID5', charge= -2) ## basic metab w/different charge

        ### compartment
        self.metab6 = MassMetabolite('TestID6', compartment="c") ## bastic metab for testing
        self.metab7 = MassMetabolite('TestID7', compartment="e") ## basic metab with differt compartment

        ### initial_condition & gibbs_formation
        self.metab8 = MassMetabolite('TestID8')

    def tearDown(self):
        del self.metab1
        del self.metab2
        del self.metab3

        del self.metab4
        del self.metab5

        del self.metab6
        del self.metab7


    ### Test IsInstance
    def test_IsIntance(self):
        self.assertIsInstance(self.metab1, MassMetabolite)
        self.assertIsInstance(self.metab2, MassMetabolite)
        self.assertIsInstance(self.metab3, MassMetabolite)
        self.assertIsInstance(self.metab4, MassMetabolite)
        self.assertIsInstance(self.metab5, MassMetabolite)
        self.assertIsInstance(self.metab6, MassMetabolite)
        self.assertIsInstance(self.metab7, MassMetabolite)
        self.assertIsInstance(self.metab8, MassMetabolite)


    ### Test Equal/NotEqual
    def test_formulas(self):

        expected_result = 'CH4'
        self.assertEqual(self.metab1.formula, expected_result)
        self.assertNotEqual(self.metab2.formula, expected_result)
        self.assertNotEqual(self.metab3.formula, expected_result)

        expected_result = 'CH3'
        self.assertNotEqual(self.metab1.formula, expected_result)
        self.assertEqual(self.metab2.formula, expected_result)
        self.assertNotEqual(self.metab3.formula, expected_result)

        expected_result = 'C6H12.5O6'
        self.assertNotEqual(self.metab1.formula, expected_result)
        self.assertNotEqual(self.metab2.formula, expected_result)
        self.assertEqual(self.metab3.formula, expected_result)

    def test_charge(self):

        expected_result = 2
        self.assertEqual(self.metab4.charge, expected_result)
        self.assertNotEqual(self.metab5.charge, expected_result)

        expected_result = -2
        self.assertNotEqual(self.metab4.charge, expected_result)
        self.assertEqual(self.metab5.charge, expected_result)

    def test_compartment(self):

        expected_result = "c"
        self.assertEqual(self.metab6.compartment, expected_result)
        self.assertNotEqual(self.metab7.compartment, expected_result)

        expected_result = "e"
        self.assertNotEqual(self.metab6.compartment, expected_result)
        self.assertEqual(self.metab7.compartment, expected_result)

    def test_elements(self):
        expected_result = {'C': 1, 'H': 4}
        self.assertEqual(self.metab1.elements, expected_result)
        self.assertNotEqual(self.metab2.elements, expected_result)
        self.assertNotEqual(self.metab3.elements, expected_result)

        expected_result = {'C': 1, 'H': 3}
        self.assertNotEqual(self.metab1.elements, expected_result)
        self.assertEqual(self.metab2.elements, expected_result)
        self.assertNotEqual(self.metab3.elements, expected_result)

        expected_result = {'C': 6, 'H': 12.5, 'O': 6}
        self.assertNotEqual(self.metab1.elements, expected_result)
        self.assertNotEqual(self.metab2.elements, expected_result)
        self.assertEqual(self.metab3.elements, expected_result)

    def test_initial_condition(self):
        self.metab8.initial_condition = 5

        expected_result = 5
        self.assertEqual(self.metab8.initial_condition, expected_result)

        expected_result = 10
        self.assertNotEqual(self.metab8.initial_condition, expected_result)

    def test_ic(self):
        self.metab8.ic = 5

        expected_result = 5
        self.assertEqual(self.metab8.ic, expected_result)

        expected_result = 10
        self.assertNotEqual(self.metab8.ic, expected_result)

    def test_gibbs_formation(self):
        self.metab8.gibbs_formation = 100

        expected_result = 100
        self.assertEqual(self.metab8.gibbs_formation, expected_result)

        expected_result = 1000
        self.assertNotEqual(self.metab8.gibbs_formation, expected_result)

    def test_gf(self):
        self.metab8.gf = 100

        expected_result = 100
        self.assertEqual(self.metab8.gf, expected_result)

        expected_result = 1000
        self.assertNotEqual(self.metab8.gf, expected_result)

    def test_formula_weight(self):

        expected_result = 16.04246
        self.assertEqual(self.metab1.formula_weight, expected_result)
        self.assertNotEqual(self.metab2.formula_weight, expected_result)
        self.assertNotEqual(self.metab3.formula_weight, expected_result)

        expected_result = 15.03452
        self.assertNotEqual(self.metab1.formula_weight, expected_result)
        self.assertEqual(self.metab2.formula_weight, expected_result)
        self.assertNotEqual(self.metab3.formula_weight, expected_result)

        expected_result = 180.65985
        self.assertNotEqual(self.metab1.formula_weight, expected_result)
        self.assertNotEqual(self.metab2.formula_weight, expected_result)
        self.assertEqual(self.metab3.formula_weight, expected_result)


    ### Test Exceptions and Warnings
    @unittest.skip("Not added yet")
    def test_formula_not_string_or_none_has_typeerror(self):
        with self.assertRaises(TypeError) as t:
            metab_temp = MassMetabolite('TestID_temp', formula=2)

        self.assertRaisesRegex(TypeError, "formula must be a string")
        del metab_temp

    @unittest.skip("Not added yet")
    def test_charge_not_float_or_none_has_typeerrory(self):
        with self.assertRaises(TypeError) as t:
            metab_temp = MassMetabolite('TestID_temp', charge="-2")

        self.assertRaisesRegex(TypeError, "charge must be a float")
        del metab_temp

    @unittest.skip("Not added yet")
    def test_compartment_not_string_or_none_has_typeerror(self):
        with self.assertRaises(TypeError) as t:
            metab_temp = MassMetabolite('TestID_temp', compartment=3)

        self.assertRaisesRegex(TypeError, "compartment must a string")
        del metab_temp

    def test_elements_asterisk_formula_has_userwarning(self):
        with self.assertWarns(UserWarning) as w:
            metab_temp = MassMetabolite('TestID_temp', formula="CoCl2*6H2O")
            metab_temp.elements

        self.assertWarnsRegex(UserWarning, "invalid character '*' found in formula")
        del metab_temp

    def test_elements_parentheses_formula_has_userwarning(self):
        """Test is a warning for now.
        May consider returning an error in future"""

        with self.assertWarns(UserWarning) as w:
            metab_temp = MassMetabolite('TestID_temp', formula="(CH3)3N+CH2CH2OCOCH3")
            metab_temp.elements

        self.assertWarnsRegex(UserWarning, "invalid formula (has parenthesis)")
        del metab_temp

    def test_elements_noninteger_formula_has_userwarning(self):
        with self.assertWarns(UserWarning) as w:
            self.metab3.elements

        self.assertWarnsRegex(UserWarning, "is not an integer")

    def test_initial_condition_not_float_or_none_has_typeerror(self):
        with self.assertRaises(TypeError) as t:
            #self.metab8.initial_condition = "1.4"
            self.metab8.initial_condition = 0.9

        self.assertRaisesRegex(TypeError, "initial condition must be a number")

    def test_initial_condition_negative_value_has_valueerror(self):
        with self.assertRaises(ValueError) as v:
            self.metab1.initial_condition = -1

        self.assertEqual(str(v.exception), "Initial condition must be a non-negative number")
