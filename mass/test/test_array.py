"""Unit test for MassModel Class

Purpose: 
    To test array methods 
    to ensure they work as intended and raise intended
    exceptions and warnings where appropriate


Intended Process for Methods:
Intended Processes for Shorthands:
"""

import unittest
import numpy as np
import pandas as pd
import mass

from scipy.sparse import dok_matrix, lil_matrix
from six import string_types, iteritems
from mass import MassMetabolite, MassReaction, MassModel
from mass.util import array
from mass.test.data import glycolysis

### Testing class begins here
class test_array(unittest.TestCase, MassModel):
    
    ### setUp/tearDown
    def setUp(self):
        self.glycolysis = glycolysis
    
    def tearDown(self):
        del self.glycolysis
    
    ### Test IsInstance
    def test_IsInstance(self):
        self.assertIsInstance(
            array.create_stoichiometric_matrix(self.glycolysis),
            np.ndarray
        )
        self.assertIsInstance(
            array.create_stoichiometric_matrix(self.glycolysis, "dense"),
            np.ndarray
        )
        self.assertIsInstance(
            array.create_stoichiometric_matrix(self.glycolysis, "dok"), 
            dok_matrix
        )
        self.assertIsInstance(
            array.create_stoichiometric_matrix(self.glycolysis, "lil"),
            lil_matrix
        )
        self.assertIsInstance(
            array.create_stoichiometric_matrix(self.glycolysis, "dataframe"),
            pd.core.frame.DataFrame
        )
    
    ### Test Equal/NotEqual
    
    ### Test Exceptions and Warnings
    
    ### ### create_stoichiometric_matrix
    def test_create_stoichiometric_matrix_model_not_MassModel_has_typeerror(self):
        with self.assertRaisesRegex(TypeError, "model must be a MassModel") as t:
            not_a_model = MassMetabolite("not_a_model")
            array.create_stoichiometric_matrix(not_a_model)
    
    def test_create_stoichiometric_matrix_update_model_not_bool_has_typeerror(self):
        with self.assertRaisesRegex(TypeError, "update_model must be a bool") as t:
            array.create_stoichiometric_matrix(glycolysis, update_model=1)
    
    def test_create_stoichiometric_matrix_matrix_type_not_proper_has_typeerror(self):
        with self.assertRaisesRegex(TypeError, "matrix_type must be a string") as t:
            array.create_stoichiometric_matrix(glycolysis, matrix_type=1)
    
    def test_create_stoichiometric_matrix_matrix_type_not_proper_has_valueerror(self):
        with self.assertRaisesRegex(ValueError, "matrix_type must be a string of one of the following"
                                    " types: {'dense', 'dok', 'lil', 'dataframe'}") as v:
            array.create_stoichiometric_matrix(glycolysis, matrix_type="lil 'sac")
    
    ### ### nullspace, left_nullspace, and matrix_rank
    def test_nullspace_A_not_numpy_or_scipy_array_type_has_typeerror(self):
        with self.assertRaisesRegex(TypeError, "Matrix must be one of the following formats: "
                                    "numpy.ndarray, dok_matrix, lil_matrix, or pandas 'DataFrame'") as t:
            array.nullspace(5)
    
    def test_left_nullspace_A_not_numpy_or_scipy_array_type_has_typeerror(self):
        with self.assertRaisesRegex(TypeError, "Matrix must be one of the following formats: "
                                    "numpy.ndarray, dok_matrix, lil_matrix, or pandas 'DataFrame'") as t:
            array.left_nullspace(5)
    
    def test_matrix_rank_A_not_numpy_or_scipy_array_type_has_typeerror(self):
        with self.assertRaisesRegex(TypeError, "Matrix must be one of the following formats: "
                                    "numpy.ndarray, dok_matrix, lil_matrix, or pandas 'DataFrame'") as t:
            array.matrix_rank(5)
    
    ### ### _update_S
    def test__update_S_model_not_MassModel_has_typeerror(self):
        with self.assertRaisesRegex(TypeError, "model must be a MassModel") as t:
            not_a_model = MassMetabolite("not_a_model")
            array._update_S(not_a_model)
    
    def test__update_S_matrix_type_not_proper_has_typeerror(self):
        with self.assertRaisesRegex(TypeError, "matrix_type must be a string") as t:
            array._update_S(glycolysis, matrix_type=1)
    
    def test__update_S_update_model_not_bool_has_typeerror(self):
        with self.assertRaisesRegex(TypeError, "update_model must be a bool") as t:
            array._update_S(glycolysis, update_model=1)
    
    def test__update_S_matrix_type_not_proper_has_valueerror(self):
        with self.assertRaisesRegex(ValueError, "matrix_type must be of one of the following"
                                    " types: {'dense', 'dok', 'lil', 'dataframe'}") as v:
            array._update_S(glycolysis, matrix_type="lil 'sac")

if __name__ == '__main__':
    unittest.main(argv=['first-arg-is-ignored'], exit=False)