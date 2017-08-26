"""Unit test for MassModel Class

Purpose: 
    To Test MassModel attributes, properties, and methods 
    to ensure they work as intended and raise intended
    exceptions and warnings where appropriate


Intended Processes for Attributes:
Intended Process for Methods:
Intended Processes for Shorthands:
"""

import unittest

import numpy as np
import pandas as pd

import mass
from mass import MassMetabolite, MassReaction, MassModel
from mass.test.data import glycolysis
from six import iteritems

import cobra
from cobra import Gene

### Testing class begins here
class test_MassModel(unittest.TestCase, MassModel):
    
    ### setUp/tearDown
    def setUp(self):
        # Reference Model
        self.glycolysis = glycolysis
        
        # Metabolites
        self.x1 = MassMetabolite("x1")
        self.x2 = MassMetabolite("x2")
        self.x3 = MassMetabolite("x3")
        self.x4 = MassMetabolite("x4")
        
        # Reactions
        self.v1 = MassReaction("v1")
        self.v1.add_metabolites({
            self.x1: -1,
            self.x2: 1
        })
        
        self.v2 = MassReaction("v2")
        self.v2.add_metabolites({
            self.x2: -1,
            self.x3: 1
        })
        
        self.v3 = MassReaction("v3")
        self.v3.add_metabolites({
            self.x3: -1,
            self.x4: 1
        })
        
        # Models
        self.s1 = MassModel("s1")
        self.s1.add_reactions([
            self.v1,
            self.v2
        ])
        
        self.s2 = MassModel("s2")
        self.s2.add_reactions([
            self.v1,
            self.v2,
            self.v3
        ])
    
    def tearDown(self):
        del self.glycolysis
        
        del self.s1, self.s2
        del self.v1, self.v2, self.v3
        del self.x1, self.x2, self.x3
    
    ### Test IsInstance
    def test_IsInstance(self):
        return None
    
    def test_IsInstance_to_cobra_reaction(self):
        return None
    
    def test_IsInstance_from_cobra_reaction(self):
        return None
    
    ### Test Equal/NotEqual
    def test_S(self):
        expected_value = np.array(
            [[-1, 0],
             [1, -1],
             [0, 1]]
        )
        self.assertTrue(np.array_equal(self.s1.S, expected_value))
        
        expected_value = np.array(
            [[-1, 0, 0],
             [1, -1, 0],
             [0, 1, -1],
             [0, 0, 1]]
        )
        self.assertTrue(np.array_equal(self.s2.S, expected_value))
    
    def test_rates(self):
        return None
    
    def test_rate_expressions(self):
        return None
    
    def test_odes(self):
        return None
    
    def test_exchanges(self):
        expected_value = []
        
        self.assertEqual(self.s1.exchanges, expected_value)
        self.assertEqual(self.s2.exchanges, expected_value)
        
        expected_value = [
            self.glycolysis.reactions.EX_glc__D,
            self.glycolysis.reactions.EX_AMPIN,
            self.glycolysis.reactions.EX_AMPOUT,
            self.glycolysis.reactions.EX_lac__L,
            self.glycolysis.reactions.EX_pyr,
            self.glycolysis.reactions.EX_h,
            self.glycolysis.reactions.EX_h2o
        ]
        
        self.assertEqual(self.glycolysis.exchanges, expected_value)
    
    @unittest.skip("leave until Zack fixes it")
    def test_get_external_metabolites(self):
        expected_value = []
        
        self.assertEqual(self.s1.get_external_metabolites, expected_value)
        self.assertEqual(self.s2.get_external_metabolites, expected_value)
        
        expected_value = [
            "glc__D_Xt",
            "amp_Xt",
            "lac__L_Xt",
            "pyr",
            "h",
            "h2o"
        ]
        
        self.assertEqual(self.glycolysis.get_external_metabolites, expected_value)
    
    def test_get_metabolite_compartments(self):
        expected_value = set()
        
        self.assertEqual(self.s1.get_metabolite_compartments, expected_value)
        self.assertEqual(self.s2.get_metabolite_compartments, expected_value)
        
        expected_value = {"c"}
        
        self.assertEqual(self.glycolysis.get_metabolite_compartments, expected_value)
    
    def test_get_irreversible_reactions(self):
        expected_value = []
        
        self.assertEqual(self.s1.get_irreversible_reactions, expected_value)
        self.assertEqual(self.s2.get_irreversible_reactions, expected_value)
        
        expected_value = [
            self.glycolysis.reactions.ATPM,
            self.glycolysis.reactions.DM_nadh,
            self.glycolysis.reactions.EX_glc__D,
            self.glycolysis.reactions.EX_AMPIN,
            self.glycolysis.reactions.EX_AMPOUT
        ]
        
        self.assertEqual(self.glycolysis.get_irreversible_reactions, expected_value)
    
    def test_steady_state_fluxes(self):
        expected_value = {}
        
        self.assertEqual(self.s1.steady_state_fluxes, expected_value)
        self.assertEqual(self.s2.steady_state_fluxes, expected_value)
        
        expected_value = {
            self.glycolysis.reactions.HEX1: 1.12,
            self.glycolysis.reactions.PGI: 1.12,
            self.glycolysis.reactions.PFK: 1.12,
            self.glycolysis.reactions.FBA: 1.12,
            self.glycolysis.reactions.TPI: 1.12,
            
            self.glycolysis.reactions.GAPD: 2.24,
            self.glycolysis.reactions.PGK: 2.24,
            self.glycolysis.reactions.PGM: 2.24,
            self.glycolysis.reactions.ENO: 2.24,
            self.glycolysis.reactions.PYK: 2.24,
            
            self.glycolysis.reactions.LDH_L: 2.016,
            
            self.glycolysis.reactions.ADK1: 0,
            self.glycolysis.reactions.ATPM: 2.24,
            
            self.glycolysis.reactions.DM_nadh: 0.224,
            self.glycolysis.reactions.EX_glc__D: 1.12,
            
            self.glycolysis.reactions.EX_AMPIN: 0.014,
            self.glycolysis.reactions.EX_AMPOUT: 0.014,
            
            self.glycolysis.reactions.EX_lac__L: 2.016,
            self.glycolysis.reactions.EX_pyr: 0.224,
            
            self.glycolysis.reactions.EX_h: 2.688,
            self.glycolysis.reactions.EX_h2o: 0
        }
        
        self.assertEqual(self.glycolysis.steady_state_fluxes, expected_value)
    
    ###############################################
    ###############################################
    ######### Integration & Testing below #########
    ### This testing depends on the copy method ###
    ###############################################
    ###############################################
    
    ### UNDER CONSTRUCTION ######
    ### HELP WANTED #############
    ### ROBOTS NEED NOT APPLY ###
    
    ### ### Attributes ### ###
    def test_copy_test_reactions(self):
        self.copy = glycolysis.copy()
        self.assertNotEqual(glycolysis.reactions, self.copy.reactions)
        
        self.assertTrue(
            [rxn.name for rxn in glycolysis.reactions] == 
            [rxn.name for rxn in self.copy.reactions]
        )
        self.assertTrue(
            [rxn.id for rxn in glycolysis.reactions] == 
            [rxn.id for rxn in self.copy.reactions]
        )
        
        del self.copy
    
    def test_copy_test_metabolites(self):
        self.copy = glycolysis.copy()
        self.assertNotEqual(glycolysis.metabolites, self.copy.metabolites)
        
        self.assertTrue(
            [metab.name for metab in glycolysis.metabolites] == 
            [metab.name for metab in self.copy.metabolites]
        )
        self.assertTrue(
            [metab.id for metab in glycolysis.metabolites] == 
            [metab.id for metab in self.copy.metabolites]
        )
        
        del self.copy
    
    def test_copy_test_genes(self):
        #glycolysis.genes = Gene("pseudo")
        #self.assertNotEqual(glycolysis.genes, self.copy.genes)
        return None
    
    def test_copy_test_initial_conditions(self):
        self.copy = glycolysis.copy()
        self.assertNotEqual(glycolysis.initial_conditions, self.copy.initial_conditions)
        
        self.assertTrue(
            {key.name: value for (key, value) in glycolysis.initial_conditions.items()} == 
            {key.name: value for (key, value) in self.copy.initial_conditions.items()}
        )
        self.assertTrue(
            {key.id: value for (key, value) in glycolysis.initial_conditions.items()} == 
            {key.id: value for (key, value) in self.copy.initial_conditions.items()}
        )
        
        del self.copy
    
    def test_copy_test_custom_rates(self):
        return None
    
    def test_copy_test_fixed_concentrations(self):
        return None
    
    def test_copy_test_compartments(self):
        self.copy = glycolysis.copy()
        self.assertEqual(glycolysis.compartments, self.copy.compartments)
        
        del self.copy
    
    def test_copy_test_units(self):
        return None
    ##########################
    
    ### ### Properties ### ###
    def test_copy_test_S(self):
        self.copy = glycolysis.copy()
        self.assertTrue(np.array_equal(glycolysis.S, self.copy.S))
        
        del self.copy
    
    def test_copy_test_rates(self):
        ### Insert rates here
        return None
    
    def test_copy_test_rate_expressions(self):
        ### Insert rate_expressions here
        return None
    
    def test_copy_test_odes(self):
        ### Insert odes here
        return None
    
    def test_copy_test_exchanges(self):
        self.copy = glycolysis.copy()
        self.assertNotEqual(glycolysis.exchanges, self.copy.exchanges)
        
        del self.copy
    
    def test_copy_test_get_external_metabolites(self):
        self.copy = glycolysis.copy()
        self.assertEqual(glycolysis.get_external_metabolites, self.copy.get_external_metabolites)
        
        del self.copy
    
    def test_copy_test_get_metabolite_compartments(self):
        self.copy = glycolysis.copy()
        self.assertEqual(glycolysis.get_metabolite_compartments, self.copy.get_metabolite_compartments)
        
        del self.copy
    
    def test_copy_test_get_irreversible_reactions(self):
        self.copy = glycolysis.copy()
        self.assertNotEqual(glycolysis.get_irreversible_reactions, self.copy.get_irreversible_reactions)
        
        del self.copy
    
    def test_copy_test_steady_state_fluxes(self):
        self.copy = glycolysis.copy()
        self.assertNotEqual(glycolysis.steady_state_fluxes, self.copy.steady_state_fluxes)
        
        self.assertTrue(
            {key.name: value for (key, value) in glycolysis.steady_state_fluxes.items()} == 
            {key.name: value for (key, value) in self.copy.steady_state_fluxes.items()}
        )
        self.assertTrue(
            {key.id: value for (key, value) in glycolysis.steady_state_fluxes.items()} == 
            {key.id: value for (key, value) in self.copy.steady_state_fluxes.items()}
        )
        
        del self.copy
    
    ##########################
    
#     def test_update_S(self):
#         return None
    
#     def test_add_metabolites(self):
#         self.temp_model = glycolysis
#         self.temp_metab = MassMetabolite("temp_metab", name="Haiman")
        
#         self.temp_model.add_metabolites(self.temp_metab)
        
#         assert "Haiman" in [metab.name for metab in self.temp_model.metabolites]
#         assert self.temp_model.S.shape == (21, 21)
        
#         del self.temp_model, self.temp_metab
        
#     def test_add_metabolites_with_initial_conditions(self):
#         self.temp_model2 = glycolysis
#         self.temp_metab2 = MassMetabolite("temp_metab2", name="Haiman")
#         self.temp_metab2.ic = 8328
        
#         self.temp_model2.add_metabolites(self.temp_metab2)
#         expected_value = 8328
        
#         assert "Haiman" in [metab.name for metab in self.temp_model2.metabolites]
#         #assert self.temp_model2.S.shape == (21, 21)
#         print(self.temp_model2.S.shape)
#         self.assertEqual(self.temp_model2.metabolites.temp_metab2.ic, expected_value)
        
    
#     def test_remove_metabolites(self):
#         return None
        
    
    ### Test Exceptions and Warnings

### Run Program ###
if __name__ == '__main__':
    unittest.main(argv=['first-arg-is-ignored'], exit=False)