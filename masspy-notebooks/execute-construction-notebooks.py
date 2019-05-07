# -*- coding: utf-8 -*-
"""TODO Module Docstrings."""
from __future__ import absolute_import

from os import path

from nbconvert.preprocessors import ExecutePreprocessor

import nbformat

ALL_CONSTRUCTION_NOTEBOOKS = [
    "SB2_Glycolysis", "SB2_Pentose Phosphate Pathway",
    "SB2_AMP_Salvage_Network", "SB2_Hemoglobin", "SB2_Textbook_Model",
    "HEX1_EnzymeModule", "PFK_EnzymeModule", "PYK_EnzymeModule",
    "G6PDH2r_EnzymeModule"]

for nb_name in ALL_CONSTRUCTION_NOTEBOOKS:
    # Create path to notebooks
    if "SB2_" in nb_name:
        construction_dir = "SB2-textbook"
    elif "_EnzymeModule" in nb_name:
        construction_dir = "EnzymeModules"
    else:
        continue
    notebook_dir = path.abspath(path.join(
        "test-model-construction", construction_dir))
    notebook_path = path.abspath(path.join(notebook_dir, nb_name + ".ipynb"))
    # Load and execute
    print("Executing construction notebook " + nb_name)
    with open(notebook_path) as f:
        nb = nbformat.read(f, as_version=4)
    ep = ExecutePreprocessor(timeout=600, kernel_name='python3')
    ep.preprocess(nb, {'metadata': {'path': notebook_dir}})
    with open(notebook_path, 'wt') as f:
        nbformat.write(nb, f)

print("Finished executing construction notebooks.")
