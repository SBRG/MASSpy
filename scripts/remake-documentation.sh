#!/bin/bash

bash_source_dir=$( cd "$(dirname "${BASH_SOURCE[0]}")" ; pwd -P )
doc_builder_dir="$bash_source_dir/../documentation_builder"

cd $doc_builder_dir
cd notebooks

DOCUMENTATION_NOTEBOOKS=("getting_started_with_masspy"\
                         "constructing_models"\
                         "reading_writing_models"\
                         "dynamic_simulation"\
                         "plot_visualization"\
                         "enzyme_modules"\
                         "thermo_concentrations"\
                         "ensemble_modeling"\
                         # "network_visualization"\
                         "quality_assurance"\
                         "cobra_to_mass"\
                         "faq"\
                         )

for doc_nb in "${DOCUMENTATION_NOTEBOOKS[@]}"; do
    jupyter nbconvert --to notebook --execute $doc_nb.ipynb  --inplace --ExecutePreprocessor.timeout=-1
done

# Execute additional example notebooks
jupyter nbconvert --to notebook --execute ./workflows/*.ipynb  --inplace --ExecutePreprocessor.timeout=-1
jupyter nbconvert --to notebook --execute ./advanced_visualization/*.ipynb  --inplace --ExecutePreprocessor.timeout=-1

rm -rf test_textbook.*
rm -rf SB2_Glycolysis.*

cd $doc_builder_dir

make clean
make html
make latex
