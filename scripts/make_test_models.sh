#!/bin/sh
parent_path=$( cd "$(dirname "${BASH_SOURCE[0]}")" ; pwd -P )

execute_notebooks="$parent_path/execute_notebooks.py"
replace_package_models="$parent_path/copy_models_to_test.py"

python $execute_notebooks -v -c all 
python $replace_package_models -v -n all