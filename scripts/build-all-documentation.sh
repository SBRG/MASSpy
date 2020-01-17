#!/bin/bash

bash_source_dir=$( cd "$(dirname "${BASH_SOURCE[0]}")" ; pwd -P )
doc_builder_dir="$bash_source_dir/../documentation_builder"

sh $bash_source_dir/clean.sh
sh $bash_source_dir/make-documentation-notebooks.sh
sh $bash_source_dir/make-additional-examples.sh
sh $bash_source_dir/make-sb2-textbook.sh

cd $doc_builder_dir

make clean
make html