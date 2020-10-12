#!/bin/sh

set -e

docs_dir=$( cd "$(dirname "$0")/../../docs" ; pwd -P )

# Default arguments if none provided
to='notebook' 
if [ "$#" -eq 0 ] ; then \
    additional_dir='additional'
    sb2_chapters_dir='education/sb2/chapters'
    sb2_models_dir='education/sb2/model_construction'
    tutorial_dir='tutorials'
    visualization_dir='gallery/visualization'
    workflow_dir='gallery/workflows'
    inplace='--inplace'
else
    # Parse  arguments
    while [ "$#" -ne 0 ] ; do
        case $1 in
            # Determine which notebooks to execute
            -a | --all )
                additional_dir='additional'
                sb2_chapters_dir='education/sb2/chapters'
                sb2_models_dir='education/sb2/model_construction'
                tutorial_dir='tutorials'
                visualization_dir='gallery/visualization'
                workflow_dir='gallery/workflows'
                ;;
            -e | --education | --edu )
                sb2_chapters_dir='education/sb2/chapters'
                sb2_models_dir='education/sb2/model_construction'                ;;
            -g | --gallery | --gal )
                visualization_dir='gallery/visualization'
                workflow_dir='gallery/workflows'
                ;;
            -t | --tutorials )
                tutorial_dir='tutorials'
                ;;
            -m | --miscellaneous | --misc )
                additional_dir='additional'
                ;;
            # Other arguments for execution
            -i | --inplace )
                inplace='--inplace'
                ;;
            --to )
                shift && to=$1
                if echo $to | grep -Ev "^(webpdf)$" ; then \
                    allow_chromium_download="--allow-chromium-download"
                elif echo $to | grep -Ev "^(html|webpdf|notebook)$" ; then \
                    echo "Bad '--to' value, defaulting to notebook" 
                    to='notebook'
                fi
                ;;
            --allow-errors )
                allow_errors='--allow-errors'
                ;;
        esac
        shift 
    done 
fi

options=$( \
    echo --to $to \
         --ExecutePreprocessor.timeout=-1 \
         --Application.log_level='CRITICAL' \
         --ExecutePreprocessor.store_widget_state=True \
         $allow_errors $output_dir $allow_chromium_download $inplace )

for directory in $additional_dir \
                 $sb2_chapters_dir \
                 $sb2_models_dir \
                 $tutorial_dir \
                 $visualization_dir \
                 $workflow_dir ; do \
    if [ $directory ] ; then \
        for notebook in $(find $docs_dir/$directory -type f -name "*.ipynb" | sort -u ); do \
            echo "Executing $( basename $notebook )"
            jupyter nbconvert --execute $notebook $options 
        done
    fi 
    done