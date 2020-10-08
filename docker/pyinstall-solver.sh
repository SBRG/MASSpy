#!/bin/sh

set -e

verbose=0
solver=
while [ "$#" -ne 0 ] ; do
    case $1 in
        --directory)
            shift && directory=$1 ;;
        --solver)
            shift && solver=$1 ;;
        -v | --verbose)
            shift && verbose=$1 ;;
    esac
    shift
done

if [ -e $directory/setup.py ] ; then
    cd $directory
    if [ $verbose -ne 0 ] ; then
        echo "Installing $solver Python API"
        python setup.py install ;
    else
        python setup.py -q install ;
    fi ;
else
    [ $verbose -ne 0 ] && echo "Skipping install for $solver Python API " || : ;
fi
