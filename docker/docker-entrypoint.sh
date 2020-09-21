#!/bin/sh
set -e

if [ "${1:0:1}" = '-' ]; then
    set -- gurobi "$@"
fi

if [[ "$VERBOSE" = "yes" ]]; then
    set -x
fi

LICENSE=/usr/src/opt/licenses/gurobi.lic
if [ -f $LICENSE ]; then
    echo "Skipping license creation"
else
    echo "Configure license $GUROBI_LICENSE"
fi