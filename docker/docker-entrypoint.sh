#!/bin/sh

set -e

# Set verbose and debug options
if [ -n '${VERBOSE}' ] ; then
    case ${VERBOSE} in
        [0-1]) set +v && set +x ;;
        # For extra verbosity/debugging purposes
        2) set -v ;;
        3) set -x ;;
    esac
fi

# Set expected license directory
license_dir=${HOME}/opt/licenses
# Let user know which optimizer licenses, and by extension,
# which optimizers are available in container.
for solver in 'CPLEX' 'Gurobi' ; do
    case $solver in
        'CPLEX')
            solver_lic_path=$license_dir/CPLEX ;;
        'Gurobi')
            solver_lic_path=$license_dir/gurobi.lic ;;
    esac
    if [ ${VERBOSE} -ne 0 ] ; then
        if [ -e $solver_lic_path ] ; then
            echo "$solver license exists at $solver_lic_path" ;
        else
            echo "No $solver license available" ;
        fi ;
    fi ;
    done

# If Gurobi avalable, set 'GRB_LICENSE_FILE' to point to license
if [ -e $license_dir/gurobi.lic ] ; then
    export GRB_LICENSE_FILE=$license_dir/gurobi.lic ;
fi

exec "$@"
