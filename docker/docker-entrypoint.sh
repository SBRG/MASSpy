#!/bin/sh

set -e

# Set verbose and debug options
if [ -n '${VERBOSE}' ] ; then 
    case ${VERBOSE} in 
        0) set +v && set +x ;;
        1) set -v ;;
        2) set -x ;;
    esac
fi

# Move client token license 

# # Configure gurobi license
# if [ ! -z ${GRB_LICENSE_FILE} ] && [ -f ${GRB_LICENSE_FILE} ]; then 
#     echo 'Gurobi license exists at '${GRB_LICENSE_FILE}
# elif [ ! -z ${GRB_LICENSE_FILE} ] ; then
#     echo 'Creating token server client license'
# else
#     echo 'No Gurobi license available.'
# fi

# # Ensures license path variable points to license volume if not set
# if [ -z ${GRB_LICENSE_FILE} ] ; then \
#     export GRB_LICENSE_FILE=${HOME}/opt/licenses/gurobi.lic
# fi

exec "$@"
