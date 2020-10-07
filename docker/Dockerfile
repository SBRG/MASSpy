ARG python_version=3.7
# Image to checkout and install a specific version of MASSpy
FROM python:${python_version} AS env_setup
ARG python_version
# Define build arguments
# User information
ARG user=mass_user
ARG group=mass_group
ARG project=mass_project
ARG uid=1000
ARG gid=1000
# Verbosity of dockerfile output when possible
ARG verbose

# Version of MASSpy to checkout
ARG mass_version=0.1.0

# User environment variables
ENV USER=${user} \
    GROUP=${group} \
    UID=${uid} \
    GID=${gid} \
    HOME=/home/${user} \
    PROJECT_VOLUME=${project}

# MASSpy environment variable, allow version to be changed at time of build
ENV MASS_VERSION=${mass_version}
ENV MASS_INSTALL /opt/MASSpy
# Verbosity of image build and entrypoint script
ENV VERBOSE=${verbose:-0}

RUN [ ${VERBOSE} -ne 0 ] \
    && echo 'Building image using the following:\n' \
    && echo 'PYTHON_VERSION '${PYTHON_VERSION} \
    && echo 'MASS_VERSION: '${MASS_VERSION} \
    && echo 'USER (UID): '${USER}' ('${UID}')' \
    && echo 'GROUP (GID): '${GROUP}' ('${GID}')' \
    && echo 'PROJECT_VOLUME: '${PROJECT_VOLUME} \
    && echo 'BIND_MOUNT: '${BIND_MOUNT}'\n' \
    || :

# CPLEX environment variables
ENV CPLEX_VERSION 12.10.0
ENV CPLEX_INSTALL /opt/CPLEX
ENV CPLEX_HOME ${CPLEX_INSTALL}/cplex/python/$python_version/x86-64_linux

# Gurobi environment variables
ENV GUROBI_VERSION 9.0.3
ENV GUROBI_INSTALL /opt/gurobi
ENV GUROBI_HOME ${GUROBI_INSTALL}/linux64
ENV LD_LIBRARY_PATH=$LD_LIBRARY_PATH:${HOME}${GUROBI_HOME}/lib

FROM env_setup AS mass_setup
# Clone MASSpy repository
WORKDIR ${HOME}/opt
RUN git clone https://github.com/SBRG/MASSpy.git $( [ ${VERBOSE} -eq 0 ] && echo '--quiet' )

# Checkout a specific version of MASSpy if desired and install
WORKDIR ${HOME}${MASS_INSTALL}
RUN echo \
    # Determine version of MASSpy to checkout from from VCS
    && if echo ${MASS_VERSION} | grep -Eq '^(latest|master)$' ; then \
        # Use latest (master) branch
        git_branch_or_tag=master ; \
    # Check whether SemVer tag given (excluding pre-release & build metabdata)
    elif echo ${MASS_VERSION} | grep -Pq '^(0|[1-9]\d*)\.(0|[1-9]\d*)\.(0|[1-9]\d*)?$' ; then \
        # Use version specified
        git_branch_or_tag=v${MASS_VERSION} ; \
    else \
        git_branch_or_tag=${MASS_VERSION} ; \
    fi \
    && git checkout ${git_branch_or_tag} $( [ ${VERBOSE} -eq 0 ] && echo '--quiet' ) \
    && pip install . $( [ ${VERBOSE} -eq 0 ] && echo '--quiet' ) \
    && pip install notebook $( [ ${VERBOSE} -eq 0 ] && echo '--quiet' )

# Copy necessary files into container from build context
WORKDIR ${HOME}
COPY *.sh cplex* gurobi* tmp/

RUN echo \
    # Make directory for docker scripts and move helper scripts into it.
    && mkdir opt/docker-scripts \
    && mv tmp/*.sh opt/docker-scripts \
    ## && cp ${HOME}${MASS_INSTALL}/docker/*.sh ./docker-scripts \
    # Make directory for licenses volume and move license files into it
    && mkdir opt/licenses \
    # Copy license files and/or templates into license directory
    && mv tmp/*.lic* opt/licenses

# Install IBM CPLEX Optimization Studio if install.binfile present in build context
WORKDIR ${HOME}
RUN echo \
    && cplex_installer=cplex_studio*.linux*.bin \
    && if [ -f tmp/$cplex_installer ] ; then \
        if [ ${VERBOSE} -ne 0 ] ; then echo "Installing CPLEX" ; fi \
        && chmod a+rwx tmp/$cplex_installer \
        && tmp/$cplex_installer \
            -f ./cplex.installer.properties \
            -DUSER_INSTALL_DIR=${HOME}${CPLEX_INSTALL} \
        # Copy CPLEX license into license directory
        && mkdir ./opt/licenses/CPLEX \
        && mv ./${CPLEX_INSTALL}/license/*.txt ./opt/licenses/CPLEX/ \
        # Clean files that aren't needed
        && for to_clean in 'concert' 'cpoptimizer' 'doc' 'license' 'opl' 'python' 'Uninstall' ; \
            do rm -rf ./${CPLEX_INSTALL}/$to_clean ; \
            done \
        && for to_clean in 'examples' 'matlab' 'readmeUNIX.html' ; \
            do rm -rf ./${CPLEX_INSTALL}/cplex/$to_clean ; \
            done ; \
    else \
        if [ ${VERBOSE} -ne 0 ] ; then echo "No installer found for CPLEX" ; fi \
    fi

# Install Gurobi Optimizer if gurobi.lic present
WORKDIR ${HOME}
RUN echo \
    # Only build if gurobi.lic file present
    && if [ -f opt/licenses/gurobi.lic ] ; then \
        gurobi_package='gurobi'${GUROBI_VERSION}'_linux64.tar.gz' \
        && if [ ! -f tmp/$gurobi_package ] ; then \
            if [ ${VERBOSE} -ne 0 ] ; then echo 'Fetching Gurobi...' ; fi \
            && wget -qP tmp/ 'https://packages.gurobi.com/'$( echo ${GUROBI_VERSION} | sed 's/\(^[0-9]*.[0-9]*\).*/\1/' )'/'$gurobi_package ; \
        fi \
        && mkdir -p ./${GUROBI_HOME} \
        && tar xzf tmp/$gurobi_package -C ./tmp \
        #Only keep what is necessary
        && for to_keep in 'lib' 'include' 'bin' 'setup.py'; \
            do mv './tmp/gurobi'$( echo ${GUROBI_VERSION} | sed -e 's/\.//g' )'/linux64/'$to_keep ./${GUROBI_HOME}/ ; \
            done ; \
    else \
        if [ ${VERBOSE} -ne 0 ] ; then echo "No Gurobi license, skipping Gurobi installation" ; fi ; \
    fi

# Remove tmp
RUN rm -rf tmp

# Install cplex module (Python API) for the IBM CPLEX Optimization Studio
RUN echo \
    && ${HOME}/opt/docker-scripts/pyinstall-solver.sh \
        --directory ${HOME}${CPLEX_HOME} \
        --solver "CPLEX" \
        --verbose ${VERBOSE} \
    && ${HOME}/opt/docker-scripts/pyinstall-solver.sh \
        --directory ${HOME}${GUROBI_HOME} \
        --solver "Gurobi" \
        --verbose ${VERBOSE}

# Make directories for volume in product workspace
WORKDIR ${HOME}/${PROJECT_VOLUME}

# Set home as final work directory
WORKDIR ${HOME}
# Add user and set user permissions for everything in home folder
RUN echo \
    && groupadd -g ${GID} ${GROUP} \
    && useradd -M -g ${GID} -G ${GROUP} -u ${UID} ${USER}
RUN echo \
    && chown -R ${USER} ${HOME}/ \
    && chmod -R 777 ${HOME}/

USER ${USER}

WORKDIR ${HOME}

ENTRYPOINT [ "./opt/docker-scripts/docker-entrypoint.sh" ]
CMD [ "sh" ]