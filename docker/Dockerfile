ARG python_version=3.7
ARG mass_version=0.1.1
# Image to checkout and install a specific version of MASSpy
FROM python:${python_version} AS env-setup
ARG python_version
ARG mass_version

# Define build arguments
# User information
ARG user=mass_user
ARG group=mass_group
ARG project=mass_project
ARG uid=1000
ARG gid=1000
# Verbosity of dockerfile output when possible
ARG verbose

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
    && echo \
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

# Build arg to determine how to install MASSpy (allows for a local version)
FROM env-setup AS mass-setup
# # Copy build context files into tmp
WORKDIR ${HOME}/tmp/MASSpy
COPY ./  ./

# Make directory for MASSpy installation
WORKDIR ${HOME}
RUN echo \
    && mkdir opt \
    && if [ "${MASS_VERSION}" != "local" ] ; then \
        cd opt \
        && git clone https://github.com/SBRG/MASSpy.git $( [ ${VERBOSE} -eq 0 ] && echo '--quiet' ) \
        && dirpath=tmp/MASSpy ; \
    else \
        # Copy local MASSpy files to expected install path
        cp -R tmp/MASSpy/ ${HOME}/opt/ \
        && dirpath=tmp/MASSpy/docker ; \
    fi \
    # Ensure workdir is home again before checking paths
    && cd ${HOME} \
    # Move solver files and entrypoint to highest level of tmp for later
    && if [ -d "$dirpath/cplex" ] ; then mv $dirpath/cplex/* tmp/ ; fi \
    && if [ -d "$dirpath/gurobi" ] ; then mv $dirpath//gurobi/* tmp/ ; fi \
    && if [ -e "$dirpath/docker-entrypoint.sh" ] ; then mv $dirpath/docker-entrypoint.sh tmp/ ; fi

# Checkout a specific version of MASSpy if desired and install
WORKDIR ${HOME}${MASS_INSTALL}
RUN echo \
    # If a local version of MASSpy, do nothing with branch checkouts
    && if [ "${MASS_VERSION}" != "local" ] ; then \
        # Determine version of MASSpy to checkout from from VCS
        if echo ${MASS_VERSION} | grep -Eq '^(latest|master)$' ; then \
            # Use latest (master) branch
            git checkout master $( [ ${VERBOSE} -eq 0 ] && echo '--quiet' ) ; \
        # Check whether SemVer tag given (excluding pre-release & build metabdata)
        elif echo ${MASS_VERSION} | grep -Pq '^(0|[1-9]\d*)\.(0|[1-9]\d*)\.(0|[1-9]\d*)?$' ; then \
            # Use version specified
            git checkout v${MASS_VERSION} $( [ ${VERBOSE} -eq 0 ] && echo '--quiet' ) ; \
        else \
            # Assume branch name given otherwise
            git checkout ${MASS_VERSION} $( [ ${VERBOSE} -eq 0 ] && echo '--quiet' ) ; \
        fi ; \
    fi \
    && pip install -r requirements.txt $( [ ${VERBOSE} -eq 0 ] && echo '--quiet' ) \
    && pip install . $( [ ${VERBOSE} -eq 0 ] && echo '--quiet' ) \
    && pip install notebook $( [ ${VERBOSE} -eq 0 ] && echo '--quiet' )

WORKDIR ${HOME}
RUN echo \
    # Make directory for docker scripts and move helper scripts into it.
    && mkdir opt/docker-scripts \
    && find ${HOME}${MASS_INSTALL}/docker -name "*.sh" -exec cp {} ./opt/docker-scripts/ ";" \
    && if [ -e "tmp/docker-entrypoint.sh" ] ; then mv tmp/docker-entrypoint.sh opt/docker-scripts/ ; fi \
    # Make directory for licenses volume and move license files into it
    && mkdir opt/licenses \
    # Copy license files and/or templates into license directory
    && if [ -e tmp/*.lic* ] ; then mv tmp/*.lic* opt/licenses/ ; fi 

# Install IBM CPLEX Optimization Studio if install.bin file present in build context
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

# Final build stage install Python API if they exist and adds new user with necessary permissions
FROM mass-setup AS mass-builder
WORKDIR ${HOME}

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