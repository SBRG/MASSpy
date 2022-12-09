ARG python_version=3.7
ARG mass_version=0.1.7
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
    PROJECT_VOLUME=/${project}

# MASSpy environment variable, allow version to be changed at time of build
ENV MASS_VERSION=${mass_version}
ENV MASS_INSTALL /opt/MASSpy
ENV LICENSE_DIR /opt/licenses
# Verbosity of image build and entrypoint script
ENV VERBOSE=${verbose:-0}

RUN [ ${VERBOSE} -ne 0 ] \
    && echo 'Building image using the following:\n' \
    && echo 'PYTHON_VERSION '${PYTHON_VERSION} \
    && echo 'MASS_VERSION: '${MASS_VERSION} \
    && echo 'USER (UID): '${USER}' ('${UID}')' \
    && echo 'GROUP (GID): '${GROUP}' ('${GID}')' \
    && echo 'PROJECT_VOLUME: '${HOME}${PROJECT_VOLUME} \
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

FROM env-setup AS mass-setup
# # Copy build context files into tmp
WORKDIR ${HOME}/tmp/MASSpy
COPY ./  ./

# Make directory for MASSpy installation
WORKDIR ${HOME}/opt
RUN echo \
    && if [ "${MASS_VERSION}" != "local" ] ; then \
        git clone https://github.com/SBRG/MASSpy.git $( [ ${VERBOSE} -eq 0 ] && echo '--quiet' ) \
        && dirpath=../tmp/MASSpy ; \
    else \
        # Copy local MASSpy files to expected install path
        cp -R ../tmp/MASSpy/ ${HOME}/opt/ \
        && dirpath=../tmp/MASSpy/docker ; \
    fi \
    # Make directory for licenses volume and move license files into it
    && mkdir ..${LICENSE_DIR} \
    # Copy license files into license directory
    && if [ -e $dirpath/**/*.lic ] ; then mv $dirpath/**/*.lic ..${LICENSE_DIR}/ ; fi

# Checkout a specific version of MASSpy if desired and install
WORKDIR ${HOME}${MASS_INSTALL}
RUN echo \
    # If a local version of MASSpy, do nothing with branch checkouts
    && if [ "${MASS_VERSION}" != "local" ] ; then \
        # Determine version of MASSpy to checkout from from VCS
        if echo ${MASS_VERSION} | grep -Eq '^(latest|main)$' ; then \
            # Use latest (main) branch
            git checkout main $( [ ${VERBOSE} -eq 0 ] && echo '--quiet' ) ; \
        else \
            # Assume SemVer given otherwise
            git checkout v${MASS_VERSION} $( [ ${VERBOSE} -eq 0 ] && echo '--quiet' ) ; \
        fi \
        && pip install . $( [ ${VERBOSE} -eq 0 ] && echo '--quiet' ) ; \
    else \
        pip install -e . $( [ ${VERBOSE} -eq 0 ] && echo '--quiet' ) ; \
    fi \
    && pip install notebook $( [ ${VERBOSE} -eq 0 ] && echo '--quiet' )

# Install IBM CPLEX Optimization Studio if install.bin file present in build context
WORKDIR ${HOME}
RUN echo \
    && cplex_installer=cplex_studio*.linux*.bin \
    && if [ -f $dirpath/cplex/$cplex_installer ] ; then \
        if [ ${VERBOSE} -ne 0 ] ; then echo "Installing CPLEX" ; fi \
        && chmod a+rwx $dirpath/cplex/$cplex_installer \
        && $dirpath/cplex/$cplex_installer \
            -f ./cplex.installer.properties \
            -DUSER_INSTALL_DIR=${HOME}${CPLEX_INSTALL} \
        # Copy CPLEX license into license directory
        && mkdir .${LICENSE_DIR}/CPLEX \
        && mv .${CPLEX_INSTALL}/license/*.txt .${LICENSE_DIR}/CPLEX/ \
        # Clean files that aren't needed
        && for to_clean in 'concert' 'cpoptimizer' 'doc' 'license' 'opl' 'python' 'Uninstall' ; \
            do rm -rf .${CPLEX_INSTALL}/$to_clean ; \
            done \
        && for to_clean in 'examples' 'matlab' 'readmeUNIX.html' ; \
            do rm -rf .${CPLEX_INSTALL}/cplex/$to_clean ; \
            done ; \
    else \
        if [ ${VERBOSE} -ne 0 ] ; then echo "No installer found for CPLEX" ; fi \
    fi

# Install Gurobi Optimizer if gurobi.lic present
WORKDIR ${HOME}
RUN echo \
    # Only build if gurobi.lic file present
    && if [ -f .${LICENSE_DIR}/gurobi.lic ] ; then \
        gurobi_package='gurobi'${GUROBI_VERSION}'_linux64.tar.gz' \
        && if [ ! -f tmp/$gurobi_package ] ; then \
            if [ ${VERBOSE} -ne 0 ] ; then echo 'Fetching Gurobi...' ; fi \
            && wget -qP tmp/ 'https://packages.gurobi.com/'$( echo ${GUROBI_VERSION} | sed 's/\(^[0-9]*.[0-9]*\).*/\1/' )'/'$gurobi_package ; \
        fi \
        && mkdir -p .${GUROBI_HOME} \
        && tar xzf tmp/$gurobi_package -C ./tmp \
        #Only keep what is necessary
        && for to_keep in 'lib' 'include' 'bin' 'setup.py'; \
            do mv './tmp/gurobi'$( echo ${GUROBI_VERSION} | sed -e 's/\.//g' )'/linux64/'$to_keep .${GUROBI_HOME}/ ; \
            done ; \
    else \
        if [ ${VERBOSE} -ne 0 ] ; then echo "No Gurobi license, skipping Gurobi installation" ; fi ; \
    fi

# Remove tmp & remove dos2unix
RUN echo \
    && rm -rf tmp \
    # Ensure scripts are Unix compatible using dos2unix
    && apt-get update && apt-get install -y dos2unix \
    && find .${MASS_INSTALL}/docker/docker-scripts/ -name "*.sh" -exec dos2unix {} ";" \
    && apt-get --purge remove -y dos2unix \
    && rm -rf /var/lib/apt/lists/*

# Final build stage install Python API if they exist and adds new user with necessary permissions
FROM mass-setup AS mass-builder
WORKDIR ${HOME}
RUN echo \
    # Install cplex module (Python API) for the IBM CPLEX Optimization Studio if present
    && .${MASS_INSTALL}/docker/docker-scripts/pyinstall-solver.sh \
        --directory ${HOME}${CPLEX_HOME} \
        --solver "CPLEX" \
        --verbose ${VERBOSE} \
    # Install gurobipy module (Python API) for the Gurobi Optimizer if present
    && .${MASS_INSTALL}/docker/docker-scripts/pyinstall-solver.sh \
        --directory ${HOME}${GUROBI_HOME} \
        --solver "Gurobi" \
        --verbose ${VERBOSE}

# Make directories for volume in product workspace
WORKDIR ${HOME}${PROJECT_VOLUME}

# Set home as final work directory
WORKDIR ${HOME}

# Add user and set user permissions for everything in home folder
RUN echo \
    && groupadd -g ${GID} ${GROUP} \
    && useradd -M -g ${GID} -G ${GROUP} -u ${UID} ${USER} \
    && chown -R ${USER} ${HOME}/ \
    && chmod -R 777 ${HOME}/

USER ${USER}
ENTRYPOINT [ "./opt/MASSpy/docker/docker-scripts/docker-entrypoint.sh" ]
CMD [ "sh" ]
