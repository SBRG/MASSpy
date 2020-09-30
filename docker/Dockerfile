ARG python_version=3.7
# Image to checkout and install a specific version of MASSpy
FROM python:${python_version} AS env_setup
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

RUN echo \
    && [ ${verbose} -ne 0 ] \
    && echo 'Building image using the following:\n' \
    && echo 'PYTHON_VERSION '${PYTHON_VERSION} \
    && echo 'MASS_VERSION: '${MASS_VERSION} \
    && echo 'USER (UID): '${USER}' ('${UID}')' \
    && echo 'GROUP (GID): '${GROUP}' ('${GID}')' \
    && echo 'PROJECT_VOLUME: '${PROJECT_VOLUME}'\n' \
    || :

# Gurobi environment variables
ENV GUROBI_MAJOR_VERSION 9.0
ENV GUROBI_VERSION 9.0.3
ENV GUROBI_INSTALL /opt/gurobi
ENV GUROBI_HOME $GUROBI_INSTALL/linux64
ENV LD_LIBRARY_PATH=$LD_LIBRARY_PATH:${HOME}$GUROBI_HOME/lib


FROM env_setup AS mass_setup
# TODO Add labels

# Clone MASSpy repository
WORKDIR ${HOME}/opt/
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

# # Image to download and install a specific version of IBM CPLEX Optimization Studio
# FROM mass_setup AS cplex_setup

# Image to download and install a specific version of Gurobi Optimizer
FROM env_setup AS gurobi_setup

WORKDIR ${HOME}${GUROBI_INSTALL}
RUN echo \
    # Get gurobi package from internet, unzip with tar and copy into the installation directory
    && wget -qO- 'https://packages.gurobi.com/'${GUROBI_MAJOR_VERSION}'/gurobi'${GUROBI_VERSION}'_linux64.tar.gz' \
    | tar xz -C . \
    && mv 'gurobi'$( echo $GUROBI_VERSION | sed -e 's/\.//g' )'/linux64/' . \
    && for to_clean in 'docs' 'examples' 'src' 'matlab' 'R' ; \
        do rm -rf 'linux64/'${to_clean} ; \
        done \
    && rm -rf 'gurobi'$( echo $GUROBI_VERSION | sed -e 's/\.//g' )

# Image to setup optimizers for MASSpy
FROM mass_setup AS mass_builder
# TODO add labels

# Install gurobipy module (Python API) for the Gurobi Optimizer
WORKDIR ${HOME}
COPY --from=gurobi_setup ${HOME} ./
WORKDIR ${HOME}${GUROBI_HOME}
RUN pwd && ls && python setup.py install

# Copy all files from local to docker 
WORKDIR ${HOME}${MASS_INSTALL}
COPY ./docker-entrypoint.sh ./docker/

# Copy token client gurobi license from local
# License created based on gurobi.lic.template file
COPY ./gurobi.lic* ./docker/

# Set home as final work directory
WORKDIR ${HOME}
RUN echo \
    # Move entrypoint script to home to avoid path issues
    && cp ${HOME}${MASS_INSTALL}/docker/docker-entrypoint.sh . \
    # Copy license files
    && mkdir /opt/licenses \
    && cp ${HOME}${MASS_INSTALL}/docker/gurobi.lic* /opt/licenses

# Add user and set user permissions for everything in home folder
RUN echo \
    && groupadd -g ${GID} ${GROUP} \
    && useradd -M -g ${GID} -G ${GROUP} -u ${UID} ${USER}
RUN echo \
    && mkdir ${HOME}${PROJECT_VOLUME} \
    && chown -R ${USER} ${HOME} \
    && chmod -R a+rxw ${HOME}

USER ${USER}

# Make volume for licenses
VOLUME licenses
# Make volume for product
VOLUME ${PROJECT_VOLUME}

ENTRYPOINT ["./docker-entrypoint.sh"]