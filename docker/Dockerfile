ARG MASS_VERSION=0.1.0
# Image to checkout and install a specific version of MASSpy
FROM python:3.7 AS env_setup
# MASSpy environment variable, allow version to be changed at time of build
ARG MASS_VERSION
ENV MASS_VERSION=$MASS_VERSION
ENV MASS_INSTALL /opt/MASSpy

# Gurobi environment variables
ENV GUROBI_MAJOR_VERSION 9.0
ENV GUROBI_VERSION=9.0.3
ENV GUROBI_INSTALL /opt/gurobi
ENV GUROBI_HOME $GUROBI_INSTALL/linux64
ENV LD_LIBRARY_PATH $GUROBI_HOME/lib

# TODO add labels

FROM env_setup AS mass_setup
# TODO add labels

WORKDIR /usr/src/opt/
RUN git clone https://github.com/SBRG/MASSpy.git

# Checkout a specific version of MASSpy if desired and install
WORKDIR /usr/src${MASS_INSTALL}
RUN echo \
    # Specify latest release, development branch, or a specific version of MASSpy, otherwise default to latest stable release
    && if [ "${MASS_VERSION}" != "latest" ] ; \
        then if [ "${MASS_VERSION}" = "development" ] ; \
            then branch_or_tag=${MASS_VERSION}; \
            else branch_or_tag="v${MASS_VERSION}"; \
            fi; \
        else branch_or_tag=master; \
        fi \
    && git checkout ${branch_or_tag} --quiet \ 
    && pip install .

# Image to download and install a specific version of Gurobi Optimizer
FROM env_setup AS gurobi_setup
# TODO add labels

WORKDIR /usr/src/${GUROBI_INSTALL}
RUN echo \
    # Get gurobi package from internet, unzip with tar and copy into the installation directory
    && wget -qO- 'https://packages.gurobi.com/'${GUROBI_MAJOR_VERSION}'/gurobi'${GUROBI_VERSION}'_linux64.tar.gz' \
    | tar xz -C . \
    && mv 'gurobi'$( echo $GUROBI_VERSION | sed -e 's/\.//g' )'/linux64/' . \
    && for to_clean in "docs" "examples" "src" "matlab" "R" ; \
        do rm -rf 'linux64/'${to_clean} ; \
        done \
    && rm -rf 'gurobi'$( echo $GUROBI_VERSION | sed -e 's/\.//g' )

# TODO Image to download and install a specific version of CPLEX Optimizer

# Image to setup optimizers for MASSpy
FROM mass_setup AS mass_builder
# TODO add labels

# Install Gurobi Optimizer
WORKDIR /usr/src/
COPY --from=gurobi_setup /usr/src/ ./
WORKDIR /usr/src/${GUROBI_HOME}
RUN python setup.py install

# TODO Install CPLEX Optimizer
# WORKDIR /usr/src/
# COPY --from=cplex_setip /usr/src/ ./
# WORKDIR /usr/src/${CPLEX_HOME}
# RUN python setup.py install

WORKDIR /usr/
COPY --from=mass_setup /usr/src/${MASS_INSTALL}/docker/docker-entrypoint.sh ./bin
RUN chmod 755 ./bin/docker-entrypoint.sh
ENTRYPOINT [ "./bin/docker-entrypoint.sh" ]