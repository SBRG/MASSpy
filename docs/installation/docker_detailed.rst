Advanced Docker Usage
=====================
This page contains additional information about the MASSpy Docker image and container.

.. _recognized-image-build-context:

Recognized Image Build Context
------------------------------
The directory structure below outlines the expected build context with all optional aspects included when building a
Docker image::

    MASSpy                   # Source directory
    └── docker               # Root directory for build context
        ├── Dockerfile
        ├── cplex
        │   ├── cplex_studio1210.linux-x86-64.bin
        │   └── cplex.install.properties
        ├── gurobi
        │   └── gurobi.lic
        └── docker-entrypoint.sh

The MASSpy image only requires the Dockerfile in its "context" to be built. Anything else is optional and will add specific funtionality
as outlined below:

**Dockerfile** :
    The `MASSpy Dockerfile <https://github.com/SBRG/MASSpy/blob/main/docker/Dockerfile>`__ required to build the image.

**cplex** :
    Directory used to install IBM CPLEX Optimization studio 12.10

    - **cplex_studio1210.linux-x86-64.bin**:
        The installer binary for CPLEX. The presence of this file triggers the CPLEX installation process.
    - **cplex.install.properties**:
        Installer properties for CPLEX. Acecpts license agreement and sets silent install. Ignored if no installer exists in build context.

**gurobi** :
    Directory used to install Gurobi Optimizer 9.0.3

    - **gurobi.lic**:
        Gurobi license file. The presence of this file triggers the Gurobi installation process.
    - **gurobi.lic.template**:
        `Template for Gurobi license <https://github.com/SBRG/MASSpy/blob/main/docker/gurobi/gurobi.lic.template>`__.
        Can be included to configure the token client license at a later point from within the container.

**docker-entrypoint.sh** :
    A shell script for the `container entrypoint <https://docs.docker.com/engine/reference/builder/#entrypoint>`__ to replace
    the customize the standard docker entrypoint behavior. Must be named ``docker-entrypoint.sh`` to work.

Build-time variables
--------------------
Certain `build-time variables <https://docs.docker.com/engine/reference/commandline/build/#set-build-time-variables---build-arg>`__ are set and passed as arguments
when building the image. Build-time variables are passed to ``--build-arg`` flag in the form of ``VARIABLE=VALUE``.
All build-args are optional and are not required to be defined at the time when the image is built.

The following build-time variables can be utilized by the MASSpy Dockerfile at the time of build:

**verbose**
    Integer 0 or 1 determining whether to include additional output as the image builds.
    Can be either the value ``0`` to disabled verbosity, or ``1`` to enabled it.
    Primarily for debugging purposes. Default value is ``0``.

**python_version**
    Indicates python base image to use. Must be Python 3.6+. Default is ``3.7``.

**mass_version**
    The branch or tagged version of MASSpy to use in the Docker container. Value will be passed to ``git checkout``. Must be one of the following:

    * A branch on the MASSpy GitHub Repository.
    * ``{MAJOR}.{MINOR}.{PATCH}`` to use a specific version of MASSpy.

    Default is ``latest`` to use the latest stable release (main branch) of MASSpy.

An example build command using all of the build-time variables::

    docker build \
        --build-arg python_version=3.7 \
        --build-arg mass_version=latest \
        --build-arg verbose=0 \
        -t sbrg/masspy:latest ./docker


Using a local installation of MASSpy
------------------------------------
To use the local installation of MASSpy when building the docker image, navigate to the directory containing the local installation of MASSpy
and run the following build command::

    docker build \
        --build-arg mass_version=local \
        -t sbrg/masspy:local \
        -f ./docker/Dockerfile ./

The resulting image ``sbrg/masspy:local`` can then be used to build a container using ``docker run``.
Note that will install the local version of **MASSpy** in `editable mode <https://pip.pypa.io/en/stable/cli/pip_install/#editable-installs>`__.
