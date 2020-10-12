Advanced Docker Usage
=====================
This page contains additional information about the MASSpy Docker image and container.

.. _recognized-image-build-context:

Recognized Image Build Context
------------------------------
The directory structure below outlines the expected build context with all optional aspects included when building a
Docker image from the `MASSpy Dockerfile <https://github.com/SBRG/MASSpy/blob/master/docker/Dockerfile>`_ ::

    MASSpy                   # Current Directory
    └── docker               # Root Directory for build context
        ├── Dockerfile
        ├── cplex 
        │   ├── cplex_studio1210.linux-x86-64.bin
        │   └── cplex.install.properties
        ├── gurobi
        │   └── gurobi.lic
        └── docker-entrypoint.sh

Build-time variables
--------------------
Certain `build-time variables <https://docs.docker.com/engine/reference/commandline/build/#set-build-time-variables---build-arg>`_ are set and passed as arguments
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

    Default is ``latest`` to use the latest stable release (master branch) of MASSpy.

An example build command using all of the build-time variables::

    docker build \
        --build-arg python_version=3.7 \
        --build-arg mass_version=latest \
        --build-arg verbose=0 \
        -t sbrg/masspy/masspy:latest ./docker
