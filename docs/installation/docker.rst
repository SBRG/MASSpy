Docker
=======
**MASSpy** comes in deployable `Docker <https://docs.docker.com/>`_ container, allowing for quick access
to an isolated Python environment prepackaged the MASSpy software, all within a virtual machine that can run
in a variety of environments.

The following guide demonstrates how to setup a Docker container for **MASSpy**. It assumes that the Docker Daemon and Client have
already been installed. The guide can be broken down into three key steps:

    1. :ref:`obtaining-the-image`: An image for the MASSpy Docker container must be either obtained from an online registry or by built.
    2. :ref:`creating-the-container`: Once obtained, a container must be created from the image. 
    3. :ref:`running-the-container`: After the container is built, the final step is to run the container and get started using MASSpy!

**Important:** In order to use the *Gurobi Optimizer* or the *IBM ILOG CPLEX Optimization Studio*, the Docker image must be built locally
from a Dockerfile and a "context" containing certain files. See the secion below on :ref:`building-the-image`

About Docker
    Interested in learning more about Docker? Read more about containerization and getting started with Docker in the 
    `Docker Quick Start <https://docs.docker.com/get-started/>`_ in the official Docker documentation.

.. _obtaining-the-image:

Obtaining the image
-------------------

An image for a MASSpy Docker container can be either be downloaded from an online registry, or it
can be built from a Dockerfile and the proper build "context". 

    * The recommended method to obtain a MASSpy image is to download an image from an online registry.
    * To enable the use of a commercial optimization solver (e.g., Gurobi, CPLEX) inside the container, the
      MASSpy image must be built locally.

.. _downloading-the-image:

Downloading the image
~~~~~~~~~~~~~~~~~~~~~
Images for the MASSpy software are be found in the following registries:

`SBRG DockerHub <https://hub.docker.com/r/sbrg/masspy>`_ : 
    * **Image Name**: ``sbrg/masspy``
    * **Valid Tags**: ``latest`` | ``{MAJOR}.{MINOR}.{PATCH}`` (following `SemVer <https://semver.org/>`_ guidelines)

To pull the MASSpy image ``sbrg/masspy``, run the following in a shell::

    docker pull sbrg/masspy

By default, the ``latest`` version of MASSpy image is pulled from the registry. A tag must be included in order to download a specific image version.
For example, to pull the ``sbrg/masspy`` image with the ``latest`` tag::

    docker pull sbrg/masspy:latest

.. _building-the-image:

Building the image
~~~~~~~~~~~~~~~~~~
**Build Context**: The following directory stucture shows the minimal requirements needed as context when building the image::

    MASSpy               # Current Directory
    └── docker           # Root Directory for build context
        └── Dockerfile   # Dockerfile from VCS (https://github.com/SBRG/MASSpy)

To build the image with tag ``latest``, navigate to the ``MASSpy`` directory and use the command line::

    docker build -t sbrg/masspy:latest ./docker

Check out the :ref:`recognized-image-build-context` section for more information about the build context. 

Including ILOG CPLEX Optimization Studio
++++++++++++++++++++++++++++++++++++++++
**Licensing**: To utilize the ILOG CPLEX Optimization Studio in a Docker container, a license must be obtained first.
**Build Context**: To include CPLEX, the build context must be modified to contain the ``cplex`` subdirectory as follows::

    MASSpy
    └── docker
        ├── Dockerfile
        └── cplex 
            ├── cplex_studio1210.linux-x86-64.bin
            └── cplex.install.properties

See :ref:`cplex-solver` for more information on obtaining academic license.

Including Gurobi Optimizer
++++++++++++++++++++++++++
To utilize the Gurobi Optimizer in a Docker container, a license must be obtained first.
**Licensing**: To utilize the ILOG CPLEX Optimization Studio in a Docker container, a license must be obtained first.
**Build Context**: To include Gurobi, the build context must be modified to contain the ``gurobi`` subdirectory as follows::

    MASSpy
    └── docker
        ├── Dockerfile
        └── gurobi
            └── gurobi.lic

See :ref:`gurobi-solver` for more information on obtaining an academic license.

.. _creating-the-container:

Creating the container
----------------------
Once the MASSpy image is obtained, the next step is to run the image as a container using the following command:

    docker run \
        --mount type=volume,src=licenses,dst=/home/masspy_user/opt/licenses \
        --publish 8888:8888 \
        --name masspy_container \
        -t sbrg/masspy/masspy:latest

To break down the above command:

    * ``--mount``
        The ``--mount``flag creates a volume to allow data to persist even after a container has been stopped. 
        In this particular example, a mount of type ``volume` called ``mass_project"``is mounted to the container at
        the location ``/home/masspy_user/mass_project``. Not required for use, but highly recommended. 
    * ``--publish``
        The ``--publish`` flag publishes the container’s port  ``8888``, binding it to the host port at ``8888``.
        Required to utilize Jupyter (iPython) notebooks from inside the container.
    * ``--name``
        An optional name for the container. In this particular example, the container is given the name ``masspy_container``.
    * ``-t`` Allocate a pseudo-TTY. Required.
    
If optimization solvers are included when building the image, it is recommended to mount the ``licenses`` volume
as well. This can be done via the following::

    docker run\
        --mount type=volume,src=licenses,dst=/home/masspy_user/opt/licenses \
        --mount type=volume,src=mass_project,dst=/home/masspy_user/mass_project \
        --publish 8888:8888 \
        --name masspy_container \
        -t sbrg/masspy/masspy:latest

.. _running-the-container:

Running the container
----------------------
To start the container at the time of creation, the ``-i`` flag can be included in the run command::

    docker run\
        --mount type=volume,src=licenses,dst=/home/masspy_user/opt/licenses \
        --mount type=volume,src=mass_project,dst=/home/masspy_user/mass_project \
        --publish 8888:8888 \
        --name masspy_container \
        -it sbrg/masspy/masspy:latest

If ``masspy_container`` was created without an interactive shell, or to resume a stopped container::

    docker start -i masspy_container

.. _stopping-the-container:


Running MASSpy from the container
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Once a container has been started with an interactive shell allocated ( the ``-it`` flag ), either a Jupyter (iPython)
notebook or Python itself can be started by running one of the following from the shell within the container

    * To start python, run ``python`` 
    * To start an Jupyter notebook, run ``jupyter notebook --ip=0.0.0.0 --port=8888``. 

To stop the inteactive shell and exit the container, run the ``exit`` command.

Stopping and removing the container
-----------------------------------
To stop the container so that it can be resumed at a later point::

    docker stop masspy_container

To remove the container entirely::

    docker rm masspy_container
