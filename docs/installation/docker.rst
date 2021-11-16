Using MASSpy with Docker
========================
**MASSpy** comes in deployable `Docker <https://docs.docker.com/>`__ container, allowing for quick access
to an isolated Python environment prepackaged the MASSpy software, all within a virtual machine that can run
in a variety of environments.

The following guide demonstrates how to setup a Docker container for **MASSpy**. It assumes that the Docker Daemon and Client have
already been installed. The guide can be broken down into three key steps:

    1. :ref:`obtaining-the-image`: An image for the MASSpy Docker container must be either obtained from an online registry or by built.
    2. :ref:`creating-the-container`: Once obtained, a container must be created from the image.
    3. :ref:`running-the-container`: After the container is built, the final step is to run the container and get started using MASSpy!

**Important:** In order to use the *Gurobi Optimizer* or the *IBM ILOG CPLEX Optimization Studio*, the Docker image must be built locally
from a Dockerfile and a "context" containing certain files. See the secion below on :ref:`building-the-image`.

About Docker
    Interested in learning more about Docker? Read more about containerization and getting started with Docker in the
    `Docker Quick Start <https://docs.docker.com/get-started/>`__ in the official Docker documentation.


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

`SBRG DockerHub <https://hub.docker.com/r/sbrg/masspy>`__ :
    * **Image Name**: ``sbrg/masspy``
    * **Tags**: A full list of tags can be found `here <https://hub.docker.com/r/sbrg/masspy/tags>`__

To pull the MASSpy image ``sbrg/masspy``, run the following in a shell::

    docker pull sbrg/masspy

A tag must be included in order to download a specific image version. For example, to pull the ``sbrg/masspy`` image with the ``latest`` tag::

    docker pull sbrg/masspy:latest

By default, the ``latest`` version of MASSpy image is pulled from the registry.

.. _building-the-image:

Building the image
~~~~~~~~~~~~~~~~~~
**Build Context**: The following directory stucture shows the minimal requirements needed as context when building the image::

    MASSpy               # Source directory
    └── docker           # Root directory for build context
        └── Dockerfile   # Dockerfile from VCS (https://github.com/SBRG/MASSpy)

To build the image with tag ``latest``, navigate to the ``MASSpy`` directory and use the command line::

    docker build -t sbrg/masspy:latest ./docker

**Windows Users:** Please note the following issue about `running linux containers using Docker for Windows <https://github.com/docker/for-win/issues/1340>`__.


.. _including-cplex-optimizer:

Including ILOG CPLEX Optimization Studio 12.10
++++++++++++++++++++++++++++++++++++++++++++++
To utilize the ILOG CPLEX Optimization Studio in a Docker container, a license must be obtained first.
See :ref:`cplex-solver` for more information on obtaining an academic license.

Once a CPLEX license has been obtained:

    1. Download the installer ``cplex_studioXXXX.linux-x86-64.bin`` from CPLEX, replacing "XXXX"
       for the version number without punctuation (e.g., 1210).
    2. Place the installer into the ``cplex`` directory in the build context as outlined below.
    3. Place the file ``cplex.install.properties`` into the build context to accept the license
       agreement and to enable silent install.

.. note::
    The CPLEX installer must be for **LINUX** to be compatible with the containers built using the
    `MASSpy Dockerfile <https://github.com/SBRG/MASSpy/blob/main/docker/Dockerfile>`__.

**Build Context**: To include CPLEX, the build context must be modified to contain the ``cplex`` subdirectory as follows::

    MASSpy
    └── docker
        ├── Dockerfile
        └── cplex
            ├── cplex_studio1210.linux-x86-64.bin
            └── cplex.install.properties


.. _including-gurobi-optimizer:

Including Gurobi Optimizer 9.0.3
++++++++++++++++++++++++++++++++
To utilize the Gurobi Optimizer in a Docker container, a `floating license <https://www.gurobi.com/documentation/9.0/quickstart_linux/setting_up_and_using_a_flo.html>`__
must be obtained first. See :ref:`gurobi-solver` for more information on obtaining a floating license.

Once a floating Gurobi license has been obtained:

    1. Copy the `gurobi.lic.template <https://github.com/SBRG/MASSpy/blob/main/docker/gurobi/gurobi.lic.template>`__ and
       rename the file ``gurobi.lic``.
    2. Modify the license file according to the
       `Gurobi documentation <https://www.gurobi.com/documentation/9.0/quickstart_linux/creating_a_token_server_cl.html>`__.
    3. Place the license file into the ``gurobi`` directory in the build context as outlined below.

**Build Context**: To include Gurobi, the build context must be modified to contain the ``gurobi`` subdirectory as follows::

    MASSpy
    └── docker
        ├── Dockerfile
        └── gurobi
            └── gurobi.lic

Additional information
++++++++++++++++++++++
For more information about the build context for the MASSpy image, see the :ref:`recognized-image-build-context` section.

.. _creating-the-container:

Creating the container
----------------------
Once the MASSpy image is obtained, the next step is to run the image as a container using the following command::

    docker run \
        --name mass-container \
        --mount type=volume,src=mass_project,dst=/home/masspy_user/mass_project \
        --publish 8888:8888 \
        -it sbrg/masspy:latest

To break down the above command:

    * --name :
        The ``--name`` flag sets an optional name for the container that can be used to reference the container
        with the Docker Client. Here, the container is named ``mass-container``.
    * --mount :
        The ``--mount`` flag creates a volume to allow data to persist even after a container has been stopped.
        In this particular example, a mount of type ``volume`` called ``mass_project`` is mounted to the container at
        the location ``/home/masspy_user/mass_project``. Not required for use, but highly recommended.
    * --publish :
        The ``--publish`` flag publishes the container’s port  ``8888``, binding it to the host port at ``8888``.
        Required to utilize Jupyter (iPython) notebooks from inside the container.
    * -it :
        Allocate a pseudo-TTY and create an interactive shell in the container.

If optimization solvers are included when building the image, it is recommended to mount the ``licenses`` volume
as well. This can be done via the following::

    docker run \
        --name mass-container \
        --mount type=volume,src=licenses,dst=/home/masspy_user/opt/licenses \
        --mount type=volume,src=mass_project,dst=/home/masspy_user/mass_project \
        --publish 8888:8888 \
        -it sbrg/masspy:latest

.. note::
    Containers names must be unique. To re-use a name for a new container, the previous container must first be removed.


.. _running-the-container:

Running MASSpy with the container
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Once a container has been started with an interactive shell allocated ( the ``-it`` flag ), either a Jupyter (iPython)
notebook or Python itself can be started by running one of the following from the shell within the container

    * To start python, run ``python``
    * To start a Jupyter notebook, run ``jupyter notebook --ip=0.0.0.0 --port=8888``.

To stop the inteactive shell and exit the container, run the ``exit`` command.


Resuming the container
~~~~~~~~~~~~~~~~~~~~~~
To resume the container ``mass-container`` after it has been stopped::

    docker start -i mass-container


Cleanup
~~~~~~~
To remove the container ``mass-container`` entirely::

    docker rm mass-container

To remove the image ``sbrg/masspy:latest`` entirely::

    docker rmi sbrg/masspy:latest


Troubleshooting
---------------
Need help trouble shooting Docker for your system? Try searching the official Docker resources:

    `Docker CE for Linux <https://github.com/docker/for-linux/>`__ |
    `Docker Desktop for Mac <https://github.com/docker/for-mac/>`__ |
    `Docker Desktop for Windows <https://github.com/docker/for-win/>`__
