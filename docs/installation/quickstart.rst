Quick Start Guide
=================
To quickly get started with the latest version of MASSpy, check out the information below!

With Python
-----------
The recommended method is to install **MASSpy** is to use ``pip`` to
install the software from the `Python Package Index <https://pypi.org/project/masspy/>`__.
It is recommended to do this inside a `virtual environment <https://docs.python-guide.org/dev/virtualenvs/>`__)::

    pip install masspy

With Docker
-----------
To quickly get started with the latest version of MASSpy using Docker, run the following commands in a shell::

    docker pull sbrg/masspy
    docker run --rm \
        --mount type=volume,src=licenses,dst=/home/masspy_user/opt/licenses \
        --mount type=volume,src=mass_project,dst=/home/masspy_user/mass_project \
        --publish 8888:8888 \
        -it sbrg/masspy

From within the container, either run ``python`` or ``jupyter notebook --ip=0.0.0.0 --port=8888`` depending on
the desired Python workspace.

Optimization in MASSpy
----------------------
By default, **MASSpy** comes with the `GLPK <https://www.gnu.org/software/glpk/>`__ solver. However, specific features of
**MASSpy** require a commercial optimization solver with additional solving capabilities. For more information, check out the
section on :ref:`/installation/solvers.rst`.
