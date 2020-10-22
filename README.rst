MASSpy - Mass Action Stoichiometric Simulation in Python
========================================================

|PyVer| |PyPiVer| |DocVer| |DocImgSz| |RTD| |LIC| 


What is MASSpy?
--------------
The **M**\ass **A**\ction **S**\toichiometric **S**\imulation **py**\thon
(**MASSpy**) package contains modules for the construction, simulation, and
analysis of kinetic models of biochemical reaction systems.

**MASSpy** is built to integrate seemlessly with **COBRApy**, a widely used
modeling software package for constraint-based reconstruction and analysis of
biochemical reaction systems. Therefore **MASSpy**  can be used seperately from
or in conjuction with **COBRApy**, thereby providing a wide range of modeling
workflows and techniques.

Additional information about **COBRApy** can be found in its
`documentation <https://cobrapy.readthedocs.io/en/latest/index.html>`_ or
`github page <https://github.com/opencobra/cobrapy>`_.


Installation
------------
To quickly get started with the latest version of MASSpy, check out the information below taken from the
`Quick Start <https://masspy.readthedocs.io/en/latest/installation/quickstart.html>`_ guide!


With Python
~~~~~~~~~~~
The recommended method is to install **MASSpy** is to use ``pip`` to
install the software from the `Python Package Index <https://pypi.python.org/pypi/masspy>`_.
It is recommended to do this inside a `virtual environment <http://docs.python-guide.org/en/latest/dev/virtualenvs/>`_)::

    pip install masspy


With Docker
~~~~~~~~~~~
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
~~~~~~~~~~~~~~~~~~~~~~
By default, **MASSpy** comes with the `GLPK <https://www.gnu.org/software/glpk/>`__ solver. However, specific features of
**MASSpy** require a commercial optimization solver with additional solving capabilities. For more information, check out the
section on `solvers <https://masspy.readthedocs.io/en/latest/installation/quickstart.html>`_.


Documentation
-------------
The documentation for **MASSpy** can be found `here <https://masspy.readthedocs.io/>`_. All documentation is generated using `Sphinx <https://www.sphinx-doc.org/>`_ and hosted by `ReadTheDocs <https://readthedocs.org/>`_. 


Cite
----
A manuscript is in preparation for publication and will be the proper reference
for citing the **MASSpy** in the future. In the meantime, feel free to cite the
`preprint at bioRxiv <https://www.biorxiv.org/content/10.1101/2020.07.31.230334v1>`_.

    MASSpy: Building, simulating, and visualizing dynamic biological models in Python using mass action kinetics
    Zachary B. Haiman, Daniel C. Zielinski, Yuko Koike, James T. Yurkovich, Bernhard O. Palsson
    bioRxiv 2020.07.31.230334; doi: https://doi.org/10.1101/2020.07.31.230334


Contributing
------------
Contributions are always welcome! Please read the `contributions
guideline <.github/CONTRIBUTING.rst>`_
to get started.


License
-------
The **MASSpy** source is released under the MIT license. However,
please note that binary packages which include GLPK (such as the binary
wheels distributed at https://pypi.python.org/pypi/cobra) and
`libRoadRunner <https://pypi.org/project/libroadrunner/>`_ will be bound
by their licenses as well.
    

.. |PyVer| image:: https://img.shields.io/pypi/pyversions/masspy?logo=Python&style=plastic
    :target: https://www.python.org/downloads/

.. |PyPiVer| image:: https://img.shields.io/pypi/v/masspy?logo=PyPi&style=plastic
    :target: https://pypi.org/project/masspy/

.. |DocVer| image:: https://img.shields.io/docker/v/sbrg/masspy?label=Docker&logo=Docker&sort=semver&style=plastic
    :target: https://hub.docker.com/r/sbrg/masspy

.. |DocImgSz| image:: https://img.shields.io/docker/image-size/sbrg/masspy?logo=docker&sort=semver&style=plastic
    :target: https://hub.docker.com/r/sbrg/masspy

.. |RTD| image:: https://img.shields.io/readthedocs/masspy/latest?logo=Read%20The%20Docs&style=plastic
    :target: https://masspy.readthedocs.io/en/latest/

.. |LIC| image:: https://img.shields.io/github/license/sbrg/masspy?logo=license&style=plastic
    :target: https://github.com/SBRG/MASSpy/blob/master/LICENSE
