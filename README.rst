MASSpy - Mass Action Stoichiometric Simulation in Python
========================================================

|PyVer| |PyPi| |RTD| |LIC| |GHDL|

What is MASSpy?
~~~~~~~~~~~~~~~
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

Cite
----
A manuscript is in preparation for publication and will be the proper reference
for citing the **MASSpy** in the future. In the meantime, feel free to cite the
`preprint at bioRxiv  <https://www.biorxiv.org/content/10.1101/2020.07.31.230334v1>`_

Installation
~~~~~~~~~~~~

The recommended method is to install **MASSpy** is to use ``pip`` to
`install MASSpy from PyPI <https://pypi.python.org/pypi/masspy>`_. It is
recommended to do this inside a `virtual environment
<http://docs.python-guide.org/en/latest/dev/virtualenvs/>`_)::

	pip install masspy

You can install all packages, including optional dependencies, directly by::

    pip install masspy[all]

If you downloaded the source code, run::

	pip install -e .

in the ``MASSpy`` source directory. For additional information, please refer to the
the `detailed installation instructions <INSTALL.rst>`_.

Contributing
~~~~~~~~~~~~

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

.. |PyPi| image:: https://img.shields.io/pypi/v/masspy?logo=PyPi&style=plastic
    :target: https://pypi.org/project/masspy/

.. |RTD| image:: https://img.shields.io/readthedocs/masspy/latest?logo=Read%20The%20Docs&style=plastic
    :target: https://masspy.readthedocs.io/en/latest/

.. |LIC| image:: https://img.shields.io/github/license/sbrg/masspy?logo=license&style=plastic
    :target: https://github.com/SBRG/MASSpy/blob/master/LICENSE

.. |GHDL| image:: https://img.shields.io/github/downloads/sbrg/masspy/total?logo=GitHub&style=social
    :target: https://github.com/SBRG/MASSpy
