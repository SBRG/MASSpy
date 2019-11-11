MASSpy - Mass Action Stoichiometric Simulation in Python
========================================================

|PyPI|

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
for citing the **MASSpy** in the future.

Installation
~~~~~~~~~~~~

The recommended method is to install **MASSpy** is to use ``pip`` to
`install masspy from PyPI <https://pypi.python.org/pypi/masspy>`_. It is
recommended to do this inside a `virtual environment
<http://docs.python-guide.org/en/latest/dev/virtualenvs/>`_)::

	pip install masspy

You can install all packages, including optional dependencies, directly by::

    pip install masspy[all]

If you downloaded the source code, run::

	pip install -e .

in the ``masspy`` source directory. For additional information, please refer to the
the `detailed installation instructions <INSTALL.rst>`_.

Contributing
~~~~~~~~~~~~

Contributions are always welcome! Please read the `contributions
guideline <.github/CONTRIBUTING.rst>`_
to get started.

License
-------

The **MASSpy** source is released under both the GPL and LGPL licenses. You
may choose which license you choose to use the software under. However,
please note that binary packages which include GLPK (such as the binary
wheels distributed at https://pypi.python.org/pypi/cobra) and
`libRoadRunner <https://pypi.org/project/libroadrunner/>`_ will be bound
by their licenses as well.

This program is free software: you can redistribute it and/or modify it
under the terms of the GNU General Public License or the Lesser GNU
General Public License as published by the Free Software Foundation,
either version 2 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General
Public License for more details.

.. |PyPI| image:: https://badge.fury.io/py/masspy.svg
    :target: https://pypi.python.org/pypi/masspy