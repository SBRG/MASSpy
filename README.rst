MASSpy - Mass Action Stoichiometric Simulation in Python
========================================================

|PyVer| |PyPiVer| |PyPiBld| 
|RTD| |DocVer| |DocBld| 
|LIC| 

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

Installation
~~~~~~~~~~~~

.. include:: docs/installation/quickstart.rst

Documentation
~~~~~~~~~~~~~
The documentation for **MASSpy** can be found `here <https://masspy.readthedocs.io/>`_.

All documentation is generated using `Sphinx <https://www.sphinx-doc.org/>`_ and hosted by `ReadTheDocs <https://readthedocs.org/>`_. 

Cite
----
A manuscript is in preparation for publication and will be the proper reference
for citing the **MASSpy** in the future. In the meantime, feel free to cite the
`preprint at bioRxiv  <https://www.biorxiv.org/content/10.1101/2020.07.31.230334v1>`_


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

.. |PyPiVer| image:: https://img.shields.io/pypi/v/masspy?logo=PyPi&style=plastic
    :target: https://pypi.org/project/masspy/

.. |PyPiBld| image:: https://img.shields.io/github/workflow/status/sbrg/masspy/python-publish
    :target: https://github.com/SBRG/MASSpy/blob/master/.github/workflows/python-publish.yml

.. |DocVer| image:: https://img.shields.io/docker/automated/sbrg/masspy
    :target: https://hub.docker.com/r/sbrg/masspy

.. |DocBld| image:: https://img.shields.io/github/workflow/status/sbrg/masspy/publish-docker-images/master
    :target: https://github.com/SBRG/MASSpy/blob/master/.github/workflows/publish-docker-images.yml

.. |RTD| image:: https://img.shields.io/readthedocs/masspy/latest?logo=Read%20The%20Docs&style=plastic
    :target: https://masspy.readthedocs.io/en/latest/

.. |LIC| image:: https://img.shields.io/github/license/sbrg/masspy?logo=license&style=plastic
    :target: https://github.com/SBRG/MASSpy/blob/master/LICENSE
