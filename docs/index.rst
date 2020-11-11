MASSpy: Modeling Dynamic Biological Processes in Python
=======================================================

|PyVer| |PyPiVer| |DocVer| |DocImgSz| |RTD| |LIC| |Black|


Welcome to MASSpy's Documentation!
----------------------------------
The **M**\ass **A**\ction **S**\toichiometric **S**\imulation **py**\thon
(**MASSpy**) package contains modules for the construction, simulation, and
analysis of kinetic models of biochemical reaction systems.

**MASSpy** is built to integrate seamlessly with **COBRApy** :cite:`ELPH13`, a widely used
modeling software package for constraint-based reconstruction and analysis of
biochemical reaction systems. **MASSpy**  can be used separately from
or in conjunction with **COBRApy**, providing a vast assortment of modeling
techniques and tools that enable different workflows. Additional information about
**COBRApy** can be found in its
`documentation <https://cobrapy.readthedocs.io/en/stable/>`__ or
`GitHub page <https://opencobra.github.io/cobrapy/>`__.

Citation
~~~~~~~~
A manuscript is in preparation for publication and will be the proper reference for citing the MASSpy software package in the future.
In the meantime, feel free to cite the preprint :cite:`HZK+20`, which can be found at bioRxiv.

   The code and instsructions to reproduce the results presented in the publication is located
   in the `MASSpy-publication GitHub Repository <https://github.com/SBRG/MASSpy-publication>`__.

Installation and Setup
----------------------
There are various ways to get started with the **MASSpy** package. The guides below provide instructions on how to set up a **MASSpy** environment best suited to your needs.

**Quick Start Guide**:
   Ready to dive into **MASSpy** right away? Check out the :ref:`/installation/quickstart.rst`.

**Optimization Solvers**:
   In order to utilize certain **MASSpy** features, additional optimization capabilities (e.g., quadratic programming) are necessary. 
   Read more at :ref:`/installation/solvers.rst`.

**Docker Containers**:
   Need a standardized, ready-to-deploy container for your project? Learn how to set up :ref:`/installation/docker.rst` for **MASSpy**.  

.. toctree::
   :maxdepth: 1
   :name: install-toc
   :caption: Installation and Setup
   :hidden:

   installation/quickstart.rst
   installation/solvers.rst
   installation/docker.rst
   installation/docker_detailed.rst


Once MASSpy is installed, check out the step-by-step tutorials below to learn how to use **MASSpy**!

.. toctree::
   :numbered:
   :maxdepth: 1
   :caption: Step-by-Step Tutorials

   tutorials/getting_started_with_masspy.ipynb
   tutorials/constructing_models.ipynb
   tutorials/reading_writing_models.ipynb
   tutorials/dynamic_simulation.ipynb
   tutorials/plot_visualization.ipynb
   tutorials/enzyme_modules.ipynb
   tutorials/thermo_concentrations.ipynb
   tutorials/ensemble_modeling.ipynb
   tutorials/network_visualization.ipynb
   tutorials/quality_assurance.ipynb
   tutorials/global_configuration.ipynb
   tutorials/cobra_to_mass.ipynb
   tutorials/compartments.ipynb
   tutorials/import_export_optimization.ipynb


Example Gallery
----------------
Interested in seeing more of **MASSpy** in action? Browse through the :ref:`Gallery </gallery/index.rst>`.

.. toctree::
   :maxdepth: 1
   :caption: Gallery
   :hidden:

   gallery/index.rst


Educational Resources
---------------------
Need to review the basic principles of dynamic simulation and analysis? Educational resources utilizing **MASSpy** are outlined below and available for academic purposes. 

   * :ref:`/education/sb2/index.rst`


.. toctree::
   :maxdepth: 1
   :caption: Educational Resources
   :hidden:

   education/sb2/index.rst


API
---
Not sure how to use a specific method or function? Try searching the :ref:`/autoapi/index.rst`!

.. toctree::
   :maxdepth: 1
   :caption: API
   :hidden:

   autoapi/index.rst


.. toctree::
   :maxdepth: 1
   :caption: Additional Resources:

   additional/faqs.ipynb
   additional/code_repositories.rst
   zreferences.rst


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`


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

.. |Black| image:: https://img.shields.io/badge/code%20style-black-000000.svg
    :target: https://github.com/psf/black
