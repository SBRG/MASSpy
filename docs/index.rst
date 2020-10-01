MASSpy: Modeling Dynamic Biological Processes in Python
=======================================================

|PyVer| |PyPi| |RTD| |LIC| |GHDL|

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
`github page <https://github.com/opencobra/cobrapy.>`__.


Installation and Setup
----------------------
There are various ways to get started with the **MASSpy** package. The guides below provide instructions on how to set up a **MASSpy** environment best suited to your needs.

**Quick Start Guide**:
   Ready to dive into **MASSpy** right away? Check out the :ref:`/installation/quickstart.rst` 

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


.. toctree::
   :numbered:
   :maxdepth: 1
   :name: tutorial-toc
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
   tutorials/exporting_optimization_problems.ipynb


.. toctree::
   :maxdepth: 1
   :name: example-gallery
   :caption: Example Gallery

   example_gallery/index.rst


.. toctree::
   :maxdepth: 1
   :name: api-toc
   :caption: API:
   :hidden:

   autoapi/index.rst


.. toctree::
   :maxdepth: 1
   :name: education-toc
   :caption: Educational Resources:

   education/sb2/index.rst


.. toctree::
   :maxdepth: 1
   :name: additional-toc
   :caption: Additional Resources:
   
   additional/faqs/faqs.rst
   additional/code_repositories.rst
   zreferences.rst


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`


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
