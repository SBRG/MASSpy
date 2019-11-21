Welcome to MASSpy's Documentation!
==================================

The **M**\ass **A**\ction **S**\toichiometric **S**\imulation **py**\thon
(**MASSpy**) package contains modules for the construction, simulation, and
analysis of kinetic models of biochemical reaction systems.

**MASSpy** is built to integrate seamlessly with **COBRApy** :cite:`ELPH13`, a widely used
modeling software package for constraint-based reconstruction and analysis of
biochemical reaction systems. **MASSpy**  can be used separately from
or in conjunction with **COBRApy**, providing a vast assortment of modeling
techniques and tools that enable different workflows. Additional information about
**COBRApy** can be found in its
`documentation <https://cobrapy.readthedocs.io/en/latest/index.html>`_ or
`github page <https://github.com/opencobra/cobrapy.>`_.

To install **MASSpy** with all features enabled::

    pip install masspy[all]

For additional information, please refer to the `detailed installation
instructions <https://github.com/SBRG/MASSpy/blob/master/INSTALL.rst>`_.


.. toctree::
   :numbered:
   :maxdepth: 1
   :caption: Contents:

   notebooks/getting_started_with_masspy
   notebooks/constructing_models
   notebooks/reading_writing_models
   notebooks/dynamic_simulation
   notebooks/plot_visualization
   notebooks/enzyme_modules
   notebooks/thermo_concentrations
   notebooks/ensemble_modeling
   notebooks/network_visualization
   notebooks/quality_assurance
   notebooks/global_configuration
   notebooks/cobra_to_mass
   notebooks/faq
   autoapi/index.rst
   zreferences.rst

Applications and  Examples
--------------------------

.. toctree::
   :maxdepth: 2

   applications.rst

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
