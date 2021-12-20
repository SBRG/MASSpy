Optimization Solvers
====================
**MASSpy** utilizes `Optlang <https://github.com/opencobra/optlang>`__ as a common interface for different optimization solvers.
By default, **MASSpy** will come with `swiglpk <https://github.com/biosustain/swiglpk>`__, an interface to the open source (mixed integer)
linear programming (LP) solver `GLPK <https://www.gnu.org/software/glpk/>`__. However, in order to utilize specific **MASSpy** features,
a (mixed integer) quadratic programming (QP) solver are necessary.

The following features require QP support:

* Concentration solution space sampling

The following optional solvers are currently supported by **Optlang**:


.. _cplex-solver:

IBM ILOG CPLEX Optimization Studio
----------------------------------
The *IBM ILOG CPLEX Optimization Studio* (CPLEX) can be utilized through the
`CPLEX Python API <https://www.ibm.com/docs/en/icos/20.1.0?topic=cplex-setting-up-python-api>`__.

* To use CPLEX, a license must be obtained. Free academic licences are available.
* To use CPLEX with Docker, an installer file must also be downloaded.

  `Homepage <https://www.ibm.com/products/ilog-cplex-optimization-studio/>`__ |
  `Documentation <https://www.ibm.com/docs/en/icos/20.1.0?topic=tutorials-python-tutorial>`__ |
  `Academic License <https://www.ibm.com/academic/faqs/agreement/>`__

.. _gurobi-solver:

Gurobi Optimizer
----------------
The *Gurobi Optimizer* (Gurobi) is utilized through the `Gurobi Python Interface <https://www.gurobi.com/documentation/9.0/quickstart_linux/the_grb_python_interface_f.html>`__.

* To use Gurobi, a license must be obtained. Free academic licences are available.
* To use Gurobi with Docker, a floating license is required.

  `Homepage <https://www.gurobi.com/products/gurobi-optimizer/>`__ |
  `Documentation <https://www.gurobi.com/documentation/>`__ |
  `Academic License <https://www.gurobi.com/academia/academic-program-and-licenses/>`__ |
  `Floating License <https://www.gurobi.com/documentation/9.0/quickstart_linux/retrieving_a_floating_lice.html>`__

Working with other solvers
--------------------------
Prefer to work with a different optimization solver? It is possible to import/export optimization problems for use in other solvers.
Read more at :ref:`/tutorials/import_export_optimization.ipynb`.
