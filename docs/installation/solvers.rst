Optimization Solvers
====================
**MASSpy** utilizes `Optlang <https://github.com/opencobra/optlang>`__ as a common interface for different optimization solvers. 
By default, **MASSpy** will come with `swiglpk <https://github.com/biosustain/swiglpk>`__, an interface to the open source (mixed integer)
linear programming (LP) solver `GLPK <https://www.gnu.org/software/glpk/>`__. However, in order to utilize specific **MASSpy** features, 
a (mixed integer) quadratic programming (QP) solver are necessary. 

The following features require QP support:

* Concentration solution space sampling

The following optional solvers are currently supported by **Optlang**:

* **Gurobi Optimizer**
  `Homepage <https://www.gurobi.com/products/gurobi-optimizer/>`__ |
  `Academic Licencing <https://www.gurobi.com/academia/academic-program-and-licenses/>`__ |
  `Documentation <https://www.gurobi.com/documentation/>`__ 

* **IBM ILOG CPLEX Optimization Studio**
  `Homepage <https://www.ibm.com/products/ilog-cplex-optimization-studio/>`__ |
  `Academic Licencing <https://www.ibm.com/academic/faqs/agreement/>`__ |
  `Documentation <https://www.ibm.com/support/knowledgecenter/SSSA5P_12.10.0/ilog.odms.cplex.help/CPLEX/UsrMan/topics/APIs/Python/01_title_synopsis.html>`__ 

Working with other solvers
~~~~~~~~~~~~~~~~~~~~~~~~~~
Prefer to work with a different optimization solver? It is possible to export the optimization problem for use in other solvers.
Read more at :ref:`/tutorials/exporting_optimization_problems.ipynb`.
