Installation of masspy
=======================

All releases require Python 3.4+ to be installed before proceeding.
Mac OS X (10.7+) and Ubuntu ship with Python. Windows users without python can
download and install python from the `python website
<https://www.python.org/ftp/python/3.6.4/python-3.6.4-amd64-webinstall.exe>`_.
Please note that though Anaconda and other python distributions may work
with masspy, they are not explicitly supported (yet!).

Stable version installation
~~~~~~~~~~~~~~~~~~~~~~~~~~~

masspy can be installed with any recent installation of pip.
Instructions for several operating systems are below:

Mac OS X or Linux
-----------------

1. `install pip <http://pip.readthedocs.org/en/latest/installing.html>`_.
2. In a terminal, run ``sudo pip install mass``

We highly recommend updating ``pip`` beforehand
(``pip install pip --upgrade``).

Microsoft Windows
-----------------

The preferred installation method on Windows is also to use pip. The
latest Windows installers for Python 3.4 and up include pip, so if you
use those you will already have pip.

1. In a terminal, run ``C:\Python36\Scripts\pip.exe install masspy``
(you may need to adjust the path accordingly).

To install without pip, you will need to download and use the
appropriate installer for your version of python from the `python
package index <https://pypi.python.org/pypi/masspy/>`_.

Installation for development
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Get the `detailed contribution instructions <CONTRIBUTING.rst>`_ for
contributing to masspy.

Installation of optional dependencies
=====================================

Optional dependencies
~~~~~~~~~~~~~~~~~~~~~

On windows, these can downloaded from `this site
<http://www.lfd.uci.edu/~gohlke/pythonlibs/>`_. On Mac/Linux, they can be
installed using pip, or from the OS package manager (e.g brew, apt, yum).

1. `libsbml <http://sbml.org>`_ >= 5.10 to read/write SBML level 2
   files

   -  `Windows libsbml installer <http://www.lfd.uci.edu/~gohlke/pythonlibs/#libsbml>`_
   -  Use ``sudo pip install python-libsbml`` on Mac/Linux

2. `lxml <http://lxml.de/>`_ to speed up read/write of SBML level 3 files.
You can install all packages directly by::

    pip install "masspy[all]"

Testing your installation
=========================

Currently masspy is still under development and testing. We apologize
for the inconvenience and provide further details on this page in a
future update.

For now, if you have any questions or concerns, please email
zhaiman@eng.ucsd.edu with "masspy concerns" in the subject line and we will
get back to you as soon as we can.
