.. VELOCIraptor documentation master file, created by
   sphinx-quickstart on Mon Jul 31 10:13:40 2017.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.
   **VELOCIraptor** is a C++ halo finder using MPI and OpenMP APIs.
   It comes with several associated analysis tools in python, example configuration files
   and analysis python scripts (and sample jupyter notebooks). There is an associated
   halo merger tree code **TreeFrog** (also C++ MPI+OpenMP)

Welcome to **VELOCIraptor**'s documentation!
#############################################

**VELOCIraptor** is a C++ halo finder using MPI and OpenMP APIs.
The repository also contains several associated analysis tools in python,
example configuration files and analysis python scripts (and sample jupyter notebooks).
The code can also be compiled as a library for on-the-fly halo finding within an
N-body/hydrodynamnical code. Currently integration is limited to swift but extensions are
in the works. There is an associated halo merger tree code **TreeFrog** (also C++ MPI+OpenMP)


.. toctree::
   :maxdepth: 2
   :numbered:
   :titlesonly:
   :glob:
   :hidden:
   :caption: Contents:
   getting.rst
   usage.rst
   output.rst



Sections
========
* :ref:`getting` : How to compile the code
* :ref:`usage` : How to run the code
* :ref:`output` : How to use output data
* :ref:`dev` : For developers.
* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
