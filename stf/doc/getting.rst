Getting *VELOCIraptor*
###################

*VELOCIraptor* is currently hosted in `GitHub <https://github.com/pelahi/VELOCIraptor-STF>`_.
To get a copy you can clone the repository::

 git clone https://github.com/pelahi/VELOCIraptor-STF


Compiling
=========

*VELOCIraptor* depends on:

* `GSL <https://www.gnu.org/software/gsl/>`_
* libNBody - a scientific library included with VELOCIraptor. VELOCIraptor needs this library for a number of structures, classes, and methods it provides.

Optional requirements are:

For parallel use may need the following non-standard libraries for compilation
depending on the compilation flags used:

* MPI - the Message Passing Interface (version 1.0 or higher). Many
  vendor supplied versions exist, in addition to excellent open source
  implementations, e.g.  MPICH
  (http://www-unix.mcs.anl.gov/mpi/mpich/) or LAM
  (http://www.lam-mpi.org/).

* `OpenMP <http://www.openmp.org/>`_ - API, generally included with many compilers

VELOCIraptor also can output in a variety of formats ASCII, binary, HDF and ADIOS.
HDF and ADIOS can be enabled and disabled, and require libaries.

* `Hiearchical Data Format (HDF) <https://www.hdfgroup.org/>`_ - self describing data format.
* `Adaptable IO System (ADIOS) <https://www.olcf.ornl.gov/center-projects/adios/>`_ - self describing data format.


*VELOCIraptor*'s compilation system is based on `GNU make <https://www.gnu.org/software/make/>`_.

To compile *VELOCIraptor* (assuming you are inside the ``VELOCIraptor-STF/stf`` directory already) simply::
 $> cp Makefile.config.template Makefile.config
Edit Makefile.config with you favourite editor and type::
 $> make
A list of compile time options is listed in .. doxygenfile:: ui.cxx
