.. _getting:

Getting |vr|
########################

**VELOCIraptor** is currently hosted in `GitHub <https://github.com/pelahi/VELOCIraptor-STF>`_.
To get a copy you can clone the repository::

  git clone https://github.com/pelahi/VELOCIraptor-STF

**VELOCIraptor**'s compilation system is based on `cmake <https://www.cmake.org/>`_. ``cmake`` will
check that you have a proper compiler (anything supporting C++11 or later should do),
and scan the system for all required dependencies.

To compile **VELOCIraptor** run (assuming you are inside the ``VELOCIraptor-STF/`` directory already)::

 $> mkdir build
 $> cd build
 $> cmake ..
 $> make

With ``cmake`` you can also specify additional compilation flags.
For example, if you want to generate the fastest possible code
you can try this::

 $> cmake .. -DCMAKE_CXX_FLAGS="-O3 -march=native"

You can also specify a different installation directory like this::

 $> cmake .. -DCMAKE_INSTALL_PREFIX=~/my/installation/directory

Other ``cmake`` options that can be given in the command-line include:

A list of compile time options is found below in :ref:`compileoptions`.

Requirments
===========

**VELOCIraptor** depends on:

* `GSL <https://www.gnu.org/software/gsl/>`_ - the GNU Scientific Library
* **NBodylib** - a internal scientific library included with **VELOCIraptor**. **VELOCIraptor** needs this library for a number of structures, classes, and methods it provides.

Optional requirements
---------------------

For parallel use may need the following libraries are required for compilation
depending on the compilation flags used:

* MPI - the Message Passing Interface (version 1.0 or higher). Many
  vendor supplied versions exist, in addition to excellent open source
  implementations, e.g. `Open MPI <https://www.open-mpi.org/>`_, `MPICH <http://www-unix.mcs.anl.gov/mpi/mpich/>`_ or
  `LAM <http://www.lam-mpi.org/>`_.

* `OpenMP <http://www.openmp.org/>`_ - API, generally included with many compilers

**VELOCIraptor** also can output in a variety of formats: ASCII, binary, HDF and ADIOS.
HDF and ADIOS can be enabled and disabled, and require libraries.

* `Hiearchical Data Format (HDF) <https://www.hdfgroup.org/>`_ - self describing data format.
* `Adaptable IO System (ADIOS) <https://www.olcf.ornl.gov/center-projects/adios/>`_ - self describing data format.

.. _compileoptions:

Compilation Options
===================

These can be passed to ``cmake``

.. topic:: External library flags

    * Parallel APIs can be enabled by setting
        * For MPI
            | ``VR_MPI``: boolean to compile with MPI support
            | ``VR_MPI_REDUCE``: boolean that reduces impact of MPI memory overhead at the cost of extra cpu cycles. Suggested this be turned on
            | ``MPI_LIBRARY``: specify library path to MPI
            | ``MPI_EXTRA_LIBRARY``: Extra MPI libraries to link against
            | ``VR_LARGE_MPI_DOMAIN`` : Enable if mpi domain is going to contain more than max 16 bit integer number of mpi processes
        * For OpenMP
            | ``NBODY_OPENMP``: boolean to compile with OpenMP support
            | ``OpenMP_CXX_FLAGS``: string, compiler flag that enables OpenMP


    * Enable input/output formats
        * For HDF
            | ``VR_HDF5``: boolean on whether to include HDF support
            | ``VR_ALLOWPARALLELHDF5``: boolean on whether to allow for parallel HDF support (if available)
            | ``VR_ALLOWPARALLELHDF5COMPRESSIONHDF5``: boolean on whether to allow for compression parallel HDF support (THIS IS UNSTABLE, USE WITH CAUTION)
            | ``HDF5_ROOT``: specify a local directory containing HDF library.
        * for XDR (nchilada) input
            | ``VR_XDR``: boolean on whether to include XDR support
            | ``XDR_DIR``: specify a local directory containing XDR library.
        * for adios output (alpha, not yet available)
            | ``VR_ADIOS``: boolean on whether to include ADIOS support
            | ``ADIOS_DIR``: specify a local directory containing ADIOS library.

    * To set directories of required libraries
        * Set the directories of the following libraries
            | ``GSL_DIR =``

.. topic:: Internal precision and data structure flags

    * Adjust precision in stored variables and calculations
        * Calculations/properties at 32 bit float precision
            ``VR_SINGLE_PRECISION=ON``
        * all integers are 64 bit integers. Enable this if dealing with more than MAXINT total number of particles
            ``VR_LONG_INT=ON``

    * Adjust **NBodylib** Particle class data precision and memory footprint
        * Do not store the mass as all particles are the same mass. :strong:`WARNING`: :emphasis:`This is not fully implement for all types of input and requires further testing, use with caution.`
            ``VR_NO_MASS=ON``
        * Use single precision to store positions,velocities, and possibly other internal properties
            ``NBODY_SINGLE_PARTICLE_PRECISION=ON``
        * Use unsigned ints (size set by whether using long int or not) to store permanent 'particle' ids
            ``NBODY_UNSIGNED_PARTICLE_PIDS=ON``
        * Use unsigned ints (size set by whether using long int or not) to store ids (index value). Note that velociraptor uses negative index values for sorting purposes so ONLY ENABLE if library to be used with other codes.
            ``NBODY_UNSIGNED_PARTICLE_IDS=ON``

    * Hydro simulations: activate extra data structures in the **NBodylib** Particle class
        * activate gas, store self-energy
            ``VR_USE_GAS=ON``
        * activate stars only, store metallicity, formation time, star foramtion rate (for gas particles)
            ``VR_USE_STARS=ON``
        * Calculate bulk black hole properties
            ``VR_USE_BH=ON``
        * stars and gas and black holes
            ``VR_USE_HYDRO=ON``

    * Adjust memory/max size of Binary KD Tree options, used to run search particles. If tree is going to be built on more than max 32 bit integer number particles then enable, memory footprint increases
            ``VR_USE_LARGE_KDTREE=ON``

.. topic:: Operation flags

    * only calculate local density distribution for particles residing in field objects (but using all particles to estimate quantity). Default.
        ``VR_STRUCTURE_DEN=ON``
    * or just use particles inside field objects, reducing cpu cycles but will bias estimates for particle in outer region of field structures, overrides STRUCTUREDEN
        ``VR_HALO_DEN=ON``
    * flag useful for zoom simulations with a high resolution region
        ``VR_ZOOM_SIM=ON``

.. topic:: Executable flags

    * Produce SWIFTSIM compatible library (executable still produced but does simply returns warning)
        | ``VR_USE_SWIFT_INTERFACE=ON``
        | ``CMAKE_CXX_FLAGS=-fPIC``
    * Enable debugging
        ``DEBUG=ON``
