.. _getting:

Getting **VELOCIraptor**
########################

**VELOCIraptor** is currently hosted in `GitHub <https://github.com/pelahi/VELOCIraptor-STF>`_.
To get a copy you can clone the repository
::
    git clone https://github.com/pelahi/VELOCIraptor-STF

**VELOCIraptor**'s compilation system is based on `GNU make <https://www.gnu.org/software/make/>`_.
To compile **VELOCIraptor** (assuming you are inside the ``VELOCIraptor-STF/stf`` directory already) first
::
    $> cp Makefile.config.template Makefile.config

Edit Makefile.config which sets the compilation flags with you favourite editor and type:
::
    $> make

A list of compile time options is found below in `Compilation Options`_.

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

Compilation Options
===================

These can be changed by editing Makefile.config (based on the :download:`template <../Makefile.config.template>`) in your favourite editor.

.. topic:: External library flags

    * Parallel APIs can be enabled by setting
        * For MPI
            | ``MPI="on"``
            | ``MPIREDUCE="on"``
        * For OpenMP
            ``OMP="on"``

    * Enable input/output formats
        * For HDF |
            ``HDFENABLE="on"``
        * for XDF (nchilada) input
            ``XDRENABLE="on"``
        * for adios output
            ``ADIOSENABLE="on"``

    * To set local libraries
        * Activate the following flags
            | ``LOCALGSL="on"``
            | ``LOCALFFTW="on"``
            | ``LOCALHDF="on"``
            | ``LOCALADIOS="on"``
            | ``LOCALXDR="on"``

        * Set the directories of the following libraries
            | ``GSL_DIR =``
            | ``FFTW_DIR =``
            | ``HDF_DIR =``
            | ``ADIOS_DIR =``

.. topic:: Internal precision and data structure flags

    * Adjust precision in stored variables and calculations
        * Calculations/properties at 32 bit float precision
            ``SINGLEPRECISION="on"``
        * all integers are 64 bit integers. Enable this if dealing with more than MAXINT total number of particles
            ``LONGINT="on"``

    * Adjust **NBodylib** Particle class data precision and memory footprint
        * Do not store the mass as all particles are the same mass. :strong:`WARNING`: :emphasis:`This is not fully implement for all types of input and requires further testing, use with caution.`
            ``NOMASS="on"``
        * Use single precision to store positions,velocities, and possibly other internal properties
            ``SINGLEPARTICLEPRECISION="on"``
        * Use unsigned ints (size set by whether using long int or not) to store permanent 'particle' ids
            ``UNSIGNEDPARTICLEPIDS="on"``
        * Use unsigned ints (size set by whether using long int or not) to store ids (index value). Note that velociraptor uses negative index values for sorting purposes so ONLY ENABLE if library to be used with other codes.
            ``UNSIGNEDPARTICLEIDS="on"``

    * Hydro simulations: activate extra data structures in the **NBodylib** Particle class
        * activate gas, store self-energy
            ``USEGAS="on"``
        * activate stars only, store metallicity, formation time, star foramtion rate (for gas particles)
            ``USESTARS="on"``
        * Calculate bulk black hole properties
            ``USEBH="on"``
        * stars and gas
            ``USEBARYONS="on"``
        * Cosmic ray quantities, currently nothing enabled
            ``USECOSMICRAYS="on"``
        * activate everything
            ``USEHYDRO="on"``

    * Adjust memory/max size of Binary KD Tree options, used to run search particles
        * if tree is going to be built on more than max 32 bit integer number particles then enable, memory footprint increases |
            ``LARGEKDTREE="on"``

.. topic:: Executable flags

    * Adjust **VELOCIraptor** operation
        * only calculate local density distribution for particles residing in field objects (but using all particles to estimate quantity). Default. |
            ``STRUCTUREDEN="on"``
        * or just use particles inside field objects, reducing cpu cycles but will bias estimates for particle in outer region of field structures, overrides STRUCTUREDEN |
            ``HALODEN="on"``
        * flag useful for zoom simulations with a high resolution region |
            ``ZOOMSIM="on"``
    * Adjust **TreeFrog** memory footprint
        * if particle ids are long integers |
            ``TREEFROGLONGIDS="on"``
        * if particle ids are unsigned |
            ``TREEFROGUNSIGNEDIDS="on"``
    * Enable if mpi domain is going to contain more than max 16 bit integer number of mpi processes
        ``LARGEMPIDOMAIN="on"``
    * Enable debugging
        ``DEBUG="on"``
