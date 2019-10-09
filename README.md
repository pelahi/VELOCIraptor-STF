```

____   _______________.____    ________  _________ .___                      __
\   \ /   /\_   _____/|    |   \_____  \ \_   ___ \|   |___________  _______/  |_  ___________
 \   Y   /  |    __)_ |    |    /   |   \/    \  \/|   \_  __ \__  \ \____ \   __\/  _ \_  __ \
  \     /   |        \|    |___/    |    \     \___|   ||  | \// __ \|  |_> >  | (  <_> )  | \/
   \___/   /_______  /|_______ \_______  /\______  /___||__|  (____  /   __/|__|  \____/|__|
                   \/         \/       \/        \/                \/|__|

```
(Formerly)
```
  ____________________                     __
 /   _____/\__    ___/______ __ __   _____/  |_ __ _________   ____
 \_____  \   |    |  \_  __ \  |  \_/ ___\   __\  |  \_  __ \_/ __ \
 /        \  |    |   |  | \/  |  /\  \___|  | |  |  /|  | \/\  ___/
/_______  /  |____|   |__|  |____/  \___  >__| |____/ |__|    \___  >
        \/                              \/                        \/
___________.__            .___
\_   _____/|__| ____    __| _/___________
 |    __)  |  |/    \  / __ |/ __ \_  __ \
 |     \   |  |   |  \/ /_/ \  ___/|  | \/
 \___  /   |__|___|  /\____ |\___  >__|
     \/            \/      \/    \/
    ___   _______________________________ ___
   /  /  /   _____/\__    ___/\_   _____/ \  \
  /  /   \_____  \   |    |    |    __)    \  \
 (  (    /        \  |    |    |     \      )  )
  \  \  /_______  /  |____|    \___  /     /  /
   \__\         \/                 \/     /__/

```

![alt text](https://github.com/pelahi/VELOCIraptor-STF/blob/master/velociraptoricon.png)

# VELOCIraptor (formerly STructure Finder)

================================================================================================
## developed by:

    Pascal Jahan Elahi (continuously)
    Additional contributors:
    Rhys Poulton
    Rodrigo Canas

================================================================================================

## Content

This is brief description of the package. For details please see online documentation at
[readthedocs](https://velociraptor-stf.readthedocs.io/)

The repo contains the following directories

    src/        contains main source code for the algorithm
    doc/        contains Doxygen generated latex and html file of code
    examples/   contains examples of configuration files and how to run the code
    NBodylib/   submodule: contains library of objects and routines used by algorithm
    tools/      submodule: contains python tools of manipulating/reading output


================================================================================================

## Compiling (see documentation for more information)

VELOCIraptor uses CMake as its build tool. cmake is used to perform system-level checks,
like looking for libraries and setting up the rules for the build, and then generates the
actual build scripts in one of the supported build systems. Among other things, cmake supports
out-of-tree builds (useful to keep more than one build with different settings, and to avoid
cluttering the original source code directories) and several build system, like make and
ninja files.

VELOCIraptor uses submodules so if you have a fresh clone you can use

    git submodule update --init --recursive

to update the submodules use

    git submodule update --recursive --remote

Note that cmake will init the submodules by default

The simplest way of building is, standing on the root your repository, run cmake to produce
Makefiles and then compile with these steps:

    mkdir build
    cd build
    cmake .. # By default will generate Makefiles
    make all

There are a variety of options that can be invoked and these can be viewed using

    cmake -LH

(though this only works after having run cmake at least once)

Although documentation is present on the readthedocs site, extra documentation can be produced
by typing

    make doc

which will produce html and latex documents using Doxygen. This will be located in
doc/html/index.html and doc/latex/refman.tex

Note that VELOCIraptor and all variants do not support non-Unix environments. (Mac OS X is fine; Windows is not).

================================================================================================

## Running (see documentation for more information)

Running is as simple as

    ./bin/stf -i input -s nsnaportype -C configfile

a sample of a configuation file is in examples.
For mpi enabled executable

    mpirun -np mutipleoftwo ./bin/stf

Note that at the moment, mpirun assumes that a single structure can fit onto the shared
memory local to the mpi thread. If larger haloes are to be analyzed, it is suggested that
the iSingleHalo option be set to 1, and the analysis is done on a shared memory machine
with enough memory. A more complete version capable of handling large structures across
mpi domains that are then searched for substructures is in the works.

## Outputs

The code will produce several files related to the configuration options (.configuration)
input data (.siminfo), units information based on the input and the configuration options (.units)
and several files containing information about the structures identified. These files can be split
into separate files containing field objects (halos) and internal structures (subhalos) if desired
(set by configuration option). The files can be in several formats: ascii, binary (not recommended, HDF5.
These files are

    i) Properties File (.properties)
    Contains a variety of properties calculate for each structure identified and also contains
    some information regarding the relationship of this object to others (if object is a substructure,
    what is its hostHaloID).

    ii) Catalog_ files (.catalog_groups, .catalog_particles, .catalog_parttypes)
    Contains particle information to extract read particles from an input file or for tracking (ie: producing halo merger trees).
    The catalog_groups contains the sizes of groups and the offsets to read the associated .catalog_particles (.catalog_particles.unbound)
    .catalog_parttypes (.catalog_parttypes.unbound) which just listed the IDS and Types of particles belonging to groups.
    For examples of how to read this information, see the python tools included. When combined with raw particle data can
    be used for extra processing (such as calculating properties/profiles not calculated by default by VELOCIraptor).
    These files are also necessary if one wishes to construct "halo merger trees" or cross match haloes between catalogues.
    TreeFrog, an associated MPI+OpenMP tool, has inbuilt readers for VELOCIraptor output.

    iii) Field Structure / Substructure relationships (.hierarchy)
    Contains the substructure hierarchy information, such as the hostID (which is -1 if it is a field structure)
    an objects ID, number of direct substructures.

The code can also output a simple list which is particle id ordered that simply has the (sub)halo of a particle
(and is zero if particle doesn't belong to a list.) These outputs are outname.fof.grp. Note that the fof.grp
format is collected from all MPI threads and is only ascii output.


================================================================================================

## Altering IO for other file types (see documentation for more information)

Naturally, not all simulations will be in the io formats already written. An example of
several implemented io routines are in the src directory. The routine needs to load all
the appropriate date into a Particle array.

Currently VELOCIraptor can read Gadget (1,2), HDF, RAMSES, TIPSY, and Nchilada (alpha)

================================================================================================

## Associated analysis:

TreeFrog (Fomerly Halotree): This code is a separate repo

    git clone https://github.com/pelahi/TreeFrog.git

but is useful for checking halo catalogs, particularly the
examples/catalogcomparisontolerancecheck.py code that runs TreeFrog to compare two catalogs.
The goal of this code is as a test to see if a new catalog produced by velociraptor matches
a reference catalog.

================================================================================================
