/*! \file TreeFrog.h
 *  \brief header file for the code
 */

#ifndef TREEFROG_H
#define TREEFROG_H

#include "allvars.h"
#include "proto.h"

#define TREEFROG 1.20

#ifdef USEOPENMP
#include <omp.h>
#endif

#endif

/*!
\mainpage  Reference documentation for TreeFrog (aka HaloMergerTree)
\author Pascal J. Elahi \n
\version 1.2
\section license License

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License as
published by the Free Software Foundation; either version 2 of
the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
General Public License for more details at
http://www.gnu.org/copyleft/gpl.html

\section Description Description

TreeFrog is a massively parallel code for building a halo (or in fact any structure from voids to filaments) temporal merger tree.

This program reads as input a file containing a list of files containing particle lists produced by structure finders. These files can
be in several formats but critically list for each structure the particles associated with it and their IDs. These IDs are assumed to be
temporally constant and unique, allowig them to be used to trace structures across time.
The program can be altered to read different formats by altering \ref proto.h, \ref allvars.h, \ref ui.cxx and \ref io.cxx, specifically \ref ReadData
and providing the appropriate interface in \ref otherio.cxx. By default it is optimised to read output from <a href="../../../../doc/html/index.html">VELOCIraptor</a>
The overall flow of the program is outlined in \ref main.cxx and a brief description of the types of linking available is in \ref linking

A full account of the numerical algorithms employed by the code is given in the halo tree comparison paper, <a href="http://adsabs.harvard.edu/abs/2013MNRAS.436..150S">Srisawat et al. (2013)</a>

\section prelim Getting started

    Getting started is as simple as compiling <a href="../../../../doc/html/index.html">VELOCIraptor</a>.

\section install Compilation

Like <a href="../../../../doc/html/index.html">VELOCIraptor</a>, TreeFrog needs the following non-standard libraries for compilation:

- \b libNBody - a scientific library included with STF (\ref libNBody).
  STF needs this library for a number of structures, classes, and methods it provides.

Also needs:

- \b GSL - the GNU scientific library. This open-source package can be
  obtained at http://www.gnu.org/software/gsl. STF
  needs this library for a few special function calls
  and for random number generation.

  For parallel use may need the following non-standard libraries for compilation
depending on the compilation flags used:

- \b MPI - the Message Passing Interface (version 1.0 or higher). Many
  vendor supplied versions exist, in addition to excellent open source
  implementations, e.g.  MPICH
  (http://www-unix.mcs.anl.gov/mpi/mpich/) or LAM
  (http://www.lam-mpi.org/).

- \b OpenMP - the OpenMP API generally included with many compilers

Note that if any of the above libraries is not installed in standard
locations on your system, the \ref STF-makeflags "Makefile.config" provided with
the code may need slight adjustments. Similarly, compiler options,
particularly with respect to optimisations, may need adjustment to the
C++-compiler that is used.

The provided makefile is compatible with GNU-make, i.e. typing \b make or
\b gmake should then build the executable <b>treefrog</b>. If your site
does not have GNU-make, get it, or write your own makefile.

\section howtorun Running the code

A typical command to start the code looks like: \n \n

 <b> ./treefrog < args > </b> \n \n

with <em>OpenMP</em>, setting the environment variable OMP_NUM_THREADS=8 will then run
the code with 8 threads in the openmp sections.

with MPI: \n \n

 <b> mpirun -np 8 ./treefrog < args > </b> \n \n

This would start the code using 8 processors, assuming that the parallel
environment uses the <em>mpirun</em> command to start MPI
applications. Depending on the operating system, other commands may be
required for this task, e.g. <em>poe</em> on IBM/AIX machines. Note that
the code can be started using an arbitrary number of
processors, however one should be careful given the memory that is required to run with many
OpenMP threads.

Note that at the moment, mpirun assumes that all the past snapshots need to generate temporal links for a given current snapshot can
be stored locally on the local memory. That is the full halo (strutcure) catalog for N+1 snapshots is processed, where N is the
number of snapshots used to identify candidate links.

\section uiargs Interface
The code is passed command line arguments. More info can be found in \ref args, \ref Options struct and \ref ui.cxx \n
It is important to note that the input file is meant to contain a list of files that are to be opened with the order
of first file earliest snapshot, last latest. Example would be \n
\b snap_000.halo_particles \n
\b snap_001.halo_particles \n
where these files are a list of particle ids (and possibly type information too) of particles in structures that you
wish to track. The only alteration to the name is when reading in <a href="../../../../doc/html/index.html">VELOCIraptor</a> output. This output is structured
into a specific set of files, .properties, .catalog_groups, .catalog_particles, .catalog_parttypes. In his case the file
list should contain a list of the base names, ex: \n
\b snap_000.VELOCIraptor \n
\b snap_001.VELOCIraptor \n
where for snap_000.VELOCIraptor, there is snap_000.VELOCIraptor.properties, snap_000.VELOCIraptor.catalog_groups, snap_000.VELOCIraptor.catalog_particles, snap_000.VELOCIraptor.catalog_particles.unbound
(snap_000.VELOCIraptor.catalog_parttypes, snap_000.VELOCIraptor.catalog_parttypes.unbound if different particle types processed in halo catalog).


\section linking Linking Haloes
The code is designed to link halo (structures) in one catalog with those listed in another using the particle IDs.
Technically the code is a particle correlator which calculates for every halo in one catalog,
a merit of how good a match each halo in another catalog is via a merit function. The standard function is


\section outputs Outputs

    \subsection mergertree Merger Tree
    The code will produce an output that contains a simple header with
    \b - number of snapshots (catalogs) processed
    \b - a description of the code opertaion
    \b - total number of haloes in all catalogs
    This is followed by info for each snapshot processed.
    \b - snapshot_value number_of_haloes_in_snapshot
    \b - halo tree data
    This tree data lists for each halo
    \b - haloID N
    \b - IdofProgenitor_0

    \b - IdofProgenitor_N
    where N is the Number_of_progenitors.
    additonal information such as the merit value for each progenitor may also be listed.

    (see \ref io.cxx)
*/
