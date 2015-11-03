/*! \file stf.h
 *  \brief header file for the code
 */

#ifndef STFINDER_H
#define STFINDER_H

#include "allvars.h"
#include "proto.h"

#define STFVERSION 1.13

#endif

///-----------------------
/* ----------------------------------------------------------------------
   The rest of this file contains documentation for compiling and
   running the code, in a format appropriate for 'doxygen'.
   ----------------------------------------------------------------------
*/

/*! 
\mainpage  Reference documentation for VELOCIraptor (aka STF)
\author Pascal J. Elahi \n
        Sydney Institute for Astronomy \n
        University of Sydney \n
        Sydney, NSW, Australia, 2006 \n
        pelahi@physics.usyd.edu.au, pascaljelahi@gmail.com \n
\version 1.13
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

VELOCIraptor/STF is a massively parallel code for finding substructure in halos from N-body simulations. 
It was designed to identify tidal debris, an impossible task for most other structures finders since most are designed to find overdensities in configuration space. 
Tidal debris is 

This program reads a N-body output such as a tipsy/gadget file and finds (sub)structures. \n
It is primarily designed to use the VELOCIraptor/STF algorithm (which generates an outlier subset) but can also
use the 3DFOF and 6DFOF algorithm (effectively finds regions of high physical or phase-space density).\n
The program can be altered to read different formats by altering \ref proto.h, \ref io.cxx, specifically \ref ReadData
and providing the appropriate io.cxx file like \ref gadgetio.cxx. The overall flow of the program is outlined in \ref main.cxx 
and a brief description of the searches available is in \ref searching

A full account of the numerical algorithms employed by the code is given
in the code paper, <a href="http://adsabs.harvard.edu/abs/2011MNRAS.418..320E">Elahi et al. (2011)</a>
and detailed instructions for usage of the code are given in the included code documentation. 
For a discussion of tidal debris in cosmological simulations and a comparison of the few codes capable in principle
of identifying physically diffuse tidal debris see <a href="http://adsabs.harvard.edu/abs/2013MNRAS.433.1537E">Elahi et al. (2013)</a>.

This html-document serves as a cross-referenced documentation of the
source code itself - in fact, using the doxygen tool, the html-pages
have been produced from comments inlined in the source code. Apart from
the source-code documentation, a brief guide to code compilation is
given below, and under <b>Related Pages (see link on top)</b> you can find
a short guide to compile-time options of the code.

\section prelim Getting started

    Getting started is as simple as editing the Makefile.config file (see \ref STF-makeflags)
    typing make and then running the code (see \ref howtorun). 

\section install Compilation 

VELOCIraptor/STF needs the following non-standard libraries for compilation:

- \b GSL - the GNU scientific library. This open-source package can be
  obtained at http://www.gnu.org/software/gsl. STF
  needs this library for a few special function calls
  and for random number generation.

- \b libNBody - a scientific library included with STF (\ref libNBody).
  STF needs this library for a number of structures, classes, and methods it provides.

STF for parallel use may need the following non-standard libraries for compilation
depending on the compilation flags used:

- \b MPI - the Message Passing Interface (version 1.0 or higher). Many
  vendor supplied versions exist, in addition to excellent open source
  implementations, e.g.  MPICH
  (http://www-unix.mcs.anl.gov/mpi/mpich/) or LAM
  (http://www.lam-mpi.org/).

- \b OpenMP - the OpenMP API generally included with many compilers

- \b CUDA - NOT IMPLEMENTED YET.

One could in principle alter the code to use high precision libraries

- \b QD - for higher precision (like double-double and quad-quad). 
  (http://crd.lbl.gov/~dhbailey/mpdist/)

- \b ARPREC - also for arbitrarily high precision. 
  (http://crd.lbl.gov/~dhbailey/mpdist/)

Note that if any of the above libraries is not installed in standard
locations on your system, the \ref STF-makeflags "Makefile.config" provided with
the code may need slight adjustments. Similarly, compiler options,
particularly with respect to optimisations, may need adjustment to the
C++-compiler that is used.

The provided makefile is compatible with GNU-make, i.e. typing \b make or
\b gmake should then build the executable <b>STF</b>. If your site
does not have GNU-make, get it, or write your own makefile.

\section howtorun Running the code

A typical command to start the code looks like: \n \n

 <b> ./stf < args > </b> \n \n

with <em>OpenMP</em>, setting the environment variable OMP_NUM_THREADS=8 will then run 
the code with 8 threads in the openmp sections. 

with MPI: \n \n

 <b> mpirun -np 8 ./stf < args > </b> \n \n

This would start the code using 8 processors, assuming that the parallel
environment uses the <em>mpirun</em> command to start MPI
applications. Depending on the operating system, other commands may be
required for this task, e.g. <em>poe</em> on IBM/AIX machines. Note that
the code can in principle be started using an arbitrary number of
processors, but the communication algorithms will be most efficient for
powers of 2. 

Note that at the moment, mpirun assumes that a single structure can fit onto the shared
memory local to the mpi thread. If larger haloes are to be analyzed, it is suggested that
the iSingleHalo option be set to 1, and the analysis is done on a shared memory machine
with enough memory. A more complete version capable of handling large structures across
mpi domains that are then searched for substructures is in the works.

\section param Parameters
The code has several parameters that can be adjusted through a configuration file. The following commands
are accepted (more info can be found in \ref Options struct and \ref ui.cxx for user interface or the sample configuration file
in the examples directory). \n \n
    \subsection configfile Preferred interface: Use a configuration file
    for details of configuration options see \ref configopt
    \arg \b \e -i < input file > \ref Options.fname \n
    \arg \b \e -s < number of files per snapshot for gadget input, 0 for tipsy [default] > \ref Options.snapshotnum \n
    \arg \b \e -Z < number of files to read in parallel (when mpi is invoked) > \ref Options.nsnapread \n
    \arg \b \e -o < output base name (this can be overwritten by a configuration option in the config file. Suggestion would be to not use this option in the config file, use explicit command> \ref Options.outname \n
    \arg \b \e -C < Config file name (see \ref configopt for discussion of what is contained in this ascii parameter file) > \ref Options.pname \n
    \n

\section searching Altering the search for substructures

The algorithm searches for substructures in a specific fashion but if the user wishes to modify the search in some way,
this can be done by altering \ref search.cxx (along with \ref proto.h, \ref allvars.h and \ref fofalgo.h as necessary).
However, before altering the search it is useful to understand what searches are available.

First, note that implemented in the code are a variety of FOF criteria. Second, also implemented is an interative search.
Given the added complexity of the iterative search, an aside is necessary. First, the iterative search finds candidate objects
using the criteria passed, then <b><i> relaxes the criteria</i></b>, thus if an iterative search is used,
<b><em>one should use it with more restrictive fof criteria than one would normally do</em></b>. The iterative search
can correct for a bias present in large subhaloes <b><i>AND</i></b> also find extended portions (typically unbound portions)
of a substructure. As candidate tidal debris substructure can be split into several "groups" using the smaller
search window and have links when using the more relaxed criteria. If enough connections are present between
candidate substructures they are merged. The iterative parameters indicating how much the initial parameters
are increased by are listed in \ref Options

\section unbinding Unbinding when searching for tidal debris
As the algorithm is designed to identify dynamically distinct but not necessarily self-bound structures residing in a background 
(that is roughly in equilibrium), the code allows the user to specific the ratio between kinetic and potential energy of particles 
AND whether particles are removed from the background potential when unbinding. Typically codes of this nature, i.e. (sub)halo finders, 
<b><em>REQUIRE</em></b> and unbinding step to function properly. This is <b>not</b> the case here. Tests show that requiring the ratio 
|Epot|/Ekin >0.2 does not throw up too many suprious objects. One can be even more relaxed (and use <em>slightly</em> cpu resources) if 
one does not require particles that do not meet the energy criterion to be ignored when estimating the potential. The idea would be for
extended tidal debris with a very small loosely unbound core, the sea of tidal debris particles contributes negligibly to the potential 
and so updating this is not cost effective. 

But for the default user interested in only completely self-bound objects, a kinetic ratio of 1.0 and ignore the bacground potential are 
the options to use. 

\section outputs Outputs

    \subsection foflists Structure ID lists
    The code will produce a list of structure ids (two different formats). These outputs are outname.fof.grp
    and outname.sublevels.fof.grp (or pglist instead of fof.grp) which will list only haloes
    (or level one hierarchy objects) and all structures found respectively. Note that the fof.grp
    format is collected from all MPI threads and is only ascii. The pglist format will be written
    in parallel, one for each thread, such that one has pglist.threadID. This format can be ascii
    or binary. This format contains several outputs, the catalog_group file, the catalog_particles file
    and the catalog_particles.unbound file. 
    
    \subsection prop Structure properties
    The code can produce a ascii/binary file outname.prop containing the properties of the structures analyzed.

    \subsection hierarchy Field Structure / Substructure relationships
    The code will produce a list of the number of levels found in the hierarchy and
    then for each object, the id, its parent id (which is -1 if it is a field structure)
    and all direct substructures (number and substructure ids). This can be in ascii and binary
    formats and is written to outname.hierarchy.

    \subsection mergertrees Merger trees
    The halomerger tree code located within the analysis directory. Work in progress.

    \subsection baryonic_analysis Analysing baryons
    Work in progress.

    \section modifications Modifications/Searching for other types of Structures

The code can be modified to search for other types of substructures by altering the definition of an outlier
particle and the linking criteria used. For instance, to search a hydrodynamical simulation for dynamically
distinct gas substructures which are also metalicity outliers would require:
    \arg altering the \ref NBody::Particle class in \ref Particle.h and \ref Particle.cxx. The particle could be
    used to store extra quantities such as temperature, metalicity, ionization fraction, etc.
    \arg Then one would need to alter the quantities calculated for the background in \ref bgfield.cxx to
    calculate for instance, the mean mass or volume weighted metalicity and the variance
    \arg Then one would have to alter \ref localfield.cxx to calculate the local velocity distribution
    <b><em>AND the local metalicity distribution</em></b>.
    \arg One would need to alter \ref localbgcomp.cxx to compare the predicted bg to the actual
    local value for each quantity independently. (Here one might assume that the metalicity distribution is lognormal.
    \arg Finally, one would need to alter the search criterion used in \ref search.cxx (and \ref fofalgo.h). 

*/

/*! \page STF-makeflags  Makefile.config

A number of features of STF are controlled with compile-time options
in the makefile rather than by arguments. The Makefile.config file contains 
a list of all available compile-time options, with most of them commented out by default. 
To activate a certain feature, the corresponding parameter should be commented in, 
and given the desired value, where appropriate. Below, a brief guide to these options is
included.

<b>Important Note:</b> Whenever one of the compile-time options
described below is modified, a full recompilation of the code may be
necessary. To guarantee that this is done when a simple <b>make</b> is
specified, all source files have been specified in the Makefile as being
dependent on the Makefile and Makefile.config files. Alternatively, one can also issue the
command <b>make clean</b>, which will erase all object files, followed
by <b>make</b>.


\n
\section secmake1 Basic operation mode of code
- \b SCALING \n Set this if you want to scale the physical and velocity 
  linking lengths by bulk values of the halo to be searched. Only valid if a single halo is passed
- \b STRUCDEN \n Set this if you want to calculate the velocity density function used to find (sub)structures \e ONLY for particles resident in structure
- \b HALOONLYDEN \n Set this if you want to calculate the velocity density function used to find (sub)structures \e ONLY for particles resident in structure \em USING 
\em ONLY \em PARTICLES in the parent structure. \b STRUCDEN is overridden by \b HALOONLYDEN (technically these are incompatible with each other as HALOONLYDEN is faster but is biased and the number of particles for which the local distribution function density is calculated for is even more incomplete).
- \b NOMASS \n Set in \ref NBody::Particle contained in file \ref Particle.h reduces memory allocation
- \b MASSVAL \n Set in \ref NBody::Particle contained in file \ref Particle.h is the value returned by \ref NBody::Particle.GetMass
- \b LOWPRECISIONPOS \n Set in \ref NBody::Particle contained in file \ref Particle.h reduces memory allocation and sets phase-space position of particles to floats
- \b LONGIDS \n If this is set, the code assumes that particle-IDs are stored as 64-bit long integers. This is only really needed if you want
     to go beyond ~2 billion particles.

\n
\section secmake2 Parallel and computational flags
Note that for practical reasons, the combination of OpenMP/MPI only works on correctly setup environments where one can explicitly state how many MPI threads to start on a given node. Otherwise, many systems will fill up a node with mpi threads and do so till all asked for mpi threads are active. Consequently, the openmp threads spawned by the MPI threads will compete with the MPI threads on the same node and the other OpenMP threads started by other MPI threads. 
- \b USEOPENMP \n Code is compiled with openmp
- \b USEMPI \n Code is compiled with mpi. must also set the appropriate compiler
- \b SINGLEPRECISION, \b QUADPRECISION \n Code is compiled with double (normal), single and quad precision. must also compile 
  \ref libNBody contained in "NBodylib" with the same precision
- \b REDUCEMPIMEM \n If this flag is set, then the amount of memory allocated for MPI routines is reduced and one does not have to worry as much about the 
MPI factors \ref MPIProcFac & \ref MPIExportFac defined in \ref mpivar.h

\n
\section secmake3 Output options
- \b OUTPUTFIELDHALOSINCLUSIVE \n Similar to \b OUTPUTFIELDHALOS but here the data is inclusive, that is file list and associated properties produced for halos including any 
    substructures present.
- \b PGLISTOUTPUT \n If this is set, the particle catalogue is a particle group list format.

\n
\section secmake4 Things for special behaviour 
- \b TWODIMS \n Projection search, not yet implemented

\n
\section secmake4 IO flags when reading gadget or other formats
- \b GLONGID \n Gadget stores 64 long ids
- \b GDPOS \n Gadget stores double precision x,v 
- \b GSMASS \n Gadget stores single precision mass
- \b GSHEAD \n Gadget has version 2 headers between each block of data


*/
