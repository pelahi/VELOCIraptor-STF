/*! \file NBody.h
 *  \brief header file for the NBody name space

    \namespace NBody
    \brief A namespace for user defined classes and routines useful for N-body simulations contained in libNBody

*/

#ifndef NBODY_H
#define NBODY_H

#include <iostream>
#include <fstream>
#include <cstring>
#include <cstdio>
#include <cmath>
#include <cstdlib>
#include <NBodyMath.h>
#include <Particle.h>
#include <System.h>

#ifdef USEOPENMP
#include <omp.h>
#endif

#endif // NBODY_H

/*! 
\mainpage  Reference documentation for libNBody
\defgroup libNBody NBody library
\author Pascal J. Elahi \n
        pascaljelahi@gmail.com \n
\version 1.0
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
    This library contains mathematical routines defined in the \ref Math namespace, useful classes for N-body simulations in the \ref NBody namespace.

\section prelim Getting started

\section install Compilation 

\b libNBody needs the following non-standard libraries for compilation:

- \b GSL - the GNU scientific library. This open-source package can be
  obtained at http://www.gnu.org/software/gsl. STF 
  needs this library for a few special function calls
  and for random number generation.

- \b QD - for higher precision (like double-double and quad-quad). 
  (http://crd.lbl.gov/~dhbailey/mpdist/)

- \b ARPREC - also for arbitrarily high precision. 
  (http://crd.lbl.gov/~dhbailey/mpdist/)

Note that if any of the above libraries is not installed in standard
locations on your system, the \ref libNBody-makeflags "Makefile.config" provided with
the code may need slight adjustments. Similarly, compiler options,
particularly with respect to optimisations, may need adjustment to the
C++-compiler that is used.

The provided makefile is compatible with GNU-make, i.e. typing \b make or
\b gmake should then build the various libraries and install them in the library directory
defined in the Makefile.config. If your site does not have GNU-make, get it, or write your own makefile.

\section howtorun Including the library 
Including the library to compile other codes is the same as including any other library.
Just call -lMath -lNBody -lKD -lAnalysis -lInitCond -lCosmology 
as necessary
and compile with the appropriate library and include directories

*/

/*! \page libNBody-makeflags  Makefile.config of libNBody

A number of features of libNBody are controlled with compile-time options
in the makefile. The Makefile.config file contains a list of all available compile-time options.
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
\n
\section secmake1 Parallel and computational flags
- \b USEOPENMP \n Code is compiled with openmp. This flag is only relevant to certain routines in \ref NBody::System and \ref NBody::KDTree

- \b QUADPRECISION \n Code is compiled with quad precision. 

*/

