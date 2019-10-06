/*! \file stf.h
 *  \brief header file for the code
 */

#ifndef STFINDER_H
#define STFINDER_H

#include "allvars.h"
#include "proto.h"

#define STFVERSION 1.30

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
        ICRAR \n
        University of Western Australia \n
        Crawley, WA, Australia, 6009 \n
        pascal.elahi@uwa.edu.au, pelahi@physics.usyd.edu.au, pascaljelahi@gmail.com \n
\version 1.30
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
Tidal debris and other forms of substructure within haloes/galaxies are identified using a 6D phase-space algorithm.
Structures are also analysed and bulk properties recorded in output files along with all particles that belong to structures.

This program runs by reading an input N-body snapshot such as a tipsy/gadget/HDF5 file and finds (sub)structures,
outputing files containing bulk properties and particles IDs belonging to (sub)structures. \n
It is primarily designed to use the VELOCIraptor/STF algorithm (which generates an outlier subset) but can also
use the 3DFOF and 6DFOF algorithm (effectively finds regions of high physical or phase-space density).\n
The program can be altered to read different formats by altering \ref proto.h, \ref io.cxx, specifically \ref ReadData
and providing the appropriate io.cxx file like \ref gadgetio.cxx. The overall flow of the program is outlined in \ref main.cxx.

A full account of the numerical algorithms employed by the code is given
in the code papers,
<a href="http://adsabs.harvard.edu/abs/2011MNRAS.418..320E">Elahi et al., (2011)</a>,
<a href="https://ui.adsabs.harvard.edu/abs/2019PASA...36...21E/abstract">Elahi et al., (2019)</a>
and detailed instructions for usage of the code are given in the included code documentation.
For a discussion of tidal debris in cosmological simulations and a comparison of the few codes capable in principle
of identifying physically diffuse tidal debris see <a href="http://adsabs.harvard.edu/abs/2013MNRAS.433.1537E">Elahi et al. (2013)</a>.

This html-document serves as a cross-referenced documentation of the
source code itself - in fact, using the doxygen tool, the html-pages
have been produced from comments inlined in the source code. Apart from
the source-code documentation, documentatino can be found
<a href="https://velociraptor-stf.readthedocs.io/en/latest/">online</a>.
*/
