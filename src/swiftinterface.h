/*! \file swiftinterface.h
 *  \brief this file contains definitions and routines for interfacing with swift
 *
 */

#ifndef SWIFTINTERFACE_H
#define SWIFTINTERFACE_H

#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>
#include <string>
#include <getopt.h>
#include <sys/stat.h>

#ifdef USEMPI
#include <mpi.h>
#endif

///include the options structure via allvars
#include "allvars.h"
#include "proto.h"

///\name Include for NBodyFramework library.
//@{
///nbody code
#include <NBody.h>
///Math code
#include <NBodyMath.h>
///Binary KD-Tree code
#include <KDTree.h>
//@}

using namespace std;
using namespace Math;
using namespace NBody;

#ifdef SWIFTINTERFACE

///structure for interfacing with swift and extract cosmological information
namespace Swift {
    struct cosmoinfo {
        double atime, littleh, Omega_m, Omega_b, Omega_Lambda, Omega_cdm, w_de;
    };
    ///structure to store unit information of swift
    struct unitinfo {
        double lengthtokpc,velocitytokms,masstosolarmass,energyperunitmass,gravity,hubbleunit;
    };

    /* Structure to hold the location of a top-level cell. */
    //struct cell_loc {
    //
    //    /* Coordinates x,y,z */
    //    double loc[3];

    //} SWIFT_STRUCT_ALIGN;

    //store simulation info
    struct siminfo {
        double period, zoomhigresolutionmass, interparticlespacing, spacedimension[3];

        /* Number of top-level cells. */
        int numcells;

        /* Number of top-level cells in each dimension. */
        int numcellsperdim;

        /* Locations of top-level cells. */
        cell_loc *cellloc;

        /*! Top-level cell width. */
        double cellwidth[3];

        /*! Inverse of the top-level cell width. */
        double icellwidth[3];

        /*! Holds the node ID of each top-level cell. */
        int *cellnodeids;

        int icosmologicalsim;
    };

}

using namespace Swift;

///\defgroup external C style interfaces that can be called in swift N-body code.
//@{
//extern "C" void InvokeVelociraptor(const int num_gravity_parts, const int num_hydro_parts, struct gpart *gravity_parts, struct part *hydro_parts, double mpi_domain [3][2]);
///initialize velociraptor
extern "C" int InitVelociraptor(char* configname, char* outputname, Swift::cosmoinfo, Swift::unitinfo, Swift::siminfo);
///actually run velociraptor
extern "C" int InvokeVelociraptor(const size_t num_gravity_parts, const size_t num_hydro_parts, struct swift_vel_part *swift_parts, const int *cell_node_ids, char* outputname);
//@}

//extern KDTree *mpimeshtree;
//extern Particle *mpimeshinfo;

///global libvelociraptorOptions structure that is used when calling library velociraptor from swift
extern Options libvelociraptorOpt;

#endif

#endif
