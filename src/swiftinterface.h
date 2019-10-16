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
    /// store basic cosmological information
    struct cosmoinfo {
  /*! Current expansion factor of the Universe. (cosmology.a) */
  double atime;

  /*! Reduced Hubble constant (H0 / (100km/s/Mpc) (cosmology.h) */
  double littleh;

  /*! Matter density parameter (cosmology.Omega_m) */
  double Omega_m;

  /*! Radiation density parameter (cosmology.Omega_r) */
  double Omega_r;

  /*! Neutrino density parameter (0 in SWIFT) */
  double Omega_nu;

  /*! Neutrino density parameter (cosmology.Omega_k) */
  double Omega_k;

  /*! Baryon density parameter (cosmology.Omega_b) */
  double Omega_b;

  /*! Radiation constant density parameter (cosmology.Omega_lambda) */
  double Omega_Lambda;

  /*! Dark matter density parameter (cosmology.Omega_m - cosmology.Omega_b) */
  double Omega_cdm;

  /*! Dark-energy equation of state at the current time (cosmology.w)*/
  double w_de;
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

    ///store simulation info
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

        /// whether run is cosmological simulation
        int icosmologicalsim;
        /// running a zoom simulation
        int izoomsim;
        /// flags indicating what types of particles are present
        int idarkmatter, igas, istar, ibh, iother;
        #ifdef NOMASS
        /// store mass as a option if saving memory by not storing masses as running uniform box
        double mass_uniform_box;
        #endif
    };

    struct groupinfo {
      int index;
      long long groupid;
    };
}


using namespace Swift;

///\defgroup external C style interfaces that can be called in swift N-body code.
//@{
//extern "C" int InitVelociraptor(char* configname, char* outputname, Swift::cosmoinfo, Swift::unitinfo, Swift::siminfo, const int numthreads);
///initialize velociraptor, check configuration options,
extern "C" int InitVelociraptor(char* configname, Swift::unitinfo, Swift::siminfo, const int numthreads);
///initialize velociraptor based on extra possible outputs requested by swift in addition to standard invocation, check configuration options,
extern "C" int InitVelociraptorExtra(const int iextra, char* configname, Swift::unitinfo, Swift::siminfo, const int numthreads);
///actually run velociraptor
extern "C" Swift::groupinfo * InvokeVelociraptor(const int snapnum, char* outputname,
    Swift::cosmoinfo, Swift::siminfo,
    const size_t num_gravity_parts, const size_t num_hydro_parts, const size_t num_star_parts,
    struct swift_vel_part *swift_parts, int *cell_node_ids,
    const int numthreads, const int ireturngroupinfoflag, int *const numingroups
);
extern "C" Swift::groupinfo * InvokeVelociraptorHydro(const int snapnum, char* outputname,
    Swift::cosmoinfo, Swift::siminfo,
    const size_t num_gravity_parts, const size_t num_hydro_parts, const size_t num_star_parts, const size_t num_bh_parts,
    struct swift_vel_part *swift_parts, int *cell_node_ids,
    const int numthreads, const int ireturngroupinfoflag, int *const numingroups,
    struct swift_vel_gas_part *swift_gas_parts = NULL,
    struct swift_vel_star_part *swift_star_parts = NULL,
    struct swift_vel_bh_part *swift_bh_parts = NULL
);
extern "C" Swift::groupinfo * InvokeVelociraptorExtra(const int iextra,
    const int snapnum, char* outputname,
    Swift::cosmoinfo, Swift::siminfo,
    const size_t num_gravity_parts, const size_t num_hydro_parts, const size_t num_star_parts,
    struct swift_vel_part *swift_parts, int *cell_node_ids,
    const int numthreads, const int ireturngroupinfoflag, int *const numingroups
);
extern "C" Swift::groupinfo * InvokeVelociraptorHydroExtra(const int iextra,
    const int snapnum, char* outputname,
    Swift::cosmoinfo, Swift::siminfo,
    const size_t num_gravity_parts, const size_t num_hydro_parts,
    const size_t num_star_parts, const size_t num_bh_parts,
    struct swift_vel_part *swift_parts, int *cell_node_ids,
    const int numthreads, const int ireturngroupinfoflag, int *const numingroups,
    struct swift_vel_gas_part *swift_gas_parts = NULL,
    struct swift_vel_star_part *swift_star_parts = NULL,
    struct swift_vel_bh_part *swift_bh_parts = NULL
);
///set simulation information
void SetVelociraptorSimulationState(Swift::cosmoinfo, Swift::siminfo);
//@}


///global libvelociraptorOptions structure that is used when calling library velociraptor from swift
extern Options libvelociraptorOpt;
extern Options libvelociraptorOptbackup;
extern Options libvelociraptorOptextra[10];

void CheckSwiftTasks(string message, const Int_t n, Particle *p);
#endif

#endif
