/*! \file allvars.h
 *  \brief declares global variables and structures.
 *
 *  This file declares all global variables and structures. Further variables should be added here, and declared as
 *  \e \b extern. The actual existence of these variables is provided by the file \ref allvars.cxx. To produce
 *  \ref allvars.cxx from \ref allvars.h, do the following:
 *
 *     \arg Erase all \#define's, typedef's, and enum's
 *     \arg add \#include "allvars.h", delete the \#ifndef ALLVARS_H conditional
 *     \arg delete all keywords 'extern'
 *     \arg delete all struct definitions enclosed in {...}, e.g.
 *        "extern struct global_data_all_processes {....} All;"
 *        becomes "struct global_data_all_processes All;"
 */

#ifndef ALLVARS_H
#define ALLVARS_H

#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>
#include <string>
#include <getopt.h>
#include <sys/stat.h> 
#include <sys/timeb.h> 

///\name Include for NBodyFramework library.
//@{
///nbody code
#include <NBody.h>
///Math code
#include <NBodyMath.h>
///Binary KD-Tree code
#include <KDTree.h>
///Extra routines that analyze a distribution of particles
#include <Analysis.h>
//@}

// ///\name for checking the endian of floats
//#include <endianutils.h>

#ifdef USEOPENMP
#include <omp.h>
#endif

#ifdef USEMPI
#include <mpi.h>
///Includes external global variables used for MPI version of code
#include "mpivar.h"
#endif

using namespace std;
using namespace Math;
using namespace NBody;

//-- Structures and external variables


/// \name Specific structure type, allow for other types beside HALO
/// \todo will have to alter this to allow for subvois possibly even subfilaments (though I'm not certain what the hell this would be)
//@{
#define HALOSTYPE 10
#define WALLSTYPE 1
#define VOIDSTYPE 2
#define FILAMENTSTYPE 3
//@}

/// \name Cross-match types
//@{

//subsets made
///call \ref NsharedN1N2
#define NsharedN1N2 1 
#define NsharedN1 2
#define Nshared 3
#define Nsharedcombo 4
//@}

/// \name Input format
//@{
#define DSUSSING 1
#define DCATALOG 2
#define DNIFTY 3
//@}

/// \name Mapping functions
//@{
#define DNOMAP 0
#define DSIMPLEMAP 1
//@}


/// \name nifty format particle type for subselecting particle matching
#define NIFTYDMTYPE 1

/// \name defining particle types for subselecting particle matching when loading velociraptor group catalogs
//@{
#define ALLTYPEMATCH -1
#define GASTYPEMATCH 0
#define DMTYPEMATCH 1
#define DMGASTYPEMATCH -2
#define STARTYPEMATCH 4
//@}

/// Options structure stores useful variables that have user determined values which are altered by \ref GetArgs in \ref ui.cxx
struct Options
{
    ///\name filenames
    //@{
    char *fname,*outname;
    //@}
    ///number of snapshots
    int numsnapshots;
    ///number of particles
    long unsigned NumPart;
    ///total number of haloes across all snapshots    
    long unsigned TotalNumberofHalos;

    ///type of cross-match search
    int matchtype;
    /// store description of code
    char *description;
    ///flag for whether a simple cross comparsion between two catalogs is run or default merger tree produced
    int icatalog;

    ///output format
    int outputformat;

    ///cross match merit limit significance, that is match only when quantity is above MERITLIM*some measure of noise
    Double_t mlsig;

    ///io format flags for basic format, whether data split into multiple files, binary, and field objects separate
    //@{
    int ioformat,nmpifiles,ibinary,ifield;
    //@}

    //for mapping particle ids to index 
    //@{
    int imapping;
    void(*mappingfunc)(long unsigned int &);
    //@}

    ///for offseting the snapshot id value when writing to a file
    Int_t snapshotvaloffset;

    ///to fix ids for nifty project
    int idcorrectflag;

    ///particle type cross matching subselection
    int itypematch;

    Options()
    {
        numsnapshots=2;

        fname=outname=NULL;
        matchtype=NsharedN1N2;
        NumPart=512*512*512;
        TotalNumberofHalos=0;

        ioformat=DCATALOG;
        nmpifiles=0;
        ibinary=1;
        ifield=1;

        mlsig=0.1;

        imapping=DNOMAP;//no mapping
        mappingfunc=NULL;

        snapshotvaloffset=0;
        icatalog=0;
        idcorrectflag=0;
        outputformat=0;

        itypematch=ALLTYPEMATCH;
    }
};

/*! 
    Halo data structure
*/
struct HaloData{
    long unsigned haloID;
    long unsigned NumberofParticles;
    long unsigned *ParticleID;
    //Coordinate         *X;
    //Coordinate         *V;
    HaloData(long unsigned nhalos=0){
        NumberofParticles=nhalos;
        if (NumberofParticles>0) {
            ParticleID=new long unsigned[NumberofParticles];
            //X=new Coordinate[NumberofParticles];
            //V=new Coordinate[NumberofParticles];
        }
    }
    void Alloc(long unsigned nhalos=0){
        if (NumberofParticles>0){
            delete[] ParticleID;
            //delete[] X;
            //delete[] V;
        }
        NumberofParticles=nhalos;
        if (NumberofParticles>0) {
            ParticleID=new long unsigned[NumberofParticles];
            //X=new Coordinate[NumberofParticles];
            //V=new Coordinate[NumberofParticles];
        }
    }
    ~HaloData(){
        if (NumberofParticles>0){
            delete[] ParticleID;
            //delete[] X;
            //delete[] V;
        }
    }
};

struct HaloTreeData{
    Int_t numhalos;
    HaloData *Halo;
    HaloTreeData(){
        numhalos=0;
        Halo=NULL;
    }
    ~HaloTreeData(){
        if (numhalos>0) delete[] Halo;
    }
};

/*! 
    Structure used to keep track of a structure's progenitor/parent structure using temporal/spatial information
*/
struct ProgenitorData
{
    ///structure type and number in current level of hierarchy 
    int stype;
    int NumberofProgenitors;
    long unsigned* ProgenitorList;
    Double_t *Merit;
    ///store the fraction of shared particles
    Double_t *nsharedfrac;
    //points to the the head pfof address of the group and parent
    ProgenitorData(){
        NumberofProgenitors=0;
        ProgenitorList=NULL;
        Merit=NULL;
        nsharedfrac=NULL;
    }
    ~ProgenitorData(){
        if (NumberofProgenitors>0) {
            delete[] ProgenitorList;
            delete[] Merit;
            if (nsharedfrac!=NULL) delete[] nsharedfrac;
        }
    }
};

/*!
    Structure used to keep track of a structure's descendant/child structure using temporal/spatial information
*/
struct DescendantData
{
    ///structure type and number in current level of hierarchy
    int stype;
    int NumberofDescendants;
    long unsigned* DescendantList;
    Double_t *Merit;
    ///store the fraction of shared particles
    Double_t *nsharedfrac;
    //points to the the head pfof address of the group and parent
    DescendantData(){
        NumberofDescendants=0;
        DescendantList=NULL;
        Merit=NULL;
        nsharedfrac=NULL;
    }
    ~DescendantData(){
        if (NumberofDescendants>0) {
            delete[] DescendantList;
            delete[] Merit;
            if (nsharedfrac!=NULL) delete[] nsharedfrac;
        }
    }
};

#endif

