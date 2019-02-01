/*! \file mpivar.h
 *  \brief this file declares variables used for mpi compilation.
 */

#ifndef MPIVAR_H
#define MPIVAR_H

#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>
#include <string>
#include <getopt.h>
#include <sys/stat.h>

#include <mpi.h>

//
// Limits the value of x to be between a and b
//
//Only wraps once. If x > 2b, the returned value will be larger than b.
//Similarly for x < -b.
//
#define box_wrap(x, a, b)                                \
    ({                                                     \
     const __typeof__(x) _x = (x);                        \
     const __typeof__(a) _a = (a);                        \
     const __typeof__(b) _b = (b);                        \
     _x < _a ? (_x + _b) : ((_x >= _b) ? (_x - _b) : _x); \
     })

/* Convert cell location to ID. */
#define cell_getid(cdim, i, j, k) \
      ((int)(k) + (cdim)[2] * ((int)(j) + (cdim)[1] * (int)(i)))

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

///default buffer size (here in number of particles to send in one go)
#define MPIPartBufSize 100000
///size of largest MPI chunck in bytes that can be sent in one go (here set by MPI count argument, which is max int)
#define LOCAL_MAX_MSGSIZE 2147483647L
///Nlocal maximum initially set to nbodies/NProc*MPProcFac, represents maximum load imbalance
#define MPIProcFac 1.0
///maximum export factor used to allocate temporary export particle arrays, thus in any given step (except io) where
///particles are sent from one MPI thread to another, assume at most this factor of particle will need to be sent
#define MPIExportFac 0
#define MAXNNEXPORT 32

///define a type to store the maxium number of mpi tasks
#ifdef HUGEMPI
typedef int short_mpi_t;
#else
typedef short short_mpi_t;
#endif

/// \name definitions for MPI types
//@{
#ifdef LONGINT
#define MPI_Int_t MPI_LONG_LONG_INT
#define MPI_UInt_t MPI_UNSIGNED_LONG_LONG
#else
#define MPI_Int_t MPI_INT
#define MPI_UInt_t MPI_UNSIGNED
#endif
#ifdef SINGLEPRECISION
#define MPI_Real_t MPI_FLOAT
#else
#define MPI_Real_t MPI_DOUBLE
#endif

//@}

/// \name definitions for MPI message flags
//@{

///flag for IO particle exchange
#define TAG_IO_A 1
#define TAG_IO_B 2

///flag for FOF particle exchange
#define TAG_FOF_A 10
#define TAG_FOF_B 11
#define TAG_FOF_C 12
#define TAG_FOF_D 13
#define TAG_FOF_E 14
#define TAG_FOF_F 15

///flag for NN particle exchange
#define TAG_NN_A 20
#define TAG_NN_B 21

///flag for Grid data exchange
#define TAG_GRID_A 31
#define TAG_GRID_B 31
#define TAG_GRID_C 32

///flags for Extended output exchange
#define TAG_EXTENDED_A 100
#define TAG_EXTENDED_B 200


///flags for swift information exchange
#define TAG_SWIFT_A 1000
//@}

/// \name for mpi tasks and domain construction
//@{
extern int ThisTask, NProcs;
extern Int_t Ntotal, NExport, NImport, Nlocal, Nlocaltot, Ngridtotal, Ngridlocal;
///array to store baryons in all, gas, star particles, BH particles, etc.
extern Int_t Ntotalbaryon[NBARYONTYPES], Nlocalbaryon[NBARYONTYPES];
///to store the amount of available memory
extern Int_t Nmemlocal,Noldlocal,Nmemlocalbaryon;

extern Double_t mpi_xlim[3][2],mpi_dxsplit[3];
extern int mpi_nxsplit[3], mpi_ideltax[3];
extern Double_t mpi_period;

///store the boundaries of the processors spatial domain
extern struct MPI_Domain {
    Double_t bnd[3][2];
} *mpi_domain;
///
//@}

/// \name for mpi FOF search
//@{
///array that stores number of particles
extern Int_t *mpi_nlocal;
///Array which is [NProcs*NProcs] that stores number of particles to be exported in FOF search from processor i to j via [i+NProcs*j]
extern Int_t *mpi_nsend;
///local array that stores a particles global id;
extern Int_t *mpi_idlist;
///local array that stores a particles global index list of input file(s);
///\todo must implement this array to be used in output produced by \ref WritePGListIndex. This index based output can be useful.
extern Int_t *mpi_indexlist;
///local array that stores the thread to which a particle's fof group belongs
extern short_mpi_t *mpi_foftask;
///array that stores number of groups, need for properly setting Ids and broadcasting data
extern Int_t *mpi_ngroups;
///array that is used by task zero to collect the group ids of every particle so that it can be written to a file.
extern Int_t *mpi_pfof;
///array that is used to indicate particle must be sent across an mpi_domain, stores the mpi thread num particle is to be sent to
///\todo issue is that particles can be sent to more than one mpi thread so have to store Nlocal*NProcs which can be too large
///since the idea is to use GETNUMEXPORT to reduce mem costs
extern int *mpi_part_send_domain;
///store the max group id across all threads and the group id offset to be used with unlinked particles
extern Int_t mpi_maxgid,mpi_gidoffset;
///structure facilitates linking across mpi threads
extern struct fofdata_in
{
    //Particle Part();
    //FLOAT Hsml;
    Int_t iGroup;
    short_mpi_t iGroupTask;
    Int_t iLen;
    Int_t Index;
    Int_t Task;
}
*FoFDataIn, *FoFDataGet;
///Particle arrays used for transmitting data between mpi threads
extern Particle *PartDataIn, *PartDataGet;
///Particle arrays that allow allocation of memory need when deallocating local particle arrays that then need to be reassigned
extern Particle *mpi_Part1,*mpi_Part2;

///structure facilitates passing ids across threads so that determine number of groups and reorder group
///ids so that reflects group size
extern struct fofid_in {
    Particle p;
    Int_t iGroup;
    Int_t ID;
    Int_t Task;
    Int_t Index;
} *FoFGroupDataLocal, *FoFGroupDataExport;

///structure facilitates NN search across mpi threads
extern struct nndata_in
{
    //Int_t Index;
    Int_t ToTask,FromTask;
    Coordinate Pos, Vel;
    Double_t R2,V2;//NNR2[MAXNNEXPORT],NNV2[MAXNNEXPORT];
}
*NNDataIn, *NNDataGet;
//extern Particle *NNPartReturn, *NNPartReturnLocal;

///For transmitting grid data
//@{
extern struct GridCell *mpi_grid;
extern Coordinate *mpi_gvel;
extern Matrix *mpi_gveldisp;
//@}

/*extern struct fofdata_out
{
    FLOAT Distance;
    int   iGroup;
    int   iGroupTask;
}
*FoFDataResult, *FoFDataPartialResult;*/
//store MinSize as when using mpi prior to stitching use min of 2;
extern Int_t MinNumMPI,MinNumOld;

//@}


#endif
