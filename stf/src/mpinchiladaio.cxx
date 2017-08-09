/*! \file mpinchiladaio.cxx
 *  \brief this file contains routines used with MPI compilation and nchilada io and domain construction.
 *  \todo this is incomplete, missing reading of the period of the cosmological volume and have also
 *  not yet implemented the determination of the number of particles in the mpi domain
 */

#if defined(USEMPI)

//-- For MPI

#include "stf.h"
#include "nchiladaitems.h"

/// \name Nchilada Domain decomposition
//@{

/*!
    Determine the domain decomposition.\n
    Here the domains are constructured in data units
    only read tasks should call this routine. It is tricky to get appropriate load balancing and correct number of particles per processor.\n

    I could use recursive binary splitting like kd-tree along most spread axis till have appropriate number of volumes corresponding  to number of processors. Or build a Peno-Hilbert space filling curve.

    BUT for identifying particles on another mpi domain, simple enough to split volume into rectangular cells.
*/
#ifdef USEXDR

///Determine Domain for Nchilada input
void MPIDomainExtentNchilada(Options &opt){
    char buf[2000];
    if (ThisTask==0) {
        if(opt.num_files>1) sprintf(buf,"%s.0.hdf5",opt.fname);
        else sprintf(buf,"%s.hdf5",opt.fname);
    }
    //There may be issues with particles exactly on the edge of a domain so before expanded limits by a small amount
    //now only done if a specific compile option passed
#ifdef MPIEXPANDLIM
    for (int i=0;i<3;i++) {
        Double_t dx=0.001*(mpi_xlim[i][1]-mpi_xlim[i][0]);
        mpi_xlim[i][0]-=dx;mpi_xlim[i][1]+=dx;
    }
#endif

    //make sure limits have been found
    MPI_Barrier(MPI_COMM_WORLD);
    if (NProcs==1) {
        for (int i=0;i<3;i++) {
            mpi_domain[ThisTask].bnd[i][0]=mpi_xlim[i][0];
            mpi_domain[ThisTask].bnd[i][1]=mpi_xlim[i][1];
        }
    }
}

void MPIDomainDecompositionNchilada(Options &opt){
    Int_t i,j,k,n,m;
    int Nsplit,isplit;

    if (ThisTask==0) {
    }
}

///reads HDF file to determine number of particles in each MPIDomain
void MPINumInDomainNchilada(Options &opt)
{
    if (NProcs>1) {
    MPIDomainExtentNchilada(opt);
    MPIInitialDomainDecomposition();
    MPIDomainDecompositionNchilada(opt);

    Int_t i,j,k,n,nchunk;
    char buf[2000];
    MPI_Status status;

    Int_t Nlocalbuf,ibuf=0,*Nbuf, *Nbaryonbuf;
    int *ireadfile,*ireadtask,*readtaskID;
    ireadtask=new int[NProcs];
    readtaskID=new int[opt.nsnapread];
    ireadfile=new int[opt.num_files];
    MPIDistributeReadTasks(opt,ireadtask,readtaskID);

    Nbuf=new Int_t[NProcs];
    Nbaryonbuf=new Int_t[NProcs];
    for (j=0;j<NProcs;j++) Nbuf[j]=0;
    for (j=0;j<NProcs;j++) Nbaryonbuf[j]=0;

    //now read position information to determine number of particles per processor

    //now having read number of particles, run all gather
    Int_t mpi_nlocal[NProcs];
    MPI_Allreduce(Nbuf,mpi_nlocal,NProcs,MPI_Int_t,MPI_SUM,MPI_COMM_WORLD);
    Nlocal=mpi_nlocal[ThisTask];
    if (opt.iBaryonSearch) {
        MPI_Allreduce(Nbaryonbuf,mpi_nlocal,NProcs,MPI_Int_t,MPI_SUM,MPI_COMM_WORLD);
        Nlocalbaryon[0]=mpi_nlocal[ThisTask];
    }
    }
}

//@}
#endif

#endif
