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
    //determine the number of splits in each dimension
    Nsplit=log((float)NProcs)/log(2.0);
    mpi_ideltax[0]=0;mpi_ideltax[1]=1;mpi_ideltax[2]=2;
    isplit=0;
    for (j=0;j<3;j++) mpi_nxsplit[j]=0;
    for (j=0;j<Nsplit;j++) {
        mpi_nxsplit[mpi_ideltax[isplit++]]++;
        if (isplit==3) isplit=0;
    }
    for (j=0;j<3;j++) mpi_nxsplit[j]=pow(2.0,mpi_nxsplit[j]);
    //for all the cells along the boundary of axis with the third split axis (smallest variance if I actually tried load balancing with a KD-Tree)
    //set the domain limits to the sims limits
    int ix=mpi_ideltax[0],iy=mpi_ideltax[1],iz=mpi_ideltax[2];
    int mpitasknum;
    for (j=0;j<mpi_nxsplit[iy];j++) {
        for (i=0;i<mpi_nxsplit[ix];i++) {
            mpitasknum=i+j*mpi_nxsplit[ix]+0*(mpi_nxsplit[ix]*mpi_nxsplit[iy]);
            mpi_domain[mpitasknum].bnd[iz][0]=mpi_xlim[iz][0];
            mpitasknum=i+j*mpi_nxsplit[ix]+(mpi_nxsplit[iz]-1)*(mpi_nxsplit[ix]*mpi_nxsplit[iy]);
            mpi_domain[mpitasknum].bnd[iz][1]=mpi_xlim[iz][1];
        }
    }
    //here for domains along second axis
    for (k=0;k<mpi_nxsplit[iz];k++) {
        for (i=0;i<mpi_nxsplit[ix];i++) {
            mpitasknum=i+0*mpi_nxsplit[ix]+k*(mpi_nxsplit[ix]*mpi_nxsplit[iy]);
            mpi_domain[mpitasknum].bnd[iy][0]=mpi_xlim[iy][0];
            mpitasknum=i+(mpi_nxsplit[iy]-1)*mpi_nxsplit[ix]+k*(mpi_nxsplit[ix]*mpi_nxsplit[iy]);
            mpi_domain[mpitasknum].bnd[iy][1]=mpi_xlim[iy][1];
        }
    }
    //finally along axis with largest variance
    for (k=0;k<mpi_nxsplit[iz];k++) {
        for (j=0;j<mpi_nxsplit[iy];j++) {
            mpitasknum=0+j*mpi_nxsplit[ix]+k*(mpi_nxsplit[ix]*mpi_nxsplit[iy]);
            mpi_domain[mpitasknum].bnd[ix][0]=mpi_xlim[ix][0];
            mpitasknum=(mpi_nxsplit[ix]-1)+j*mpi_nxsplit[ix]+k*(mpi_nxsplit[ix]*mpi_nxsplit[iy]);
            mpi_domain[mpitasknum].bnd[ix][1]=mpi_xlim[ix][1];
        }
    }
    //here use the three different histograms to define the boundary
    int start[3],end[3];
    Double_t bndval[3],binsum[3],lastbin;
    start[0]=start[1]=start[2]=0;
    for (i=0;i<mpi_nxsplit[ix];i++) {
        bndval[0]=(mpi_xlim[ix][1]-mpi_xlim[ix][0])*(Double_t)(i+1)/(Double_t)mpi_nxsplit[ix];
        if(i<mpi_nxsplit[ix]-1) {
        for (j=0;j<mpi_nxsplit[iy];j++) {
            for (k=0;k<mpi_nxsplit[iz];k++) {
                //define upper limit
                mpitasknum=i+j*mpi_nxsplit[ix]+k*(mpi_nxsplit[ix]*mpi_nxsplit[iy]);
                mpi_domain[mpitasknum].bnd[ix][1]=bndval[0];
                //define lower limit
                mpitasknum=(i+1)+j*mpi_nxsplit[ix]+k*(mpi_nxsplit[ix]*mpi_nxsplit[iy]);
                mpi_domain[mpitasknum].bnd[ix][0]=bndval[0];
            }
        }
        }
        //now for secondary splitting
        if (mpi_nxsplit[iy]>1)
        for (j=0;j<mpi_nxsplit[iy];j++) {
            bndval[1]=(mpi_xlim[iy][1]-mpi_xlim[iy][0])*(Double_t)(j+1)/(Double_t)mpi_nxsplit[iy];
            if(j<mpi_nxsplit[iy]-1) {
            for (k=0;k<mpi_nxsplit[iz];k++) {
                mpitasknum=i+j*mpi_nxsplit[ix]+k*(mpi_nxsplit[ix]*mpi_nxsplit[iy]);
                mpi_domain[mpitasknum].bnd[iy][1]=bndval[1];
                mpitasknum=i+(j+1)*mpi_nxsplit[ix]+k*(mpi_nxsplit[ix]*mpi_nxsplit[iy]);
                mpi_domain[mpitasknum].bnd[iy][0]=bndval[1];
            }
            }
            if (mpi_nxsplit[iz]>1)
            for (k=0;k<mpi_nxsplit[iz];k++) {
                bndval[2]=(mpi_xlim[iz][1]-mpi_xlim[iz][0])*(Double_t)(k+1)/(Double_t)mpi_nxsplit[iz];
                if (k<mpi_nxsplit[iz]-1){
                mpitasknum=i+j*mpi_nxsplit[ix]+k*(mpi_nxsplit[ix]*mpi_nxsplit[iy]);
                mpi_domain[mpitasknum].bnd[iz][1]=bndval[2];
                mpitasknum=i+j*mpi_nxsplit[ix]+(k+1)*(mpi_nxsplit[ix]*mpi_nxsplit[iy]);
                mpi_domain[mpitasknum].bnd[iz][0]=bndval[2];
                }
            }
        }
    }
    }
}

///reads HDF file to determine number of particles in each MPIDomain
void MPINumInDomainNchilada(Options &opt)
{
    if (NProcs>1) {
    MPIDomainExtentNchilada(opt);
    MPIDomainDecompositionNchilada(opt);
    MPIInitialDomainDecomposition();

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
