/*! \file mpiroutines.cxx
 *  \brief this file contains routines used with MPI compilation. 

    MPI routines generally pertain to domain decomposition or to specific MPI tasks that determine what needs to be broadcast between various threads.

 */

#ifdef USEMPI

//-- For MPI

#include "halomergertree.h"

///load balance how snapshots are split across mpi domains
void MPILoadBalanceSnapshots(Options &opt){
    /*
    //now if there are too many tasks for the number of snapshots read, exit
    //must decided how best to do this
    if (NProcs>(0.5*opt.numsnapshots/opt.numsteps)) {
        if (ThisTask==0) {
            cerr<<"Too many mpi threads for number of snapshots and number of linking steps desired. "<<endl;
            cerr<<"Suggested number is less than "<<floor(0.5*opt.numsnapshots/opt.numsteps)<<endl;
            cerr<<"Exiting"<<endl;
        }
        MPI_Finalize();
        exit(9);
    }*/
    //if there is only one mpi thread no splitting to be done.
    if (NProcs==1) {
        StartSnap=0;
        EndSnap=opt.numsnapshots;
        NSnap=opt.numsnapshots;
        return;
    }

    ///determine the total number of haloes and ideal splitting
    int startsnap[NProcs],endsnap[NProcs];
    if (ThisTask==0) {
        //there are two ways to load balance the data, halos or particles in haloes 
        //load balancing via particles is the defaul option but could compile with halo load balancing
        Int_t *numinfo=new Int_t[opt.numsnapshots];
#ifdef MPIHALOBALANCE
        cout<<"Halo based splitting"<<endl;
        Int_t totpart=ReadNumberofHalos(opt,numinfo);
#else
        cout<<"Particle in halos based splitting"<<endl;
        Int_t totpart=ReadNumberofParticlesInHalos(opt,numinfo);
#endif
        //now number of particles per mpi thread
        if (opt.numpermpi==0) opt.numpermpi=totpart/NProcs;
        Int_t sum=0, itask=NProcs-1, inumsteps=-1;
        endsnap[itask]=opt.numsnapshots;
        cout<<" Total number of items is "<<totpart<<" and should have "<<opt.numpermpi<<" per mpi"<<endl;
        for (int i=opt.numsnapshots-1;i>=0;i--) {
            sum+=numinfo[i];
            inumsteps++;
            if (sum>opt.numpermpi && inumsteps>=opt.numsteps) {
                startsnap[itask]=i;
                endsnap[itask-1]=i+opt.numsteps;
                itask--;
                //inumsteps=-1;
                //sum=0;
                inumsteps=opt.numsteps-1;
                sum=0;for (int j=0;j<=opt.numsteps;j++) sum+=numinfo[i+j];
                if (itask==0) break;
            }
        }
        startsnap[itask]=0;
        int itaskworkload=0;
        for (itask=0;itask<NProcs;itask++) 
        {
            sum=0;
            for (int i=startsnap[itask];i<endsnap[itask];i++) sum+=numinfo[i];
            cout<<itask<<" will use snaps "<<startsnap[itask]<<" to "<<endsnap[itask]<<" with overlap of "<<opt.numsteps<<" that contains "<<sum<<" item"<<endl;
            itaskworkload+=(sum>0);
            if (sum==0) {
                cerr<<"WARNING : task "<<itask<<" will not analyze any snapshots with current number of MPI threads"<<endl;
            }
        }
        delete[] numinfo;
        if (itaskworkload<NProcs) {
            cerr<<"MPI theads number + number of desired steps less than optimal, only "<<itaskworkload<<" of "<<NProcs<<" running "<<endl;
            if (itaskworkload<0.75*NProcs) {cerr<<"Exiting, adjust to "<<itaskworkload<<" number of mpi theads "<<endl;MPI_Abort(MPI_COMM_WORLD,1);}
        }
    }
    //task 0 has finished determining which snapshots a given mpi thread will process

    MPI_Bcast(startsnap, NProcs, MPI_Int_t, 0, MPI_COMM_WORLD);
    MPI_Bcast(endsnap, NProcs, MPI_Int_t, 0, MPI_COMM_WORLD);
    StartSnap=startsnap[ThisTask];
    EndSnap=endsnap[ThisTask];
    NSnap=EndSnap-StartSnap;
}

#endif
