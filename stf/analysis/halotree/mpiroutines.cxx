/*! \file mpiroutines.cxx
 *  \brief this file contains routines used with MPI compilation. 

    MPI routines generally pertain to domain decomposition or to specific MPI tasks that determine what needs to be broadcast between various threads.

 */

#ifdef USEMPI

//-- For MPI

#include "halomergertree.h"

///load balance how snapshots are split across mpi domains
void MPILoadBalanceSnapshots(Options &opt){
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

///mpi routine to broadcast the progenitor based descendant list for snapshots that overlap between mpi tasks. Ensures that descendent lists is complete and 
///can be used to clean progenitor list
void MPIUpdateProgenitorsUsingDescendants(Options &opt, HaloTreeData *&pht, DescendantDataProgenBased **&pprogendescen, ProgenitorData **&pprogen)
{
    ///broadcast overlaping snapshots
    ///this is done by broadcasting start and end snapshots 
    int mpi_startsnap[NProcs],mpi_endsnap[NProcs];
    int sendtask,recvtask, isnap;
    MPI_Status *status;
    MPI_Allgather(&StartSnap,1,MPI_INT,mpi_startsnap,NProcs,MPI_INT,MPI_COMM_WORLD);
    MPI_Allgather(&EndSnap,1,MPI_INT,mpi_endsnap,NProcs,MPI_INT,MPI_COMM_WORLD);

    ///vector data structures used to store the data that is parsed using an offset pointer
    Int_t totitems;
    vector<unsigned long> noffset;
    vector<Int_t>numdescen;
    vector<long unsigned> tothaloindex;
    vector<long unsigned> tothalotemporalindex;
    vector<Double_t> totmerit;
    vector<int> totMPITask;
    //first send information from task itask to itask+1
    //note that as the DescendantDataProgenBased is an array have to send size then each element which itself contains arrays
    //ideally we could use boost to specify the mpi data format but for now lets just generate vectors containing all the data that can then be parsed
    //first send info up to the next task
    ///\todo need to clean up send style so that processors are waiting to do nothing
    ///plus there are issues with the buffer with this type of send style. Might need to be more clever
    for (int itask=0;itask<NProcs-1;itask++) {
    for (int i=1;i<opt.numsteps;i++) {
        //sending information
        if (itask==ThisTask) {
            recvtask=ThisTask+1;
            isnap=EndSnap-i;
            //build data structures to be sent
            totitems=0;
            numdescen.resize(pht[isnap].numhalos);
            noffset.resize(pht[isnap].numhalos);
            for (Int_t j=0;j<pht[isnap].numhalos;j++) {
                numdescen[j]=pprogendescen[isnap][j].NumberofDescendants;
                noffset[j]=totitems;
                totitems+=pprogendescen[isnap][j].NumberofDescendants;
            }
            tothaloindex.resize(totitems);
            tothalotemporalindex.resize(totitems);
            totmerit.resize(totitems);
            totMPITask.resize(totitems);
            totitems=0;
            for (Int_t j=0;j<pht[isnap].numhalos;j++) {
                for (Int_t k=0;i<numdescen[j];k++) {
                    tothaloindex[noffset[j]+k]=pprogendescen[isnap][j].haloindex[k];
                    tothalotemporalindex[noffset[j]+k]=pprogendescen[isnap][j].halotemporalindex[k];
                    totmerit[noffset[j]+k]=pprogendescen[isnap][j].Merit[k];
                    totMPITask[noffset[j]+k]=pprogendescen[isnap][j].MPITask[k];
                }
                totitems+=pprogendescen[isnap][j].NumberofDescendants;
            }
            //then send number of halos, total number of items and then the offset data so that the combined data can be parsed
            MPI_Send(&totitems,1, MPI_Int_t, recvtask, ThisTask+NProcs, MPI_COMM_WORLD);
            MPI_Send(&noffset[0],pht[isnap].numhalos, MPI_LONG, recvtask, ThisTask+2*NProcs, MPI_COMM_WORLD);
            MPI_Send(&numdescen[0],pht[isnap].numhalos, MPI_INT, recvtask, ThisTask+3*NProcs, MPI_COMM_WORLD);
            MPI_Send(&tothaloindex[0],totitems, MPI_LONG, recvtask, ThisTask+4*NProcs, MPI_COMM_WORLD);
            MPI_Send(&tothalotemporalindex[0],totitems, MPI_LONG, recvtask, ThisTask+5*NProcs, MPI_COMM_WORLD);
            MPI_Send(&totmerit[0],totitems, MPI_LONG, recvtask, ThisTask+6*NProcs, MPI_COMM_WORLD);
            MPI_Send(&totMPITask[0],totitems, MPI_LONG, recvtask, ThisTask+7*NProcs, MPI_COMM_WORLD);
        }
        //receiving information
        else if (itask==ThisTask-1) {
            isnap=mpi_endsnap[itask]-i;
            //store the number of items, allocate memory and then receive
            MPI_Recv(&totitems,1, MPI_Int_t, itask, itask+NProcs, MPI_COMM_WORLD,status);
            numdescen.resize(pht[isnap].numhalos);
            noffset.resize(pht[isnap].numhalos);
            tothaloindex.resize(totitems);
            tothalotemporalindex.resize(totitems);
            totmerit.resize(totitems);
            totMPITask.resize(totitems);
            MPI_Recv(&noffset[0],pht[isnap].numhalos, MPI_LONG, itask, itask+2*NProcs, MPI_COMM_WORLD,status);
            MPI_Recv(&numdescen[0],pht[isnap].numhalos, MPI_INT, itask, itask+3*NProcs, MPI_COMM_WORLD,status);
            MPI_Recv(&tothaloindex[0],totitems, MPI_LONG, itask, itask+4*NProcs, MPI_COMM_WORLD,status);
            MPI_Recv(&tothalotemporalindex[0],totitems, MPI_LONG, itask, itask+5*NProcs, MPI_COMM_WORLD,status);
            MPI_Recv(&totmerit[0],totitems, MPI_LONG, itask, itask+6*NProcs, MPI_COMM_WORLD,status);
            MPI_Recv(&totMPITask[0],totitems, MPI_LONG, itask, itask+7*NProcs, MPI_COMM_WORLD,status);

            //then merge the mpi data together
            for (Int_t j=0;j<pht[isnap].numhalos;j++)
                pprogendescen[isnap][j].Merge(numdescen[j],&tothaloindex[noffset[j]],&tothalotemporalindex[noffset[j]],&totmerit[noffset[j]],&totMPITask[noffset[j]]);
        }
    }
    MPI_Barrier(MPI_COMM_WORLD);
    }
    
    //now repeat but going the other way
    for (int itask=NProcs-1;itask>0;itask++) {
    for (int i=0;i<opt.numsteps-1;i++) {
        //sending information
        if (itask==ThisTask) {
            recvtask=ThisTask-1;
            //build data structures to be sent
            totitems=0;
            isnap=StartSnap+i;
            numdescen.resize(pht[isnap].numhalos);
            noffset.resize(pht[isnap].numhalos);
            for (Int_t j=0;j<pht[isnap].numhalos;j++) {
                numdescen[j]=pprogendescen[isnap][j].NumberofDescendants;
                noffset[j]=totitems;
                totitems+=pprogendescen[isnap][j].NumberofDescendants;
            }
            tothaloindex.resize(totitems);
            tothalotemporalindex.resize(totitems);
            totmerit.resize(totitems);
            totMPITask.resize(totitems);
            totitems=0;
            for (Int_t j=0;j<pht[isnap].numhalos;j++) {
                for (Int_t k=0;i<numdescen[j];k++) {
                    tothaloindex[noffset[j]+k]=pprogendescen[isnap][j].haloindex[k];
                    tothalotemporalindex[noffset[j]+k]=pprogendescen[isnap][j].halotemporalindex[k];
                    totmerit[noffset[j]+k]=pprogendescen[isnap][j].Merit[k];
                    totMPITask[noffset[j]+k]=pprogendescen[isnap][j].MPITask[k];
                }
                totitems+=pprogendescen[isnap][j].NumberofDescendants;
            }
            //then send number of halos, total number of items and then the offset data so that the combined data can be parsed
            MPI_Send(&totitems,1, MPI_Int_t, recvtask, ThisTask+NProcs, MPI_COMM_WORLD);
            MPI_Send(&noffset[0],pht[isnap].numhalos, MPI_LONG, recvtask, ThisTask+2*NProcs, MPI_COMM_WORLD);
            MPI_Send(&numdescen[0],pht[isnap].numhalos, MPI_INT, recvtask, ThisTask+3*NProcs, MPI_COMM_WORLD);
            MPI_Send(&tothaloindex[0],totitems, MPI_LONG, recvtask, ThisTask+4*NProcs, MPI_COMM_WORLD);
            MPI_Send(&tothalotemporalindex[0],totitems, MPI_LONG, recvtask, ThisTask+5*NProcs, MPI_COMM_WORLD);
            MPI_Send(&totmerit[0],totitems, MPI_LONG, recvtask, ThisTask+6*NProcs, MPI_COMM_WORLD);
            MPI_Send(&totMPITask[0],totitems, MPI_LONG, recvtask, ThisTask+7*NProcs, MPI_COMM_WORLD);
        }
        //receiving information
        else if (itask==ThisTask+1) {
            isnap=mpi_startsnap[itask]+i;
            //store the number of items, allocate memory and then receive
            MPI_Recv(&totitems,1, MPI_Int_t, itask, itask+NProcs, MPI_COMM_WORLD,status);
            numdescen.resize(pht[isnap].numhalos);
            noffset.resize(pht[isnap].numhalos);
            tothaloindex.resize(totitems);
            tothalotemporalindex.resize(totitems);
            totmerit.resize(totitems);
            totMPITask.resize(totitems);
            MPI_Recv(&noffset[0],pht[isnap].numhalos, MPI_LONG, itask, itask+2*NProcs, MPI_COMM_WORLD,status);
            MPI_Recv(&numdescen[0],pht[isnap].numhalos, MPI_INT, itask, itask+3*NProcs, MPI_COMM_WORLD,status);
            MPI_Recv(&tothaloindex[0],totitems, MPI_LONG, itask, itask+4*NProcs, MPI_COMM_WORLD,status);
            MPI_Recv(&tothalotemporalindex[0],totitems, MPI_LONG, itask, itask+5*NProcs, MPI_COMM_WORLD,status);
            MPI_Recv(&totmerit[0],totitems, MPI_LONG, itask, itask+6*NProcs, MPI_COMM_WORLD,status);
            MPI_Recv(&totMPITask[0],totitems, MPI_LONG, itask, itask+7*NProcs, MPI_COMM_WORLD,status);

            //then merge the mpi data together
            for (Int_t j=0;j<pht[isnap].numhalos;j++)
                pprogendescen[isnap][j].Merge(numdescen[j],&tothaloindex[noffset[j]],&tothalotemporalindex[noffset[j]],&totmerit[noffset[j]],&totMPITask[noffset[j]]);
        }
    }
    MPI_Barrier(MPI_COMM_WORLD);
    }

}

#endif
