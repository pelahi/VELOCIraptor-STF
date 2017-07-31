/*! \file idroutines.cxx
 *  \brief this file contains routines for manipulating/updating ids, particles and halos
 */

#include "TreeFrog.h"

/// \name if halo ids need to be adjusted
//@{
void UpdateHaloIDs(Options &opt, HaloTreeData *&pht) {
    Int_t i,j;
    if (opt.haloidval>0) {
        for (i=opt.numsnapshots-1;i>=0;i--) {
#ifdef USEMPI
        if (i>=StartSnap && i<EndSnap) {
#endif
            if(pht[i].numhalos>0) for (j=0;j<pht[i].numhalos;j++) pht[i].Halo[j].haloID+=opt.haloidval*(i+opt.snapshotvaloffset)+opt.haloidoffset;
#ifdef USEMPI
        }
#endif
        }
    }
}
//@}

/// \name if particle ids need to be mapped to indices
//@{
///builds a map by identifying the ids of particles in structure across snapshots
///which is memory efficient as only needs array size of maximum number of particles in structures
map<IDTYPE, IDTYPE> ConstructMemoryEfficientPIDStoIndexMap(Options &opt, HaloTreeData *&pht) {
    Int_t i,j,k;
    Int_t index,indexoffset,offsetsnap;
    double time1,time2;
#ifndef USEMPI
    int ThisTask=0,NProcs=1,NSnap=opt.numsnapshots,StartSnap=0,EndSnap=opt.numsnapshots;
#endif
    //to handle overlapping snapshots between mpi domains
    if (ThisTask==0) offsetsnap=0;
    else offsetsnap=opt.numsteps;

    if (ThisTask==0) cout<<"Mapping PIDS to index "<<endl;
    //place ids in a set so have unique ordered set of ids
    unordered_set<IDTYPE> idset;
    vector<IDTYPE> idvec;
    map<IDTYPE, IDTYPE> idmap;
    //index=0;
    time1=MyGetTime();
    time2=MyGetTime();
    Int_t reservesize=0;
    for (j=0;j<pht[EndSnap-1].numhalos;j++) reservesize+=pht[EndSnap-1].Halo[j].NumberofParticles;
    idset.reserve(reservesize);
    cout<<ThisTask<<" initial reserve of memory for construction of id map is "<<reservesize*sizeof(IDTYPE)/1024./1024./1024.<<"GB"<<endl;
    for (i=StartSnap+offsetsnap;i<EndSnap;i++)
        for (j=0;j<pht[i].numhalos;j++)
            for (k=0;k<pht[i].Halo[j].NumberofParticles;k++) idset.insert(pht[i].Halo[j].ParticleID[k]);
    time1=MyGetTime()-time1;
    if (opt.iverbose) cout<<ThisTask<<" finished getting unique ids of "<<idset.size()<<" in "<<time1<<endl;
#ifdef USEMPI

    //now if using mpi then must broadcast the sets that a final unique set can be made
    //this is done each task sending information to its neighbouring task, the lower task
    //number then building a set of unique ids in the receving task.
    //once a task has sent information, it nolonger needs to send info
    //all the receiving tasks in the previous round make up the new set of sending and receiving tasks.
    //keep going till all info cascaded to task 0 and a unique set of ids are generated across all snapshots and mpi domains

    //used to transmit info
    IDTYPE *idarray;
    //int nmpiloops=(int)(floor(exp(log((Double_t)NProcs)-log(2.0))))
    int nrecvtasks,nsendtasks, numloops=0;
    MPI_Status status;
    int sendtask[NProcs],recvtask[NProcs], numhavesent=0;
    vector<int> sendset,recvset;
    Int_t commsize,oldsize;

    //for each task determine whether sending/receiving and to/from who
    //initialise the send/recv
    for (i=0;i<NProcs;i++) sendtask[i]=recvtask[i]=-1;
    //task zero is explicitly never a send task
    for (i=NProcs-1;i>=1;i-=2) {sendtask[i]=i-1;sendset.push_back(i);}
    for (i=NProcs-2;i>=0;i-=2) {recvtask[i]=i+1;recvset.push_back(i);}
    //if odd number of tasks, then need to explicity add 0 to the recv task set
    if (NProcs%2==1) recvset.push_back(0);

    //keep sending till number of sending tasks == NProcs-1 and task 0 has all the info
    do {
        //if process no longer sending or receiving do nothing in this round
        if (!(sendtask[ThisTask]==-1 && recvtask[ThisTask]==-1)) {
            //if a send task
            if (sendtask[ThisTask]>=0 && recvtask[ThisTask]==-1) {
                commsize=idset.size();
                idarray=new IDTYPE[commsize];
                j=0; for (auto idval: idset) idarray[j++]=idval;
                //send size
                MPI_Send(&commsize,1, MPI_Int_t, sendtask[ThisTask], numloops*NProcs*NProcs+ThisTask, MPI_COMM_WORLD);
                //send info
                //MPI_Send(idset.data,commsize, MPI_Id_type, sendtask[ThisTask], numloops*NProcs*NProcs+NProcs+ThisTask, MPI_COMM_WORLD);
                MPI_Send(idarray,commsize, MPI_Id_type, sendtask[ThisTask], numloops*NProcs*NProcs+NProcs+ThisTask, MPI_COMM_WORLD);
                delete[] idarray;
            }
            //if a recvtask
            if (recvtask[ThisTask]>=0 && sendtask[ThisTask]==-1) {
                //recv size
                MPI_Recv(&commsize,1, MPI_Int_t, recvtask[ThisTask], numloops*NProcs*NProcs+recvtask[ThisTask], MPI_COMM_WORLD,&status);
                idarray=new IDTYPE[commsize];
                //recv info
                MPI_Recv(idarray,commsize, MPI_Id_type, recvtask[ThisTask], numloops*NProcs*NProcs+NProcs+recvtask[ThisTask], MPI_COMM_WORLD,&status);
                for (i=0;i<commsize;i++) idset.insert(idarray[i]);
                delete[] idarray;
            }
        }
        //update the send receive taks for the next time around
        //all tasks that were send tasks now must not be included.
        MPI_Barrier(MPI_COMM_WORLD);
        nsendtasks=sendset.size();
        numhavesent+=nsendtasks;
        nrecvtasks=recvset.size();
        //reset tasks
        for (j=0;j<nsendtasks;j++) sendtask[sendset[j]]=-1;
        for (j=0;j<nrecvtasks;j++) recvtask[recvset[j]]=-1;
        //only tasks that were receiving before are part of the new send/recv set
        for (j=0;j<nrecvtasks-1;j+=2) sendtask[recvset[j]]=recvset[j+1];
        for (j=1;j<nrecvtasks;j+=2) recvtask[recvset[j]]=recvset[j-1];
        sendset.clear();
        recvset.clear();
        for (j=NProcs-1;j>=0;j--) if (sendtask[j]>=0) sendset.push_back(j);
        for (j=NProcs-1;j>=0;j--) if (recvtask[j]>=0) recvset.push_back(j);
        //if current set of communicating tasks set by the old recv task set is odd must explicity include 0
        if (nrecvtasks%2==1)recvset.push_back(0);
        numloops++;
    } while(numhavesent<NProcs-1);

    //now broadcast from root to all processors
    if (ThisTask==0) {
        commsize=idset.size();
        idvec=vector<IDTYPE>(idset.begin(), idset.end());
    }
    MPI_Bcast(&commsize,1, MPI_Int_t, 0, MPI_COMM_WORLD);
    if (ThisTask!=0) idvec=vector<IDTYPE>(commsize);
    MPI_Bcast(idvec.data(),commsize, MPI_Id_type, 0, MPI_COMM_WORLD);
#else
    idvec=vector<IDTYPE>(idset.begin(), idset.end());
#endif
    opt.MaxIDValue=idvec.size();
    for (i=0;i<opt.MaxIDValue;i++) idmap.insert(pair<IDTYPE, IDTYPE>(idvec[i],(IDTYPE)i));
    time2=MyGetTime()-time2;
#ifdef USEMPI
    if (opt.iverbose && ThisTask==0) cout<<ThisTask<<" finished getting GLOBAL unique ids of "<<idvec.size()<<" in "<<time2<<endl;
#endif
    return idmap;
}

void MapPIDStoIndex(Options &opt, HaloTreeData *&pht, map<IDTYPE, IDTYPE> &idmap) {
#ifndef USEMPI
    int ThisTask=0;
#endif
    Int_t i,j,k;
    if (ThisTask==0) cout<<"Mapping PIDS to index "<<endl;
#ifdef USEOPENMP
#pragma omp parallel default(shared) \
private(i,j,k)
{
#pragma omp for schedule(dynamic) nowait
#endif
    for (i=0;i<opt.numsnapshots;i++) {
#ifdef USEMPI
    if (i>=StartSnap && i<EndSnap) {
#endif
        for (j=0;j<pht[i].numhalos;j++) {
            for (k=0;k<pht[i].Halo[j].NumberofParticles;k++){
                pht[i].Halo[j].ParticleID[k]=idmap[pht[i].Halo[j].ParticleID[k]];
            }
        }
#ifdef USEMPI
    }
#endif
    }
#ifdef USEOPENMP
}
#endif
}

void MapPIDStoIndex(Options &opt, HaloTreeData *&pht) {
#ifndef USEMPI
    int ThisTask=0;
#endif
    Int_t i,j,k;
    if (ThisTask==0) cout<<"Mapping PIDS to index "<<endl;
#ifdef USEOPENMP
#pragma omp parallel default(shared) \
private(i,j,k)
{
#pragma omp for schedule(dynamic) nowait
#endif
    for (i=0;i<opt.numsnapshots;i++) {
#ifdef USEMPI
    if (i>=StartSnap && i<EndSnap) {
#endif
        for (j=0;j<pht[i].numhalos;j++)
            for (k=0;k<pht[i].Halo[j].NumberofParticles;k++)
                opt.mappingfunc(pht[i].Halo[j].ParticleID[k]);
#ifdef USEMPI
    }
#endif
    }
#ifdef USEOPENMP
}
#endif
}

void simplemap(IDTYPE &i) {}

//@}

//check to see if ID data compatible with accessing index array allocated
void IDcheck(Options &opt, HaloTreeData *&pht){
    int ierrorflag=0,ierrorsumflag=0;
#ifndef USEMPI
    int ThisTask=0,NProcs=1,NSnap=opt.numsnapshots,StartSnap=0,EndSnap=opt.numsnapshots;
#endif
    for (int i=opt.numsnapshots-1;i>=0;i--) {
#ifdef USEMPI
    if (i>=StartSnap && i<EndSnap) {
#endif
        ierrorflag=0;
        for (Int_t j=0;j<pht[i].numhalos;j++) {
            for (Int_t k=0;k<pht[i].Halo[j].NumberofParticles;k++)
                if (pht[i].Halo[j].ParticleID[k]<0||pht[i].Halo[j].ParticleID[k]>opt.MaxIDValue) {
                    cout<<ThisTask<<" snapshot "<<i<<" particle id out of range "<<pht[i].Halo[j].ParticleID[k]<<" not in [0,"<<opt.MaxIDValue<<")"<<endl;
                    ierrorflag=1;
                    break;
                }
            if (ierrorflag==1) {ierrorsumflag+=1;break;}
        }
#ifdef USEMPI
    }
#endif
    }
#ifdef USEMPI
    int mpierrorflag;
    MPI_Allreduce(&ierrorsumflag,&mpierrorflag,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
    ierrorsumflag=mpierrorflag;
#endif
    if (ierrorsumflag>0) {
#ifdef USEMPI
        if (ThisTask==0)
#endif
        cout<<"Error in particle ids, outside of range. Exiting"<<endl;
#ifdef USEMPI
        MPI_Abort(MPI_COMM_WORLD,9);
#else
        exit(9);
#endif
    }
}
