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
set<long long> ConstructMemoryEfficientPIDStoIndexMap(Options &opt, HaloTreeData *&pht) {
    cout<<"Mapping PIDS to index "<<endl;
    Int_t i,j,k;
    Int_t index,indexoffset;
#ifndef USEMPI
    int ThisTask=0,NProcs=1,NSnap=opt.numsnapshots,StartSnap=0,EndSnap=opt.numsnapshots;
#endif
    //we first create vector to store all particle ids of particles in snapshots
    Int_t totnumparts=0,*numparts,*noffset;
    numparts=new Int_t[opt.numsnapshots];
    noffset=new Int_t[opt.numsnapshots];
    noffset[StartSnap]=0;
    for (i=StartSnap;i<EndSnap;i++) {
        numparts[StartSnap]=0;
        for (j=0;j<pht[i].numhalos;j++) numparts[i]+=pht[i].Halo[j].NumberofParticles;
        totnumparts+=numparts[i];
        if (i>StartSnap) noffset[i]=noffset[i-1]+numparts[i-1];
    }
    
    //place ids in a set so have unique ordered set of ids
    set<long long> idset;
    //index=0;
    for (i=0;i<opt.numsnapshots;i++) 
        for (j=0;j<pht[i].numhalos;j++) 
            for (k=0;k<pht[i].Halo[j].NumberofParticles;k++) idset.insert(pht[i].Halo[j].ParticleID[k]);
//            for (k=0;k<pht[i].Halo[j].NumberofParticles;k++) idvec[index++]=pht[i].Halo[j].ParticleID[k];
    //now make unique set
    //set<long long> idset(idvec.begin(), idvec.end());
#ifdef USEMPI

    //now if using mpi then must broadcast the sets that a final unique set can be made
    //this is done each task sending information to its neighbouring task, the lower task 
    //number then building a set of unique ids in the receving task. 
    //once a task has sent information, it nolonger needs to send info
    //all the receiving tasks in the previous round make up the new set of sending and receiving tasks. 
    //keep going till all info cascaded to task 0 and a unique set of ids are generated across all snapshots and mpi domains

    //store uniqe set of ids into vector to send info
    vector<long long> idvec;
    idvec=vector<long long>(idset.begin(), idset.end()); 

    //long long *idarray;
    //int nmpiloops=(int)(floor(exp(log((Double_t)NProcs)-log(2.0))))
    int nrecvtasks,nsendtasks;
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
                //idarray=new long long[commsize];
                //for (j=0;j<commsize;j++) idarray[j]=idvec[j];
                idvec=vector<long long>(idset.begin(), idset.end()); 
                //send size
                MPI_Send(&commsize,1, MPI_Int_t, sendtask[ThisTask], i*NProcs+sendtask[ThisTask], MPI_COMM_WORLD);
                //send info
                MPI_Send(idvec.data(),commsize, MPI_LONG, sendtask[ThisTask], NProcs*NProcs+i*NProcs+sendtask[ThisTask], MPI_COMM_WORLD);
                //delete[] idarray;
            }
            //if a recvtask
            if (recvtask[ThisTask]>=0 && sendtask[ThisTask]==-1) {
                //recv size
                MPI_Recv(&commsize,1, MPI_Int_t, recvtask[ThisTask], i*NProcs+sendtask[ThisTask], MPI_COMM_WORLD,&status);
                //idarray=new long long[commsize];
                idvec.reserve(commsize); 
                //recv info
                MPI_Recv(idvec.data(),commsize, MPI_LONG, recvtask[ThisTask], NProcs*NProcs+i*NProcs+sendtask[ThisTask], MPI_COMM_WORLD,&status);
                //MPI_Recv(idarray,commsize, MPI_LONG, recvtask[ThisTask], NProcs*NProcs+i*NProcs+sendtask[ThisTask], MPI_COMM_WORLD,&status);
                //now generate unique set
                //oldsize=idvec.size();
                //idvec.resize(oldsize+commsize);
                //for (j=0;j<commsize;j++)idvec[j+oldsize]=idarray[j];
                //delete[] idarray;
                //make updated set
                idset.insert(idvec.begin(), idvec.end());
                //store uniqe set of ids into vector 
                //idvec=vector<long long>(idset.begin(), idset.end()); 
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
    } while(numhavesent<NProcs-1);

    //now broadcast from root to all processors
    commsize=idvec.size();
    MPI_Bcast(&commsize,1, MPI_Int_t, 0, MPI_COMM_WORLD);
    idvec.reserve(commsize); 
    MPI_Bcast(idvec.data(),commsize, MPI_LONG, 0, MPI_COMM_WORLD);

    //idarray=new long long[commsize];
    //if (ThisTask==0) for (i=0;i<commsize;i++) idarray[i]=idvec[i];
    //MPI_Bcast(idarray,commsize, MPI_LONG, 0, MPI_COMM_WORLD);
    /*if (ThisTask!=0) {
        idvec.resize(commsize);
        for (i=0;i<commsize;i++) idvec[i]=idarray[i];
    }
    delete idarray;*/
    //from this vector generate the set 
    if (ThisTask!=0) idset=set<long long>(idvec.begin(), idvec.end());
#endif
    opt.MaxIDValue=idset.size();
    return idset;
}

void MapPIDStoIndex(Options &opt, HaloTreeData *&pht, set<long long> &idset) {
    cout<<"Mapping PIDS to index "<<endl;
    Int_t i,j,k;
#ifdef USEOPENMP
#pragma omp parallel default(shared) \
private(i,j,k)
{
#pragma omp for schedule(dynamic) nowait
#endif
    for (i=0;i<opt.numsnapshots;i++) {
            for (j=0;j<pht[i].numhalos;j++) {
                for (k=0;k<pht[i].Halo[j].NumberofParticles;k++){
                    //the set returns an interator that is then used to adjust particle ids
                    pht[i].Halo[j].ParticleID[k]=*idset.find(pht[i].Halo[j].ParticleID[k]);
                }
            }
    }
#ifdef USEOPENMP
}
#endif
}

void MapPIDStoIndex(Options &opt, HaloTreeData *&pht) {
    cout<<"Mapping PIDS to index "<<endl;
    Int_t i,j,k;
#ifdef USEOPENMP
#pragma omp parallel default(shared) \
private(i,j,k)
{
#pragma omp for schedule(dynamic) nowait
#endif
    for (i=0;i<opt.numsnapshots;i++) {
            for (j=0;j<pht[i].numhalos;j++) 
                for (k=0;k<pht[i].Halo[j].NumberofParticles;k++) 
                    opt.mappingfunc(pht[i].Halo[j].ParticleID[k]);
    }
#ifdef USEOPENMP
}
#endif
}

void simplemap(long unsigned int &i) {}

//@}

//check to see if ID data compatible with accessing index array allocated
void IDcheck(Options &opt, HaloTreeData *&pht){
    int ierrorflag=0,ierrorsumflag=0;
    for (int i=opt.numsnapshots-1;i>=0;i--) {
#ifdef USEMPI
    if (i>=StartSnap && i<EndSnap) {
#endif
        ierrorflag=0;
        for (Int_t j=0;j<pht[i].numhalos;j++) {
            for (Int_t k=0;k<pht[i].Halo[j].NumberofParticles;k++) 
                if (pht[i].Halo[j].ParticleID[k]<0||pht[i].Halo[j].ParticleID[k]>opt.MaxIDValue) {
                    cout<<"snapshot "<<i<<" particle id out of range "<<pht[i].Halo[j].ParticleID[k]<<" not in [0,"<<opt.MaxIDValue<<")"<<endl;
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
