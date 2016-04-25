/*! \file main.cxx
 *  \brief Main program

*/

#include "halomergertree.h"

using namespace std;
using namespace Math;
using namespace NBody;

int main(int argc,char **argv)
{
#ifdef USEMPI
    //start MPI
    MPI_Init(&argc,&argv);
    //find out how big the SPMD world is
    MPI_Comm_size(MPI_COMM_WORLD,&NProcs);
    //and this processes' rank is
    MPI_Comm_rank(MPI_COMM_WORLD,&ThisTask);
#else
    int ThisTask=0,NProcs=1,Nlocal,Ntotal;
#endif

    HaloTreeData *pht;
    ProgenitorData **pprogen;
    DescendantData **pdescen;
    ProgenitorData *pprogentemp;
    DescendantData *pdescentemp;
    long unsigned *pfofp,*pfofd;
    long unsigned *noffset;
    long unsigned *pglist;
    //long unsigned **pglist1;
    long int i,j;
    long unsigned nh,nhp,nhd;
    double time1,time2;
    Options opt;
    char fname[1000];

    GetArgs(argc, argv, opt);
    if(opt.matchtype==NsharedN1N2)
        opt.description=(char*)"VELOCIraptor halo merger tree constructed by identifying the main progenitor with the highest value of Nshared^2/Nh/Np";
    else if(opt.matchtype==NsharedN1)
        opt.description=(char*)"VELOCIraptor halo merger tree constructed by identifying the main progenitor with the highest value of Nshared/Nh";
    else if(opt.matchtype==Nshared)
        opt.description=(char*)"VELOCIraptor halo merger tree constructed by identifying the main progenitor with the highest value of Nshared";
    else if (opt.matchtype==Nsharedcombo)
        opt.description=(char*)"VELOCIraptor halo merger tree constructed by identifying the main progenitor with the highest value of Nshared/Nh+(Nshared^2/Nh/Np) so as to weight progenitors that contribute similar amounts by how much of their mass contributes to the new object";

#ifdef USEMPI
    //now if there are too many tasks for the number of snapshots read, exit
    if (NProcs>(0.5*opt.numsnapshots/opt.numsteps)) {
        if (ThisTask==0) {
            cerr<<"Too many mpi threads for number of snapshots and number of linking steps desired. "<<endl;
            cerr<<"Suggested number is less than "<<floor(0.5*opt.numsnapshots/opt.numsteps)<<endl;
            cerr<<"Exiting"<<endl;
        }
        MPI_Finalize();
        exit(9);
    }
    if (NProcs>1) {
    ///determine the total number of haloes and ideal splitting
    int startsnap[NProcs],endsnap[NProcs];
    if (ThisTask==0) {
    Int_t *halonum=new Int_t[opt.numsnapshots];
    Int_t totnum=ReadNumberofHalos(opt,halonum);
    //MPI_Bcast(totnum, 1, MPI_Int_t, 0, MPI_COMM_WORLD);
    //MPI_Bcast(halonum, opt.numsnapshots, MPI_Int_t, 0, MPI_COMM_WORLD);
    //now number of halos per mpi thread
    Int_t halonumpermpi=totnum/NProcs;
    Int_t sum=0, itask=NProcs-1, inumsteps=0;
    endsnap[itask]=opt.numsnapshots;
    cout<<ThisTask<<" total is "<<totnum<<" and should have "<<halonumpermpi<<endl;
    for (int i=opt.numsnapshots-1;i>=0;i--) {
        sum+=halonum[i];
        if (sum>halonumpermpi && inumsteps>opt.numsteps) {
            startsnap[itask]=i;
            endsnap[itask-1]=i+opt.numsteps;
            itask--;
            inumsteps=0;
            sum=0;
            if (itask==0) break;
        }
        inumsteps++;
    }
    startsnap[itask]=0;
    for (itask=0;itask<NProcs;itask++) 
    {
        sum=0;
        for (int i=startsnap[itask];i<endsnap[itask];i++) sum+=halonum[i];
        cout<<itask<<" will use snaps "<<startsnap[itask]<<" to "<<endsnap[itask]<<" with overlap of "<<opt.numsteps<<" that contains "<<sum<<" halos"<<endl;
    }
    delete[] halonum;
    }
    MPI_Bcast(startsnap, NProcs, MPI_Int_t, 0, MPI_COMM_WORLD);
    MPI_Bcast(endsnap, NProcs, MPI_Int_t, 0, MPI_COMM_WORLD);
    StartSnap=startsnap[ThisTask];
    EndSnap=endsnap[ThisTask];
    NSnap=EndSnap-StartSnap;
    
    MPI_Barrier(MPI_COMM_WORLD);
    }
    else {
        StartSnap=0;
        EndSnap=opt.numsnapshots;
        NSnap=opt.numsnapshots;
	}
#endif 
    
    //read particle information and allocate memory
    if (ThisTask==0) cout<<"Read Data ... "<<endl;
    pht=ReadData(opt);
    if (ThisTask==0) {
        cout<<"Done Loading"<<endl;
        cout<<"Found "<<opt.TotalNumberofHalos<<" halos "<<endl;
        cout<<"Memory needed to store addressing for ids "<<sizeof(long unsigned)*opt.NumPart/1024./1024.0/1024.0<<" "<<opt.NumPart<<endl;
    }
    //adjust ids if necessary
    if (opt.imapping>DNOMAP) MapPIDStoIndex(opt,pht);
    //check to see if ID data compatible with accessing index array allocated
    int ierrorflag=0,ierrorsumflag=0;
    for (i=opt.numsnapshots-1;i>=0;i--) {
#ifdef USEMPI
    if (i>=StartSnap && i<EndSnap) {
#endif
        ierrorflag=0;
        for (j=0;j<pht[i].numhalos;j++) {
            for (Int_t k=0;k<pht[i].Halo[j].NumberofParticles;k++) 
                if (pht[i].Halo[j].ParticleID[k]<0||pht[i].Halo[j].ParticleID[k]>opt.NumPart) {
                    cout<<ThisTask<<" for snapshot "<<i<<" particle id out of range "<<pht[i].Halo[j].ParticleID[k]<<" not in [0,"<<opt.NumPart<<")"<<endl;
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

    pfofp=new long unsigned[opt.NumPart];
    for (i=0;i<opt.NumPart;i++) pfofp[i]=0;
    pprogen=new ProgenitorData*[opt.numsnapshots];
    //if particle ids need to be mapped to index

    
    for (i=opt.numsnapshots-1;i>=0;i--) {
#ifdef USEMPI
    //check if data is load making sure i is in appropriate range
    if (i>=StartSnap && i<EndSnap) {
#endif
        //if snapshot contains halos/structures then search for progenitors 
        if(pht[i].numhalos>0){
            cout<<i<<" "<<pht[i].numhalos<<" cross matching objects in standard merger tree direction (progenitors)"<<endl;

            //allocate offset to easily access particle ID/index list
            noffset=new long unsigned[pht[i].numhalos];
            noffset[0]=0;
            nh=pht[i].Halo[0].NumberofParticles;
            for (j=1;j<pht[i].numhalos;j++){nh+=pht[i].Halo[j].NumberofParticles;noffset[j]=noffset[j-1]+pht[i].Halo[j-1].NumberofParticles;}
            pglist=new long unsigned[nh];
            for (j=0;j<pht[i].numhalos;j++){
                for (Int_t k=0;k<pht[i].Halo[j].NumberofParticles;k++)
                    pglist[noffset[j]+k]=pht[i].Halo[j].ParticleID[k];
            }

            //if not last snapshot then can look back in time and produce links
            if (i>0) {
            //loop over all possible steps
            for (Int_t istep=1;istep<=opt.numsteps;istep++) if (i-istep>=0) {
#ifdef USEMPI
            //when using mpi also check that data is local 
            if (i-istep-StartSnap>=0) {
#endif
                //set pfof progenitor data structure, used to produce links. Only produced IF snapshot not first one
                for (j=0;j<pht[i-istep].numhalos;j++) 
                    for (Int_t k=0;k<pht[i-istep].Halo[j].NumberofParticles;k++) 
                        pfofp[pht[i-istep].Halo[j].ParticleID[k]]=j+1;

                //begin cross matching with previous snapshot(s)
                //for first linking, cross match and allocate memory 
                if (istep==1) {
                    pprogen[i]=CrossMatch(opt, opt.NumPart, pht[i].numhalos, pht[i-istep].numhalos, pht[i].Halo, pht[i-istep].Halo, pglist, noffset, pfofp);
                    CleanCrossMatch(pht[i].numhalos, pht[i-istep].numhalos, pht[i].Halo, pht[i-istep].Halo, pprogen[i]);
                }
                //otherwise only care about objects with no links
                else {
                    pprogentemp=CrossMatch(opt, opt.NumPart, pht[i].numhalos, pht[i-istep].numhalos, pht[i].Halo, pht[i-istep].Halo, pglist, noffset, pfofp, istep, pprogen[i]);
                    CleanCrossMatch(pht[i].numhalos, pht[i-istep].numhalos, pht[i].Halo, pht[i-istep].Halo, pprogentemp);
                    UpdateRefProgenitors(opt.imultsteplinkcrit,pht[i].numhalos, pprogen[i], pprogentemp);
                    delete[] pprogentemp;
                }

                //reset pfof2 values 
                for (j=0;j<pht[i-istep].numhalos;j++)
                    for (int k=0;k<pht[i-istep].Halo[j].NumberofParticles;k++) 
                        pfofp[pht[i-istep].Halo[j].ParticleID[k]]=0;
#ifdef USEMPI
            }
#endif
            }
            }
            //otherwise allocate memory but do nothing with it
            else pprogen[i]=new ProgenitorData[pht[i].numhalos];

            //free memory associated with accessing current haloes
            delete[] pglist;
            delete[] noffset;
        }
        //if no haloes/structures then cannot have an progenitors
        else {
            pprogen[i]=NULL;
        }
#ifdef USEMPI
    }
    else pprogen[i]=NULL;
#endif
    }
    delete[] pfofp;
    //Produce a reverse cross comparison or descendant tree if a full graph is requested. 
    if(opt.icatalog==DGRAPH) {
    pdescen=new DescendantData*[opt.numsnapshots];
    pfofd=new long unsigned[opt.NumPart];
    for (i=0;i<opt.NumPart;i++) {pfofd[i]=0;}
    //if particle ids need to be mapped to index
    if (opt.imapping>DNOMAP) MapPIDStoIndex(opt,pht);

    for (i=opt.numsnapshots-1;i>=0;i--) {
#ifdef USEMPI
    if (i>=StartSnap && i<EndSnap) {
#endif
        if(pht[i].numhalos>0){
            cout<<i<<" "<<pht[i].numhalos<<" cross matching objects in other direction (descendants)"<<endl;

            noffset=new long unsigned[pht[i].numhalos];
            noffset[0]=0;
            nh=pht[i].Halo[0].NumberofParticles;
            for (j=1;j<pht[i].numhalos;j++){nh+=pht[i].Halo[j].NumberofParticles;noffset[j]=noffset[j-1]+pht[i].Halo[j-1].NumberofParticles;}
            pglist=new long unsigned[nh];
            for (j=0;j<pht[i].numhalos;j++){
                for (Int_t k=0;k<pht[i].Halo[j].NumberofParticles;k++)
                    pglist[noffset[j]+k]=pht[i].Halo[j].ParticleID[k];
            }


            if (i<opt.numsnapshots-1) {
            for (Int_t istep=1;istep<=opt.numsteps;istep++) if (i+istep<=opt.numsnapshots-1) {
#ifdef USEMPI
            if (i+istep<EndSnap) {
#endif
                //set pfof progenitor data structure, used to produce links. Only produced IF snapshot not first one
                for (j=0;j<pht[i+istep].numhalos;j++)
                    for (int k=0;k<pht[i+istep].Halo[j].NumberofParticles;k++) 
                        pfofd[pht[i+istep].Halo[j].ParticleID[k]]=j+1;

                //begin cross matching with  snapshot(s)
                //for first linking, cross match and allocate memory 
                if (istep==1) {
                    pdescen[i]=CrossMatchDescendant(opt, opt.NumPart, pht[i].numhalos, pht[i+istep].numhalos, pht[i].Halo, pht[i+istep].Halo, pglist, noffset, pfofd);
                    CleanCrossMatchDescendant(pht[i].numhalos, pht[i+istep].numhalos, pht[i].Halo, pht[i+istep].Halo, pdescen[i]);
                }
                //otherwise only care about objects with no links
                else {
                    pdescentemp=CrossMatchDescendant(opt, opt.NumPart, pht[i].numhalos, pht[i+istep].numhalos, pht[i].Halo, pht[i+istep].Halo, pglist, noffset, pfofd,istep, pdescen[i]);
                    CleanCrossMatchDescendant(pht[i].numhalos, pht[i+istep].numhalos, pht[i].Halo, pht[i+istep].Halo, pdescen[i]);
                    UpdateRefDescendants(opt.imultsteplinkcrit,pht[i].numhalos, pdescen[i], pdescentemp);
                    delete[] pdescentemp;
                }

                for (j=0;j<pht[i+istep].numhalos;j++)
                    for (int k=0;k<pht[i+istep].Halo[j].NumberofParticles;k++) 
                        pfofd[pht[i+istep].Halo[j].ParticleID[k]]=0;
#ifdef USEMPI
            }
#endif
            }
            }
            //otherwise allocate memory but do nothing with it
            else pdescen[i]=new DescendantData[pht[i].numhalos];

            delete[] pglist;
            delete[] noffset;
        }
        else {
            pdescen[i]=NULL;
        }
#ifdef USEMPI
    }
    else pdescen[i]=NULL;
#endif
    }
    delete[] pfofd;
    }
    cout<<ThisTask<<" Done"<<endl;

    if (opt.haloidval>0) {
        for (i=opt.numsnapshots-1;i>=0;i--) {
#ifdef USEMPI
        if (i>=StartSnap && i<EndSnap) {
#endif
            if(pht[i].numhalos>0) for (j=0;j<pht[i].numhalos;j++) pht[i].Halo[j].haloID+=opt.haloidval*(i+opt.snapshotvaloffset);
#ifdef USEMPI
        }
#endif
        }
    }
#ifdef USEMPI
    MPI_Barrier(MPI_COMM_WORLD);
#endif
    if (opt.icatalog!=DCROSSCAT) {
        WriteHaloMergerTree(opt,pprogen,pht);
        if (opt.icatalog==DGRAPH) {
            sprintf(fname,"%s.graph",opt.outname);
            sprintf(opt.outname,"%s",fname);
            WriteHaloGraph(opt,pprogen,pdescen,pht);
        }
    }
    else WriteCrossComp(opt,pprogen,pht);

    //free up memory
    for (i=0;i<opt.numsnapshots;i++) if (pprogen[i]!=NULL) delete[] pprogen[i];
    delete[] pprogen;
    //if producing graph as well
    if(opt.icatalog==DGRAPH) {for (i=0;i<opt.numsnapshots;i++) if (pdescen[i]!=NULL) delete[] pdescen[i];delete[] pdescen;}
#ifdef USEMPI
    MPI_Barrier(MPI_COMM_WORLD);
#endif
    delete[] pht;

#ifdef USEMPI
    MPI_Finalize();
#endif
    return 0;
}
