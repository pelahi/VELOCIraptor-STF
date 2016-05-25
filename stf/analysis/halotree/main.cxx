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
#ifdef USEMPI
    MPILoadBalanceSnapshots(opt);
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
    //check that ids are within allowed range for linking
    IDcheck(opt,pht);

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
