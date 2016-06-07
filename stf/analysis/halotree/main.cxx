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

    //store the configuration settings in the options object
    Options opt;

    //store the halo tree
    HaloTreeData *pht;
    //store the progenitors of a given set of objects
    ProgenitorData **pprogen;
    //if multiple snapshots are used in constructing the progenitor list, must store the possible set of descendents 
    //which can then be cleaned to ensure one descendent
    ///\todo eventually must clean up graph construction so that the graph is completely reversible
    DescendantDataProgenBased **pprogendescen;
    
    //stoe the descendants of a given set of objects
    DescendantData **pdescen;
    //these are used to store temporary sets of posssible candidate progenitors/descendants when code links across multiple snaps
    ProgenitorData *pprogentemp;
    DescendantData *pdescentemp;
    ///store the halo ids of particles in objects, with mapping of ids, makes it easy to look up the haloid of a given particle
    long unsigned *pfofp,*pfofd;
    long unsigned *noffset;
    long unsigned *pglist;
    long int i,j;
    long unsigned nh,nhp,nhd;

    //store the time taken by segments of the code
    double time1,time2;

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
    //adjust ids if particle ids need to be mapped to index
    if (opt.imapping>DNOMAP) MapPIDStoIndex(opt,pht);
    //check that ids are within allowed range for linking
    IDcheck(opt,pht);
    //then allocate simple array used for accessing halo ids of particles through their IDs
    pfofp=new long unsigned[opt.NumPart];
    for (i=0;i<opt.NumPart;i++) pfofp[i]=0;
    //allocate memory associated with progenitors
    pprogen=new ProgenitorData*[opt.numsnapshots];
    //if more than a single snapshot is used to identify possible progenitors then must store the descendants information
    //so the list can later on be cleaned
    if (opt.numsnapshots==1) pprogendescen=NULL;
    else {
        pprogendescen=new DescendantDataProgenBased*[opt.numsnapshots];
        for (i=opt.numsnapshots-1;i>=0;i--) {
#ifdef USEMPI
        //check if data is load making sure i is in appropriate range
        if (i>=StartSnap && i<EndSnap) {
#endif
            if (pht[i].numhalos>0) pprogendescen[i]=new DescendantDataProgenBased[pht[i].numhalos];
            else pprogendescen[i]=NULL;
#ifdef USEMPI
        }
#endif
        }
    }
    
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
                    if (opt.numsteps>1) BuildProgenitorBasedDescendantList(i, i-istep, pht[i].numhalos, pprogen[i], pprogendescen);
                }
                //otherwise only care about objects with no links
                else {
                    pprogentemp=CrossMatch(opt, opt.NumPart, pht[i].numhalos, pht[i-istep].numhalos, pht[i].Halo, pht[i-istep].Halo, pglist, noffset, pfofp, istep, pprogen[i]);
                    CleanCrossMatch(pht[i].numhalos, pht[i-istep].numhalos, pht[i].Halo, pht[i-istep].Halo, pprogentemp);
                    UpdateRefProgenitors(opt.imultsteplinkcrit,pht[i].numhalos, pprogen[i], pprogentemp);
                    delete[] pprogentemp;
                    BuildProgenitorBasedDescendantList(i, i-istep, pht[i].numhalos, pprogen[i], pprogendescen, istep);

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
    
    //now if more than one snapshot is used to generate links, must clean up across multiple snapshots to ensure that an object is the progenitor of a single other object
    //this requires parsing the descendent data list produced by identifying progenitors
    if (opt.numsteps>1) {
    for (i=opt.numsnapshots-1;i>=0;i--) {
#ifdef USEMPI
    //check if data is load making sure i is in appropriate range
    if (i>=StartSnap && i<EndSnap) {
#endif
        CleanProgenitorsUsingDescendants(i, pht, pprogendescen, pprogen);
#ifdef USEMPI
    }
#endif
    }
    }

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

    //now adjust the halo ids as prior, halo ids simple represent read order, allowing for addressing of memory by halo id value
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
    
    //now start writing the results
    if (opt.icatalog!=DCROSSCAT) {
        WriteHaloMergerTree(opt,pprogen,pht);
        if (opt.icatalog==DGRAPH) {
            char fname[1000];
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
