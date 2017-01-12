/*! \file main.cxx
 *  \brief Main program

*/

#include "TreeFrog.h"

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
    int StartSnap,EndSnap;
#endif
    int nthreads;
#ifdef USEOPENMP
#pragma omp parallel
    {
    if (omp_get_thread_num()==0) nthreads=omp_get_num_threads();
    }
#else
    nthreads=1;
#endif

#ifdef USEMPI
    if (ThisTask==0) cout<<"VELOCIraptor/Tree running with MPI. Number of mpi threads: "<<NProcs<<endl;
#endif
#ifdef USEOPENMP
    if (ThisTask==0) cout<<"VELOCIraptor/Tree running with OpenMP. Number of openmp threads: "<<nthreads<<endl;
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
    //store the halo ids of particles in objects, with mapping of ids, makes it easy to look up the haloid of a given particle
    unsigned int *pfofp,*pfofd;
    long long i,j;
    long unsigned nh,nhp,nhd;
    //flags indicating whether updates to list are necessary
    int ilistupdated;

    //store the time taken by segments of the code
    double time1,time2;

    GetArgs(argc, argv, opt);
#ifdef USEMPI
    MPILoadBalanceSnapshots(opt);
#else
    StartSnap=0;
    EndSnap=opt.numsnapshots;
#endif
    //read particle information and allocate memory
    if (ThisTask==0) cout<<"Read Data ... "<<endl;
    pht=ReadData(opt);
    if (ThisTask==0) {
        cout<<"Done Loading"<<endl;
        cout<<"Found "<<opt.TotalNumberofHalos<<" halos "<<endl;
        cout<<"Memory needed to store addressing for ids "<<sizeof(long unsigned)*opt.MaxIDValue/1024./1024.0/1024.0<<", maximum ID of "<<opt.MaxIDValue<<endl;
    }
    //adjust ids if particle ids need to be mapped to index
    if (opt.imapping>DNOMAP) MapPIDStoIndex(opt,pht);
    //check that ids are within allowed range for linking
    IDcheck(opt,pht);
    //then allocate simple array used for accessing halo ids of particles through their IDs
    pfofp=new unsigned int[opt.MaxIDValue];
    for (i=0;i<opt.MaxIDValue;i++) pfofp[i]=0;
    //allocate memory associated with progenitors
    pprogen=new ProgenitorData*[opt.numsnapshots];
    //if more than a single snapshot is used to identify possible progenitors then must store the descendants information
    //so the list can later on be cleaned
    if (opt.numsteps==1) pprogendescen=NULL;
    else {
        pprogendescen=new DescendantDataProgenBased*[opt.numsnapshots];
        //initialize all to null
        for (i=opt.numsnapshots-1;i>=0;i--) {
            pprogendescen[i]=NULL;
        }
        //save the last step
        pprogendescen[EndSnap-1]=new DescendantDataProgenBased[pht[EndSnap-1].numhalos];
    }

    time1=MyGetTime();
    if (opt.iverbose) cout<<" Starting the cross matching "<<endl;
    //beginning loop over snapshots and building tree
    for (i=opt.numsnapshots-1;i>=0;i--) {
    //check if data is load making sure i is in appropriate range
    if (i>=StartSnap && i<EndSnap) {
        time2=MyGetTime();
        //if snapshot contains halos/structures then search for progenitors 
        if(pht[i].numhalos>0){
            cout<<i<<" "<<pht[i].numhalos<<" cross matching objects in standard merger tree direction (progenitors)"<<endl;
            //if need to store DescendantDataProgenBased as multiple snapshots, allocate any unallocated pprogendescen needed to process
            //current step if this has not already been allocated
            if (opt.numsteps>1) {
                for (j=1;j<=opt.numsteps;j++) if (i-j>=StartSnap) 
                    if (pprogendescen[i-j]==NULL && pht[i-j].numhalos>0) pprogendescen[i-j]=new DescendantDataProgenBased[pht[i-j].numhalos];
            }
            //if not last snapshot then can look back in time and produce links
            if (i>StartSnap) {
            //loop over all possible steps
            for (Int_t istep=1;istep<=opt.numsteps;istep++) if (i-istep>=0) {
            //when using mpi also check that data is local 
            if (i-istep-StartSnap>=0) {
                //set pfof progenitor data structure, used to produce links. Only produced IF snapshot not first one
                for (j=0;j<pht[i-istep].numhalos;j++) 
                    for (Int_t k=0;k<pht[i-istep].Halo[j].NumberofParticles;k++) 
                        pfofp[pht[i-istep].Halo[j].ParticleID[k]]=j+1;

                //begin cross matching with previous snapshot(s)
                //for first linking, cross match and allocate memory 
                if (istep==1) {
                    pprogen[i]=CrossMatch(opt, pht[i].numhalos, pht[i-istep].numhalos, pht[i].Halo, pht[i-istep].Halo, pfofp, ilistupdated);
                    CleanCrossMatch(pht[i].numhalos, pht[i-istep].numhalos, pht[i].Halo, pht[i-istep].Halo, pprogen[i]);
                    if (opt.numsteps>1) BuildProgenitorBasedDescendantList(i, i-istep, pht[i].numhalos, pprogen[i], pprogendescen);
                }
                //otherwise only care about objects with no links
                else {
                    pprogentemp=CrossMatch(opt, pht[i].numhalos, pht[i-istep].numhalos, pht[i].Halo, pht[i-istep].Halo, pfofp, ilistupdated, istep, pprogen[i]);
                    //if new candidate progenitors have been found then need to update the reference list
                    if (ilistupdated>0) {
                        //if merit limit was applied when deciding to search earlier snapshots
                        CleanCrossMatch(pht[i].numhalos, pht[i-istep].numhalos, pht[i].Halo, pht[i-istep].Halo, pprogentemp);
                        UpdateRefProgenitors(opt,pht[i].numhalos, pprogen[i], pprogentemp, pprogendescen,i);
                        BuildProgenitorBasedDescendantList(i, i-istep, pht[i].numhalos, pprogen[i], pprogendescen, istep);
                    }
                    delete[] pprogentemp;

                }

                //reset pfof2 values 
                for (j=0;j<pht[i-istep].numhalos;j++)
                    for (int k=0;k<pht[i-istep].Halo[j].NumberofParticles;k++) 
                        pfofp[pht[i-istep].Halo[j].ParticleID[k]]=0;
            }
            }
            }
            //otherwise allocate memory but do nothing with it
            else pprogen[i]=new ProgenitorData[pht[i].numhalos];

            //free memory associated with accessing current haloes
            //delete[] pglist;
            //delete[] noffset;
        }
        //if no haloes/structures then cannot have an progenitors
        else {
            pprogen[i]=NULL;
        }
        //now if using multiple snapshots, 
        //if enough non-overlapping (mpi wise) snapshots have been processed, one can cleanup progenitor list using the DescendantDataProgenBased data
        //then free this data
        //this occurs if current snapshot is at least Endsnap-opt.numsteps*2 or lower as then Endsnap-opt.numsteps have had progenitor list processed
        if (opt.numsteps>1 && pht[i].numhalos>0 && (i<EndSnap-2*opt.numsteps && i>StartSnap+2*opt.numsteps)) {
            if (opt.iverbose) cout<<"Cleaning Progenitor list using descendant information for "<<i<<endl;
            CleanProgenitorsUsingDescendants(i, pht, pprogendescen, pprogen);
            delete[] pprogendescen[i];
            pprogendescen[i]=NULL;
        }
        if (opt.iverbose) cout<<ThisTask<<" finished Progenitor processing for snapshot "<<i<<" in "<<MyGetTime()-time2<<endl;
    }
    else pprogen[i]=NULL;
    }
    delete[] pfofp;


    //now if more than one snapshot is used to generate links, must clean up across multiple snapshots to ensure that an object is the progenitor of a single other object
    //this requires parsing the descendent data list produced by identifying progenitors
    if (opt.numsteps>1) {
        //if mpi is used then must broadcast this descendant based progenitor data for snapshots that overlap so that the halos can be appropriately cleaned
#ifdef USEMPI
        if (NProcs>1) MPIUpdateProgenitorsUsingDescendants(opt, pht, pprogendescen, pprogen);
#endif
        for (i=opt.numsnapshots-1;i>=0;i--) {
        //check if data is load making sure i is in appropriate range (note that only look above StartSnap (as first snap can't have progenitors)
        //not that last snapshot locally has halos with no descendants so no need to clean this list
        if (i>=StartSnap && i<EndSnap-1) {
            if (opt.iverbose) cout<<"Cleaning Progenitor list using descendant information for "<<i<<endl;
            if (pprogendescen[i]!=NULL) {
                CleanProgenitorsUsingDescendants(i, pht, pprogendescen, pprogen);
                delete[] pprogendescen[i];
                pprogendescen[i]=NULL;
            }
        }
        }
    }
    if (opt.iverbose) cout<<"Finished the Progenitor cross matching "<<MyGetTime()-time1<<endl;

    //Produce a reverse cross comparison or descendant tree if a full graph is requested. 
    if(opt.icatalog==DGRAPH) {
    time1=MyGetTime();
    if (opt.iverbose) cout<<"Starting descendant cross matching "<<endl;
    pdescen=new DescendantData*[opt.numsnapshots];
    pfofd=new unsigned int[opt.MaxIDValue];
    for (i=0;i<opt.MaxIDValue;i++) {pfofd[i]=0;}
    //if particle ids need to be mapped to index
    if (opt.imapping>DNOMAP) MapPIDStoIndex(opt,pht);

    for (i=opt.numsnapshots-1;i>=0;i--) {
    if (i>=StartSnap && i<EndSnap) {
        time2=MyGetTime();
        if(pht[i].numhalos>0){
            cout<<i<<" "<<pht[i].numhalos<<" cross matching objects in other direction (descendants)"<<endl;

            if (i<=EndSnap) {
            for (Int_t istep=1;istep<=opt.numsteps;istep++) if (i+istep<=opt.numsnapshots-1) {
            if (i+istep<EndSnap) {
                //set pfof progenitor data structure, used to produce links. Only produced IF snapshot not first one
                for (j=0;j<pht[i+istep].numhalos;j++)
                    for (int k=0;k<pht[i+istep].Halo[j].NumberofParticles;k++) 
                        pfofd[pht[i+istep].Halo[j].ParticleID[k]]=j+1;

                //begin cross matching with  snapshot(s)
                //for first linking, cross match and allocate memory 
                if (istep==1) {
                    pdescen[i]=CrossMatchDescendant(opt,  pht[i].numhalos, pht[i+istep].numhalos, pht[i].Halo, pht[i+istep].Halo, pfofd, ilistupdated);
                    CleanCrossMatchDescendant(pht[i].numhalos, pht[i+istep].numhalos, pht[i].Halo, pht[i+istep].Halo, pdescen[i]);
                }
                //otherwise only care about objects with no links
                else {
                    pdescentemp=CrossMatchDescendant(opt, pht[i].numhalos, pht[i+istep].numhalos, pht[i].Halo, pht[i+istep].Halo, pfofd, ilistupdated, istep, pdescen[i]);
                    CleanCrossMatchDescendant(pht[i].numhalos, pht[i+istep].numhalos, pht[i].Halo, pht[i+istep].Halo, pdescen[i]);
                    UpdateRefDescendants(opt,pht[i].numhalos, pdescen[i], pdescentemp);
                    delete[] pdescentemp;
                }

                for (j=0;j<pht[i+istep].numhalos;j++)
                    for (int k=0;k<pht[i+istep].Halo[j].NumberofParticles;k++) 
                        pfofd[pht[i+istep].Halo[j].ParticleID[k]]=0;
            }
            }
            }
            //otherwise allocate memory but do nothing with it
            else pdescen[i]=new DescendantData[pht[i].numhalos];

            //delete[] pglist;
            //delete[] noffset;
        }
        else {
            pdescen[i]=NULL;
        }
        if (opt.iverbose) cout<<ThisTask<<" finished Progenitor processing for snapshot "<<i<<" in "<<MyGetTime()-time2<<endl;
    }
    else pdescen[i]=NULL;
    }
    delete[] pfofd;
    if (opt.iverbose) cout<<"Finished Descendant cross matching "<<MyGetTime()-time1<<endl;
    }

    //now adjust the halo ids as prior, halo ids simple represent read order, allowing for addressing of memory by halo id value
    UpdateHaloIDs(opt,pht);
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
    //if used multiple time steps
    if (opt.numsnapshots>1) {
        for (i=opt.numsnapshots-1;i>=0;i--) {
#ifdef USEMPI
        //check if data is load making sure i is in appropriate range
        if (i>=StartSnap && i<EndSnap) {
#endif
            if (pprogendescen[i]!=NULL) delete[] pprogendescen[i];
#ifdef USEMPI
        }
#endif
        }
    }
#ifdef USEMPI
    delete[] mpi_startsnap;
    delete[] mpi_endsnap;
    MPI_Barrier(MPI_COMM_WORLD);
#endif
    delete[] pht;

#ifdef USEMPI
    MPI_Finalize();
#endif
    return 0;
}
