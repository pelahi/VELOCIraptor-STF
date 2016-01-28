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
    ///\todo need to implement mpi routines. At the moment, just ThisTask zero does anything
    if (ThisTask==0) {

    HaloTreeData *pht;
    ProgenitorData **pprogen;
    DescendantData **pdescen;
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

    //read particle information and allocate memory
    cout<<"Read Data ... "<<endl;
    pht=ReadData(opt);
    cout<<"Done Loading"<<endl;
    cout<<"Memory needed to store addressing for ids "<<sizeof(long unsigned)*opt.NumPart/1024./1024.0/1024.0<<" "<<opt.NumPart<<endl;
    pprogen=new ProgenitorData*[opt.numsnapshots];
    pfofp=new long unsigned[opt.NumPart];

    for (i=0;i<opt.NumPart;i++) pfofp[i]=0;
    //if particle ids need to be mapped to index
    if (opt.imapping>DNOMAP) MapPIDStoIndex(opt,pht);

    for (i=opt.numsnapshots-1;i>=0;i--) {
        if(pht[i].numhalos>0){
            cout<<i<<" "<<pht[i].numhalos<<" cross matching objects in standard merger tree direction (progenitors)"<<endl;

            //need a pfof list, id to index list and pglist
            //set pfof2 values
            if (i>0) {
            for (j=0;j<pht[i-1].numhalos;j++) 
                for (Int_t k=0;k<pht[i-1].Halo[j].NumberofParticles;k++) 
                    pfofp[pht[i-1].Halo[j].ParticleID[k]]=j+1;
            }

            noffset=new long unsigned[pht[i].numhalos];
            noffset[0]=0;
            nh=pht[i].Halo[0].NumberofParticles;
            for (j=1;j<pht[i].numhalos;j++){nh+=pht[i].Halo[j].NumberofParticles;noffset[j]=noffset[j-1]+pht[i].Halo[j-1].NumberofParticles;}
            pglist=new long unsigned[nh];
            for (j=0;j<pht[i].numhalos;j++){
                for (Int_t k=0;k<pht[i].Halo[j].NumberofParticles;k++)
                    pglist[noffset[j]+k]=pht[i].Halo[j].ParticleID[k];
            }
            if (i>0) {
                pprogen[i]=CrossMatch(opt, opt.NumPart, pht[i].numhalos, pht[i-1].numhalos, pht[i].Halo, pht[i-1].Halo, pglist, noffset, pfofp);
                CleanCrossMatch(pht[i].numhalos, pht[i-1].numhalos, pht[i].Halo, pht[i-1].Halo, pprogen[i]);
            }
            else pprogen[i]=new ProgenitorData[pht[i].numhalos];
            if (i>0) {
            for (j=0;j<pht[i-1].numhalos;j++)
                for (int k=0;k<pht[i-1].Halo[j].NumberofParticles;k++) 
                    pfofp[pht[i-1].Halo[j].ParticleID[k]]=0;
            }

            delete[] pglist;
            delete[] noffset;
        }
        else {
            pprogen[i]=NULL;
        }
    }
    delete[] pfofp;
    //if only runing a cross comparison between two different catalogues don't need to determine "descendents"
    if(opt.icatalog==0) {
    pdescen=new DescendantData*[opt.numsnapshots];
    pfofd=new long unsigned[opt.NumPart];
    for (i=0;i<opt.NumPart;i++) {pfofd[i]=0;}
    //if particle ids need to be mapped to index
    if (opt.imapping>DNOMAP) MapPIDStoIndex(opt,pht);

    for (i=opt.numsnapshots-1;i>=0;i--) {
        if(pht[i].numhalos>0){
            cout<<i<<" "<<pht[i].numhalos<<" cross matching objects in other direction (descendants)"<<endl;
            //time1=omp_get_wtime();
            if (i<opt.numsnapshots-1) {
            for (j=0;j<pht[i+1].numhalos;j++)
                for (Int_t k=0;k<pht[i+1].Halo[j].NumberofParticles;k++) 
                    pfofd[pht[i+1].Halo[j].ParticleID[k]]=j+1;
            }

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
                pdescen[i]=CrossMatchDescendant(opt, opt.NumPart, pht[i].numhalos, pht[i+1].numhalos, pht[i].Halo, pht[i+1].Halo, pglist, noffset, pfofd);
                CleanCrossMatchDescendant(pht[i].numhalos, pht[i+1].numhalos, pht[i].Halo, pht[i+1].Halo, pdescen[i]);
            }
            else pdescen[i]=new DescendantData[pht[i].numhalos];

            if (i<opt.numsnapshots-1) {
            for (j=0;j<pht[i+1].numhalos;j++)
                for (int k=0;k<pht[i+1].Halo[j].NumberofParticles;k++) 
                    pfofd[pht[i+1].Halo[j].ParticleID[k]]=0;
            }

            delete[] pglist;
            delete[] noffset;
        }
        else {
            pdescen[i]=NULL;
        }
    }
    delete[] pfofd;
    }
    cout<<"Done"<<endl;
    if (opt.haloidval>0) {
        for (i=opt.numsnapshots-1;i>=0;i--) 
            if(pht[i].numhalos>0) for (j=0;j<pht[i].numhalos;j++) pht[i].Halo[j].haloID+=opt.haloidval*i;
    }
    if (opt.icatalog==0) {
    WriteHaloMergerTree(opt,pprogen,pht);
    sprintf(fname,"%s.graph",opt.outname);
    sprintf(opt.outname,"%s",fname);
    WriteHaloGraph(opt,pprogen,pdescen,pht);
    }   
    else {
    WriteCrossComp(opt,pprogen,pht);
    }
    for (i=0;i<opt.numsnapshots;i++) if (pprogen[i]!=NULL) delete[] pprogen[i];
    delete[] pprogen;
    if(opt.icatalog==0) {for (i=0;i<opt.numsnapshots;i++) if (pdescen[i]!=NULL) delete[] pdescen[i];delete[] pdescen;}
    delete[] pht;
    }

#ifdef USEMPI
    MPI_Finalize();
#endif
    return 0;
}
