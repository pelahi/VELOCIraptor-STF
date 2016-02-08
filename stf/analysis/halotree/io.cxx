/*! \file io.cxx
 *  \brief this file contains routines for io
 */

#include "halomergertree.h"

///Read data from a number of snapshot files
HaloTreeData *ReadData(Options &opt)
{
    HaloTreeData *HaloTree;
    fstream Fin;//file is list of halo data files
    long unsigned j,nparts,haloid;
    HaloTree=new HaloTreeData[opt.numsnapshots];
    char buf[1000*opt.numsnapshots];
    Int_t tothalos=0;
    int i,nthreads;
    Fin.open(opt.fname);
    if (!Fin.is_open()) {cerr<<"file containing snapshot list can't be opened"<<endl;exit(9);}
#ifndef USEOPENMP
    nthreads=1;
#else
#pragma omp parallel 
    {
            if (omp_get_thread_num()==0) nthreads=omp_get_num_threads();
    }
#endif

#ifdef USEOPENMP
  for(i=0; i<opt.numsnapshots; i++) Fin>>&(buf[i*1000]);
#pragma omp parallel default(shared) \
private(i)
{
#pragma omp for reduction(+:tothalos)
  for(i=0; i<opt.numsnapshots; i++)
  {
    if (opt.ioformat==DSUSSING) HaloTree[i].Halo=ReadHaloData(&buf[i*1000],HaloTree[i].numhalos);
    else if (opt.ioformat==DCATALOG) HaloTree[i].Halo=ReadHaloGroupCatalogData(&buf[i*1000],HaloTree[i].numhalos, opt.nmpifiles, opt.ibinary,opt.ifield, opt.itypematch);
    else if (opt.ioformat==DNIFTY) HaloTree[i].Halo=ReadNIFTYData(&buf[i*1000],HaloTree[i].numhalos, opt.idcorrectflag);
    tothalos+=HaloTree[i].numhalos;
  }
}
  opt.TotalNumberofHalos=tothalos;
#else
  for(i=0; i<opt.numsnapshots; i++)
  {
    Fin>>buf;
    if (opt.ioformat==DSUSSING) HaloTree[i].Halo=ReadHaloData(buf,HaloTree[i].numhalos);
    else if (opt.ioformat==DCATALOG) HaloTree[i].Halo=ReadHaloGroupCatalogData(buf,HaloTree[i].numhalos, opt.nmpifiles, opt.ibinary, opt.ifield, opt.itypematch);
    else if (opt.ioformat==DNIFTY) HaloTree[i].Halo=ReadNIFTYData(buf,HaloTree[i].numhalos, opt.idcorrectflag);
    opt.TotalNumberofHalos+=HaloTree[i].numhalos;
  }
#endif
  Fin.close();
  return HaloTree;
}


void WriteHaloMergerTree(Options &opt, ProgenitorData **p, HaloTreeData *h) {
    fstream Fout;
    char fname[1000];
    int istep;
    sprintf(fname,"%s",opt.outname);
    cout<<"saving halo merger tree to "<<fname<<endl;
    Fout.open(fname,ios::out);
    Fout<<1<<endl;
    Fout<<opt.description<<endl;
    Fout<<opt.TotalNumberofHalos<<endl;
    if (opt.outputformat==0) {
    for (int i=opt.numsnapshots-1;i>0;i--) {
        Fout<<i+opt.snapshotvaloffset<<"\t"<<h[i].numhalos<<endl;
        for (int j=0;j<h[i].numhalos;j++) {
            Fout<<h[i].Halo[j].haloID<<"\t"<<p[i][j].NumberofProgenitors<<endl;
            for (int k=0;k<p[i][j].NumberofProgenitors;k++) {
                Fout<<h[i-p[i][j].istep].Halo[p[i][j].ProgenitorList[k]-1].haloID<<endl;
                //Fout<<h[i-1].Halo[p[i][j].ProgenitorList[k]-1].haloID<<endl;
            }
        }
    }
    }
    else {
    for (int i=opt.numsnapshots-1;i>0;i--) {
        Fout<<i+opt.snapshotvaloffset<<"\t"<<h[i].numhalos<<endl;
        for (int j=0;j<h[i].numhalos;j++) {
            Fout<<h[i].Halo[j].haloID<<"\t"<<p[i][j].NumberofProgenitors<<endl;
            for (int k=0;k<p[i][j].NumberofProgenitors;k++) {
                Fout<<h[i-p[i][j].istep].Halo[p[i][j].ProgenitorList[k]-1].haloID<<" "<<p[i][j].Merit[k]<<endl;
                //Fout<<h[i-1].Halo[p[i][j].ProgenitorList[k]-1].haloID<<" "<<p[i][j].Merit[k]<<endl;
            }
        }
    }
    }
    ///last file has no connections
    for (int j=0;j<h[0].numhalos;j++) {
        Fout<<h[0].Halo[j].haloID<<"\t"<<0<<endl;
    }
    Fout<<"END"<<endl;
    Fout.close();
}

void WriteHaloGraph(Options &opt, ProgenitorData **p, DescendantData **d, HaloTreeData *h) {
    fstream Fout;
    char fname[1000];
    sprintf(fname,"%s",opt.outname);
    cout<<"saving halo merger tree to "<<fname<<endl;
    Fout.open(fname,ios::out);
    Fout<<1<<endl;
    Fout<<opt.description<<endl;
    Fout<<opt.TotalNumberofHalos<<endl;
    if (opt.outputformat==0) {
    for (int i=opt.numsnapshots-1;i>=0;i--) {
        Fout<<i+opt.snapshotvaloffset<<"\t"<<h[i].numhalos<<endl;
        for (int j=0;j<h[i].numhalos;j++) {
            Fout<<h[i].Halo[j].haloID<<"\t"<<p[i][j].NumberofProgenitors<<"\t"<<d[i][j].NumberofDescendants<<endl;
            for (int k=0;k<p[i][j].NumberofProgenitors;k++) {
                Fout<<h[i-p[i][j].istep].Halo[p[i][j].ProgenitorList[k]-1].haloID<<endl;
            }
            for (int k=0;k<d[i][j].NumberofDescendants;k++) {
                Fout<<h[i+d[i][j].istep].Halo[d[i][j].DescendantList[k]-1].haloID<<endl;
            }
        }
    }
    }
    else {
    for (int i=opt.numsnapshots-1;i>=0;i--) {
        Fout<<i+opt.snapshotvaloffset<<"\t"<<h[i].numhalos<<endl;
        for (int j=0;j<h[i].numhalos;j++) {
            Fout<<h[i].Halo[j].haloID<<"\t"<<p[i][j].NumberofProgenitors<<"\t"<<d[i][j].NumberofDescendants<<endl;
            for (int k=0;k<p[i][j].NumberofProgenitors;k++) {
                Fout<<h[i-p[i][j].istep].Halo[p[i][j].ProgenitorList[k]-1].haloID<<" "<<p[i][j].Merit[k]<<endl;
            }
            for (int k=0;k<d[i][j].NumberofDescendants;k++) {
                Fout<<h[i+d[i][j].istep].Halo[d[i][j].DescendantList[k]-1].haloID<<" "<<d[i][j].Merit[k]<<endl;
            }
        }
    }
    }
    Fout<<"END"<<endl;
    Fout.close();
}

void WriteCrossComp(Options &opt, ProgenitorData **p, HaloTreeData *h) {
    fstream Fout;
    char fname[1000];
    sprintf(fname,"%s",opt.outname);
    cout<<"saving halo merger tree to "<<fname<<endl;
    Fout.open(fname,ios::out);
    Fout<<opt.description<<endl;
    Fout<<opt.TotalNumberofHalos<<endl;
    if (opt.outputformat==0) {
    for (int i=opt.numsnapshots-1;i>0;i--) {
        Fout<<i<<"\t"<<h[i].numhalos<<endl;
        for (int j=0;j<h[i].numhalos;j++) {
            Fout<<h[i].Halo[j].haloID<<"\t"<<p[i][j].NumberofProgenitors<<endl;
            for (int k=0;k<p[i][j].NumberofProgenitors;k++) {
                Fout<<h[i-1].Halo[p[i][j].ProgenitorList[k]-1].haloID<<" "<<p[i][j].Merit[k]<<endl;
            }
        }
    }
    }
    else if (opt.outputformat==1) {
    for (int i=opt.numsnapshots-1;i>0;i--) {
        Fout<<i<<"\t"<<h[i].numhalos<<endl;
        for (int j=0;j<h[i].numhalos;j++) {
            Fout<<h[i].Halo[j].haloID<<"\t"<<p[i][j].NumberofProgenitors<<endl;
            for (int k=0;k<p[i][j].NumberofProgenitors;k++) {
                Fout<<h[i-1].Halo[p[i][j].ProgenitorList[k]-1].haloID<<" "<<p[i][j].Merit[k]<<" "<<p[i][j].nsharedfrac[k]<<" ";
            }
        }
    }
    }
    Fout<<"END"<<endl;
    Fout.close();
}

