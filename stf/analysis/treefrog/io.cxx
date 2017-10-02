/*! \file io.cxx
 *  \brief this file contains routines for io
 */

#include "TreeFrog.h"

///Checks if file exits by attempting to get the file attributes
///If success file obviously exists.
///If failure may mean that we don't have permission to access the folder which contains this file or doesn't exist.
///If we need to do that level of checking, lookup return values of stat which will give you more details on why stat failed.
bool FileExists(const char *fname) {
  struct stat stFileInfo;
  bool blnReturn;
  int intStat;
  intStat = stat(fname,&stFileInfo);
  if(intStat == 0) return  true;
  else return false;
}

///\name Reading routines
//@{
///Reads number of halos at each snapshot, useful for mpi decomposition
#ifdef USEMPI
Int_t ReadNumberofHalos(Options &opt, Int_t *numhalos)
{
    fstream Fin;//file is list of halo data files
    string *buf=new string[opt.numsnapshots];
    Int_t tothalos=0;

    Fin.open(opt.fname);
    if (!Fin.is_open()) {
        cerr<<"file containing snapshot list can't be opened"<<endl;
#ifdef USEMPI
        MPI_Abort(MPI_COMM_WORLD,9);
#else
        exit(9);
#endif
    }
    for(int i=0; i<opt.numsnapshots; i++)
    {
        Fin>>buf[i];
            //if (opt.ioformat==DSUSSING) HaloTree[i].Halo=ReadHaloData(&buf[i*1000],HaloTree[i].numhalos);
            //else if (opt.ioformat==DCATALOG)
                numhalos[i]=MPIReadHaloGroupCatalogDataNum(buf[i],opt.nmpifiles, opt.ibinary,opt.ifield, opt.itypematch);
            //else if (opt.ioformat==DNIFTY) HaloTree[i].Halo=ReadNIFTYData(&buf[i*1000],HaloTree[i].numhalos, opt.idcorrectflag);
            tothalos+=numhalos[i];
    }
    Fin.close();
    delete[] buf;
    return tothalos;
}

Int_t ReadNumberofParticlesInHalos(Options &opt, Int_t *numpartinhalos)
{
    fstream Fin;//file is list of halo data files
    string *buf=new string[opt.numsnapshots];
    Int_t tothalos=0,totpart=0;

    Fin.open(opt.fname);
    if (!Fin.is_open()) {
        cerr<<"file containing snapshot list can't be opened"<<endl;
#ifdef USEMPI
        MPI_Abort(MPI_COMM_WORLD,9);
#else
        exit(9);
#endif
    }
    for(int i=0; i<opt.numsnapshots; i++)
    {
        Fin>>buf[i];
            //if (opt.ioformat==DSUSSING) HaloTree[i].Halo=ReadHaloData(&buf[i*1000],HaloTree[i].numhalos);
            //else if (opt.ioformat==DCATALOG)
            numpartinhalos[i]=MPIReadHaloGroupCatalogDataParticleNum(buf[i],opt.nmpifiles, opt.ibinary,opt.ifield, opt.itypematch);
            //else if (opt.ioformat==DNIFTY) HaloTree[i].Halo=ReadNIFTYData(&buf[i*1000],HaloTree[i].numhalos, opt.idcorrectflag);
            //tothalos+=numhalos[i];
            totpart+=numpartinhalos[i];
    }
    Fin.close();
    delete[] buf;
    return totpart;
}
#endif


///Read data from a number of snapshot files
///\todo need to figure out the best way to optimize the openmp reading as it can be unstable as it stands
HaloTreeData *ReadData(Options &opt)
{
    HaloTreeData *HaloTree;
    fstream Fin;//file is list of halo data files
    long unsigned j,nparts,haloid;
    HaloTree=new HaloTreeData[opt.numsnapshots];
    string *buf=new string[opt.numsnapshots];
    Int_t tothalos=0;
    Double_t t0;
    t0=MyGetTime();
#ifdef USEMPI
    Int_t mpi_tothalos;
#endif
    int i,nthreads;
    nthreads=1;
#ifdef USEOPENMP
#pragma omp parallel
    {
            if (omp_get_thread_num()==0) nthreads=omp_get_num_threads();
    }
#endif

#ifdef USEMPI
    if (ThisTask==0) {
        Fin.open(opt.fname);
        if (!Fin.is_open()) {
            cerr<<"file containing snapshot list can't be opened"<<endl;
            MPI_Abort(MPI_COMM_WORLD,9);
        }
        Fin.close();
    }
    MPI_Barrier(MPI_COMM_WORLD);
    for (int itask=0;itask<NProcs;itask++) {
        if (ThisTask==itask) {
            Fin.open(opt.fname);
            for(i=0; i<opt.numsnapshots; i++) Fin>>buf[i];
            Fin.close();
        }
        MPI_Barrier(MPI_COMM_WORLD);
    }
    if (opt.iverbose==1&&ThisTask==0) cout<<"Reading data"<<endl;
#else
    Fin.open(opt.fname);
    if (!Fin.is_open()) {
        cerr<<"file containing snapshot list can't be opened"<<endl;
        exit(9);
    }
    for(i=0; i<opt.numsnapshots; i++) Fin>>buf[i];
    Fin.close();
    if (opt.iverbose==1) cout<<"Reading data"<<endl;
#endif

#if (defined(USEOPENMP) && !defined(USEMPI))
#pragma omp parallel default(shared) \
private(i)
{
#pragma omp for reduction(+:tothalos)
#endif
    for(i=0; i<opt.numsnapshots; i++)
    {
#ifdef USEMPI
        //if mpi only read relavant data
        if (i>=StartSnap && i<EndSnap) {
#endif
            if (opt.ioformat==DSUSSING) HaloTree[i].Halo=ReadHaloData(buf[i],HaloTree[i].numhalos);
            else if (opt.ioformat==DCATALOG) HaloTree[i].Halo=ReadHaloGroupCatalogData(buf[i],HaloTree[i].numhalos, opt.nmpifiles, opt.ibinary,opt.ifield, opt.itypematch,opt.iverbose);
            else if (opt.ioformat==DNIFTY) HaloTree[i].Halo=ReadNIFTYData(buf[i],HaloTree[i].numhalos, opt.idcorrectflag);
            else if (opt.ioformat==DVOID) HaloTree[i].Halo=ReadVoidData(buf[i],HaloTree[i].numhalos, opt.idcorrectflag);
#ifdef USEMPI
            //if mpi then there is data overlap so only add to total if no overlap
            if (ThisTask<NProcs-1 && NProcs>1) {
                if (i<EndSnap-opt.numsteps) tothalos+=HaloTree[i].numhalos;
            }
            else tothalos+=HaloTree[i].numhalos;
#else
            tothalos+=HaloTree[i].numhalos;
#endif


#ifdef USEMPI
        }
#endif

    }
#if (defined(USEOPENMP) && !defined(USEMPI))
}
#endif

#ifdef USEMPI
    cout<<ThisTask<<" has ---- "<<tothalos<<endl;
    mpi_tothalos=tothalos;
    MPI_Allreduce(&mpi_tothalos, &tothalos, 1, MPI_Int_t, MPI_SUM, MPI_COMM_WORLD);
    if (ThisTask==0) cout<<" all tasks have "<<tothalos<<endl;
#endif
    opt.TotalNumberofHalos=tothalos;
    if (opt.iverbose==1) cout<<"Finished reading data "<<MyGetTime()-t0<<endl;

    return HaloTree;
}
//@}

///\name Write routines
//@{
void WriteHaloMergerTree(Options &opt, ProgenitorData **p, HaloTreeData *h) {
    fstream Fout;
    char fname[1000];
    char fnamempi[2000];
    int istep;
    sprintf(fname,"%s",opt.outname);
    Double_t time1=MyGetTime();
#ifndef USEMPI
    int ThisTask=0, NProcs=1;
    int StartSnap=0,EndSnap=opt.numsnapshots;
#endif
#ifdef USEHDF
    H5File Fhdf;
    H5std_string datasetname;
    DataSpace dataspace;
    DataSet dataset;
    Attribute attr;
    DataSpace attrspace;
    hsize_t *dims;
    int rank;
    DataSpace *propdataspace;
    DataSet *propdataset;
    HDFCatalogNames hdfnames;
    int itemp=0;
#endif
    int istart,iend;

    if (ThisTask==0) cout<<"Writing to "<<fname<<endl;
    if (opt.outputformat==OUTASCII)
    {
        //if mpi then can have a single aggregate file or each thread write their own
        if (NProcs>1) {
#ifdef USEMPI
        if (opt.iwriteparallel==1 && ThisTask==0) cout<<"Writing files in parallel "<<endl;
        //now if mpi then last task writes header
        //all tasks starting from last and moves backwards till it reaches its
        //startpoint+numofsteps used to produce links
        //save first task which goes all the way to its StartSnap
        iend=EndSnap;
        istart=StartSnap+opt.numsteps;
        if (opt.iwriteparallel==1) sprintf(fnamempi,"%s.mpi_task-%d.isnap-%d.fsnap-%d",opt.outname,ThisTask,istart,iend);
        else sprintf(fnamempi,"%s",fname);
        if (ThisTask==0) istart=0;
        for (int itask=NProcs-1;itask>=0;itask--) {
            if (ThisTask==itask) {
                if (ThisTask==NProcs-1) {
                    Fout.open(fnamempi,ios::out);
                    Fout<<opt.numsnapshots<<endl;
                    Fout<<opt.description<<endl;
                    Fout<<opt.TotalNumberofHalos<<endl;
                }
                else {
                    Fout.open(fnamempi,ios::out | ios::app);
                }
                if (opt.iverbose)cout<<ThisTask<<" starting to write "<<fnamempi<<" for "<<iend<<" down to "<<istart<<flush<<endl;
                for (int i=opt.numsnapshots-1;i>0;i--) if (i>=istart && i<iend) {
                    Fout<<i+opt.snapshotvaloffset<<"\t"<<h[i].numhalos<<endl;
                    for (int j=0;j<h[i].numhalos;j++) {
                        Fout<<h[i].Halo[j].haloID<<"\t"<<p[i][j].NumberofProgenitors;
                        if (opt.outdataformat>=2) {
                            Fout<<"\t"<<h[i].Halo[j].NumberofParticles;
                        }
                        Fout<<endl;
                        for (int k=0;k<p[i][j].NumberofProgenitors;k++) {
                            Fout<<h[i-p[i][j].istep].Halo[p[i][j].ProgenitorList[k]-1].haloID<<" ";
                            if (opt.outdataformat>=1) {
                                Fout<<p[i][j].Merit[k]<<" ";
                            }
                            if (opt.outdataformat>=2) {
                                Fout<<h[i-p[i][j].istep].Halo[p[i][j].ProgenitorList[k]-1].NumberofParticles<<" ";
                            }
                            Fout<<endl;
                        }
                    }
                }
                Fout.close();
            }
            if (opt.iwriteparallel==0) MPI_Barrier(MPI_COMM_WORLD);
        }
        if (ThisTask==0) {
            Fout.open(fnamempi,ios::out | ios::app);
            ///last file has no connections
            Fout<<0+opt.snapshotvaloffset<<"\t"<<h[0].numhalos<<endl;
            for (int j=0;j<h[0].numhalos;j++) {
                Fout<<h[0].Halo[j].haloID<<"\t"<<0;
                if (opt.outdataformat>=2) {
                    Fout<<"\t"<<h[0].Halo[j].NumberofParticles;
                }
                Fout<<endl;
            }
            Fout<<"END"<<endl;
            Fout.close();
        }
#endif
        }
        //if not mpi (ie: NProcs==1)
        else{
            Fout.open(fname,ios::out);
            Fout<<opt.numsnapshots<<endl;
            Fout<<opt.description<<endl;
            Fout<<opt.TotalNumberofHalos<<endl;
            for (int i=opt.numsnapshots-1;i>0;i--) {
                Fout<<i+opt.snapshotvaloffset<<"\t"<<h[i].numhalos<<endl;
                for (int j=0;j<h[i].numhalos;j++) {
                    Fout<<h[i].Halo[j].haloID<<"\t"<<p[i][j].NumberofProgenitors;
                    if (opt.outdataformat>=2) {
                        Fout<<"\t"<<h[i].Halo[j].NumberofParticles;
                    }
                    Fout<<endl;
                    for (int k=0;k<p[i][j].NumberofProgenitors;k++) {
                        Fout<<h[i-p[i][j].istep].Halo[p[i][j].ProgenitorList[k]-1].haloID<<" ";
                        if (opt.outdataformat>=1) {
                            Fout<<p[i][j].Merit[k]<<" ";
                        }
                        if (opt.outdataformat>=2) {
                            Fout<<h[i-p[i][j].istep].Halo[p[i][j].ProgenitorList[k]-1].NumberofParticles<<" ";
                        }
                        Fout<<endl;
                    }
                }
            }
            Fout<<0+opt.snapshotvaloffset<<"\t"<<h[0].numhalos<<endl;
            ///last file has no connections
            for (int j=0;j<h[0].numhalos;j++) {
                Fout<<h[0].Halo[j].haloID<<"\t"<<0;
                if (opt.outdataformat>=2) {
                    Fout<<"\t"<<h[0].Halo[j].NumberofParticles;
                }
                Fout<<endl;
            }
            Fout<<"END"<<endl;
            Fout.close();
        }
    }
    //end of ascii output
#ifdef USEHDF
    else if (opt.outputformat==OUTHDF)
    {
#ifdef USEMPI
        if (opt.iwriteparallel==1 && ThisTask==0) cout<<"Writing files in parallel "<<endl;
        if (opt.iwriteparallel==1) sprintf(fname,"%s.mpi_task-%d.isnap-%d.fsnap-%d",opt.outname,ThisTask,istart,iend);
        else sprintf(fname,"%s",fname);
#else
        sprintf(fname,"%s",fname);
#endif
        //write header information
#ifdef USEMPI
        if ((opt.iwriteparallel==0 && ThisTask==0) || (opt.iwriteparallel==1))
#endif
        {
            Fhdf=H5File(fname,H5F_ACC_TRUNC);
            attrspace=DataSpace(H5S_SCALAR);
            attr=Fhdf.createAttribute("Number_of_snapshots", PredType::STD_U32LE, attrspace);
            attr.write(PredType::STD_U32LE,&opt.numsnapshots);
            attrspace=DataSpace(H5S_SCALAR);
            attr=Fhdf.createAttribute("Total_number_of_halos", PredType::STD_U64LE, attrspace);
            attr.write(PredType::STD_U64LE,&opt.TotalNumberofHalos);
            attrspace=DataSpace(H5S_SCALAR);
            attr=Fhdf.createAttribute("Merit_limit", PredType::NATIVE_DOUBLE, attrspace);
            attr.write(PredType::NATIVE_DOUBLE,&opt.mlsig);
            //for multistep info
            attrspace=DataSpace(H5S_SCALAR);
            attr=Fhdf.createAttribute("Number_of_steps", PredType::STD_U32LE, attrspace);
            attr.write(PredType::STD_U32LE,&opt.numsteps);
            attrspace=DataSpace(H5S_SCALAR);
            attr=Fhdf.createAttribute("Search_next_step_criterion", PredType::STD_U32LE, attrspace);
            attr.write(PredType::STD_U32LE,&opt.imultsteplinkcrit);
            attrspace=DataSpace(H5S_SCALAR);
            attr=Fhdf.createAttribute("Merit_limit_for_next_step", PredType::NATIVE_DOUBLE, attrspace);
            attr.write(PredType::NATIVE_DOUBLE,&opt.meritlimit);
            //for core matching info
            attrspace=DataSpace(H5S_SCALAR);
            attr=Fhdf.createAttribute("Core_fraction", PredType::NATIVE_DOUBLE, attrspace);
            attr.write(PredType::NATIVE_DOUBLE,&opt.particle_frac);
            attrspace=DataSpace(H5S_SCALAR);
            attr=Fhdf.createAttribute("Core_min_number_of_particles", PredType::STD_U32LE, attrspace);
            attr.write(PredType::STD_U32LE,&opt.min_numpart);
            // general description
            // Create new string datatype for attribute
            StrType strdatatype(PredType::C_S1, 1000);
            // Set up write buffer for attribute
            const H5std_string strwritebuf (opt.description);
            attr = Fhdf.createAttribute("Description", strdatatype, attrspace);
            attr.write(strdatatype, strwritebuf);
            Fhdf.close();
        }

        iend=EndSnap;
        istart=StartSnap+opt.numsteps;
        if (ThisTask==0) istart=0;
        //having written header if necessary, reopen and append to hdf file
        for (int itask=NProcs-1;itask>=0;itask--) {
            Fhdf=H5File(fname,H5F_ACC_RDWR);
            /*
            if (opt.outputformat==0) {
            for (int i=opt.numsnapshots-1;i>0;i--) if (i>=istart && i<iend) {
                Fout<<i+opt.snapshotvaloffset<<"\t"<<h[i].numhalos<<endl;
                for (int j=0;j<h[i].numhalos;j++) {
                    Fout<<h[i].Halo[j].haloID<<"\t"<<p[i][j].NumberofProgenitors<<endl;
                    for (int k=0;k<p[i][j].NumberofProgenitors;k++) {
                        Fout<<h[i-p[i][j].istep].Halo[p[i][j].ProgenitorList[k]-1].haloID<<endl;
                    }
                }
            }
            }
            else {
            for (int i=opt.numsnapshots-1;i>0;i--) if (i>=istart && i<iend) {
                Fout<<i+opt.snapshotvaloffset<<"\t"<<h[i].numhalos<<endl;
                for (int j=0;j<h[i].numhalos;j++) {
                    Fout<<h[i].Halo[j].haloID<<"\t"<<p[i][j].NumberofProgenitors<<endl;
                    for (int k=0;k<p[i][j].NumberofProgenitors;k++) {
                        Fout<<h[i-p[i][j].istep].Halo[p[i][j].ProgenitorList[k]-1].haloID<<" "<<p[i][j].Merit[k]<<endl;
                    }
                }
            }
            }
            */
            //might need to be careful about stagging writes to the file
            Fhdf.close();
        }
        if (ThisTask==0) {
            Fhdf=H5File(fname,H5F_ACC_RDWR);
            /*
            ///last file has no connections
            Fout<<0+opt.snapshotvaloffset<<"\t"<<h[0].numhalos<<endl;
            for (int j=0;j<h[0].numhalos;j++) {
                Fout<<h[0].Halo[j].haloID<<"\t"<<0<<endl;
            }
            */
            Fhdf.close();
        }
    }
#endif
    if (ThisTask==0) cout<<"Done writing to "<<fname<<" "<<MyGetTime()-time1<<endl;
}

/// same as \ref WriteHaloMergerTree  but going reverse direction using DescendantData
void WriteHaloMergerTree(Options &opt, DescendantData **p, HaloTreeData *h) {
    fstream Fout;
    char fname[1000];
    char fnamempi[2000];
    int istep;
    sprintf(fname,"%s",opt.outname);
    Double_t time1=MyGetTime();
#ifndef USEMPI
    int ThisTask=0, NProcs=1;
    int StartSnap=0,EndSnap=opt.numsnapshots;
#endif
#ifdef USEHDF
    H5File Fhdf;
    H5std_string datasetname;
    DataSpace dataspace;
    DataSet dataset;
    Attribute attr;
    DataSpace attrspace;
    hsize_t *dims;
    int rank;
    DataSpace *propdataspace;
    DataSet *propdataset;
    HDFCatalogNames hdfnames;
    int itemp=0;
#endif
    int istart,iend;

    if (ThisTask==0) cout<<"Writing descedant based halo merger tree to "<<fname<<endl;
    if (opt.outputformat==OUTASCII)
    {
        //if mpi then can have a single aggregate file or each thread write their own
        if (NProcs>1) {
#ifdef USEMPI
        if (opt.iwriteparallel==1 && ThisTask==0) cout<<"Writing files in parallel "<<endl;
        //now if mpi then last task writes header
        //all tasks starting from last and moves backwards till it reaches its
        //startpoint+numofsteps used to produce links
        //save first task which goes all the way to its StartSnap
        iend=EndSnap-opt.numsteps;
        istart=StartSnap;
        if (opt.iwriteparallel==1) sprintf(fnamempi,"%s.mpi_task-%d.isnap-%d.fsnap-%d",opt.outname,ThisTask,istart,iend);
        else sprintf(fnamempi,"%s",fname);
        if (ThisTask==0) istart=0;
        if (ThisTask==NProcs-1) iend=opt.numsnapshots-1;
        for (int itask=0;itask<NProcs;itask++) {
            if (ThisTask==itask) {
                if (ThisTask==0) {
                    Fout.open(fnamempi,ios::out);
                    Fout<<opt.numsnapshots<<endl;
                    Fout<<opt.description<<endl;
                    Fout<<opt.TotalNumberofHalos<<endl;
                }
                else {
                    Fout.open(fnamempi,ios::out | ios::app);
                }
                if (opt.iverbose)cout<<ThisTask<<" starting to write "<<fnamempi<<" for "<<iend<<" down to "<<istart<<flush<<endl;
                for (int i=0;i<opt.numsnapshots-1;i++) if (i>=istart && i<iend) {
                    Fout<<i+opt.snapshotvaloffset<<"\t"<<h[i].numhalos<<endl;
                    for (int j=0;j<h[i].numhalos;j++) {
                        Fout<<h[i].Halo[j].haloID<<"\t"<<p[i][j].NumberofDescendants;
                        if (opt.outdataformat>=2) {
                            Fout<<"\t"<<h[i].Halo[j].NumberofParticles;
                        }
                        Fout<<endl;
                        for (int k=0;k<p[i][j].NumberofDescendants;k++) {
                            Fout<<h[i+p[i][j].istep].Halo[p[i][j].DescendantList[k]-1].haloID<<" ";
                            Fout<<p[i][j].dtoptype[k]<<" ";
                            if (opt.outdataformat>=1) {
                                Fout<<p[i][j].Merit[k]<<" ";
                            }
                            if (opt.outdataformat>=2) {
                                Fout<<h[i+p[i][j].istep].Halo[p[i][j].DescendantList[k]-1].NumberofParticles<<" ";
                            }
                            Fout<<endl;
                        }
                    }
                }
                Fout.close();
            }
            if (opt.iwriteparallel==0) MPI_Barrier(MPI_COMM_WORLD);
        }
        if (ThisTask==NProcs-1) {
            Fout.open(fnamempi,ios::out | ios::app);
            ///last file has no connections
            Fout<<opt.numsnapshots-1+opt.snapshotvaloffset<<"\t"<<h[opt.numsnapshots-1].numhalos<<endl;
            for (int j=0;j<h[opt.numsnapshots-1].numhalos;j++) {
                Fout<<h[opt.numsnapshots-1].Halo[j].haloID<<"\t"<<0;
                if (opt.outdataformat>=2) {
                    Fout<<"\t"<<h[opt.numsnapshots-1].Halo[j].NumberofParticles;
                }
                Fout<<endl;
            }
            Fout<<"END"<<endl;
            Fout.close();
        }
#endif
        }
        //if not mpi (ie: NProcs==1)
        else{
            Fout.open(fname,ios::out);
            Fout<<opt.numsnapshots<<endl;
            Fout<<opt.description<<endl;
            Fout<<opt.TotalNumberofHalos<<endl;
            for (int i=0;i<opt.numsnapshots-1;i++) {
                Fout<<i+opt.snapshotvaloffset<<"\t"<<h[i].numhalos<<endl;
                for (int j=0;j<h[i].numhalos;j++) {
                    Fout<<h[i].Halo[j].haloID<<"\t"<<p[i][j].NumberofDescendants;
                    if (opt.outdataformat>=2) {
                        Fout<<"\t"<<h[i].Halo[j].NumberofParticles;
                    }
                    Fout<<endl;
                    for (int k=0;k<p[i][j].NumberofDescendants;k++) {
                        Fout<<h[i+p[i][j].istep].Halo[p[i][j].DescendantList[k]-1].haloID<<" ";
                        Fout<<p[i][j].dtoptype[k]<<" ";
                        if (opt.outdataformat>=1) {
                            Fout<<p[i][j].Merit[k]<<" ";
                        }
                        if (opt.outdataformat>=2) {
                            Fout<<h[i+p[i][j].istep].Halo[p[i][j].DescendantList[k]-1].NumberofParticles<<" ";
                        }
                        Fout<<endl;
                    }
                }
            }
            Fout<<opt.numsnapshots-1+opt.snapshotvaloffset<<"\t"<<h[opt.numsnapshots-1].numhalos<<endl;
            ///last file has no connections
            for (int j=0;j<h[opt.numsnapshots-1].numhalos;j++) {
                Fout<<h[opt.numsnapshots-1].Halo[j].haloID<<"\t"<<0;
                if (opt.outdataformat>=2) {
                    Fout<<"\t"<<h[opt.numsnapshots-1].Halo[j].NumberofParticles;
                }
                Fout<<endl;
            }
            Fout<<"END"<<endl;
            Fout.close();
        }
    }
    //end of ascii output
#ifdef USEHDF
    else if (opt.outputformat==OUTHDF)
    {

    }
#endif
    if (ThisTask==0) cout<<"Done writing to "<<fname<<" "<<MyGetTime()-time1<<endl;
}

void WriteHaloGraph(Options &opt, ProgenitorData **p, DescendantData **d, HaloTreeData *h) {
    fstream Fout;
    char fname[1000];
    sprintf(fname,"%s",opt.outname);
    cout<<"saving halo merger tree to "<<fname<<endl;
    Fout.open(fname,ios::out);
    Fout<<opt.numsnapshots<<endl;
    Fout<<opt.description<<endl;
    Fout<<opt.TotalNumberofHalos<<endl;
    for (int i=opt.numsnapshots-1;i>=0;i--) {
        Fout<<i+opt.snapshotvaloffset<<"\t"<<h[i].numhalos<<endl;
        for (int j=0;j<h[i].numhalos;j++) {
            if (i==opt.numsnapshots-1) Fout<<h[i].Halo[j].haloID<<"\t"<<p[i][j].NumberofProgenitors<<"\t"<<0;
            else if (i==0)Fout<<h[i].Halo[j].haloID<<"\t"<<0<<"\t"<<d[i][j].NumberofDescendants;
            else Fout<<h[i].Halo[j].haloID<<"\t"<<p[i][j].NumberofProgenitors<<"\t"<<d[i][j].NumberofDescendants;
            if (opt.outdataformat>=2) {
                Fout<<"\t"<<h[i].Halo[j].NumberofParticles;
            }
            Fout<<endl;
            if (i>0) {
                for (int k=0;k<p[i][j].NumberofProgenitors;k++) {
                    Fout<<h[i-p[i][j].istep].Halo[p[i][j].ProgenitorList[k]-1].haloID<<" ";
                    if (opt.outdataformat>=1) {
                        Fout<<p[i][j].Merit[k]<<" ";
                    }
                    if (opt.outdataformat>=2) {
                        Fout<<h[i-p[i][j].istep].Halo[p[i][j].ProgenitorList[k]-1].NumberofParticles<<" ";
                    }
                    Fout<<endl;
                }
            }
            if (i<opt.numsnapshots-1) {
                for (int k=0;k<d[i][j].NumberofDescendants;k++) {
                    Fout<<h[i+d[i][j].istep].Halo[d[i][j].DescendantList[k]-1].haloID<<" ";
                    Fout<<d[i][j].dtoptype[k]<<" ";
                    if (opt.outdataformat>=1) {
                        Fout<<d[i][j].Merit[k]<<" ";
                    }
                    if (opt.outdataformat>=2) {
                        Fout<<h[i+d[i][j].istep].Halo[d[i][j].DescendantList[k]-1].NumberofParticles<<" ";
                    }
                    Fout<<endl;
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
    for (int i=opt.numsnapshots-1;i>0;i--) {
        Fout<<i<<"\t"<<h[i].numhalos<<endl;
        for (int j=0;j<h[i].numhalos;j++) {
            Fout<<h[i].Halo[j].haloID<<"\t"<<p[i][j].NumberofProgenitors<<endl;
            for (int k=0;k<p[i][j].NumberofProgenitors;k++) {
                Fout<<h[i-1].Halo[p[i][j].ProgenitorList[k]-1].haloID<<" "<<p[i][j].Merit[k]<<endl;
            }
        }
    }
    Fout<<"END"<<endl;
    Fout.close();
}
//@}

///\name write associated information of treefrog run
//@{
///if a memory efficient particle id to index map was produced save the file
void SavePIDStoIndexMap(Options &opt,map<IDTYPE, IDTYPE>&idmap)
{
    char fname[1000];
    fstream Fout;
    IDTYPE *keys,*indices;
    Int_t i;
    size_t idsize,mapsize;
#ifndef USEMPI
    int ThisTask=0;
#endif
    if (ThisTask==0)
    {
        sprintf(fname,"%s.pidtoindexmap.dat",opt.outname);
        Fout.open(fname,ios::out|ios::binary);
        cout<<"Writing unique memory efficent mapping for particle IDS to index to "<<fname<<endl;
        //write header information
        Fout.write((char*)&opt.numsnapshots,sizeof(int));
        Fout.write((char*)&opt.TotalNumberofHalos,sizeof(long unsigned));
        idsize=sizeof(IDTYPE);
        Fout.write((char*)&idsize,sizeof(size_t));
        mapsize=idmap.size();
        Fout.write((char*)&mapsize,sizeof(size_t));
        keys=new IDTYPE[mapsize];
        indices=new IDTYPE[mapsize];
        i=0;
        for (auto it: idmap) {
            keys[i]=it.first;
            indices[i]=it.second;
            i++;
        }
        Fout.write((char*)keys,sizeof(IDTYPE)*mapsize);
        Fout.write((char*)indices,sizeof(IDTYPE)*mapsize);
        Fout.close();
        delete[] keys;
        delete[] indices;
    }
}
///if a memory efficient particle id to index map file exists read it
int ReadPIDStoIndexMap(Options &opt,map<IDTYPE, IDTYPE>&idmap)
{
    char fname[1000];
    fstream Fout;
    IDTYPE *keys,*indices;
    Int_t i;
    size_t idsize,mapsize;
    int numsnap;
    long unsigned tothalo;
    int iflag=0;
    double time1;
#ifndef USEMPI
    int ThisTask=0;
#endif
    if (ThisTask==0) {
        sprintf(fname,"%s.pidtoindexmap.dat",opt.outname);
        Fout.open(fname,ios::in|ios::binary);
        cout<<"Attempting to read unique memory efficent mapping for particle IDS to index to "<<fname<<endl;
        //write header information
        Fout.read((char*)&numsnap,sizeof(int));
        Fout.read((char*)&tothalo,sizeof(long unsigned));
        Fout.read((char*)&idsize,sizeof(size_t));
        //if all is well then keep reading
        if (opt.numsnapshots==numsnap && opt.TotalNumberofHalos && idsize==sizeof(IDTYPE)) {
            iflag=1;
            cout<<"Reading information"<<endl;
            Fout.read((char*)&mapsize,sizeof(size_t));
            keys=new IDTYPE[mapsize];
            indices=new IDTYPE[mapsize];
            Fout.read((char*)keys,sizeof(IDTYPE)*mapsize);
            Fout.read((char*)indices,sizeof(IDTYPE)*mapsize);
        }
        Fout.close();
    }
#ifdef USEMPI
    MPI_Bcast(&iflag, 1, MPI_INT, 0, MPI_COMM_WORLD);
#endif
    if (iflag==0) {
        if (ThisTask==0) cout<<"Unable to read data, must generate map "<<endl;
        return iflag;
    }
#ifdef USEMPI
    //communicate information
#endif
    if (ThisTask==0) cout<<"And producing internal map "<<endl;
    time1=MyGetTime();
    for (i=0;i<mapsize;i++) idmap.insert(pair<IDTYPE, IDTYPE>(keys[i],indices[i]));
    if (ThisTask==0) cout<<"Took "<<MyGetTime()-time1<<endl;
    return iflag;
}
//@}
