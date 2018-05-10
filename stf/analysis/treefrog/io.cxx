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
unsigned long ReadNumberofHalos(Options &opt, unsigned long *numhalos)
{
    fstream Fin;//file is list of halo data files
    string *buf=new string[opt.numsnapshots];
    unsigned long tothalos=0;

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

unsigned long ReadNumberofParticlesInHalos(Options &opt, unsigned long *numpartinhalos)
{
    fstream Fin;//file is list of halo data files
    string *buf=new string[opt.numsnapshots];
    unsigned long tothalos=0,totpart=0;

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
    hsize_t dims[1], chunk_dims[1];
    int rank;
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
                        if (opt.outdataformat>=DATAOUTMERITNPART) {
                            Fout<<"\t"<<h[i].Halo[j].NumberofParticles;
                        }
                        Fout<<endl;
                        for (int k=0;k<p[i][j].NumberofProgenitors;k++) {
                            Fout<<h[i-p[i][j].istep].Halo[p[i][j].ProgenitorList[k]-1].haloID<<" ";
                            if (opt.outdataformat>=DATAOUTMERIT) {
                                Fout<<p[i][j].Merit[k]<<" ";
                            }
                            if (opt.outdataformat>=DATAOUTMERITNPART) {
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
                if (opt.outdataformat>=DATAOUTMERITNPART) {
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
                    if (opt.outdataformat>=DATAOUTMERITNPART) {
                        Fout<<"\t"<<h[i].Halo[j].NumberofParticles;
                    }
                    Fout<<endl;
                    for (int k=0;k<p[i][j].NumberofProgenitors;k++) {
                        Fout<<h[i-p[i][j].istep].Halo[p[i][j].ProgenitorList[k]-1].haloID<<" ";
                        if (opt.outdataformat>=DATAOUTMERIT) {
                            Fout<<p[i][j].Merit[k]<<" ";
                        }
                        if (opt.outdataformat>=DATAOUTMERITNPART) {
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
                if (opt.outdataformat>=DATAOUTMERITNPART) {
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
        //If hdf5 then write a tree file per snapshot meaning it can be written out in parallel
        iend=EndSnap;
        istart=StartSnap+opt.numsteps;
        if (ThisTask==0) istart=0;

        for (int i=opt.numsnapshots-1;i>0;i--) if (i>=istart && i<iend) {


            sprintf(fname,"%s/snapshot_%03d.VELOCIraptor.tree",opt.outname,i+opt.snapshotvaloffset);
            cout<<ThisTask<<" is writing to "<<fname<<endl;

            //Header information
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


            //Set the datasets properties
            dims[0]=h[i].numhalos;
            rank=1;
            // Set the minmum chunk size
            chunk_dims[0]=min((unsigned long)HDFOUTPUTCHUNKSIZE,(unsigned long)h[i].numhalos);
            //Create dataset proplist
            DSetCreatPropList hdfdatasetproplist1;

            // ID
            dataspace = DataSpace(rank,dims);
            datasetname=H5std_string("ID");

            // Check if there are halos to output so it can be compressed
            if (chunk_dims[0]>0) {
                // Modify dataset creation property to enable chunking
                hdfdatasetproplist1.setChunk(rank, chunk_dims);
                // Set ZLIB (DEFLATE) Compression using level 6.
                hdfdatasetproplist1.setDeflate(6);
                dataset = Fhdf.createDataSet(datasetname,PredType::STD_U64LE,dataspace,hdfdatasetproplist1);
            }
            else {
                dataset = Fhdf.createDataSet(datasetname,PredType::STD_U64LE,dataspace);
            }

            //Write out the dataset
            long unsigned *data1 = new long unsigned[h[i].numhalos];
            for(Int_t j=0;j<h[i].numhalos;j++) data1[j]=h[i].Halo[j].haloID;
            dataset.write(data1,PredType::STD_U64LE);

            delete[] data1;

            //Num progenitors
            dataspace = DataSpace(rank,dims);
            datasetname=H5std_string("NumProgen");

            // Check if there are halos to output so it can be compressed
            if (chunk_dims[0]>0) {
                // Modify dataset creation property to enable chunking
                hdfdatasetproplist1.setChunk(rank, chunk_dims);
                // Set ZLIB (DEFLATE) Compression using level 6.
                hdfdatasetproplist1.setDeflate(6);
                dataset = Fhdf.createDataSet(datasetname,PredType::STD_I32LE,dataspace,hdfdatasetproplist1);
            }
            else {
                dataset = Fhdf.createDataSet(datasetname,PredType::STD_I32LE,dataspace);
            }

            //Write out the dataset
            int *data2 = new int[h[i].numhalos];
            for(Int_t j=0;j<h[i].numhalos;j++) data2[j]=p[i][j].NumberofProgenitors;
            dataset.write(data2,PredType::STD_I32LE);

            delete[] data2;

            // Number of particles
            if (opt.outdataformat>=DATAOUTMERITNPART) {
                dataspace = DataSpace(rank,dims);
                datasetname=H5std_string("Npart");

                // Check if there are halos to output so it can be compressed
                if (chunk_dims[0]>0) {
                    // Modify dataset creation property to enable chunking
                    hdfdatasetproplist1.setChunk(rank, chunk_dims);
                    // Set ZLIB (DEFLATE) Compression using level 6.
                    hdfdatasetproplist1.setDeflate(6);
                    dataset = Fhdf.createDataSet(datasetname,PredType::STD_U64LE,dataspace,hdfdatasetproplist1);
                }
                else {
                    dataset = Fhdf.createDataSet(datasetname,PredType::STD_U64LE,dataspace);
                }

                //Write out the dataset
                long unsigned *data3 = new long unsigned[h[i].numhalos];
                for(Int_t j=0;j<h[i].numhalos;j++)  data3[j]=h[i].Halo[j].NumberofParticles;
                dataset.write(data3,PredType::STD_U64LE);

                delete[] data3;
            }

            // Keep track of the number of progenitors
            long unsigned totnprogen=0;

            //Progenitors offsets
            dataspace = DataSpace(rank,dims);
            datasetname=H5std_string("ProgenOffsets");

            // Check if there are halos to output so it can be compressed
            if (chunk_dims[0]>0) {
                // Modify dataset creation property to enable chunking
                hdfdatasetproplist1.setChunk(rank, chunk_dims);
                // Set ZLIB (DEFLATE) Compression using level 6.
                hdfdatasetproplist1.setDeflate(6);
                dataset = Fhdf.createDataSet(datasetname,PredType::STD_U64LE,dataspace,hdfdatasetproplist1);
            }
            else {
                dataset = Fhdf.createDataSet(datasetname,PredType::STD_U64LE,dataspace);
            }

            //Write out the dataset
            long unsigned *data4 = new long unsigned[h[i].numhalos];
            for(Int_t j=0;j<h[i].numhalos;j++){
                data4[j]=totnprogen;
                totnprogen+=p[i][j].NumberofProgenitors;
            }
            dataset.write(data4,PredType::STD_U64LE);

            delete[] data4;


            //Set the new dataset parameters
            long unsigned itemp=0;
            dims[0]=totnprogen;
            rank=1;
            //Set minimum chunk size
            chunk_dims[0]=min((unsigned long)HDFOUTPUTCHUNKSIZE,totnprogen);
            //Create dataset proplist
            DSetCreatPropList hdfdatasetproplist2;

            //Progenitor IDs
            dataspace = DataSpace(rank,dims);
            datasetname=H5std_string("Progenitors");


            // Check if there are halos to output so it can be compressed
            if (chunk_dims[0]>0) {
                // Modify dataset creation property to enable chunking
                hdfdatasetproplist2.setChunk(rank, chunk_dims);
                // Set ZLIB (DEFLATE) Compression using level 6.
                hdfdatasetproplist2.setDeflate(6);
                dataset = Fhdf.createDataSet(datasetname,PredType::STD_U64LE,dataspace,hdfdatasetproplist2);
            }
            else {
                dataset = Fhdf.createDataSet(datasetname,PredType::STD_U64LE,dataspace);
            }

            //Write out the dataset
            long unsigned *data5= new long unsigned[totnprogen];
            for(Int_t j=0;j<h[i].numhalos;j++)
                for (Int_t k=0;k<p[i][j].NumberofProgenitors;k++)
                    data5[itemp++] = h[i-p[i][j].istep].Halo[p[i][j].ProgenitorList[k]-1].haloID;
            dataset.write(data5,PredType::STD_U64LE);

            delete[] data5;

            //Merits
            if (opt.outdataformat>=DATAOUTMERIT) {
                itemp=0;

                dataspace = DataSpace(rank,dims);
                datasetname=H5std_string("Merits");
                // Check if there are halos to output so it can be compressed
                if (chunk_dims[0]>0) {
                    // Modify dataset creation property to enable chunking
                    hdfdatasetproplist2.setChunk(rank, chunk_dims);
                    // Set ZLIB (DEFLATE) Compression using level 6.
                    hdfdatasetproplist2.setDeflate(6);
                    dataset = Fhdf.createDataSet(datasetname,PredType::NATIVE_FLOAT,dataspace,hdfdatasetproplist2);
                }
                else {
                    dataset = Fhdf.createDataSet(datasetname,PredType::NATIVE_FLOAT,dataspace);
                }

                //Write out the dataset
                float *data6 = new float[totnprogen];
                for(Int_t j=0;j<h[i].numhalos;j++)
                    for (Int_t k=0;k<p[i][j].NumberofProgenitors;k++)
                        data6[itemp++] = p[i][j].Merit[k];
                dataset.write(data6,PredType::NATIVE_FLOAT);

                delete[] data6;
            }

            //Progenitor number of particles
            if (opt.outdataformat>=DATAOUTMERITNPART) {
                itemp=0;

                dataspace = DataSpace(rank,dims);
                datasetname=H5std_string("ProgenNpart");

                // Check if there are halos to output so it can be compressed
                if (chunk_dims[0]>0) {
                    // Modify dataset creation property to enable chunking
                    hdfdatasetproplist2.setChunk(rank, chunk_dims);
                    // Set ZLIB (DEFLATE) Compression using level 6.
                    hdfdatasetproplist2.setDeflate(6);
                    dataset = Fhdf.createDataSet(datasetname,PredType::STD_U64LE,dataspace,hdfdatasetproplist2);
                }
                else {
                    dataset = Fhdf.createDataSet(datasetname,PredType::STD_U64LE,dataspace);
                }

                //Write out the dataset
                long unsigned *data7 = new long unsigned[totnprogen];
                for(Int_t j=0;j<h[i].numhalos;j++)
                    for (Int_t k=0;k<p[i][j].NumberofProgenitors;k++)
                        data7[itemp++] = h[i-p[i][j].istep].Halo[p[i][j].ProgenitorList[k]-1].NumberofParticles;
                dataset.write(data7,PredType::STD_U64LE);

                delete[] data7;
            }

            Fhdf.close();

            cout<<ThisTask<<" Done writing to "<<fname<<" "<<MyGetTime()-time1<<endl;

        }

        ///last file has no connections
        if(ThisTask==0){

            sprintf(fname,"%s/snapshot_%03d.VELOCIraptor.tree",opt.outname,0+opt.snapshotvaloffset);
            cout<<ThisTask<<" is writing to "<<fname<<endl;

            //Header information
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


            ///last file has no connections
            dims[0]=h[0].numhalos;
            rank=1;
            //Set minimum chunk size
            chunk_dims[0]=min((unsigned long)HDFOUTPUTCHUNKSIZE,(unsigned long)h[0].numhalos);
            //Create dataset proplist
            DSetCreatPropList hdfdatasetproplist;

            //ID
            dataspace = DataSpace(rank,dims);
            datasetname=H5std_string("ID");

            // Check if there are halos to output so it can be compressed
            if (chunk_dims[0]>0) {
                // Modify dataset creation property to enable chunking
                hdfdatasetproplist.setChunk(rank, chunk_dims);
                // Set ZLIB (DEFLATE) Compression using level 6.
                hdfdatasetproplist.setDeflate(6);
                dataset = Fhdf.createDataSet(datasetname,PredType::STD_U64LE,dataspace,hdfdatasetproplist);
            }
            else {
                dataset = Fhdf.createDataSet(datasetname,PredType::STD_U64LE,dataspace);
            }

            //Write out the dataset
            long unsigned *data8 = new long unsigned[h[0].numhalos];
            for(Int_t j=0;j<h[0].numhalos;j++) data8[j]=h[0].Halo[j].haloID;
            dataset.write(data8,PredType::STD_U64LE);

            delete[] data8;

            //Num Progenitors
            dataspace = DataSpace(rank,dims);
            datasetname=H5std_string("NumProgen");

            // Check if there are halos to output so it can be compressed
            if (chunk_dims[0]>0) {
                // Modify dataset creation property to enable chunking
                hdfdatasetproplist.setChunk(rank, chunk_dims);
                // Set ZLIB (DEFLATE) Compression using level 6.
                hdfdatasetproplist.setDeflate(6);
                dataset = Fhdf.createDataSet(datasetname,PredType::STD_I32LE,dataspace,hdfdatasetproplist);
            }
            else {
                dataset = Fhdf.createDataSet(datasetname,PredType::STD_I32LE,dataspace);
            }

            //Write out the dataset
            int *data9 = new int[h[0].numhalos];
            for(Int_t j=0;j<h[0].numhalos;j++) data9[j]=0;
            dataset.write(data9,PredType::STD_I32LE);

            delete[] data9;

            // Number of particles
            if (opt.outdataformat>=DATAOUTMERITNPART) {
                dataspace = DataSpace(rank,dims);
                datasetname=H5std_string("Npart");

                // Check if there are halos to output so it can be compressed
                if (chunk_dims[0]>0) {
                    // Modify dataset creation property to enable chunking
                    hdfdatasetproplist.setChunk(rank, chunk_dims);
                    // Set ZLIB (DEFLATE) Compression using level 6.
                    hdfdatasetproplist.setDeflate(6);
                    dataset = Fhdf.createDataSet(datasetname,PredType::STD_U64LE,dataspace,hdfdatasetproplist);
                }
                else {
                    dataset = Fhdf.createDataSet(datasetname,PredType::STD_U64LE,dataspace);
                }

                //Write out the dataset
                long unsigned *data10 = new long unsigned[h[0].numhalos];
                for(Int_t j=0;j<h[0].numhalos;j++)  data10[j]=h[0].Halo[j].NumberofParticles;
                dataset.write(data10,PredType::STD_U64LE);

                delete[] data10;
            }

            Fhdf.close();

            cout<<ThisTask<<" Done writing to "<<fname<<" "<<MyGetTime()-time1<<endl;

        }
    } //End of hdf5 output
#endif
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
    hsize_t dims[1], chunk_dims[1];
    int rank;
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
                        if (opt.outdataformat>=DATAOUTMERITNPART) {
                            Fout<<"\t"<<h[i].Halo[j].NumberofParticles;
                        }
                        Fout<<endl;
                        for (int k=0;k<p[i][j].NumberofDescendants;k++) {
                            Fout<<h[i+p[i][j].istep].Halo[p[i][j].DescendantList[k]-1].haloID<<" ";
                            Fout<<p[i][j].dtoptype[k]<<" ";
                            if (opt.outdataformat>=DATAOUTMERIT) {
                                Fout<<p[i][j].Merit[k]<<" ";
                            }
                            if (opt.outdataformat>=DATAOUTMERITNPART) {
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
                if (opt.outdataformat>=DATAOUTMERITNPART) {
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
                    if (opt.outdataformat>=DATAOUTMERITNPART) {
                        Fout<<"\t"<<h[i].Halo[j].NumberofParticles;
                    }
                    Fout<<endl;
                    for (int k=0;k<p[i][j].NumberofDescendants;k++) {
                        Fout<<h[i+p[i][j].istep].Halo[p[i][j].DescendantList[k]-1].haloID<<" ";
                        Fout<<p[i][j].dtoptype[k]<<" ";
                        if (opt.outdataformat>=DATAOUTMERIT) {
                            Fout<<p[i][j].Merit[k]<<" ";
                        }
                        if (opt.outdataformat>=DATAOUTMERITNPART) {
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
                if (opt.outdataformat>=DATAOUTMERITNPART) {
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
        //If hdf5 then write a tree file per snapshot meaning it can be written out in parallel
        iend=EndSnap-opt.numsteps;
        istart=StartSnap;
        if (ThisTask==0) istart=0;
        if (ThisTask==NProcs-1) iend=opt.numsnapshots-1;

        for (int i=istart;i<iend;i++) {

            sprintf(fname,"%s/snapshot_%03d.VELOCIraptor.tree",opt.outname,i+opt.snapshotvaloffset);
            cout<<ThisTask<<" is writing to "<<fname<<endl;

            //Header information
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


            //Setup the datasets
            dims[0]=h[i].numhalos;
            rank=1;
            // Set the minmum chunk size
            chunk_dims[0]=min((unsigned long)HDFOUTPUTCHUNKSIZE,(unsigned long)h[i].numhalos);
            //Create dataset proplist
            DSetCreatPropList hdfdatasetproplist1;

            //ID
            dataspace = DataSpace(rank,dims);
            datasetname=H5std_string("ID");

            // Check if there are halos to output so it can be compressed
            if (chunk_dims[0]>0) {
                // Modify dataset creation property to enable chunking
                hdfdatasetproplist1.setChunk(rank, chunk_dims);
                // Set ZLIB (DEFLATE) Compression using level 6.
                hdfdatasetproplist1.setDeflate(6);
                dataset = Fhdf.createDataSet(datasetname,PredType::STD_U64LE,dataspace,hdfdatasetproplist1);
            }
            else {
                dataset = Fhdf.createDataSet(datasetname,PredType::STD_U64LE,dataspace);
            }

            //Write out the dataset
            long unsigned *data1 = new long unsigned[h[i].numhalos];
            for(Int_t j=0;j<h[i].numhalos;j++) data1[j]=h[i].Halo[j].haloID;
            dataset.write(data1,PredType::STD_U64LE);

            delete[] data1;

            //Num Descendants
            dataspace = DataSpace(rank,dims);
            datasetname=H5std_string("NumDesc");

            // Check if there are halos to output so it can be compressed
            if (chunk_dims[0]>0) {
                // Modify dataset creation property to enable chunking
                hdfdatasetproplist1.setChunk(rank, chunk_dims);
                // Set ZLIB (DEFLATE) Compression using level 6.
                hdfdatasetproplist1.setDeflate(6);
                dataset = Fhdf.createDataSet(datasetname,PredType::STD_I32LE,dataspace,hdfdatasetproplist1);
            }
            else {
                dataset = Fhdf.createDataSet(datasetname,PredType::STD_I32LE,dataspace);
            }

            //Write out the dataset
            int *data2 = new int[h[i].numhalos];
            for(Int_t j=0;j<h[i].numhalos;j++) data2[j]=p[i][j].NumberofDescendants;
            dataset.write(data2,PredType::STD_I32LE);

            delete[] data2;


            // If want to output the number of particles
            if (opt.outdataformat>=DATAOUTMERITNPART) {
                dataspace = DataSpace(rank,dims);
                datasetname=H5std_string("Npart");

                // Check if there are halos to output so it can be compressed
                if (chunk_dims[0]>0) {
                    // Modify dataset creation property to enable chunking
                    hdfdatasetproplist1.setChunk(rank, chunk_dims);
                    // Set ZLIB (DEFLATE) Compression using level 6.
                    hdfdatasetproplist1.setDeflate(6);
                    dataset = Fhdf.createDataSet(datasetname,PredType::STD_U64LE,dataspace,hdfdatasetproplist1);
                }
                else {
                    dataset = Fhdf.createDataSet(datasetname,PredType::STD_U64LE,dataspace);
                }

                //Write out the dataset
                long unsigned *data3 = new long unsigned[h[i].numhalos];
                for(Int_t j=0;j<h[i].numhalos;j++)  data3[j]=h[i].Halo[j].NumberofParticles;
                dataset.write(data3,PredType::STD_U64LE);

                delete[] data3;
            }


            // Store the total number of decendants
            long unsigned totndesc=0;

            //Decendants offsets
            dataspace = DataSpace(rank,dims);
            datasetname=H5std_string("DescOffsets");

            // Check if there are halos to output so it can be compressed
            if (chunk_dims[0]>0) {
                // Modify dataset creation property to enable chunking
                hdfdatasetproplist1.setChunk(rank, chunk_dims);
                // Set ZLIB (DEFLATE) Compression using level 6.
                hdfdatasetproplist1.setDeflate(6);
                dataset = Fhdf.createDataSet(datasetname,PredType::STD_U64LE,dataspace,hdfdatasetproplist1);
            }
            else {
                dataset = Fhdf.createDataSet(datasetname,PredType::STD_U64LE,dataspace);
            }

            //Write out the dataset
            long unsigned *data4 = new long unsigned[h[i].numhalos];
            for(Int_t j=0;j<h[i].numhalos;j++){
                data4[j]=totndesc;
                totndesc+=p[i][j].NumberofDescendants;
            }
            dataset.write(data4,PredType::STD_U64LE);

            delete[] data4;

            //Set the new dataset parameters
            long unsigned itemp=0;
            dims[0]=totndesc;
            rank=1;
            // Set the minmum chunk size
            chunk_dims[0]=min((unsigned long)HDFOUTPUTCHUNKSIZE,totndesc);
            //Create dataset proplist
            DSetCreatPropList hdfdatasetproplist2;

            // Descendant IDs
            dataspace = DataSpace(rank,dims);
            datasetname=H5std_string("Descendants");

            // Check if there are halos to output so it can be compressed
            if (chunk_dims[0]>0) {
                // Modify dataset creation property to enable chunking
                hdfdatasetproplist2.setChunk(rank, chunk_dims);
                // Set ZLIB (DEFLATE) Compression using level 6.
                hdfdatasetproplist2.setDeflate(6);
                dataset = Fhdf.createDataSet(datasetname,PredType::STD_U64LE,dataspace,hdfdatasetproplist2);
            }
            else {
                dataset = Fhdf.createDataSet(datasetname,PredType::STD_U64LE,dataspace);
            }

            //Write out the dataset
            long unsigned *data5 = new long unsigned[totndesc];
            for(Int_t j=0;j<h[i].numhalos;j++)
                for (Int_t k=0;k<p[i][j].NumberofDescendants;k++)
                    data5[itemp++] = h[i+p[i][j].istep].Halo[p[i][j].DescendantList[k]-1].haloID;

            dataset.write(data5,PredType::STD_U64LE);

            delete[] data5;

            // Descendant rank
            itemp=0;
            dataspace = DataSpace(rank,dims);
            datasetname=H5std_string("Ranks");

            // Check if there are halos to output so it can be compressed
            if (chunk_dims[0]>0) {
                // Modify dataset creation property to enable chunking
                hdfdatasetproplist2.setChunk(rank, chunk_dims);
                // Set ZLIB (DEFLATE) Compression using level 6.
                hdfdatasetproplist2.setDeflate(6);
                dataset = Fhdf.createDataSet(datasetname,PredType::STD_I32LE,dataspace,hdfdatasetproplist2);
            }
            else {
                dataset = Fhdf.createDataSet(datasetname,PredType::STD_I32LE,dataspace);
            }

            //Write out the dataset
            int *data6 = new int[totndesc];
            for(Int_t j=0;j<h[i].numhalos;j++)
                for (Int_t k=0;k<p[i][j].NumberofDescendants;k++)
                    data6[itemp++] = p[i][j].dtoptype[k];
            dataset.write(data6,PredType::STD_I32LE);

            delete[] data6;

            //Merits
            if (opt.outdataformat>=DATAOUTMERIT) {
                itemp=0;
                datasetname=H5std_string("Merits");
                dataspace = DataSpace(rank,dims);

                // Check if there are halos to output so it can be compressed
                if (chunk_dims[0]>0) {
                    // Modify dataset creation property to enable chunking
                    hdfdatasetproplist2.setChunk(rank, chunk_dims);
                    // Set ZLIB (DEFLATE) Compression using level 6.
                    hdfdatasetproplist2.setDeflate(6);
                    dataset = Fhdf.createDataSet(datasetname,PredType::NATIVE_FLOAT,dataspace,hdfdatasetproplist2);
                }
                else {
                    dataset = Fhdf.createDataSet(datasetname,PredType::NATIVE_FLOAT,dataspace);
                }

                //Write out the dataset
                float *data7 = new float[totndesc];
                for(Int_t j=0;j<h[i].numhalos;j++)
                    for (Int_t k=0;k<p[i][j].NumberofDescendants;k++)
                        data7[itemp++] = p[i][j].Merit[k];
                dataset.write(data7,PredType::NATIVE_FLOAT);

                delete[] data7;
            }

            //Descendant number of particles
            if (opt.outdataformat>=DATAOUTMERITNPART) {
                itemp=0;
                dataspace = DataSpace(rank,dims);
                datasetname=H5std_string("DescNpart");

                // Check if there are halos to output so it can be compressed
                if (chunk_dims[0]>0) {
                    // Modify dataset creation property to enable chunking
                    hdfdatasetproplist2.setChunk(rank, chunk_dims);
                    // Set ZLIB (DEFLATE) Compression using level 6.
                    hdfdatasetproplist2.setDeflate(6);
                    dataset = Fhdf.createDataSet(datasetname,PredType::STD_U64LE,dataspace,hdfdatasetproplist2);
                }
                else {
                    dataset = Fhdf.createDataSet(datasetname,PredType::STD_U64LE,dataspace);
                }

                //Write out the dataset
                long unsigned *data8 = new long unsigned[totndesc];
                for(Int_t j=0;j<h[i].numhalos;j++)
                    for (Int_t k=0;k<p[i][j].NumberofDescendants;k++)
                        data8[itemp++] = h[i+p[i][j].istep].Halo[p[i][j].DescendantList[k]-1].NumberofParticles;
                dataset.write(data8,PredType::STD_U64LE);

                delete[] data8;
            }

            Fhdf.close();

            cout<<ThisTask<<" Done writing to "<<fname<<" "<<MyGetTime()-time1<<endl;

        }


        ///last file has no connections
        if(ThisTask==NProcs-1){

            sprintf(fname,"%s/snapshot_%03d.VELOCIraptor.tree",opt.outname,opt.numsnapshots-1+opt.snapshotvaloffset);
            cout<<ThisTask<<" is writing to "<<fname<<endl;

            //Header information
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




            ///last file has no connections
            dims[0]=h[opt.numsnapshots-1].numhalos;
            rank=1;
            // Set the minmum chunk size
            chunk_dims[0]=min((unsigned long)HDFOUTPUTCHUNKSIZE,(unsigned long)h[opt.numsnapshots-1].numhalos);
            //Create dataset proplist
            DSetCreatPropList hdfdatasetproplist;

            //ID
            dataspace = DataSpace(rank,dims);
            datasetname=H5std_string("ID");

            // Check if there are halos to output so it can be compressed
            if (chunk_dims[0]>0) {
                // Modify dataset creation property to enable chunking
                hdfdatasetproplist.setChunk(rank, chunk_dims);
                // Set ZLIB (DEFLATE) Compression using level 6.
                hdfdatasetproplist.setDeflate(6);
                dataset = Fhdf.createDataSet(datasetname,PredType::STD_U64LE,dataspace,hdfdatasetproplist);
            }
            else {
                dataset = Fhdf.createDataSet(datasetname,PredType::STD_U64LE,dataspace);
            }

            //Write out the dataset
            long unsigned *data9 = new long unsigned[h[opt.numsnapshots-1].numhalos];
            for(Int_t j=0;j<h[opt.numsnapshots-1].numhalos;j++) data9[j]=h[opt.numsnapshots-1].Halo[j].haloID;
            dataset.write(data9,PredType::STD_U64LE);

            delete[] data9;

            //Num Descendants
            dataspace = DataSpace(rank,dims);
            datasetname=H5std_string("NumDesc");

            // Check if there are halos to output so it can be compressed
            if (chunk_dims[0]>0) {
                // Modify dataset creation property to enable chunking
                hdfdatasetproplist.setChunk(rank, chunk_dims);
                // Set ZLIB (DEFLATE) Compression using level 6.
                hdfdatasetproplist.setDeflate(6);
                dataset = Fhdf.createDataSet(datasetname,PredType::STD_I32LE,dataspace,hdfdatasetproplist);
            }
            else {
                dataset = Fhdf.createDataSet(datasetname,PredType::STD_I32LE,dataspace);
            }

            //Write out the dataset
            int *data10 = new int[h[opt.numsnapshots-1].numhalos];
            for(Int_t j=0;j<h[opt.numsnapshots-1].numhalos;j++) data10[j]=0;
            dataset.write(data10,PredType::STD_I32LE);

            delete[] data10;

            // If want to output the number of particles
            if (opt.outdataformat>=DATAOUTMERITNPART) {
                dataspace = DataSpace(rank,dims);
                datasetname=H5std_string("Npart");

                // Check if there are halos to output so it can be compressed
                if (chunk_dims[0]>0) {
                    // Modify dataset creation property to enable chunking
                    hdfdatasetproplist.setChunk(rank, chunk_dims);
                    // Set ZLIB (DEFLATE) Compression using level 6.
                    hdfdatasetproplist.setDeflate(6);
                    dataset = Fhdf.createDataSet(datasetname,PredType::STD_U64LE,dataspace,hdfdatasetproplist);
                }
                else {
                    dataset = Fhdf.createDataSet(datasetname,PredType::STD_U64LE,dataspace);
                }

                //Write out the dataset
                long unsigned *data11 = new long unsigned[h[opt.numsnapshots-1].numhalos];
                for(Int_t j=0;j<h[opt.numsnapshots-1].numhalos;j++)  data11[j]=h[opt.numsnapshots-1].Halo[j].NumberofParticles;
                dataset.write(data11,PredType::STD_U64LE);

                delete[] data11;
            }

            Fhdf.close();


            cout<<ThisTask<<" Done writing to "<<fname<<" "<<MyGetTime()-time1<<endl;

        }
    }//End of hdf5 output
#endif
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
            if (opt.outdataformat>=DATAOUTMERITNPART) {
                Fout<<"\t"<<h[i].Halo[j].NumberofParticles;
            }
            Fout<<endl;
            if (i>0) {
                for (int k=0;k<p[i][j].NumberofProgenitors;k++) {
                    Fout<<h[i-p[i][j].istep].Halo[p[i][j].ProgenitorList[k]-1].haloID<<" ";
                    if (opt.outdataformat>=DATAOUTMERIT) {
                        Fout<<p[i][j].Merit[k]<<" ";
                    }
                    if (opt.outdataformat>=DATAOUTMERITNPART) {
                        Fout<<h[i-p[i][j].istep].Halo[p[i][j].ProgenitorList[k]-1].NumberofParticles<<" ";
                    }
                    Fout<<endl;
                }
            }
            if (i<opt.numsnapshots-1) {
                for (int k=0;k<d[i][j].NumberofDescendants;k++) {
                    Fout<<h[i+d[i][j].istep].Halo[d[i][j].DescendantList[k]-1].haloID<<" ";
                    Fout<<d[i][j].dtoptype[k]<<" ";
                    if (opt.outdataformat>=DATAOUTMERIT) {
                        Fout<<d[i][j].Merit[k]<<" ";
                    }
                    if (opt.outdataformat>=DATAOUTMERITNPART) {
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
            Fout<<h[i].Halo[j].haloID<<"\t"<<p[i][j].NumberofProgenitors;
            if (opt.outdataformat>=DATAOUTMERITNPART) {
                Fout<<"\t"<<h[i].Halo[j].NumberofParticles;
            }
            Fout<<endl;
            for (int k=0;k<p[i][j].NumberofProgenitors;k++) {
                Fout<<h[i-1].Halo[p[i][j].ProgenitorList[k]-1].haloID<<" ";
                if (opt.outdataformat>=DATAOUTMERIT) {
                    Fout<<p[i][j].Merit[k]<<" ";
                }
                if (opt.outdataformat>=DATAOUTMERITNPART) {
                    Fout<<h[i-p[i][j].istep].Halo[p[i][j].ProgenitorList[k]-1].NumberofParticles<<" ";
                }
                Fout<<endl;
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
    unsigned long chunksize,offset,oldchunksize;
    unsigned int nchunks;
#ifndef USEMPI
    int ThisTask=0, NProcs=1;
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
        chunksize=10000;
        nchunks=ceil(mapsize/(Double_t)chunksize);
        if (chunksize>mapsize) {
            chunksize=mapsize;
            nchunks=1;
        }
        oldchunksize=chunksize;
        keys=new IDTYPE[chunksize];
        offset=0;
        i=0;
        for (auto it: idmap) {
            keys[i++]=it.first;
            if (i==chunksize) {
                i=0;
                Fout.write((char*)keys,sizeof(IDTYPE)*chunksize);
                offset+=chunksize;
                chunksize=min(chunksize,mapsize-offset);
            }
        }
        delete[] keys;
        chunksize=oldchunksize;
        indices=new IDTYPE[chunksize];
        offset=0;
        i=0;
        for (auto it: idmap) {
            indices[i++]=it.second;
            if (i==chunksize) {
                i=0;
                Fout.write((char*)indices,sizeof(IDTYPE)*chunksize);
                offset+=chunksize;
                chunksize=min(chunksize,mapsize-offset);
            }
        }
        delete[] indices;
        Fout.close();
    }
}
///if a memory efficient particle id to index map file exists read it
int ReadPIDStoIndexMap(Options &opt,map<IDTYPE, IDTYPE>&idmap)
{
    char fname[1000];
    fstream Fin, Fin2;
    IDTYPE *keys,*indices;
    Int_t i;
    size_t idsize,mapsize;
    unsigned long chunksize,offset;
    unsigned int nchunks;
    int numsnap;
    long unsigned tothalo;
    int iflag=0;
    double time1;
#ifndef USEMPI
    int ThisTask=0, NProcs=1;
#endif
    if (ThisTask==0) {
        sprintf(fname,"%s.pidtoindexmap.dat",opt.outname);
        Fin.open(fname,ios::in|ios::binary);
        //used to get to desired offset in file to read indices
        Fin2.open(fname,ios::in|ios::binary);
        cout<<"Attempting to read unique memory efficent mapping for particle IDS to index to "<<fname<<endl;
        //write header information
        Fin.read((char*)&numsnap,sizeof(int));
        Fin.read((char*)&tothalo,sizeof(long unsigned));
        Fin.read((char*)&idsize,sizeof(size_t));
        //if all is well then keep reading
        if (opt.numsnapshots==numsnap && opt.TotalNumberofHalos && idsize==sizeof(IDTYPE)) {
            iflag=1;
        }
    }
#ifdef USEMPI
    MPI_Bcast(&iflag, 1, MPI_INT, 0, MPI_COMM_WORLD);
#endif
    if (iflag==0) {
        if (ThisTask==0) cout<<"Unable to read data, must generate map "<<endl;
        Fin.close();
        Fin2.close();
        return iflag;
    }
    time1=MyGetTime();
    idmap.clear();
    if (ThisTask==0) {
        cout<<"Reading information"<<endl;
        Fin.read((char*)&mapsize,sizeof(size_t));
    }
    //communicate information
#ifdef USEMPI
    MPI_Bcast(&mapsize, sizeof(size_t), MPI_BYTE, 0, MPI_COMM_WORLD);
#endif
    opt.MaxIDValue=mapsize;
    if (ThisTask==0) {
        cout<<"Map will have "<<opt.MaxIDValue<<" elements and need roughly "<<opt.MaxIDValue*(sizeof(IDTYPE)*2+3*4)/1024./1024./1024.<<"GB of mem "<<endl;
        //offset to start of indices
        Fin2.seekg(mapsize*sizeof(IDTYPE)+sizeof(int)+sizeof(long unsigned)+sizeof(size_t)+sizeof(size_t));
    }
    //and send information in chunks
    chunksize=floor(2147483648/((Double_t)NProcs*sizeof(IDTYPE)));
    nchunks=ceil(mapsize/(Double_t)chunksize);
    if (chunksize>mapsize) {
        chunksize=mapsize;
        nchunks=1;
    }
    keys=new IDTYPE[chunksize];
    indices=new IDTYPE[chunksize];
    offset=0;
    for (auto ichunk=0;ichunk<nchunks;ichunk++)
    {
        if (ThisTask==0) {
            Fin.read((char*)keys,sizeof(IDTYPE)*chunksize);
            Fin2.read((char*)indices,sizeof(IDTYPE)*chunksize);
        }
#ifdef USEMPI
        if (NProcs>1) {
        MPI_Bcast(keys, sizeof(IDTYPE)*chunksize, MPI_BYTE, 0, MPI_COMM_WORLD);
        MPI_Bcast(indices, sizeof(IDTYPE)*chunksize, MPI_BYTE, 0, MPI_COMM_WORLD);
        }
#endif
        for (i=0;i<chunksize;i++) idmap.insert(pair<IDTYPE, IDTYPE>(keys[i],indices[i]));
        offset+=chunksize;
        chunksize=min(chunksize,mapsize-offset);
    }
    if (ThisTask==0) {Fin.close(); Fin2.close();}
    if (ThisTask==0) cout<<"Took "<<MyGetTime()-time1<<endl;
    delete[] indices;
    delete[] keys;
    return iflag;
}
//@}
