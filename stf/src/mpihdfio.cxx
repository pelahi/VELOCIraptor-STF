/*! \file mpigadgetio.cxx
 *  \brief this file contains routines used with MPI compilation and gadget io and domain construction. 
 */

#if defined(USEMPI) && defined(USEHDF)

//-- For MPI

#include "stf.h"
#include "hdfitems.h"

/// \name HDF Domain decomposition
//@{

/*! 
    Determine the domain decomposition.\n
    Here the domains are constructured in data units
    only read tasks should call this routine. It is tricky to get appropriate load balancing and correct number of particles per processor.\n
    
    I could use recursive binary splitting like kd-tree along most spread axis till have appropriate number of volumes corresponding 
    to number of processors. Or build a Peno-Hilbert space filling curve. 

    BUT for identifying particles on another mpi domain, simple enough to split volume into rectangular cells.
*/

///Determine Domain for HDF input
void MPIDomainExtentHDF(Options &opt){
    char buf[2000];
    H5File Fhdf;
    HDF_Group_Names hdf_gnames;
    HDF_Header hdf_header_info;
    Group headergroup;
    Attribute headerattribs;
    DataSpace headerdataspace;
    float floatbuff;
    double doublebuff;
    FloatType floattype;
    PredType HDFREALTYPE(PredType::NATIVE_FLOAT);
    int ireaderror = 0;

    if (ThisTask==0) {
        if(opt.num_files>1) sprintf(buf,"%s.0.hdf5",opt.fname);
        else sprintf(buf,"%s.hdf5",opt.fname);

        //Try block to detect exceptions raised by any of the calls inside it
        try
        {
            //turn off the auto-printing when failure occurs so that we can
            //handle the errors appropriately
            Exception::dontPrint();

            //Open the specified file and the specified dataset in the file.
            Fhdf.openFile(buf, H5F_ACC_RDONLY);
            cout<<"Loading HDF header info in header group: "<<hdf_gnames.Header_name<<endl;
            //get header group 
            headergroup=Fhdf.openGroup(hdf_gnames.Header_name);

            //start reading attributes
            headerattribs=headergroup.openAttribute(hdf_header_info.names[hdf_header_info.IBoxSize]);
            headerdataspace=headerattribs.getSpace();
            floattype=headerattribs.getFloatType();
            if (floattype.getSize()==sizeof(float)) {
                headerattribs.read(PredType::NATIVE_FLOAT,&floatbuff);
                hdf_header_info.BoxSize=floatbuff;
            }
            if (floattype.getSize()==sizeof(double)) {
                headerattribs.read(PredType::NATIVE_DOUBLE,&doublebuff);
                hdf_header_info.BoxSize=doublebuff;
            }
        }
        catch(GroupIException error)
        {
            error.printError();
        }
        // catch failure caused by the H5File operations
        catch( FileIException error )
        {
            error.printError();

        }
        // catch failure caused by the DataSet operations
        catch( DataSetIException error )
        {
            error.printError();
            ireaderror=1;
        }
        // catch failure caused by the DataSpace operations
        catch( DataSpaceIException error )
        {
            error.printError();
            ireaderror=1;
        }
        // catch failure caused by the DataSpace operations
        catch( DataTypeIException error )
        {
            error.printError();
            ireaderror=1;
        }
        Fhdf.close();
        for (int i=0;i<3;i++) {mpi_xlim[i][0]=0;mpi_xlim[i][1]=hdf_header_info.BoxSize;}
    }
    //There may be issues with particles exactly on the edge of a domain so before expanded limits by a small amount
    //now only done if a specific compile option passed
#ifdef MPIEXPANDLIM
    for (int i=0;i<3;i++) {
        Double_t dx=0.001*(mpi_xlim[i][1]-mpi_xlim[i][0]);
        mpi_xlim[i][0]-=dx;mpi_xlim[i][1]+=dx;
    }
#endif

    //make sure limits have been found
    MPI_Barrier(MPI_COMM_WORLD);
    if (NProcs==1) {
        for (int i=0;i<3;i++) {
            mpi_domain[ThisTask].bnd[i][0]=mpi_xlim[i][0];
            mpi_domain[ThisTask].bnd[i][1]=mpi_xlim[i][1];
        }
    }
}

void MPIDomainDecompositionHDF(Options &opt){
    Int_t i,j,k,n,m;
    int Nsplit,isplit;

    if (ThisTask==0) {
    //determine the number of splits in each dimension
    Nsplit=log((float)NProcs)/log(2.0);
    mpi_ideltax[0]=0;mpi_ideltax[1]=1;mpi_ideltax[2]=2;
    isplit=0;
    for (j=0;j<3;j++) mpi_nxsplit[j]=0;
    for (j=0;j<Nsplit;j++) {
        mpi_nxsplit[mpi_ideltax[isplit++]]++;
        if (isplit==3) isplit=0;
    }
    for (j=0;j<3;j++) mpi_nxsplit[j]=pow(2.0,mpi_nxsplit[j]);
    //for all the cells along the boundary of axis with the third split axis (smallest variance if I actually tried load balancing with a KD-Tree)
    //set the domain limits to the sims limits
    int ix=mpi_ideltax[0],iy=mpi_ideltax[1],iz=mpi_ideltax[2];
    int mpitasknum;
    for (j=0;j<mpi_nxsplit[iy];j++) {
        for (i=0;i<mpi_nxsplit[ix];i++) {
            mpitasknum=i+j*mpi_nxsplit[ix]+0*(mpi_nxsplit[ix]*mpi_nxsplit[iy]);
            mpi_domain[mpitasknum].bnd[iz][0]=mpi_xlim[iz][0];
            mpitasknum=i+j*mpi_nxsplit[ix]+(mpi_nxsplit[iz]-1)*(mpi_nxsplit[ix]*mpi_nxsplit[iy]);
            mpi_domain[mpitasknum].bnd[iz][1]=mpi_xlim[iz][1];
        }
    }
    //here for domains along second axis
    for (k=0;k<mpi_nxsplit[iz];k++) {
        for (i=0;i<mpi_nxsplit[ix];i++) {
            mpitasknum=i+0*mpi_nxsplit[ix]+k*(mpi_nxsplit[ix]*mpi_nxsplit[iy]);
            mpi_domain[mpitasknum].bnd[iy][0]=mpi_xlim[iy][0];
            mpitasknum=i+(mpi_nxsplit[iy]-1)*mpi_nxsplit[ix]+k*(mpi_nxsplit[ix]*mpi_nxsplit[iy]);
            mpi_domain[mpitasknum].bnd[iy][1]=mpi_xlim[iy][1];
        }
    }
    //finally along axis with largest variance
    for (k=0;k<mpi_nxsplit[iz];k++) {
        for (j=0;j<mpi_nxsplit[iy];j++) {
            mpitasknum=0+j*mpi_nxsplit[ix]+k*(mpi_nxsplit[ix]*mpi_nxsplit[iy]);
            mpi_domain[mpitasknum].bnd[ix][0]=mpi_xlim[ix][0];
            mpitasknum=(mpi_nxsplit[ix]-1)+j*mpi_nxsplit[ix]+k*(mpi_nxsplit[ix]*mpi_nxsplit[iy]);
            mpi_domain[mpitasknum].bnd[ix][1]=mpi_xlim[ix][1];
        }
    }
    //here use the three different histograms to define the boundary
    int start[3],end[3];
    Double_t bndval[3],binsum[3],lastbin;
    start[0]=start[1]=start[2]=0;
    for (i=0;i<mpi_nxsplit[ix];i++) {
        bndval[0]=(mpi_xlim[ix][1]-mpi_xlim[ix][0])*(Double_t)(i+1)/(Double_t)mpi_nxsplit[ix];
        if(i<mpi_nxsplit[ix]-1) {
        for (j=0;j<mpi_nxsplit[iy];j++) {
            for (k=0;k<mpi_nxsplit[iz];k++) {
                //define upper limit 
                mpitasknum=i+j*mpi_nxsplit[ix]+k*(mpi_nxsplit[ix]*mpi_nxsplit[iy]);
                mpi_domain[mpitasknum].bnd[ix][1]=bndval[0];
                //define lower limit
                mpitasknum=(i+1)+j*mpi_nxsplit[ix]+k*(mpi_nxsplit[ix]*mpi_nxsplit[iy]);
                mpi_domain[mpitasknum].bnd[ix][0]=bndval[0];
            }
        }
        }
        //now for secondary splitting
        if (mpi_nxsplit[iy]>1) 
        for (j=0;j<mpi_nxsplit[iy];j++) {
            bndval[1]=(mpi_xlim[iy][1]-mpi_xlim[iy][0])*(Double_t)(j+1)/(Double_t)mpi_nxsplit[iy];
            if(j<mpi_nxsplit[iy]-1) {
            for (k=0;k<mpi_nxsplit[iz];k++) {
                mpitasknum=i+j*mpi_nxsplit[ix]+k*(mpi_nxsplit[ix]*mpi_nxsplit[iy]);
                mpi_domain[mpitasknum].bnd[iy][1]=bndval[1];
                mpitasknum=i+(j+1)*mpi_nxsplit[ix]+k*(mpi_nxsplit[ix]*mpi_nxsplit[iy]);
                mpi_domain[mpitasknum].bnd[iy][0]=bndval[1];
            }
            }
            if (mpi_nxsplit[iz]>1)
            for (k=0;k<mpi_nxsplit[iz];k++) {
                bndval[2]=(mpi_xlim[iz][1]-mpi_xlim[iz][0])*(Double_t)(k+1)/(Double_t)mpi_nxsplit[iz];
                if (k<mpi_nxsplit[iz]-1){
                mpitasknum=i+j*mpi_nxsplit[ix]+k*(mpi_nxsplit[ix]*mpi_nxsplit[iy]);
                mpi_domain[mpitasknum].bnd[iz][1]=bndval[2];
                mpitasknum=i+j*mpi_nxsplit[ix]+(k+1)*(mpi_nxsplit[ix]*mpi_nxsplit[iy]);
                mpi_domain[mpitasknum].bnd[iz][0]=bndval[2];
                }
            }
        }
    }
    }
}

///reads HDF file to determine number of particles in each MPIDomain
void MPINumInDomainHDF(Options &opt)
{
    if (NProcs>1) {
    MPIDomainExtentHDF(opt);
    MPIDomainDecompositionHDF(opt);
    MPIInitialDomainDecomposition();

    Int_t i,j,k,n,nchunk;
    char buf[2000];
    MPI_Status status;

    //structure stores the names of the groups in the hdf input
    HDF_Group_Names hdf_gnames;
    HDF_Part_Info hdf_gas_info(HDFGASTYPE);
    HDF_Part_Info hdf_dm_info(HDFDMTYPE);
    HDF_Part_Info hdf_tracer_info(HDFTRACERTYPE);
    HDF_Part_Info hdf_star_info(HDFSTARTYPE);
    HDF_Part_Info hdf_bh_info(HDFBHTYPE);

    HDF_Part_Info *hdf_parts[NHDFTYPE];
    hdf_parts[0]=&hdf_gas_info;
    hdf_parts[1]=&hdf_dm_info;
    //hdf_parts[2]=(void*)&hdf_extra_info;
    hdf_parts[3]=&hdf_tracer_info;
    hdf_parts[4]=&hdf_star_info;
    hdf_parts[5]=&hdf_bh_info;

    H5File *Fhdf;
    //to store the groups, data sets and their associated data spaces
    HDF_Header *hdf_header_info;
    Group *headergroup;
    Attribute *headerattribs;
    DataSpace *headerdataspace;
    Group *partsgroup;
    DataSet *partsdataset;
    DataSpace *partsdataspace;
    DataSpace chunkspace;
    //buffers to load data
    int *intbuff=new int[NHDFTYPE];
    long long *longbuff=new long long[NHDFTYPE];
    float *floatbuff=new float[HDFCHUNKSIZE*3];
    double *doublebuff=new double[HDFCHUNKSIZE*3];
    void *realbuff;
    //arrays to store number of items to read and offsets when selecting hyperslabs
    //at most one needs a dimensionality of 13 for the tracer particles in Illustris
    hsize_t filespacecount[13],filespaceoffset[13];
    //to determine types 
    IntType inttype;
    FloatType floattype;
    PredType HDFREALTYPE(PredType::NATIVE_FLOAT);
    PredType HDFINTEGERTYPE(PredType::NATIVE_LONG);
    int ifloat,iint;
    int datarank;
    hsize_t datadim[5];
    Int_t Nlocalbuf,ibuf=0,*Nbuf, *Nbaryonbuf;
    int *ireadfile,*ireadtask,*readtaskID;
    ireadtask=new int[NProcs];
    readtaskID=new int[opt.nsnapread];
    ireadfile=new int[opt.num_files];
    MPIDistributeReadTasks(opt,ireadtask,readtaskID);

    Nbuf=new Int_t[NProcs];
    Nbaryonbuf=new Int_t[NProcs];
    for (j=0;j<NProcs;j++) Nbuf[j]=0;
    for (j=0;j<NProcs;j++) Nbaryonbuf[j]=0;

        ///array listing number of particle types used.
    ///Since Illustris contains an unused type of particles (2) and tracer particles (3) really not useful to iterate over all particle types in loops
    int nusetypes,nbusetypes;
    int usetypes[NHDFTYPE];
    if (ireadtask[ThisTask]>=0) {
        if (opt.partsearchtype==PSTALL) {
            //lets assume there are dm/stars/gas.
            nusetypes=3;
            usetypes[0]=HDFGASTYPE;usetypes[1]=HDFDMTYPE;usetypes[2]=HDFSTARTYPE;
            //now if also blackholes/sink particles increase number of types
            if (opt.iusesinkparticles) usetypes[nusetypes++]=HDFBHTYPE;
            if (opt.iusewindparticles) usetypes[nusetypes++]=HDFWINDTYPE;
        }
        else if (opt.partsearchtype==PSTDARK) {nusetypes=1;usetypes[0]=HDFDMTYPE;}
        else if (opt.partsearchtype==PSTGAS) {nusetypes=1;usetypes[0]=HDFGASTYPE;}
        else if (opt.partsearchtype==PSTSTAR) {nusetypes=1;usetypes[0]=HDFSTARTYPE;}
        else if (opt.partsearchtype==PSTBH) {nusetypes=1;usetypes[0]=HDFBHTYPE;}

        Fhdf=new H5File[opt.num_files];
        hdf_header_info=new HDF_Header[opt.num_files];
        headergroup=new Group[opt.num_files];
        headerdataspace=new DataSpace[opt.num_files];
        headerattribs=new Attribute[opt.num_files];
        partsgroup=new Group[opt.num_files*NHDFTYPE];
        partsdataset=new DataSet[opt.num_files*NHDFTYPE];
        partsdataspace=new DataSpace[opt.num_files*NHDFTYPE];
        MPISetFilesRead(opt,ireadfile,ireadtask);
        for(i=0; i<opt.num_files; i++) {
            if(opt.num_files>1) sprintf(buf,"%s.%d.hdf5",opt.fname,i);
            else sprintf(buf,"%s.hdf5",opt.fname);
            //Open the specified file and the specified dataset in the file.
            Fhdf[i].openFile(buf, H5F_ACC_RDONLY);
            headergroup[i]=Fhdf[i].openGroup(hdf_gnames.Header_name);
            //get number in file
            headerattribs[i]=headergroup[i].openAttribute(hdf_header_info[i].names[hdf_header_info[i].INuminFile]);
            headerdataspace[i]=headerattribs[i].getSpace();
            inttype=headerattribs[i].getIntType();
            if (inttype.getSize()==sizeof(int)) {
                headerattribs[i].read(PredType::NATIVE_INT,&intbuff[0]);
                for (k=0;k<NHDFTYPE;k++) hdf_header_info[i].npart[k]=intbuff[k];
            }
            if (inttype.getSize()==sizeof(long long)) {
                headerattribs[i].read(PredType::NATIVE_LONG,&longbuff[0]);
                for (k=0;k<NHDFTYPE;k++) hdf_header_info[i].npart[k]=longbuff[k];
            }

            //open particle group structures 
            for (j=0;j<nusetypes;j++) {k=usetypes[j]; partsgroup[i*NHDFTYPE+k]=Fhdf[i].openGroup(hdf_gnames.part_names[k]);}
            if (opt.partsearchtype==PSTDARK && opt.iBaryonSearch) for (j=1;j<=nbusetypes;j++) {k=usetypes[j];partsgroup[i*NHDFTYPE+k]=Fhdf[i].openGroup(hdf_gnames.part_names[k]);}

            //get positions
            for (j=0;j<nusetypes;j++) {
                k=usetypes[j]; 
                partsdataset[i*NHDFTYPE+k]=partsgroup[i*NHDFTYPE+k].openDataSet(hdf_parts[k]->names[0]);
                partsdataspace[i*NHDFTYPE+k]=partsdataset[i*NHDFTYPE+k].getSpace();
                //assuming all particles use the same float type for shared property structures
                floattype=partsdataset[i*NHDFTYPE+k].getFloatType();
            }
            if (opt.partsearchtype==PSTDARK && opt.iBaryonSearch) for (j=1;j<=nbusetypes;j++) {
                k=usetypes[j]; 
                partsdataset[i*NHDFTYPE+k]=partsgroup[i*NHDFTYPE+k].openDataSet(hdf_parts[k]->names[0]);
                partsdataspace[i*NHDFTYPE+k]=partsdataset[i*NHDFTYPE+k].getSpace();
            }
            if (floattype.getSize()==sizeof(float)) {HDFREALTYPE=PredType::NATIVE_FLOAT;realbuff=floatbuff;ifloat=1;}
            else {HDFREALTYPE=PredType::NATIVE_DOUBLE ;realbuff=doublebuff;ifloat=0;}
            for (j=0;j<nusetypes;j++) {
                k=usetypes[j];
                //data loaded into memory in chunks
                if (hdf_header_info[i].npart[k]<HDFCHUNKSIZE)nchunk=hdf_header_info[i].npart[k];
                else nchunk=HDFCHUNKSIZE;
                for(n=0;n<hdf_header_info[i].npart[k];n+=nchunk)
                {
                    if (hdf_header_info[i].npart[k]-n<HDFCHUNKSIZE&&hdf_header_info[i].npart[k]-n>0)nchunk=hdf_header_info[i].npart[k]-n;
                    //setup hyperslab so that it is loaded into the buffer
                    datarank=1;
                    datadim[0]=nchunk*3;
                    chunkspace=DataSpace(datarank,datadim);
                    filespacecount[0]=nchunk;filespacecount[1]=3;
                    filespaceoffset[0]=n;filespaceoffset[1]=0;
                    partsdataspace[i*NHDFTYPE+k].selectHyperslab(H5S_SELECT_SET, filespacecount, filespaceoffset);
                    partsdataset[i*NHDFTYPE+k].read(realbuff,HDFREALTYPE,chunkspace,partsdataspace[i*NHDFTYPE+k]);

                    if (ifloat) {
                        for (int nn=0;nn<nchunk;nn++) {
                            ibuf=MPIGetParticlesProcessor(floatbuff[nn*3],floatbuff[nn*3+1],floatbuff[nn*3+2]);
                            Nbuf[ibuf]++;
                        }
                    }
                    else {
                        for (int nn=0;nn<nchunk;nn++) {
                            ibuf=MPIGetParticlesProcessor(doublebuff[nn*3],doublebuff[nn*3+1],doublebuff[nn*3+2]);
                            Nbuf[ibuf]++;
                        }
                    }
                }
            }
            if (opt.partsearchtype==PSTDARK && opt.iBaryonSearch) {
                for (j=1;j<=nbusetypes;j++) {
                    k=usetypes[j];
                    //data loaded into memory in chunks
                    if (hdf_header_info[i].npart[k]<HDFCHUNKSIZE)nchunk=hdf_header_info[i].npart[k];
                    else nchunk=HDFCHUNKSIZE;
                    for(n=0;n<hdf_header_info[i].npart[k];n+=nchunk)
                    {
                        if (hdf_header_info[i].npart[k]-n<HDFCHUNKSIZE&&hdf_header_info[i].npart[k]-n>0)nchunk=hdf_header_info[i].npart[k]-n;
                        //setup hyperslab so that it is loaded into the buffer
                        datarank=1;
                        datadim[0]=nchunk*3;
                        chunkspace=DataSpace(datarank,datadim);
                        filespacecount[0]=nchunk;filespacecount[1]=3;
                        filespaceoffset[0]=n;filespaceoffset[1]=0;
                        partsdataspace[i*NHDFTYPE+k].selectHyperslab(H5S_SELECT_SET, filespacecount, filespaceoffset);
                        partsdataset[i*NHDFTYPE+k].read(realbuff,HDFREALTYPE,chunkspace,partsdataspace[i*NHDFTYPE+k]);

                        if (ifloat) {
                            for (int nn=0;nn<nchunk;nn++) {
                                ibuf=MPIGetParticlesProcessor(floatbuff[nn*3],floatbuff[nn*3+1],floatbuff[nn*3+2]);
                                Nbaryonbuf[ibuf]++;
                            }
                        }
                        else {
                            for (int nn=0;nn<nchunk;nn++) {
                                ibuf=MPIGetParticlesProcessor(doublebuff[nn*3],doublebuff[nn*3+1],doublebuff[nn*3+2]);
                                Nbaryonbuf[ibuf]++;
                            }
                        }
                    }
                }
            }
            Fhdf[i].close();
        }
    }
    //now having read number of particles, run all gather
    Int_t mpi_nlocal[NProcs];
    MPI_Allreduce(Nbuf,mpi_nlocal,NProcs,MPI_Int_t,MPI_SUM,MPI_COMM_WORLD);
    Nlocal=mpi_nlocal[ThisTask];
    if (opt.iBaryonSearch) {
        MPI_Allreduce(Nbaryonbuf,mpi_nlocal,NProcs,MPI_Int_t,MPI_SUM,MPI_COMM_WORLD);
        Nlocalbaryon[0]=mpi_nlocal[ThisTask];
    }
    Nmemlocal=Nlocal*(1.0+MPIExportFac);
    if (opt.iBaryonSearch) Nmemlocalbaryon=Nlocalbaryon[0]*(1.0+MPIExportFac);
    }
}

//@}

#endif
