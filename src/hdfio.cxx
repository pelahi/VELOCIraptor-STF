/*! \file hdfio.cxx
 *  \brief this file contains routines for hdf snapshot file io

 Note that the code is partly based on HDF examples which have the following liscence

 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * Copyright by The HDF Group.                                               *
 * Copyright by the Board of Trustees of the University of Illinois.         *
 * All rights reserved.                                                      *
 *                                                                           *
 * This file is part of HDF5.  The full HDF5 copyright notice, including     *
 * terms governing use, modification, and redistribution, is contained in    *
 * the files COPYING and Copyright.html.  COPYING can be found at the root   *
 * of the source code distribution tree; Copyright.html can be found at the  *
 * root level of an installed copy of the electronic HDF5 document set and   *
 * is linked from the top-level documents page.  It can also be found at     *
 * hdfgroup.org/HDF5/doc/Copyright.html.  If you do not have                 *
 * access to either file, you may request a copy from help@hdfgroup.org.     *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

 */
#ifdef USEHDF

//-- HDF5 SPECIFIC IO

#include "stf.h"

#include "hdfitems.h"

extern "C" herr_t file_attrib_info(hid_t loc_id, const char *name, const H5L_info_t *linfo, void *opdata)
{
    hid_t attrib_id;
    //Open the group using its name via the c interface
    attrib_id = H5Aopen(loc_id, name, H5P_DEFAULT);
    //Display group name.
    cout<<"Attribute Name : " <<name<<" "<<attrib_id<<endl;
    //close the group via the c interface
    H5Aclose(attrib_id);
    return 0;
}

extern "C" herr_t file_data_info(hid_t loc_id, const char *name, const H5L_info_t *linfo, void *opdata)
{
    hid_t dataset_id;
    //Open the group using its name via the c interface
    dataset_id = H5Dopen2(loc_id, name, H5P_DEFAULT);
    //Display group name.
    cout<<"Data Name : " <<name<<" "<<dataset_id<<endl;
    //close the group via the c interface
    H5Dclose(dataset_id);
    return 0;
}

// operator function to interface with HDF c code and print out the Group names in the hdf file
extern "C" herr_t file_info(hid_t loc_id, const char *name, const H5L_info_t *linfo, void *opdata)
{
    hid_t group_id;
    //Open the group using its name via the c interface
    group_id = H5Gopen2(loc_id, name, H5P_DEFAULT);
    //Display group name.
    cout<<"Group Name : " <<name<<endl;
    //H5Literate(group_id, H5_INDEX_NAME, H5_ITER_INC, NULL, file_data_info, NULL);
    //close the group via the c interface
    H5Gclose(group_id);
    return 0;
}

///reads an hdf5 formatted file.
void ReadHDF(Options &opt, vector<Particle> &Part, const Int_t nbodies,Particle *&Pbaryons, Int_t nbaryons)
{
    //structure stores the names of the groups in the hdf input
    char buf[2000];
    HDF_Group_Names hdf_gnames (opt.ihdfnameconvention);
    //structures store names in groups
    HDF_Header *hdf_header_info;
    HDF_Part_Info hdf_gas_info(HDFGASTYPE,opt.ihdfnameconvention);
    HDF_Part_Info hdf_dm_info(HDFDMTYPE,opt.ihdfnameconvention);
    HDF_Part_Info hdf_tracer_info(HDFTRACERTYPE,opt.ihdfnameconvention);
    HDF_Part_Info hdf_star_info(HDFSTARTYPE,opt.ihdfnameconvention);
    HDF_Part_Info hdf_bh_info(HDFBHTYPE,opt.ihdfnameconvention);

    HDF_Part_Info *hdf_parts[NHDFTYPE];
    hdf_parts[0]=&hdf_gas_info;
    hdf_parts[1]=&hdf_dm_info;
    //hdf_parts[2]=(void*)&hdf_extra_info;
    hdf_parts[3]=&hdf_tracer_info;
    hdf_parts[4]=&hdf_star_info;
    hdf_parts[5]=&hdf_bh_info;

    H5File *Fhdf;
    //to store the groups, data sets and their associated data spaces
    Group *partsgroup;
    Attribute *headerattribs;
    DataSpace *headerdataspace;
    DataSet *partsdataset;
    DataSpace *partsdataspace;
    DataSpace chunkspace;
    int chunksize=opt.inputbufsize;
    //buffers to load data
    int *intbuff=new int[chunksize];
    long long *longbuff=new long long[chunksize];
    unsigned int *uintbuff=new unsigned int[chunksize];
    float *floatbuff=new float[chunksize*3];
    double *doublebuff=new double[chunksize*3];
    void *integerbuff,*realbuff;
    //arrays to store number of items to read and offsets when selecting hyperslabs
    hsize_t filespacecount[HDFMAXPROPDIM],filespaceoffset[HDFMAXPROPDIM];
    //to determine types
    IntType inttype;
    FloatType floattype;
    PredType HDFREALTYPE(PredType::NATIVE_FLOAT);
    PredType HDFINTEGERTYPE(PredType::NATIVE_LONG);
    int ifloat,ifloat_pos, iint;
    int datarank;
    hsize_t datadim[5];
    //if any conversion is need for metallicity
    float zmetconversion=1;
    if (opt.ihdfnameconvention == HDFILLUSTISNAMES) zmetconversion=ILLUSTRISZMET;

    ///array listing number of particle types used.
    ///Since Illustris contains an unused type of particles (2) and tracer particles (3) really not useful to iterate over all particle types in loops
    int nusetypes,nbusetypes;
    int usetypes[NHDFTYPE];
    if (opt.partsearchtype==PSTALL) {
        nusetypes=0;
        //assume existance of dark matter and gas
        usetypes[nusetypes++]=HDFGASTYPE;usetypes[nusetypes++]=HDFDMTYPE;
        if (opt.iuseextradarkparticles) {
            usetypes[nusetypes++]=HDFDM1TYPE;
            usetypes[nusetypes++]=HDFDM2TYPE;
    	}
        if (opt.iusestarparticles) usetypes[nusetypes++]=HDFSTARTYPE;
        if (opt.iusesinkparticles) usetypes[nusetypes++]=HDFBHTYPE;
        if (opt.iusewindparticles) usetypes[nusetypes++]=HDFWINDTYPE;
    }
    else if (opt.partsearchtype==PSTDARK) {
        nusetypes=1;usetypes[0]=HDFDMTYPE;
        if (opt.iuseextradarkparticles) {
            usetypes[nusetypes++]=HDFDM1TYPE;
            usetypes[nusetypes++]=HDFDM2TYPE;
        }
        if (opt.iBaryonSearch) {
            nbusetypes=1;usetypes[nusetypes+nbusetypes++]=HDFGASTYPE;
            if (opt.iusestarparticles) usetypes[nusetypes+nbusetypes++]=HDFSTARTYPE;
            if (opt.iusesinkparticles) usetypes[nusetypes+nbusetypes++]=HDFBHTYPE;
        }
    }
    else if (opt.partsearchtype==PSTGAS) {nusetypes=1;usetypes[0]=HDFGASTYPE;}
    else if (opt.partsearchtype==PSTSTAR) {nusetypes=1;usetypes[0]=HDFSTARTYPE;}
    else if (opt.partsearchtype==PSTBH) {
        nusetypes=1;usetypes[0]=HDFBHTYPE;
    }

    Int_t i,j,k,n,nchunk,count,bcount,itemp,count2,bcount2;

    //store cosmology
    double z,aadjust,Hubble,Hubbleflow;

    Double_t mscale,lscale,lvscale;
    Double_t MP_DM=MAXVALUE,LN,N_DM,MP_B=0;
    int ifirstfile=0,*ireadfile,ireaderror=0;
    int *ireadtask,*readtaskID;
    Int_t ninputoffset;

#ifdef USEMPI
    if (ThisTask == 0)
#endif
    opt.num_files = HDF_get_nfiles (opt.fname, opt.partsearchtype);

#ifndef USEMPI
    Int_t Ntotal;
    int ThisTask=0,NProcs=1;
    ireadfile=new int[opt.num_files];
    for (i=0;i<opt.num_files;i++) ireadfile[i]=1;
#else
    MPI_Bcast(&(opt.num_files), sizeof(opt.num_files), MPI_BYTE, 0, MPI_COMM_WORLD);
#endif

    //if verbose spit out the types of particles that are going to be searched for
    if (ThisTask==0 && opt.iverbose>1) {
        cout<<" --------------- "<<endl;
        cout<<"Expecting "<<nusetypes<<" types of particles to be read "<<endl;
        for (i=0;i<nusetypes;i++) cout<<"Particle "<<usetypes[i]<<" with name "<<hdf_gnames.part_names[usetypes[i]]<<endl;
        if (opt.partsearchtype==PSTDARK && opt.iBaryonSearch) {
            cout<<"Additionally, as full separate baryon search , expecting "<<nbusetypes<<" baryon particles"<<endl;
            for (i=1;i<=nbusetypes;i++) cout<<"Particle "<<usetypes[i]<<" with name "<<hdf_gnames.part_names[usetypes[i]]<<endl;
        }
    }

    //if MPI is used, read processors (all tasks with task numbers less than the number of snapshots) opens the file and loads the data into a particle buffer
    //this particle buffer is used to broadcast data to the appropriate processor
#ifdef USEMPI
    //since positions, velocities, masses are all at different points in the file,
    //to correctly assign particle to proccessor with correct velocities and mass must have several file pointers
    Particle *Pbuf;
    int mpi_ireaderror;

    //for parallel input
    MPI_Comm mpi_comm_read;
    vector<Particle> *Preadbuf;
    Int_t BufSize=opt.mpiparticlebufsize;
    Int_t *Nbuf, *Nreadbuf,*nreadoffset;
    int ibuf=0;
    Int_t ibufindex;
    Int_t *Nlocalthreadbuf;
    int *irecv, *mpi_irecvflag;
    MPI_Request *mpi_request;
    Int_t *mpi_nsend_baryon;
    if (opt.iBaryonSearch && opt.partsearchtype!=PSTALL) mpi_nsend_baryon=new Int_t[NProcs*NProcs];
    Int_t inreadsend,totreadsend;
    Int_t *mpi_nsend_readthread;
    Int_t *mpi_nsend_readthread_baryon;
    if (opt.iBaryonSearch) mpi_nsend_baryon=new Int_t[NProcs*NProcs];
    if (opt.nsnapread>1) {
        mpi_nsend_readthread=new Int_t[opt.nsnapread*opt.nsnapread];
        if (opt.iBaryonSearch) mpi_nsend_readthread_baryon=new Int_t[opt.nsnapread*opt.nsnapread];
    }

    //used in mpi to load access to all the data blocks of interest
    DataSet *partsdatasetall;
    DataSpace *partsdataspaceall;

    //extra blocks to store info
    float *velfloatbuff=new float[chunksize*3];
    double *veldoublebuff=new double[chunksize*3];
    float *massfloatbuff=new float[chunksize];
    double *massdoublebuff=new double[chunksize];
#ifdef GASON
    float *ufloatbuff=new float[chunksize];
    double *udoublebuff=new double[chunksize];
#endif
#if defined(GASON)&&defined(STARON)
    float *Zfloatbuff=new float[chunksize];
    double *Zdoublebuff=new double[chunksize];
    float *SFRfloatbuff=new float[chunksize];
    double *SFRdoublebuff=new double[chunksize];
#endif
#ifdef STARON
    float *Tagefloatbuff=new float[chunksize];
    double *Tagedoublebuff=new double[chunksize];
#endif
    Pbuf = NULL; /* Keep Pbuf NULL or allocated so we can check its status later */

    Nbuf=new Int_t[NProcs];
    for (int j=0;j<NProcs;j++) Nbuf[j]=0;
    nreadoffset=new Int_t[opt.nsnapread];
    ireadtask=new int[NProcs];
    readtaskID=new int[opt.nsnapread];
    MPIDistributeReadTasks(opt,ireadtask,readtaskID);
    MPI_Comm_split(MPI_COMM_WORLD, (ireadtask[ThisTask]>=0), ThisTask, &mpi_comm_read);

    if (ThisTask==0) cout<<"There are "<<opt.nsnapread<<" threads reading "<<opt.num_files<<" files "<<endl;
    if (ireadtask[ThisTask]>=0)
    {
        //to temporarily store data from gadget file
        Pbuf=new Particle[BufSize*NProcs];
        Nreadbuf=new Int_t[opt.num_files];
        for (int j=0;j<opt.num_files;j++) Nreadbuf[j]=0;
        if (opt.nsnapread>1){
            Preadbuf=new vector<Particle>[opt.nsnapread];
            for (int j=0;j<opt.nsnapread;j++) Preadbuf[j].reserve(BufSize);
        }
        //to determine which files the thread should read
        ireadfile=new int[opt.num_files];
        ifirstfile=MPISetFilesRead(opt,ireadfile,ireadtask);
        inreadsend=0;
        for (int j=0;j<opt.num_files;j++) inreadsend+=ireadfile[j];
        MPI_Allreduce(&inreadsend,&totreadsend,1,MPI_Int_t,MPI_MIN,mpi_comm_read);

    }
    else {
        Nlocalthreadbuf=new Int_t[opt.nsnapread];
        irecv=new int[opt.nsnapread];
        mpi_irecvflag=new int[opt.nsnapread];
        for (i=0;i<opt.nsnapread;i++) irecv[i]=1;
        mpi_request=new MPI_Request[opt.nsnapread];
    }
    Nlocal=0;
    if (opt.iBaryonSearch) Nlocalbaryon[0]=0;

    if (ireadtask[ThisTask]>=0) {
#endif
    //read the header
    Fhdf=new H5File[opt.num_files];
    hdf_header_info=new HDF_Header[opt.num_files];
    for (i=0; i<opt.num_files; i++) hdf_header_info[i] = HDF_Header(opt.ihdfnameconvention);
    headerdataspace=new DataSpace[opt.num_files];
    headerattribs=new Attribute[opt.num_files];
    partsgroup=new Group[opt.num_files*NHDFTYPE];
    partsdataset=new DataSet[opt.num_files*NHDFTYPE];
    partsdataspace=new DataSpace[opt.num_files*NHDFTYPE];
#ifdef USEMPI
    partsdatasetall=new DataSet[opt.num_files*NHDFTYPE*NHDFDATABLOCK];
    partsdataspaceall=new DataSpace[opt.num_files*NHDFTYPE*NHDFDATABLOCK];
#endif
    for(i=0; i<opt.num_files; i++) if(ireadfile[i]) {
        if(opt.num_files>1) sprintf(buf,"%s.%d.hdf5",opt.fname,(int)i);
        else sprintf(buf,"%s.hdf5",opt.fname);
        //Try block to detect exceptions raised by any of the calls inside it
        try
        {
            //turn off the auto-printing when failure occurs so that we can
            //handle the errors appropriately
            Exception::dontPrint();

            //Open the specified file and the specified dataset in the file.
            Fhdf[i].openFile(buf, H5F_ACC_RDONLY);
            if (ThisTask==0 && i==0) {
                cout<<buf<<endl;
                cout<<"HDF file contains the following group structures "<<endl;
                H5Literate(Fhdf[i].getId(), H5_INDEX_NAME, H5_ITER_INC, NULL, file_info, NULL);
                cout<<" Expecting "<<endl;
                for (j=0;j<NHDFTYPE+1;j++) cout<<hdf_gnames.names[j]<<endl;
            }

            //start reading attributes
            headerattribs[i]=get_attribute(Fhdf[i], hdf_header_info[i].names[hdf_header_info[i].IBoxSize]);
            headerdataspace[i]=headerattribs[i].getSpace();
            floattype=headerattribs[i].getFloatType();
            if (floattype.getSize()==sizeof(float)) {
                headerattribs[i].read(PredType::NATIVE_FLOAT,&floatbuff[0]);
                hdf_header_info[i].BoxSize=floatbuff[0];
            }
            if (floattype.getSize()==sizeof(double)) {
                headerattribs[i].read(PredType::NATIVE_DOUBLE,&doublebuff[0]);
                hdf_header_info[i].BoxSize=doublebuff[0];
            }

            headerattribs[i]=get_attribute(Fhdf[i], hdf_header_info[i].names[hdf_header_info[i].IMass]);
            headerdataspace[i]=headerattribs[i].getSpace();
            if (headerdataspace[i].getSimpleExtentNdims()!=1) ireaderror=1;
            floattype=headerattribs[i].getFloatType();
            if (floattype.getSize()==sizeof(float)) {
                headerattribs[i].read(PredType::NATIVE_FLOAT,&floatbuff[0]);
                for (k=0;k<NHDFTYPE;k++)hdf_header_info[i].mass[k]=floatbuff[k];
            }
            if (floattype.getSize()==sizeof(double)) {
                headerattribs[i].read(PredType::NATIVE_DOUBLE,&doublebuff[0]);
                for (k=0;k<NHDFTYPE;k++)hdf_header_info[i].mass[k]=doublebuff[k];
            }

            headerattribs[i]=get_attribute(Fhdf[i], hdf_header_info[i].names[hdf_header_info[i].INuminFile]);
            headerdataspace[i]=headerattribs[i].getSpace();
            if (headerdataspace[i].getSimpleExtentNdims()!=1) ireaderror=1;
            inttype=headerattribs[i].getIntType();
            if (inttype.getSize()==sizeof(int)) {
                headerattribs[i].read(PredType::NATIVE_INT,&intbuff[0]);
                for (k=0;k<NHDFTYPE;k++) hdf_header_info[i].npart[k]=intbuff[k];
            }
            if (inttype.getSize()==sizeof(long long)) {
                headerattribs[i].read(PredType::NATIVE_LONG,&longbuff[0]);
                for (k=0;k<NHDFTYPE;k++) hdf_header_info[i].npart[k]=longbuff[k];
            }

            headerattribs[i]=get_attribute(Fhdf[i], hdf_header_info[i].names[hdf_header_info[i].INumTot]);
            headerdataspace[i]=headerattribs[i].getSpace();
            headerattribs[i].read(PredType::NATIVE_UINT,&uintbuff[0]);
            for (k=0;k<NHDFTYPE;k++) hdf_header_info[i].npartTotal[k]=uintbuff[k];

            headerattribs[i]=get_attribute(Fhdf[i], hdf_header_info[i].names[hdf_header_info[i].INumTotHW]);
            headerdataspace[i]=headerattribs[i].getSpace();
            headerattribs[i].read(PredType::NATIVE_UINT,&uintbuff[0]);
            for (k=0;k<NHDFTYPE;k++) hdf_header_info[i].npartTotalHW[k]=uintbuff[k];

            headerattribs[i]=get_attribute(Fhdf[i], hdf_header_info[i].names[hdf_header_info[i].IOmega0]);
            headerdataspace[i]=headerattribs[i].getSpace();
            floattype=headerattribs[i].getFloatType();

            if (floattype.getSize()==sizeof(float)) {
                headerattribs[i].read(PredType::NATIVE_FLOAT,&floatbuff[0]);
                hdf_header_info[i].Omega0=floatbuff[0];
            }
            if (floattype.getSize()==sizeof(double)) {
                headerattribs[i].read(PredType::NATIVE_DOUBLE,&doublebuff[0]);
                hdf_header_info[i].Omega0=doublebuff[0];
            }

            headerattribs[i]=get_attribute(Fhdf[i], hdf_header_info[i].names[hdf_header_info[i].IOmegaL]);
            headerdataspace[i]=headerattribs[i].getSpace();
            floattype=headerattribs[i].getFloatType();
            if (floattype.getSize()==sizeof(float)) {
                headerattribs[i].read(PredType::NATIVE_FLOAT,&floatbuff[0]);
                hdf_header_info[i].OmegaLambda=floatbuff[0];
            }
            if (floattype.getSize()==sizeof(double)) {
                headerattribs[i].read(PredType::NATIVE_DOUBLE,&doublebuff[0]);
                hdf_header_info[i].OmegaLambda=doublebuff[0];
            }

            headerattribs[i]=get_attribute(Fhdf[i], hdf_header_info[i].names[hdf_header_info[i].IRedshift]);
            headerdataspace[i]=headerattribs[i].getSpace();
            floattype=headerattribs[i].getFloatType();
            if (floattype.getSize()==sizeof(float)) {
                headerattribs[i].read(PredType::NATIVE_FLOAT,&floatbuff[0]);
                hdf_header_info[i].redshift=floatbuff[0];
            }
            if (floattype.getSize()==sizeof(double)) {
                headerattribs[i].read(PredType::NATIVE_DOUBLE,&doublebuff[0]);
                hdf_header_info[i].redshift=doublebuff[0];
            }

            headerattribs[i]=get_attribute(Fhdf[i], hdf_header_info[i].names[hdf_header_info[i].ITime]);
            headerdataspace[i]=headerattribs[i].getSpace();
            floattype=headerattribs[i].getFloatType();

            if (floattype.getSize()==sizeof(float)) {
                headerattribs[i].read(PredType::NATIVE_FLOAT,&floatbuff[0]);
                hdf_header_info[i].time=floatbuff[0];
            }
            if (floattype.getSize()==sizeof(double)) {
                headerattribs[i].read(PredType::NATIVE_DOUBLE,&doublebuff[0]);
                hdf_header_info[i].time=doublebuff[0];
            }

            headerattribs[i]=get_attribute(Fhdf[i], hdf_header_info[i].names[hdf_header_info[i].IHubbleParam]);
            headerdataspace[i]=headerattribs[i].getSpace();
            floattype=headerattribs[i].getFloatType();

            if (floattype.getSize()==sizeof(float)) {
                headerattribs[i].read(PredType::NATIVE_FLOAT,&floatbuff[0]);
                hdf_header_info[i].HubbleParam=floatbuff[0];
            }
            if (floattype.getSize()==sizeof(double)) {
                headerattribs[i].read(PredType::NATIVE_DOUBLE,&doublebuff[0]);
                hdf_header_info[i].HubbleParam=doublebuff[0];
            }

            headerattribs[i]=get_attribute(Fhdf[i], hdf_header_info[i].names[hdf_header_info[i].INumFiles]);
            headerdataspace[i]=headerattribs[i].getSpace();
            inttype=headerattribs[i].getIntType();
            if (inttype.getSize()==sizeof(int)) {
                headerattribs[i].read(PredType::NATIVE_INT,&intbuff[0]);
                hdf_header_info[i].num_files=intbuff[0];
            }
            if (inttype.getSize()==sizeof(long long)) {
                headerattribs[i].read(PredType::NATIVE_LONG,&longbuff[0]);
                hdf_header_info[i].num_files=longbuff[0];
            }
        }
        catch(GroupIException &error)
        {
          HDF5PrintError(error);
          cerr<<"Error in group might suggest config file has the incorrect HDF naming convention. ";
          cerr<<"Check HDF_name_convetion or add new naming convention updating hdfitems.h in the source code. "<<endl;
          Fhdf[i].close();
#ifdef USEMPI
          MPI_Abort(MPI_COMM_WORLD,8);
#else
          exit(8);
#endif
        }
        // catch failure caused by the H5File operations
        catch( FileIException &error )
        {
          HDF5PrintError(error);
          cerr<<"Error reading file. Exiting "<<endl;
          Fhdf[i].close();
#ifdef USEMPI
          MPI_Abort(MPI_COMM_WORLD,8);
#else
          exit(8);
#endif
        }
        // catch failure caused by the DataSet operations
        catch( DataSetIException &error )
        {
          HDF5PrintError(error);
          cerr<<"Error in data set might suggest config file has the incorrect HDF naming convention. ";
          cerr<<"Check HDF_name_convetion or update hdfio.cxx in the source code to read correct format"<<endl;
          Fhdf[i].close();
#ifdef USEMPI
          MPI_Abort(MPI_COMM_WORLD,8);
#else
          exit(8);
#endif
        }
        // catch failure caused by the DataSpace operations
        catch( DataSpaceIException &error )
        {
          HDF5PrintError(error);
          cerr<<"Error in data space might suggest config file has the incorrect HDF naming convention. ";
          cerr<<"Check HDF_name_convetion or update hdfio.cxx in the source code to read correct format"<<endl;
          Fhdf[i].close();
#ifdef USEMPI
          MPI_Abort(MPI_COMM_WORLD,8);
#else
          exit(8);
#endif
        }
        // catch failure caused by the DataSpace operations
        catch( DataTypeIException &error )
        {
          HDF5PrintError(error);
          cerr<<"Error in data type might suggest need to update hdfio.cxx in the source code to read correct format"<<endl;
          Fhdf[i].close();
#ifdef USEMPI
          MPI_Abort(MPI_COMM_WORLD,8);
#else
          exit(8);
#endif
        }
    }
    //after info read, initialise cosmological parameters
    opt.p=hdf_header_info[ifirstfile].BoxSize;
    if (opt.icosmologicalin) {
        z=hdf_header_info[ifirstfile].redshift;
        opt.a=1./(1.+z);
        opt.Omega_m=hdf_header_info[ifirstfile].Omega0;
        opt.Omega_Lambda=hdf_header_info[ifirstfile].OmegaLambda;
        opt.h=hdf_header_info[ifirstfile].HubbleParam;
        opt.Omega_cdm=opt.Omega_m-opt.Omega_b;
        CalcOmegak(opt);
        //Hubble flow
        if (opt.comove) aadjust=1.0;
        else aadjust=opt.a;
        Hubble=GetHubble(opt, aadjust);
        CalcCriticalDensity(opt, aadjust);
        CalcBackgroundDensity(opt, aadjust);
        CalcVirBN98(opt,aadjust);
        //if opt.virlevel<0, then use virial overdensity based on Bryan and Norman 1997 virialization level is given by
        if (opt.virlevel<0) opt.virlevel=opt.virBN98;
        PrintCosmology(opt);
    }
    else {
      opt.a=1.0;
      aadjust=opt.a;
      Hubbleflow=0.;
      cout<<"Non-cosmological input, using h = "<< opt.h<<endl;
    }

    // SWIFT snapshots already include the 1/h factor factor,
    // so there is no need to include it.
    if(opt.ihdfnameconvention == HDFSWIFTEAGLENAMES) {
      mscale=opt.M;lscale=opt.L*aadjust;lvscale=opt.L*opt.a;
    }
    else {
      mscale=opt.M/opt.h;lscale=opt.L/opt.h*aadjust;lvscale=opt.L/opt.h*opt.a;
    }

    //ignore hubble flow
    Hubbleflow=0.;
    Ntotal=0;
    for (int j=0;j<NHDFTYPE;j++)
    {
      opt.numpart[j]=hdf_header_info[ifirstfile].npartTotal[j];
      Ntotal+=hdf_header_info[ifirstfile].npartTotal[j];
    }
    for (int j=0;j<NHDFTYPE;j++)
    {
      opt.numpart[j]+=((long long)hdf_header_info[ifirstfile].npartTotalHW[j]<<32);
      Ntotal+=((long long)hdf_header_info[ifirstfile].npartTotalHW[j]<<32);
    }
    if (ThisTask==0) {
      cout<<"File contains "<<Ntotal<<" particles and is at time "<<opt.a<<endl;
      cout<<"Particle system contains "<<nbodies<<" particles and is at time "<<opt.a<<" in a box of size "<<opt.p<<endl;
      cout<<"Cosmology (h,Omega_m,Omega_cdm,Omega_b,Omega_L) = ("<< opt.h<<","<<opt.Omega_m<<","<<opt.Omega_cdm<<","<<opt.Omega_b<<","<<opt.Omega_Lambda<<")"<<endl;
    }
    //by default the interparticle spacing is determined using GDMTYPE
    //which is particle of type 1
    N_DM=hdf_header_info[ifirstfile].npartTotal[HDFDMTYPE];
    N_DM+=((long long)hdf_header_info[ifirstfile].npartTotalHW[HDFDMTYPE]<<32);
    LN=(opt.p*lscale/pow(N_DM,1.0/3.0));
#ifdef USEMPI
    }
#endif
    //after finished reading the header, start on the actual particle information

#ifndef USEMPI
    //init counters
    count2=bcount2=0;
    //start loding particle data
    for(i=0; i<opt.num_files; i++) {
        cout<<ThisTask<<" is reading file "<<i<<endl;
        ///\todo should be more rigorous with try/catch stuff
        try {
            //open particle group structures
            for (j=0;j<nusetypes;j++) {
              k=usetypes[j];
              if (ThisTask==0 && opt.iverbose>1) cout<<"Opening group "<<hdf_gnames.part_names[k]<<endl;
              partsgroup[i*NHDFTYPE+k]=Fhdf[i].openGroup(hdf_gnames.part_names[k]);
            }
            if (opt.partsearchtype==PSTDARK && opt.iBaryonSearch) for (j=1;j<=nbusetypes;j++) {
              k=usetypes[j];
              if (ThisTask==0 && opt.iverbose>1) cout<<"Opening group "<<hdf_gnames.part_names[k]<<endl;
              partsgroup[i*NHDFTYPE+k]=Fhdf[i].openGroup(hdf_gnames.part_names[k]);
            }
            itemp=0;
            //get positions
            for (j=0;j<nusetypes;j++) {
              k=usetypes[j];
              if (ThisTask==0 && opt.iverbose>1) cout<<"Opening group "<<hdf_gnames.part_names[k]<<": Data set "<<hdf_parts[k]->names[itemp]<<endl;
              partsdataset[i*NHDFTYPE+k]=partsgroup[i*NHDFTYPE+k].openDataSet(hdf_parts[k]->names[itemp]);
              partsdataspace[i*NHDFTYPE+k]=partsdataset[i*NHDFTYPE+k].getSpace();
              //assuming all particles use the same float type for shared property structures
              floattype=partsdataset[i*NHDFTYPE+k].getFloatType();
            }
            if (opt.partsearchtype==PSTDARK && opt.iBaryonSearch) for (j=1;j<=nbusetypes;j++) {
              k=usetypes[j];
              partsdataset[i*NHDFTYPE+k]=partsgroup[i*NHDFTYPE+k].openDataSet(hdf_parts[k]->names[itemp]);
              partsdataspace[i*NHDFTYPE+k]=partsdataset[i*NHDFTYPE+k].getSpace();
            }
            if (floattype.getSize()==sizeof(float)) {HDFREALTYPE=PredType::NATIVE_FLOAT;realbuff=floatbuff;ifloat=1;}
            else {HDFREALTYPE=PredType::NATIVE_DOUBLE ;realbuff=doublebuff;ifloat=0;}
            count=count2;
            bcount=bcount2;
            for (j=0;j<nusetypes;j++) {
              k=usetypes[j];
              //data loaded into memory in chunks
              if (hdf_header_info[i].npart[k]<chunksize)nchunk=hdf_header_info[i].npart[k];
              else nchunk=chunksize;
              for(n=0;n<hdf_header_info[i].npart[k];n+=nchunk)
              {
                if (hdf_header_info[i].npart[k]-n<chunksize&&hdf_header_info[i].npart[k]-n>0)nchunk=hdf_header_info[i].npart[k]-n;
                //setup hyperslab so that it is loaded into the buffer
                datarank=1;
                datadim[0]=nchunk*3;
                chunkspace=DataSpace(datarank,datadim);
                filespacecount[0]=nchunk;filespacecount[1]=3;
                filespaceoffset[0]=n;filespaceoffset[1]=0;
                partsdataspace[i*NHDFTYPE+k].selectHyperslab(H5S_SELECT_SET, filespacecount, filespaceoffset);
                partsdataset[i*NHDFTYPE+k].read(realbuff,HDFREALTYPE,chunkspace,partsdataspace[i*NHDFTYPE+k]);

                if (ifloat) for (int nn=0;nn<nchunk;nn++) Part[count++].SetPosition(floatbuff[nn*3],floatbuff[nn*3+1],floatbuff[nn*3+2]);
                else for (int nn=0;nn<nchunk;nn++) Part[count++].SetPosition(doublebuff[nn*3],doublebuff[nn*3+1],doublebuff[nn*3+2]);
              }
            }
            if (opt.partsearchtype==PSTDARK && opt.iBaryonSearch) {
              for (j=1;j<=nbusetypes;j++) {
                k=usetypes[j];
                //data loaded into memory in chunks
                if (hdf_header_info[i].npart[k]<chunksize)nchunk=hdf_header_info[i].npart[k];
                else nchunk=chunksize;
                for(n=0;n<hdf_header_info[i].npart[k];n+=nchunk)
                {
                  if (hdf_header_info[i].npart[k]-n<chunksize&&hdf_header_info[i].npart[k]-n>0)nchunk=hdf_header_info[i].npart[k]-n;
                  //setup hyperslab so that it is loaded into the buffer
                  datarank=1;
                  datadim[0]=nchunk*3;
                  chunkspace=DataSpace(datarank,datadim);
                  filespacecount[0]=nchunk;filespacecount[1]=3;
                  filespaceoffset[0]=n;filespaceoffset[1]=0;
                  partsdataspace[i*NHDFTYPE+k].selectHyperslab(H5S_SELECT_SET, filespacecount, filespaceoffset);
                  partsdataset[i*NHDFTYPE+k].read(realbuff,HDFREALTYPE,chunkspace,partsdataspace[i*NHDFTYPE+k]);

                  if (ifloat) for (int nn=0;nn<nchunk;nn++) Pbaryons[bcount++].SetPosition(floatbuff[nn*3],floatbuff[nn*3+1],floatbuff[nn*3+2]);
                  else for (int nn=0;nn<nchunk;nn++) Pbaryons[bcount++].SetPosition(doublebuff[nn*3],doublebuff[nn*3+1],doublebuff[nn*3+2]);
                }
              }
            }
            //get velocities
            itemp++;
            for (j=0;j<nusetypes;j++) {
              k=usetypes[j];
              if (ThisTask==0 && opt.iverbose>1) cout<<"Opening group "<<hdf_gnames.part_names[k]<<": Data set "<<hdf_parts[k]->names[itemp]<<endl;
              partsdataset[i*NHDFTYPE+k]=partsgroup[i*NHDFTYPE+k].openDataSet(hdf_parts[k]->names[itemp]);
              partsdataspace[i*NHDFTYPE+k]=partsdataset[i*NHDFTYPE+k].getSpace();
              //assuming all particles use the same float type for shared property structures
              floattype=partsdataset[i*NHDFTYPE+k].getFloatType();
            }
            if (opt.partsearchtype==PSTDARK && opt.iBaryonSearch) for (j=1;j<=nbusetypes;j++) {
              k=usetypes[j];
              partsdataset[i*NHDFTYPE+k]=partsgroup[i*NHDFTYPE+k].openDataSet(hdf_parts[k]->names[itemp]);
              partsdataspace[i*NHDFTYPE+k]=partsdataset[i*NHDFTYPE+k].getSpace();
            }
            if (floattype.getSize()==sizeof(float)) {HDFREALTYPE=PredType::NATIVE_FLOAT;realbuff=floatbuff;ifloat=1;}
            else {HDFREALTYPE=PredType::NATIVE_DOUBLE ;realbuff=doublebuff;ifloat=0;}
            count=count2;
            bcount=bcount2;
            for (j=0;j<nusetypes;j++) {
              k=usetypes[j];
              //data loaded into memory in chunks
              if (hdf_header_info[i].npart[k]<chunksize)nchunk=hdf_header_info[i].npart[k];
              else nchunk=chunksize;
              for(n=0;n<hdf_header_info[i].npart[k];n+=nchunk)
              {
                if (hdf_header_info[i].npart[k]-n<chunksize&&hdf_header_info[i].npart[k]-n>0)nchunk=hdf_header_info[i].npart[k]-n;
                //setup hyperslab so that it is loaded into the buffer
                datarank=1;
                datadim[0]=nchunk*3;
                chunkspace=DataSpace(datarank,datadim);
                filespacecount[0]=nchunk;filespacecount[1]=3;
                filespaceoffset[0]=n;filespaceoffset[1]=0;
                partsdataspace[i*NHDFTYPE+k].selectHyperslab(H5S_SELECT_SET, filespacecount, filespaceoffset);
                partsdataset[i*NHDFTYPE+k].read(realbuff,HDFREALTYPE,chunkspace,partsdataspace[i*NHDFTYPE+k]);

                if (ifloat) for (int nn=0;nn<nchunk;nn++) Part[count++].SetVelocity(floatbuff[nn*3],floatbuff[nn*3+1],floatbuff[nn*3+2]);
                else for (int nn=0;nn<nchunk;nn++) Part[count++].SetVelocity(doublebuff[nn*3],doublebuff[nn*3+1],doublebuff[nn*3+2]);
              }
            }
            if (opt.partsearchtype==PSTDARK && opt.iBaryonSearch) {
              for (j=1;j<=nbusetypes;j++) {
                k=usetypes[j];
                //data loaded into memory in chunks
                if (hdf_header_info[i].npart[k]<chunksize)nchunk=hdf_header_info[i].npart[k];
                else nchunk=chunksize;
                for(n=0;n<hdf_header_info[i].npart[k];n+=nchunk)
                {
                  if (hdf_header_info[i].npart[k]-n<chunksize&&hdf_header_info[i].npart[k]-n>0)nchunk=hdf_header_info[i].npart[k]-n;
                  //setup hyperslab so that it is loaded into the buffer
                  datarank=1;
                  datadim[0]=nchunk*3;
                  chunkspace=DataSpace(datarank,datadim);
                  filespacecount[0]=nchunk;filespacecount[1]=3;
                  filespaceoffset[0]=n;filespaceoffset[1]=0;
                  partsdataspace[i*NHDFTYPE+k].selectHyperslab(H5S_SELECT_SET, filespacecount, filespaceoffset);
                  partsdataset[i*NHDFTYPE+k].read(realbuff,HDFREALTYPE,chunkspace,partsdataspace[i*NHDFTYPE+k]);

                  if (ifloat) for (int nn=0;nn<nchunk;nn++) Pbaryons[bcount++].SetVelocity(floatbuff[nn*3],floatbuff[nn*3+1],floatbuff[nn*3+2]);
                  else for (int nn=0;nn<nchunk;nn++) Pbaryons[bcount++].SetVelocity(doublebuff[nn*3],doublebuff[nn*3+1],doublebuff[nn*3+2]);
                }
              }
            }
            //get ids
            itemp++;
            for (j=0;j<nusetypes;j++) {
              k=usetypes[j];
              if (ThisTask==0 && opt.iverbose>1) cout<<"Opening group "<<hdf_gnames.part_names[k]<<": Data set "<<hdf_parts[k]->names[itemp]<<endl;
              partsdataset[i*NHDFTYPE+k]=partsgroup[i*NHDFTYPE+k].openDataSet(hdf_parts[k]->names[itemp]);
              partsdataspace[i*NHDFTYPE+k]=partsdataset[i*NHDFTYPE+k].getSpace();
              //assuming all particles use the same float type for shared property structures
              inttype=partsdataset[i*NHDFTYPE+k].getIntType();
            }
            if (opt.partsearchtype==PSTDARK && opt.iBaryonSearch) for (j=1;j<=nbusetypes;j++) {
              k=usetypes[j];
              partsdataset[i*NHDFTYPE+k]=partsgroup[i*NHDFTYPE+k].openDataSet(hdf_parts[k]->names[itemp]);
              partsdataspace[i*NHDFTYPE+k]=partsdataset[i*NHDFTYPE+k].getSpace();
            }
            if (inttype.getSize()==sizeof(int)) {HDFINTEGERTYPE=PredType::NATIVE_INT;integerbuff=intbuff;iint=1;}
            else {HDFINTEGERTYPE=PredType::NATIVE_LONG;integerbuff=longbuff;iint=0;}
            count=count2;
            bcount=bcount2;
            for (j=0;j<nusetypes;j++) {
              k=usetypes[j];
              //data loaded into memory in chunks
              if (hdf_header_info[i].npart[k]<chunksize)nchunk=hdf_header_info[i].npart[k];
              else nchunk=chunksize;
              ninputoffset=0;
              for(n=0;n<hdf_header_info[i].npart[k];n+=nchunk)
              {
                if (hdf_header_info[i].npart[k]-n<chunksize&&hdf_header_info[i].npart[k]-n>0)nchunk=hdf_header_info[i].npart[k]-n;
                //setup hyperslab so that it is loaded into the buffer
                datarank=1;
                datadim[0]=nchunk;
                chunkspace=DataSpace(datarank,datadim);
                filespacecount[0]=nchunk;
                filespaceoffset[0]=n;
                partsdataspace[i*NHDFTYPE+k].selectHyperslab(H5S_SELECT_SET, filespacecount, filespaceoffset);
                partsdataset[i*NHDFTYPE+k].read(integerbuff,HDFINTEGERTYPE,chunkspace,partsdataspace[i*NHDFTYPE+k]);

                for (int nn=0;nn<nchunk;nn++) {
                    if (iint) Part[count].SetPID(intbuff[nn]);
                    else Part[count].SetPID(longbuff[nn]);
                    Part[count].SetID(count);
                    if (k==HDFGASTYPE) Part[count].SetType(GASTYPE);
                    else if (k==HDFDMTYPE) Part[count].SetType(DARKTYPE);
                    else if (k==HDFSTARTYPE) Part[count].SetType(STARTYPE);
                    else if (k==HDFBHTYPE) Part[count].SetType(BHTYPE);
#ifdef EXTRAINPUTINFO
                    if (opt.iextendedoutput)
                    {
                        Part[count].SetInputFileID(i);
                        Part[count].SetInputIndexInFile(nn+ninputoffset);
                    }
#endif
                    count++;
                }
                ninputoffset += nchunk;
              }
            }
            if (opt.partsearchtype==PSTDARK && opt.iBaryonSearch) {
              for (j=1;j<=nbusetypes;j++) {
                k=usetypes[j];
                //data loaded into memory in chunks
                if (hdf_header_info[i].npart[k]<chunksize)nchunk=hdf_header_info[i].npart[k];
                else nchunk=chunksize;
                ninputoffset=0;
                for(n=0;n<hdf_header_info[i].npart[k];n+=nchunk)
                {
                  if (hdf_header_info[i].npart[k]-n<chunksize&&hdf_header_info[i].npart[k]-n>0)nchunk=hdf_header_info[i].npart[k]-n;
                  datarank=1;
                  datadim[0]=nchunk;
                  chunkspace=DataSpace(datarank,datadim);
                  filespacecount[0]=nchunk;
                  filespaceoffset[0]=n;
                  partsdataspace[i*NHDFTYPE+k].selectHyperslab(H5S_SELECT_SET, filespacecount, filespaceoffset);
                  partsdataset[i*NHDFTYPE+k].read(integerbuff,HDFINTEGERTYPE,chunkspace,partsdataspace[i*NHDFTYPE+k]);

                  for (int nn=0;nn<nchunk;nn++) {
                    if (iint) Pbaryons[bcount].SetPID(intbuff[nn]);
                    else Pbaryons[bcount].SetPID(longbuff[nn]);
                    Pbaryons[bcount].SetID(bcount);
                    if (k==HDFGASTYPE) Pbaryons[bcount].SetType(GASTYPE);
                    else if (k==HDFSTARTYPE) Pbaryons[bcount].SetType(STARTYPE);
                    else if (k==HDFBHTYPE) Pbaryons[bcount].SetType(BHTYPE);
#ifdef EXTRAINPUTINFO
                    if (opt.iextendedoutput)
                    {
                        Pbaryons[bcount].SetInputFileID(i);
                        Pbaryons[bcount].SetInputIndexInFile(nn+ninputoffset);
                    }
#endif
                    bcount++;
                  }
                  ninputoffset+=nchunk;
                }
              }
            }

            //get masses, note that DM do not contain a mass field
            itemp++;
            for (j=0;j<nusetypes;j++) {
              k=usetypes[j];
              if (ThisTask==0 && opt.iverbose>1) cout<<"Opening group "<<hdf_gnames.part_names[k]<<": Data set "<<hdf_parts[k]->names[itemp]<<endl;
              if (hdf_header_info[i].mass[k]==0){
                partsdataset[i*NHDFTYPE+k]=partsgroup[i*NHDFTYPE+k].openDataSet(hdf_parts[k]->names[itemp]);
                partsdataspace[i*NHDFTYPE+k]=partsdataset[i*NHDFTYPE+k].getSpace();
                floattype=partsdataset[i*NHDFTYPE+k].getFloatType();
              }
            }
            if (opt.partsearchtype==PSTDARK && opt.iBaryonSearch) for (j=1;j<=nbusetypes;j++) {
              k=usetypes[j];
              if (hdf_header_info[i].mass[k]==0){
                partsdataset[i*NHDFTYPE+k]=partsgroup[i*NHDFTYPE+k].openDataSet(hdf_parts[k]->names[itemp]);
                partsdataspace[i*NHDFTYPE+k]=partsdataset[i*NHDFTYPE+k].getSpace();
                floattype=partsdataset[i*NHDFTYPE+k].getFloatType();
              }
            }
            if (floattype.getSize()==sizeof(float)) {HDFREALTYPE=PredType::NATIVE_FLOAT;realbuff=floatbuff;ifloat=1;}
            else {HDFREALTYPE=PredType::NATIVE_DOUBLE ;realbuff=doublebuff;ifloat=0;}
            count=count2;
            bcount=bcount2;
            for (j=0;j<nusetypes;j++) {
              k=usetypes[j];
              if (hdf_header_info[i].mass[k]==0) {
                //data loaded into memory in chunks
                if (hdf_header_info[i].npart[k]<chunksize)nchunk=hdf_header_info[i].npart[k];
                else nchunk=chunksize;
                for(n=0;n<hdf_header_info[i].npart[k];n+=nchunk)
                {
                  if (hdf_header_info[i].npart[k]-n<chunksize&&hdf_header_info[i].npart[k]-n>0)nchunk=hdf_header_info[i].npart[k]-n;
                  //setup hyperslab so that it is loaded into the buffer
                  datarank=1;
                  datadim[0]=nchunk;
                  chunkspace=DataSpace(datarank,datadim);
                  filespacecount[0]=nchunk;filespacecount[1]=1;
                  filespaceoffset[0]=n;filespaceoffset[1]=0;
                  partsdataspace[i*NHDFTYPE+k].selectHyperslab(H5S_SELECT_SET, filespacecount, filespaceoffset);
                  partsdataset[i*NHDFTYPE+k].read(realbuff,HDFREALTYPE,chunkspace,partsdataspace[i*NHDFTYPE+k]);

                  if (ifloat) for (int nn=0;nn<nchunk;nn++) Part[count++].SetMass(floatbuff[nn]);
                  else for (int nn=0;nn<nchunk;nn++) Part[count++].SetMass(doublebuff[nn]);
                }
              }
              else {
                for (int nn=0;nn<hdf_header_info[i].npart[k];nn++) Part[count++].SetMass(hdf_header_info[i].mass[k]);
              }
            }
            if (opt.partsearchtype==PSTDARK && opt.iBaryonSearch) {
              for (j=1;j<=nbusetypes;j++) {
                k=usetypes[j];
                if (hdf_header_info[i].mass[k]==0) {
                  //data loaded into memory in chunks
                  if (hdf_header_info[i].npart[k]<chunksize)nchunk=hdf_header_info[i].npart[k];
                  else nchunk=chunksize;
                  for(n=0;n<hdf_header_info[i].npart[k];n+=nchunk)
                  {
                    if (hdf_header_info[i].npart[k]-n<chunksize&&hdf_header_info[i].npart[k]-n>0)nchunk=hdf_header_info[i].npart[k]-n;
                    //setup hyperslab so that it is loaded into the buffer
                    datarank=1;
                    datadim[0]=nchunk;
                    chunkspace=DataSpace(datarank,datadim);
                    filespacecount[0]=nchunk;filespacecount[1]=1;
                    filespaceoffset[0]=n;filespaceoffset[1]=0;
                    partsdataspace[i*NHDFTYPE+k].selectHyperslab(H5S_SELECT_SET, filespacecount, filespaceoffset);
                    partsdataset[i*NHDFTYPE+k].read(realbuff,HDFREALTYPE,chunkspace,partsdataspace[i*NHDFTYPE+k]);

                    if (ifloat) for (int nn=0;nn<nchunk;nn++) Pbaryons[bcount++].SetMass(floatbuff[nn]);
                    else for (int nn=0;nn<nchunk;nn++) Pbaryons[bcount++].SetMass(doublebuff[nn]);
                  }
                }
                else {
                  for (int nn=0;nn<hdf_header_info[i].npart[k];nn++) Pbaryons[bcount++].SetMass(hdf_header_info[i].mass[k]);
                }
              }
            }

            //and if not just searching DM, load other parameters
            if (!(opt.partsearchtype==PSTDARK && opt.iBaryonSearch==0)) {
#ifdef GASON
              //first gas internal energy
              for (j=0;j<nusetypes;j++) {
                k=usetypes[j];
                if (k==HDFGASTYPE){
                  if (ThisTask==0 && opt.iverbose>1) cout<<"Opening group "<<hdf_gnames.part_names[k]<<": Data set "<<hdf_parts[k]->names[5]<<endl;
                  partsdataset[i*NHDFTYPE+k]=partsgroup[i*NHDFTYPE+k].openDataSet(hdf_parts[k]->names[5]);
                  partsdataspace[i*NHDFTYPE+k]=partsdataset[i*NHDFTYPE+k].getSpace();
                  floattype=partsdataset[i*NHDFTYPE+k].getFloatType();
                }
              }
              if (opt.partsearchtype==PSTDARK && opt.iBaryonSearch) for (j=1;j<=nbusetypes;j++) {
                k=usetypes[j];
                if (k==HDFGASTYPE){
                  partsdataset[i*NHDFTYPE+k]=partsgroup[i*NHDFTYPE+k].openDataSet(hdf_parts[k]->names[5]);
                  partsdataspace[i*NHDFTYPE+k]=partsdataset[i*NHDFTYPE+k].getSpace();
                  floattype=partsdataset[i*NHDFTYPE+k].getFloatType();
                }
              }
              count=count2;
              bcount=bcount2;
              for (j=0;j<nusetypes;j++) {
                k=usetypes[j];
                if (k==HDFGASTYPE) {
                  //data loaded into memory in chunks
                  if (hdf_header_info[i].npart[k]<chunksize)nchunk=hdf_header_info[i].npart[k];
                  else nchunk=chunksize;
                  for(n=0;n<hdf_header_info[i].npart[k];n+=nchunk)
                  {
                    if (hdf_header_info[i].npart[k]-n<chunksize&&hdf_header_info[i].npart[k]-n>0)nchunk=hdf_header_info[i].npart[k]-n;
                    //setup hyperslab so that it is loaded into the buffer
                    datarank=1;
                    datadim[0]=nchunk;
                    chunkspace=DataSpace(datarank,datadim);
                    filespacecount[0]=nchunk;filespacecount[1]=1;
                    filespaceoffset[0]=n;filespaceoffset[1]=0;
                    partsdataspace[i*NHDFTYPE+k].selectHyperslab(H5S_SELECT_SET, filespacecount, filespaceoffset);
                    partsdataset[i*NHDFTYPE+k].read(realbuff,HDFREALTYPE,chunkspace,partsdataspace[i*NHDFTYPE+k]);

                    if (ifloat) for (int nn=0;nn<nchunk;nn++) Part[count++].SetU(floatbuff[nn]);
                    else for (int nn=0;nn<nchunk;nn++) Part[count++].SetU(doublebuff[nn]);
                  }
                }
                else {
                  count+=hdf_header_info[i].npart[k];
                }
              }
              if (opt.partsearchtype==PSTDARK && opt.iBaryonSearch) {
                for (j=1;j<=nbusetypes;j++) {
                  k=usetypes[j];
                  if (k==HDFGASTYPE) {
                    //data loaded into memory in chunks
                    if (hdf_header_info[i].npart[k]<chunksize)nchunk=hdf_header_info[i].npart[k];
                    else nchunk=chunksize;
                    for(n=0;n<hdf_header_info[i].npart[k];n+=nchunk)
                    {
                      if (hdf_header_info[i].npart[k]-n<chunksize&&hdf_header_info[i].npart[k]-n>0)nchunk=hdf_header_info[i].npart[k]-n;
                      //setup hyperslab so that it is loaded into the buffer
                      datarank=1;
                      datadim[0]=nchunk;
                      chunkspace=DataSpace(datarank,datadim);
                      filespacecount[0]=nchunk;filespacecount[1]=1;
                      filespaceoffset[0]=n;filespaceoffset[1]=0;
                      partsdataspace[i*NHDFTYPE+k].selectHyperslab(H5S_SELECT_SET, filespacecount, filespaceoffset);
                      partsdataset[i*NHDFTYPE+k].read(realbuff,HDFREALTYPE,chunkspace,partsdataspace[i*NHDFTYPE+k]);

                      if (ifloat) for (int nn=0;nn<nchunk;nn++) Pbaryons[bcount++].SetU(floatbuff[nn]);
                      else for (int nn=0;nn<nchunk;nn++) Pbaryons[bcount++].SetU(doublebuff[nn]);
                    }
                  }
                  else {
                    count+=hdf_header_info[i].npart[k];
                  }
                }
              }
#ifdef STARON
              //if star forming get star formation rate
              for (j=0;j<nusetypes;j++) {
                k=usetypes[j];
                if (k==HDFGASTYPE){
                  if (ThisTask==0 && opt.iverbose>1) cout<<"Opening group "<<hdf_gnames.part_names[k]<<": Data set "<<hdf_parts[k]->names[6]<<endl;
                  partsdataset[i*NHDFTYPE+k]=partsgroup[i*NHDFTYPE+k].openDataSet(hdf_parts[k]->names[6]);
                  partsdataspace[i*NHDFTYPE+k]=partsdataset[i*NHDFTYPE+k].getSpace();
                  floattype=partsdataset[i*NHDFTYPE+k].getFloatType();
                }
              }
              if (opt.partsearchtype==PSTDARK && opt.iBaryonSearch) for (j=1;j<=nbusetypes;j++) {
                k=usetypes[j];
                if (k==HDFGASTYPE){
                  partsdataset[i*NHDFTYPE+k]=partsgroup[i*NHDFTYPE+k].openDataSet(hdf_parts[k]->names[6]);
                  partsdataspace[i*NHDFTYPE+k]=partsdataset[i*NHDFTYPE+k].getSpace();
                  floattype=partsdataset[i*NHDFTYPE+k].getFloatType();
                }
              }
              count=count2;
              bcount=bcount2;
              for (j=0;j<nusetypes;j++) {
                k=usetypes[j];
                if (k==HDFGASTYPE) {
                  //data loaded into memory in chunks
                  if (hdf_header_info[i].npart[k]<chunksize)nchunk=hdf_header_info[i].npart[k];
                  else nchunk=chunksize;
                  for(n=0;n<hdf_header_info[i].npart[k];n+=nchunk)
                  {
                    if (hdf_header_info[i].npart[k]-n<chunksize&&hdf_header_info[i].npart[k]-n>0)nchunk=hdf_header_info[i].npart[k]-n;
                    //setup hyperslab so that it is loaded into the buffer
                    datarank=1;
                    datadim[0]=nchunk;
                    chunkspace=DataSpace(datarank,datadim);
                    filespacecount[0]=nchunk;filespacecount[1]=1;
                    filespaceoffset[0]=n;filespaceoffset[1]=0;
                    partsdataspace[i*NHDFTYPE+k].selectHyperslab(H5S_SELECT_SET, filespacecount, filespaceoffset);
                    partsdataset[i*NHDFTYPE+k].read(realbuff,HDFREALTYPE,chunkspace,partsdataspace[i*NHDFTYPE+k]);

                    if (ifloat) for (int nn=0;nn<nchunk;nn++) Part[count++].SetSFR(floatbuff[nn]);
                    else for (int nn=0;nn<nchunk;nn++) Part[count++].SetSFR(doublebuff[nn]);
                  }
                }
                else {
                  count+=hdf_header_info[i].npart[k];
                }
              }
              if (opt.partsearchtype==PSTDARK && opt.iBaryonSearch) {
                for (j=1;j<=nbusetypes;j++) {
                  k=usetypes[j];
                  if (k==HDFGASTYPE) {
                    //data loaded into memory in chunks
                    if (hdf_header_info[i].npart[k]<chunksize)nchunk=hdf_header_info[i].npart[k];
                    else nchunk=chunksize;
                    for(n=0;n<hdf_header_info[i].npart[k];n+=nchunk)
                    {
                      if (hdf_header_info[i].npart[k]-n<chunksize&&hdf_header_info[i].npart[k]-n>0)nchunk=hdf_header_info[i].npart[k]-n;
                      //setup hyperslab so that it is loaded into the buffer
                      datarank=1;
                      datadim[0]=nchunk;
                      chunkspace=DataSpace(datarank,datadim);
                      filespacecount[0]=nchunk;filespacecount[1]=1;
                      filespaceoffset[0]=n;filespaceoffset[1]=0;
                      partsdataspace[i*NHDFTYPE+k].selectHyperslab(H5S_SELECT_SET, filespacecount, filespaceoffset);
                      partsdataset[i*NHDFTYPE+k].read(realbuff,HDFREALTYPE,chunkspace,partsdataspace[i*NHDFTYPE+k]);

                      if (ifloat) for (int nn=0;nn<nchunk;nn++) Pbaryons[bcount++].SetSFR(floatbuff[nn]);
                      else for (int nn=0;nn<nchunk;nn++) Pbaryons[bcount++].SetSFR(doublebuff[nn]);
                    }
                  }
                  else {
                    count+=hdf_header_info[i].npart[k];
                  }
                }
              }
              //then metallicity
              for (j=0;j<nusetypes;j++) {
                k=usetypes[j];
                if (k==HDFGASTYPE){
                  if (ThisTask==0 && opt.iverbose>1) cout<<"Opening group "<<hdf_gnames.part_names[k]<<": Data set "<<hdf_parts[k]->names[hdf_parts[k]->propindex[HDFGASIMETAL]]<<endl;
                  partsdataset[i*NHDFTYPE+k]=partsgroup[i*NHDFTYPE+k].openDataSet(hdf_parts[k]->names[hdf_parts[k]->propindex[HDFGASIMETAL]]);
                  partsdataspace[i*NHDFTYPE+k]=partsdataset[i*NHDFTYPE+k].getSpace();
                  floattype=partsdataset[i*NHDFTYPE+k].getFloatType();
                }
                if (k==HDFSTARTYPE){
                  if (ThisTask==0 && opt.iverbose>1) cout<<"Opening group "<<hdf_gnames.part_names[k]<<": Data set "<<hdf_parts[k]->names[hdf_parts[k]->propindex[HDFSTARIMETAL]]<<endl;
                  partsdataset[i*NHDFTYPE+k]=partsgroup[i*NHDFTYPE+k].openDataSet(hdf_parts[k]->names[hdf_parts[k]->propindex[HDFSTARIMETAL]]);
                  partsdataspace[i*NHDFTYPE+k]=partsdataset[i*NHDFTYPE+k].getSpace();
                  floattype=partsdataset[i*NHDFTYPE+k].getFloatType();
                }
              }
              if (opt.partsearchtype==PSTDARK && opt.iBaryonSearch) for (j=1;j<=nbusetypes;j++) {
                k=usetypes[j];
                if (k==HDFGASTYPE){
                  partsdataset[i*NHDFTYPE+k]=partsgroup[i*NHDFTYPE+k].openDataSet(hdf_parts[k]->names[hdf_parts[k]->propindex[HDFGASIMETAL]]);
                  partsdataspace[i*NHDFTYPE+k]=partsdataset[i*NHDFTYPE+k].getSpace();
                  floattype=partsdataset[i*NHDFTYPE+k].getFloatType();
                }
                if (k==HDFSTARTYPE){
                  partsdataset[i*NHDFTYPE+k]=partsgroup[i*NHDFTYPE+k].openDataSet(hdf_parts[k]->names[hdf_parts[k]->propindex[HDFSTARIMETAL]]);
                  partsdataspace[i*NHDFTYPE+k]=partsdataset[i*NHDFTYPE+k].getSpace();
                  floattype=partsdataset[i*NHDFTYPE+k].getFloatType();
                }
              }
              count=count2;
              bcount=bcount2;
              for (j=0;j<nusetypes;j++) {
                k=usetypes[j];
                if (k==HDFGASTYPE||k==HDFSTARTYPE) {
                  //data loaded into memory in chunks
                  if (hdf_header_info[i].npart[k]<chunksize)nchunk=hdf_header_info[i].npart[k];
                  else nchunk=chunksize;
                  for(n=0;n<hdf_header_info[i].npart[k];n+=nchunk)
                  {
                    if (hdf_header_info[i].npart[k]-n<chunksize&&hdf_header_info[i].npart[k]-n>0)nchunk=hdf_header_info[i].npart[k]-n;
                    //setup hyperslab so that it is loaded into the buffer
                    datarank=1;
                    datadim[0]=nchunk;
                    chunkspace=DataSpace(datarank,datadim);
                    filespacecount[0]=nchunk;filespacecount[1]=1;
                    filespaceoffset[0]=n;filespaceoffset[1]=0;
                    partsdataspace[i*NHDFTYPE+k].selectHyperslab(H5S_SELECT_SET, filespacecount, filespaceoffset);
                    partsdataset[i*NHDFTYPE+k].read(realbuff,HDFREALTYPE,chunkspace,partsdataspace[i*NHDFTYPE+k]);

                    if (ifloat) for (int nn=0;nn<nchunk;nn++) Part[count++].SetZmet(floatbuff[nn]*zmetconversion);
                    else for (int nn=0;nn<nchunk;nn++) Part[count++].SetZmet(doublebuff[nn]*zmetconversion);
                  }
                }
                else {
                  count+=hdf_header_info[i].npart[k];
                }
              }
              if (opt.partsearchtype==PSTDARK && opt.iBaryonSearch) {
                for (j=1;j<=nbusetypes;j++) {
                  k=usetypes[j];
                  if (k==HDFGASTYPE||k==HDFSTARTYPE) {
                    //data loaded into memory in chunks
                    if (hdf_header_info[i].npart[k]<chunksize)nchunk=hdf_header_info[i].npart[k];
                    else nchunk=chunksize;
                    for(n=0;n<hdf_header_info[i].npart[k];n+=nchunk)
                    {
                      if (hdf_header_info[i].npart[k]-n<chunksize&&hdf_header_info[i].npart[k]-n>0)nchunk=hdf_header_info[i].npart[k]-n;
                      //setup hyperslab so that it is loaded into the buffer
                      datarank=1;
                      datadim[0]=nchunk;
                      chunkspace=DataSpace(datarank,datadim);
                      filespacecount[0]=nchunk;filespacecount[1]=1;
                      filespaceoffset[0]=n;filespaceoffset[1]=0;
                      partsdataspace[i*NHDFTYPE+k].selectHyperslab(H5S_SELECT_SET, filespacecount, filespaceoffset);
                      partsdataset[i*NHDFTYPE+k].read(realbuff,HDFREALTYPE,chunkspace,partsdataspace[i*NHDFTYPE+k]);

                      if (ifloat) for (int nn=0;nn<nchunk;nn++) Pbaryons[bcount++].SetZmet(floatbuff[nn]*zmetconversion);
                      else for (int nn=0;nn<nchunk;nn++) Pbaryons[bcount++].SetZmet(doublebuff[nn]*zmetconversion);
                    }
                  }
                  else {
                    count+=hdf_header_info[i].npart[k];
                  }
                }
              }
              //then get star formation time, must also adjust so that if tage<0 this is a wind particle in Illustris so change particle type
              for (j=0;j<nusetypes;j++) {
                k=usetypes[j];
                if (k==HDFSTARTYPE){
                  if (ThisTask==0 && opt.iverbose>1) cout<<"Opening group "<<hdf_gnames.part_names[k]<<": Data set "<<hdf_parts[k]->names[hdf_parts[k]->propindex[HDFSTARIAGE]]<<endl;
                  partsdataset[i*NHDFTYPE+k]=partsgroup[i*NHDFTYPE+k].openDataSet(hdf_parts[k]->names[hdf_parts[k]->propindex[HDFSTARIAGE]]);
                  partsdataspace[i*NHDFTYPE+k]=partsdataset[i*NHDFTYPE+k].getSpace();
                  floattype=partsdataset[i*NHDFTYPE+k].getFloatType();
                }
              }
              if (opt.partsearchtype==PSTDARK && opt.iBaryonSearch) for (j=1;j<=nbusetypes;j++) {
                k=usetypes[j];
                if (k==HDFSTARTYPE){
                  partsdataset[i*NHDFTYPE+k]=partsgroup[i*NHDFTYPE+k].openDataSet(hdf_parts[k]->names[hdf_parts[k]->propindex[HDFSTARIAGE]]);
                  partsdataspace[i*NHDFTYPE+k]=partsdataset[i*NHDFTYPE+k].getSpace();
                  floattype=partsdataset[i*NHDFTYPE+k].getFloatType();
                }
              }
              count=count2;
              bcount=bcount2;
              for (j=0;j<nusetypes;j++) {
                k=usetypes[j];
                if (k==HDFSTARTYPE) {
                  //data loaded into memory in chunks
                  if (hdf_header_info[i].npart[k]<chunksize)nchunk=hdf_header_info[i].npart[k];
                  else nchunk=chunksize;
                  for(n=0;n<hdf_header_info[i].npart[k];n+=nchunk)
                  {
                    if (hdf_header_info[i].npart[k]-n<chunksize&&hdf_header_info[i].npart[k]-n>0)nchunk=hdf_header_info[i].npart[k]-n;
                    //setup hyperslab so that it is loaded into the buffer
                    datarank=1;
                    datadim[0]=nchunk;
                    chunkspace=DataSpace(datarank,datadim);
                    filespacecount[0]=nchunk;filespacecount[1]=1;
                    filespaceoffset[0]=n;filespaceoffset[1]=0;
                    partsdataspace[i*NHDFTYPE+k].selectHyperslab(H5S_SELECT_SET, filespacecount, filespaceoffset);
                    partsdataset[i*NHDFTYPE+k].read(realbuff,HDFREALTYPE,chunkspace,partsdataspace[i*NHDFTYPE+k]);

                    if (ifloat) for (int nn=0;nn<nchunk;nn++) {if (floatbuff[nn]<0) Part[count].SetType(WINDTYPE);Part[count++].SetTage(floatbuff[nn]);}
                    else for (int nn=0;nn<nchunk;nn++) {if (doublebuff[nn]<0) Part[count].SetType(WINDTYPE);Part[count++].SetTage(doublebuff[nn]);}
                  }
                }
                else {
                  count+=hdf_header_info[i].npart[k];
                }
              }
              if (opt.partsearchtype==PSTDARK && opt.iBaryonSearch) {
                for (j=1;j<=nbusetypes;j++) {
                  k=usetypes[j];
                  if (k==HDFGASTYPE||k==HDFSTARTYPE) {
                    //data loaded into memory in chunks
                    if (hdf_header_info[i].npart[k]<chunksize)nchunk=hdf_header_info[i].npart[k];
                    else nchunk=chunksize;
                    for(n=0;n<hdf_header_info[i].npart[k];n+=nchunk)
                    {
                      if (hdf_header_info[i].npart[k]-n<chunksize&&hdf_header_info[i].npart[k]-n>0)nchunk=hdf_header_info[i].npart[k]-n;
                      //setup hyperslab so that it is loaded into the buffer
                      datarank=1;
                      datadim[0]=nchunk;
                      chunkspace=DataSpace(datarank,datadim);
                      filespacecount[0]=nchunk;filespacecount[1]=1;
                      filespaceoffset[0]=n;filespaceoffset[1]=0;
                      partsdataspace[i*NHDFTYPE+k].selectHyperslab(H5S_SELECT_SET, filespacecount, filespaceoffset);
                      partsdataset[i*NHDFTYPE+k].read(realbuff,HDFREALTYPE,chunkspace,partsdataspace[i*NHDFTYPE+k]);

                      if (ifloat) for (int nn=0;nn<nchunk;nn++) {if (floatbuff[nn]<0) Pbaryons[bcount].SetType(WINDTYPE);Pbaryons[bcount++].SetTage(floatbuff[nn]);}
                      else for (int nn=0;nn<nchunk;nn++) {if (doublebuff[nn]<0) Pbaryons[bcount].SetType(WINDTYPE);Pbaryons[bcount++].SetTage(doublebuff[nn]);}
                    }
                  }
                  else {
                    count+=hdf_header_info[i].npart[k];
                  }
                }
              }
#endif
#endif
            }//end of if not dark matter then baryon search
            count2=count;
            bcount2=bcount;
        }//end of try
        catch(GroupIException error)
        {
            HDF5PrintError(error);
    		cerr<<"Error in group might suggest config file has the incorrect HDF naming convention. ";
    		cerr<<"Check HDF_name_convetion or add new naming convention updating hdfitems.h in the source code. "<<endl;
    		Fhdf[i].close();
    		exit(8);
    	}
        // catch failure caused by the H5File operations
        catch( FileIException error )
        {
            HDF5PrintError(error);
    		cerr<<"Error reading file. Exiting "<<endl;
    		Fhdf[i].close();
    		exit(8);
        }
        // catch failure caused by the DataSet operations
        catch( DataSetIException error )
        {
            HDF5PrintError(error);
    		cerr<<"Error in data set might suggest config file has the incorrect HDF naming convention. ";
    		cerr<<"Check HDF_name_convetion or update hdfio.cxx in the source code to read correct format"<<endl;
    		Fhdf[i].close();
    		exit(8);
        }
        // catch failure caused by the DataSpace operations
        catch( DataSpaceIException error )
        {
            HDF5PrintError(error);
    		cerr<<"Error in data space might suggest config file has the incorrect HDF naming convention. ";
    		cerr<<"Check HDF_name_convetion or update hdfio.cxx in the source code to read correct format"<<endl;
    		Fhdf[i].close();
    		exit(8);
        }
        // catch failure caused by the DataSpace operations
        catch( DataTypeIException error )
        {
            HDF5PrintError(error);
    		cerr<<"Error in data type might suggest need to update hdfio.cxx in the source code to read correct format"<<endl;
    		Fhdf[i].close();
    		exit(8);
        }
        Fhdf[i].close();
    }

    double vscale = 0.0;

    // SWIFT snapshot velocities already contain the sqrt(a) factor,
    // so there is no need to include it.
    if(opt.ihdfnameconvention == HDFSWIFTEAGLENAMES) vscale = opt.V;
    else vscale = opt.V*sqrt(opt.a);

    //finally adjust to appropriate units
    for (i=0;i<nbodies;i++)
    {
      Part[i].SetMass(Part[i].GetMass()*mscale);
      for (int j=0;j<3;j++) Part[i].SetVelocity(j,Part[i].GetVelocity(j)*vscale+Hubbleflow*Part[i].GetPosition(j));
      for (int j=0;j<3;j++) Part[i].SetPosition(j,Part[i].GetPosition(j)*lscale);
#ifdef GASON
      if (Part[i].GetType()==GASTYPE) Part[i].SetU(Part[i].GetU()*opt.V*opt.V);
#endif
    }
    if (Pbaryons!=NULL && opt.iBaryonSearch==1) {
      for (i=0;i<nbaryons;i++)
      {
        Pbaryons[i].SetMass(Pbaryons[i].GetMass()*mscale);
        for (int j=0;j<3;j++) Pbaryons[i].SetVelocity(j,Pbaryons[i].GetVelocity(j)*vscale+Hubbleflow*Pbaryons[i].GetPosition(j));
        for (int j=0;j<3;j++) Pbaryons[i].SetPosition(j,Pbaryons[i].GetPosition(j)*lscale);
#ifdef GASON
        Pbaryons[i].SetU(Pbaryons[i].GetU()*opt.V*opt.V);
#endif
      }
    }

#else
    //for all mpi threads that are reading input data, open file load access to data structures and begin loading into either local buffer or temporary buffer to be send to
    //non-read threads
    if (ireadtask[ThisTask]>=0) {
        inreadsend=0;
        count2=bcount2=0;
        for(i=0; i<opt.num_files; i++) if(ireadfile[i])
        {
            cout<<ThisTask<<" is reading file "<<i<<endl;
            ///\todo should be more rigorous with try/catch stuff
            try
            {
                //open particle group structures
                for (j=0;j<nusetypes;j++) {k=usetypes[j]; partsgroup[i*NHDFTYPE+k]=Fhdf[i].openGroup(hdf_gnames.part_names[k]);}
                if (opt.partsearchtype==PSTDARK && opt.iBaryonSearch) for (j=1;j<=nbusetypes;j++) {k=usetypes[j];partsgroup[i*NHDFTYPE+k]=Fhdf[i].openGroup(hdf_gnames.part_names[k]);}
                //open data structures that exist for all data blocks
                for (itemp=0;itemp<NHDFDATABLOCKALL;itemp++) {
                  //for everything but mass no header check needed.
                  if (itemp!=3) {
                    for (j=0;j<nusetypes;j++) {
                      k=usetypes[j];
                      if (ThisTask==0 && opt.iverbose>1) cout<<"Opening group "<<hdf_gnames.part_names[k]<<": Data set "<<hdf_parts[k]->names[itemp]<<endl;
                      partsdatasetall[i*NHDFTYPE*NHDFDATABLOCK+k*NHDFDATABLOCK+itemp]=partsgroup[i*NHDFTYPE+k].openDataSet(hdf_parts[k]->names[itemp]);
                      partsdataspaceall[i*NHDFTYPE*NHDFDATABLOCK+k*NHDFDATABLOCK+itemp]=partsdatasetall[i*NHDFTYPE*NHDFDATABLOCK+k*NHDFDATABLOCK+itemp].getSpace();
                    }
                    if (opt.partsearchtype==PSTDARK && opt.iBaryonSearch) for (j=1;j<=nbusetypes;j++) {
                      k=usetypes[j];
                      partsdatasetall[i*NHDFTYPE*NHDFDATABLOCK+k*NHDFDATABLOCK+itemp]=partsgroup[i*NHDFTYPE+k].openDataSet(hdf_parts[k]->names[itemp]);
                      partsdataspaceall[i*NHDFTYPE*NHDFDATABLOCK+k*NHDFDATABLOCK+itemp]=partsdatasetall[i*NHDFTYPE*NHDFDATABLOCK+k*NHDFDATABLOCK+itemp].getSpace();
                    }
                  }
                  else {
                    for (j=0;j<nusetypes;j++) {
                      k=usetypes[j];
                      if (hdf_header_info[i].mass[k]==0){
                        if (ThisTask==0 && opt.iverbose>1) cout<<"Opening group "<<hdf_gnames.part_names[k]<<": Data set "<<hdf_parts[k]->names[itemp]<<endl;
                        partsdatasetall[i*NHDFTYPE*NHDFDATABLOCK+k*NHDFDATABLOCK+itemp]=partsgroup[i*NHDFTYPE+k].openDataSet(hdf_parts[k]->names[itemp]);
                        partsdataspaceall[i*NHDFTYPE*NHDFDATABLOCK+k*NHDFDATABLOCK+itemp]=partsdatasetall[i*NHDFTYPE*NHDFDATABLOCK+k*NHDFDATABLOCK+itemp].getSpace();
                      }
                    }
                    if (opt.partsearchtype==PSTDARK && opt.iBaryonSearch) for (j=1;j<=nbusetypes;j++) {
                      k=usetypes[j];
                      partsdatasetall[i*NHDFTYPE*NHDFDATABLOCK+k*NHDFDATABLOCK+itemp]=partsgroup[i*NHDFTYPE+k].openDataSet(hdf_parts[k]->names[itemp]);
                      partsdataspaceall[i*NHDFTYPE*NHDFDATABLOCK+k*NHDFDATABLOCK+itemp]=partsdatasetall[i*NHDFTYPE*NHDFDATABLOCK+k*NHDFDATABLOCK+itemp].getSpace();
                    }
                  }
                }
                //now for extra data blocks
                //and if not just searching DM, load other parameters
                if (!(opt.partsearchtype==PSTDARK && opt.iBaryonSearch==0)) {
                  itemp=4;
#ifdef GASON
                  //first gas internal energy
                  for (j=0;j<nusetypes;j++) {
                    k=usetypes[j];
                    if (k==HDFGASTYPE){
                      if (ThisTask==0 && opt.iverbose>1) cout<<"Opening group "<<hdf_gnames.part_names[k]<<": Data set "<<hdf_parts[k]->names[5]<<endl;
                      partsdatasetall[i*NHDFTYPE*NHDFDATABLOCK+k*NHDFDATABLOCK+itemp]=partsgroup[i*NHDFTYPE+k].openDataSet(hdf_parts[k]->names[5]);
                      partsdataspaceall[i*NHDFTYPE*NHDFDATABLOCK+k*NHDFDATABLOCK+itemp]=partsdatasetall[i*NHDFTYPE*NHDFDATABLOCK+k*NHDFDATABLOCK+itemp].getSpace();
                    }
                  }
                  if (opt.partsearchtype==PSTDARK && opt.iBaryonSearch) for (j=1;j<=nbusetypes;j++) {
                    k=usetypes[j];
                    if (k==HDFGASTYPE){
                      partsdatasetall[i*NHDFTYPE*NHDFDATABLOCK+k*NHDFDATABLOCK+itemp]=partsgroup[i*NHDFTYPE+k].openDataSet(hdf_parts[k]->names[5]);
                      partsdataspaceall[i*NHDFTYPE*NHDFDATABLOCK+k*NHDFDATABLOCK+itemp]=partsdatasetall[i*NHDFTYPE*NHDFDATABLOCK+k*NHDFDATABLOCK+itemp].getSpace();
                    }
                  }
#ifdef STARON
                  //if star forming get star formation rate
                  itemp++;
                  for (j=0;j<nusetypes;j++) {
                    k=usetypes[j];
                    if (k==HDFGASTYPE){
                      if (ThisTask==0 && opt.iverbose>1) cout<<"Opening group "<<hdf_gnames.part_names[k]<<": Data set "<<hdf_parts[k]->names[6]<<endl;
                      partsdatasetall[i*NHDFTYPE*NHDFDATABLOCK+k*NHDFDATABLOCK+itemp]=partsgroup[i*NHDFTYPE+k].openDataSet(hdf_parts[k]->names[6]);
                      partsdataspaceall[i*NHDFTYPE*NHDFDATABLOCK+k*NHDFDATABLOCK+itemp]=partsdatasetall[i*NHDFTYPE*NHDFDATABLOCK+k*NHDFDATABLOCK+itemp].getSpace();
                    }
                  }
                  if (opt.partsearchtype==PSTDARK && opt.iBaryonSearch) for (j=1;j<=nbusetypes;j++) {
                    k=usetypes[j];
                    if (k==HDFGASTYPE){
                      partsdatasetall[i*NHDFTYPE*NHDFDATABLOCK+k*NHDFDATABLOCK+itemp]=partsgroup[i*NHDFTYPE+k].openDataSet(hdf_parts[k]->names[6]);
                      partsdataspaceall[i*NHDFTYPE*NHDFDATABLOCK+k*NHDFDATABLOCK+itemp]=partsdatasetall[i*NHDFTYPE*NHDFDATABLOCK+k*NHDFDATABLOCK+itemp].getSpace();
                    }
                  }
                  //then metallicity
                  itemp++;
                  for (j=0;j<nusetypes;j++) {
                    k=usetypes[j];
                    if (k==HDFGASTYPE){
                      if (ThisTask==0 && opt.iverbose>1) cout<<"Opening group "<<hdf_gnames.part_names[k]<<": Data set "<<hdf_parts[k]->names[hdf_parts[k]->propindex[HDFGASIMETAL]]<<endl;
                      partsdatasetall[i*NHDFTYPE*NHDFDATABLOCK+k*NHDFDATABLOCK+itemp]=partsgroup[i*NHDFTYPE+k].openDataSet(hdf_parts[k]->names[hdf_parts[k]->propindex[HDFGASIMETAL]]);
                      partsdataspaceall[i*NHDFTYPE*NHDFDATABLOCK+k*NHDFDATABLOCK+itemp]=partsdatasetall[i*NHDFTYPE*NHDFDATABLOCK+k*NHDFDATABLOCK+itemp].getSpace();
                    }
                    if (k==HDFSTARTYPE){
                      if (ThisTask==0 && opt.iverbose>1) cout<<"Opening group "<<hdf_gnames.part_names[k]<<": Data set "<<hdf_parts[k]->names[hdf_parts[k]->propindex[HDFSTARIMETAL]]<<endl;
                      partsdatasetall[i*NHDFTYPE*NHDFDATABLOCK+k*NHDFDATABLOCK+itemp]=partsgroup[i*NHDFTYPE+k].openDataSet(hdf_parts[k]->names[hdf_parts[k]->propindex[HDFSTARIMETAL]]);
                      partsdataspaceall[i*NHDFTYPE*NHDFDATABLOCK+k*NHDFDATABLOCK+itemp]=partsdatasetall[i*NHDFTYPE*NHDFDATABLOCK+k*NHDFDATABLOCK+itemp].getSpace();
                    }
                  }
                  if (opt.partsearchtype==PSTDARK && opt.iBaryonSearch) for (j=1;j<=nbusetypes;j++) {
                    k=usetypes[j];
                    if (k==HDFGASTYPE){
                      partsdatasetall[i*NHDFTYPE*NHDFDATABLOCK+k*NHDFDATABLOCK+itemp]=partsgroup[i*NHDFTYPE+k].openDataSet(hdf_parts[k]->names[hdf_parts[k]->propindex[HDFGASIMETAL]]);
                      partsdataspaceall[i*NHDFTYPE*NHDFDATABLOCK+k*NHDFDATABLOCK+itemp]=partsdatasetall[i*NHDFTYPE*NHDFDATABLOCK+k*NHDFDATABLOCK+itemp].getSpace();
                    }
                    if (k==HDFSTARTYPE){
                      partsdatasetall[i*NHDFTYPE*NHDFDATABLOCK+k*NHDFDATABLOCK+itemp]=partsgroup[i*NHDFTYPE+k].openDataSet(hdf_parts[k]->names[hdf_parts[k]->propindex[HDFSTARIMETAL]]);
                      partsdataspaceall[i*NHDFTYPE*NHDFDATABLOCK+k*NHDFDATABLOCK+itemp]=partsdatasetall[i*NHDFTYPE*NHDFDATABLOCK+k*NHDFDATABLOCK+itemp].getSpace();
                    }
                  }
                  //then get star formation time, must also adjust so that if tage<0 this is a wind particle in Illustris so change particle type
                  itemp++;
                  for (j=0;j<nusetypes;j++) {
                    k=usetypes[j];
                    if (k==HDFSTARTYPE){
                      if (ThisTask==0 && opt.iverbose>1) cout<<"Opening group "<<hdf_gnames.part_names[k]<<": Data set "<<hdf_parts[k]->names[hdf_parts[k]->propindex[HDFSTARIAGE]]<<endl;
                      partsdatasetall[i*NHDFTYPE*NHDFDATABLOCK+k*NHDFDATABLOCK+itemp]=partsgroup[i*NHDFTYPE+k].openDataSet(hdf_parts[k]->names[hdf_parts[k]->propindex[HDFSTARIAGE]]);
                      partsdataspaceall[i*NHDFTYPE*NHDFDATABLOCK+k*NHDFDATABLOCK+itemp]=partsdatasetall[i*NHDFTYPE*NHDFDATABLOCK+k*NHDFDATABLOCK+itemp].getSpace();
                    }
                  }
                  if (opt.partsearchtype==PSTDARK && opt.iBaryonSearch) for (j=1;j<=nbusetypes;j++) {
                    k=usetypes[j];
                    if (k==HDFSTARTYPE){
                      partsdatasetall[i*NHDFTYPE*NHDFDATABLOCK+k*NHDFDATABLOCK+itemp]=partsgroup[i*NHDFTYPE+k].openDataSet(hdf_parts[k]->names[hdf_parts[k]->propindex[HDFSTARIAGE]]);
                      partsdataspaceall[i*NHDFTYPE*NHDFDATABLOCK+k*NHDFDATABLOCK+itemp]=partsdatasetall[i*NHDFTYPE*NHDFDATABLOCK+k*NHDFDATABLOCK+itemp].getSpace();
                    }
                  }
#endif

#endif
                } //end of baryon read if not running search dm then baryons


                for (j=0;j<nusetypes;j++) {
                  k=usetypes[j];
                  //data loaded into memory in chunks
                  if (hdf_header_info[i].npart[k]<chunksize)nchunk=hdf_header_info[i].npart[k];
                  else nchunk=chunksize;
                  ninputoffset = 0;
                  for(n=0;n<hdf_header_info[i].npart[k];n+=nchunk)
                  {
                    if (hdf_header_info[i].npart[k]-n<chunksize&&hdf_header_info[i].npart[k]-n>0)nchunk=hdf_header_info[i].npart[k]-n;
                    //setup hyperslab so that it is loaded into the buffer
                    //load positions
                    itemp=0;
                    //set hyperslab
                    datarank=1;
                    datadim[0]=nchunk*3;
                    chunkspace=DataSpace(datarank,datadim);
                    filespacecount[0]=nchunk;filespacecount[1]=3;
                    filespaceoffset[0]=n;filespaceoffset[1]=0;
                    //set type
                    floattype=partsdatasetall[i*NHDFTYPE*NHDFDATABLOCK+k*NHDFDATABLOCK+itemp].getFloatType();
                    if (floattype.getSize()==sizeof(float)) {HDFREALTYPE=PredType::NATIVE_FLOAT;realbuff=floatbuff;ifloat_pos=1;}
                    else {HDFREALTYPE=PredType::NATIVE_DOUBLE ;realbuff=doublebuff;ifloat_pos=0;}
                    //read hyperslab into local buffer
                    partsdataspaceall[i*NHDFTYPE*NHDFDATABLOCK+k*NHDFDATABLOCK+itemp].selectHyperslab(H5S_SELECT_SET, filespacecount, filespaceoffset);
                    partsdatasetall[i*NHDFTYPE*NHDFDATABLOCK+k*NHDFDATABLOCK+itemp].read(realbuff,HDFREALTYPE,chunkspace,partsdataspaceall[i*NHDFTYPE*NHDFDATABLOCK+k*NHDFDATABLOCK+itemp]);
                    //velocities
                    itemp++;
                    datarank=1;
                    datadim[0]=nchunk*3;
                    chunkspace=DataSpace(datarank,datadim);
                    filespacecount[0]=nchunk;filespacecount[1]=3;
                    filespaceoffset[0]=n;filespaceoffset[1]=0;
                    floattype=partsdatasetall[i*NHDFTYPE*NHDFDATABLOCK+k*NHDFDATABLOCK+itemp].getFloatType();
                    if (floattype.getSize()==sizeof(float)) {HDFREALTYPE=PredType::NATIVE_FLOAT;realbuff=velfloatbuff;ifloat=1;}
                    else {HDFREALTYPE=PredType::NATIVE_DOUBLE ;realbuff=veldoublebuff;ifloat=0;}
                    partsdataspaceall[i*NHDFTYPE*NHDFDATABLOCK+k*NHDFDATABLOCK+itemp].selectHyperslab(H5S_SELECT_SET, filespacecount, filespaceoffset);
                    partsdatasetall[i*NHDFTYPE*NHDFDATABLOCK+k*NHDFDATABLOCK+itemp].read(realbuff,HDFREALTYPE,chunkspace,partsdataspaceall[i*NHDFTYPE*NHDFDATABLOCK+k*NHDFDATABLOCK+itemp]);
                    //ids
                    itemp++;
                    datarank=1;
                    datadim[0]=nchunk;
                    chunkspace=DataSpace(datarank,datadim);
                    filespacecount[0]=nchunk;filespacecount[1]=1;
                    filespaceoffset[0]=n;filespaceoffset[1]=0;
                    inttype=partsdatasetall[i*NHDFTYPE*NHDFDATABLOCK+k*NHDFDATABLOCK+itemp].getIntType();
                    if (inttype.getSize()==sizeof(int)) {HDFINTEGERTYPE=PredType::NATIVE_INT;integerbuff=intbuff;iint=1;}
                    else {HDFINTEGERTYPE=PredType::NATIVE_LONG;integerbuff=longbuff;iint=0;}
                    partsdataspaceall[i*NHDFTYPE*NHDFDATABLOCK+k*NHDFDATABLOCK+itemp].selectHyperslab(H5S_SELECT_SET, filespacecount, filespaceoffset);
                    partsdatasetall[i*NHDFTYPE*NHDFDATABLOCK+k*NHDFDATABLOCK+itemp].read(integerbuff,HDFINTEGERTYPE,chunkspace,partsdataspaceall[i*NHDFTYPE*NHDFDATABLOCK+k*NHDFDATABLOCK+itemp]);
                    //masses
                    itemp++;
                    datarank=1;
                    datadim[0]=nchunk;
                    chunkspace=DataSpace(datarank,datadim);
                    filespacecount[0]=nchunk;filespacecount[1]=1;
                    filespaceoffset[0]=n;filespaceoffset[1]=0;
                    if (hdf_header_info[i].mass[k]==0){
                      floattype=partsdatasetall[i*NHDFTYPE*NHDFDATABLOCK+k*NHDFDATABLOCK+itemp].getFloatType();
                      if (floattype.getSize()==sizeof(float)) {HDFREALTYPE=PredType::NATIVE_FLOAT;realbuff=massfloatbuff;ifloat=1;}
                      else {HDFREALTYPE=PredType::NATIVE_DOUBLE;realbuff=massdoublebuff;ifloat=0;}
                      partsdataspaceall[i*NHDFTYPE*NHDFDATABLOCK+k*NHDFDATABLOCK+itemp].selectHyperslab(H5S_SELECT_SET, filespacecount, filespaceoffset);
                      partsdatasetall[i*NHDFTYPE*NHDFDATABLOCK+k*NHDFDATABLOCK+itemp].read(realbuff,HDFREALTYPE,chunkspace,partsdataspaceall[i*NHDFTYPE*NHDFDATABLOCK+k*NHDFDATABLOCK+itemp]);
                    }
#ifdef GASON
                    //self-energy
                    itemp++;
                    datarank=1;
                    datadim[0]=nchunk;
                    chunkspace=DataSpace(datarank,datadim);
                    filespacecount[0]=nchunk;filespacecount[1]=1;
                    filespaceoffset[0]=n;filespaceoffset[1]=0;
                    if (k == HDFGASTYPE) {
                      floattype=partsdatasetall[i*NHDFTYPE*NHDFDATABLOCK+k*NHDFDATABLOCK+itemp].getFloatType();
                      if (floattype.getSize()==sizeof(float)) {HDFREALTYPE=PredType::NATIVE_FLOAT;realbuff=ufloatbuff;ifloat=1;}
                      else {HDFREALTYPE=PredType::NATIVE_DOUBLE;realbuff=udoublebuff;ifloat=0;}
                      partsdataspaceall[i*NHDFTYPE*NHDFDATABLOCK+k*NHDFDATABLOCK+itemp].selectHyperslab(H5S_SELECT_SET, filespacecount, filespaceoffset);
                      partsdatasetall[i*NHDFTYPE*NHDFDATABLOCK+k*NHDFDATABLOCK+itemp].read(realbuff,HDFREALTYPE,chunkspace,partsdataspaceall[i*NHDFTYPE*NHDFDATABLOCK+k*NHDFDATABLOCK+itemp]);
                    }
#ifdef STARON
                    //star formation rate
                    itemp++;
                    datarank=1;
                    datadim[0]=nchunk;
                    chunkspace=DataSpace(datarank,datadim);
                    filespacecount[0]=nchunk;filespacecount[1]=1;
                    filespaceoffset[0]=n;filespaceoffset[1]=0;
                    if (k == HDFGASTYPE) {
                      floattype=partsdatasetall[i*NHDFTYPE*NHDFDATABLOCK+k*NHDFDATABLOCK+itemp].getFloatType();
                      if (floattype.getSize()==sizeof(float)) {HDFREALTYPE=PredType::NATIVE_FLOAT;realbuff=SFRfloatbuff;ifloat=1;}
                      else {HDFREALTYPE=PredType::NATIVE_DOUBLE;realbuff=SFRdoublebuff;ifloat=0;}
                      partsdataspaceall[i*NHDFTYPE*NHDFDATABLOCK+k*NHDFDATABLOCK+itemp].selectHyperslab(H5S_SELECT_SET, filespacecount, filespaceoffset);
                      partsdatasetall[i*NHDFTYPE*NHDFDATABLOCK+k*NHDFDATABLOCK+itemp].read(realbuff,HDFREALTYPE,chunkspace,partsdataspaceall[i*NHDFTYPE*NHDFDATABLOCK+k*NHDFDATABLOCK+itemp]);
                    }

                    //metallicity
                    itemp++;
                    datarank=1;
                    datadim[0]=nchunk;
                    chunkspace=DataSpace(datarank,datadim);
                    filespacecount[0]=nchunk;filespacecount[1]=1;
                    filespaceoffset[0]=n;filespaceoffset[1]=0;
                    if (k == HDFGASTYPE || k == HDFSTARTYPE) {
                      floattype=partsdatasetall[i*NHDFTYPE*NHDFDATABLOCK+k*NHDFDATABLOCK+itemp].getFloatType();
                      if (floattype.getSize()==sizeof(float)) {HDFREALTYPE=PredType::NATIVE_FLOAT;realbuff=Zfloatbuff;ifloat=1;}
                      else {HDFREALTYPE=PredType::NATIVE_DOUBLE;realbuff=Zdoublebuff;ifloat=0;}
                      partsdataspaceall[i*NHDFTYPE*NHDFDATABLOCK+k*NHDFDATABLOCK+itemp].selectHyperslab(H5S_SELECT_SET, filespacecount, filespaceoffset);
                      partsdatasetall[i*NHDFTYPE*NHDFDATABLOCK+k*NHDFDATABLOCK+itemp].read(realbuff,HDFREALTYPE,chunkspace,partsdataspaceall[i*NHDFTYPE*NHDFDATABLOCK+k*NHDFDATABLOCK+itemp]);
                    }

                    //stellar age
                    itemp++;
                    datarank=1;
                    datadim[0]=nchunk;
                    chunkspace=DataSpace(datarank,datadim);
                    filespacecount[0]=nchunk;filespacecount[1]=1;
                    filespaceoffset[0]=n;filespaceoffset[1]=0;
                    if (k == HDFSTARTYPE) {
                      floattype=partsdatasetall[i*NHDFTYPE*NHDFDATABLOCK+k*NHDFDATABLOCK+itemp].getFloatType();
                      if (floattype.getSize()==sizeof(float)) {HDFREALTYPE=PredType::NATIVE_FLOAT;realbuff=Tagefloatbuff;ifloat=1;}
                      else {HDFREALTYPE=PredType::NATIVE_DOUBLE;realbuff=Tagedoublebuff;ifloat=0;}
                      partsdataspaceall[i*NHDFTYPE*NHDFDATABLOCK+k*NHDFDATABLOCK+itemp].selectHyperslab(H5S_SELECT_SET, filespacecount, filespaceoffset);
                      partsdatasetall[i*NHDFTYPE*NHDFDATABLOCK+k*NHDFDATABLOCK+itemp].read(realbuff,HDFREALTYPE,chunkspace,partsdataspaceall[i*NHDFTYPE*NHDFDATABLOCK+k*NHDFDATABLOCK+itemp]);
                    }
#endif
#endif

                    for (int nn=0;nn<nchunk;nn++) {
                        if (ifloat_pos) ibuf=MPIGetParticlesProcessor(floatbuff[nn*3],floatbuff[nn*3+1],floatbuff[nn*3+2]);
                        else ibuf=MPIGetParticlesProcessor(doublebuff[nn*3],doublebuff[nn*3+1],doublebuff[nn*3+2]);
                        ibufindex=ibuf*BufSize+Nbuf[ibuf];
                        //reset hydro quantities of buffer
#ifdef GASON
                        Pbuf[ibufindex].SetU(0);
#ifdef STARON
                        Pbuf[ibufindex].SetSFR(0);
                        Pbuf[ibufindex].SetZmet(0);
#endif
#endif
#ifdef STARON
                        Pbuf[ibufindex].SetZmet(0);
                        Pbuf[ibufindex].SetTage(0);
#endif
#ifdef BHON
#endif
                        //store particle info in Ptemp;
                        if (ifloat_pos)
                            Pbuf[ibufindex].SetPosition(floatbuff[nn*3],floatbuff[nn*3+1],floatbuff[nn*3+2]);
                        else
                            Pbuf[ibufindex].SetPosition(doublebuff[nn*3],doublebuff[nn*3+1],doublebuff[nn*3+2]);
                        if (ifloat) {
                            Pbuf[ibufindex].SetVelocity(velfloatbuff[nn*3],velfloatbuff[nn*3+1],velfloatbuff[nn*3+2]);
                            if (hdf_header_info[i].mass[k]==0)Pbuf[ibufindex].SetMass(massfloatbuff[nn]);
                            else Pbuf[ibufindex].SetMass(hdf_header_info[i].mass[k]);
                        }
                        else {
                            Pbuf[ibufindex].SetVelocity(veldoublebuff[nn*3],veldoublebuff[nn*3+1],veldoublebuff[nn*3+2]);
                            if (hdf_header_info[i].mass[k]==0)Pbuf[ibufindex].SetMass(massdoublebuff[nn]);
                            else Pbuf[ibufindex].SetMass(hdf_header_info[i].mass[k]);
                        }
                        if (iint) Pbuf[ibufindex].SetPID(intbuff[nn]);
                        else Pbuf[ibufindex].SetPID(longbuff[nn]);
                        Pbuf[ibufindex].SetID(nn);
                        if (k==HDFGASTYPE) Pbuf[ibufindex].SetType(GASTYPE);
                        else if (k==HDFDMTYPE) Pbuf[ibufindex].SetType(DARKTYPE);
                        else if (k==HDFSTARTYPE) Pbuf[ibufindex].SetType(STARTYPE);
                        else if (k==HDFBHTYPE) Pbuf[ibufindex].SetType(BHTYPE);

#ifdef GASON
                      if (k==HDFGASTYPE) {
                        if (ifloat) Pbuf[ibufindex].SetU(ufloatbuff[nn]);
                        else Pbuf[ibufindex].SetU(udoublebuff[nn]);
#ifdef STARON
                        if (ifloat) Pbuf[ibufindex].SetSFR(SFRfloatbuff[nn]);
                        else Pbuf[ibufindex].SetSFR(SFRdoublebuff[nn]);
                        if (ifloat) Pbuf[ibufindex].SetZmet(Zfloatbuff[nn]);
                        else Pbuf[ibufindex].SetZmet(Zdoublebuff[nn]);
#endif
                    }
#endif
#ifdef STARON
                      if (k==HDFSTARTYPE) {
                        if (ifloat) Pbuf[ibufindex].SetZmet(Zfloatbuff[nn]);
                        else Pbuf[ibufindex].SetZmet(Zdoublebuff[nn]);
                        if (ifloat) {if (Tagefloatbuff[nn]<0) Pbuf[ibufindex].SetType(WINDTYPE); Pbuf[ibufindex].SetTage(Tagefloatbuff[nn]);}
                        else {if (Tagedoublebuff[nn]<0) Pbuf[ibufindex].SetType(WINDTYPE);Pbuf[ibufindex].SetTage(Tagedoublebuff[nn]);}
                      }
#endif
#ifdef EXTRAINPUTINFO
                        if (opt.iextendedoutput)
                        {
                            Pbuf[ibufindex].SetInputFileID(i);
                            Pbuf[ibufindex].SetInputIndexInFile(nn+ninputoffset);
                        }
#endif
                      Nbuf[ibuf]++;
                      MPIAddParticletoAppropriateBuffer(ibuf, ibufindex, ireadtask, BufSize, Nbuf, Pbuf, Nlocal, Part.data(), Nreadbuf, Preadbuf);
                    }
                    ninputoffset += nchunk;
                  }
                }
                if (opt.partsearchtype==PSTDARK && opt.iBaryonSearch) {
                  for (j=1;j<=nbusetypes;j++) {
                    k=usetypes[j];
                    if (hdf_header_info[i].npart[k]<chunksize)nchunk=hdf_header_info[i].npart[k];
                    else nchunk=chunksize;
                    for(n=0;n<hdf_header_info[i].npart[k];n+=nchunk)
                    {
                      if (hdf_header_info[i].npart[k]-n<chunksize&&hdf_header_info[i].npart[k]-n>0)nchunk=hdf_header_info[i].npart[k]-n;
                      //setup hyperslab so that it is loaded into the buffer
                      //load positions
                      itemp=0;
                      //set hyperslab
                      datarank=1;
                      datadim[0]=nchunk*3;
                      chunkspace=DataSpace(datarank,datadim);
                      filespacecount[0]=nchunk;filespacecount[1]=3;
                      filespaceoffset[0]=n;filespaceoffset[1]=0;
                      //set type
                      floattype=partsdatasetall[i*NHDFTYPE*NHDFDATABLOCK+k*NHDFDATABLOCK+itemp].getFloatType();
                      if (floattype.getSize()==sizeof(float)) {HDFREALTYPE=PredType::NATIVE_FLOAT;realbuff=floatbuff;ifloat_pos=1;}
                      else {HDFREALTYPE=PredType::NATIVE_DOUBLE ;realbuff=doublebuff;ifloat_pos=0;}
                      //read hyperslab into local buffer
                      partsdataspaceall[i*NHDFTYPE*NHDFDATABLOCK+k*NHDFDATABLOCK+itemp].selectHyperslab(H5S_SELECT_SET, filespacecount, filespaceoffset);
                      partsdatasetall[i*NHDFTYPE*NHDFDATABLOCK+k*NHDFDATABLOCK+itemp].read(realbuff,HDFREALTYPE,chunkspace,partsdataspaceall[i*NHDFTYPE*NHDFDATABLOCK+k*NHDFDATABLOCK+itemp]);
                      //velocities
                      itemp++;
                      datarank=1;
                      datadim[0]=nchunk*3;
                      chunkspace=DataSpace(datarank,datadim);
                      filespacecount[0]=nchunk;filespacecount[1]=3;
                      filespaceoffset[0]=n;filespaceoffset[1]=0;
                      floattype=partsdatasetall[i*NHDFTYPE*NHDFDATABLOCK+k*NHDFDATABLOCK+itemp].getFloatType();
                      if (floattype.getSize()==sizeof(float)) {HDFREALTYPE=PredType::NATIVE_FLOAT;realbuff=velfloatbuff;ifloat=1;}
                      else {HDFREALTYPE=PredType::NATIVE_DOUBLE ;realbuff=veldoublebuff;ifloat=0;}
                      partsdataspaceall[i*NHDFTYPE*NHDFDATABLOCK+k*NHDFDATABLOCK+itemp].selectHyperslab(H5S_SELECT_SET, filespacecount, filespaceoffset);
                      partsdatasetall[i*NHDFTYPE*NHDFDATABLOCK+k*NHDFDATABLOCK+itemp].read(realbuff,HDFREALTYPE,chunkspace,partsdataspaceall[i*NHDFTYPE*NHDFDATABLOCK+k*NHDFDATABLOCK+itemp]);
                      //ids
                      itemp++;
                      datarank=1;
                      datadim[0]=nchunk;
                      chunkspace=DataSpace(datarank,datadim);
                      filespacecount[0]=nchunk;filespacecount[1]=1;
                      filespaceoffset[0]=n;filespaceoffset[1]=0;
                      inttype=partsdatasetall[i*NHDFTYPE*NHDFDATABLOCK+k*NHDFDATABLOCK+itemp].getIntType();
                      if (inttype.getSize()==sizeof(int)) {HDFINTEGERTYPE=PredType::NATIVE_INT;integerbuff=intbuff;iint=1;}
                      else {HDFINTEGERTYPE=PredType::NATIVE_LONG;integerbuff=longbuff;iint=0;}
                      partsdataspaceall[i*NHDFTYPE*NHDFDATABLOCK+k*NHDFDATABLOCK+itemp].selectHyperslab(H5S_SELECT_SET, filespacecount, filespaceoffset);
                      partsdatasetall[i*NHDFTYPE*NHDFDATABLOCK+k*NHDFDATABLOCK+itemp].read(integerbuff,HDFINTEGERTYPE,chunkspace,partsdataspaceall[i*NHDFTYPE*NHDFDATABLOCK+k*NHDFDATABLOCK+itemp]);
                      //masses
                      itemp++;
                      datarank=1;
                      datadim[0]=nchunk;
                      chunkspace=DataSpace(datarank,datadim);
                      filespacecount[0]=nchunk;filespacecount[1]=1;
                      filespaceoffset[0]=n;filespaceoffset[1]=0;
                      if (hdf_header_info[i].mass[k]==0){
                        floattype=partsdatasetall[i*NHDFTYPE*NHDFDATABLOCK+k*NHDFDATABLOCK+itemp].getFloatType();
                        if (floattype.getSize()==sizeof(float)) {HDFREALTYPE=PredType::NATIVE_FLOAT;realbuff=massfloatbuff;ifloat=1;}
                        else {HDFREALTYPE=PredType::NATIVE_DOUBLE;realbuff=massdoublebuff;ifloat=0;}
                        partsdataspaceall[i*NHDFTYPE*NHDFDATABLOCK+k*NHDFDATABLOCK+itemp].selectHyperslab(H5S_SELECT_SET, filespacecount, filespaceoffset);
                        partsdatasetall[i*NHDFTYPE*NHDFDATABLOCK+k*NHDFDATABLOCK+itemp].read(realbuff,HDFREALTYPE,chunkspace,partsdataspaceall[i*NHDFTYPE*NHDFDATABLOCK+k*NHDFDATABLOCK+itemp]);
                      }
#ifdef GASON
                      //self-energy
                      itemp++;
                      datarank=1;
                      datadim[0]=nchunk;
                      chunkspace=DataSpace(datarank,datadim);
                      filespacecount[0]=nchunk;filespacecount[1]=1;
                      filespaceoffset[0]=n;filespaceoffset[1]=0;
                      if (k==HDFGASTYPE) {
                        floattype=partsdatasetall[i*NHDFTYPE*NHDFDATABLOCK+k*NHDFDATABLOCK+itemp].getFloatType();
                        if (floattype.getSize()==sizeof(float)) {HDFREALTYPE=PredType::NATIVE_FLOAT;realbuff=ufloatbuff;ifloat=1;}
                        else {HDFREALTYPE=PredType::NATIVE_DOUBLE;realbuff=udoublebuff;ifloat=0;}
                        partsdataspaceall[i*NHDFTYPE*NHDFDATABLOCK+k*NHDFDATABLOCK+itemp].selectHyperslab(H5S_SELECT_SET, filespacecount, filespaceoffset);
                        partsdatasetall[i*NHDFTYPE*NHDFDATABLOCK+k*NHDFDATABLOCK+itemp].read(realbuff,HDFREALTYPE,chunkspace,partsdataspaceall[i*NHDFTYPE*NHDFDATABLOCK+k*NHDFDATABLOCK+itemp]);
                      }
#ifdef STARON
                      //star formation rate
                      itemp++;
                      datarank=1;
                      datadim[0]=nchunk;
                      chunkspace=DataSpace(datarank,datadim);
                      filespacecount[0]=nchunk;filespacecount[1]=1;
                      filespaceoffset[0]=n;filespaceoffset[1]=0;
                      if (k==HDFGASTYPE) {
                        floattype=partsdatasetall[i*NHDFTYPE*NHDFDATABLOCK+k*NHDFDATABLOCK+itemp].getFloatType();
                        if (floattype.getSize()==sizeof(float)) {HDFREALTYPE=PredType::NATIVE_FLOAT;realbuff=SFRfloatbuff;ifloat=1;}
                        else {HDFREALTYPE=PredType::NATIVE_DOUBLE;realbuff=SFRdoublebuff;ifloat=0;}
                        partsdataspaceall[i*NHDFTYPE*NHDFDATABLOCK+k*NHDFDATABLOCK+itemp].selectHyperslab(H5S_SELECT_SET, filespacecount, filespaceoffset);
                        partsdatasetall[i*NHDFTYPE*NHDFDATABLOCK+k*NHDFDATABLOCK+itemp].read(realbuff,HDFREALTYPE,chunkspace,partsdataspaceall[i*NHDFTYPE*NHDFDATABLOCK+k*NHDFDATABLOCK+itemp]);
                      }

                      //metallicity
                      itemp++;
                      datarank=1;
                      datadim[0]=nchunk;
                      chunkspace=DataSpace(datarank,datadim);
                      filespacecount[0]=nchunk;filespacecount[1]=1;
                      filespaceoffset[0]=n;filespaceoffset[1]=0;
                      if (k==HDFGASTYPE || k==HDFSTARTYPE) {
                        floattype=partsdatasetall[i*NHDFTYPE*NHDFDATABLOCK+k*NHDFDATABLOCK+itemp].getFloatType();
                        if (floattype.getSize()==sizeof(float)) {HDFREALTYPE=PredType::NATIVE_FLOAT;realbuff=Zfloatbuff;ifloat=1;}
                        else {HDFREALTYPE=PredType::NATIVE_DOUBLE;realbuff=Zdoublebuff;ifloat=0;}
                        partsdataspaceall[i*NHDFTYPE*NHDFDATABLOCK+k*NHDFDATABLOCK+itemp].selectHyperslab(H5S_SELECT_SET, filespacecount, filespaceoffset);
                        partsdatasetall[i*NHDFTYPE*NHDFDATABLOCK+k*NHDFDATABLOCK+itemp].read(realbuff,HDFREALTYPE,chunkspace,partsdataspaceall[i*NHDFTYPE*NHDFDATABLOCK+k*NHDFDATABLOCK+itemp]);
                      }

                      //stellar age
                      itemp++;
                      datarank=1;
                      datadim[0]=nchunk;
                      chunkspace=DataSpace(datarank,datadim);
                      filespacecount[0]=nchunk;filespacecount[1]=1;
                      filespaceoffset[0]=n;filespaceoffset[1]=0;
                      if (k==HDFSTARTYPE) {
                        floattype=partsdatasetall[i*NHDFTYPE*NHDFDATABLOCK+k*NHDFDATABLOCK+itemp].getFloatType();
                        if (floattype.getSize()==sizeof(float)) {HDFREALTYPE=PredType::NATIVE_FLOAT;realbuff=Tagefloatbuff;ifloat=1;}
                        else {HDFREALTYPE=PredType::NATIVE_DOUBLE;realbuff=Tagedoublebuff;ifloat=0;}
                        partsdataspaceall[i*NHDFTYPE*NHDFDATABLOCK+k*NHDFDATABLOCK+itemp].selectHyperslab(H5S_SELECT_SET, filespacecount, filespaceoffset);
                        partsdatasetall[i*NHDFTYPE*NHDFDATABLOCK+k*NHDFDATABLOCK+itemp].read(realbuff,HDFREALTYPE,chunkspace,partsdataspaceall[i*NHDFTYPE*NHDFDATABLOCK+k*NHDFDATABLOCK+itemp]);
                      }
#endif
#endif
                      for (int nn=0;nn<nchunk;nn++) {
                        if (ifloat_pos) ibuf=MPIGetParticlesProcessor(floatbuff[nn*3],floatbuff[nn*3+1],floatbuff[nn*3+2]);
                        else ibuf=MPIGetParticlesProcessor(doublebuff[nn*3],doublebuff[nn*3+1],doublebuff[nn*3+2]);
                        ibufindex=ibuf*BufSize+Nbuf[ibuf];
                        //reset hydro quantities of buffer
#ifdef GASON
                        Pbuf[ibufindex].SetU(0);
#ifdef STARON
                        Pbuf[ibufindex].SetSFR(0);
                        Pbuf[ibufindex].SetZmet(0);
#endif
#endif
#ifdef STARON
                        Pbuf[ibufindex].SetZmet(0);
                        Pbuf[ibufindex].SetTage(0);
#endif
#ifdef BHON
#endif
                        //store particle info in Ptemp;
                        if(ifloat_pos) {
                          Pbuf[ibufindex].SetPosition(floatbuff[nn*3],floatbuff[nn*3+1],floatbuff[nn*3+2]);
                        }
                        else {
                          Pbuf[ibufindex].SetPosition(doublebuff[nn*3],doublebuff[nn*3+1],doublebuff[nn*3+2]);
                        }
                        if (ifloat) {
                          Pbuf[ibufindex].SetVelocity(velfloatbuff[nn*3],velfloatbuff[nn*3+1],velfloatbuff[nn*3+2]);
                          if (hdf_header_info[i].mass[k]==0)Pbuf[ibufindex].SetMass(massfloatbuff[nn]);
                          else Pbuf[ibufindex].SetMass(hdf_header_info[i].mass[k]);
                        }
                        else {
                          Pbuf[ibufindex].SetVelocity(veldoublebuff[nn*3],veldoublebuff[nn*3+1],veldoublebuff[nn*3+2]);
                          if (hdf_header_info[i].mass[k]==0)Pbuf[ibufindex].SetMass(massdoublebuff[nn]);
                          else Pbuf[ibufindex].SetMass(hdf_header_info[i].mass[k]);
                        }
                        if (iint) Pbuf[ibufindex].SetPID(intbuff[nn]);
                        else Pbuf[ibufindex].SetPID(longbuff[nn]);
                        Pbuf[ibufindex].SetID(nn);
                        if (k==HDFGASTYPE) Pbuf[ibufindex].SetType(GASTYPE);
                        else if (k==HDFDMTYPE) Pbuf[ibufindex].SetType(DARKTYPE);
                        else if (k==HDFSTARTYPE) Pbuf[ibufindex].SetType(STARTYPE);
                        else if (k==HDFBHTYPE) Pbuf[ibufindex].SetType(BHTYPE);
#ifdef GASON
                        if (k==HDFGASTYPE) {
                          if (ifloat) Pbuf[ibufindex].SetU(ufloatbuff[nn]);
                          else Pbuf[ibufindex].SetU(udoublebuff[nn]);
#ifdef STARON
                          if (ifloat) Pbuf[ibufindex].SetSFR(SFRfloatbuff[nn]);
                          else Pbuf[ibufindex].SetSFR(SFRdoublebuff[nn]);
                          if (ifloat) Pbuf[ibufindex].SetZmet(Zfloatbuff[nn]);
                          else Pbuf[ibufindex].SetZmet(Zdoublebuff[nn]);
#endif
                        }
#endif
#ifdef STARON
                        if (k==HDFSTARTYPE) {
                          if (ifloat) Pbuf[ibufindex].SetZmet(Zfloatbuff[nn]);
                          else Pbuf[ibufindex].SetZmet(Zdoublebuff[nn]);
                          if (ifloat) {if (Tagefloatbuff[nn]<0) Pbuf[ibufindex].SetType(WINDTYPE); Pbuf[ibufindex].SetTage(Tagefloatbuff[nn]);}
                          else {if (Tagedoublebuff[nn]<0) Pbuf[ibufindex].SetType(WINDTYPE);Pbuf[ibufindex].SetTage(Tagedoublebuff[nn]);}
                        }
#endif
#ifdef EXTRAINPUTINFO
                        if (opt.iextendedoutput)
                        {
                            Pbuf[ibufindex].SetInputFileID(i);
                            Pbuf[ibufindex].SetInputIndexInFile(nn+ninputoffset);
                        }
#endif
                        Nbuf[ibuf]++;
                        MPIAddParticletoAppropriateBuffer(ibuf, ibufindex, ireadtask, BufSize, Nbuf, Pbuf, Nlocalbaryon[0], Pbaryons, Nreadbuf, Preadbuf);
                      }
                      ninputoffset+=nchunk;
                    }//end of chunk
                  }//end of party type
                }//end of baryon if
            }//end of try block
            catch(GroupIException error)
            {
                HDF5PrintError(error);
        		cerr<<"Error in group might suggest config file has the incorrect HDF naming convention. ";
        		cerr<<"Check HDF_name_convetion or add new naming convention updating hdfitems.h in the source code. "<<endl;
        		Fhdf[i].close();
        		MPI_Abort(MPI_COMM_WORLD,8);
        	}
            // catch failure caused by the H5File operations
            catch( FileIException error )
            {
                HDF5PrintError(error);
        		cerr<<"Error reading file. Exiting "<<endl;
        		Fhdf[i].close();
        		MPI_Abort(MPI_COMM_WORLD,8);
            }
            // catch failure caused by the DataSet operations
            catch( DataSetIException error )
            {
                HDF5PrintError(error);
        		cerr<<"Error in data set might suggest config file has the incorrect HDF naming convention. ";
        		cerr<<"Check HDF_name_convetion or update hdfio.cxx in the source code to read correct format"<<endl;
        		Fhdf[i].close();
        		MPI_Abort(MPI_COMM_WORLD,8);
            }
            // catch failure caused by the DataSpace operations
            catch( DataSpaceIException error )
            {
                HDF5PrintError(error);
        		cerr<<"Error in data space might suggest config file has the incorrect HDF naming convention. ";
        		cerr<<"Check HDF_name_convetion or update hdfio.cxx in the source code to read correct format"<<endl;
        		Fhdf[i].close();
        		MPI_Abort(MPI_COMM_WORLD,8);
            }
            // catch failure caused by the DataSpace operations
            catch( DataTypeIException error )
            {
                HDF5PrintError(error);
        		cerr<<"Error in data type might suggest need to update hdfio.cxx in the source code to read correct format"<<endl;
        		Fhdf[i].close();
        		MPI_Abort(MPI_COMM_WORLD,8);
            }
            Fhdf[i].close();
            //send info between read threads
            if (opt.nsnapread>1&&inreadsend<totreadsend){
                MPI_Allgather(Nreadbuf, opt.nsnapread, MPI_Int_t, mpi_nsend_readthread, opt.nsnapread, MPI_Int_t, mpi_comm_read);
                MPISendParticlesBetweenReadThreads(opt, Preadbuf, Part.data(), ireadtask, readtaskID, Pbaryons, mpi_comm_read, mpi_nsend_readthread, mpi_nsend_readthread_baryon);
                inreadsend++;
                for(ibuf = 0; ibuf < opt.nsnapread; ibuf++) Nreadbuf[ibuf]=0;
            }
        }//end of file if read
        //once finished reading the file if there are any particles left in the buffer broadcast them
        for(ibuf = 0; ibuf < NProcs; ibuf++) if (ireadtask[ibuf]<0)
        {
            MPI_Ssend(&Nbuf[ibuf],1,MPI_Int_t, ibuf, ibuf+NProcs, MPI_COMM_WORLD);
            if (Nbuf[ibuf]>0) {
                MPI_Ssend(&Pbuf[ibuf*BufSize], sizeof(Particle)*Nbuf[ibuf], MPI_BYTE, ibuf, ibuf, MPI_COMM_WORLD);
                Nbuf[ibuf]=0;
                //last broadcast with Nbuf[ibuf]=0 so that receiver knows no more particles are to be broadcast
                MPI_Ssend(&Nbuf[ibuf],1,MPI_Int_t,ibuf,ibuf+NProcs,MPI_COMM_WORLD);
            }
        }
        //do final send between read threads
        if (opt.nsnapread>1){
            MPI_Allgather(Nreadbuf, opt.nsnapread, MPI_Int_t, mpi_nsend_readthread, opt.nsnapread, MPI_Int_t, mpi_comm_read);
            MPISendParticlesBetweenReadThreads(opt, Preadbuf, Part.data(), ireadtask, readtaskID, Pbaryons, mpi_comm_read, mpi_nsend_readthread, mpi_nsend_readthread_baryon);
            inreadsend++;
            for(ibuf = 0; ibuf < opt.nsnapread; ibuf++) Nreadbuf[ibuf]=0;
        }
    }
    //if not reading information than waiting to receive information
    else {
        MPIReceiveParticlesFromReadThreads(opt,Pbuf,Part.data(),readtaskID, irecv, mpi_irecvflag, Nlocalthreadbuf, mpi_request,Pbaryons);
    }
#endif


#ifdef USEMPI
    if (ThisTask==0) {
#endif
      ///if gas found and Omega_b not set correctly (ie: ==0), assumes that
      ///lowest mass gas particle found corresponds to Omega_b
      ///Note that if there is mass evolution this WILL NOT WORK!
      if (opt.Omega_b==0 && MP_B==MAXVALUE){
        opt.Omega_b=MP_B/(MP_DM+MP_B)*opt.Omega_m;
        opt.Omega_cdm=opt.Omega_m-opt.Omega_b;
      }

      // SWIFT snapshots already include the 1/h factor factor,
      // so there is no need to include it.
      if(opt.ihdfnameconvention == HDFSWIFTEAGLENAMES) {
        //adjust period
        if (opt.comove) opt.p*=opt.L;
        else opt.p*=opt.L*opt.a;
      }
      else {
        //adjust period
        if (opt.comove) opt.p*=opt.L/opt.h;
        else opt.p*=opt.L/opt.h*opt.a;
      }

#ifdef USEMPI
    }
#endif
#ifdef USEMPI
    MPI_Barrier(MPI_COMM_WORLD);
    //update cosmological data and boundary in code units
    MPI_Bcast(&(opt.p),sizeof(opt.p),MPI_BYTE,0,MPI_COMM_WORLD);
    MPI_Bcast(&(opt.a),sizeof(opt.a),MPI_BYTE,0,MPI_COMM_WORLD);
    MPI_Bcast(&(opt.Omega_cdm),sizeof(opt.Omega_cdm),MPI_BYTE,0,MPI_COMM_WORLD);
    MPI_Bcast(&(opt.Omega_b),sizeof(opt.Omega_b),MPI_BYTE,0,MPI_COMM_WORLD);
    MPI_Bcast(&(opt.Omega_m),sizeof(opt.Omega_m),MPI_BYTE,0,MPI_COMM_WORLD);
    MPI_Bcast(&(opt.Omega_Lambda),sizeof(opt.Omega_Lambda),MPI_BYTE,0,MPI_COMM_WORLD);
    MPI_Bcast(&(opt.h),sizeof(opt.h),MPI_BYTE,0,MPI_COMM_WORLD);
    MPI_Bcast(&(opt.rhocrit),sizeof(opt.rhocrit),MPI_BYTE,0,MPI_COMM_WORLD);
    MPI_Bcast(&(opt.rhobg),sizeof(opt.rhobg),MPI_BYTE,0,MPI_COMM_WORLD);
    MPI_Bcast(&(opt.virlevel),sizeof(opt.virlevel),MPI_BYTE,0,MPI_COMM_WORLD);
    MPI_Bcast(&(opt.virBN98),sizeof(opt.virBN98),MPI_BYTE,0,MPI_COMM_WORLD);
#ifdef NOMASS
    MPI_Bcast(&(opt.MassValue),sizeof(opt.MassValue),MPI_BYTE,0,MPI_COMM_WORLD);
#endif
    MPI_Bcast(&(Ntotal),sizeof(Ntotal),MPI_BYTE,0,MPI_COMM_WORLD);

    MPI_Bcast(&(lscale),sizeof(lscale),MPI_BYTE,0,MPI_COMM_WORLD);
    MPI_Bcast(&(lvscale),sizeof(lvscale),MPI_BYTE,0,MPI_COMM_WORLD);
    MPI_Bcast(&(mscale),sizeof(mscale),MPI_BYTE,0,MPI_COMM_WORLD);
    MPI_Bcast(&(Hubbleflow),sizeof(Hubbleflow),MPI_BYTE,0,MPI_COMM_WORLD);
#endif
    ///If compiled with HIGHRES, the code assumes that the gadget data is a multi-resolution simulation
    ///with the lowest mass dark matter particle corresponding to the highest resolution and
    ///thus the physical linking length is assumed to be in fraction of interparticle spacing
    ///and is adjusted to a physical distance. Note that if high res and neff is not passed
    ///code assumes that lowest mass gas particle can be used to determine Omega_b and thus can be used
    ///to calculate the mean interparticle spacing.
    ///but one can also pass opt.Neff to adjust what the code thinks is the average inter particle spacing
#ifdef HIGHRES
    if (opt.Neff==-1) {
      //Once smallest mass particle is found (which should correspond to highest resolution area,
      if (opt.Omega_b==0) MP_B=0;
      LN=pow(((MP_DM+MP_B)*opt.M/opt.h)/(opt.Omega_m*3.0*opt.H*opt.h*opt.H*opt.h/(8.0*M_PI*opt.G)),1./3.)*opt.a;
    }
    else {
      LN=opt.p/(Double_t)opt.Neff;
    }
#endif
#ifdef USEMPI
    MPI_Bcast(&LN, 1, MPI_Real_t, 0, MPI_COMM_WORLD);
#endif
    ///if not an individual halo and cosmological and store scale of the highest resolution interparticle spacing to scale the physical FOF linking length
    //if (opt.iSingleHalo==0 && opt.icosmologicalin==1)
    // Set linking length when using a SWIFT snapshot
    if (opt.iSingleHalo==0)
    {
      opt.ellxscale=LN;
      opt.uinfo.eps*=LN;
    }
    //a bit of clean up
#ifdef USEMPI
    MPI_Comm_free(&mpi_comm_read);
    if (opt.iBaryonSearch) delete[] mpi_nsend_baryon;
    if (opt.nsnapread>1) {
      delete[] mpi_nsend_readthread;
      if (opt.iBaryonSearch) delete[] mpi_nsend_readthread_baryon;
      if (ireadtask[ThisTask]>=0) delete[] Preadbuf;
    }
    delete[] Nbuf;
    if (ireadtask[ThisTask]>=0) {
      delete[] Nreadbuf;
      delete[] Pbuf;
      delete[] ireadfile;
    }
    delete[] ireadtask;
    delete[] readtaskID;
#endif

#ifdef USEMPI

    double vscale = 0.0;

    // SWIFT snapshot velocities already contain the sqrt(a) factor,
    // so there is no need to include it.
    if(opt.ihdfnameconvention == HDFSWIFTEAGLENAMES) vscale = opt.V;
    else vscale = opt.V*sqrt(opt.a);

    //finally adjust to appropriate units
    for (i=0;i<Nlocal;i++)
    {
      Part[i].SetMass(Part[i].GetMass()*mscale);
      for (int j=0;j<3;j++) Part[i].SetVelocity(j,Part[i].GetVelocity(j)*vscale+Hubbleflow*Part[i].GetPosition(j));
      for (int j=0;j<3;j++) Part[i].SetPosition(j,Part[i].GetPosition(j)*lscale);
#ifdef GASON
      if (Part[i].GetType()==GASTYPE) Part[i].SetU(Part[i].GetU()*opt.V*opt.V);
#endif
    }
    if (Pbaryons!=NULL && opt.iBaryonSearch==1) {
      for (i=0;i<Nlocalbaryon[0];i++)
      {
        Pbaryons[i].SetMass(Pbaryons[i].GetMass()*mscale);
        for (int j=0;j<3;j++) Pbaryons[i].SetVelocity(j,Pbaryons[i].GetVelocity(j)*vscale+Hubbleflow*Pbaryons[i].GetPosition(j));
        for (int j=0;j<3;j++) Pbaryons[i].SetPosition(j,Pbaryons[i].GetPosition(j)*lscale);
#ifdef GASON
        Pbaryons[i].SetU(Pbaryons[i].GetU()*opt.V*opt.V);
#endif
    }
    }
#endif

}

#endif
