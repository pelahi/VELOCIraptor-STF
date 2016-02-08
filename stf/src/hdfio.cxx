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
void ReadHDF(Options &opt, Particle *&Part, const Int_t nbodies,Particle *&Pbaryons, Int_t nbaryons) 
{
    //structure stores the names of the groups in the hdf input
    char buf[2000];
    HDF_Group_Names hdf_gnames;
    //structures store names in groups
    HDF_Header *hdf_header_info;
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
    Group *headergroup;
    Group *partsgroup;
    Attribute *headerattribs;
    DataSpace *headerdataspace;
    DataSet *partsdataset;
    DataSpace *partsdataspace;
    DataSpace chunkspace;
    //buffers to load data
    int *intbuff=new int[HDFCHUNKSIZE];
    long long *longbuff=new long long[HDFCHUNKSIZE];
    float *floatbuff=new float[HDFCHUNKSIZE*3];
    double *doublebuff=new double[HDFCHUNKSIZE*3];
    void *integerbuff,*realbuff;
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

    ///array listing number of particle types used.
    ///Since Illustris contains an unused type of particles (2) and tracer particles (3) really not useful to iterate over all particle types in loops
    int nusetypes,nbusetypes;
    int usetypes[NHDFTYPE];
    if (opt.partsearchtype==PSTALL) {nusetypes=4;usetypes[0]=0;usetypes[1]=1;usetypes[2]=4;usetypes[3]=5;}
    else if (opt.partsearchtype==PSTDARK) {nusetypes=1;usetypes[0]=1;if (opt.iBaryonSearch) {nbusetypes=3;usetypes[1]=0;usetypes[2]=4;usetypes[3]=5;}}
    else if (opt.partsearchtype==PSTGAS) {nusetypes=1;usetypes[0]=0;}
    else if (opt.partsearchtype==PSTSTAR) {nusetypes=1;usetypes[0]=4;}
    else if (opt.partsearchtype==PSTBH) {nusetypes=1;usetypes[0]=5;}

    Int_t i,j,k,n,nchunk,count,bcount,itemp,count2,bcount2;

    //store cosmology
    double z,aadjust,Hubble,Hubbleflow;

    Double_t mscale,lscale,lvscale;
    Double_t MP_DM=MAXVALUE,LN,N_DM,MP_B=0;
    int ifirstfile=0,*ireadfile,ireaderror=0;
#ifndef USEMPI
    Int_t Ntotal;
    int ThisTask=0,NProcs=1;
    ireadfile=new int[opt.num_files];
    for (i=0;i<opt.num_files;i++) ireadfile[i]=1;
#endif

    //if MPI is used, read processors (all tasks with task numbers less than the number of snapshots) opens the file and loads the data into a particle buffer
    //this particle buffer is used to broadcast data to the appropriate processor
#ifdef USEMPI
    //since positions, velocities, masses are all at different points in the file,
    //to correctly assign particle to proccessor with correct velocities and mass must have several file pointers
    MPI_Status status;
    Particle *Pbuf;
    int mpi_ireaderror;

    //for parallel io
    Int_t Nlocalbuf,ibuf=0,*Nbuf, *Nreadbuf,*nreadoffset;
    Int_t *Nlocalthreadbuf,Nlocaltotalbuf;
    int *irecv, sendTask,recvTask,irecvflag, *mpi_irecvflag;
    MPI_Request *mpi_request;
    Int_t *mpi_nsend_baryon;
    if (opt.iBaryonSearch) mpi_nsend_baryon=new Int_t[NProcs*NProcs];
    //used in mpi to load access to all the data blocks of interest
    DataSet *partsdatasetall;
    DataSpace *partsdataspaceall;
    DataSpace *chunkspaceall;

    //extra blocks to store info
    float *velfloatbuff=new float[HDFCHUNKSIZE*3];
    double *veldoublebuff=new double[HDFCHUNKSIZE*3];
    float *massfloatbuff=new float[HDFCHUNKSIZE];
    double *massdoublebuff=new double[HDFCHUNKSIZE];
    float *ufloatbuff=new float[HDFCHUNKSIZE];
    double *udoublebuff=new double[HDFCHUNKSIZE];

    Nbuf=new Int_t[NProcs];
    for (int j=0;j<NProcs;j++) Nbuf[j]=0;
    nreadoffset=new Int_t[opt.nsnapread];
    if (ThisTask==0) cout<<"There are "<<opt.nsnapread<<" threads reading "<<opt.num_files<<" files "<<endl;
    if (ThisTask<opt.nsnapread)
    {
        //to temporarily store data from gadget file
        Pbuf=new Particle[BufSize*NProcs];
        Nreadbuf=new Int_t[opt.num_files];
    	for (int j=0;j<opt.num_files;j++) Nreadbuf[j]=0;

        //to determine which files the thread should read
        ireadfile=new int[opt.num_files];
        for (i=0;i<opt.num_files;i++) ireadfile[i]=0;
        int nread=opt.num_files/opt.nsnapread;
        int niread=ThisTask*nread,nfread=(ThisTask+1)*nread;
        if (ThisTask==opt.nsnapread-1) nfread=opt.num_files;
        for (i=niread;i<nfread;i++) ireadfile[i]=1;
        ifirstfile=niread;
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

#ifndef MPIREDUCEMEM
    MPIDomainExtentHDF(opt);
    if (NProcs>1) {
    MPIDomainDecompositionHDF(opt);
    MPIInitialDomainDecomposition();
    }
    MPI_Barrier(MPI_COMM_WORLD);
#endif
    if (ThisTask<opt.nsnapread) {
#endif
    Fhdf=new H5File[opt.num_files];
    hdf_header_info=new HDF_Header[opt.num_files];
    headergroup=new Group[opt.num_files];
    headerdataspace=new DataSpace[opt.num_files];
    headerattribs=new Attribute[opt.num_files];
    partsgroup=new Group[opt.num_files*NHDFTYPE];
    partsdataset=new DataSet[opt.num_files*NHDFTYPE];
    partsdataspace=new DataSpace[opt.num_files*NHDFTYPE];
#ifdef USEMPI
    partsdatasetall=new DataSet[opt.num_files*NHDFTYPE*NHDFDATABLOCK];
    partsdataspaceall=new DataSpace[opt.num_files*NHDFTYPE*NHDFDATABLOCK];
#endif
    for(i=0; i<opt.num_files; i++) {
    if(ireadfile[i])
    {
        if(opt.num_files>1) sprintf(buf,"%s.%d.hdf5",opt.fname,i);
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
            headergroup[i]=Fhdf[i].openGroup(hdf_gnames.Header_name);

            //start reading attributes
            headerattribs[i]=headergroup[i].openAttribute(hdf_header_info[i].names[hdf_header_info[i].IBoxSize]);
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

            headerattribs[i]=headergroup[i].openAttribute(hdf_header_info[i].names[hdf_header_info[i].IMass]);
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

            headerattribs[i]=headergroup[i].openAttribute(hdf_header_info[i].names[hdf_header_info[i].INuminFile]);
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

            headerattribs[i]=headergroup[i].openAttribute(hdf_header_info[i].names[hdf_header_info[i].INumTot]);
            headerdataspace[i]=headerattribs[i].getSpace();
            inttype=headerattribs[i].getIntType();
            if (inttype.getSize()==sizeof(int)) {
                headerattribs[i].read(PredType::NATIVE_INT,&intbuff[0]);
                for (k=0;k<NHDFTYPE;k++) hdf_header_info[i].npartTotal[k]=intbuff[k];
            }
            if (inttype.getSize()==sizeof(long long)) {
                headerattribs[i].read(PredType::NATIVE_LONG,&longbuff[0]);
                for (k=0;k<NHDFTYPE;k++) hdf_header_info[i].npartTotal[k]=longbuff[k];
            }

            headerattribs[i]=headergroup[i].openAttribute(hdf_header_info[i].names[hdf_header_info[i].INumTotHW]);
            headerdataspace[i]=headerattribs[i].getSpace();
            inttype=headerattribs[i].getIntType();
            if (inttype.getSize()==sizeof(int)) {
                headerattribs[i].read(PredType::NATIVE_INT,&intbuff[0]);
                for (k=0;k<NHDFTYPE;k++) hdf_header_info[i].npartTotalHW[k]=intbuff[k];
            }
            if (inttype.getSize()==sizeof(long long)) {
                headerattribs[i].read(PredType::NATIVE_LONG,&longbuff[0]);
                for (k=0;k<NHDFTYPE;k++) hdf_header_info[i].npartTotalHW[k]=longbuff[k];
            }

            headerattribs[i]=headergroup[i].openAttribute(hdf_header_info[i].names[hdf_header_info[i].IOmega0]);
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

            headerattribs[i]=headergroup[i].openAttribute(hdf_header_info[i].names[hdf_header_info[i].IOmegaL]);
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

            headerattribs[i]=headergroup[i].openAttribute(hdf_header_info[i].names[hdf_header_info[i].IRedshift]);
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

            headerattribs[i]=headergroup[i].openAttribute(hdf_header_info[i].names[hdf_header_info[i].ITime]);
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

            headerattribs[i]=headergroup[i].openAttribute(hdf_header_info[i].names[hdf_header_info[i].IHubbleParam]);
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

            headerattribs[i]=headergroup[i].openAttribute(hdf_header_info[i].names[hdf_header_info[i].INumFiles]);
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
    }
    }
    //after info read, initialise cosmological parameters
    opt.p=hdf_header_info[ifirstfile].BoxSize;
    z=hdf_header_info[ifirstfile].redshift;
    opt.a=1./(1.+z);
    opt.Omega_m=hdf_header_info[ifirstfile].Omega0;
    opt.Omega_Lambda=hdf_header_info[ifirstfile].OmegaLambda;
    opt.h=hdf_header_info[ifirstfile].HubbleParam;
    opt.Omega_cdm=opt.Omega_m-opt.Omega_b;
    //Hubble flow
    if (opt.comove) aadjust=1.0;
    else aadjust=opt.a;
    Hubble=opt.h*opt.H*sqrt((1-opt.Omega_m-opt.Omega_Lambda)*pow(aadjust,-2.0)+opt.Omega_m*pow(aadjust,-3.0)+opt.Omega_Lambda);
    opt.rhobg=3.*Hubble*Hubble/(8.0*M_PI*opt.G)*opt.Omega_m;
    //if opt.virlevel<0, then use virial overdensity based on Bryan and Norman 1998 virialization level is given by
    if (opt.virlevel<0) 
    {
        Double_t bnx=-((1-opt.Omega_m-opt.Omega_Lambda)*pow(aadjust,-2.0)+opt.Omega_Lambda)/((1-opt.Omega_m-opt.Omega_Lambda)*pow(aadjust,-2.0)+opt.Omega_m*pow(aadjust,-3.0)+opt.Omega_Lambda);
        opt.virlevel=(18.0*M_PI*M_PI+82.0*bnx-39*bnx*bnx)/opt.Omega_m;
    }
    mscale=opt.M/opt.h;lscale=opt.L/opt.h*aadjust;lvscale=opt.L/opt.h*opt.a;
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
    cout<<"File contains "<<Ntotal<<" particles at is at time "<<opt.a<<endl;
    cout<<"Particle system contains "<<nbodies<<" particles at is at time "<<opt.a<<" in a box of size "<<opt.p<<endl;
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
#ifdef USEMPI
    MPI_Allreduce(&ireaderror, &mpi_ireaderror, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    if (mpi_ireaderror) {
        MPI_Finalize();
        exit(9);
    }
#else
    if (ireaderror) exit(9);
#endif

#ifndef USEMPI
    //start loding particle data
    for(i=0; i<opt.num_files; i++) {
    if(ireadfile[i])
    {
        cout<<ThisTask<<" is reading file "<<i<<endl;
        ///\todo should be more rigorous with try/catch stuff
        //open particle group structures 
        for (j=0;j<nusetypes;j++) {k=usetypes[j]; partsgroup[i*NHDFTYPE+k]=Fhdf[i].openGroup(hdf_gnames.part_names[k]);}
        if (opt.partsearchtype==PSTDARK && opt.iBaryonSearch) for (j=1;j<=nbusetypes;j++) {k=usetypes[j];partsgroup[i*NHDFTYPE+k]=Fhdf[i].openGroup(hdf_gnames.part_names[k]);}
        itemp=0;
        //get positions
        for (j=0;j<nusetypes;j++) {
            k=usetypes[j]; 
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

                if (ifloat) for (int nn=0;nn<nchunk;nn++) Part[count++].SetPosition(floatbuff[nn*3],floatbuff[nn*3+1],floatbuff[nn*3+2]);
                else for (int nn=0;nn<nchunk;nn++) Part[count++].SetPosition(doublebuff[nn*3],doublebuff[nn*3+1],doublebuff[nn*3+2]);
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

                if (ifloat) for (int nn=0;nn<nchunk;nn++) Pbaryons[bcount++].SetPosition(floatbuff[nn*3],floatbuff[nn*3+1],floatbuff[nn*3+2]);
                else for (int nn=0;nn<nchunk;nn++) Pbaryons[bcount++].SetPosition(doublebuff[nn*3],doublebuff[nn*3+1],doublebuff[nn*3+2]);
            }
        }
        }
        //get velocities
        itemp++;
        for (j=0;j<nusetypes;j++) {
            k=usetypes[j]; 
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

                if (ifloat) for (int nn=0;nn<nchunk;nn++) Part[count++].SetVelocity(floatbuff[nn*3],floatbuff[nn*3+1],floatbuff[nn*3+2]);
                else for (int nn=0;nn<nchunk;nn++) Part[count++].SetVelocity(doublebuff[nn*3],doublebuff[nn*3+1],doublebuff[nn*3+2]);
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

                if (ifloat) for (int nn=0;nn<nchunk;nn++) Pbaryons[bcount++].SetVelocity(floatbuff[nn*3],floatbuff[nn*3+1],floatbuff[nn*3+2]);
                else for (int nn=0;nn<nchunk;nn++) Pbaryons[bcount++].SetVelocity(doublebuff[nn*3],doublebuff[nn*3+1],doublebuff[nn*3+2]);
            }
        }
        }
        //get ids
        itemp++;
        for (j=0;j<nusetypes;j++) {
            k=usetypes[j]; 
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
            if (hdf_header_info[i].npart[k]<HDFCHUNKSIZE)nchunk=hdf_header_info[i].npart[k];
            else nchunk=HDFCHUNKSIZE;
            for(n=0;n<hdf_header_info[i].npart[k];n+=nchunk)
            {
                if (hdf_header_info[i].npart[k]-n<HDFCHUNKSIZE&&hdf_header_info[i].npart[k]-n>0)nchunk=hdf_header_info[i].npart[k]-n;
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
                    count++;
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
                    bcount++;
                }
            }
        }
        }

        //get masses, note that DM do not contain a mass field
        itemp++;
        for (j=0;j<nusetypes;j++) {
            k=usetypes[j]; 
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
            if (hdf_header_info[i].npart[k]<HDFCHUNKSIZE)nchunk=hdf_header_info[i].npart[k];
            else nchunk=HDFCHUNKSIZE;
            for(n=0;n<hdf_header_info[i].npart[k];n+=nchunk)
            {
                if (hdf_header_info[i].npart[k]-n<HDFCHUNKSIZE&&hdf_header_info[i].npart[k]-n>0)nchunk=hdf_header_info[i].npart[k]-n;
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
            if (hdf_header_info[i].npart[k]<HDFCHUNKSIZE)nchunk=hdf_header_info[i].npart[k];
            else nchunk=HDFCHUNKSIZE;
            for(n=0;n<hdf_header_info[i].npart[k];n+=nchunk)
            {
                if (hdf_header_info[i].npart[k]-n<HDFCHUNKSIZE&&hdf_header_info[i].npart[k]-n>0)nchunk=hdf_header_info[i].npart[k]-n;
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
            }
        }
        count=count2;
        bcount=bcount2;
        for (j=0;j<nusetypes;j++) {
            k=usetypes[j]; 
            if (k==HDFGASTYPE) {
            //data loaded into memory in chunks
            if (hdf_header_info[i].npart[k]<HDFCHUNKSIZE)nchunk=hdf_header_info[i].npart[k];
            else nchunk=HDFCHUNKSIZE;
            for(n=0;n<hdf_header_info[i].npart[k];n+=nchunk)
            {
                if (hdf_header_info[i].npart[k]-n<HDFCHUNKSIZE&&hdf_header_info[i].npart[k]-n>0)nchunk=hdf_header_info[i].npart[k]-n;
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
            if (hdf_header_info[i].npart[k]<HDFCHUNKSIZE)nchunk=hdf_header_info[i].npart[k];
            else nchunk=HDFCHUNKSIZE;
            for(n=0;n<hdf_header_info[i].npart[k];n+=nchunk)
            {
                if (hdf_header_info[i].npart[k]-n<HDFCHUNKSIZE&&hdf_header_info[i].npart[k]-n>0)nchunk=hdf_header_info[i].npart[k]-n;
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
#endif
        }

        count2=count;
        bcount2=bcount;
    }
    }

    //finally adjust to appropriate units
    for (i=0;i<nbodies;i++)
    {
        Part[i].SetMass(Part[i].GetMass()*mscale);
        for (int j=0;j<3;j++) Part[i].SetVelocity(j,Part[i].GetVelocity(j)*opt.V*sqrt(opt.a)+Hubbleflow*Part[i].GetPosition(j));
        for (int j=0;j<3;j++) Part[i].SetPosition(j,Part[i].GetPosition(j)*lscale);
#ifdef GASON
        if (Part[i].GetType()==GASTYPE) Part[i].SetU(Part[i].GetU()*opt.V*opt.V);
#endif
    }
    if (Pbaryons!=NULL && opt.iBaryonSearch==1) {
    for (i=0;i<nbaryons;i++)
    {
        Pbaryons[i].SetMass(Pbaryons[i].GetMass()*mscale);
        for (int j=0;j<3;j++) Pbaryons[i].SetVelocity(j,Pbaryons[i].GetVelocity(j)*opt.V*sqrt(opt.a)+Hubbleflow*Pbaryons[i].GetPosition(j));
        for (int j=0;j<3;j++) Pbaryons[i].SetPosition(j,Pbaryons[i].GetPosition(j)*lscale);
#ifdef GASON
        Pbaryons[i].SetU(Pbaryons[i].GetU()*opt.V*opt.V);
#endif
    }
    }

#else
    //for all mpi threads that are reading input data, open file load access to data structures and begin loading into either local buffer or temporary buffer to be send to 
    //non-read threads
    if (ThisTask<opt.nsnapread) {
    count2=bcount2=0;
    for(i=0; i<opt.num_files; i++) {
    if(ireadfile[i])
    {
        cout<<ThisTask<<" is reading file "<<i<<endl;
        ///\todo should be more rigorous with try/catch stuff
        //open particle group structures 
        for (j=0;j<nusetypes;j++) {k=usetypes[j]; partsgroup[i*NHDFTYPE+k]=Fhdf[i].openGroup(hdf_gnames.part_names[k]);}
        if (opt.partsearchtype==PSTDARK && opt.iBaryonSearch) for (j=1;j<=nbusetypes;j++) {k=usetypes[j];partsgroup[i*NHDFTYPE+k]=Fhdf[i].openGroup(hdf_gnames.part_names[k]);}
        //open data structures that exist for all data blocks
        for (itemp=0;itemp<NHDFDATABLOCKALL;itemp++) {
            //for everything but mass no header check needed. 
            if (itemp!=3) {
                for (j=0;j<nusetypes;j++) {
                    k=usetypes[j]; 
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
#endif
        }

        for (j=0;j<nusetypes;j++) {
            k=usetypes[j];
            //data loaded into memory in chunks
            if (hdf_header_info[i].npart[k]<HDFCHUNKSIZE)nchunk=hdf_header_info[i].npart[k];
            else nchunk=HDFCHUNKSIZE;
            for(n=0;n<hdf_header_info[i].npart[k];n+=nchunk)
            {
                if (hdf_header_info[i].npart[k]-n<HDFCHUNKSIZE&&hdf_header_info[i].npart[k]-n>0)nchunk=hdf_header_info[i].npart[k]-n;
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
                if (floattype.getSize()==sizeof(float)) {HDFREALTYPE=PredType::NATIVE_FLOAT;realbuff=floatbuff;ifloat=1;}
                else {HDFREALTYPE=PredType::NATIVE_DOUBLE ;realbuff=doublebuff;ifloat=0;}
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

                for (int nn=0;nn<nchunk;nn++) {
                    if (ifloat) ibuf=MPIGetParticlesProcessor(floatbuff[nn*3],floatbuff[nn*3+1],floatbuff[nn*3+2]);
                    else ibuf=MPIGetParticlesProcessor(doublebuff[nn*3],doublebuff[nn*3+1],doublebuff[nn*3+2]);
                    //store particle info in Ptemp;
                    if (ifloat) {
                        Pbuf[ibuf*BufSize+Nbuf[ibuf]].SetPosition(floatbuff[nn*3],floatbuff[nn*3+1],floatbuff[nn*3+2]);
                        Pbuf[ibuf*BufSize+Nbuf[ibuf]].SetVelocity(velfloatbuff[nn*3],velfloatbuff[nn*3+1],velfloatbuff[nn*3+2]);
                        if (hdf_header_info[i].mass[k]==0)Pbuf[ibuf*BufSize+Nbuf[ibuf]].SetMass(massfloatbuff[nn]);
                        else Pbuf[ibuf*BufSize+Nbuf[ibuf]].SetMass(hdf_header_info[i].mass[k]);
                    }
                    else {
                        Pbuf[ibuf*BufSize+Nbuf[ibuf]].SetPosition(doublebuff[nn*3],doublebuff[nn*3+1],doublebuff[nn*3+2]);
                        Pbuf[ibuf*BufSize+Nbuf[ibuf]].SetVelocity(veldoublebuff[nn*3],veldoublebuff[nn*3+1],veldoublebuff[nn*3+2]);
                        if (hdf_header_info[i].mass[k]==0)Pbuf[ibuf*BufSize+Nbuf[ibuf]].SetMass(massdoublebuff[nn]);
                        else Pbuf[ibuf*BufSize+Nbuf[ibuf]].SetMass(hdf_header_info[i].mass[k]);
                    }
                    if (iint) Pbuf[ibuf*BufSize+Nbuf[ibuf]].SetPID(intbuff[nn]);
                    else Pbuf[ibuf*BufSize+Nbuf[ibuf]].SetPID(longbuff[nn]);
                    Pbuf[ibuf*BufSize+Nbuf[ibuf]].SetID(nn);
                    if (k==HDFGASTYPE) Pbuf[ibuf*BufSize+Nbuf[ibuf]].SetType(GASTYPE);
                    else if (k==HDFDMTYPE) Pbuf[ibuf*BufSize+Nbuf[ibuf]].SetType(DARKTYPE);
                    else if (k==HDFSTARTYPE) Pbuf[ibuf*BufSize+Nbuf[ibuf]].SetType(STARTYPE);
                    else if (k==HDFBHTYPE) Pbuf[ibuf*BufSize+Nbuf[ibuf]].SetType(BHTYPE);
#ifdef GASON
                    if (k==HDFGASTYPE) {
                        if (ifloat) Pbuf[ibuf*BufSize+Nbuf[ibuf]].SetU(ufloatbuff[nn]); 
                        else Pbuf[ibuf*BufSize+Nbuf[ibuf]].SetU(udoublebuff[nn]);
#ifdef STARON
#endif

                    }
#endif
#ifdef STARON
                    if (k==HDFSTARTYPE) {
                    }
#endif
                    if(ibuf<opt.nsnapread&&ibuf!=ThisTask) Nreadbuf[ibuf]++;
                    Nbuf[ibuf]++;
                    if (ibuf==ThisTask) {
                        Nbuf[ibuf]--;
                        Part[Nlocal++]=Pbuf[ibuf*BufSize+Nbuf[ibuf]];
                    }
                    else {
                        //before a simple send was done because only Task zero was reading the data
                        //but now if ibuf<opt.nsnapread, care must be taken.
                        //blocking sends that are matched by non-blocking receives
                        if(Nbuf[ibuf]==BufSize&&ibuf>=opt.nsnapread) {
                            MPI_Send(&Nbuf[ibuf], 1, MPI_Int_t, ibuf, ibuf+NProcs, MPI_COMM_WORLD);
                            MPI_Send(&Pbuf[ibuf*BufSize],sizeof(Particle)*Nbuf[ibuf],MPI_BYTE,ibuf,ibuf,MPI_COMM_WORLD);
                            Nbuf[ibuf]=0;
                        }
                        else if (Nbuf[ibuf]==BufSize&&ibuf<opt.nsnapread) {
                            Nbuf[ibuf]=0;
                        }
                    }
                }
            }
        }
        if (opt.partsearchtype==PSTDARK && opt.iBaryonSearch) {
        for (j=1;j<=nbusetypes;j++) {
            k=usetypes[j];
            if (hdf_header_info[i].npart[k]<HDFCHUNKSIZE)nchunk=hdf_header_info[i].npart[k];
            else nchunk=HDFCHUNKSIZE;
            for(n=0;n<hdf_header_info[i].npart[k];n+=nchunk)
            {
                if (hdf_header_info[i].npart[k]-n<HDFCHUNKSIZE&&hdf_header_info[i].npart[k]-n>0)nchunk=hdf_header_info[i].npart[k]-n;
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
                if (floattype.getSize()==sizeof(float)) {HDFREALTYPE=PredType::NATIVE_FLOAT;realbuff=floatbuff;ifloat=1;}
                else {HDFREALTYPE=PredType::NATIVE_DOUBLE ;realbuff=doublebuff;ifloat=0;}
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

                for (int nn=0;nn<nchunk;nn++) {
                    if (ifloat) ibuf=MPIGetParticlesProcessor(floatbuff[nn*3],floatbuff[nn*3+1],floatbuff[nn*3+2]);
                    else ibuf=MPIGetParticlesProcessor(doublebuff[nn*3],doublebuff[nn*3+1],doublebuff[nn*3+2]);
                    //store particle info in Ptemp;
                    if (ifloat) {
                        Pbuf[ibuf*BufSize+Nbuf[ibuf]].SetPosition(floatbuff[nn*3],floatbuff[nn*3+1],floatbuff[nn*3+2]);
                        Pbuf[ibuf*BufSize+Nbuf[ibuf]].SetVelocity(velfloatbuff[nn*3],velfloatbuff[nn*3+1],velfloatbuff[nn*3+2]);
                        if (hdf_header_info[i].mass[k]==0)Pbuf[ibuf*BufSize+Nbuf[ibuf]].SetMass(massfloatbuff[nn]);
                        else Pbuf[ibuf*BufSize+Nbuf[ibuf]].SetMass(hdf_header_info[i].mass[k]);
                    }
                    else {
                        Pbuf[ibuf*BufSize+Nbuf[ibuf]].SetPosition(doublebuff[nn*3],doublebuff[nn*3+1],doublebuff[nn*3+2]);
                        Pbuf[ibuf*BufSize+Nbuf[ibuf]].SetVelocity(veldoublebuff[nn*3],veldoublebuff[nn*3+1],veldoublebuff[nn*3+2]);
                        if (hdf_header_info[i].mass[k]==0)Pbuf[ibuf*BufSize+Nbuf[ibuf]].SetMass(massdoublebuff[nn]);
                        else Pbuf[ibuf*BufSize+Nbuf[ibuf]].SetMass(hdf_header_info[i].mass[k]);
                    }
                    if (iint) Pbuf[ibuf*BufSize+Nbuf[ibuf]].SetPID(intbuff[nn]);
                    else Pbuf[ibuf*BufSize+Nbuf[ibuf]].SetPID(longbuff[nn]);
                    Pbuf[ibuf*BufSize+Nbuf[ibuf]].SetID(nn);
                    if (k==HDFGASTYPE) Pbuf[ibuf*BufSize+Nbuf[ibuf]].SetType(GASTYPE);
                    else if (k==HDFDMTYPE) Pbuf[ibuf*BufSize+Nbuf[ibuf]].SetType(DARKTYPE);
                    else if (k==HDFSTARTYPE) Pbuf[ibuf*BufSize+Nbuf[ibuf]].SetType(STARTYPE);
                    else if (k==HDFBHTYPE) Pbuf[ibuf*BufSize+Nbuf[ibuf]].SetType(BHTYPE);
#ifdef GASON
                    if (k==HDFGASTYPE) {
                        if (ifloat) Pbuf[ibuf*BufSize+Nbuf[ibuf]].SetU(ufloatbuff[nn]); 
                        else Pbuf[ibuf*BufSize+Nbuf[ibuf]].SetU(udoublebuff[nn]);
#ifdef STARON
#endif

                    }
#endif
#ifdef STARON
                    if (k==HDFSTARTYPE) {
                    }
#endif
                    if(ibuf<opt.nsnapread&&ibuf!=ThisTask) Nreadbuf[ibuf]++;
                    Nbuf[ibuf]++;
                    if (ibuf==ThisTask) {
                        Nbuf[ibuf]--;
                        //Pbaryons[Nlocalbaryon[0]++]=Pbuf[ibuf*BufSize+Nbuf[ibuf]];
                        Part[Nlocal++]=Pbuf[ibuf*BufSize+Nbuf[ibuf]];
                        //if (k==HDFGASTYPE) Nlocalbaryon[1]++;
                        //else if (k==HDFSTARTYPE) Nlocalbaryon[2]++;
                        //else if (k==HDFBHTYPE) Nlocalbaryon[3]++;
                    }
                    else {
                        //before a simple send was done because only Task zero was reading the data
                        //but now if ibuf<opt.nsnapread, care must be taken.
                        //blocking sends that are matched by non-blocking receives
                        if(Nbuf[ibuf]==BufSize&&ibuf>=opt.nsnapread) {
                            MPI_Send(&Nbuf[ibuf], 1, MPI_Int_t, ibuf, ibuf+NProcs, MPI_COMM_WORLD);
                            MPI_Send(&Pbuf[ibuf*BufSize],sizeof(Particle)*Nbuf[ibuf],MPI_BYTE,ibuf,ibuf,MPI_COMM_WORLD);
                            Nbuf[ibuf]=0;
                        }
                        else if (Nbuf[ibuf]==BufSize&&ibuf<opt.nsnapread) {
                            Nbuf[ibuf]=0;
                        }
                    }
                }
            }
        }
        }

    }
    }
    //once finished reading the file if there are any particles left in the buffer broadcast them
    for(ibuf = opt.nsnapread; ibuf < NProcs; ibuf++)
    {
        MPI_Ssend(&Nbuf[ibuf],1,MPI_Int_t, ibuf, ibuf+NProcs, MPI_COMM_WORLD);
        if (Nbuf[ibuf]>0) {
            MPI_Ssend(&Pbuf[ibuf*BufSize], sizeof(Particle)*Nbuf[ibuf], MPI_BYTE, ibuf, ibuf, MPI_COMM_WORLD);
            Nbuf[ibuf]=0;
            //last broadcast with Nbuf[ibuf]=0 so that receiver knows no more particles are to be broadcast
            MPI_Ssend(&Nbuf[ibuf],1,MPI_Int_t,ibuf,ibuf+NProcs,MPI_COMM_WORLD);
        }
    }
    }
    //if not reading information than waiting to receive information
    else {
        //for all threads not reading snapshots, simply receive particles as necessary from all threads involved with reading the data
        //first determine which threads are going to send information to this thread.
        for (i=0;i<opt.nsnapread;i++) if (irecv[i]) {
            mpi_irecvflag[i]=0;
            MPI_Irecv(&Nlocalthreadbuf[i], 1, MPI_Int_t, i, ThisTask+NProcs, MPI_COMM_WORLD, &mpi_request[i]);
        }
        Nlocaltotalbuf=0;
        //non-blocking receives for the number of particles one expects to receive
        do {
            irecvflag=0;
            for (i=0;i<opt.nsnapread;i++) if (irecv[i]) {
                if (mpi_irecvflag[i]==0) {
                    //test if a request has been sent for a Recv call by one of the read threads
                    MPI_Test(&mpi_request[i], &mpi_irecvflag[i], &status);
                    if (mpi_irecvflag[i]) {
                        if (Nlocalthreadbuf[i]>0) {
                            MPI_Recv(&Part[Nlocal],sizeof(Particle)*Nlocalthreadbuf[i],MPI_BYTE,i,ThisTask, MPI_COMM_WORLD,&status);
                            Nlocal+=Nlocalthreadbuf[i];
                            Nlocaltotalbuf+=Nlocalthreadbuf[i];
                            mpi_irecvflag[i]=0;
                            MPI_Irecv(&Nlocalthreadbuf[i], 1, MPI_Int_t, i, ThisTask+NProcs, MPI_COMM_WORLD, &mpi_request[i]);
                        }
                        else {
                            irecv[i]=0;
                        }
                    }
                }
            }
            for (i=0;i<opt.nsnapread;i++) irecvflag+=irecv[i];
        } while(irecvflag>0);
        //now that data is local, must adjust data iff a separate baryon search is required. 
        if (opt.partsearchtype==PSTDARK && opt.iBaryonSearch) {
            for (i=0;i<Nlocal;i++) {
                k=Part[i].GetType();
                if (!(k==GASTYPE||k==STARTYPE||k==BHTYPE)) Part[i].SetID(0);
                else {
                    Nlocalbaryon[0]++;
                    if  (k==GASTYPE) {Part[i].SetID(1);Nlocalbaryon[1]++;}
                    else if  (k==STARTYPE) {Part[i].SetID(2);Nlocalbaryon[2]++;}
                    else if  (k==BHTYPE) {Part[i].SetID(3);Nlocalbaryon[3]++;}
                }
            }
            //sorted so that dark matter particles first, baryons after
            qsort(Part,Nlocal, sizeof(Particle), IDCompare);
            Nlocal-=Nlocalbaryon[0];
            //index type separated
            for (i=0;i<Nlocal;i++) Part[i].SetID(i);
            for (i=0;i<Nlocalbaryon[0];i++) Part[i+Nlocal].SetID(i+Nlocal);
            //finally, need to move baryons forward by the Export Factor * Nlocal as need that extra buffer to copy data two and from mpi threads
//#ifndef MPIREDUCE
//            for (i=Nlocalbaryon[0]-1;i>=0;i--) Part[i+(Int_t)(Nlocal*MPIExportFac)]=Part[i+Nlocal];
//#endif
        }
    }

    //finally need to send info between read threads once all threads reading data have broadcasted the data appropriately to all other threads
    //must deallocate Pbuf, reallocate it for the local Nreadbuf amount and read the files again. This ensures little memory overhead and files are only read twice

    //since Nbuf is used to determine what is going to be sent between threads in point-to-point communication
    //via an allgather, reset Nbuf
    for (i=0;i<NProcs;i++) Nbuf[i]=0;
    if (ThisTask<opt.nsnapread && opt.nsnapread>1) {
    delete[] Pbuf;
    Nlocalbuf=0;
    for (i=0;i<opt.nsnapread;i++) Nlocalbuf+=Nreadbuf[i];
    if (Nlocalbuf>0)
    {
    Pbuf=new Particle[Nlocalbuf];
    //determine offsets
    nreadoffset[0]=0;for (i=1;i<opt.nsnapread;i++)nreadoffset[i]=nreadoffset[i-1]+Nreadbuf[i-1];
    for(i=0;i<opt.num_files; i++)
    if (ireadfile[i])
    {
        for (j=0;j<nusetypes;j++) {
            k=usetypes[j];
            //data loaded into memory in chunks
            if (hdf_header_info[i].npart[k]<HDFCHUNKSIZE)nchunk=hdf_header_info[i].npart[k];
            else nchunk=HDFCHUNKSIZE;
            for(n=0;n<hdf_header_info[i].npart[k];n+=nchunk)
            {
                if (hdf_header_info[i].npart[k]-n<HDFCHUNKSIZE&&hdf_header_info[i].npart[k]-n>0)nchunk=hdf_header_info[i].npart[k]-n;
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
                if (floattype.getSize()==sizeof(float)) {HDFREALTYPE=PredType::NATIVE_FLOAT;realbuff=floatbuff;ifloat=1;}
                else {HDFREALTYPE=PredType::NATIVE_DOUBLE ;realbuff=doublebuff;ifloat=0;}
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

                for (int nn=0;nn<nchunk;nn++) {
                    if (ifloat) ibuf=MPIGetParticlesProcessor(floatbuff[nn*3],floatbuff[nn*3+1],floatbuff[nn*3+2]);
                    else ibuf=MPIGetParticlesProcessor(doublebuff[nn*3],doublebuff[nn*3+1],doublebuff[nn*3+2]);
                    if (ibuf<opt.nsnapread&&ibuf!=ThisTask) {
                    //store particle info in Ptemp;
                    if (ifloat) {
                        Pbuf[nreadoffset[ibuf]+Nbuf[ibuf]].SetPosition(floatbuff[nn*3],floatbuff[nn*3+1],floatbuff[nn*3+2]);
                        Pbuf[nreadoffset[ibuf]+Nbuf[ibuf]].SetVelocity(velfloatbuff[nn*3],velfloatbuff[nn*3+1],velfloatbuff[nn*3+2]);
                        if (hdf_header_info[i].mass[k]==0)Pbuf[nreadoffset[ibuf]+Nbuf[ibuf]].SetMass(massfloatbuff[nn]);
                        else Pbuf[nreadoffset[ibuf]+Nbuf[ibuf]].SetMass(hdf_header_info[i].mass[k]);
                    }
                    else {
                        Pbuf[nreadoffset[ibuf]+Nbuf[ibuf]].SetPosition(doublebuff[nn*3],doublebuff[nn*3+1],doublebuff[nn*3+2]);
                        Pbuf[nreadoffset[ibuf]+Nbuf[ibuf]].SetVelocity(veldoublebuff[nn*3],veldoublebuff[nn*3+1],veldoublebuff[nn*3+2]);
                        if (hdf_header_info[i].mass[k]==0)Pbuf[nreadoffset[ibuf]+Nbuf[ibuf]].SetMass(massdoublebuff[nn]);
                        else Pbuf[nreadoffset[ibuf]+Nbuf[ibuf]].SetMass(hdf_header_info[i].mass[k]);
                    }
                    if (iint) Pbuf[nreadoffset[ibuf]+Nbuf[ibuf]].SetPID(intbuff[nn]);
                    else Pbuf[nreadoffset[ibuf]+Nbuf[ibuf]].SetPID(longbuff[nn]);
                    Pbuf[nreadoffset[ibuf]+Nbuf[ibuf]].SetID(nn);
                    if (k==HDFGASTYPE) Pbuf[nreadoffset[ibuf]+Nbuf[ibuf]].SetType(GASTYPE);
                    else if (k==HDFDMTYPE) Pbuf[nreadoffset[ibuf]+Nbuf[ibuf]].SetType(DARKTYPE);
                    else if (k==HDFSTARTYPE) Pbuf[nreadoffset[ibuf]+Nbuf[ibuf]].SetType(STARTYPE);
                    else if (k==HDFBHTYPE) Pbuf[nreadoffset[ibuf]+Nbuf[ibuf]].SetType(BHTYPE);
#ifdef GASON
                    if (k==HDFGASTYPE) {
                        if (ifloat) Pbuf[nreadoffset[ibuf]+Nbuf[ibuf]].SetU(ufloatbuff[nn]); 
                        else Pbuf[nreadoffset[ibuf]+Nbuf[ibuf]].SetU(udoublebuff[nn]);
#ifdef STARON
#endif

                    }
#endif
#ifdef STARON
                    if (k==HDFSTARTYPE) {
                    }
#endif
                    Nbuf[ibuf]++;
                    }
                }
            }
        }
        if (opt.partsearchtype==PSTDARK && opt.iBaryonSearch) {
        for (j=1;j<=nbusetypes;j++) {
            k=usetypes[j];
            if (hdf_header_info[i].npart[k]<HDFCHUNKSIZE)nchunk=hdf_header_info[i].npart[k];
            else nchunk=HDFCHUNKSIZE;
            for(n=0;n<hdf_header_info[i].npart[k];n+=nchunk)
            {
                if (hdf_header_info[i].npart[k]-n<HDFCHUNKSIZE&&hdf_header_info[i].npart[k]-n>0)nchunk=hdf_header_info[i].npart[k]-n;
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
                if (floattype.getSize()==sizeof(float)) {HDFREALTYPE=PredType::NATIVE_FLOAT;realbuff=floatbuff;ifloat=1;}
                else {HDFREALTYPE=PredType::NATIVE_DOUBLE ;realbuff=doublebuff;ifloat=0;}
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

                for (int nn=0;nn<nchunk;nn++) {
                    if (ifloat) ibuf=MPIGetParticlesProcessor(floatbuff[nn*3],floatbuff[nn*3+1],floatbuff[nn*3+2]);
                    else ibuf=MPIGetParticlesProcessor(doublebuff[nn*3],doublebuff[nn*3+1],doublebuff[nn*3+2]);
                    if (ibuf<opt.nsnapread&&ibuf!=ThisTask) {
                    //store particle info in Ptemp;
                    if (ifloat) {
                        Pbuf[nreadoffset[ibuf]+Nbuf[ibuf]].SetPosition(floatbuff[nn*3],floatbuff[nn*3+1],floatbuff[nn*3+2]);
                        Pbuf[nreadoffset[ibuf]+Nbuf[ibuf]].SetVelocity(velfloatbuff[nn*3],velfloatbuff[nn*3+1],velfloatbuff[nn*3+2]);
                        if (hdf_header_info[i].mass[k]==0)Pbuf[nreadoffset[ibuf]+Nbuf[ibuf]].SetMass(massfloatbuff[nn]);
                        else Pbuf[nreadoffset[ibuf]+Nbuf[ibuf]].SetMass(hdf_header_info[i].mass[k]);
                    }
                    else {
                        Pbuf[nreadoffset[ibuf]+Nbuf[ibuf]].SetPosition(doublebuff[nn*3],doublebuff[nn*3+1],doublebuff[nn*3+2]);
                        Pbuf[nreadoffset[ibuf]+Nbuf[ibuf]].SetVelocity(veldoublebuff[nn*3],veldoublebuff[nn*3+1],veldoublebuff[nn*3+2]);
                        if (hdf_header_info[i].mass[k]==0)Pbuf[nreadoffset[ibuf]+Nbuf[ibuf]].SetMass(massdoublebuff[nn]);
                        else Pbuf[nreadoffset[ibuf]+Nbuf[ibuf]].SetMass(hdf_header_info[i].mass[k]);
                    }
                    if (iint) Pbuf[nreadoffset[ibuf]+Nbuf[ibuf]].SetPID(intbuff[nn]);
                    else Pbuf[nreadoffset[ibuf]+Nbuf[ibuf]].SetPID(longbuff[nn]);
                    Pbuf[nreadoffset[ibuf]+Nbuf[ibuf]].SetID(nn);
                    if (k==HDFGASTYPE) Pbuf[nreadoffset[ibuf]+Nbuf[ibuf]].SetType(GASTYPE);
                    else if (k==HDFDMTYPE) Pbuf[nreadoffset[ibuf]+Nbuf[ibuf]].SetType(DARKTYPE);
                    else if (k==HDFSTARTYPE) Pbuf[nreadoffset[ibuf]+Nbuf[ibuf]].SetType(STARTYPE);
                    else if (k==HDFBHTYPE) Pbuf[nreadoffset[ibuf]+Nbuf[ibuf]].SetType(BHTYPE);
#ifdef GASON
                    if (k==HDFGASTYPE) {
                        if (ifloat) Pbuf[nreadoffset[ibuf]+Nbuf[ibuf]].SetU(ufloatbuff[nn]); 
                        else Pbuf[nreadoffset[ibuf]+Nbuf[ibuf]].SetU(udoublebuff[nn]);
#ifdef STARON
#endif

                    }
#endif
#ifdef STARON
                    if (k==HDFSTARTYPE) {
                    }
#endif
                    Nbuf[ibuf]++;
                    }
                }

            }
        }
        }
        //more information contained in sph particles and if there is sf feed back but for the moment, ignore
        Fhdf[i].close();
    }
    }
    }

    //gather all the items that must be sent.
    MPI_Allgather(Nbuf, NProcs, MPI_Int_t, mpi_nsend, NProcs, MPI_Int_t, MPI_COMM_WORLD);
    //if separate baryon search then sort the Pbuf array so that it is separated by type 
    if (opt.partsearchtype==PSTDARK && opt.iBaryonSearch) {
        if (ThisTask<opt.nsnapread) {
        for(ibuf = 0; ibuf < opt.nsnapread; ibuf++) if (mpi_nsend[ThisTask * NProcs + ibuf] > 0)
        {
            Nbuf[ibuf]=0;
            for (i=0;i<mpi_nsend[ThisTask * NProcs + ibuf];i++) {
                k=Pbuf[nreadoffset[ibuf]+i].GetType();
                if (!(k==GASTYPE||k==STARTYPE||k==BHTYPE)) Pbuf[nreadoffset[ibuf]+i].SetID(0);
                else {
                    if  (k==GASTYPE) Pbuf[nreadoffset[ibuf]+i].SetID(1);
                    else if  (k==STARTYPE) Pbuf[nreadoffset[ibuf]+i].SetID(2);
                    else if  (k==BHTYPE) Pbuf[nreadoffset[ibuf]+i].SetID(3);
                    Nbuf[ibuf]++;
                }
            }
            qsort(&Pbuf[nreadoffset[ibuf]],mpi_nsend[ThisTask*NProcs+ibuf], sizeof(Particle), IDCompare);
        }
	    }
        MPI_Allgather(Nbuf, NProcs, MPI_Int_t, mpi_nsend_baryon, NProcs, MPI_Int_t, MPI_COMM_WORLD);
        for (ibuf=0;ibuf<NProcs*NProcs;ibuf++) mpi_nsend[ibuf]-=mpi_nsend_baryon[ibuf];
    }
    //and then send all the data between the read threads
    if (ThisTask<opt.nsnapread) {
    for(ibuf = 0; ibuf < opt.nsnapread; ibuf++)
    {
        if (ibuf!=ThisTask)
        {
            sendTask = ThisTask;
            recvTask = ibuf;
            if(mpi_nsend[ThisTask * NProcs + recvTask] > 0 || mpi_nsend[recvTask * NProcs + ThisTask] > 0)
            {
                //blocking point-to-point send and receive. Here must determine the appropriate offset point in the local export buffer
                //for sending data and also the local appropriate offset in the local the receive buffer for information sent from the local receiving buffer
                MPI_Sendrecv(&Pbuf[nreadoffset[recvTask]],sizeof(Particle)*mpi_nsend[ThisTask * NProcs + recvTask], MPI_BYTE, recvTask, TAG_IO_A,
                    &Part[Nlocal],sizeof(Particle)*mpi_nsend[recvTask * NProcs + ThisTask], MPI_BYTE, recvTask, TAG_IO_A, MPI_COMM_WORLD, &status);
                Nlocal+=mpi_nsend[recvTask * NProcs + ThisTask];
            }
        }
    }
    }
    if (ThisTask<opt.nsnapread) {
    if (opt.partsearchtype==PSTDARK && opt.iBaryonSearch) {
    for(ibuf = 0; ibuf < opt.nsnapread; ibuf++)
    {
        if (ibuf!=ThisTask)
        {
            sendTask = ThisTask;
            recvTask = ibuf;
            //if separate baryon search, send baryons too
                if(mpi_nsend_baryon[ThisTask * NProcs + recvTask] > 0 || mpi_nsend_baryon[recvTask * NProcs + ThisTask] > 0) 
                MPI_Sendrecv(&Pbuf[nreadoffset[recvTask]+mpi_nsend[ThisTask * NProcs + recvTask]],sizeof(Particle)*mpi_nsend_baryon[ThisTask * NProcs + recvTask], MPI_BYTE, recvTask, TAG_IO_B,
                    &Part[Nlocal],sizeof(Particle)*mpi_nsend_baryon[recvTask * NProcs + ThisTask], MPI_BYTE, recvTask, TAG_IO_B, MPI_COMM_WORLD, &status);
                Nlocal+=mpi_nsend_baryon[recvTask * NProcs + ThisTask];
                //MPI_Sendrecv(&Pbuf[nreadoffset[recvTask]+mpi_nsend[ThisTask * NProcs + recvTask]],sizeof(Particle)*mpi_nsend_baryon[ThisTask * NProcs + recvTask], MPI_BYTE, recvTask, TAG_IO_B,
                //    &Pbaryons[Nlocalbaryon[0]],sizeof(Particle)*mpi_nsend_baryon[recvTask * NProcs + ThisTask], MPI_BYTE, recvTask, TAG_IO_B, MPI_COMM_WORLD, &status);
                //Nlocalbaryon[0]+=mpi_nsend_baryon[recvTask * NProcs + ThisTask];
        }
    }

    for (i=0;i<Nlocal;i++) {
        k=Part[i].GetType();
        if (!(k==GASTYPE||k==STARTYPE||k==BHTYPE)) Part[i].SetID(0);
        else {
            Nlocalbaryon[0]++;
            if  (k==GASTYPE) {Part[i].SetID(1);Nlocalbaryon[1]++;}
            else if  (k==STARTYPE) {Part[i].SetID(2);Nlocalbaryon[2]++;}
            else if  (k==BHTYPE) {Part[i].SetID(3);Nlocalbaryon[3]++;}
        }
    }
            //sorted so that dark matter particles first, baryons after
            qsort(Part,Nlocal, sizeof(Particle), IDCompare);
            Nlocal-=Nlocalbaryon[0];
            //index type separated
            for (i=0;i<Nlocal;i++) Part[i].SetID(i);
            for (i=0;i<Nlocalbaryon[0];i++) Part[i+Nlocal].SetID(i+Nlocal);
            //finally, need to move baryons forward by the Export Factor * Nlocal as need that extra buffer to copy data two and from mpi threads
//#ifndef MPIREDUCE
            //for (i=Nlocalbaryon[0]-1;i>=0;i--) Part[i+(Int_t)(Nlocal*MPIExportFac)]=Part[i+Nlocal];
//#endif
    delete[] mpi_nsend_baryon;
    }
    //set IDS
    for (i=0;i<Nlocal;i++) Part[i].SetID(i);
    if (opt.iBaryonSearch) for (i=0;i<Nlocalbaryon[0];i++) Pbaryons[i].SetID(i+Nlocal);
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
    //adjust period
    if (opt.comove) opt.p*=opt.L/opt.h;
    else opt.p*=opt.L/opt.h*opt.a;
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
    MPI_Bcast(&(opt.rhobg),sizeof(opt.rhobg),MPI_BYTE,0,MPI_COMM_WORLD);
    MPI_Bcast(&(opt.virlevel),sizeof(opt.virlevel),MPI_BYTE,0,MPI_COMM_WORLD);
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
    ///if not an individual halo, assume cosmological and store scale of the highest resolution interparticle spacing to scale the physical FOF linking length 
    if (opt.iSingleHalo==0) 
    {
        opt.ellxscale=LN;
        opt.uinfo.eps*=LN;
    }

#ifdef USEMPI
    //finally adjust to appropriate units
    for (i=0;i<Nlocal;i++)
    {
        Part[i].SetMass(Part[i].GetMass()*mscale);
        for (int j=0;j<3;j++) Part[i].SetVelocity(j,Part[i].GetVelocity(j)*opt.V*sqrt(opt.a)+Hubbleflow*Part[i].GetPosition(j));
        for (int j=0;j<3;j++) Part[i].SetPosition(j,Part[i].GetPosition(j)*lscale);
#ifdef GASON
        if (Part[i].GetType()==GASTYPE) Part[i].SetU(Part[i].GetU()*opt.V*opt.V);
#endif
    }
    if (Pbaryons!=NULL && opt.iBaryonSearch==1) {
    for (i=0;i<Nlocalbaryon[0];i++)
    {
        Pbaryons[i].SetMass(Pbaryons[i].GetMass()*mscale);
        for (int j=0;j<3;j++) Pbaryons[i].SetVelocity(j,Pbaryons[i].GetVelocity(j)*opt.V*sqrt(opt.a)+Hubbleflow*Pbaryons[i].GetPosition(j));
        for (int j=0;j<3;j++) Pbaryons[i].SetPosition(j,Pbaryons[i].GetPosition(j)*lscale);
#ifdef GASON
        Pbaryons[i].SetU(Pbaryons[i].GetU()*opt.V*opt.V);
#endif
    }
    }
#endif

}


#endif
