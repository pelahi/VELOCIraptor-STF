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
    HDF_Part_Info hdf_extradm_info(HDFDM1TYPE,opt.ihdfnameconvention);
    HDF_Part_Info hdf_tracer_info(HDFTRACERTYPE,opt.ihdfnameconvention);
    HDF_Part_Info hdf_star_info(HDFSTARTYPE,opt.ihdfnameconvention);
    HDF_Part_Info hdf_bh_info(HDFBHTYPE,opt.ihdfnameconvention);

    HDF_Part_Info *hdf_parts[NHDFTYPE];
    hdf_parts[0]=&hdf_gas_info;
    hdf_parts[1]=&hdf_dm_info;
    #ifdef HIGHRES
    hdf_parts[2]=&hdf_extradm_info;
    hdf_parts[3]=&hdf_extradm_info;
    #else
    hdf_parts[2]=&hdf_extradm_info;
    hdf_parts[3]=&hdf_tracer_info;
    #endif
    hdf_parts[4]=&hdf_star_info;
    hdf_parts[5]=&hdf_bh_info;

    //to store the groups, data sets and their associated data spaces
    // H5File *Fhdf;
    // Group *partsgroup;
    // Attribute *headerattribs;
    // DataSpace *headerdataspace;
    // DataSet *partsdataset;
    // DataSpace *partsdataspace;
    // DataSpace chunkspace;
    vector<hid_t> Fhdf;
    vector<hid_t> partsgroup;
    vector<hid_t> headerattribs;
    vector<hid_t> headerdataspace;
    vector<hid_t> partsdataset;
    vector<hid_t> partsdataspace;
    hid_t chunkspace;
    int chunksize=opt.inputbufsize;
    //buffers to load data
    int *intbuff=new int[chunksize];
    long long *longbuff=new long long[chunksize];
    unsigned int *uintbuff=new unsigned int[chunksize];
    float *floatbuff=new float[chunksize*3];
    double *doublebuff=new double[chunksize*3];
    void *integerbuff,*realbuff;
    vector<double> vdoublebuff;
    vector<int> vintbuff;
    vector<unsigned int> vuintbuff;
    vector<long long> vlongbuff;
    vector<unsigned long long> vulongbuff;
    //arrays to store number of items to read and offsets when selecting hyperslabs
    hsize_t filespacecount[HDFMAXPROPDIM],filespaceoffset[HDFMAXPROPDIM];
    //to determine types
    //IntType inttype;
    //FloatType floattype;
    //PredType HDFREALTYPE(PredType::NATIVE_FLOAT);
    //PredType HDFINTEGERTYPE(PredType::NATIVE_LONG);
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
    Int_t i,j,k,n,nchunk,count,bcount,itemp,count2,bcount2;

    //store cosmology
    double z,aadjust,Hubble,Hubbleflow;

    double mscale,lscale,lvscale;
    double MP_DM=MAXVALUE,LN,N_DM,MP_B=0;
    int ifirstfile=0,*ireadfile,ireaderror=0;
    int *ireadtask,*readtaskID;
    Int_t ninputoffset;

#ifdef USEMPI
    if (ThisTask == 0)
#endif
    opt.num_files = HDF_get_nfiles (opt.fname, opt.partsearchtype);
    HDFSetUsedParticleTypes(opt,nusetypes,nbusetypes,usetypes);

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
    // DataSet *partsdatasetall;
    // DataSpace *partsdataspaceall;
    vector<hid_t> partsdatasetall;
    vector<hid_t> partsdataspaceall;

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
    hdf_header_info=new HDF_Header[opt.num_files];
    for (i=0; i<opt.num_files; i++) hdf_header_info[i] = HDF_Header(opt.ihdfnameconvention);
    // Fhdf=new H5File[opt.num_files];
    // headerdataspace=new DataSpace[opt.num_files];
    // headerattribs=new Attribute[opt.num_files];
    // partsgroup=new Group[opt.num_files*NHDFTYPE];
    // partsdataset=new DataSet[opt.num_files*NHDFTYPE];
    // partsdataspace=new DataSpace[opt.num_files*NHDFTYPE];
    Fhdf.resize(opt.num_files);
    headerdataspace.resize(opt.num_files);
    headerattribs.resize(opt.num_files);
    partsgroup.resize(opt.num_files*NHDFTYPE);
    partsdataset.resize(opt.num_files*NHDFTYPE);
    partsdataspace.resize(opt.num_files*NHDFTYPE);
    //init the hid_t to negative values
    for (auto &x:Fhdf) x=-1;
    for (auto &x:headerdataspace) x=-1;
    for (auto &x:headerattribs) x=-1;
    for (auto &x:partsgroup) x=-1;
    for (auto &x:partsdataset) x=-1;
    for (auto &x:partsdataspace) x=-1;
#ifdef USEMPI
    partsdatasetall.resize(opt.num_files*NHDFTYPE*NHDFDATABLOCK);
    partsdataspaceall.resize(opt.num_files*NHDFTYPE*NHDFDATABLOCK);
    for (auto &x:partsdatasetall) x=-1;
    for (auto &x:partsdataspaceall) x=-1;
#endif
    for(i=0; i<opt.num_files; i++) if(ireadfile[i]) {
        if(opt.num_files>1) sprintf(buf,"%s.%d.hdf5",opt.fname,(int)i);
        else sprintf(buf,"%s.hdf5",opt.fname);
        //Try block to detect exceptions raised by any of the calls inside it
        //try
        {
            //turn off the auto-printing when failure occurs so that we can
            //handle the errors appropriately
            //Exception::dontPrint();

            //Open the specified file and the specified dataset in the file.
            Fhdf[i] = H5Fopen(buf, H5F_ACC_RDONLY, H5P_DEFAULT);
            if (ThisTask==0 && i==0) {
                cout<<buf<<endl;
                cout<<"HDF file contains the following group structures "<<endl;
                //H5Literate(Fhdf[i].getId(), H5_INDEX_NAME, H5_ITER_INC, NULL, file_info, NULL);
                H5Literate(Fhdf[i], H5_INDEX_NAME, H5_ITER_INC, NULL, file_info, NULL);
                cout<<" Expecting "<<endl;
                for (j=0;j<NHDFTYPE+1;j++) cout<<hdf_gnames.names[j]<<endl;
            }
            hdf_header_info[i].BoxSize = read_attribute<double>(Fhdf[i], hdf_header_info[i].names[hdf_header_info[i].IBoxSize]);
            vdoublebuff=read_attribute_v<double>(Fhdf[i], hdf_header_info[i].names[hdf_header_info[i].IMass]);
            for (k=0;k<NHDFTYPE;k++)hdf_header_info[i].mass[k]=vdoublebuff[k];
            if (opt.ihdfnameconvention==HDFSWIFTEAGLENAMES) {
                vlongbuff = read_attribute_v<long long>(Fhdf[i], hdf_header_info[i].names[hdf_header_info[i].INuminFile]);
                for (k=0;k<NHDFTYPE;k++) hdf_header_info[i].npart[k]=vlongbuff[k];
            }
            else{
                vuintbuff = read_attribute_v<unsigned int>(Fhdf[i], hdf_header_info[i].names[hdf_header_info[i].INuminFile]);
                for (k=0;k<NHDFTYPE;k++) hdf_header_info[i].npart[k]=vuintbuff[k];
            }
            vuintbuff=read_attribute_v<unsigned int>(Fhdf[i], hdf_header_info[i].names[hdf_header_info[i].INumTot]);
            for (k=0;k<NHDFTYPE;k++)hdf_header_info[i].npartTotal[k]=vuintbuff[k];
            vuintbuff=read_attribute_v<unsigned int>(Fhdf[i], hdf_header_info[i].names[hdf_header_info[i].INumTotHW]);
            for (k=0;k<NHDFTYPE;k++)hdf_header_info[i].npartTotalHW[k]=vuintbuff[k];
            hdf_header_info[i].Omega0 = read_attribute<double>(Fhdf[i], hdf_header_info[i].names[hdf_header_info[i].IOmega0]);
            hdf_header_info[i].OmegaLambda = read_attribute<double>(Fhdf[i], hdf_header_info[i].names[hdf_header_info[i].IOmegaL]);
            hdf_header_info[i].redshift = read_attribute<double>(Fhdf[i], hdf_header_info[i].names[hdf_header_info[i].IRedshift]);
            hdf_header_info[i].time = read_attribute<double>(Fhdf[i], hdf_header_info[i].names[hdf_header_info[i].ITime]);
            hdf_header_info[i].HubbleParam = read_attribute<double>(Fhdf[i], hdf_header_info[i].names[hdf_header_info[i].IHubbleParam]);
            hdf_header_info[i].num_files = read_attribute<int>(Fhdf[i], hdf_header_info[i].names[hdf_header_info[i].INumFiles]);
        }
        /*catch(GroupIException &error)
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
        */
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
      mscale=opt.massinputconversion;lscale=opt.lengthinputconversion*aadjust;lvscale=opt.lengthinputconversion*opt.a;
    }
    else {
      mscale=opt.massinputconversion/opt.h;lscale=opt.lengthinputconversion/opt.h*aadjust;lvscale=opt.lengthinputconversion/opt.h*opt.a;
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
        //try
        {
            //open particle group structures
            for (j=0;j<nusetypes;j++) {
              k=usetypes[j];
              if (ThisTask==0 && opt.iverbose>1) cout<<"Opening group "<<hdf_gnames.part_names[k]<<endl;
              // partsgroup[i*NHDFTYPE+k]=Fhdf[i].openGroup(hdf_gnames.part_names[k]);
              partsgroup[i*NHDFTYPE+k]=HDF5OpenGroup(Fhdf[i],hdf_gnames.part_names[k]);
            }
            if (opt.partsearchtype==PSTDARK && opt.iBaryonSearch) for (j=1;j<=nbusetypes;j++) {
              k=usetypes[j];
              if (ThisTask==0 && opt.iverbose>1) cout<<"Opening group "<<hdf_gnames.part_names[k]<<endl;
              // partsgroup[i*NHDFTYPE+k]=Fhdf[i].openGroup(hdf_gnames.part_names[k]);
              partsgroup[i*NHDFTYPE+k]=HDF5OpenGroup(Fhdf[i],hdf_gnames.part_names[k]);
            }
            itemp=0;
            //get positions
            for (j=0;j<nusetypes;j++) {
              k=usetypes[j];
              if (ThisTask==0 && opt.iverbose>1) cout<<"Opening group "<<hdf_gnames.part_names[k]<<": Data set "<<hdf_parts[k]->names[itemp]<<endl;
              partsdataset[i*NHDFTYPE+k]=HDF5OpenDataSet(partsgroup[i*NHDFTYPE+k],hdf_parts[k]->names[0]);
              partsdataspace[i*NHDFTYPE+k]=HDF5OpenDataSpace(partsdataset[i*NHDFTYPE+k]);
            }
            if (opt.partsearchtype==PSTDARK && opt.iBaryonSearch) for (j=1;j<=nbusetypes;j++) {
              k=usetypes[j];
              partsdataset[i*NHDFTYPE+k]=HDF5OpenDataSet(partsgroup[i*NHDFTYPE+k],hdf_parts[k]->names[0]);
              partsdataspace[i*NHDFTYPE+k]=HDF5OpenDataSpace(partsdataset[i*NHDFTYPE+k]);
            }
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
                HDF5ReadHyperSlabReal(doublebuff,partsdataset[i*NHDFTYPE+k], partsdataspace[i*NHDFTYPE+k], 1, 3, nchunk, n);
                for (int nn=0;nn<nchunk;nn++) Part[count++].SetPosition(doublebuff[nn*3],doublebuff[nn*3+1],doublebuff[nn*3+2]);
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
                  // //setup hyperslab so that it is loaded into the buffer
                  HDF5ReadHyperSlabReal(doublebuff,partsdataset[i*NHDFTYPE+k], partsdataspace[i*NHDFTYPE+k], 1, 3, nchunk, n);
                  for (int nn=0;nn<nchunk;nn++) Pbaryons[bcount++].SetPosition(doublebuff[nn*3],doublebuff[nn*3+1],doublebuff[nn*3+2]);
                }
              }
            }
            //close data spaces
            for (auto &hidval:partsdataspace) HDF5CloseDataSpace(hidval);
            for (auto &hidval:partsdataset) HDF5CloseDataSet(hidval);
            //get velocities
            itemp++;
            for (j=0;j<nusetypes;j++) {
              k=usetypes[j];
              if (ThisTask==0 && opt.iverbose>1) cout<<"Opening group "<<hdf_gnames.part_names[k]<<": Data set "<<hdf_parts[k]->names[itemp]<<endl;
              partsdataset[i*NHDFTYPE+k]=HDF5OpenDataSet(partsgroup[i*NHDFTYPE+k],hdf_parts[k]->names[itemp]);
              partsdataspace[i*NHDFTYPE+k]=HDF5OpenDataSpace(partsdataset[i*NHDFTYPE+k]);
            }
            if (opt.partsearchtype==PSTDARK && opt.iBaryonSearch) for (j=1;j<=nbusetypes;j++) {
              k=usetypes[j];
              partsdataset[i*NHDFTYPE+k]=HDF5OpenDataSet(partsgroup[i*NHDFTYPE+k],hdf_parts[k]->names[itemp]);
              partsdataspace[i*NHDFTYPE+k]=HDF5OpenDataSpace(partsdataset[i*NHDFTYPE+k]);
            }
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
                HDF5ReadHyperSlabReal(doublebuff,partsdataset[i*NHDFTYPE+k], partsdataspace[i*NHDFTYPE+k], 1, 3, nchunk, n);
                for (int nn=0;nn<nchunk;nn++) Part[count++].SetVelocity(doublebuff[nn*3],doublebuff[nn*3+1],doublebuff[nn*3+2]);
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
                  // //setup hyperslab so that it is loaded into the buffer
                  HDF5ReadHyperSlabReal(doublebuff,partsdataset[i*NHDFTYPE+k], partsdataspace[i*NHDFTYPE+k], 1, 3, nchunk, n);
                  for (int nn=0;nn<nchunk;nn++) Pbaryons[bcount++].SetVelocity(doublebuff[nn*3],doublebuff[nn*3+1],doublebuff[nn*3+2]);
                }
              }
            }
            //close data spaces
            for (auto &hidval:partsdataspace) HDF5CloseDataSpace(hidval);
            for (auto &hidval:partsdataset) HDF5CloseDataSet(hidval);
            //get ids
            itemp++;
            for (j=0;j<nusetypes;j++) {
              k=usetypes[j];
              if (ThisTask==0 && opt.iverbose>1) cout<<"Opening group "<<hdf_gnames.part_names[k]<<": Data set "<<hdf_parts[k]->names[itemp]<<endl;
              partsdataset[i*NHDFTYPE+k]=HDF5OpenDataSet(partsgroup[i*NHDFTYPE+k],hdf_parts[k]->names[itemp]);
              partsdataspace[i*NHDFTYPE+k]=HDF5OpenDataSpace(partsdataset[i*NHDFTYPE+k]);
            }
            if (opt.partsearchtype==PSTDARK && opt.iBaryonSearch) for (j=1;j<=nbusetypes;j++) {
              k=usetypes[j];
              partsdataset[i*NHDFTYPE+k]=HDF5OpenDataSet(partsgroup[i*NHDFTYPE+k],hdf_parts[k]->names[itemp]);
              partsdataspace[i*NHDFTYPE+k]=HDF5OpenDataSpace(partsdataset[i*NHDFTYPE+k]);
            }
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
                HDF5ReadHyperSlabInteger(longbuff,partsdataset[i*NHDFTYPE+k], partsdataspace[i*NHDFTYPE+k], 1, 1, nchunk, n);

                for (int nn=0;nn<nchunk;nn++) {
                    Part[count].SetPID(longbuff[nn]);
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
                  HDF5ReadHyperSlabInteger(longbuff,partsdataset[i*NHDFTYPE+k], partsdataspace[i*NHDFTYPE+k], 1, 1, nchunk, n);

                  for (int nn=0;nn<nchunk;nn++) {
                    Pbaryons[bcount].SetPID(longbuff[nn]);
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
            //close data spaces
            for (auto &hidval:partsdataspace) HDF5CloseDataSpace(hidval);
            for (auto &hidval:partsdataset) HDF5CloseDataSet(hidval);

            //get masses, note that DM do not contain a mass field
            itemp++;
            for (j=0;j<nusetypes;j++) {
              k=usetypes[j];
              if (ThisTask==0 && opt.iverbose>1) cout<<"Opening group "<<hdf_gnames.part_names[k]<<": Data set "<<hdf_parts[k]->names[itemp]<<endl;
              if (hdf_header_info[i].mass[k]==0){
                partsdataset[i*NHDFTYPE+k]=HDF5OpenDataSet(partsgroup[i*NHDFTYPE+k],hdf_parts[k]->names[itemp]);
                partsdataspace[i*NHDFTYPE+k]=HDF5OpenDataSpace(partsdataset[i*NHDFTYPE+k]);
              }
            }
            if (opt.partsearchtype==PSTDARK && opt.iBaryonSearch) for (j=1;j<=nbusetypes;j++) {
              k=usetypes[j];
              if (hdf_header_info[i].mass[k]==0){
                partsdataset[i*NHDFTYPE+k]=HDF5OpenDataSet(partsgroup[i*NHDFTYPE+k],hdf_parts[k]->names[itemp]);
                partsdataspace[i*NHDFTYPE+k]=HDF5OpenDataSpace(partsdataset[i*NHDFTYPE+k]);
              }
            }
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
                  HDF5ReadHyperSlabReal(doublebuff,partsdataset[i*NHDFTYPE+k], partsdataspace[i*NHDFTYPE+k], 1, 1, nchunk, n);
                  for (int nn=0;nn<nchunk;nn++) Part[count++].SetMass(doublebuff[nn]);
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
                    HDF5ReadHyperSlabReal(doublebuff,partsdataset[i*NHDFTYPE+k], partsdataspace[i*NHDFTYPE+k], 1, 1, nchunk, n);
                    for (int nn=0;nn<nchunk;nn++) Pbaryons[bcount++].SetMass(doublebuff[nn]);
                  }
                }
                else {
                  for (int nn=0;nn<hdf_header_info[i].npart[k];nn++) Pbaryons[bcount++].SetMass(hdf_header_info[i].mass[k]);
                }
              }
            }
            //close data spaces
            for (auto &hidval:partsdataspace) HDF5CloseDataSpace(hidval);
            for (auto &hidval:partsdataset) HDF5CloseDataSet(hidval);

            //and if not just searching DM, load other parameters
            if (!(opt.partsearchtype==PSTDARK && opt.iBaryonSearch==0)) {
#ifdef GASON
              //first gas internal energy
              for (j=0;j<nusetypes;j++) {
                k=usetypes[j];
                if (k==HDFGASTYPE){
                  if (ThisTask==0 && opt.iverbose>1) cout<<"Opening group "<<hdf_gnames.part_names[k]<<": Data set "<<hdf_parts[k]->names[5]<<endl;
                  partsdataset[i*NHDFTYPE+k]=HDF5OpenDataSet(partsgroup[i*NHDFTYPE+k],hdf_parts[k]->names[5]);
                  partsdataspace[i*NHDFTYPE+k]=HDF5OpenDataSpace(partsdataset[i*NHDFTYPE+k]);
                }
              }
              if (opt.partsearchtype==PSTDARK && opt.iBaryonSearch) for (j=1;j<=nbusetypes;j++) {
                k=usetypes[j];
                if (k==HDFGASTYPE){
                  partsdataset[i*NHDFTYPE+k]=HDF5OpenDataSet(partsgroup[i*NHDFTYPE+k],hdf_parts[k]->names[5]);
                  partsdataspace[i*NHDFTYPE+k]=HDF5OpenDataSpace(partsdataset[i*NHDFTYPE+k]);
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
                    // datarank=1;
                    // datadim[0]=nchunk;
                    // chunkspace=DataSpace(datarank,datadim);
                    // filespacecount[0]=nchunk;filespacecount[1]=1;
                    // filespaceoffset[0]=n;filespaceoffset[1]=0;
                    // partsdataspace[i*NHDFTYPE+k].selectHyperslab(H5S_SELECT_SET, filespacecount, filespaceoffset);
                    // partsdataset[i*NHDFTYPE+k].read(realbuff,HDFREALTYPE,chunkspace,partsdataspace[i*NHDFTYPE+k]);
                    //
                    // if (ifloat) for (int nn=0;nn<nchunk;nn++) Part[count++].SetU(floatbuff[nn]);
                    // else for (int nn=0;nn<nchunk;nn++) Part[count++].SetU(doublebuff[nn]);
                    HDF5ReadHyperSlabReal(doublebuff,partsdataset[i*NHDFTYPE+k], partsdataspace[i*NHDFTYPE+k], 1, 1, nchunk, n);
                    for (int nn=0;nn<nchunk;nn++) Part[count++].SetU(doublebuff[nn]);
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
                      // datarank=1;
                      // datadim[0]=nchunk;
                      // chunkspace=DataSpace(datarank,datadim);
                      // filespacecount[0]=nchunk;filespacecount[1]=1;
                      // filespaceoffset[0]=n;filespaceoffset[1]=0;
                      // partsdataspace[i*NHDFTYPE+k].selectHyperslab(H5S_SELECT_SET, filespacecount, filespaceoffset);
                      // partsdataset[i*NHDFTYPE+k].read(realbuff,HDFREALTYPE,chunkspace,partsdataspace[i*NHDFTYPE+k]);
                      //
                      // if (ifloat) for (int nn=0;nn<nchunk;nn++) Pbaryons[bcount++].SetU(floatbuff[nn]);
                      // else for (int nn=0;nn<nchunk;nn++) Pbaryons[bcount++].SetU(doublebuff[nn]);
                      HDF5ReadHyperSlabReal(doublebuff,partsdataset[i*NHDFTYPE+k], partsdataspace[i*NHDFTYPE+k], 1, 1, nchunk, n);
                      for (int nn=0;nn<nchunk;nn++) Pbaryons[bcount++].SetU(doublebuff[nn]);
                    }
                  }
                  else {
                    count+=hdf_header_info[i].npart[k];
                  }
                }
              }
              //close data spaces
              for (auto &hidval:partsdataspace) HDF5CloseDataSpace(hidval);
              for (auto &hidval:partsdataset) HDF5CloseDataSet(hidval);
#ifdef STARON
              //if star forming get star formation rate
              for (j=0;j<nusetypes;j++) {
                k=usetypes[j];
                if (k==HDFGASTYPE){
                  if (ThisTask==0 && opt.iverbose>1) cout<<"Opening group "<<hdf_gnames.part_names[k]<<": Data set "<<hdf_parts[k]->names[6]<<endl;
                  partsdataset[i*NHDFTYPE+k]=HDF5OpenDataSet(partsgroup[i*NHDFTYPE+k],hdf_parts[k]->names[6]);
                  partsdataspace[i*NHDFTYPE+k]=HDF5OpenDataSpace(partsdataset[i*NHDFTYPE+k]);
                }
              }
              if (opt.partsearchtype==PSTDARK && opt.iBaryonSearch) for (j=1;j<=nbusetypes;j++) {
                k=usetypes[j];
                if (k==HDFGASTYPE){
                  partsdataset[i*NHDFTYPE+k]=HDF5OpenDataSet(partsgroup[i*NHDFTYPE+k],hdf_parts[k]->names[6]);
                  partsdataspace[i*NHDFTYPE+k]=HDF5OpenDataSpace(partsdataset[i*NHDFTYPE+k]);
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
                    HDF5ReadHyperSlabReal(doublebuff,partsdataset[i*NHDFTYPE+k], partsdataspace[i*NHDFTYPE+k], 1, 1, nchunk, n);
                    for (int nn=0;nn<nchunk;nn++) Part[count++].SetSFR(doublebuff[nn]);
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
                      HDF5ReadHyperSlabReal(doublebuff,partsdataset[i*NHDFTYPE+k], partsdataspace[i*NHDFTYPE+k], 1, 1, nchunk, n);
                      for (int nn=0;nn<nchunk;nn++) Pbaryons[bcount++].SetU(doublebuff[nn]);
                    }
                  }
                  else {
                    count+=hdf_header_info[i].npart[k];
                  }
                }
              }
              //close data spaces
              for (auto &hidval:partsdataspace) HDF5CloseDataSpace(hidval);
              for (auto &hidval:partsdataset) HDF5CloseDataSet(hidval);
              //then metallicity
              for (j=0;j<nusetypes;j++) {
                k=usetypes[j];
                if (k==HDFGASTYPE){
                  if (ThisTask==0 && opt.iverbose>1) cout<<"Opening group "<<hdf_gnames.part_names[k]<<": Data set "<<hdf_parts[k]->names[hdf_parts[k]->propindex[HDFGASIMETAL]]<<endl;
                  partsdataset[i*NHDFTYPE+k]=HDF5OpenDataSet(partsgroup[i*NHDFTYPE+k],hdf_parts[k]->names[hdf_parts[k]->propindex[HDFGASIMETAL]]);
                  partsdataspace[i*NHDFTYPE+k]=HDF5OpenDataSpace(partsdataset[i*NHDFTYPE+k]);
                }
                if (k==HDFSTARTYPE){
                  if (ThisTask==0 && opt.iverbose>1) cout<<"Opening group "<<hdf_gnames.part_names[k]<<": Data set "<<hdf_parts[k]->names[hdf_parts[k]->propindex[HDFSTARIMETAL]]<<endl;
                  partsdataset[i*NHDFTYPE+k]=HDF5OpenDataSet(partsgroup[i*NHDFTYPE+k],hdf_parts[k]->names[hdf_parts[k]->propindex[HDFSTARIMETAL]]);
                  partsdataspace[i*NHDFTYPE+k]=HDF5OpenDataSpace(partsdataset[i*NHDFTYPE+k]);
                }
              }
              if (opt.partsearchtype==PSTDARK && opt.iBaryonSearch) for (j=1;j<=nbusetypes;j++) {
                k=usetypes[j];
                if (k==HDFGASTYPE){
                  partsdataset[i*NHDFTYPE+k]=HDF5OpenDataSet(partsgroup[i*NHDFTYPE+k],hdf_parts[k]->names[hdf_parts[k]->propindex[HDFGASIMETAL]]);
                  partsdataspace[i*NHDFTYPE+k]=HDF5OpenDataSpace(partsdataset[i*NHDFTYPE+k]);
                }
                if (k==HDFSTARTYPE){
                  partsdataset[i*NHDFTYPE+k]=HDF5OpenDataSet(partsgroup[i*NHDFTYPE+k],hdf_parts[k]->names[hdf_parts[k]->propindex[HDFSTARIMETAL]]);
                  partsdataspace[i*NHDFTYPE+k]=HDF5OpenDataSpace(partsdataset[i*NHDFTYPE+k]);
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
                    HDF5ReadHyperSlabReal(doublebuff,partsdataset[i*NHDFTYPE+k], partsdataspace[i*NHDFTYPE+k], 1, 1, nchunk, n);
                    for (int nn=0;nn<nchunk;nn++) Part[count++].SetZmet(doublebuff[nn]*zmetconversion);
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
                      HDF5ReadHyperSlabReal(doublebuff,partsdataset[i*NHDFTYPE+k], partsdataspace[i*NHDFTYPE+k], 1, 1, nchunk, n);
                      Pbaryons[bcount++].SetZmet(doublebuff[nn]*zmetconversion);
                    }
                  }
                  else {
                    count+=hdf_header_info[i].npart[k];
                  }
                }
              }
              //close data spaces
              for (auto &hidval:partsdataspace) HDF5CloseDataSpace(hidval);
              for (auto &hidval:partsdataset) HDF5CloseDataSet(hidval);
              //then get star formation time, must also adjust so that if tage<0 this is a wind particle in Illustris so change particle type
              for (j=0;j<nusetypes;j++) {
                k=usetypes[j];
                if (k==HDFSTARTYPE){
                  if (ThisTask==0 && opt.iverbose>1) cout<<"Opening group "<<hdf_gnames.part_names[k]<<": Data set "<<hdf_parts[k]->names[hdf_parts[k]->propindex[HDFSTARIAGE]]<<endl;
                  partsdataset[i*NHDFTYPE+k]=HDF5OpenDataSet(partsgroup[i*NHDFTYPE+k],hdf_parts[k]->names[hdf_parts[k]->propindex[HDFSTARIAGE]]);
                  partsdataspace[i*NHDFTYPE+k]=HDF5OpenDataSpace(partsdataset[i*NHDFTYPE+k]);
                }
              }
              if (opt.partsearchtype==PSTDARK && opt.iBaryonSearch) for (j=1;j<=nbusetypes;j++) {
                k=usetypes[j];
                if (k==HDFSTARTYPE){
                  partsdataset[i*NHDFTYPE+k]=HDF5OpenDataSet(partsgroup[i*NHDFTYPE+k],hdf_parts[k]->names[hdf_parts[k]->propindex[HDFSTARIAGE]]);
                  partsdataspace[i*NHDFTYPE+k]=HDF5OpenDataSpace(partsdataset[i*NHDFTYPE+k]);
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
                    HDF5ReadHyperSlabReal(doublebuff,partsdataset[i*NHDFTYPE+k], partsdataspace[i*NHDFTYPE+k], 1, 1, nchunk, n);
                    for (int nn=0;nn<nchunk;nn++) {if (doublebuff[nn]<0) Part[count].SetType(WINDTYPE);Part[count++].SetTage(doublebuff[nn]);}
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
                      HDF5ReadHyperSlabReal(doublebuff,partsdataset[i*NHDFTYPE+k], partsdataspace[i*NHDFTYPE+k], 1, 1, nchunk, n);
                      for (int nn=0;nn<nchunk;nn++) {if (doublebuff[nn]<0) Pbaryons[bcount].SetType(WINDTYPE);Pbaryons[bcount++].SetTage(doublebuff[nn]);}
                    }
                  }
                  else {
                    count+=hdf_header_info[i].npart[k];
                  }
                }
              }
              //close data spaces
              for (auto &hidval:partsdataspace) HDF5CloseDataSpace(hidval);
              for (auto &hidval:partsdataset) HDF5CloseDataSet(hidval);
#endif
#endif
            }//end of if not dark matter then baryon search
            //close groups
            for (auto &hidval:partsgroup) HDF5CloseGroup(hidval);
            count2=count;
            bcount2=bcount;
        }//end of try
        /*
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
        */
        //Fhdf[i].close();
        HDF5CloseFile(Fhdf[i]);
    }

    double vscale = 0.0;

    // SWIFT snapshot velocities already contain the sqrt(a) factor,
    // so there is no need to include it.
    if(opt.ihdfnameconvention == HDFSWIFTEAGLENAMES) vscale = opt.velocityinputconversion;
    else vscale = opt.velocityinputconversion*sqrt(opt.a);

    //finally adjust to appropriate units
    for (i=0;i<nbodies;i++)
    {
#ifdef HIGHRES
      if (Part[i].GetType()==DARKTYPE && Part[i].GetMass()<MP_DM) MP_DM=Part[i].GetMass();
#endif
      Part[i].SetMass(Part[i].GetMass()*mscale);
      for (int j=0;j<3;j++) Part[i].SetVelocity(j,Part[i].GetVelocity(j)*vscale+Hubbleflow*Part[i].GetPosition(j));
      for (int j=0;j<3;j++) Part[i].SetPosition(j,Part[i].GetPosition(j)*lscale);
#ifdef GASON
      if (Part[i].GetType()==GASTYPE) Part[i].SetU(Part[i].GetU()*opt.velocityinputconversion*opt.velocityinputconversion);
#endif
    }
    if (Pbaryons!=NULL && opt.iBaryonSearch==1) {
      for (i=0;i<nbaryons;i++)
      {
        Pbaryons[i].SetMass(Pbaryons[i].GetMass()*mscale);
        for (int j=0;j<3;j++) Pbaryons[i].SetVelocity(j,Pbaryons[i].GetVelocity(j)*vscale+Hubbleflow*Pbaryons[i].GetPosition(j));
        for (int j=0;j<3;j++) Pbaryons[i].SetPosition(j,Pbaryons[i].GetPosition(j)*lscale);
#ifdef GASON
        Pbaryons[i].SetU(Pbaryons[i].GetU()*opt.velocityinputconversion*opt.velocityinputconversion);
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
            //try
            {
                //open particle group structures
                for (j=0;j<nusetypes;j++) {
                    k=usetypes[j];
                    // partsgroup[i*NHDFTYPE+k]=Fhdf[i].openGroup(hdf_gnames.part_names[k]);
                    partsgroup[i*NHDFTYPE+k]=HDF5OpenGroup(Fhdf[i],hdf_gnames.part_names[k]);
                }
                if (opt.partsearchtype==PSTDARK && opt.iBaryonSearch) {
                    for (j=1;j<=nbusetypes;j++) {
                        k=usetypes[j];
                        partsgroup[i*NHDFTYPE+k]=HDF5OpenGroup(Fhdf[i],hdf_gnames.part_names[k]);
                    }
                }
                //open data structures that exist for all data blocks
                for (itemp=0;itemp<NHDFDATABLOCKALL;itemp++) {
                  //for everything but mass no header check needed.
                  if (itemp!=3) {
                    for (j=0;j<nusetypes;j++) {
                      k=usetypes[j];
                      if (ThisTask==0 && opt.iverbose>1) cout<<"Opening group "<<hdf_gnames.part_names[k]<<": Data set "<<hdf_parts[k]->names[itemp]<<endl;
                      partsdatasetall[i*NHDFTYPE*NHDFDATABLOCK+k*NHDFDATABLOCK+itemp]=HDF5OpenDataSet(partsgroup[i*NHDFTYPE+k],hdf_parts[k]->names[itemp]);
                      partsdataspaceall[i*NHDFTYPE*NHDFDATABLOCK+k*NHDFDATABLOCK+itemp]=HDF5OpenDataSpace(partsdatasetall[i*NHDFTYPE*NHDFDATABLOCK+k*NHDFDATABLOCK+itemp]);
                    }
                    if (opt.partsearchtype==PSTDARK && opt.iBaryonSearch) for (j=1;j<=nbusetypes;j++) {
                      k=usetypes[j];
                      partsdatasetall[i*NHDFTYPE*NHDFDATABLOCK+k*NHDFDATABLOCK+itemp]=HDF5OpenDataSet(partsgroup[i*NHDFTYPE+k],hdf_parts[k]->names[itemp]);
                      partsdataspaceall[i*NHDFTYPE*NHDFDATABLOCK+k*NHDFDATABLOCK+itemp]=HDF5OpenDataSpace(partsdatasetall[i*NHDFTYPE*NHDFDATABLOCK+k*NHDFDATABLOCK+itemp]);
                    }
                  }
                  else {
                    for (j=0;j<nusetypes;j++) {
                      k=usetypes[j];
                      if (hdf_header_info[i].mass[k]==0){
                        if (ThisTask==0 && opt.iverbose>1) cout<<"Opening group "<<hdf_gnames.part_names[k]<<": Data set "<<hdf_parts[k]->names[itemp]<<endl;
                        partsdatasetall[i*NHDFTYPE*NHDFDATABLOCK+k*NHDFDATABLOCK+itemp]=HDF5OpenDataSet(partsgroup[i*NHDFTYPE+k],hdf_parts[k]->names[itemp]);
                        partsdataspaceall[i*NHDFTYPE*NHDFDATABLOCK+k*NHDFDATABLOCK+itemp]=HDF5OpenDataSpace(partsdatasetall[i*NHDFTYPE*NHDFDATABLOCK+k*NHDFDATABLOCK+itemp]);
                      }
                    }
                    if (opt.partsearchtype==PSTDARK && opt.iBaryonSearch) for (j=1;j<=nbusetypes;j++) {
                      k=usetypes[j];
                      partsdatasetall[i*NHDFTYPE*NHDFDATABLOCK+k*NHDFDATABLOCK+itemp]=HDF5OpenDataSet(partsgroup[i*NHDFTYPE+k],hdf_parts[k]->names[itemp]);
                      partsdataspaceall[i*NHDFTYPE*NHDFDATABLOCK+k*NHDFDATABLOCK+itemp]=HDF5OpenDataSpace(partsdatasetall[i*NHDFTYPE*NHDFDATABLOCK+k*NHDFDATABLOCK+itemp]);
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
                      partsdatasetall[i*NHDFTYPE*NHDFDATABLOCK+k*NHDFDATABLOCK+itemp]=HDF5OpenDataSet(partsgroup[i*NHDFTYPE+k],hdf_parts[k]->names[5]);
                      partsdataspaceall[i*NHDFTYPE*NHDFDATABLOCK+k*NHDFDATABLOCK+itemp]=HDF5OpenDataSpace(partsdatasetall[i*NHDFTYPE*NHDFDATABLOCK+k*NHDFDATABLOCK+itemp]);
                    }
                  }
                  if (opt.partsearchtype==PSTDARK && opt.iBaryonSearch) for (j=1;j<=nbusetypes;j++) {
                    k=usetypes[j];
                    if (k==HDFGASTYPE){
                      partsdatasetall[i*NHDFTYPE*NHDFDATABLOCK+k*NHDFDATABLOCK+itemp]=HDF5OpenDataSet(partsgroup[i*NHDFTYPE+k],hdf_parts[k]->names[5]);
                      partsdataspaceall[i*NHDFTYPE*NHDFDATABLOCK+k*NHDFDATABLOCK+itemp]=HDF5OpenDataSpace(partsdatasetall[i*NHDFTYPE*NHDFDATABLOCK+k*NHDFDATABLOCK+itemp]);
                    }
                  }
#ifdef STARON
                  //if star forming get star formation rate
                  itemp++;
                  for (j=0;j<nusetypes;j++) {
                    k=usetypes[j];
                    if (k==HDFGASTYPE){
                      if (ThisTask==0 && opt.iverbose>1) cout<<"Opening group "<<hdf_gnames.part_names[k]<<": Data set "<<hdf_parts[k]->names[6]<<endl;
                      partsdatasetall[i*NHDFTYPE*NHDFDATABLOCK+k*NHDFDATABLOCK+itemp]=HDF5OpenDataSet(partsgroup[i*NHDFTYPE+k],hdf_parts[k]->names[6]);
                      partsdataspaceall[i*NHDFTYPE*NHDFDATABLOCK+k*NHDFDATABLOCK+itemp]=HDF5OpenDataSpace(partsdatasetall[i*NHDFTYPE*NHDFDATABLOCK+k*NHDFDATABLOCK+itemp]);
                    }
                  }
                  if (opt.partsearchtype==PSTDARK && opt.iBaryonSearch) for (j=1;j<=nbusetypes;j++) {
                    k=usetypes[j];
                    if (k==HDFGASTYPE){
                      partsdatasetall[i*NHDFTYPE*NHDFDATABLOCK+k*NHDFDATABLOCK+itemp]=HDF5OpenDataSet(partsgroup[i*NHDFTYPE+k],hdf_parts[k]->names[6]);
                      partsdataspaceall[i*NHDFTYPE*NHDFDATABLOCK+k*NHDFDATABLOCK+itemp]=HDF5OpenDataSpace(partsdatasetall[i*NHDFTYPE*NHDFDATABLOCK+k*NHDFDATABLOCK+itemp]);
                    }
                  }
                  //then metallicity
                  itemp++;
                  for (j=0;j<nusetypes;j++) {
                    k=usetypes[j];
                    if (k==HDFGASTYPE){
                      if (ThisTask==0 && opt.iverbose>1) cout<<"Opening group "<<hdf_gnames.part_names[k]<<": Data set "<<hdf_parts[k]->names[hdf_parts[k]->propindex[HDFGASIMETAL]]<<endl;
                      partsdatasetall[i*NHDFTYPE*NHDFDATABLOCK+k*NHDFDATABLOCK+itemp]=HDF5OpenDataSet(partsgroup[i*NHDFTYPE+k],hdf_parts[k]->names[hdf_parts[k]->propindex[HDFGASIMETAL]]);
                      partsdataspaceall[i*NHDFTYPE*NHDFDATABLOCK+k*NHDFDATABLOCK+itemp]=HDF5OpenDataSpace(partsdatasetall[i*NHDFTYPE*NHDFDATABLOCK+k*NHDFDATABLOCK+itemp]);
                    }
                    if (k==HDFSTARTYPE){
                      if (ThisTask==0 && opt.iverbose>1) cout<<"Opening group "<<hdf_gnames.part_names[k]<<": Data set "<<hdf_parts[k]->names[hdf_parts[k]->propindex[HDFSTARIMETAL]]<<endl;
                      partsdatasetall[i*NHDFTYPE*NHDFDATABLOCK+k*NHDFDATABLOCK+itemp]=HDF5OpenDataSet(partsgroup[i*NHDFTYPE+k],hdf_parts[k]->names[hdf_parts[k]->propindex[HDFSTARIMETAL]]);
                      partsdataspaceall[i*NHDFTYPE*NHDFDATABLOCK+k*NHDFDATABLOCK+itemp]=HDF5OpenDataSpace(partsdatasetall[i*NHDFTYPE*NHDFDATABLOCK+k*NHDFDATABLOCK+itemp]);
                    }
                  }
                  if (opt.partsearchtype==PSTDARK && opt.iBaryonSearch) for (j=1;j<=nbusetypes;j++) {
                    k=usetypes[j];
                    if (k==HDFGASTYPE){
                      partsdatasetall[i*NHDFTYPE*NHDFDATABLOCK+k*NHDFDATABLOCK+itemp]=HDF5OpenDataSet(partsgroup[i*NHDFTYPE+k],hdf_parts[k]->names[hdf_parts[k]->propindex[HDFGASIMETAL]]);
                      partsdataspaceall[i*NHDFTYPE*NHDFDATABLOCK+k*NHDFDATABLOCK+itemp]=HDF5OpenDataSpace(partsdatasetall[i*NHDFTYPE*NHDFDATABLOCK+k*NHDFDATABLOCK+itemp]);
                    }
                    if (k==HDFSTARTYPE){
                      partsdatasetall[i*NHDFTYPE*NHDFDATABLOCK+k*NHDFDATABLOCK+itemp]=HDF5OpenDataSet(partsgroup[i*NHDFTYPE+k],hdf_parts[k]->names[hdf_parts[k]->propindex[HDFSTARIMETAL]]);
                      partsdataspaceall[i*NHDFTYPE*NHDFDATABLOCK+k*NHDFDATABLOCK+itemp]=HDF5OpenDataSpace(partsdatasetall[i*NHDFTYPE*NHDFDATABLOCK+k*NHDFDATABLOCK+itemp]);
                    }
                  }
                  //then get star formation time, must also adjust so that if tage<0 this is a wind particle in Illustris so change particle type
                  itemp++;
                  for (j=0;j<nusetypes;j++) {
                    k=usetypes[j];
                    if (k==HDFSTARTYPE){
                      if (ThisTask==0 && opt.iverbose>1) cout<<"Opening group "<<hdf_gnames.part_names[k]<<": Data set "<<hdf_parts[k]->names[hdf_parts[k]->propindex[HDFSTARIAGE]]<<endl;
                      partsdatasetall[i*NHDFTYPE*NHDFDATABLOCK+k*NHDFDATABLOCK+itemp]=HDF5OpenDataSet(partsgroup[i*NHDFTYPE+k],hdf_parts[k]->names[hdf_parts[k]->propindex[HDFSTARIAGE]]);
                      partsdataspaceall[i*NHDFTYPE*NHDFDATABLOCK+k*NHDFDATABLOCK+itemp]=HDF5OpenDataSpace(partsdatasetall[i*NHDFTYPE*NHDFDATABLOCK+k*NHDFDATABLOCK+itemp]);
                    }
                  }
                  if (opt.partsearchtype==PSTDARK && opt.iBaryonSearch) for (j=1;j<=nbusetypes;j++) {
                    k=usetypes[j];
                    if (k==HDFSTARTYPE){
                      partsdatasetall[i*NHDFTYPE*NHDFDATABLOCK+k*NHDFDATABLOCK+itemp]=HDF5OpenDataSet(partsgroup[i*NHDFTYPE+k],hdf_parts[k]->names[hdf_parts[k]->propindex[HDFSTARIAGE]]);
                      partsdataspaceall[i*NHDFTYPE*NHDFDATABLOCK+k*NHDFDATABLOCK+itemp]=HDF5OpenDataSpace(partsdatasetall[i*NHDFTYPE*NHDFDATABLOCK+k*NHDFDATABLOCK+itemp]);
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
                    HDF5ReadHyperSlabReal(doublebuff,partsdatasetall[i*NHDFTYPE*NHDFDATABLOCK+k*NHDFDATABLOCK+itemp], partsdataspaceall[i*NHDFTYPE*NHDFDATABLOCK+k*NHDFDATABLOCK+itemp], 1, 3, nchunk, n);
                    //velocities
                    itemp++;
                    HDF5ReadHyperSlabReal(veldoublebuff,partsdatasetall[i*NHDFTYPE*NHDFDATABLOCK+k*NHDFDATABLOCK+itemp], partsdataspaceall[i*NHDFTYPE*NHDFDATABLOCK+k*NHDFDATABLOCK+itemp], 1, 3, nchunk, n);
                    //ids
                    itemp++;
                    HDF5ReadHyperSlabInteger(longbuff,partsdatasetall[i*NHDFTYPE*NHDFDATABLOCK+k*NHDFDATABLOCK+itemp], partsdataspaceall[i*NHDFTYPE*NHDFDATABLOCK+k*NHDFDATABLOCK+itemp], 1, 1, nchunk, n);

                    //masses
                    itemp++;
                    if (hdf_header_info[i].mass[k]==0) {
                        HDF5ReadHyperSlabReal(massdoublebuff,partsdatasetall[i*NHDFTYPE*NHDFDATABLOCK+k*NHDFDATABLOCK+itemp], partsdataspaceall[i*NHDFTYPE*NHDFDATABLOCK+k*NHDFDATABLOCK+itemp], 1, 1, nchunk, n);
                    }
#ifdef GASON
                    //self-energy
                    itemp++;
                    if (k == HDFGASTYPE) {
                        HDF5ReadHyperSlabReal(udoublebuff,partsdatasetall[i*NHDFTYPE*NHDFDATABLOCK+k*NHDFDATABLOCK+itemp], partsdataspaceall[i*NHDFTYPE*NHDFDATABLOCK+k*NHDFDATABLOCK+itemp], 1, 1, nchunk, n);
                    }
#ifdef STARON
                    //star formation rate
                    itemp++;
                    if (k == HDFGASTYPE) {
                        HDF5ReadHyperSlabReal(SFRdoublebuff,partsdatasetall[i*NHDFTYPE*NHDFDATABLOCK+k*NHDFDATABLOCK+itemp], partsdataspaceall[i*NHDFTYPE*NHDFDATABLOCK+k*NHDFDATABLOCK+itemp], 1, 1, nchunk, n);
                    }

                    //metallicity
                    itemp++;
                    if (k == HDFGASTYPE || k == HDFSTARTYPE) {
                        HDF5ReadHyperSlabReal(Zdoublebuff,partsdatasetall[i*NHDFTYPE*NHDFDATABLOCK+k*NHDFDATABLOCK+itemp], partsdataspaceall[i*NHDFTYPE*NHDFDATABLOCK+k*NHDFDATABLOCK+itemp], 1, 1, nchunk, n);
                    }

                    //stellar age
                    itemp++;
                    if (k == HDFSTARTYPE) {
                        HDF5ReadHyperSlabReal(Tagedoublebuff,partsdatasetall[i*NHDFTYPE*NHDFDATABLOCK+k*NHDFDATABLOCK+itemp], partsdataspaceall[i*NHDFTYPE*NHDFDATABLOCK+k*NHDFDATABLOCK+itemp], 1, 1, nchunk, n);
                    }
#endif
#endif

                    for (int nn=0;nn<nchunk;nn++) {
                        ibuf=MPIGetParticlesProcessor(doublebuff[nn*3],doublebuff[nn*3+1],doublebuff[nn*3+2]);
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
                        Pbuf[ibufindex].SetPosition(doublebuff[nn*3],doublebuff[nn*3+1],doublebuff[nn*3+2]);
                        Pbuf[ibufindex].SetVelocity(veldoublebuff[nn*3],veldoublebuff[nn*3+1],veldoublebuff[nn*3+2]);
                        if (hdf_header_info[i].mass[k]==0)Pbuf[ibufindex].SetMass(massdoublebuff[nn]);
                        else Pbuf[ibufindex].SetMass(hdf_header_info[i].mass[k]);
                        Pbuf[ibufindex].SetPID(longbuff[nn]);
                        Pbuf[ibufindex].SetID(nn);
                        if (k==HDFGASTYPE) Pbuf[ibufindex].SetType(GASTYPE);
                        else if (k==HDFDMTYPE) Pbuf[ibufindex].SetType(DARKTYPE);
#ifdef HIGHRES
                        else if (k==HDFDM1TYPE) Pbuf[ibufindex].SetType(DARKTYPE);
                        else if (k==HDFDM2TYPE) Pbuf[ibufindex].SetType(DARKTYPE);
#endif
                        else if (k==HDFSTARTYPE) Pbuf[ibufindex].SetType(STARTYPE);
                        else if (k==HDFBHTYPE) Pbuf[ibufindex].SetType(BHTYPE);

#ifdef HIGHRES
                        if (k==HDFDMTYPE && MP_DM>Pbuf[ibufindex].GetMass()) MP_DM=Pbuf[ibufindex].GetMass();
#endif
#ifdef GASON
                      if (k==HDFGASTYPE) {
                        Pbuf[ibufindex].SetU(udoublebuff[nn]);
#ifdef STARON
                        Pbuf[ibufindex].SetSFR(SFRdoublebuff[nn]);
                        Pbuf[ibufindex].SetZmet(Zdoublebuff[nn]);
#endif
                    }
#endif
#ifdef STARON
                      if (k==HDFSTARTYPE) {
                        Pbuf[ibufindex].SetZmet(Zdoublebuff[nn]);
                        if (Tagedoublebuff[nn]<0) Pbuf[ibufindex].SetType(WINDTYPE);
                        Pbuf[ibufindex].SetTage(Tagedoublebuff[nn]);
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
                      HDF5ReadHyperSlabReal(doublebuff,partsdatasetall[i*NHDFTYPE*NHDFDATABLOCK+k*NHDFDATABLOCK+itemp], partsdataspaceall[i*NHDFTYPE*NHDFDATABLOCK+k*NHDFDATABLOCK+itemp], 1, 3, nchunk, n);
                      //velocities
                      itemp++;
                      HDF5ReadHyperSlabReal(veldoublebuff,partsdatasetall[i*NHDFTYPE*NHDFDATABLOCK+k*NHDFDATABLOCK+itemp], partsdataspaceall[i*NHDFTYPE*NHDFDATABLOCK+k*NHDFDATABLOCK+itemp], 1, 3, nchunk, n);
                      //ids
                      itemp++;
                      HDF5ReadHyperSlabInteger(longbuff,partsdatasetall[i*NHDFTYPE*NHDFDATABLOCK+k*NHDFDATABLOCK+itemp], partsdataspaceall[i*NHDFTYPE*NHDFDATABLOCK+k*NHDFDATABLOCK+itemp], 1, 3, nchunk, n);
                      //masses
                      itemp++;
                      if (hdf_header_info[i].mass[k]==0) {
                          HDF5ReadHyperSlabReal(massdoublebuff,partsdatasetall[i*NHDFTYPE*NHDFDATABLOCK+k*NHDFDATABLOCK+itemp], partsdataspaceall[i*NHDFTYPE*NHDFDATABLOCK+k*NHDFDATABLOCK+itemp], 1, 1, nchunk, n);
                      }
#ifdef GASON
                      //self-energy
                      itemp++;
                      if (k == HDFGASTYPE) {
                          HDF5ReadHyperSlabReal(udoublebuff,partsdatasetall[i*NHDFTYPE*NHDFDATABLOCK+k*NHDFDATABLOCK+itemp], partsdataspaceall[i*NHDFTYPE*NHDFDATABLOCK+k*NHDFDATABLOCK+itemp], 1, 1, nchunk, n);
                      }
#ifdef STARON
                      //star formation rate
                      itemp++;
                      if (k == HDFGASTYPE) {
                          HDF5ReadHyperSlabReal(SFRdoublebuff,partsdatasetall[i*NHDFTYPE*NHDFDATABLOCK+k*NHDFDATABLOCK+itemp], partsdataspaceall[i*NHDFTYPE*NHDFDATABLOCK+k*NHDFDATABLOCK+itemp], 1, 1, nchunk, n);
                      }

                      //metallicity
                      itemp++;
                      if (k == HDFGASTYPE || k == HDFSTARTYPE) {
                          HDF5ReadHyperSlabReal(Zdoublebuff,partsdatasetall[i*NHDFTYPE*NHDFDATABLOCK+k*NHDFDATABLOCK+itemp], partsdataspaceall[i*NHDFTYPE*NHDFDATABLOCK+k*NHDFDATABLOCK+itemp], 1, 1, nchunk, n);
                      }

                      //stellar age
                      itemp++;
                      if (k == HDFSTARTYPE) {
                          HDF5ReadHyperSlabReal(Tagedoublebuff,partsdatasetall[i*NHDFTYPE*NHDFDATABLOCK+k*NHDFDATABLOCK+itemp], partsdataspaceall[i*NHDFTYPE*NHDFDATABLOCK+k*NHDFDATABLOCK+itemp], 1, 1, nchunk, n);
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
                        Pbuf[ibufindex].SetPosition(doublebuff[nn*3],doublebuff[nn*3+1],doublebuff[nn*3+2]);
                        Pbuf[ibufindex].SetVelocity(veldoublebuff[nn*3],veldoublebuff[nn*3+1],veldoublebuff[nn*3+2]);
                        if (hdf_header_info[i].mass[k]==0)Pbuf[ibufindex].SetMass(massdoublebuff[nn]);
                        else Pbuf[ibufindex].SetMass(hdf_header_info[i].mass[k]);
                        Pbuf[ibufindex].SetPID(longbuff[nn]);
                        Pbuf[ibufindex].SetID(nn);
                        if (k==HDFGASTYPE) Pbuf[ibufindex].SetType(GASTYPE);
                        else if (k==HDFDMTYPE) Pbuf[ibufindex].SetType(DARKTYPE);
#ifdef HIGHRES
                        else if (k==HDFDM1TYPE) Pbuf[ibufindex].SetType(DARKTYPE);
                        else if (k==HDFDM2TYPE) Pbuf[ibufindex].SetType(DARKTYPE);
#endif
                        else if (k==HDFSTARTYPE) Pbuf[ibufindex].SetType(STARTYPE);
                        else if (k==HDFBHTYPE) Pbuf[ibufindex].SetType(BHTYPE);
#ifdef HIGHRES
                        if (k==HDFDMTYPE && MP_DM>Pbuf[ibufindex].GetMass()) MP_DM=Pbuf[ibufindex].GetMass();
#endif
#ifdef GASON
                        if (k==HDFGASTYPE) {
                          Pbuf[ibufindex].SetU(udoublebuff[nn]);
#ifdef STARON
                          Pbuf[ibufindex].SetSFR(SFRdoublebuff[nn]);
                          Pbuf[ibufindex].SetZmet(Zdoublebuff[nn]);
#endif
                        }
#endif
#ifdef STARON
                        if (k==HDFSTARTYPE) {
                          Pbuf[ibufindex].SetZmet(Zdoublebuff[nn]);
                          if (Tagedoublebuff[nn]<0) Pbuf[ibufindex].SetType(WINDTYPE);
                          Pbuf[ibufindex].SetTage(Tagedoublebuff[nn]);
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
                  }//end of part type
                }//end of baryon if
                //close data spaces
                for (auto &hidval:partsdataspace) HDF5CloseDataSpace(hidval);
                for (auto &hidval:partsdataset) HDF5CloseDataSet(hidval);
                for (auto &hidval:partsgroup) HDF5CloseGroup(hidval);
            }//end of try block
            /*
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
            */
            HDF5CloseFile(Fhdf[i]);
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
        if (opt.comove) opt.p*=opt.lengthinputconversion;
        else opt.p*=opt.lengthinputconversion*opt.a;
      }
      else {
        //adjust period
        if (opt.comove) opt.p*=opt.lengthinputconversion/opt.h;
        else opt.p*=opt.lengthinputconversion/opt.h*opt.a;
      }

#ifdef USEMPI
    }
#endif
#ifdef USEMPI
    if (opt.nsnapread>1) {
        MPI_Allreduce(&MP_DM,&MP_DM, 1, MPI_DOUBLE, MPI_MIN,mpi_comm_read);
    }
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
    opt.zoomlowmassdm=MP_DM*mscale;
    if (opt.Neff==-1) {
      //Once smallest mass particle is found (which should correspond to highest resolution area,
      if (opt.Omega_b==0) MP_B=0;
      LN=pow((MP_DM+MP_B)*mscale/opt.rhobg,1.0/3.0);
    }
    else {
      LN=opt.p/(Double_t)opt.Neff;
    }
    #ifdef USEMPI
    MPI_Bcast(&opt.zoomlowmassdm,sizeof(opt.zoomlowmassdm),MPI_BYTE,0,MPI_COMM_WORLD);
    #endif
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
    if(opt.ihdfnameconvention == HDFSWIFTEAGLENAMES) vscale = opt.velocityinputconversion;
    else vscale = opt.velocityinputconversion*sqrt(opt.a);

    //finally adjust to appropriate units
    for (i=0;i<Nlocal;i++)
    {
      Part[i].SetMass(Part[i].GetMass()*mscale);
      for (int j=0;j<3;j++) Part[i].SetVelocity(j,Part[i].GetVelocity(j)*vscale+Hubbleflow*Part[i].GetPosition(j));
      for (int j=0;j<3;j++) Part[i].SetPosition(j,Part[i].GetPosition(j)*lscale);
#ifdef GASON
      if (Part[i].GetType()==GASTYPE) Part[i].SetU(Part[i].GetU()*opt.velocityinputconversion*opt.velocityinputconversion);
#endif
    }
    if (Pbaryons!=NULL && opt.iBaryonSearch==1) {
      for (i=0;i<Nlocalbaryon[0];i++)
      {
        Pbaryons[i].SetMass(Pbaryons[i].GetMass()*mscale);
        for (int j=0;j<3;j++) Pbaryons[i].SetVelocity(j,Pbaryons[i].GetVelocity(j)*vscale+Hubbleflow*Pbaryons[i].GetPosition(j));
        for (int j=0;j<3;j++) Pbaryons[i].SetPosition(j,Pbaryons[i].GetPosition(j)*lscale);
#ifdef GASON
        Pbaryons[i].SetU(Pbaryons[i].GetU()*opt.velocityinputconversion*opt.velocityinputconversion);
#endif
    }
    }
#endif


    delete[] intbuff;
    delete[] longbuff;
    delete[] uintbuff;
    delete[] floatbuff;
    delete[] doublebuff;
#ifdef USEMPI
    delete[] velfloatbuff;
    delete[] veldoublebuff;
    delete[] massfloatbuff;
    delete[] massdoublebuff;
#ifdef GASON
    delete[] ufloatbuff;
    delete[] udoublebuff;
#endif
#if defined(GASON)&&defined(STARON)
    delete[] Zfloatbuff;
    delete[] Zdoublebuff;
    delete[] SFRfloatbuff;
    delete[] SFRdoublebuff;
#endif
#ifdef STARON
    delete[] Tagefloatbuff;
    delete[] Tagedoublebuff;
#endif
#endif

}

#endif
